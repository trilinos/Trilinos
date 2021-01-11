// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid2_IntegrationToolsDef.hpp
    \brief  Definition file for the Intrepid2::IntegrationTools class.
    \author Created by Nathan V. Roberts.
*/

#ifndef __INTREPID2_INTEGRATIONTOOLS_DEF_HPP__
#define __INTREPID2_INTEGRATIONTOOLS_DEF_HPP__

#include "Intrepid2_FunctorIterator.hpp"
#include "Intrepid2_TensorArgumentIterator.hpp"

#include "Teuchos_TimeMonitor.hpp"

namespace Intrepid2 {

  namespace Impl
  {
    /**
      \brief Implementation of a general sum factorization algorithm, abstracted from the algorithm described by Mora and Demkowicz, for integration.  Uses hierarchical parallelism.
     */
    template<class Scalar, class ExecSpaceType, int integralViewRank>
    class F_Integrate
    {
      using TeamPolicy = Kokkos::TeamPolicy<ExecSpaceType>;
      using TeamMember = typename TeamPolicy::member_type;
      
      using IntegralViewType = Kokkos::View<typename RankExpander<Scalar, integralViewRank>::value_type, ExecSpaceType>;
      IntegralViewType integralView_;
      TensorData<Scalar,ExecSpaceType> leftComponent_;
      Data<Scalar,ExecSpaceType> composedTransform_;
      TensorData<Scalar,ExecSpaceType> rightComponent_;
      TensorData<Scalar,ExecSpaceType> cellMeasures_;
      int a_offset_;
      int b_offset_;
      int leftComponentSpan_;   //  leftComponentSpan tracks the dimensions spanned by the left component
      int rightComponentSpan_;  // rightComponentSpan tracks the dimensions spanned by the right component
      int numTensorComponents_;
      int  leftFieldOrdinalOffset_;
      int rightFieldOrdinalOffset_;
      bool forceNonSpecialized_; // if true, don't use the specialized (more manual) implementation(s) for particular component counts.  Primary use case is for testing.
      
      size_t fad_size_output_ = 0; // 0 if not a fad type
      
      Kokkos::Array<int, 7> offsetsForComponentOrdinal_;
      
      // as an optimization, we do all the bounds and argument iteration within the functor rather than relying on TensorArgumentIterator
      // (this also makes it easier to reorder loops, etc., for further optimizations)
      Kokkos::Array<int,Parameters::MaxTensorComponents>  leftFieldBounds_;
      Kokkos::Array<int,Parameters::MaxTensorComponents> rightFieldBounds_;
      Kokkos::Array<int,Parameters::MaxTensorComponents> pointBounds_;
      
      Kokkos::Array<int,Parameters::MaxTensorComponents>  leftFieldRelativeEnumerationSpans_; // total number of enumeration indices with arguments prior to the startingComponent fixed
      Kokkos::Array<int,Parameters::MaxTensorComponents> rightFieldRelativeEnumerationSpans_;
      
      int maxFieldsLeft_;
      int maxFieldsRight_;
      int maxPointCount_;
    public:
      F_Integrate(Data<Scalar,ExecSpaceType> integralData,
                  TensorData<Scalar,ExecSpaceType> leftComponent,
                  Data<Scalar,ExecSpaceType> composedTransform,
                  TensorData<Scalar,ExecSpaceType> rightComponent,
                  TensorData<Scalar,ExecSpaceType> cellMeasures,
                  int a_offset,
                  int b_offset,
                  int leftFieldOrdinalOffset,
                  int rightFieldOrdinalOffset,
                  bool forceNonSpecialized)
      :
      integralView_(integralData.template getUnderlyingView<integralViewRank>()),
      leftComponent_(leftComponent),
      composedTransform_(composedTransform),
      rightComponent_(rightComponent),
      cellMeasures_(cellMeasures),
      a_offset_(a_offset),
      b_offset_(b_offset),
      leftComponentSpan_(leftComponent.extent_int(2)),
      rightComponentSpan_(rightComponent.extent_int(2)),
      numTensorComponents_(leftComponent.numTensorComponents()),
      leftFieldOrdinalOffset_(leftFieldOrdinalOffset),
      rightFieldOrdinalOffset_(rightFieldOrdinalOffset),
      forceNonSpecialized_(forceNonSpecialized)
      {
        INTREPID2_TEST_FOR_EXCEPTION(numTensorComponents_ != rightComponent_.numTensorComponents(), std::invalid_argument, "Left and right components must have matching number of tensorial components");
        
        // set up bounds containers
        const int FIELD_DIM = 0;
        const int POINT_DIM = 1;
        maxFieldsLeft_  = 0;
        maxFieldsRight_ = 0;
        maxPointCount_  = 0;
        for (int r=0; r<numTensorComponents_; r++)
        {
          leftFieldBounds_[r]  = leftComponent_.getTensorComponent(r).extent_int(FIELD_DIM);
          maxFieldsLeft_       = std::max(maxFieldsLeft_, leftFieldBounds_[r]);
          rightFieldBounds_[r] = rightComponent_.getTensorComponent(r).extent_int(FIELD_DIM);
          maxFieldsRight_      = std::max(maxFieldsRight_, rightFieldBounds_[r]);
          pointBounds_[r]      = leftComponent_.getTensorComponent(r).extent_int(POINT_DIM);
          maxPointCount_       = std::max(maxPointCount_, pointBounds_[r]);
        }
        
        // set up relative enumeration spans: total number of enumeration indices with arguments prior to the startingComponent fixed.  These are for *truncated* iterators; hence the -2 rather than -1 for the first startingComponent value.
        int  leftRelativeEnumerationSpan = 1;
        int rightRelativeEnumerationSpan = 1;
        for (int startingComponent=numTensorComponents_-2; startingComponent>=0; startingComponent--)
        {
          leftRelativeEnumerationSpan  *=  leftFieldBounds_[startingComponent];
          rightRelativeEnumerationSpan *= rightFieldBounds_[startingComponent];
          leftFieldRelativeEnumerationSpans_ [startingComponent] =  leftRelativeEnumerationSpan;
          rightFieldRelativeEnumerationSpans_[startingComponent] = rightRelativeEnumerationSpan;
        }
        
        // prepare for allocation of temporary storage
        // note: tempStorage goes "backward", starting from the final component, which needs just one entry
        
        const bool allocateFadStorage = !std::is_pod<Scalar>::value;
        if (allocateFadStorage)
        {
          fad_size_output_ = dimension_scalar(integralView_);
        }
        
        const int R = numTensorComponents_ - 1;
        
        int num_ij = 1; // this counts how many entries there are corresponding to components from r to R-1.
        int allocationSoFar = 0;
        offsetsForComponentOrdinal_[R] = allocationSoFar;
        allocationSoFar++; // we store one entry corresponding to R, the final component
        
        for (int r=R-1; r>0; r--) // last component is innermost in the for loops (requires least storage)
        {
          const int leftFields  =  leftComponent.getTensorComponent(r).extent_int(0);
          const int rightFields = rightComponent.getTensorComponent(r).extent_int(0);
          
          num_ij *= leftFields * rightFields;
          
          offsetsForComponentOrdinal_[r] = allocationSoFar;
          allocationSoFar += num_ij;
        }
        offsetsForComponentOrdinal_[0] = allocationSoFar; // first component stores directly to final integralView.
      }
      
      template<size_t maxComponents, size_t numComponents = maxComponents>
      KOKKOS_INLINE_FUNCTION
      int incrementArgument(      Kokkos::Array<int,maxComponents> &arguments,
                            const Kokkos::Array<int,maxComponents> &bounds) const
      {
        if (numComponents == 0)
        {
          return -1;
        }
        else
        {
          int r = static_cast<int>(numComponents - 1);
          while (arguments[r] + 1 >= bounds[r])
          {
            arguments[r] = 0; // reset
            r--;
            if (r < 0) break;
          }
          if (r >= 0) ++arguments[r];
          return r;
        }
      }
      
      //! runtime-sized variant of incrementArgument; gets used by approximate flop count.
      KOKKOS_INLINE_FUNCTION
      int incrementArgument(      Kokkos::Array<int,Parameters::MaxTensorComponents> &arguments,
                            const Kokkos::Array<int,Parameters::MaxTensorComponents> &bounds,
                            const int &numComponents) const
      {
        if (numComponents == 0) return -1;
        int r = static_cast<int>(numComponents - 1);
        while (arguments[r] + 1 >= bounds[r])
        {
          arguments[r] = 0; // reset
          r--;
          if (r < 0) break;
        }
        if (r >= 0) ++arguments[r];
        return r;
      }
      
      template<size_t maxComponents, size_t numComponents = maxComponents>
      KOKKOS_INLINE_FUNCTION
      int nextIncrementResult(const Kokkos::Array<int,maxComponents> &arguments,
                              const Kokkos::Array<int,maxComponents> &bounds) const
      {
        if (numComponents == 0)
        {
          return -1;
        }
        else
        {
          int r = static_cast<int>(numComponents - 1);
          while (arguments[r] + 1 >= bounds[r])
          {
            r--;
            if (r < 0) break;
          }
          return r;
        }
      }
      
      //! runtime-sized variant of nextIncrementResult; gets used by approximate flop count.
      KOKKOS_INLINE_FUNCTION
      int nextIncrementResult(const Kokkos::Array<int,Parameters::MaxTensorComponents> &arguments,
                              const Kokkos::Array<int,Parameters::MaxTensorComponents> &bounds,
                              const int &numComponents) const
      {
        if (numComponents == 0) return -1;
        int r = numComponents - 1;
        while (arguments[r] + 1 >= bounds[r])
        {
          r--;
          if (r < 0) break;
        }
        return r;
      }
      
      template<size_t maxComponents, size_t numComponents = maxComponents>
      KOKKOS_INLINE_FUNCTION
      int relativeEnumerationIndex(const Kokkos::Array<int,maxComponents> &arguments,
                                   const Kokkos::Array<int,maxComponents> &bounds,
                                   const int startIndex) const
      {
        // the following mirrors what is done in TensorData
        if (numComponents == 0)
        {
          return 0;
        }
        else
        {
          int enumerationIndex = 0;
          for (size_t r=numComponents-1; r>static_cast<size_t>(startIndex); r--)
          {
            enumerationIndex += arguments[r];
            enumerationIndex *= bounds[r-1];
          }
          enumerationIndex += arguments[startIndex];
          return enumerationIndex;
        }
      }
      
      //! runSpecialized implementations are hand-coded variants of run() for a particular number of components.  To allow comparisons with the generic implementation (both in terms of performance and for verification), we use the member variable forceNonSpecialized_ to determine whether runSpecialized is selected when a specialized implementation is available.
      
      // nvcc refuses to compile the below with error, "explicit specialization is not allowed in the current scope".  Clang is OK with it.  We just do a non-templated version below.
//      template<size_t numTensorComponents>
//      KOKKOS_INLINE_FUNCTION
//      void runSpecialized( const TeamMember & teamMember ) const;
      
//      template<>
//      KOKKOS_INLINE_FUNCTION
//      void runSpecialized<3>( const TeamMember & teamMember ) const
      KOKKOS_INLINE_FUNCTION
      void runSpecialized3( const TeamMember & teamMember ) const
      {
        constexpr int numTensorComponents = 3;
        
        Kokkos::Array<int,numTensorComponents>  pointBounds;
        Kokkos::Array<int,numTensorComponents>  leftFieldBounds;
        Kokkos::Array<int,numTensorComponents> rightFieldBounds;
        for (unsigned r=0; r<numTensorComponents; r++)
        {
          pointBounds[r] = pointBounds_[r];
          leftFieldBounds[r]  =  leftFieldBounds_[r];
          rightFieldBounds[r] = rightFieldBounds_[r];
        }
        
        const int cellDataOrdinal = teamMember.league_rank();
        const int numThreads      = teamMember.team_size(); // num threads
        const int scratchViewSize = offsetsForComponentOrdinal_[0]; // per thread
        
        Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged> scratchView; // for caching partial integration values
        Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged> pointWeights; // indexed by (expanded) point; stores M_ab * cell measure
        Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged> leftFields_x, leftFields_y, leftFields_z, rightFields_x, rightFields_y, rightFields_z; // cache the field values for faster access
        if (fad_size_output_ > 0) {
          scratchView = Kokkos::View<Scalar*,   ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(),    scratchViewSize * numThreads,      fad_size_output_);
          pointWeights   = Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged>  (teamMember.team_shmem(), composedTransform_.extent_int(1),  fad_size_output_);
          leftFields_x  = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(),  leftFieldBounds[0], pointBounds[0], fad_size_output_);
          rightFields_x = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), rightFieldBounds[0], pointBounds[0], fad_size_output_);
          leftFields_y  = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(),  leftFieldBounds[1], pointBounds[1], fad_size_output_);
          rightFields_y = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), rightFieldBounds[1], pointBounds[1], fad_size_output_);
          leftFields_z  = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(),  leftFieldBounds[2], pointBounds[2], fad_size_output_);
          rightFields_z = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), rightFieldBounds[2], pointBounds[2], fad_size_output_);
        }
        else {
          scratchView = Kokkos::View<Scalar*,   ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(),    scratchViewSize * numThreads);
          pointWeights   = Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged>  (teamMember.team_shmem(), composedTransform_.extent_int(1));
          leftFields_x  = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(),  leftFieldBounds[0], pointBounds[0]);
          rightFields_x = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), rightFieldBounds[0], pointBounds[0]);
          leftFields_y  = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(),  leftFieldBounds[1], pointBounds[1]);
          rightFields_y = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), rightFieldBounds[1], pointBounds[1]);
          leftFields_z  = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(),  leftFieldBounds[2], pointBounds[2]);
          rightFields_z = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), rightFieldBounds[2], pointBounds[2]);
        }
        
//        int approximateFlopCount = 0;
//        int flopsPerCellMeasuresAccess = cellMeasures_.numTensorComponents() - 1;
        
        constexpr int R = numTensorComponents - 1;
        
        for (int a_component=0; a_component < leftComponentSpan_; a_component++)
        {
          const int a = a_offset_ + a_component;
          for (int b_component=0; b_component < rightComponentSpan_; b_component++)
          {
            const int b = b_offset_ + b_component;
            
            const Data<Scalar,ExecSpaceType> &  leftFinalComponent =  leftComponent_.getTensorComponent(R);
            const Data<Scalar,ExecSpaceType> & rightFinalComponent = rightComponent_.getTensorComponent(R);
            
            const int numLeftFieldsFinal  =  leftFinalComponent.extent_int(0); // shape (F,P[,D])
            const int numRightFieldsFinal = rightFinalComponent.extent_int(0); // shape (F,P[,D])
            
            const Data<Scalar,ExecSpaceType> &  leftTensorComponent_x =  leftComponent_.getTensorComponent(0);
            const Data<Scalar,ExecSpaceType> & rightTensorComponent_x = rightComponent_.getTensorComponent(0);
            const Data<Scalar,ExecSpaceType> &  leftTensorComponent_y =  leftComponent_.getTensorComponent(1);
            const Data<Scalar,ExecSpaceType> & rightTensorComponent_y = rightComponent_.getTensorComponent(1);
            const Data<Scalar,ExecSpaceType> &  leftTensorComponent_z =  leftComponent_.getTensorComponent(2);
            const Data<Scalar,ExecSpaceType> & rightTensorComponent_z = rightComponent_.getTensorComponent(2);
            
            const int maxFields = (maxFieldsLeft_ > maxFieldsRight_) ? maxFieldsLeft_ : maxFieldsRight_;
            Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,maxFields * maxPointCount_), [&] (const int& fieldOrdinalPointOrdinal) {
              const int fieldOrdinal = fieldOrdinalPointOrdinal % maxFields;
              const int pointOrdinal = fieldOrdinalPointOrdinal / maxFields;
              if ((fieldOrdinal < leftTensorComponent_x.extent_int(0)) && (pointOrdinal < leftTensorComponent_x.extent_int(1)))
              {
                const int   leftRank = leftTensorComponent_x.rank();
                leftFields_x(fieldOrdinal,pointOrdinal) = (leftRank == 2) ? leftTensorComponent_x(fieldOrdinal,pointOrdinal) : leftTensorComponent_x(fieldOrdinal,pointOrdinal,a_component);
              }
              if ((fieldOrdinal < leftTensorComponent_y.extent_int(0)) && (pointOrdinal < leftTensorComponent_y.extent_int(1)))
              {
                const int   leftRank = leftTensorComponent_y.rank();
                leftFields_y(fieldOrdinal,pointOrdinal) = (leftRank == 2) ? leftTensorComponent_y(fieldOrdinal,pointOrdinal) : leftTensorComponent_y(fieldOrdinal,pointOrdinal,a_component);
              }
              if ((fieldOrdinal < leftTensorComponent_z.extent_int(0)) && (pointOrdinal < leftTensorComponent_z.extent_int(1)))
              {
                const int   leftRank = leftTensorComponent_z.rank();
                leftFields_z(fieldOrdinal,pointOrdinal) = (leftRank == 2) ? leftTensorComponent_z(fieldOrdinal,pointOrdinal) : leftTensorComponent_z(fieldOrdinal,pointOrdinal,a_component);
              }
              if ((fieldOrdinal < rightTensorComponent_x.extent_int(0)) && (pointOrdinal < rightTensorComponent_x.extent_int(1)))
              {
                const int   rightRank = rightTensorComponent_x.rank();
                rightFields_x(fieldOrdinal,pointOrdinal) = (rightRank == 2) ? rightTensorComponent_x(fieldOrdinal,pointOrdinal) : rightTensorComponent_x(fieldOrdinal,pointOrdinal,b_component);
              }
              if ((fieldOrdinal < rightTensorComponent_y.extent_int(0)) && (pointOrdinal < rightTensorComponent_y.extent_int(1)))
              {
                const int   rightRank = rightTensorComponent_y.rank();
                rightFields_y(fieldOrdinal,pointOrdinal) = (rightRank == 2) ? rightTensorComponent_y(fieldOrdinal,pointOrdinal) : rightTensorComponent_y(fieldOrdinal,pointOrdinal,b_component);
              }
              if ((fieldOrdinal < rightTensorComponent_z.extent_int(0)) && (pointOrdinal < rightTensorComponent_z.extent_int(1)))
              {
                const int   rightRank = rightTensorComponent_z.rank();
                rightFields_z(fieldOrdinal,pointOrdinal) = (rightRank == 2) ? rightTensorComponent_z(fieldOrdinal,pointOrdinal) : rightTensorComponent_z(fieldOrdinal,pointOrdinal,b_component);
              }
            });
            
            if (composedTransform_.underlyingMatchesNotional())
            {
              const auto & composedTransformView = composedTransform_.getUnderlyingView4();
              Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,composedTransformView.extent_int(1)), [&] (const int& pointOrdinal) {
                pointWeights(pointOrdinal) = composedTransformView(cellDataOrdinal,pointOrdinal,a,b) * cellMeasures_(cellDataOrdinal,pointOrdinal);
              });
            }
            else
            {
              Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,composedTransform_.extent_int(1)), [&] (const int& pointOrdinal) {
                pointWeights(pointOrdinal) = composedTransform_(cellDataOrdinal,pointOrdinal,a,b) * cellMeasures_(cellDataOrdinal,pointOrdinal);
              });
            }
            
            // approximateFlopCount += composedTransform_.extent_int(1) * cellMeasures.numTensorComponents(); // cellMeasures does numTensorComponents - 1 multiplies on each access
            
            // synchronize threads
            teamMember.team_barrier();
            
            // Setting scratchView to 0 here is not necessary from an algorithmic point of view, but *might* help with performance (due to a first-touch policy)
            const int scratchOffsetForThread = teamMember.team_rank() * scratchViewSize;
            for (int i=scratchOffsetForThread; i<scratchOffsetForThread+scratchViewSize; i++)
            {
              scratchView(i) = 0.0;
            }
            
            // TODO: consider adding an innerLoopRange that is sized to be the maximum of the size of the inner loops we'd like to parallelize over.  (note that this means we do the work in the outer loop redundantly that many times…)
            Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numLeftFieldsFinal * numRightFieldsFinal), [&] (const int& leftRightFieldEnumeration) {
              const int iz = leftRightFieldEnumeration % numLeftFieldsFinal;
              const int jz = leftRightFieldEnumeration / numLeftFieldsFinal;
              
              // as an optimization, we are using the same argument array for both the truncated field that spans s to R-1 and the full one that goes to R
              // the "R" argument here is ignored by the methods that treat the truncated field (it's beyond their bounds)
              Kokkos::Array<int,numTensorComponents>  leftFieldArguments;
              Kokkos::Array<int,numTensorComponents> rightFieldArguments;
              rightFieldArguments[R] = jz;
              leftFieldArguments[R]  = iz;
              
              Kokkos::Array<int,numTensorComponents> pointArguments;
              for (int i=0; i<numTensorComponents; i++)
              {
                pointArguments[i] = 0;
              }
              
              for (int lx=0; lx<pointBounds[0]; lx++)
              {
                pointArguments[0] = lx;
                
                // clear Gy scratch:
                // in scratch, Gz (1 entry) comes first, then Gy entries.
                const int Gy_start_index = scratchOffsetForThread + offsetsForComponentOrdinal_[1];
                const int Gy_end_index   = scratchOffsetForThread + offsetsForComponentOrdinal_[0];
                
                for (int Gy_index=Gy_start_index; Gy_index < Gy_end_index; Gy_index++)
                {
                  scratchView(Gy_index) = 0;
                }
                
                for (int ly=0; ly<pointBounds[1]; ly++)
                {
                  pointArguments[1] = ly;
                  
                  Scalar * Gz = &scratchView(scratchOffsetForThread + offsetsForComponentOrdinal_[R]);
                  *Gz = 0;
                  
                  for (int lz=0; lz < pointBounds[2]; lz++)
                  {
                    const Scalar &  leftValue =  leftFields_z(iz,lz);
                    const Scalar & rightValue = rightFields_z(jz,lz);
                    
                    pointArguments[2] = lz;
                    int pointEnumerationIndex = relativeEnumerationIndex<numTensorComponents>(pointArguments, pointBounds, 0);
                    
                    *Gz += leftValue * pointWeights(pointEnumerationIndex) * rightValue;
                    
  //                approximateFlopCount += 3; // 2 multiplies, 1 sum
                  } // lz
                  
                  for (int iy=0; iy<leftFieldBounds_[1]; iy++)
                  {
                    leftFieldArguments[1] = iy;
                    const int leftEnumerationIndex_y = relativeEnumerationIndex<numTensorComponents,R>(leftFieldArguments, leftFieldBounds, 1);
                    
                    const Scalar & leftValue = leftFields_y(iy,ly);
                    
                    for (int jy=0; jy<rightFieldBounds_[1]; jy++)
                    {
                      rightFieldArguments[1] = jy;
                      
                      const int rightEnumerationIndex_y = relativeEnumerationIndex<numTensorComponents,R>(rightFieldArguments, rightFieldBounds, 1);
                      const Scalar & rightValue = rightFields_y(jy,ly);
                      
                      const int & rightEnumerationSpan_y = rightFieldRelativeEnumerationSpans_[1];
                      const int Gy_index = leftEnumerationIndex_y  * rightEnumerationSpan_y + rightEnumerationIndex_y;
                      
                      const int Gz_index = 0;
                      const Scalar & Gz = scratchView(scratchOffsetForThread + offsetsForComponentOrdinal_[2] + Gz_index);
                      
                      Scalar & Gy = scratchView(scratchOffsetForThread + offsetsForComponentOrdinal_[1] + Gy_index);
                                            
                      Gy += leftValue * Gz * rightValue;
    //                approximateFlopCount += 3; // 2 multiplies, 1 sum
                    }
                  }
                } // ly
                for (int ix=0; ix<leftFieldBounds_[0]; ix++)
                {
                  leftFieldArguments[0] = ix;
                  const Scalar & leftValue = leftFields_x(ix,lx);
                  
                  for (int iy=0; iy<leftFieldBounds_[1]; iy++)
                  {
                    leftFieldArguments[1] = iy;
                    
                    const int leftEnumerationIndex_y = relativeEnumerationIndex<numTensorComponents,R>(leftFieldArguments, leftFieldBounds, 1);
                    
                    for (int jx=0; jx<rightFieldBounds_[0]; jx++)
                    {
                      rightFieldArguments[0] = jx;
                      const Scalar & rightValue = rightFields_x(jx,lx);
                      
                      for (int jy=0; jy<rightFieldBounds_[1]; jy++)
                      {
                        rightFieldArguments[1] = jy;
                        const int rightEnumerationIndex_y = relativeEnumerationIndex<numTensorComponents,R>(rightFieldArguments, rightFieldBounds, 1);
                        
                        const int rightEnumerationSpan_y = rightFieldRelativeEnumerationSpans_[1];
                        
                        const int Gy_index = leftEnumerationIndex_y  * rightEnumerationSpan_y + rightEnumerationIndex_y;
                        Scalar & Gy = scratchView(scratchOffsetForThread + offsetsForComponentOrdinal_[1] + Gy_index);
                        
                        // compute enumeration indices to get field indices into output view
                        int  leftEnumerationIndex = relativeEnumerationIndex<numTensorComponents>( leftFieldArguments,  leftFieldBounds, 0);
                        int rightEnumerationIndex = relativeEnumerationIndex<numTensorComponents>(rightFieldArguments, rightFieldBounds, 0);
                        const int  leftFieldIndex =  leftEnumerationIndex +  leftFieldOrdinalOffset_;
                        const int rightFieldIndex = rightEnumerationIndex + rightFieldOrdinalOffset_;

                        if (integralViewRank == 3)
                        {
//                          if ((leftFieldIndex==0) && (rightFieldIndex==2))
//                          {
//                            using std::cout;
//                            using std::endl;
//                            cout << "***** Contribution to (0,0,2) *****\n";
//                            cout << "lx = " << lx << endl;
//                            cout << "ix = " << ix << endl;
//                            cout << "iy = " << iy << endl;
//                            cout << "jx = " << jx << endl;
//                            cout << "jy = " << jy << endl;
//                            cout << "iz = " << iz << endl;
//                            cout << "jz = " << jz << endl;
//                            cout << " leftValue = " << leftValue << endl;
//                            cout << "rightValue = " << rightValue << endl;
//                            cout << "Gy = " <<  Gy << endl;
//
//                            cout << endl;
//                          }
                          
                          // shape (C,F1,F2)
                          integralView_.access(cellDataOrdinal,leftFieldIndex,rightFieldIndex) += leftValue * Gy * rightValue;
                        }
                        else
                        {
                          // shape (F1,F2)
                          integralView_.access(leftFieldIndex,rightFieldIndex,0) += leftValue * Gy * rightValue;
                        }
      //                approximateFlopCount += 3; // 2 multiplies, 1 sum
                      } // jy
                    } // ix
                  } // iy
                } // ix
              } // lx
            }); // TeamThreadRange parallel_for - (iz,jz) loop
          }
        }
//        std::cout << "flop count per cell (within operator()) : " << approximateFlopCount << std::endl;
      }
      
      template<size_t numTensorComponents>
      KOKKOS_INLINE_FUNCTION
      void run( const TeamMember & teamMember ) const
      {
        Kokkos::Array<int,numTensorComponents>  pointBounds;
        Kokkos::Array<int,numTensorComponents>  leftFieldBounds;
        Kokkos::Array<int,numTensorComponents> rightFieldBounds;
        for (unsigned r=0; r<numTensorComponents; r++)
        {
          pointBounds[r] = pointBounds_[r];
          leftFieldBounds[r]  =  leftFieldBounds_[r];
          rightFieldBounds[r] = rightFieldBounds_[r];
        }
        
        const int cellDataOrdinal = teamMember.league_rank();
        const int numThreads      = teamMember.team_size(); // num threads
        const int scratchViewSize = offsetsForComponentOrdinal_[0]; // per thread
        Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged> scratchView; // for caching partial integration values
        Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged> pointWeights; // indexed by (expanded) point; stores M_ab * cell measure
        Kokkos::View<Scalar***, ExecSpaceType, Kokkos::MemoryUnmanaged> leftFields, rightFields; // cache the field values for faster access
        if (fad_size_output_ > 0) {
          scratchView = Kokkos::View<Scalar*,   ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), scratchViewSize * numThreads,     fad_size_output_);
          pointWeights   = Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged>  (teamMember.team_shmem(), composedTransform_.extent_int(1), fad_size_output_);
          leftFields  = Kokkos::View<Scalar***, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), numTensorComponents, maxFieldsLeft_,  maxPointCount_, fad_size_output_);
          rightFields = Kokkos::View<Scalar***, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), numTensorComponents, maxFieldsRight_, maxPointCount_, fad_size_output_);
        }
        else {
          scratchView = Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), scratchViewSize*numThreads);
          pointWeights   = Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged>  (teamMember.team_shmem(), composedTransform_.extent_int(1));
          leftFields  = Kokkos::View<Scalar***, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), numTensorComponents, maxFieldsLeft_,  maxPointCount_);
          rightFields = Kokkos::View<Scalar***, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), numTensorComponents, maxFieldsRight_, maxPointCount_);
        }
        
//        int approximateFlopCount = 0;
//        int flopsPerCellMeasuresAccess = cellMeasures_.numTensorComponents() - 1;
        
        constexpr int R = numTensorComponents - 1;
        
        for (int a_component=0; a_component < leftComponentSpan_; a_component++)
        {
          const int a = a_offset_ + a_component;
          for (int b_component=0; b_component < rightComponentSpan_; b_component++)
          {
            const int b = b_offset_ + b_component;
            
            const Data<Scalar,ExecSpaceType> &  leftFinalComponent =  leftComponent_.getTensorComponent(R);
            const Data<Scalar,ExecSpaceType> & rightFinalComponent = rightComponent_.getTensorComponent(R);
            
            const int numLeftFieldsFinal  =  leftFinalComponent.extent_int(0); // shape (F,P[,D])
            const int numRightFieldsFinal = rightFinalComponent.extent_int(0); // shape (F,P[,D])
            
            Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numTensorComponents), [&] (const int& r) {
              const Data<Scalar,ExecSpaceType> &  leftTensorComponent =  leftComponent_.getTensorComponent(r);
              const Data<Scalar,ExecSpaceType> & rightTensorComponent = rightComponent_.getTensorComponent(r);
              const int  leftFieldCount =  leftTensorComponent.extent_int(0);
              const int      pointCount =  leftTensorComponent.extent_int(1);
              const int        leftRank =  leftTensorComponent.rank();
              const int rightFieldCount = rightTensorComponent.extent_int(0);
              const int       rightRank = rightTensorComponent.rank();
              for (int fieldOrdinal=0; fieldOrdinal<maxFieldsLeft_; fieldOrdinal++)
              {
                // slightly weird logic here in effort to avoid branch divergence
                const int fieldAddress = (fieldOrdinal < leftFieldCount) ? fieldOrdinal : leftFieldCount - 1;
                for (int pointOrdinal=0; pointOrdinal<maxPointCount_; pointOrdinal++)
                {
                  const int pointAddress = (pointOrdinal < pointCount) ? pointOrdinal : pointCount - 1;
                  leftFields(r,fieldAddress,pointAddress) = (leftRank == 2) ? leftTensorComponent(fieldAddress,pointAddress) : leftTensorComponent(fieldAddress,pointAddress,a_component);
                }
              }
              for (int fieldOrdinal=0; fieldOrdinal<maxFieldsRight_; fieldOrdinal++)
              {
                // slightly weird logic here in effort to avoid branch divergence
                const int fieldAddress = (fieldOrdinal < rightFieldCount) ? fieldOrdinal : rightFieldCount - 1;
                for (int pointOrdinal=0; pointOrdinal<maxPointCount_; pointOrdinal++)
                {
                  const int pointAddress = (pointOrdinal < pointCount) ? pointOrdinal : pointCount - 1;
                  rightFields(r,fieldAddress,pointAddress) = (rightRank == 2) ? rightTensorComponent(fieldAddress,pointAddress) : rightTensorComponent(fieldAddress,pointAddress,b_component);
                }
              }
            });
            
            Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,composedTransform_.extent_int(1)), [&] (const int& pointOrdinal) {
              pointWeights(pointOrdinal) = composedTransform_(cellDataOrdinal,pointOrdinal,a,b) * cellMeasures_(cellDataOrdinal,pointOrdinal);
            });
            // approximateFlopCount += composedTransform_.extent_int(1) * cellMeasures.numTensorComponents(); // cellMeasures does numTensorComponents - 1 multiplies on each access
            
            // synchronize threads
            teamMember.team_barrier();
            
            Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numLeftFieldsFinal * numRightFieldsFinal), [&] (const int& leftRightFieldEnumeration) {
              const int scratchOffsetForThread = teamMember.team_rank() * scratchViewSize;
              const int i_R = leftRightFieldEnumeration % numLeftFieldsFinal;
              const int j_R = leftRightFieldEnumeration / numLeftFieldsFinal;
              
              // as an optimization, we are using the same argument array for both the truncated field that spans s to R-1 and the full one that goes to R
              // the "R" argument here is ignored by the methods that treat the truncated field (it's beyond their bounds)
              Kokkos::Array<int,numTensorComponents>  leftFieldArguments;
              Kokkos::Array<int,numTensorComponents> rightFieldArguments;
              rightFieldArguments[R] = j_R;
              leftFieldArguments[R]  = i_R;
              
              //TODO: I believe that this can be moved outside the thread parallel_for
              for (int i=scratchOffsetForThread; i<scratchOffsetForThread+scratchViewSize; i++)
              {
                scratchView(i) = 0.0;
              }
              Kokkos::Array<int,numTensorComponents> pointArguments;
              for (unsigned i=0; i<numTensorComponents; i++)
              {
                pointArguments[i] = 0;
              }
              
              int r = R;
              while (r >= 0)
              {
                // integrate in final component dimension; this is where we need the M weight, as well as the weighted measure
                const int pointBounds_R = pointBounds[R];
                int & pointArgument_R   = pointArguments[R];
                for (pointArgument_R=0; pointArgument_R < pointBounds_R; pointArgument_R++)
                {
                  Scalar * G_R;
                  if (R != 0)
                  {
                    G_R = &scratchView(scratchOffsetForThread + offsetsForComponentOrdinal_[R]);
                  }
                  else
                  {
                    const int leftFieldIndex  = i_R + leftFieldOrdinalOffset_;
                    const int rightFieldIndex = j_R + rightFieldOrdinalOffset_;
                    
                    if (integralViewRank == 3)
                    {
                      // shape (C,F1,F2)
                      G_R = &integralView_.access(cellDataOrdinal,leftFieldIndex,rightFieldIndex);
                    }
                    else
                    {
                      // shape (F1,F2)
                      G_R = &integralView_.access(leftFieldIndex,rightFieldIndex,0);
                    }
                  }
                
                  const int & pointIndex_R = pointArguments[R];
                  
                  const Scalar &  leftValue =  leftFields(R,i_R,pointIndex_R);
                  const Scalar & rightValue = rightFields(R,j_R,pointIndex_R);
                  
                  int pointEnumerationIndex = relativeEnumerationIndex<numTensorComponents>(pointArguments, pointBounds, 0);
                  
                  *G_R += leftValue * pointWeights(pointEnumerationIndex) * rightValue;
                  
//                approximateFlopCount += 3; // 2 multiplies, 1 sum
                } // pointArgument_R

                const int r_next = nextIncrementResult<numTensorComponents,numTensorComponents-1>(pointArguments, pointBounds); // numTensorComponents_-1 means that we won't touch the [R] argument, which is treated in for loop above
                const int r_min  = (r_next >= 0) ? r_next : 0;
                
                for (int s=R-1; s>=r_min; s--)
                {
                  const int & pointIndex_s = pointArguments[s];
                  
                  // want to cover all the multi-indices from s to R-1
                  for (int s1=s; s1<R; s1++)
                  {
                    leftFieldArguments[s1] = 0;
                  }
                  
                  // i_s, j_s are the indices into the "current" tensor component; these are references, so they actually vary as the arguments are incremented.
                  const int & i_s = leftFieldArguments[s];
                  const int & j_s = rightFieldArguments[s];
                
                  int sLeft = s; // hereafter, sLeft is the return from the left field increment: the lowest rank that was incremented.  If this is less than s, we've gotten through all the multi-indexes from rank s through R-1 (inclusive).
                  const int & rightEnumerationSpan_s  = rightFieldRelativeEnumerationSpans_[s];
                  const int & rightEnumerationSpan_sp = rightFieldRelativeEnumerationSpans_[s+1];
                                    
                  while (sLeft >= s)
                  {
                    const int leftEnumerationIndex_s  = relativeEnumerationIndex<numTensorComponents,R>(leftFieldArguments, leftFieldBounds, s);
                    
                    // for final tensor component, the indices i_R, j_R are fixed, so we only have one slot for temporary storage (hence the "0" index possibility here)
                    const int leftEnumerationIndex_sp = (s+1 == R) ? 0 : relativeEnumerationIndex<numTensorComponents,R>(leftFieldArguments, leftFieldBounds, s+1);
                    
                    const Scalar &  leftValue = leftFields(s,i_s,pointIndex_s);
                    
                    for (int s1=s; s1<R; s1++)
                    {
                      rightFieldArguments[s1] = 0;
                    }
                    int sRight = s; // hereafter, sRight is the return from the leftFieldTensorIterator: the lowest rank that was incremented.  If this is less than s, we've gotten through all the multi-indexes from rank s through R-1 (inclusive).
                    while (sRight >= s)
                    {
                      const int rightEnumerationIndex_s  = relativeEnumerationIndex<numTensorComponents,R>(rightFieldArguments, rightFieldBounds, s);
                      
                      // for final tensor component, the indices i_R, j_R are fixed, so we only have one slot for temporary storage (hence the "0" index possibility here)
                      const int rightEnumerationIndex_sp = (s+1 == R) ? 0 : relativeEnumerationIndex<numTensorComponents,R>(rightFieldArguments, rightFieldBounds, s+1);
                                        
                      const Scalar & rightValue = rightFields(s,j_s,pointIndex_s);
                      
                      const int G_s_index  = leftEnumerationIndex_s  * rightEnumerationSpan_s + rightEnumerationIndex_s;
                      
                      Scalar* G_s;
                      
                      // for final tensor component, the indices i_R, j_R are fixed, so we only have one slot for temporary storage
                      // (above, we have set the leftEnumerationIndex_sp, rightEnumerationIndex_sp to be 0 in this case, so G_sp_index will then also be 0)
                      const int G_sp_index = leftEnumerationIndex_sp * rightEnumerationSpan_sp + rightEnumerationIndex_sp;

                      const Scalar & G_sp = scratchView(scratchOffsetForThread + offsetsForComponentOrdinal_[s+1] + G_sp_index);

                      
                      if (s != 0)
                      {
                        G_s = &scratchView(scratchOffsetForThread + offsetsForComponentOrdinal_[s] + G_s_index);
                      }
                      else
                      {
                        // compute enumeration indices
                        int  leftEnumerationIndex = relativeEnumerationIndex<numTensorComponents>( leftFieldArguments,  leftFieldBounds, 0);
                        int rightEnumerationIndex = relativeEnumerationIndex<numTensorComponents>(rightFieldArguments, rightFieldBounds, 0);
                        const int  leftFieldIndex =  leftEnumerationIndex +  leftFieldOrdinalOffset_;
                        const int rightFieldIndex = rightEnumerationIndex + rightFieldOrdinalOffset_;
                        
                        if (integralViewRank == 3)
                        {
                          // shape (C,F1,F2)
                          G_s = &integralView_.access(cellDataOrdinal,leftFieldIndex,rightFieldIndex);
                        }
                        else
                        {
                          // shape (F1,F2)
                          G_s = &integralView_.access(leftFieldIndex,rightFieldIndex,0);
                        }
                      }
                                            
                      *G_s += leftValue * G_sp * rightValue;
                      
//                    approximateFlopCount += 3; // 2 multiplies, 1 sum
                      
                      // increment rightField
                      sRight = incrementArgument<numTensorComponents,R>(rightFieldArguments, rightFieldBounds);
                    }
                    
                    // increment leftField
                    sLeft = incrementArgument<numTensorComponents,R>(leftFieldArguments, leftFieldBounds);
                  }
                }
                
                // clear tempStorage for r_next+1 to R
                if (r_min+1 <= R)
                {
                  const int endIndex = scratchOffsetForThread + offsetsForComponentOrdinal_[r_min];
                  for (int i=scratchOffsetForThread; i<endIndex; i++)
                  {
                    scratchView(i) = 0.0;
                  }
//                      auto tempStorageSubview = Kokkos::subview(scratchView, Kokkos::pair<int,int>{0,offsetsForComponentOrdinal_[r_min]});
//                      Kokkos::deep_copy(tempStorageSubview, 0.0);
                }
                
                // proceed to next point
                r = incrementArgument<numTensorComponents,numTensorComponents-1>(pointArguments, pointBounds);  // numTensorComponents_-1 means that we won't touch the [R] argument, which is treated in the G_R integration for loop above
              }
            }); // TeamThreadRange parallel_for
          }
        }
//        std::cout << "flop count per cell (within operator()) : " << approximateFlopCount << std::endl;
      }
      
      KOKKOS_INLINE_FUNCTION
      void operator()( const TeamMember & teamMember ) const
      {
        switch (numTensorComponents_)
        {
          case 1: run<1>(teamMember); break;
          case 2: run<2>(teamMember); break;
          case 3:
            if (forceNonSpecialized_)
              run<3>(teamMember);
            else
              runSpecialized3(teamMember);
            break;
          case 4: run<4>(teamMember); break;
          case 5: run<5>(teamMember); break;
          case 6: run<6>(teamMember); break;
          case 7: run<7>(teamMember); break;
#ifdef INTREPID2_HAVE_DEBUG
          default:
            INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true,std::invalid_argument,"Unsupported component count");
#endif
        }
      }
      
      //! returns an estimate of the number of floating point operations per cell (counting sums, subtractions, divisions, and multiplies, each of which counts as one operation).
      long approximateFlopCountPerCell() const
      {
        // compute flop count on a single cell, then multiply by the number of cells
        int flopCount = 0;
        
        const int R = numTensorComponents_ - 1;

        // we cache the value of M_ab * cellMeasure at each point.
        // access to cellMeasures involves cellMeasures.numTensorComponents() - 1 multiplies, so total is the below:
        const int flopsPerPoint_ab = cellMeasures_.numTensorComponents(); // the access involves multiplying all the components together
        
        for (int a_component=0; a_component < leftComponentSpan_; a_component++)
        {
          for (int b_component=0; b_component < rightComponentSpan_; b_component++)
          {
            const Data<Scalar,ExecSpaceType> &  leftFinalComponent =  leftComponent_.getTensorComponent(R);
            const Data<Scalar,ExecSpaceType> & rightFinalComponent = rightComponent_.getTensorComponent(R);
            
            const int numLeftFieldsFinal  =  leftFinalComponent.extent_int(0); // shape (F,P[,D])
            const int numRightFieldsFinal = rightFinalComponent.extent_int(0); // shape (F,P[,D])

            flopCount += flopsPerPoint_ab * cellMeasures_.extent_int(1);
            
            int flopsPer_i_R_j_R = 0;
            
            // as an optimization, we are using the same argument array for both the truncated field that spans s to R-1 and the full one that goes to R
            // the "R" argument here is ignored by the methods that treat the truncated field (it's beyond their bounds)
            Kokkos::Array<int,Parameters::MaxTensorComponents>  leftFieldArguments;
            leftFieldArguments[R] = 0;

            Kokkos::Array<int,Parameters::MaxTensorComponents> pointArguments;
            for (int i=0; i<numTensorComponents_; i++)
            {
              pointArguments[i] = 0;
            }
            
            // as an optimization, we are using the same argument array for both the truncated field that spans s to R-1 and the full one that goes to R
            // the "R" argument here is ignored by the methods that treat the truncated field (it's beyond their bounds)
            Kokkos::Array<int,Parameters::MaxTensorComponents> rightFieldArguments;
            rightFieldArguments[R] = 0;
            
            int r = R;
            while (r >= 0)
            {
              // integrate in final component dimension
              for (pointArguments[R]=0; pointArguments[R] < pointBounds_[R]; pointArguments[R]++)
              {
                flopsPer_i_R_j_R += 4;
              }
              // TODO: figure out why the below is not the same thing as the above -- the below overcounts, somehow
//                  if (0 < pointBounds_[R])
//                  {
//                    flopsPer_i_R_j_R += pointBounds_[R] * 4;
//                  }

              const int r_next = nextIncrementResult(pointArguments, pointBounds_, numTensorComponents_);
              const int r_min  = (r_next >= 0) ? r_next : 0;

              for (int s=R-1; s>=r_min; s--)
              {
                // want to cover all the multi-indices from s to R-1: for each we have 2 multiplies and one add (3 flops)
                int numLeftIterations  =  leftFieldRelativeEnumerationSpans_[s];
                int numRightIterations = rightFieldRelativeEnumerationSpans_[s];

                flopsPer_i_R_j_R += numLeftIterations * numRightIterations * 3;
              }

              // proceed to next point
              r = incrementArgument(pointArguments, pointBounds_, numTensorComponents_);
            }
            
            flopCount += flopsPer_i_R_j_R * numLeftFieldsFinal * numRightFieldsFinal;
          }
        }
//        std::cout << "flop count per cell: " << flopCount << std::endl;
        return flopCount;
      }
      
      //! returns the team size that should be provided to the policy constructor, based on the Kokkos maximum and the amount of thread parallelism we have available.
      int teamSize(const int &maxTeamSizeFromKokkos) const
      {
        const int R = numTensorComponents_ - 1;
        const int threadParallelismExpressed = leftFieldBounds_[R] * rightFieldBounds_[R];
        return std::min(maxTeamSizeFromKokkos, threadParallelismExpressed);
      }
      
      //! Provide the shared memory capacity.
      size_t team_shmem_size (int team_size) const
      {
        // we use shared memory to create a fast buffer for intermediate values, as well as fast access to the current-cell's field values
        size_t shmem_size = 0;
        
        if (fad_size_output_ > 0)
        {
          shmem_size += Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size(offsetsForComponentOrdinal_[0] * team_size, fad_size_output_);
          shmem_size += Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size(composedTransform_.extent_int(1),           fad_size_output_);
          shmem_size += Kokkos::View<Scalar***, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size(numTensorComponents_, maxFieldsLeft_,  maxPointCount_, fad_size_output_);
          shmem_size += Kokkos::View<Scalar***, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size(numTensorComponents_, maxFieldsRight_, maxPointCount_, fad_size_output_);
        }
        else
        {
          shmem_size += Kokkos::View<Scalar*, ExecSpaceType>::shmem_size(offsetsForComponentOrdinal_[0] * team_size);
          shmem_size += Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size(composedTransform_.extent_int(1));
          shmem_size += Kokkos::View<Scalar***, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size(numTensorComponents_, maxFieldsLeft_,  maxPointCount_);
          shmem_size += Kokkos::View<Scalar***, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size(numTensorComponents_, maxFieldsRight_, maxPointCount_);
        }
        
        return shmem_size;
      }
    };
  
    /**
        \brief Implementation of a general sum factorization algorithm, using a novel approach developed by Roberts, for integration.  Uses hierarchical parallelism.
     
     Whereas F_Integrate, and Mora and Demkowicz, and all others we are aware of, cache partial sums at intermediate component levels — the cached values are indexed
     by component basis ordinals — we integrate the first component in its dimension(s) and store values for integration points in the remaining dimensions, so that our caches
     are indexed by point ordinals.  If there are L_x, L_y, and L_z quadrature points in dimensions x,y,z, we require a cache of size L_y * L_z +1 for a 3D, 3-component integral.  The standard approach
     requires a cache of size (p_x+1)*(p_y+1).  So long as one is not over-integrating by too much, these sizes are about the same.  The real advantage of our approach here is (we expect)
     that it improves data locality.
       */
      template<class Scalar, class ExecSpaceType, int integralViewRank>
      class F_IntegratePointValueCache
      {
      using TeamPolicy = Kokkos::TeamPolicy<ExecSpaceType>;
      using TeamMember = typename TeamPolicy::member_type;

      using IntegralViewType = Kokkos::View<typename RankExpander<Scalar, integralViewRank>::value_type, ExecSpaceType>;
      IntegralViewType integralView_;
      TensorData<Scalar,ExecSpaceType> leftComponent_;
      Data<Scalar,ExecSpaceType> composedTransform_;
      TensorData<Scalar,ExecSpaceType> rightComponent_;
      TensorData<Scalar,ExecSpaceType> cellMeasures_;
      int a_offset_;
      int b_offset_;
      int leftComponentSpan_;   //  leftComponentSpan tracks the dimensions spanned by the left component
      int rightComponentSpan_;  // rightComponentSpan tracks the dimensions spanned by the right component
      int numTensorComponents_;
      int  leftFieldOrdinalOffset_;
      int rightFieldOrdinalOffset_;

      size_t fad_size_output_ = 0; // 0 if not a fad type

      // as an optimization, we do all the bounds and argument iteration within the functor rather than relying on TensorArgumentIterator
      // (this also makes it easier to reorder loops, etc., for further optimizations)
        Kokkos::Array<int,Parameters::MaxTensorComponents>  leftFieldBounds_;
        Kokkos::Array<int,Parameters::MaxTensorComponents> rightFieldBounds_;
      Kokkos::Array<int,Parameters::MaxTensorComponents> pointBounds_;

      int maxFieldsLeft_;
      int maxFieldsRight_;
      int maxPointCount_;
    public:
      F_IntegratePointValueCache(Data<Scalar,ExecSpaceType> integralData,
                                 TensorData<Scalar,ExecSpaceType> leftComponent,
                                 Data<Scalar,ExecSpaceType> composedTransform,
                                 TensorData<Scalar,ExecSpaceType> rightComponent,
                                 TensorData<Scalar,ExecSpaceType> cellMeasures,
                                 int a_offset,
                                 int b_offset,
                                 int leftFieldOrdinalOffset,
                                 int rightFieldOrdinalOffset)
      :
      integralView_(integralData.template getUnderlyingView<integralViewRank>()),
      leftComponent_(leftComponent),
      composedTransform_(composedTransform),
      rightComponent_(rightComponent),
      cellMeasures_(cellMeasures),
      a_offset_(a_offset),
      b_offset_(b_offset),
      leftComponentSpan_(leftComponent.extent_int(2)),
      rightComponentSpan_(rightComponent.extent_int(2)),
      numTensorComponents_(leftComponent.numTensorComponents()),
      leftFieldOrdinalOffset_(leftFieldOrdinalOffset),
      rightFieldOrdinalOffset_(rightFieldOrdinalOffset)
      {
        INTREPID2_TEST_FOR_EXCEPTION(numTensorComponents_ != rightComponent_.numTensorComponents(), std::invalid_argument, "Left and right components must have matching number of tensorial components");

        const int FIELD_DIM = 0;
        const int POINT_DIM = 1;
        maxFieldsLeft_  = 0;
        maxFieldsRight_ = 0;
        maxPointCount_  = 0;
        for (int r=0; r<numTensorComponents_; r++)
        {
          leftFieldBounds_[r]  = leftComponent_.getTensorComponent(r).extent_int(FIELD_DIM);
          maxFieldsLeft_       = std::max(maxFieldsLeft_, leftFieldBounds_[r]);
          rightFieldBounds_[r] = rightComponent_.getTensorComponent(r).extent_int(FIELD_DIM);
          maxFieldsRight_      = std::max(maxFieldsRight_, rightFieldBounds_[r]);
          pointBounds_[r]      = leftComponent_.getTensorComponent(r).extent_int(POINT_DIM);
          maxPointCount_       = std::max(maxPointCount_, pointBounds_[r]);
        }

        // prepare for allocation of temporary storage
        // note: tempStorage goes "backward", starting from the final component, which needs just one entry

        const bool allocateFadStorage = !std::is_pod<Scalar>::value;
        if (allocateFadStorage)
        {
          fad_size_output_ = dimension_scalar(integralView_);
        }
      }

      template<size_t maxComponents, size_t numComponents = maxComponents>
      KOKKOS_INLINE_FUNCTION
      int incrementArgument(      Kokkos::Array<int,maxComponents> &arguments,
                            const Kokkos::Array<int,maxComponents> &bounds) const
      {
        if (numComponents == 0) return -1;
        int r = static_cast<int>(numComponents - 1);
        while (arguments[r] + 1 >= bounds[r])
        {
          arguments[r] = 0; // reset
          r--;
          if (r < 0) break;
        }
        if (r >= 0) ++arguments[r];
        return r;
      }

      //! runtime-sized variant of incrementArgument; gets used by approximate flop count.
      KOKKOS_INLINE_FUNCTION
      int incrementArgument(      Kokkos::Array<int,Parameters::MaxTensorComponents> &arguments,
                            const Kokkos::Array<int,Parameters::MaxTensorComponents> &bounds,
                            const int &numComponents) const
      {
        if (numComponents == 0) return -1;
        int r = static_cast<int>(numComponents - 1);
        while (arguments[r] + 1 >= bounds[r])
        {
          arguments[r] = 0; // reset
          r--;
          if (r < 0) break;
        }
        if (r >= 0) ++arguments[r];
        return r;
      }

      template<size_t maxComponents, size_t numComponents = maxComponents>
      KOKKOS_INLINE_FUNCTION
      int nextIncrementResult(const Kokkos::Array<int,maxComponents> &arguments,
                              const Kokkos::Array<int,maxComponents> &bounds) const
      {
        if (numComponents == 0) return -1;
        int r = static_cast<int>(numComponents - 1);
        while (arguments[r] + 1 >= bounds[r])
        {
          r--;
          if (r < 0) break;
        }
        return r;
      }

      //! runtime-sized variant of nextIncrementResult; gets used by approximate flop count.
      KOKKOS_INLINE_FUNCTION
      int nextIncrementResult(const Kokkos::Array<int,Parameters::MaxTensorComponents> &arguments,
                              const Kokkos::Array<int,Parameters::MaxTensorComponents> &bounds,
                              const int &numComponents) const
      {
        if (numComponents == 0) return -1;
        int r = numComponents - 1;
        while (arguments[r] + 1 >= bounds[r])
        {
          r--;
          if (r < 0) break;
        }
        return r;
      }

      template<size_t maxComponents, size_t numComponents = maxComponents>
      KOKKOS_INLINE_FUNCTION
      int relativeEnumerationIndex(const Kokkos::Array<int,maxComponents> &arguments,
                                   const Kokkos::Array<int,maxComponents> &bounds,
                                   const int startIndex) const
      {
        // the following mirrors what is done in TensorData
        if (numComponents == 0)
        {
          return 0;
        }
        int enumerationIndex = 0;
        for (size_t r=numComponents-1; r>static_cast<size_t>(startIndex); r--)
        {
          enumerationIndex += arguments[r];
          enumerationIndex *= bounds[r-1];
        }
        enumerationIndex += arguments[startIndex];
        return enumerationIndex;
      }

      template<int rank>
      KOKKOS_INLINE_FUNCTION
      enable_if_t<rank==3 && rank==integralViewRank, Scalar &>
      integralViewEntry(const IntegralViewType& integralView, const int &cellDataOrdinal, const int &i, const int &j) const
      {
        return integralView(cellDataOrdinal,i,j);
      }
        
      template<int rank>
      KOKKOS_INLINE_FUNCTION
      enable_if_t<rank==2 && rank==integralViewRank, Scalar &>
      integralViewEntry(const IntegralViewType& integralView, const int &cellDataOrdinal, const int &i, const int &j) const
      {
        return integralView(i,j);
      }
        
        //! Hand-coded 3-component version
      KOKKOS_INLINE_FUNCTION
      void runSpecialized3( const TeamMember & teamMember ) const
      {
        constexpr int numTensorComponents = 3;
        
        const int pointBounds_x = pointBounds_[0];
        const int pointBounds_y = pointBounds_[1];
        const int pointBounds_z = pointBounds_[2];
        const int pointsInNonzeroComponentDimensions = pointBounds_y * pointBounds_z;
        
        const int  leftFieldBounds_x =  leftFieldBounds_[0];
        const int rightFieldBounds_x = rightFieldBounds_[0];
        const int  leftFieldBounds_y =  leftFieldBounds_[1];
        const int rightFieldBounds_y = rightFieldBounds_[1];
        const int  leftFieldBounds_z =  leftFieldBounds_[2];
        const int rightFieldBounds_z = rightFieldBounds_[2];
        
        Kokkos::Array<int,numTensorComponents>  leftFieldBounds;
        Kokkos::Array<int,numTensorComponents> rightFieldBounds;
        for (unsigned r=0; r<numTensorComponents; r++)
        {
          leftFieldBounds[r]  =  leftFieldBounds_[r];
          rightFieldBounds[r] = rightFieldBounds_[r];
        }
        
        const auto integralView = integralView_;
        const auto  leftFieldOrdinalOffset =  leftFieldOrdinalOffset_;
        const auto rightFieldOrdinalOffset = rightFieldOrdinalOffset_;

        const int cellDataOrdinal  = teamMember.league_rank();
        const int threadNumber     = teamMember.team_rank();
        
        const int numThreads       = teamMember.team_size(); // num threads
        const int GyEntryCount     = pointBounds_z; // for each thread: store one Gy value per z coordinate
        Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged> GxIntegrals; // for caching Gx values: we integrate out the first component dimension for each coordinate in the remaining dimensios
        Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged> GyIntegrals; // for caching Gy values (each thread gets a stack, of the same height as tensorComponents - 1)
        Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged> GzIntegral;  // for one Gz value that we sum into before summing into the destination matrix
        Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged> pointWeights; // indexed by (expanded) point; stores M_ab * cell measure; shared by team
        
        Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged> leftFields_x, rightFields_x;
        Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged> leftFields_y, rightFields_y;
        Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged> leftFields_z, rightFields_z;
        if (fad_size_output_ > 0) {
          GxIntegrals   = Kokkos::View<Scalar*,   ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(),  pointsInNonzeroComponentDimensions, fad_size_output_);
          GyIntegrals   = Kokkos::View<Scalar*,   ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(),  GyEntryCount * numThreads,          fad_size_output_);
          GzIntegral    = Kokkos::View<Scalar*,   ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(),  numThreads,                         fad_size_output_);
          pointWeights  = Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged>  (teamMember.team_shmem(),  composedTransform_.extent_int(1),   fad_size_output_);
          
          leftFields_x  = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(),  leftFieldBounds_x, pointBounds_x, fad_size_output_);
          rightFields_x = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), rightFieldBounds_x, pointBounds_x, fad_size_output_);
          leftFields_y  = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(),  leftFieldBounds_y, pointBounds_y, fad_size_output_);
          rightFields_y = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), rightFieldBounds_y, pointBounds_y, fad_size_output_);
          leftFields_z  = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(),  leftFieldBounds_z, pointBounds_z, fad_size_output_);
          rightFields_z = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), rightFieldBounds_z, pointBounds_z, fad_size_output_);
        }
        else {
          GxIntegrals   = Kokkos::View<Scalar*,   ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(),  pointsInNonzeroComponentDimensions);
          GyIntegrals   = Kokkos::View<Scalar*,   ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(),  GyEntryCount * numThreads);
          GzIntegral    = Kokkos::View<Scalar*,   ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(),  numThreads);
          pointWeights  = Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged>  (teamMember.team_shmem(),  composedTransform_.extent_int(1));
        
          leftFields_x  = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(),  leftFieldBounds_x, pointBounds_x);
          rightFields_x = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), rightFieldBounds_x, pointBounds_x);
          leftFields_y  = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(),  leftFieldBounds_y, pointBounds_y);
          rightFields_y = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), rightFieldBounds_y, pointBounds_y);
          leftFields_z  = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(),  leftFieldBounds_z, pointBounds_z);
          rightFields_z = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), rightFieldBounds_z, pointBounds_z);
        }

//        int approximateFlopCount = 0;
//        int flopsPerCellMeasuresAccess = cellMeasures_.numTensorComponents() - 1;

        // approximateFlopCount += composedTransform_.extent_int(1) * cellMeasures.numTensorComponents(); // cellMeasures does numTensorComponents - 1 multiplies on each access

        // synchronize threads
        teamMember.team_barrier();

        for (int a_component=0; a_component < leftComponentSpan_; a_component++)
        {
          const int a = a_offset_ + a_component;
          for (int b_component=0; b_component < rightComponentSpan_; b_component++)
          {
            const int b = b_offset_ + b_component;

            const Data<Scalar,ExecSpaceType> &  leftTensorComponent_x =  leftComponent_.getTensorComponent(0);
            const Data<Scalar,ExecSpaceType> & rightTensorComponent_x = rightComponent_.getTensorComponent(0);
            const Data<Scalar,ExecSpaceType> &  leftTensorComponent_y =  leftComponent_.getTensorComponent(1);
            const Data<Scalar,ExecSpaceType> & rightTensorComponent_y = rightComponent_.getTensorComponent(1);
            const Data<Scalar,ExecSpaceType> &  leftTensorComponent_z =  leftComponent_.getTensorComponent(2);
            const Data<Scalar,ExecSpaceType> & rightTensorComponent_z = rightComponent_.getTensorComponent(2);
            
            const int maxFields = (maxFieldsLeft_ > maxFieldsRight_) ? maxFieldsLeft_ : maxFieldsRight_;
            Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,maxFields), [&] (const int& fieldOrdinal) {
              if (fieldOrdinal < leftTensorComponent_x.extent_int(0))
              {
                const int pointCount = leftTensorComponent_x.extent_int(1);
                const int   leftRank = leftTensorComponent_x.rank();
                for (int pointOrdinal=0; pointOrdinal<pointCount; pointOrdinal++)
                {
                  leftFields_x(fieldOrdinal,pointOrdinal) = (leftRank == 2) ? leftTensorComponent_x(fieldOrdinal,pointOrdinal) : leftTensorComponent_x(fieldOrdinal,pointOrdinal,a_component);
                }
              }
              if (fieldOrdinal < leftTensorComponent_y.extent_int(0))
              {
                const int pointCount = leftTensorComponent_y.extent_int(1);
                const int   leftRank = leftTensorComponent_y.rank();
                for (int pointOrdinal=0; pointOrdinal<pointCount; pointOrdinal++)
                {
                  leftFields_y(fieldOrdinal,pointOrdinal) = (leftRank == 2) ? leftTensorComponent_y(fieldOrdinal,pointOrdinal) : leftTensorComponent_y(fieldOrdinal,pointOrdinal,a_component);
                }
              }
              if (fieldOrdinal < leftTensorComponent_z.extent_int(0))
              {
                const int pointCount = leftTensorComponent_z.extent_int(1);
                const int   leftRank = leftTensorComponent_z.rank();
                for (int pointOrdinal=0; pointOrdinal<pointCount; pointOrdinal++)
                {
                  leftFields_z(fieldOrdinal,pointOrdinal) = (leftRank == 2) ? leftTensorComponent_z(fieldOrdinal,pointOrdinal) : leftTensorComponent_z(fieldOrdinal,pointOrdinal,a_component);
                }
              }
              if (fieldOrdinal < rightTensorComponent_x.extent_int(0))
              {
                const int pointCount = rightTensorComponent_x.extent_int(1);
                const int   rightRank = rightTensorComponent_x.rank();
                for (int pointOrdinal=0; pointOrdinal<pointCount; pointOrdinal++)
                {
                  rightFields_x(fieldOrdinal,pointOrdinal) = (rightRank == 2) ? rightTensorComponent_x(fieldOrdinal,pointOrdinal) : rightTensorComponent_x(fieldOrdinal,pointOrdinal,a_component);
                }
              }
              if (fieldOrdinal < rightTensorComponent_y.extent_int(0))
              {
                const int pointCount = rightTensorComponent_y.extent_int(1);
                const int   rightRank = rightTensorComponent_y.rank();
                for (int pointOrdinal=0; pointOrdinal<pointCount; pointOrdinal++)
                {
                  rightFields_y(fieldOrdinal,pointOrdinal) = (rightRank == 2) ? rightTensorComponent_y(fieldOrdinal,pointOrdinal) : rightTensorComponent_y(fieldOrdinal,pointOrdinal,a_component);
                }
              }
              if (fieldOrdinal < rightTensorComponent_z.extent_int(0))
              {
                const int pointCount = rightTensorComponent_z.extent_int(1);
                const int   rightRank = rightTensorComponent_z.rank();
                for (int pointOrdinal=0; pointOrdinal<pointCount; pointOrdinal++)
                {
                  rightFields_z(fieldOrdinal,pointOrdinal) = (rightRank == 2) ? rightTensorComponent_z(fieldOrdinal,pointOrdinal) : rightTensorComponent_z(fieldOrdinal,pointOrdinal,a_component);
                }
              }
            });
                        
            Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,composedTransform_.extent_int(1)), [&] (const int& pointOrdinal) {
              pointWeights(pointOrdinal) = composedTransform_(cellDataOrdinal,pointOrdinal,a,b) * cellMeasures_(cellDataOrdinal,pointOrdinal);
            });

            // synchronize threads
            teamMember.team_barrier();
            
            for (int i0=0; i0<leftFieldBounds_x; i0++)
            {
              for (int j0=0; j0<rightFieldBounds_x; j0++)
              {
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,pointsInNonzeroComponentDimensions), [&] (const int& pointEnumerationIndexLaterDimensions) {
                  // first component is fastest-moving; we can get a tensorPointEnumerationOffset just by multiplying by the pointBounds in x
                  const int tensorPointEnumerationOffset = pointBounds_x * pointEnumerationIndexLaterDimensions; // compute offset for pointWeights container, for which x is the fastest-moving
                  
                  Scalar & Gx = GxIntegrals(pointEnumerationIndexLaterDimensions);
                  
                  Gx = 0;
                  if (fad_size_output_ == 0)
                  {
                    // not a Fad type; we're allow to have a vector range
                    Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(teamMember, pointBounds_x), [&] (const int &x_pointOrdinal, Scalar &integralThusFar)
                    {
                      integralThusFar += leftFields_x(i0,x_pointOrdinal) * rightFields_x(j0,x_pointOrdinal) * pointWeights(tensorPointEnumerationOffset + x_pointOrdinal);
                    }, Gx);
                  }
                  else
                  {
                    for (int x_pointOrdinal=0; x_pointOrdinal<pointBounds_x; x_pointOrdinal++)
                    {
                      Gx += leftFields_x(i0,x_pointOrdinal) * rightFields_x(j0,x_pointOrdinal) * pointWeights(tensorPointEnumerationOffset + x_pointOrdinal);
                    }
                  }
                });

                // synchronize threads
                teamMember.team_barrier();
                
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,leftFieldBounds_y * rightFieldBounds_y), [&] (const int& i1j1) {
                  const int i1 = i1j1 % leftFieldBounds_y;
                  const int j1 = i1j1 / leftFieldBounds_y;
                  
                  int Gy_index = GyEntryCount * threadNumber; // thread-relative index into GyIntegrals container; store one value per z coordinate
                  
                  int pointEnumerationIndex = 0; // incremented at bottom of lz loop below.
                  for (int lz=0; lz<pointBounds_z; lz++)
                  {
                    Scalar & Gy = GyIntegrals(Gy_index);
                    Gy = 0.0;
                    
                    for (int ly=0; ly<pointBounds_y; ly++)
                    {
                      const Scalar &  leftValue =  leftFields_y(i1,ly);
                      const Scalar & rightValue = rightFields_y(j1,ly);
                    
                      Gy += leftValue * rightValue * GxIntegrals(pointEnumerationIndex);
                      
                      pointEnumerationIndex++;
                    }
                    Gy_index++;
                  }
                      
                  Scalar & Gz = GzIntegral(threadNumber); // one entry per thread
                  for (int i2=0; i2<leftFieldBounds_z; i2++)
                  {
                    for (int j2=0; j2<rightFieldBounds_z; j2++)
                    {
                      Gz = 0.0;
                      
                      int Gy_index = GyEntryCount * threadNumber; // thread-relative index into GyIntegrals container; store one value per z coordinate
                      
                      for (int lz=0; lz<pointBounds_z; lz++)
                      {
                        const Scalar &  leftValue =  leftFields_z(i2,lz);
                        const Scalar & rightValue = rightFields_z(j2,lz);
                        
                        Gz += leftValue * rightValue * GyIntegrals(Gy_index);
                        
                        Gy_index++;
                      }
                      
                      const int i =  leftFieldOrdinalOffset + i0 + (i1 + i2 *  leftFieldBounds_y) *  leftFieldBounds_x;
                      const int j = rightFieldOrdinalOffset + j0 + (j1 + j2 * rightFieldBounds_y) * rightFieldBounds_x;
                      // the above are an optimization of the below, undertaken on the occasion of a weird Intel compiler segfault, possibly a compiler bug.
//                      const int i = relativeEnumerationIndex( leftArguments,  leftFieldBounds, 0) +  leftFieldOrdinalOffset;
//                      const int j = relativeEnumerationIndex(rightArguments, rightFieldBounds, 0) + rightFieldOrdinalOffset;
                      
                      integralViewEntry<integralViewRank>(integralView, cellDataOrdinal, i, j) += Gz;
                    }
                  }
                });
                // synchronize threads
                teamMember.team_barrier();
              }
            }
          }
        }
      }
        
      template<size_t numTensorComponents>
      KOKKOS_INLINE_FUNCTION
      void run( const TeamMember & teamMember ) const
      {
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::logic_error, "implementation incomplete");
//        Kokkos::Array<int,numTensorComponents>  pointBounds;
//        Kokkos::Array<int,numTensorComponents>  leftFieldBounds;
//        Kokkos::Array<int,numTensorComponents> rightFieldBounds;
//
//        int pointsInNonzeroComponentDimensions = 1;
//        for (unsigned r=0; r<numTensorComponents; r++)
//        {
//          pointBounds[r] = pointBounds_[r];
//          if (r > 0) pointsInNonzeroComponentDimensions *= pointBounds[r];
//          leftFieldBounds[r]  =  leftFieldBounds_[r];
//          rightFieldBounds[r] = rightFieldBounds_[r];
//        }
//
//        const int cellDataOrdinal  = teamMember.league_rank();
//        const int numThreads       = teamMember.team_size(); // num threads
//        const int G_k_StackHeight = numTensorComponents - 1; // per thread
//        Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged> G_0_IntegralsView; // for caching G0 values: we integrate out the first component dimension for each coordinate in the remaining dimensios
//        Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged> G_k_StackView; // for caching G_k values (each thread gets a stack, of the same height as tensorComponents - 1)
//        Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged> pointWeights; // indexed by (expanded) point; stores M_ab * cell measure; shared by team
//        Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged> leftFields, rightFields; // cache the field values at each level for faster access
//        if (fad_size_output_ > 0) {
//          G_k_StackView     = Kokkos::View<Scalar*,   ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), G_k_StackHeight * numThreads,       fad_size_output_);
//          G_0_IntegralsView = Kokkos::View<Scalar*,   ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), pointsInNonzeroComponentDimensions, fad_size_output_);
//          pointWeights   = Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged>  (teamMember.team_shmem(), composedTransform_.extent_int(1), fad_size_output_);
//          leftFields     = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), maxPointCount_, fad_size_output_);
//          rightFields    = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), maxPointCount_, fad_size_output_);
//        }
//        else {
//          G_k_StackView     = Kokkos::View<Scalar*,   ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), G_k_StackHeight * numThreads);
//          G_0_IntegralsView = Kokkos::View<Scalar*,   ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), pointsInNonzeroComponentDimensions);
//          pointWeights   = Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged>  (teamMember.team_shmem(), composedTransform_.extent_int(1));
//          leftFields     = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), maxPointCount_);
//          rightFields    = Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>(teamMember.team_shmem(), maxPointCount_);
//        }
//
////        int approximateFlopCount = 0;
////        int flopsPerCellMeasuresAccess = cellMeasures_.numTensorComponents() - 1;
//
//        constexpr int R = numTensorComponents - 1;
//
//        // approximateFlopCount += composedTransform_.extent_int(1) * cellMeasures.numTensorComponents(); // cellMeasures does numTensorComponents - 1 multiplies on each access
//
//        // synchronize threads
//        teamMember.team_barrier();
//
//        for (int a_component=0; a_component < leftComponentSpan_; a_component++)
//        {
//          const int a = a_offset_ + a_component;
//          for (int b_component=0; b_component < rightComponentSpan_; b_component++)
//          {
//            const int b = b_offset_ + b_component;
//
//            const Data<Scalar,ExecSpaceType> &  leftFirstComponent =  leftComponent_.getTensorComponent(0);
//            const Data<Scalar,ExecSpaceType> & rightFirstComponent = rightComponent_.getTensorComponent(0);
//
//            const int numLeftFieldsFirst  =  leftFirstComponent.extent_int(0); // shape (F,P[,D])
//            const int numRightFieldsFirst = rightFirstComponent.extent_int(0); // shape (F,P[,D])
//
//            const int numPointsFirst = leftFirstComponent.extent_int(1);
//
//            Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,composedTransform_.extent_int(1)), [&] (const int& pointOrdinal) {
//              pointWeights(pointOrdinal) = composedTransform_(cellDataOrdinal,pointOrdinal,a,b) * cellMeasures_(cellDataOrdinal,pointOrdinal);
//            });
//
//            // synchronize threads
//            teamMember.team_barrier();
//
//            for (int i0=0; i0<numLeftFieldsFirst; i0++)
//            {
//              for (int j0=0; j0<numRightFieldsFirst; j0++)
//              {
//                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numLeftFieldsFirst*numPointsFirst), [&] (const int& fieldOrdinalByPointOrdinal) {
//                  const int fieldOrdinal = fieldOrdinalByPointOrdinal % numPointsFirst;
//                  const int pointOrdinal = fieldOrdinalByPointOrdinal / numPointsFirst;
//                  leftFields(pointOrdinal) = leftFirstComponent(fieldOrdinal,pointOrdinal);
//                });
//
//                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numRightFieldsFirst*numPointsFirst), [&] (const int& fieldOrdinalByPointOrdinal) {
//                  const int fieldOrdinal = fieldOrdinalByPointOrdinal % numPointsFirst;
//                  const int pointOrdinal = fieldOrdinalByPointOrdinal / numPointsFirst;
//                  rightFields(pointOrdinal) = rightFirstComponent(fieldOrdinal,pointOrdinal);
//                });
//
//                // synchronize threads
//                teamMember.team_barrier();
//
//                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,pointsInNonzeroComponentDimensions), [&] (const int& pointEnumerationIndexLaterDimensions) {
//                  Kokkos::Array<int,numTensorComponents-1> pointArgumentsInLaterDimensions;
//                  int remainingIndex = pointEnumerationIndexLaterDimensions;
//
//                  for (int d=R-1; d>0; d--) // last component (z in 3D hypercube) is fastest-moving // TODO: consider doing first component as fastest-moving.  That would make indexing into pointWeights simpler
//                  {
//                    pointArgumentsInLaterDimensions[d] = pointEnumerationIndexLaterDimensions % pointBounds[d+1];
//                    remainingIndex /= pointBounds[d+1];
//                  }
//                  pointArgumentsInLaterDimensions[0] = remainingIndex;
//
//                  int tensorPointEnumerationOffset = 0; // compute offset for pointWeights container, for which x is the fastest-moving
//                  for (int d=R; d>0; d--)
//                  {
//                    tensorPointEnumerationOffset += pointArgumentsInLaterDimensions[d-1]; // pointArgumentsInLaterDimensions does not have an x component, hence d-1 here
//                    tensorPointEnumerationOffset *= pointBounds[d-1];
//                  }
//
//                  Scalar integralValue = 0;
//                  if (fad_size_output_ == 0)
//                  {
//                    // not a Fad type; we're allow to have a vector range
//                    Kokkos::parallel_reduce("first component integral", Kokkos::ThreadVectorRange(teamMember, numPointsFirst), [&] (const int &x_pointOrdinal, Scalar &integralThusFar)
//                    {
//                      integralThusFar += leftFields(x_pointOrdinal) * rightFields(x_pointOrdinal) * pointWeights(tensorPointEnumerationOffset);
//                    }, integralValue);
//                  }
//                  else
//                  {
//                    for (int pointOrdinal=0; pointOrdinal<numPointsFirst; pointOrdinal++)
//                    {
//                      integralValue += leftFields(pointOrdinal) * rightFields(pointOrdinal) * pointWeights(tensorPointEnumerationOffset);
//                    }
//                  }
//
//                  G_0_IntegralsView(pointEnumerationIndexLaterDimensions) = integralValue;
//                });
//
//                // synchronize threads
//                teamMember.team_barrier();
//
//                // TODO: finish this, probably after having written up the algorithm for arbitrary component count.  (I have it written down for 3D.)
//              }
//            }
//          }
//        }
////        std::cout << "flop count per cell (within operator()) : " << approximateFlopCount << std::endl;
      }

      KOKKOS_INLINE_FUNCTION
      void operator()( const TeamMember & teamMember ) const
      {
        switch (numTensorComponents_)
        {
          case 1: run<1>(teamMember); break;
          case 2: run<2>(teamMember); break;
          case 3: runSpecialized3(teamMember); break;
          case 4: run<4>(teamMember); break;
          case 5: run<5>(teamMember); break;
          case 6: run<6>(teamMember); break;
          case 7: run<7>(teamMember); break;
#ifdef INTREPID2_HAVE_DEBUG
          default:
            INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true,std::invalid_argument,"Unsupported component count");
#endif
        }
      }

      //! returns an estimate of the number of floating point operations per cell (counting sums, subtractions, divisions, and multiplies, each of which counts as one operation).
      long approximateFlopCountPerCell() const
      {
        // compute flop count on a single cell
        int flopCount = 0;

        constexpr int numTensorComponents = 3;
        Kokkos::Array<int,numTensorComponents>  pointBounds;
        Kokkos::Array<int,numTensorComponents>  leftFieldBounds;
        Kokkos::Array<int,numTensorComponents> rightFieldBounds;

        int pointsInNonzeroComponentDimensions = 1;
        for (unsigned r=0; r<numTensorComponents; r++)
        {
          pointBounds[r] = pointBounds_[r];
          if (r > 0) pointsInNonzeroComponentDimensions *= pointBounds[r];
          leftFieldBounds[r]  =  leftFieldBounds_[r];
          rightFieldBounds[r] = rightFieldBounds_[r];
        }

        for (int a_component=0; a_component < leftComponentSpan_; a_component++)
        {
          for (int b_component=0; b_component < rightComponentSpan_; b_component++)
          {
//            Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,composedTransform_.extent_int(1)), [&] (const int& pointOrdinal) {
//              pointWeights(pointOrdinal) = composedTransform_(cellDataOrdinal,pointOrdinal,a,b) * cellMeasures_(cellDataOrdinal,pointOrdinal);
//            });
            flopCount += composedTransform_.extent_int(1) * cellMeasures_.numTensorComponents(); // cellMeasures does numTensorComponents - 1 multiplies on each access

            for (int i0=0; i0<leftFieldBounds[0]; i0++)
            {
              for (int j0=0; j0<rightFieldBounds[0]; j0++)
              {
                flopCount += pointsInNonzeroComponentDimensions * pointBounds[0] * 3; // 3 flops per integration point in the loop commented out below
//                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,pointsInNonzeroComponentDimensions), [&] (const int& pointEnumerationIndexLaterDimensions) {
//                  Kokkos::Array<int,numTensorComponents-1> pointArgumentsInLaterDimensions;
//                  int remainingIndex = pointEnumerationIndexLaterDimensions;
//
//                  for (int d=0; d<R-1; d++) // first component is fastest-moving; this is determined by order of access in the lz/ly loop to compute Gy (integrals in y dimension)
//                  {
//                    pointArgumentsInLaterDimensions[d] = pointEnumerationIndexLaterDimensions % pointBounds[d+1]; // d+1 because x dimension is being integrated away
//                    remainingIndex /= pointBounds[d+1];
//                  }
//                  pointArgumentsInLaterDimensions[R-1] = remainingIndex;
//
//                  int tensorPointEnumerationOffset = 0; // compute offset for pointWeights container, for which x is the fastest-moving
//                  for (int d=R; d>0; d--)
//                  {
//                    tensorPointEnumerationOffset += pointArgumentsInLaterDimensions[d-1]; // pointArgumentsInLaterDimensions does not have an x component, hence d-1 here
//                    tensorPointEnumerationOffset *= pointBounds[d-1];
//                  }
//
//                  Scalar integralValue = 0;
//                  if (fad_size_output_ == 0)
//                  {
//                    // not a Fad type; we're allow to have a vector range
//                    Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(teamMember, numPointsFirst), [&] (const int &x_pointOrdinal, Scalar &integralThusFar)
//                    {
//                      integralThusFar += leftFields(x_pointOrdinal) * rightFields(x_pointOrdinal) * pointWeights(tensorPointEnumerationOffset + x_pointOrdinal);
//                    }, integralValue);
//                  }
//                  else
//                  {
//                    for (int x_pointOrdinal=0; x_pointOrdinal<numPointsFirst; x_pointOrdinal++)
//                    {
//                      integralValue += leftFields_x(x_pointOrdinal) * rightFields_x(x_pointOrdinal) * pointWeights(tensorPointEnumerationOffset + x_pointOrdinal);
//                    }
//                  }
//
//                  GxIntegrals(pointEnumerationIndexLaterDimensions) = integralValue;
//                });

                
                flopCount += leftFieldBounds[1] * rightFieldBounds[1] * pointBounds[1] * pointBounds[2] * 3; // 3 flops for each Gy += line in the below
                flopCount += leftFieldBounds[1] * rightFieldBounds[1] * leftFieldBounds[2] * rightFieldBounds[2] * pointBounds[2] * 3; // 3 flops for each Gz += line in the below
                flopCount += leftFieldBounds[1] * rightFieldBounds[1] * leftFieldBounds[2] * rightFieldBounds[2] * 1; // 1 flops for the integralView += line below
                
//                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,leftFieldBounds[1] * rightFieldBounds[1]), [&] (const int& i1j1) {
//                  const int i1 = i1j1 % leftFieldBounds[1];
//                  const int j1 = i1j1 / leftFieldBounds[1];
//
////                  int Gy_index = GyEntryCount * threadNumber; // thread-relative index into GyIntegrals container; store one value per z coordinate
//
//                  int pointEnumerationIndex = 0; // incremented at bottom of lz loop below.
//                  for (int lz=0; lz<pointBounds[2]; lz++)
//                  {
//                    Scalar & Gy = GyIntegrals(Gy_index);
//                    Gy = 0.0;
//
//                    const bool  leftRankIs3 = ( leftFields_y.rank() == 3);
//                    const bool rightRankIs3 = (rightFields_y.rank() == 3);
//                    for (int ly=0; ly<pointBounds[1]; ly++)
//                    {
//                      const Scalar &  leftValue =  leftRankIs3 ?  leftFields_y(i1,ly,a_component) :  leftFields_y(i1,ly);
//                      const Scalar & rightValue = rightRankIs3 ? rightFields_y(j1,ly,b_component) : rightFields_y(j1,ly);
//
//                      Gy += leftValue * rightValue * GxIntegrals(pointEnumerationIndex);
//
//                      pointEnumerationIndex++;
//                    }
//                    Gy_index++;
//                  }
//
//                  for (int i2=0; i2<leftFieldBounds[2]; i2++)
//                  {
//                    for (int j2=0; j2<rightFieldBounds[2]; j2++)
//                    {
//                      Scalar Gz = 0.0;
//
//                      int Gy_index = GyEntryCount * threadNumber; // thread-relative index into GyIntegrals container; store one value per z coordinate
//
//                      const bool  leftRankIs3 = ( leftFields_z.rank() == 3);
//                      const bool rightRankIs3 = (rightFields_z.rank() == 3);
//                      for (int lz=0; lz<pointBounds[2]; lz++)
//                      {
//                        const Scalar &  leftValue =  leftRankIs3 ?  leftFields_z(i2,lz,a_component) :  leftFields_z(i2,lz);
//                        const Scalar & rightValue = rightRankIs3 ? rightFields_z(j2,lz,b_component) : rightFields_z(j2,lz);
//
//                        Gz += leftValue * rightValue * GyIntegrals(Gy_index);
//
//                        Gy_index++;
//                      }
//
//                      Kokkos::Array<int,3>  leftArguments {i0,i1,i2};
//                      Kokkos::Array<int,3> rightArguments {j0,j1,j2};
//
//                      const int i = relativeEnumerationIndex( leftArguments,  leftFieldBounds, 0) +  leftFieldOrdinalOffset_;
//                      const int j = relativeEnumerationIndex(rightArguments, rightFieldBounds, 0) + rightFieldOrdinalOffset_;
//
//                      if (integralViewRank == 2)
//                      {
//                        integralView_.access(i,j,0) += Gz;
//                      }
//                      else
//                      {
//                        integralView_.access(cellDataOrdinal,i,j) += Gz;
//                      }
//                    }
//                  }
//                });
              }
            }
          }
        }
        return flopCount;
      }

      //! returns the team size that should be provided to the policy constructor, based on the Kokkos maximum and the amount of thread parallelism we have available.
      int teamSize(const int &maxTeamSizeFromKokkos) const
      {
        // TODO: fix this to match the actual parallelism expressed
        const int R = numTensorComponents_ - 1;
        const int threadParallelismExpressed = leftFieldBounds_[R] * rightFieldBounds_[R];
        return std::min(maxTeamSizeFromKokkos, threadParallelismExpressed);
      }

      //! Provide the shared memory capacity.
      size_t team_shmem_size (int numThreads) const
      {
        // we use shared memory to create a fast buffer for intermediate values, as well as fast access to the current-cell's field values
        size_t shmem_size = 0;
        
        const int GyEntryCount = pointBounds_[2]; // for each thread: store one Gy value per z coordinate

        int pointsInNonzeroComponentDimensions = 1;
        for (int d=1; d<numTensorComponents_; d++)
        {
          pointsInNonzeroComponentDimensions *= pointBounds_[d];
        }
        
        if (fad_size_output_ > 0)
        {
          shmem_size += Kokkos::View<Scalar*,   ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size(pointsInNonzeroComponentDimensions, fad_size_output_); // GxIntegrals: entries with x integrated away
          shmem_size += Kokkos::View<Scalar*,   ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size(GyEntryCount * numThreads,          fad_size_output_); // GyIntegrals: entries with x,y integrated away
          shmem_size += Kokkos::View<Scalar*,   ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size(           1 * numThreads,          fad_size_output_); // GzIntegral:  entry   with x,y,z integrated away
          shmem_size += Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size  (composedTransform_.extent_int(1), fad_size_output_); // pointWeights
          
          shmem_size += Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size(  leftFieldBounds_[0], pointBounds_[0], fad_size_output_); // leftFields_x
          shmem_size += Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size( rightFieldBounds_[0], pointBounds_[0], fad_size_output_); // rightFields_x
          shmem_size += Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size(  leftFieldBounds_[1], pointBounds_[1], fad_size_output_); // leftFields_y
          shmem_size += Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size( rightFieldBounds_[1], pointBounds_[1], fad_size_output_); // rightFields_y
          shmem_size += Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size(  leftFieldBounds_[2], pointBounds_[2], fad_size_output_); // leftFields_z
          shmem_size += Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size( rightFieldBounds_[2], pointBounds_[2], fad_size_output_); // rightFields_z
        }
        else
        {
          shmem_size += Kokkos::View<Scalar*,   ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size(pointsInNonzeroComponentDimensions);  // GxIntegrals: entries with x integrated away
          shmem_size += Kokkos::View<Scalar*,   ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size(GyEntryCount * numThreads);           // GyIntegrals: entries with x,y integrated away
          shmem_size += Kokkos::View<Scalar*,   ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size( 1 * numThreads);                     // GzIntegral:  entry   with x,y,z integrated away
          shmem_size += Kokkos::View<Scalar*, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size  (composedTransform_.extent_int(1)); // pointWeights
          
          shmem_size += Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size(  leftFieldBounds_[0], pointBounds_[0]); // leftFields_x
          shmem_size += Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size( rightFieldBounds_[0], pointBounds_[0]); // rightFields_x
          shmem_size += Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size(  leftFieldBounds_[1], pointBounds_[1]); // leftFields_y
          shmem_size += Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size( rightFieldBounds_[1], pointBounds_[1]); // rightFields_y
          shmem_size += Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size(  leftFieldBounds_[2], pointBounds_[2]); // leftFields_z
          shmem_size += Kokkos::View<Scalar**, ExecSpaceType, Kokkos::MemoryUnmanaged>::shmem_size( rightFieldBounds_[2], pointBounds_[2]); // rightFields_z
        }

        return shmem_size;
      }
    };
  }

template<typename ExecSpaceType>
template<class Scalar>
Data<Scalar,ExecSpaceType> IntegrationTools<ExecSpaceType>::allocateIntegralData(const TransformedVectorData<Scalar,ExecSpaceType> vectorDataLeft,
                                                                                 const TensorData<Scalar,ExecSpaceType> cellMeasures,
                                                                                 const TransformedVectorData<Scalar,ExecSpaceType> vectorDataRight)
{
  // allocates a (C,F,F) container for storing integral data
  
  // determine cellDataExtent and variation type.  We currently support CONSTANT, MODULAR, and GENERAL as possible output variation types, depending on the inputs.
  // If cellMeasures has non-trivial tensor structure, the rank-1 cell Data object is the first component.
  // If cellMeasures has trivial tensor structure, then the first and only component has the cell index in its first dimension.
  // I.e., either way the relevant Data object is cellMeasures.getTensorComponent(0)
  const int CELL_DIM = 0;
  const auto cellMeasureData = cellMeasures.getTensorComponent(0);
  const auto leftTransform = vectorDataLeft.transform();
  
  DimensionInfo combinedCellDimInfo = cellMeasureData.getDimensionInfo(CELL_DIM);
  combinedCellDimInfo = combinedDimensionInfo(combinedCellDimInfo,  vectorDataLeft.transform().getDimensionInfo(CELL_DIM));
  combinedCellDimInfo = combinedDimensionInfo(combinedCellDimInfo, vectorDataRight.transform().getDimensionInfo(CELL_DIM));

  DataVariationType cellVariationType = combinedCellDimInfo.variationType;
  int cellDataExtent                  = combinedCellDimInfo.dataExtent;
  
  const int numCells       = vectorDataLeft.numCells();
  const int numFieldsLeft  = vectorDataLeft.numFields();
  const int numFieldsRight = vectorDataRight.numFields();
  
  Kokkos::Array<int,3> extents {numCells, numFieldsLeft, numFieldsRight};
  Kokkos::Array<DataVariationType,3> variationTypes {cellVariationType,GENERAL,GENERAL};
  
  if (cellVariationType != CONSTANT)
  {
    Kokkos::View<Scalar***,ExecSpaceType> data("Intrepid2 integral data",cellDataExtent,numFieldsLeft,numFieldsRight);
    return Data<Scalar,ExecSpaceType>(data, extents, variationTypes);
  }
  else
  {
    Kokkos::View<Scalar**,ExecSpaceType> data("Intrepid2 integral data",numFieldsLeft,numFieldsRight);
    return Data<Scalar,ExecSpaceType>(data, extents, variationTypes);
  }
}

//! Two use cases:
//! 1. affine tensor-topology mesh: cellMeasures is a simple tensor product in this case.
//! 2. arbitrary mesh: cellMeasures has trivial tensor product structure (one tensorial component).
template<typename ExecSpaceType>
template<class Scalar>
void IntegrationTools<ExecSpaceType>::integrate(Data<Scalar,ExecSpaceType> integrals, const TransformedVectorData<Scalar,ExecSpaceType> & vectorDataLeft,
                                                const TensorData<Scalar,ExecSpaceType> & cellMeasures,
                                                const TransformedVectorData<Scalar,ExecSpaceType> & vectorDataRight, const bool sumInto,
                                                double* approximateFlops)
{
  if (approximateFlops != NULL)
  {
    *approximateFlops = 0;
  }
    
  // integral data may have shape (C,F1,F2) or (if the variation type is CONSTANT in the cell dimension) shape (F1,F2)
  const int integralViewRank = integrals.getUnderlyingViewRank();
  
  if (!sumInto)
  {
    integrals.clear();
  }
  
  const int cellDataExtent = integrals.getDataExtent(0);
  const int numFieldsLeft  = vectorDataLeft.numFields();
  const int numFieldsRight = vectorDataRight.numFields();
  const int spaceDim       = vectorDataLeft.spaceDim();
  
  INTREPID2_TEST_FOR_EXCEPTION(vectorDataLeft.spaceDim() != vectorDataRight.spaceDim(), std::invalid_argument, "vectorDataLeft and vectorDataRight must agree on the space dimension");
  
  // TODO: add more checks against validity of inputs
  
  // we require that the number of tensor components in the vectors are the same for each vector entry
  // this is not strictly necessary, but it makes implementation easier, and we don't at present anticipate other use cases
  int numTensorComponentsLeft = -1;
  const auto &refVectorLeft   = vectorDataLeft.vectorData();
  int numFamiliesLeft         = refVectorLeft.numFamilies();
  int numVectorComponentsLeft = refVectorLeft.numComponents();
  Kokkos::Array<int,7> maxFieldsForComponentLeft  {0,0,0,0,0,0,0};
  Kokkos::Array<int,7> maxFieldsForComponentRight {0,0,0,0,0,0,0};
  for (int familyOrdinal=0; familyOrdinal<numFamiliesLeft; familyOrdinal++)
  {
    for (int vectorComponent=0; vectorComponent<numVectorComponentsLeft; vectorComponent++)
    {
      const TensorData<Scalar,ExecSpaceType> &tensorData = refVectorLeft.getComponent(familyOrdinal,vectorComponent);
      if (tensorData.numTensorComponents() > 0)
      {
        if (numTensorComponentsLeft == -1)
        {
          numTensorComponentsLeft = tensorData.numTensorComponents();
        }
        INTREPID2_TEST_FOR_EXCEPTION(numVectorComponentsLeft != tensorData.numTensorComponents(), std::invalid_argument, "Each valid entry in vectorDataLeft must have the same number of tensor components as every other");
        for (int r=0; r<numTensorComponentsLeft; r++)
        {
          maxFieldsForComponentLeft[r] = std::max(tensorData.getTensorComponent(r).extent_int(0), maxFieldsForComponentLeft[r]);
        }
      }
    }
  }
  int numTensorComponentsRight = -1;
  const auto &refVectorRight   = vectorDataRight.vectorData();
  int numFamiliesRight         = refVectorRight.numFamilies();
  int numVectorComponentsRight = refVectorRight.numComponents();
  for (int familyOrdinal=0; familyOrdinal<numFamiliesRight; familyOrdinal++)
  {
    for (int vectorComponent=0; vectorComponent<numVectorComponentsRight; vectorComponent++)
    {
      const auto &tensorData = refVectorRight.getComponent(familyOrdinal,vectorComponent);
      if (tensorData.numTensorComponents() > 0)
      {
        if (numTensorComponentsRight == -1)
        {
          numTensorComponentsRight = tensorData.numTensorComponents();
        }
        INTREPID2_TEST_FOR_EXCEPTION(numVectorComponentsRight != tensorData.numTensorComponents(), std::invalid_argument, "Each valid entry in vectorDataRight must have the same number of tensor components as every other");
        for (int r=0; r<numTensorComponentsRight; r++)
        {
          maxFieldsForComponentRight[r] = std::max(tensorData.getTensorComponent(r).extent_int(0), maxFieldsForComponentRight[r]);
        }
      }
    }
  }
  INTREPID2_TEST_FOR_EXCEPTION(numVectorComponentsLeft != numVectorComponentsRight, std::invalid_argument, "Left and right vector entries must have the same number of tensorial components");
  const int numPointTensorComponents = cellMeasures.numTensorComponents() - 1;
    
  if ((numPointTensorComponents == numTensorComponentsLeft) && vectorDataLeft.axisAligned() && vectorDataRight.axisAligned())
  {
    // cellMeasures is a non-trivial tensor product, and the pullbacks are all diagonals.
    
    // in this case, the integrals in each tensorial direction are entirely separable
    // allocate some temporary storage for the integrals in each tensorial direction
    
    // if cellMeasures is a nontrivial tensor product, that means that all cells have the same shape, up to scaling.

    Kokkos::Array<int,Parameters::MaxTensorComponents> pointDimensions;
    for (int r=0; r<numPointTensorComponents; r++)
    {
      // first tensorial component of cell measures is the cell dimension; after that we have (P1,P2,…)
      pointDimensions[r] = cellMeasures.getTensorComponent(r+1).extent_int(0);
    }

    Kokkos::Array< Kokkos::Array<ScalarView<Scalar,ExecSpaceType>, Parameters::MaxTensorComponents>, Parameters::MaxTensorComponents> componentIntegrals; // these are reference-space integrals; no Jacobians or cell measures applied as yet
    
    for (int r=0; r<numPointTensorComponents; r++)
    {
      for (int d=0; d<spaceDim; d++)
      {
        /*
         Four cases to consider:
         1. Both vector data containers are filled with axial components - first component in 3D has form (f,0,0), second (0,f,0), third (0,0,f).
         2. Both vector data containers have arbitrary components - in 3D: (f1,f2,f3) where f1 is given by the first component, f2 by the second, f3 by the third.
         3. First container is axial, second arbitrary.
         4. First is arbitrary, second axial.
         
         But note that in all four cases, the structure of the integral is the same: you have three vector component integrals that get summed.  The actual difference between
         the cases does not show up in the reference-space integrals here, but in the accumulation in physical space below, where the tensor field numbering comes into play.
         
         The choice between axial and arbitrary affects the way the fields are numbered; the arbitrary components' indices refer to the same vector function, so they correspond,
         while the axial components refer to distinct scalar functions, so their numbering in the data container is cumulative.
        */
        
        const Data<Scalar,ExecSpaceType> &  leftComponent = vectorDataLeft.vectorData().getComponent(d).getTensorComponent(r);
        const Data<Scalar,ExecSpaceType> & rightComponent = vectorDataRight.vectorData().getComponent(d).getTensorComponent(r);
        
        // It may be worth considering the possibility that some of these components point to the same data -- if so, we could possibly get better data locality by making the corresponding componentIntegral entries point to the same location as well.  (And we could avoid some computations here.)
        
        const Data<Scalar,ExecSpaceType> & quadratureWeights = cellMeasures.getTensorComponent(r+1);
        
        int leftComponentCount  =  leftComponent.extent_int(0);
        int rightComponentCount = rightComponent.extent_int(0);
        
        const bool allocateFadStorage = !std::is_pod<Scalar>::value;
        if (allocateFadStorage)
        {
          auto fad_size_output = dimension_scalar(integrals.getUnderlyingView());
          componentIntegrals[r][d] = ScalarView<Scalar,ExecSpaceType>("componentIntegrals for tensor component " + std::to_string(r) + ", in dimension " + std::to_string(d), leftComponentCount, rightComponentCount, fad_size_output);
        }
        else
        {
          componentIntegrals[r][d] = ScalarView<Scalar,ExecSpaceType>("componentIntegrals for tensor component " + std::to_string(r) + ", in dimension " + std::to_string(d), leftComponentCount, rightComponentCount);
        }
        
        Kokkos::deep_copy(componentIntegrals[r][d], 0.0); // initialize
        
        const int & numPoints = pointDimensions[r];
        const int  leftComponentRank =  leftComponent.rank();
        const int rightComponentRank = rightComponent.rank();
        
        const int  leftComponentDimSpan =  leftComponent.extent_int(2);
        const int rightComponentDimSpan = rightComponent.extent_int(2);
        INTREPID2_TEST_FOR_EXCEPTION(( leftComponentDimSpan != rightComponentDimSpan), std::invalid_argument, "left and right components must span the same number of dimensions.");
        
//        // TODO: add support for cases where the component ranks are 3, and spaceDim > 1.
//        INTREPID2_TEST_FOR_EXCEPTION(( leftComponentRank == 3) && ( leftComponent.extent_int(2) > 1), std::invalid_argument, "spaceDim > 1 not yet supported by this integrate() use case.");
//        INTREPID2_TEST_FOR_EXCEPTION((rightComponentRank == 3) && (rightComponent.extent_int(2) > 1), std::invalid_argument, "spaceDim > 1 not yet supported by this integrate() use case.");
        
        auto policy = Kokkos::MDRangePolicy<ExecSpaceType,Kokkos::Rank<2>>({0,0},{leftComponentCount,rightComponentCount});
        
        for (int a=0; a<leftComponentDimSpan; a++)
        {
          Kokkos::parallel_for("compute componentIntegrals", policy,
          KOKKOS_LAMBDA (const int &i, const int &j) {
            Scalar refSpaceIntegral = 0.0;
            for (int ptOrdinal=0; ptOrdinal<numPoints; ptOrdinal++)
            {
              const Scalar &  leftValue = ( leftComponentRank == 2) ?  leftComponent(i,ptOrdinal) :  leftComponent(i,ptOrdinal,a);
              const Scalar & rightValue = (rightComponentRank == 2) ? rightComponent(j,ptOrdinal) : rightComponent(j,ptOrdinal,a);
              refSpaceIntegral += leftValue * rightValue * quadratureWeights(ptOrdinal);
            }
            componentIntegrals[r][d](i,j) = refSpaceIntegral;
          });
        }
        
        if (approximateFlops != NULL)
        {
          *approximateFlops += leftComponentCount*rightComponentCount*numPoints*(3); // two multiplies, one add in innermost loop
        }
      }
    }
    ExecSpaceType().fence(); // ensure that kernels launched above have completed before the kernels below launch.

    // only one of these will be a valid container:
    Kokkos::View<Scalar**,  ExecSpaceType> integralView2;
    Kokkos::View<Scalar***, ExecSpaceType> integralView3;
    if (integralViewRank == 3)
    {
      integralView3 = integrals.getUnderlyingView3();
    }
    else
    {
      integralView2 = integrals.getUnderlyingView2();
    }
    
    const int leftFamilyCount  =  vectorDataLeft.vectorData().numFamilies();
    const int rightFamilyCount = vectorDataRight.vectorData().numFamilies();
    
    for (int leftFamilyOrdinal=0; leftFamilyOrdinal<leftFamilyCount; leftFamilyOrdinal++)
    {
      bool haveLaunchedContributionToLeftFamily = false; // for tracking whether we need a fence() below
      
      int a_offset = 0; // left vector component offset
      int leftFieldOffset = vectorDataLeft.vectorData().familyFieldOrdinalOffset(leftFamilyOrdinal);
      
      const int leftVectorComponentCount = vectorDataLeft.vectorData().numComponents();
      for (int leftVectorComponentOrdinal = 0; leftVectorComponentOrdinal < leftVectorComponentCount; leftVectorComponentOrdinal++)
      {
        const auto & leftComponent = vectorDataLeft.vectorData().getComponent(leftFamilyOrdinal, leftVectorComponentOrdinal);
        if (!leftComponent.isValid())
        {
          a_offset++; // empty components are understood to take up one dimension
          continue;
        }
        const int leftDimSpan = leftComponent.extent_int(2);
          
        const int leftComponentFieldCount = leftComponent.extent_int(0);
        
        for (int rightFamilyOrdinal=0; rightFamilyOrdinal<rightFamilyCount; rightFamilyOrdinal++)
        {
          bool haveLaunchedContributionToRightFamily = false; // for tracking whether we need a fence() below
          
          int b_offset = 0; // right vector component offset
          int rightFieldOffset = vectorDataRight.vectorData().familyFieldOrdinalOffset(rightFamilyOrdinal);

          const int rightVectorComponentCount = vectorDataRight.vectorData().numComponents();
          for (int rightVectorComponentOrdinal = 0; rightVectorComponentOrdinal < rightVectorComponentCount; rightVectorComponentOrdinal++)
          {
            const auto & rightComponent = vectorDataRight.vectorData().getComponent(rightFamilyOrdinal, rightVectorComponentOrdinal);
            if (!rightComponent.isValid())
            {
              b_offset++; // empty components are understood to take up one dimension
              continue;
            }
            const int rightDimSpan = rightComponent.extent_int(2);
              
            const int rightComponentFieldCount = rightComponent.extent_int(0);
                    
            // we only accumulate for a == b (since this is the axis-aligned case).  Do the a, b spans intersect?
            if ((a_offset + leftDimSpan <= b_offset) || (b_offset + rightDimSpan <= a_offset)) // no intersection
            {
              b_offset += rightDimSpan;
              
              continue;
            }
            
            // if the a, b spans intersect, by assumption they should align exactly
            INTREPID2_TEST_FOR_EXCEPTION(( a_offset != b_offset), std::logic_error, "left and right dimension offsets should match.");
            INTREPID2_TEST_FOR_EXCEPTION(( leftDimSpan != rightDimSpan), std::invalid_argument, "left and right components must span the same number of dimensions.");
            
            const int d_start = a_offset;
            const int d_end   = d_start + leftDimSpan;
            
            if (haveLaunchedContributionToLeftFamily && haveLaunchedContributionToRightFamily)
            {
              // accumulation to the same block in integrals container; fence to ensure that that has completed before the next kernel is launched
              ExecSpaceType().fence();
            }
            haveLaunchedContributionToLeftFamily  = true;
            haveLaunchedContributionToRightFamily = true;
            
            std::initializer_list<int> upperBounds {cellDataExtent,leftComponentFieldCount,rightComponentFieldCount}; // separately declared in effort to get around Intel 17.0.1 compiler weirdness.
            auto policy = Kokkos::MDRangePolicy<ExecSpaceType,Kokkos::Rank<3>>({0,0,0}, upperBounds);
            // TODO: note that for best performance, especially with Fad types, we should replace this parallel for with a Functor and use hierarchical parallelism
            Kokkos::parallel_for("compute field integrals", policy,
                                 KOKKOS_LAMBDA (const int &cellDataOrdinal, const int &leftFieldOrdinal, const int &rightFieldOrdinal) {
              const Scalar & cellMeasureWeight = cellMeasures.getTensorComponent(0)(cellDataOrdinal);
          
              TensorArgumentIterator leftTensorIterator(leftComponent, 0); // shape is (F,P), and we walk the F dimension of the container
              leftTensorIterator.setEnumerationIndex(leftFieldOrdinal);
          
              TensorArgumentIterator rightTensorIterator(rightComponent, 0); // shape is (F,P), and we walk the F dimension of the container
              rightTensorIterator.setEnumerationIndex(rightFieldOrdinal);
          
              Scalar integralSum = 0.0;
              for (int d=d_start; d<d_end; d++)
              {
                const Scalar & transformLeft_d  =  vectorDataLeft.transformWeight(cellDataOrdinal,0,d,d);
                const Scalar & transformRight_d = vectorDataRight.transformWeight(cellDataOrdinal,0,d,d);
                
                const Scalar & leftRightTransform_d = transformLeft_d * transformRight_d;
        //            approximateFlopCount++;
                
                Scalar integral_d = 1.0;
                
                for (int r=0; r<numPointTensorComponents; r++)
                {
                  integral_d *= componentIntegrals[r][d](leftTensorIterator.argument(r),rightTensorIterator.argument(r));
        //              approximateFlopCount++; // product
                }
                integralSum += leftRightTransform_d * integral_d;
        //            approximateFlopCount += 2; // multiply and sum
                
                const int i =  leftFieldOrdinal +  leftFieldOffset;
                const int j = rightFieldOrdinal + rightFieldOffset;
                
                if (integralViewRank == 3)
                {
                  integralView3(cellDataOrdinal,i,j) += cellMeasureWeight * integralSum;
                }
                else
                {
                  integralView2(i,j) += cellMeasureWeight * integralSum;
                }
              }
            });
            
            b_offset += rightDimSpan;
          } // rightVectorComponentOrdinal loop
        } // rightFamilyOrdinal loop
        
        a_offset += leftDimSpan;
      } // leftVectorComponentOrdinal
    } // leftFamilyOrdinal
    
    if (approximateFlops != NULL)
    {
      // TODO: check the accuracy of this
      *approximateFlops += (2 + spaceDim * (3 + numPointTensorComponents)) * cellDataExtent * numFieldsLeft * numFieldsRight;
    }
  }
  else // general case (not axis-aligned + affine tensor-product structure)
  {
    // prepare composed transformation matrices
    const Data<Scalar,ExecSpaceType> & leftTransform  = vectorDataLeft.transform();
    const Data<Scalar,ExecSpaceType> & rightTransform = vectorDataRight.transform();
    const bool transposeLeft  = true;
    const bool transposeRight = false;
//    auto timer = Teuchos::TimeMonitor::getNewTimer("mat-mat");
//    timer->start();
    Data<Scalar,ExecSpaceType> composedTransform = leftTransform.allocateMatMatResult(transposeLeft, leftTransform, transposeRight, rightTransform);
    composedTransform.storeMatMat(transposeLeft, leftTransform, transposeRight, rightTransform);
//    timer->stop();
//    std::cout << "Completed mat-mat in " << timer->totalElapsedTime() << " seconds.\n";
//    timer->reset();
    
    // if the composedTransform matrices are full, the following is a good estimate.  If they have some diagonal portions, this will overcount.
    if (approximateFlops != NULL)
    {
      *approximateFlops += composedTransform.getUnderlyingViewSize() * (spaceDim - 1) * 2;
    }
    
    const int leftFamilyCount     = vectorDataLeft. vectorData().numFamilies();
    const int rightFamilyCount    = vectorDataRight.vectorData().numFamilies();
    const int leftComponentCount  = vectorDataLeft. vectorData().numComponents();
    const int rightComponentCount = vectorDataRight.vectorData().numComponents();
    
    int leftFieldOrdinalOffset = 0; // keeps track of the number of fields in prior families
    for (int leftFamilyOrdinal=0; leftFamilyOrdinal<leftFamilyCount; leftFamilyOrdinal++)
    {
      // "a" keeps track of the spatial dimension over which we are integrating in the left vector
      // components are allowed to span several dimensions; we keep track of the offset for the component in a_offset
      int a_offset = 0;
      bool haveLaunchedContributionToCurrentFamilyLeft = false; // helps to track whether we need a Kokkos::fence before launching a kernel.
      for (int leftComponentOrdinal=0; leftComponentOrdinal<leftComponentCount; leftComponentOrdinal++)
      {
        const TensorData<Scalar,ExecSpaceType> & leftComponent = vectorDataLeft.vectorData().getComponent(leftFamilyOrdinal, leftComponentOrdinal);
        if (!leftComponent.isValid())
        {
           // represents zero
          a_offset += vectorDataLeft.vectorData().numDimsForComponent(leftComponentOrdinal);
          continue;
        }
           
        int rightFieldOrdinalOffset = 0; // keeps track of the number of fields in prior families
        for (int rightFamilyOrdinal=0; rightFamilyOrdinal<rightFamilyCount; rightFamilyOrdinal++)
        {
          // "b" keeps track of the spatial dimension over which we are integrating in the right vector
          // components are allowed to span several dimensions; we keep track of the offset for the component in b_offset
          bool haveLaunchedContributionToCurrentFamilyRight = false; // helps to track whether we need a Kokkos::fence before launching a kernel.
          int b_offset = 0;
          for (int rightComponentOrdinal=0; rightComponentOrdinal<rightComponentCount; rightComponentOrdinal++)
          {
            const TensorData<Scalar,ExecSpaceType> & rightComponent = vectorDataRight.vectorData().getComponent(rightFamilyOrdinal, rightComponentOrdinal);
            if (!rightComponent.isValid())
            {
               // represents zero
              b_offset += vectorDataRight.vectorData().numDimsForComponent(rightComponentOrdinal);
              continue;
            }
            
            INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(leftComponent.numTensorComponents() != rightComponent.numTensorComponents(), std::invalid_argument, "left TensorData and right TensorData have different number of tensor components.  This is not supported.");
            
            const int vectorSize = getVectorSizeForHierarchicalParallelism<Scalar>();
            Kokkos::TeamPolicy<ExecSpaceType> policy = Kokkos::TeamPolicy<ExecSpaceType>(cellDataExtent,Kokkos::AUTO(),vectorSize);
            
            // TODO: expose the options for forceNonSpecialized and usePointCacheForRank3Tensor through an IntegrationAlgorithm enumeration.
            // AUTOMATIC: let Intrepid2 choose an algorithm based on the inputs (and, perhaps, the execution space)
            // STANDARD: don't use sum factorization or axis alignment -- just do the simple contraction, a (p+1)^9 algorithm in 3D
            // SUM_FACTORIZATION                              // (p+1)^7 algorithm in 3D
            // SUM_FACTORIZATION_AXIS_ALIGNED                 // (p+1)^6 algorithm in 3D
            // SUM_FACTORIZATION_FORCE_GENERIC_IMPLEMENTATION // mainly intended for testing purposes (specialized implementations perform better when they are provided)
            // SUM_FACTORIZATION_WITH_POINT_CACHE             // novel (p+1)^7 (in 3D) algorithm in Intrepid2; unclear as yet when and whether this may be a superior approach
            bool forceNonSpecialized = false; // We might expose this in the integrate() arguments in the future.  We *should* default to false in the future.
            bool usePointCacheForRank3Tensor = true; // EXPERIMENTAL; has better performance under CUDA, but slightly worse performance than standard on serial CPU
            
            // in one branch or another below, we will launch a parallel kernel that contributes to (leftFamily, rightFamily) field ordinal pairs.
            // if we have already launched something that contributes to that part of the integral container, we need a Kokkos fence() to ensure that these do not interfere with each other.
            if (haveLaunchedContributionToCurrentFamilyLeft && haveLaunchedContributionToCurrentFamilyRight)
            {
              ExecSpaceType().fence();
            }
            haveLaunchedContributionToCurrentFamilyLeft  = true;
            haveLaunchedContributionToCurrentFamilyRight = true;
            
            if (integralViewRank == 2)
            {
              if (usePointCacheForRank3Tensor && (leftComponent.numTensorComponents() == 3))
              {
                auto functor = Impl::F_IntegratePointValueCache<Scalar, ExecSpaceType, 2>(integrals, leftComponent, composedTransform, rightComponent, cellMeasures, a_offset, b_offset, leftFieldOrdinalOffset, rightFieldOrdinalOffset);
                
                const int maxTeamSize = policy.team_size_max(functor,Kokkos::ParallelForTag());
                const int teamSize    = functor.teamSize(maxTeamSize);
                
                policy = Kokkos::TeamPolicy<ExecSpaceType>(cellDataExtent,teamSize,vectorSize);
                
                Kokkos::parallel_for( policy , functor, "F_IntegratePointValueCache rank 2");
                
                if (approximateFlops != NULL)
                {
                  *approximateFlops += functor.approximateFlopCountPerCell() * integrals.getDataExtent(0);
                }
              }
              else
              {
                auto functor = Impl::F_Integrate<Scalar, ExecSpaceType, 2>(integrals, leftComponent, composedTransform, rightComponent, cellMeasures, a_offset, b_offset, leftFieldOrdinalOffset, rightFieldOrdinalOffset, forceNonSpecialized);
                
                const int maxTeamSize = policy.team_size_max(functor,Kokkos::ParallelForTag());
                const int teamSize    = functor.teamSize(maxTeamSize);
                
                policy = Kokkos::TeamPolicy<ExecSpaceType>(cellDataExtent,teamSize,vectorSize);
                
                Kokkos::parallel_for( policy , functor, "F_Integrate rank 2");
                
                if (approximateFlops != NULL)
                {
                  *approximateFlops += functor.approximateFlopCountPerCell() * integrals.getDataExtent(0);
                }
              }
            }
            else if (integralViewRank == 3)
            {
              if (usePointCacheForRank3Tensor && (leftComponent.numTensorComponents() == 3))
              {
                auto functor = Impl::F_IntegratePointValueCache<Scalar, ExecSpaceType, 3>(integrals, leftComponent, composedTransform, rightComponent, cellMeasures, a_offset, b_offset, leftFieldOrdinalOffset, rightFieldOrdinalOffset);
                
                const int maxTeamSize = policy.team_size_max(functor,Kokkos::ParallelForTag());
                const int teamSize    = functor.teamSize(maxTeamSize);
                
                policy = Kokkos::TeamPolicy<ExecSpaceType>(cellDataExtent,teamSize,vectorSize);
                
                Kokkos::parallel_for( policy , functor, "F_IntegratePointValueCache rank 3");
                
                if (approximateFlops != NULL)
                {
                  *approximateFlops += functor.approximateFlopCountPerCell() * integrals.getDataExtent(0);
                }
              }
              else
              {
                auto functor = Impl::F_Integrate<Scalar, ExecSpaceType, 3>(integrals, leftComponent, composedTransform, rightComponent, cellMeasures, a_offset, b_offset, leftFieldOrdinalOffset, rightFieldOrdinalOffset, forceNonSpecialized);
                
                const int maxTeamSize = policy.team_size_max(functor,Kokkos::ParallelForTag());
                const int teamSize    = functor.teamSize(maxTeamSize);
                
                policy = Kokkos::TeamPolicy<ExecSpaceType>(cellDataExtent,teamSize,vectorSize);
                
                Kokkos::parallel_for( policy , functor, "F_Integrate rank 3");
                
                if (approximateFlops != NULL)
                {
                  *approximateFlops += functor.approximateFlopCountPerCell() * integrals.getDataExtent(0);
                }
              }
            }
            b_offset += vectorDataRight.vectorData().numDimsForComponent(rightComponentOrdinal);
          }
          rightFieldOrdinalOffset += vectorDataRight.vectorData().numFieldsInFamily(rightFamilyOrdinal);
        }
        a_offset += vectorDataLeft.vectorData().numDimsForComponent(leftComponentOrdinal);
      }
      leftFieldOrdinalOffset += vectorDataLeft.vectorData().numFieldsInFamily(leftFamilyOrdinal);
    }
  }
//  if (approximateFlops != NULL)
//  {
//    std::cout << "Approximate flop count (new): " << *approximateFlops << std::endl;
//  }
  ExecSpaceType().fence(); // make sure we've finished writing to integrals container before exiting…
}

} // end namespace Intrepid2

#endif
