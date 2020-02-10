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

/** \file   Intrepid2_TensorBasis.hpp
 \brief  Implementation of bases that are tensor products of two or three component bases.
 \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_TensorBasis_h
#define Intrepid2_TensorBasis_h

#include <Kokkos_View.hpp>
#include <Kokkos_DynRankView.hpp>

#include <Intrepid2_config.h>

#include <map>
#include <set>
#include <vector>

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_DeviceAssert.hpp"
#include "Intrepid2_TensorTopologyMap.hpp"
#include "Intrepid2_TensorViewIterator.hpp"
#include "Intrepid2_Utils.hpp" // defines FAD_VECTOR_SIZE, VECTOR_SIZE

namespace Intrepid2
{
  template<unsigned spaceDim>
  KOKKOS_INLINE_FUNCTION
  void getDkEnumerationInverse(Kokkos::Array<int,spaceDim> &entries, const ordinal_type dkEnum, const ordinal_type operatorOrder);
  
  template<>
  KOKKOS_INLINE_FUNCTION
  void getDkEnumerationInverse<1>(Kokkos::Array<int,1> &entries, const ordinal_type dkEnum, const ordinal_type operatorOrder)
  {
    entries[0] = operatorOrder;
  }
  
  template<>
  KOKKOS_INLINE_FUNCTION
  void getDkEnumerationInverse<2>(Kokkos::Array<int,2> &entries, const ordinal_type dkEnum, const ordinal_type operatorOrder)
  {
    entries[0] = operatorOrder - dkEnum;
    entries[1] = dkEnum;
  }
  
  template<>
  KOKKOS_INLINE_FUNCTION
  void getDkEnumerationInverse<3>(Kokkos::Array<int,3> &entries, const ordinal_type dkEnum, const ordinal_type operatorOrder)
  {
    // formula is zMult + (yMult+zMult)*(yMult+zMult+1)/2; where xMult+yMult+zMult = operatorOrder
    // it seems complicated to work out a formula that will invert this.  For the present we just take a brute force approach,
    // using getDkEnumeration() to check each possibility
    for (ordinal_type yMult=0; yMult<=operatorOrder; yMult++)
    {
      for (ordinal_type zMult=0; zMult<=operatorOrder-yMult; zMult++)
      {
        const ordinal_type xMult = operatorOrder-(zMult+yMult);
        if (dkEnum == getDkEnumeration<3>(xMult,yMult,zMult))
        {
          entries[0] = xMult;
          entries[1] = yMult;
          entries[2] = zMult;
        }
      }
    }
  }
  
  template<unsigned spaceDim>
  ordinal_type getDkEnumeration(Kokkos::Array<int,spaceDim> &entries);
  
  template<>
  inline ordinal_type getDkEnumeration<1>(Kokkos::Array<int,1> &entries)
  {
    return getDkEnumeration<1>(entries[0]);
  }
  
  template<>
  inline ordinal_type getDkEnumeration<2>(Kokkos::Array<int,2> &entries)
  {
    return getDkEnumeration<2>(entries[0],entries[1]);
  }
  
  template<>
  inline ordinal_type getDkEnumeration<3>(Kokkos::Array<int,3> &entries)
  {
    return getDkEnumeration<3>(entries[0],entries[1],entries[2]);
  }
  
  template<unsigned spaceDim1, unsigned spaceDim2>
  inline ordinal_type getDkTensorIndex(const ordinal_type dkEnum1, const ordinal_type operatorOrder1,
                                       const ordinal_type dkEnum2, const ordinal_type operatorOrder2)
  {
    Kokkos::Array<int,spaceDim1> entries1;
    getDkEnumerationInverse<spaceDim1>(entries1, dkEnum1, operatorOrder1);
    
    Kokkos::Array<int,spaceDim2> entries2;
    getDkEnumerationInverse<spaceDim2>(entries2, dkEnum2, operatorOrder2);
    
    const int spaceDim = spaceDim1 + spaceDim2;
    Kokkos::Array<int,spaceDim> entries;
    
    for (unsigned d=0; d<spaceDim1; d++)
    {
      entries[d] = entries1[d];
    }
    
    for (unsigned d=0; d<spaceDim2; d++)
    {
      entries[d+spaceDim1] = entries2[d];
    }
    
    return getDkEnumeration<spaceDim>(entries);
  }

  /** \class  Intrepid2::TensorViewFunctor
      \brief  Functor for computing values for the TensorBasis class.
   
   This functor is not intended for use outside of \ref Intrepid2::Basis_TensorBasis.
  */
  template<class ExecutionSpace, class OutputScalar, class OutputFieldType>
  class TensorViewFunctor
  {
    using ScratchSpace       = Kokkos::DefaultExecutionSpace::scratch_memory_space;
    using OutputScratchView  = Kokkos::View<OutputScalar*,ScratchSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    
    using TeamPolicy = Kokkos::TeamPolicy<>;
    using TeamMember = TeamPolicy::member_type;
    
    using TensorViewIteratorType = ::Intrepid2::TensorViewIterator<OutputFieldType, OutputFieldType, OutputFieldType, OutputScalar>;
    using RankCombinationType = typename TensorViewIteratorType::RankCombinationType;
    
    OutputFieldType  output_; // F,P[,D…]
    OutputFieldType  input1_; // F1,P[,D…] or F1,P1[,D…]
    OutputFieldType  input2_; // F2,P[,D…] or F2,P2[,D…]
    
    int numFields_, numPoints_;
    int numFields1_, numPoints1_;
    int numFields2_, numPoints2_;
    
    bool tensorPoints_; // if true, input1 and input2 refer to values at decomposed points, and P = P1 * P2.  If false, then the two inputs refer to points in the full-dimensional space, and their point lengths are the same as that of the final output.
    
    Kokkos::vector<RankCombinationType> rank_combinations_; // indicates the policy by which the input views will be combined in output view
    
    double weight_;
    
  public:
    
    TensorViewFunctor(OutputFieldType output, OutputFieldType inputValues1, OutputFieldType inputValues2,
                      bool tensorPoints, double weight)
    : output_(output), input1_(inputValues1), input2_(inputValues2), tensorPoints_(tensorPoints), weight_(weight)
    {
      numFields_ = output.extent_int(0);
      numPoints_ = output.extent_int(1);
      
      numFields1_ = inputValues1.extent_int(0);
      numPoints1_ = inputValues1.extent_int(1);
      
      numFields2_ = inputValues2.extent_int(0);
      numPoints2_ = inputValues2.extent_int(1);
      
      if (!tensorPoints_)
      {
        // then the point counts should all match
        INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != numPoints1_, std::invalid_argument, "incompatible point counts");
        INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != numPoints2_, std::invalid_argument, "incompatible point counts");
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != numPoints1_ * numPoints2_, std::invalid_argument, "incompatible point counts");
      }
      
      INTREPID2_TEST_FOR_EXCEPTION(numFields_ != numFields1_ * numFields2_, std::invalid_argument, "incompatible field sizes");
      
      unsigned max_rank = std::max(inputValues1.rank(),inputValues2.rank());
      // at present, no supported case will result in an output rank greater than both input ranks
      INTREPID2_TEST_FOR_EXCEPTION(output.rank() > max_rank, std::invalid_argument, "Unsupported view combination.");
      rank_combinations_ = Kokkos::vector<RankCombinationType>(max_rank);
      
      rank_combinations_[0] = TensorViewIteratorType::TENSOR_PRODUCT; // field combination is always tensor product
      rank_combinations_[1] = tensorPoints ? TensorViewIteratorType::TENSOR_PRODUCT : TensorViewIteratorType::DIMENSION_MATCH; // tensorPoints controls interpretation of the point dimension
      for (unsigned d=2; d<max_rank; d++)
      {
        // d >= 2 have the interpretation of spatial dimensions (gradients, etc.)
        // we let the extents of the containers determine what we're doing here
        if ((inputValues1.extent_int(d) == inputValues2.extent_int(d)) && (output.extent_int(d) == 1))
        {
          rank_combinations_[d] = TensorViewIteratorType::TENSOR_CONTRACTION;
        }
        else if (((inputValues1.extent_int(d) == output.extent_int(d)) && (inputValues2.extent_int(d) == 1))
                 || ((inputValues2.extent_int(d) == output.extent_int(d)) && (inputValues1.extent_int(d) == 1))
                 )
        {
          // this looks like multiplication of a vector by a scalar, resulting in a vector
          // this can be understood as a tensor product
          rank_combinations_[d] = TensorViewIteratorType::TENSOR_PRODUCT;
        }
        else if ((inputValues1.extent_int(d) == inputValues2.extent_int(d)) && (output.extent_int(d) == inputValues1.extent_int(d) * inputValues2.extent_int(d)))
        {
          // this is actually a generalization of the above case: a tensor product, something like a vector outer product
          rank_combinations_[d] = TensorViewIteratorType::TENSOR_PRODUCT;
        }
        else if ((inputValues1.extent_int(d) == inputValues2.extent_int(d)) && (output.extent_int(d) == inputValues1.extent_int(d)))
        {
          // it's a bit weird (I'm not aware of the use case, in the present context), but we can handle this case by adopting DIMENSION_MATCH here
          // this is something like MATLAB's .* and .+ operators, which operate entry-wise
          rank_combinations_[d] = TensorViewIteratorType::DIMENSION_MATCH;
        }
        else
        {
          std::cout << "inputValues1.extent_int(" << d << ") = " << inputValues1.extent_int(d) << std::endl;
          std::cout << "inputValues2.extent_int(" << d << ") = " << inputValues2.extent_int(d) << std::endl;
          std::cout << "output.extent_int("       << d << ") = " << output.extent_int(d) << std::endl;
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "unable to find an interpretation for this combination of views");
        }
      }
    }
    
    KOKKOS_INLINE_FUNCTION
    void operator()( const TeamMember & teamMember ) const
    {
      auto fieldOrdinal1 = teamMember.league_rank();
      
      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numFields2_), [&] (const int& fieldOrdinal2) {
        TensorViewIteratorType it(output_,input1_,input2_,rank_combinations_);
        const int FIELD_ORDINAL_DIMENSION = 0;
        it.setLocation({fieldOrdinal1,0,0,0,0,0,0},{fieldOrdinal2,0,0,0,0,0,0});
        int next_increment_rank = FIELD_ORDINAL_DIMENSION; // used to initialize prev_increment_rank at the start of the do/while loop.  Notionally, we last incremented in the field ordinal rank to get to the {fieldOrdinal1,0,0,0,0,0,0},{fieldOrdinal2,0,0,0,0,0,0} location.
        OutputScalar accumulator = 0;
        
        do
        {
          accumulator += weight_ * it.getView1Entry() * it.getView2Entry();
          next_increment_rank = it.nextIncrementRank();
          
          if ((next_increment_rank < 0) || (rank_combinations_[next_increment_rank] != TensorViewIteratorType::TENSOR_CONTRACTION))
          {
            // then we've finished the accumulation and should set the value
            it.set(accumulator);
            // reset the accumulator:
            accumulator = 0;
          }
        } while (it.increment() > FIELD_ORDINAL_DIMENSION);
      });
    }
  };
  
  /** \class  Intrepid2::Basis_TensorBasis
      \brief  Basis defined as the tensor product of two component bases.
   
   The cell topology for the tensor basis is the tensor product of the cell topologies on which the component bases are defined;
   \see Intrepid2::TensorTopologyMap.
   
   The basis is ordered such that the Basis1 field ordinals are the fastest-moving index; the formula for the composite field ordinal is:
     compositeFieldOrdinal = componentFieldOrdinal2 * basis1Cardinality + componentFieldOrdinal1
   This is done so that we can consider, e.g. Basis1 as the "x" dimension and Basis2 as the "y" dimension, and have the basis ordered in
   the same way that existing Intrepid2 bases on the quadrilateral are ordered, namely, one moves along the x dimension first, moving
   across the quadrilateral dofs "row-wise".
   
  */
  template<typename Basis1, typename Basis2>
  class Basis_TensorBasis
  :
  public Basis<typename Basis1::ExecutionSpace,typename Basis1::OutputValueType,typename Basis1::PointValueType>
  {
  protected:
    Basis1 basis1_;
    Basis2 basis2_;
  public:
    using BasisSuper = ::Intrepid2::Basis<typename Basis1::ExecutionSpace,typename Basis1::OutputValueType,typename Basis1::PointValueType>;
    
    using ExecutionSpace  = typename BasisSuper::ExecutionSpace;
    using OutputValueType = typename BasisSuper::OutputValueType;
    using PointValueType  = typename BasisSuper::PointValueType;
    
    using OrdinalTypeArray1DHost = typename BasisSuper::OrdinalTypeArray1DHost;
    using OrdinalTypeArray2DHost = typename BasisSuper::OrdinalTypeArray2DHost;
    using OutputViewType         = typename BasisSuper::OutputViewType;
    using PointViewType          = typename BasisSuper::PointViewType;
  public:
    /** \brief  Constructor.
        \param [in] basis1 - the first component basis
        \param [in] basis2 - the second component basis
     */
    Basis_TensorBasis(Basis1 basis1, Basis2 basis2)
    :
    basis1_(basis1),basis2_(basis2)
    {
      this->basisCardinality_  = basis1.getCardinality() * basis2.getCardinality();
      this->basisDegree_       = std::max(basis1.getDegree(), basis2.getDegree());
      
      // set cell topology
      shards::CellTopology cellTopo1 = basis1.getBaseCellTopology();
      shards::CellTopology cellTopo2 = basis2.getBaseCellTopology();
      
      auto cellKey1 = basis1.getBaseCellTopology().getKey();
      auto cellKey2 = basis2.getBaseCellTopology().getKey();
      if ((cellKey1 == shards::Line<2>::key) && (cellKey2 == shards::Line<2>::key))
      {
        this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
      }
      else if (((cellKey1 == shards::Quadrilateral<4>::key) && (cellKey2 == shards::Line<2>::key))
               ||     ((cellKey2 == shards::Quadrilateral<4>::key) && (cellKey1 == shards::Line<2>::key)))
      {
        this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() );
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Cell topology combination not yet supported");
      }
      
      this->basisType_         = basis1.getBasisType();
      this->basisCoordinates_  = COORDINATES_CARTESIAN;
      
      // initialize tags
      {
        const auto & cardinality = this->basisCardinality_;
        
        // Basis-dependent initializations
        const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
        const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
        const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
        const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
        
        OrdinalTypeArray1DHost tagView("tag view", cardinality*tagSize);
        
        shards::CellTopology cellTopo = this->basisCellTopology_;
        
        int tensorSpaceDim  = cellTopo.getDimension();
        int spaceDim1       = cellTopo1.getDimension();
        int spaceDim2       = cellTopo2.getDimension();
        
        if (this->getBasisType() == BASIS_FEM_HIERARCHICAL)
        {
          int degreeSize = basis1_.getPolynomialDegreeLength() + basis2_.getPolynomialDegreeLength();
          this->fieldOrdinalPolynomialDegree_ = OrdinalTypeArray2DHost("TensorBasis - field ordinal polynomial degree", this->basisCardinality_, degreeSize);
        }
        
        TensorTopologyMap topoMap(cellTopo1, cellTopo2);
        
        for (int d=0; d<=tensorSpaceDim; d++) // d: tensorial dimension
        {
          int d2_max = std::min(spaceDim2,d);
          int subcellOffset = 0; // for this dimension of tensor subcells, how many subcells have we already counted with other d2/d1 combos?
          for (int d2=0; d2<=d2_max; d2++)
          {
            int d1 = d-d2;
            if (d1 > spaceDim1) continue;
            
            unsigned subcellCount2 = cellTopo2.getSubcellCount(d2);
            unsigned subcellCount1 = cellTopo1.getSubcellCount(d1);
            for (unsigned subcellOrdinal2=0; subcellOrdinal2<subcellCount2; subcellOrdinal2++)
            {
              ordinal_type subcellDofCount2 = basis2_.getDofCount(d2, subcellOrdinal2);
              for (unsigned subcellOrdinal1=0; subcellOrdinal1<subcellCount1; subcellOrdinal1++)
              {
                ordinal_type subcellDofCount1 = basis1_.getDofCount(d1, subcellOrdinal1);
                ordinal_type tensorLocalDofCount = subcellDofCount1 * subcellDofCount2;
                for (ordinal_type localDofID2 = 0; localDofID2<subcellDofCount2; localDofID2++)
                {
                  ordinal_type fieldOrdinal2 = basis2_.getDofOrdinal(d2, subcellOrdinal2, localDofID2);
                  OrdinalTypeArray1DHost degreesField2;
                  if (this->basisType_ == BASIS_FEM_HIERARCHICAL) degreesField2 = basis2_.getPolynomialDegreeOfField(fieldOrdinal2);
                  for (ordinal_type localDofID1 = 0; localDofID1<subcellDofCount1; localDofID1++)
                  {
                    ordinal_type fieldOrdinal1 = basis1_.getDofOrdinal(d1, subcellOrdinal1, localDofID1);
                    ordinal_type tensorLocalDofID = localDofID2 * subcellDofCount1 + localDofID1;
                    ordinal_type tensorFieldOrdinal = fieldOrdinal2 * basis1_.getCardinality() + fieldOrdinal1;
                    tagView(tensorFieldOrdinal*tagSize+0) = d; // subcell dimension
                    tagView(tensorFieldOrdinal*tagSize+1) = topoMap.getCompositeSubcellOrdinal(d1, subcellOrdinal1, d2, subcellOrdinal2);
                    tagView(tensorFieldOrdinal*tagSize+2) = tensorLocalDofID;
                    tagView(tensorFieldOrdinal*tagSize+3) = tensorLocalDofCount;
                    
                    if (this->basisType_ == BASIS_FEM_HIERARCHICAL)
                    {
                      // fill in degree lookup:
                      OrdinalTypeArray1DHost degreesField1 = basis1_.getPolynomialDegreeOfField(fieldOrdinal1);
                      
                      int degreeLengthField1 = degreesField1.extent_int(0);
                      int degreeLengthField2 = degreesField2.extent_int(0);
                      for (int d3=0; d3<degreeLengthField1; d3++)
                      {
                        this->fieldOrdinalPolynomialDegree_(tensorFieldOrdinal,d3) = degreesField1(d3);
                      }
                      for (int d3=0; d3<degreeLengthField2; d3++)
                      {
                        this->fieldOrdinalPolynomialDegree_(tensorFieldOrdinal,d3+degreeLengthField1) = degreesField2(d3);
                      }
                    }
                  } // localDofID1
                } // localDofID2
              } // subcellOrdinal1
            } // subcellOrdinal2
            subcellOffset += subcellCount1 * subcellCount2;
          }
        }
        
        //        // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
        //        // tags are constructed on host
        this->setOrdinalTagData(this->tagToOrdinal_,
                                this->ordinalToTag_,
                                tagView,
                                this->basisCardinality_,
                                tagSize,
                                posScDim,
                                posScOrd,
                                posDfOrd);
      }
    }
    
    // since the getValues() below only overrides the FEM variant, we specify that
    // we use the base class's getValues(), which implements the FVD variant by throwing an exception.
    // (It's an error to use the FVD variant on this basis.)
    using BasisSuper::getValues;
    
    /** \brief  Method to extract component points from composite points.
        \param [in]  inputPoints                  - points defined on the composite cell topology
        \param [in]  attemptTensorDecomposition   - if true, attempt to find a tensor decomposition.
        \param [out] inputPoints1                 - points defined on the first component cell topology
        \param [out] inputPoints2                 - points defined on the second component cell topology
        \param [out] tensorDecompositionSucceeded - if true, the attempt to find a tensor decomposition succeeded.
     
     At present, attemptTensorDecomposition is ignored, and tensorDecompositionSucceeded will always return false.
     However, we intend to support the tensor decomposition in the future, which will allow substantial optimizations
     in computation of tensor bases.
     */
    void getComponentPoints(const PointViewType inputPoints, const bool attemptTensorDecomposition,
                            PointViewType & inputPoints1, PointViewType & inputPoints2, bool &tensorDecompositionSucceeded) const
    {
      INTREPID2_TEST_FOR_EXCEPTION(attemptTensorDecomposition, std::invalid_argument, "tensor decomposition not yet supported");
      
      // for inputPoints that are actually tensor-product of component quadrature points (say),
      // having just the one input (which will have a lot of redundant point data) is suboptimal
      // The general case can have unique x/y/z coordinates at every point, though, so we have to support that
      // when this interface is used.  But we may try detecting that the data is tensor-product and compressing
      // from there...  Ultimately, we should also add a getValues() variant that takes multiple input point containers,
      // one for each tensorial dimension.
      
      // this initial implementation is intended to simplify development of 2D and 3D bases, while also opening
      // the possibility of higher-dimensional bases.  It is not necessarily optimized for speed/memory.  There
      // are things we can do in this regard, which may become important for matrix-free computations wherein
      // basis values don't get stored but are computed dynamically.
      
      int spaceDim1 = basis1_.getBaseCellTopology().getDimension();
      int spaceDim2 = basis2_.getBaseCellTopology().getDimension();
      
      int totalSpaceDim   = inputPoints.extent_int(1);
      
      TEUCHOS_ASSERT(spaceDim1 + spaceDim2 == totalSpaceDim);
      
      // first pass: just take subviews to get input points -- this will result in redundant computations when points are themselves tensor product (i.e., inputPoints itself contains redundant data)
      
      inputPoints1 = Kokkos::subview(inputPoints,Kokkos::ALL(),std::make_pair(0,spaceDim1));
      inputPoints2 = Kokkos::subview(inputPoints,Kokkos::ALL(),std::make_pair(spaceDim1,totalSpaceDim));
      
      //      std::cout << "inputPoints : " << inputPoints.extent(0) << " x " << inputPoints.extent(1) << std::endl;
      //      std::cout << "inputPoints1 : " << inputPoints1.extent(0) << " x " << inputPoints1.extent(1) << std::endl;
      //      std::cout << "inputPoints2 : " << inputPoints2.extent(0) << " x " << inputPoints2.extent(1) << std::endl;
      
      tensorDecompositionSucceeded = false;
    }
    
    /** \brief  Fills in spatial locations (coordinates) of degrees of freedom (nodes) on the reference cell
        \param [out] dofCoords - the container into which to place the degrees of freedom.
     
     dofCoords should have shape (F,D), where the field dimension matches the cardinality of the basis, and D is the
     spatial dimension of the topology on which the basis is defined.
     
     Note that getDofCoords() is not supported by all bases; in particular, hierarchical bases do not generally support this.
     */
    virtual void getDofCoords( typename BasisSuper::ScalarViewType dofCoords ) const override
    {
      int spaceDim1 = basis1_.getBaseCellTopology().getDimension();
      int spaceDim2 = basis2_.getBaseCellTopology().getDimension();
      
      using ValueType    = typename BasisSuper::ScalarViewType::value_type;
      using ResultLayout = typename DeduceLayout< typename BasisSuper::ScalarViewType >::result_layout;
      using DeviceType   = typename BasisSuper::ScalarViewType::device_type;
      using ViewType     = Kokkos::DynRankView<ValueType, ResultLayout, DeviceType >;
      
      ViewType dofCoords1("dofCoords1",basis1_.getCardinality(),spaceDim1);
      ViewType dofCoords2("dofCoords2",basis2_.getCardinality(),spaceDim2);
      
      basis1_.getDofCoords(dofCoords1);
      basis2_.getDofCoords(dofCoords2);
      
      const ordinal_type basisCardinality1 = basis1_.getCardinality();
      const ordinal_type basisCardinality2 = basis2_.getCardinality();
      
      Kokkos::parallel_for(basisCardinality2, KOKKOS_LAMBDA (const int fieldOrdinal2)
                           {
                             for (int fieldOrdinal1=0; fieldOrdinal1<basisCardinality1; fieldOrdinal1++)
                             {
                               const ordinal_type fieldOrdinal = fieldOrdinal1 + fieldOrdinal2 * basisCardinality1;
                               for (int d1=0; d1<spaceDim1; d1++)
                               {
                                 dofCoords(fieldOrdinal,d1) = dofCoords1(fieldOrdinal1,d1);
                               }
                               for (int d2=0; d2<spaceDim2; d2++)
                               {
                                 dofCoords(fieldOrdinal,spaceDim1+d2) = dofCoords2(fieldOrdinal2,d2);
                               }
                             }
                           });
    }
    
    /** \brief  Given "Dk" enumeration indices for the component bases, returns a Dk enumeration index for the composite basis.
        \param [in] dkEnum1         - Dk enumeration index for first component basis
        \param [in] operatorOrder1  - operator order for the first component basis
        \param [in] dkEnum2         - Dk enumeration index for second component basis
        \param [in] operatorOrder2  - operator order for the second component basis
     
        \return Dk enumeration index for the composite basis, corresponding to operator order operatorOrder1 + operatorOrder2.
     */
    ordinal_type getTensorDkEnumeration(ordinal_type dkEnum1, ordinal_type operatorOrder1,
                                        ordinal_type dkEnum2, ordinal_type operatorOrder2) const
    {
      unsigned spaceDim1 = basis1_.getBaseCellTopology().getDimension();
      unsigned spaceDim2 = basis2_.getBaseCellTopology().getDimension();
      
      // for now, we only support total spaceDim <= 3.  It would not be too hard to extend to support higher dimensions,
      // but the support needs to be built out in e.g. shards::CellTopology for this, as well as our DkEnumeration, etc.
      switch (spaceDim1)
      {
        case 1:
          switch (spaceDim2)
        {
          case 1:
            return getDkTensorIndex<1, 1>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          case 2:
            return getDkTensorIndex<1, 2>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported dimension combination");
        }
        case 2:
          switch (spaceDim2)
        {
          case 1:
            return getDkTensorIndex<2, 1>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported dimension combination");
        }
          //        case 3:
          //          switch (spaceDim2)
          //        {
          //          case 1:
          //            return getDkTensorIndex<3, 1>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          //          case 2:
          //            return getDkTensorIndex<3, 2>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          //          case 3:
          //            return getDkTensorIndex<3, 3>(dkEnum1, operatorOrder1, dkEnum2, operatorOrder2);
          //          default:
          //            INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported dimension combination");
          //        }
        default:
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported dimension combination");
      }
    }
    
    /** \brief  Evaluation of a FEM basis on a <strong>reference cell</strong>.

        Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
        points in the <strong>reference cell</strong> for which the basis is defined.

        \param  outputValues      [out] - variable rank array with the basis values
        \param  inputPoints       [in]  - rank-2 array (P,D) with the evaluation points
        \param  operatorType      [in]  - the operator acting on the basis functions

        \remark For rank and dimension specifications of the output array see Section
        \ref basis_md_array_sec.  Dimensions of <var>ArrayScalar</var> arguments are checked
        at runtime if HAVE_INTREPID2_DEBUG is defined.

        \remark A FEM basis spans a COMPLETE or INCOMPLETE polynomial space on the reference cell
        which is a smooth function space. Thus, all operator types that are meaningful for the
        approximated function space are admissible. When the order of the operator exceeds the
        degree of the basis, the output array is filled with the appropriate number of zeros.
    */
    void getValues( OutputViewType outputValues, const PointViewType  inputPoints,
                   const EOperator operatorType = OPERATOR_VALUE ) const override
    {
      bool tensorPoints;  // true would mean that we take the tensor product of inputPoints1 and inputPoints2 (and that this would be equivalent to inputPoints as given -- i.e., inputPoints1 and inputPoints2 would be a tensor decomposition of inputPoints)
      bool attemptTensorDecomposition = false; // support for this not yet implemented
      PointViewType inputPoints1, inputPoints2;
      getComponentPoints(inputPoints, attemptTensorDecomposition, inputPoints1, inputPoints2, tensorPoints);
      
      switch (operatorType)
      {
        case OPERATOR_D1:
        case OPERATOR_D2:
        case OPERATOR_D3:
        case OPERATOR_D4:
        case OPERATOR_D5:
        case OPERATOR_D6:
        case OPERATOR_D7:
        case OPERATOR_D8:
        case OPERATOR_D9:
        case OPERATOR_D10:
        {
          auto opOrder = getOperatorOrder(operatorType); // number of derivatives that we take in total
          // the Dk enumeration happens in lexicographic order (reading from left to right: x, y, z, etc.)
          // this governs the nesting order of the dkEnum1, dkEnum2 for loops below: dkEnum2 should increment fastest.
          for (int derivativeCountComp2=0; derivativeCountComp2<=opOrder; derivativeCountComp2++)
          {
            int derivativeCountComp1=opOrder-derivativeCountComp2;
            EOperator op1 = (derivativeCountComp1 == 0) ? OPERATOR_VALUE : EOperator(OPERATOR_D1 + (derivativeCountComp1 - 1));
            EOperator op2 = (derivativeCountComp2 == 0) ? OPERATOR_VALUE : EOperator(OPERATOR_D1 + (derivativeCountComp2 - 1));
            
            int spaceDim1 = inputPoints1.extent_int(1);
            int spaceDim2 = inputPoints2.extent_int(1);
            
            int dkCardinality1 = (op1 != OPERATOR_VALUE) ? getDkCardinality(op1, spaceDim1) : 1;
            int dkCardinality2 = (op2 != OPERATOR_VALUE) ? getDkCardinality(op2, spaceDim2) : 1;
            
            int basisCardinality1 = basis1_.getCardinality();
            int basisCardinality2 = basis2_.getCardinality();
            
            int totalPointCount = tensorPoints ? inputPoints1.extent_int(0) * inputPoints2.extent_int(0) : inputPoints1.extent_int(0);
            
            int pointCount1, pointCount2;
            if (tensorPoints)
            {
              pointCount1 = inputPoints1.extent_int(0);
              pointCount2 = inputPoints2.extent_int(0);
            }
            else
            {
              pointCount1 = totalPointCount;
              pointCount2 = totalPointCount;
            }
            
            OutputViewType outputValues1, outputValues2;
            if (op1 == OPERATOR_VALUE)
              outputValues1 = getMatchingViewWithLabel(outputValues, "output values - basis 1",basisCardinality1,pointCount1);
            else
              outputValues1 = getMatchingViewWithLabel(outputValues, "output values - basis 1",basisCardinality1,pointCount1,dkCardinality1);
            
            if (op2 == OPERATOR_VALUE)
              outputValues2 = getMatchingViewWithLabel(outputValues, "output values - basis 2",basisCardinality2,pointCount2);
            else
              outputValues2 = getMatchingViewWithLabel(outputValues, "output values - basis 2",basisCardinality2,pointCount2,dkCardinality2);
              
            basis1_.getValues(outputValues1,inputPoints1,op1);
            basis2_.getValues(outputValues2,inputPoints2,op2);
            
            const int outputVectorSize = getVectorSizeForHierarchicalParallelism<OutputValueType>();
            const int pointVectorSize  = getVectorSizeForHierarchicalParallelism<PointValueType>();
            const int vectorSize = std::max(outputVectorSize,pointVectorSize);
            
            auto policy = Kokkos::TeamPolicy<ExecutionSpace>(basisCardinality1,Kokkos::AUTO(),vectorSize);
            
            double weight = 1.0;
            using FunctorType = TensorViewFunctor<ExecutionSpace, OutputValueType, OutputViewType>;
            
            for (int dkEnum1=0; dkEnum1<dkCardinality1; dkEnum1++)
            {
              auto outputValues1_dkEnum1 = (op1 != OPERATOR_VALUE) ? Kokkos::subview(outputValues1,Kokkos::ALL(),Kokkos::ALL(),dkEnum1)
              : Kokkos::subview(outputValues1,Kokkos::ALL(),Kokkos::ALL());
              for (int dkEnum2=0; dkEnum2<dkCardinality2; dkEnum2++)
              {
                auto outputValues2_dkEnum2 = (op2 != OPERATOR_VALUE) ? Kokkos::subview(outputValues2,Kokkos::ALL(),Kokkos::ALL(),dkEnum2)
                : Kokkos::subview(outputValues2,Kokkos::ALL(),Kokkos::ALL());
                
                ordinal_type dkTensorIndex = getTensorDkEnumeration(dkEnum1, derivativeCountComp1, dkEnum2, derivativeCountComp2);
                auto outputValues_dkTensor = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),dkTensorIndex);
                // Note that there may be performance optimizations available here:
                // - could eliminate interior for loop in favor of having a vector-valued outputValues1_dk
                // - could add support to TensorViewFunctor (and probably TensorViewIterator) for this kind of tensor Dk type of traversal
                //   (this would allow us to eliminate both for loops here)
                // At the moment, we defer such optimizations on the idea that this may not ever become a performance bottleneck.
                FunctorType functor(outputValues_dkTensor, outputValues1_dkEnum1, outputValues2_dkEnum2, tensorPoints, weight);
                Kokkos::parallel_for( policy , functor, "TensorViewFunctor");
              }
            }
          }
        }
          break;
        default: // non-OPERATOR_Dn case must be handled by subclass.
          this->getValues(outputValues, operatorType, inputPoints1, inputPoints2, tensorPoints);
      }
    }
    
    /** \brief  Evaluation of a tensor FEM basis on a <strong>reference cell</strong>; subclasses should override this.

        Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
        points in the <strong>reference cell</strong> for which the basis is defined.

        \param  outputValues      [out] - variable rank array with the basis values
        \param  operatorType      [in]  - the operator acting on the basis functions
        \param  inputPoints1      [in]  - rank-2 array (P1,D1) with the evaluation points for basis1
        \param  inputPoints2      [in]  - rank-2 array (P2,D2) with the evaluation points for basis2
        \param  tensorPoints      [in]  - whether the points should be interpreted as tensor components of the evaluation points, or in a one-to-one correspondence

       Subclasses should override this method; this gives them an opportunity to specify how operatorType should be decomposed into operators on the component bases.
     
       If tensorPoints is true, then the points dimension of outputValues should be (P1*P2).
       If tensorPoints is false, then P1 should equal P2, and these should match the points dimension of outputValues.
     
     There are three variants of getValues:
     1. The three-argument version defined by Intrepid2::Basis.  TensorBasis provides an implementation of this, which calls the five-argument version (this one).
     2. The five-argument version (this method), which provides separate point sets for the component bases, and must be specified by subclasses.  Typical implementations call the seven-argument version.
     3. The seven-argument version (below), implemented by TensorBasis, which provides separate point sets and operators for the component bases, as well as an optional weight.
     
     The intent is that subclasses implement this five-argument version; in that implementation, they need to do little else than call the seven-argument version below.
     
     Note that the three-argument implementation handles the OPERATOR_Dn operators directly; that is, subclasses can omit any consideration of OPERATOR_Dn operators in their implementation of the five-argument version.
    */
    virtual void getValues(OutputViewType outputValues, const EOperator operatorType,
                           const PointViewType  inputPoints1, const PointViewType  inputPoints2,
                           bool tensorPoints) const
    {
      INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "one-operator, two-inputPoints getValues should be overridden by TensorBasis subclasses");
    }
    
    /** \brief  Evaluation of a tensor FEM basis on a <strong>reference cell</strong>; subclasses should override this.

        Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
        points in the <strong>reference cell</strong> for which the basis is defined.

        \param  outputValues      [out] - variable rank array with the basis values
        \param  inputPoints1      [in]  - rank-2 array (P1,D1) with the evaluation points for basis1
        \param  operatorType1     [in]  - the operator acting on basis1
        \param  inputPoints2      [in]  - rank-2 array (P2,D2) with the evaluation points for basis2
        \param  operatorType2     [in]  - the operator acting on basis2
        \param  tensorPoints      [in]  - whether the points should be interpreted as tensor components of the evaluation points, or in a one-to-one correspondence
        \param  weight            [in]  - optional weight (typically 1.0 or -1.0)
     
       If tensorPoints is true, then the points dimension of outputValues should be (P1*P2).
       If tensorPoints is false, then P1 should equal P2, and these should match the points dimension of outputValues.
     
     There are three variants of getValues:
     1. The three-argument version defined by Intrepid2::Basis.  TensorBasis provides an implementation of this, which calls the five-argument version (this one).
     2. The five-argument version (above), which provides separate point sets for the component bases, and must be specified by subclasses.  Typical implementations call the seven-argument version.
     3. The seven-argument version (this method), implemented by TensorBasis, which provides separate point sets and operators for the component bases, as well as an optional weight.
     
     Subclasses should override the five-argument version above; in their implementation, they need to do little else than call this seven-argument version.
    */
    void getValues( OutputViewType outputValues,
                   const PointViewType  inputPoints1, const EOperator operatorType1,
                   const PointViewType  inputPoints2, const EOperator operatorType2,
                   bool tensorPoints, double weight=1.0) const
    {
      int basisCardinality1 = basis1_.getCardinality();
      int basisCardinality2 = basis2_.getCardinality();
      
      int totalPointCount = tensorPoints ? inputPoints1.extent_int(0) * inputPoints2.extent_int(0) : inputPoints1.extent_int(0);
      
      int pointCount1, pointCount2;
      if (tensorPoints)
      {
        pointCount1 = inputPoints1.extent_int(0);
        pointCount2 = inputPoints2.extent_int(0);
      }
      else
      {
        pointCount1 = totalPointCount;
        pointCount2 = totalPointCount;
      }
      
      int spaceDim1 = inputPoints1.extent_int(1);
      int spaceDim2 = inputPoints2.extent_int(1);
      
      INTREPID2_TEST_FOR_EXCEPTION(!tensorPoints && (totalPointCount != inputPoints2.extent_int(0)),
                                   std::invalid_argument, "If tensorPoints is false, the point counts must match!");
            
      int opRank1 = getOperatorRank(basis1_.getFunctionSpace(), operatorType1, spaceDim1);
      int opRank2 = getOperatorRank(basis2_.getFunctionSpace(), operatorType2, spaceDim2);
      
      OutputViewType outputValues1, outputValues2;
      if (opRank1 == 0)
      {
        outputValues1 = getMatchingViewWithLabel(outputValues,"output values - basis 1",basisCardinality1,pointCount1);
      }
      else if (opRank1 == 1)
      {
        outputValues1 = getMatchingViewWithLabel(outputValues,"output values - basis 1",basisCardinality1,pointCount1,spaceDim1);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported opRank1");
      }
      
      if (opRank2 == 0)
      {
        outputValues2 = getMatchingViewWithLabel(outputValues,"output values - basis 2",basisCardinality2,pointCount2);
      }
      else if (opRank2 == 1)
      {
        outputValues2 = getMatchingViewWithLabel(outputValues,"output values - basis 2",basisCardinality2,pointCount2,spaceDim2);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported opRank2");
      }
      
      basis1_.getValues(outputValues1,inputPoints1,operatorType1);
      basis2_.getValues(outputValues2,inputPoints2,operatorType2);
      
      const int outputVectorSize = getVectorSizeForHierarchicalParallelism<OutputValueType>();
      const int pointVectorSize  = getVectorSizeForHierarchicalParallelism<PointValueType>();
      const int vectorSize = std::max(outputVectorSize,pointVectorSize);
      
      auto policy = Kokkos::TeamPolicy<ExecutionSpace>(basisCardinality1,Kokkos::AUTO(),vectorSize);
      
      using FunctorType = TensorViewFunctor<ExecutionSpace, OutputValueType, OutputViewType>;
      
      FunctorType functor(outputValues, outputValues1, outputValues2, tensorPoints, weight);
      Kokkos::parallel_for( policy , functor, "TensorViewFunctor");
    }
  };
  
  /** \struct Intrepid2::TensorBasis3_Functor
      \brief  Functor for computing values for the TensorBasis3 class.
   
   This functor is not intended for use outside of \ref Intrepid2::Basis_TensorBasis3.
   
   We may replace usage of this functor with TensorViewFunctor in the future.  This would likely allow more TensorBasis3 use cases.
  */
  template<class ExecutionSpace, class OutputScalar, class OutputFieldType>
  struct TensorBasis3_Functor
  {
    using ScratchSpace       = Kokkos::DefaultExecutionSpace::scratch_memory_space;
    using OutputScratchView  = Kokkos::View<OutputScalar*,ScratchSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    
    using TeamPolicy = Kokkos::TeamPolicy<>;
    using TeamMember = TeamPolicy::member_type;
    
    OutputFieldType  output_; // F,P
    OutputFieldType  input1_; // F1,P[,D] or F1,P1[,D]
    OutputFieldType  input2_; // F2,P[,D] or F2,P2[,D]
    OutputFieldType  input3_; // F2,P[,D] or F2,P2[,D]
    
    int numFields_, numPoints_;
    int numFields1_, numPoints1_;
    int numFields2_, numPoints2_;
    int numFields3_, numPoints3_;
    
    bool tensorPoints_; // if true, input1, input2, input3 refer to values at decomposed points, and P = P1 * P2 * P3.  If false, then the three inputs refer to points in the full-dimensional space, and their point lengths are the same as that of the final output.
    
    double weight_;
    
    TensorBasis3_Functor(OutputFieldType output, OutputFieldType inputValues1, OutputFieldType inputValues2, OutputFieldType inputValues3,
                         bool tensorPoints, double weight)
    : output_(output), input1_(inputValues1), input2_(inputValues2), input3_(inputValues3), tensorPoints_(tensorPoints), weight_(weight)
    {
      numFields_ = output.extent_int(0);
      numPoints_ = output.extent_int(1);
      
      numFields1_ = inputValues1.extent_int(0);
      numPoints1_ = inputValues1.extent_int(1);
      
      numFields2_ = inputValues2.extent_int(0);
      numPoints2_ = inputValues2.extent_int(1);
      
      numFields3_ = inputValues3.extent_int(0);
      numPoints3_ = inputValues3.extent_int(1);
      /*
       We don't yet support tensor-valued bases here (only vector and scalar).  The main design question is how the layouts
       of the input containers relates to the layout of the output container.  The work we've done in TensorViewIterator basically
       shows the choices that can be made.  It does appear that in most cases (at least (most of?) those supported by TensorViewIterator),
       we can infer from the dimensions of input/output containers what choice should be made in each dimension.
       */
      INTREPID2_TEST_FOR_EXCEPTION(inputValues1.rank() > 3, std::invalid_argument, "ranks greater than 3 not yet supported");
      INTREPID2_TEST_FOR_EXCEPTION(inputValues2.rank() > 3, std::invalid_argument, "ranks greater than 3 not yet supported");
      INTREPID2_TEST_FOR_EXCEPTION(inputValues3.rank() > 3, std::invalid_argument, "ranks greater than 3 not yet supported");
      INTREPID2_TEST_FOR_EXCEPTION((inputValues1.rank() == 3) && (inputValues2.rank() == 3), std::invalid_argument, "two vector-valued input ranks not yet supported");
      INTREPID2_TEST_FOR_EXCEPTION((inputValues1.rank() == 3) && (inputValues3.rank() == 3), std::invalid_argument, "two vector-valued input ranks not yet supported");
      INTREPID2_TEST_FOR_EXCEPTION((inputValues2.rank() == 3) && (inputValues3.rank() == 3), std::invalid_argument, "two vector-valued input ranks not yet supported");
      
      if (!tensorPoints_)
      {
        // then the point counts should all match
        INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != numPoints1_, std::invalid_argument, "incompatible point counts");
        INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != numPoints2_, std::invalid_argument, "incompatible point counts");
        INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != numPoints3_, std::invalid_argument, "incompatible point counts");
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != numPoints1_ * numPoints2_ * numPoints3_, std::invalid_argument, "incompatible point counts");
      }
      
      INTREPID2_TEST_FOR_EXCEPTION(numFields_ != numFields1_ * numFields2_ * numFields3_, std::invalid_argument, "incompatible field sizes");
    }
    
    KOKKOS_INLINE_FUNCTION
    void operator()( const TeamMember & teamMember ) const
    {
      auto fieldOrdinal1 = teamMember.league_rank();
      
      if (!tensorPoints_)
      {
        if ((input1_.rank() == 2) && (input2_.rank() == 2))
        {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numFields2_), [&] (const int& fieldOrdinal2) {
            for (int fieldOrdinal3=0; fieldOrdinal3 < numFields3_; fieldOrdinal3++)
            {
              int fieldOrdinal = (fieldOrdinal3 * numFields2_ + fieldOrdinal2) * numFields1_ + fieldOrdinal1;
              for (int pointOrdinal=0; pointOrdinal<numPoints_; pointOrdinal++)
              {
                output_(fieldOrdinal,pointOrdinal) = weight_ * input1_(fieldOrdinal1,pointOrdinal) * input2_(fieldOrdinal2,pointOrdinal) * input3_(fieldOrdinal3,pointOrdinal);
              }
            }
          });
        }
        else if (input1_.rank() == 3)
        {
          int spaceDim = input1_.extent_int(2);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numFields2_), [&] (const int& fieldOrdinal2) {
            for (int fieldOrdinal3=0; fieldOrdinal3 < numFields3_; fieldOrdinal3++)
            {
              int fieldOrdinal = (fieldOrdinal3 * numFields2_ + fieldOrdinal2) * numFields1_ + fieldOrdinal1;
              for (int pointOrdinal=0; pointOrdinal<numPoints_; pointOrdinal++)
              {
                for (int d=0; d<spaceDim; d++)
                {
                  output_(fieldOrdinal,pointOrdinal,d) = weight_ * input1_(fieldOrdinal1,pointOrdinal,d) * input2_(fieldOrdinal2,pointOrdinal) * input3_(fieldOrdinal3,pointOrdinal);
                }
              }
            }
          });
        }
        else if (input2_.rank() == 3)
        {
          int spaceDim = input2_.extent_int(2);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numFields2_), [&] (const int& fieldOrdinal2) {
            for (int fieldOrdinal3=0; fieldOrdinal3 < numFields3_; fieldOrdinal3++)
            {
              int fieldOrdinal = (fieldOrdinal3 * numFields2_ + fieldOrdinal2) * numFields1_ + fieldOrdinal1;
              for (int pointOrdinal=0; pointOrdinal<numPoints_; pointOrdinal++)
              {
                for (int d=0; d<spaceDim; d++)
                {
                  output_(fieldOrdinal,pointOrdinal,d) = weight_ * input1_(fieldOrdinal1,pointOrdinal) * input2_(fieldOrdinal2,pointOrdinal) * input3_(fieldOrdinal3,pointOrdinal,d);
                }
              }
            }
          });
        }
        else if (input3_.rank() == 3)
        {
          int spaceDim = input3_.extent_int(2);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numFields2_), [&] (const int& fieldOrdinal2) {
            for (int fieldOrdinal3=0; fieldOrdinal3 < numFields3_; fieldOrdinal3++)
            {
              int fieldOrdinal = (fieldOrdinal3 * numFields2_ + fieldOrdinal2) * numFields1_ + fieldOrdinal1;
              for (int pointOrdinal=0; pointOrdinal<numPoints_; pointOrdinal++)
              {
                for (int d=0; d<spaceDim; d++)
                {
                  output_(fieldOrdinal,pointOrdinal,d) = weight_ * input1_(fieldOrdinal1,pointOrdinal) * input2_(fieldOrdinal2,pointOrdinal) * input3_(fieldOrdinal3,pointOrdinal,d);
                }
              }
            }
          });
        }
        else
        {
          // unsupported rank combination -- enforced in constructor
        }
      }
      else
      {
        if ((input1_.rank() == 2) && (input2_.rank() == 2))
        {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numFields2_), [&] (const int& fieldOrdinal2) {
            for (int fieldOrdinal3=0; fieldOrdinal3 < numFields3_; fieldOrdinal3++)
            {
              int fieldOrdinal = (fieldOrdinal3 * numFields2_ + fieldOrdinal2) * numFields1_ + fieldOrdinal1;
              for (int pointOrdinal3=0; pointOrdinal3<numPoints3_; pointOrdinal3++)
              {
                for (int pointOrdinal2=0; pointOrdinal2<numPoints2_; pointOrdinal2++)
                {
                  for (int pointOrdinal1=0; pointOrdinal1<numPoints1_; pointOrdinal1++)
                  {
                    int pointOrdinal = (pointOrdinal3 * numPoints2_ + pointOrdinal2) * numPoints1_ + pointOrdinal1;
                    output_(fieldOrdinal,pointOrdinal) = weight_ * input1_(fieldOrdinal1,pointOrdinal1) * input2_(fieldOrdinal2,pointOrdinal2) * input3_(fieldOrdinal3,pointOrdinal3);
                  }
                }
              }
            }
          });
        }
        else if (input1_.rank() == 3) // based on constructor requirements, this means the others are rank 2
        {
          int spaceDim = input1_.extent_int(2);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numFields2_), [&] (const int& fieldOrdinal2) {
            for (int fieldOrdinal3=0; fieldOrdinal3 < numFields3_; fieldOrdinal3++)
            {
              int fieldOrdinal = (fieldOrdinal3 * numFields2_ + fieldOrdinal2) * numFields1_ + fieldOrdinal1;
              for (int pointOrdinal3=0; pointOrdinal3<numPoints3_; pointOrdinal3++)
              {
                for (int pointOrdinal2=0; pointOrdinal2<numPoints2_; pointOrdinal2++)
                {
                  for (int pointOrdinal1=0; pointOrdinal1<numPoints1_; pointOrdinal1++)
                  {
                    int pointOrdinal = (pointOrdinal3 * numPoints2_ + pointOrdinal2) * numPoints1_ + pointOrdinal1;
                    for (int d=0; d<spaceDim; d++)
                    {
                      output_(fieldOrdinal,pointOrdinal,d) = weight_ * input1_(fieldOrdinal1,pointOrdinal1,d) * input2_(fieldOrdinal2,pointOrdinal2) * input3_(fieldOrdinal3,pointOrdinal3);
                    }
                  }
                }
              }
            }
          });
        }
        else if (input2_.rank() == 3) // based on constructor requirements, this means the others are rank 2
        {
          int spaceDim = input2_.extent_int(2);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numFields2_), [&] (const int& fieldOrdinal2) {
            for (int fieldOrdinal3=0; fieldOrdinal3 < numFields3_; fieldOrdinal3++)
            {
              int fieldOrdinal = (fieldOrdinal3 * numFields2_ + fieldOrdinal2) * numFields1_ + fieldOrdinal1;
              for (int pointOrdinal3=0; pointOrdinal3<numPoints3_; pointOrdinal3++)
              {
                for (int pointOrdinal2=0; pointOrdinal2<numPoints2_; pointOrdinal2++)
                {
                  for (int pointOrdinal1=0; pointOrdinal1<numPoints1_; pointOrdinal1++)
                  {
                    int pointOrdinal = (pointOrdinal3 * numPoints2_ + pointOrdinal2) * numPoints1_ + pointOrdinal1;
                    for (int d=0; d<spaceDim; d++)
                    {
                      output_(fieldOrdinal,pointOrdinal,d) = weight_ * input1_(fieldOrdinal1,pointOrdinal1) * input2_(fieldOrdinal2,pointOrdinal2,d) * input3_(fieldOrdinal3,pointOrdinal3);
                    }
                  }
                }
              }
            }
          });
        }
        else if (input3_.rank() == 3) // based on constructor requirements, this means the others are rank 2
        {
          int spaceDim = input3_.extent_int(2);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,0,numFields2_), [&] (const int& fieldOrdinal2) {
            for (int fieldOrdinal3=0; fieldOrdinal3 < numFields3_; fieldOrdinal3++)
            {
              int fieldOrdinal = (fieldOrdinal3 * numFields2_ + fieldOrdinal2) * numFields1_ + fieldOrdinal1;
              for (int pointOrdinal3=0; pointOrdinal3<numPoints3_; pointOrdinal3++)
              {
                for (int pointOrdinal2=0; pointOrdinal2<numPoints2_; pointOrdinal2++)
                {
                  for (int pointOrdinal1=0; pointOrdinal1<numPoints1_; pointOrdinal1++)
                  {
                    int pointOrdinal = (pointOrdinal3 * numPoints2_ + pointOrdinal2) * numPoints1_ + pointOrdinal1;
                    for (int d=0; d<spaceDim; d++)
                    {
                      output_(fieldOrdinal,pointOrdinal,d) = weight_ * input1_(fieldOrdinal1,pointOrdinal1) * input2_(fieldOrdinal2,pointOrdinal2) * input3_(fieldOrdinal3,pointOrdinal3,d);
                    }
                  }
                }
              }
            }
          });
        }
        else
        {
          // unsupported rank combination -- enforced in constructor
        }
      }
    }
  };
  
  
  template<typename Basis1, typename Basis2, typename Basis3,
  typename ExecutionSpace=Kokkos::DefaultExecutionSpace,
  typename OutputScalar = double,
  typename PointScalar  = double>
  class Basis_TensorBasis3
  : public Basis_TensorBasis< Basis_TensorBasis<Basis1,Basis2>,
  Basis3>
  {
    using Basis12 = Basis_TensorBasis<Basis1,Basis2>;
    using TensorBasis123 = Basis_TensorBasis<Basis12,Basis3>;
    
  public:
    using OutputViewType = typename TensorBasis123::OutputViewType;
    using PointViewType  = typename TensorBasis123::PointViewType;
    using ScalarViewType = typename TensorBasis123::ScalarViewType;
    
    using OutputValueType = typename TensorBasis123::OutputValueType;
    using PointValueType  = typename TensorBasis123::PointValueType;
  protected:
    Basis1 basis1_;
    Basis2 basis2_;
    Basis3 basis3_;
  public:
    Basis_TensorBasis3(Basis1 basis1, Basis2 basis2, Basis3 basis3)
    :
    TensorBasis123(Basis12(basis1,basis2),basis3),
    basis1_(basis1),
    basis2_(basis2),
    basis3_(basis3)
    {}
    
    using TensorBasis123::getValues;
    
    /** \brief  Evaluation of a tensor FEM basis on a <strong>reference cell</strong>.

        Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
        points in the <strong>reference cell</strong> for which the basis is defined.

        \param  outputValues      [out] - variable rank array with the basis values
        \param  operatorType      [in]  - the operator acting on the basis functions
        \param  inputPoints12     [in]  - rank-2 array (P12,D12) with the evaluation points for basis12
        \param  inputPoints3      [in]  - rank-2 array (P3,D3) with the evaluation points for basis3
        \param  tensorPoints      [in]  - whether the points should be interpreted as tensor components of the evaluation points, or in a one-to-one correspondence
     
       If tensorPoints is true, then the points dimension of outputValues should be (P12*P3).
       If tensorPoints is false, then P12 should equal P3, and these should match the points dimension of outputValues.
     
     There are four variants of getValues:
     1. The three-argument version defined by Intrepid2::Basis.  TensorBasis provides an implementation of this, which calls the five-argument version (this one).
     2. The five-argument version (this method), which provides partially separated point sets for the component bases.  TensorBasis3 provides an implementation of this, which calls the six-argument version.
     3. The six-argument version, which fully separates the point sets for the component bases.  Subclasses should implement this; they essentially specify the decomposition of the operator.
     4. The nine-argument version (below), implemented by TensorBasis3, which provides separate point sets and operators for the component bases, as well as an optional weight.
     
     The intent is that subclasses implement the six-argument version; in that implementation, they need to do little else than call the nine-argument version below.
     
     Note that the three-argument implementation handles the OPERATOR_Dn operators directly; that is, subclasses can omit any consideration of OPERATOR_Dn operators in their implementation of the six-argument version.
    */
    virtual void getValues(OutputViewType outputValues, const EOperator operatorType,
                           const PointViewType inputPoints12, const PointViewType  inputPoints3,
                           bool tensorPoints) const override
    {
      // TODO: rework this to use superclass's getComponentPoints.
      
      int spaceDim1 = basis1_.getBaseCellTopology().getDimension();
      int spaceDim2 = basis2_.getBaseCellTopology().getDimension();
      
      int totalSpaceDim12 = inputPoints12.extent_int(1);
      
      TEUCHOS_ASSERT(spaceDim1 + spaceDim2 == totalSpaceDim12);
      
      if (!tensorPoints)
      {
        auto inputPoints1 = Kokkos::subview(inputPoints12,Kokkos::ALL(),std::make_pair(0,spaceDim1));
        auto inputPoints2 = Kokkos::subview(inputPoints12,Kokkos::ALL(),std::make_pair(spaceDim1,totalSpaceDim12));
        
        this->getValues(outputValues, operatorType, inputPoints1, inputPoints2, inputPoints3, tensorPoints);
      }
      else
      {
        // superclass doesn't (yet) have a clever way to detect tensor points in a single container
        // we'd need something along those lines here to detect them in inputPoints12.
        // if we do add such a mechanism to superclass, it should be simple enough to call that from here
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "This method does not yet handle tensorPoints=true");
      }
    }
    
    /** \brief  Evaluation of a tensor FEM basis on a <strong>reference cell</strong>; subclasses should override this.

        Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
        points in the <strong>reference cell</strong> for which the basis is defined.

        \param  outputValues      [out] - variable rank array with the basis values
        \param  operatorType      [in]  - the operator acting on the basis functions
        \param  inputPoints1      [in]  - rank-2 array (P1,D1) with the evaluation points for basis1
        \param  inputPoints1      [in]  - rank-2 array (P2,D2) with the evaluation points for basis2
        \param  inputPoints3      [in]  - rank-2 array (P3,D3) with the evaluation points for basis3
        \param  tensorPoints      [in]  - whether the points should be interpreted as tensor components of the evaluation points, or in a one-to-one correspondence

       Subclasses should override this method; this gives them an opportunity to specify how operatorType should be decomposed into operators on the component bases.
     
       If tensorPoints is true, then the points dimension of outputValues should be (P1*P2*P3).
       If tensorPoints is false, then P1 should equal P2 and P2 should equal P3, and these should match the points dimension of outputValues.
     
     There are four variants of getValues:
     1. The three-argument version defined by Intrepid2::Basis.  TensorBasis provides an implementation of this, which calls the five-argument version (this one).
     2. The five-argument version (above), which provides partially separated point sets for the component bases.  TensorBasis3 provides an implementation of this, which calls the six-argument version.
     3. The six-argument version (this method), which fully separates the point sets for the component bases.  Subclasses should implement this; they essentially specify the decomposition of the operator.
     4. The nine-argument version (below), implemented by TensorBasis3, which provides separate point sets and operators for the component bases, as well as an optional weight.
     
     The intent is that subclasses implement this six-argument version; in that implementation, they need to do little else than call the nine-argument version below.
     
     Note that the three-argument implementation handles the OPERATOR_Dn operators directly; that is, subclasses can omit any consideration of OPERATOR_Dn operators in their implementation of the five-argument version.
    */
    virtual void getValues(OutputViewType outputValues, const EOperator operatorType,
                           const PointViewType  inputPoints1, const PointViewType  inputPoints2, const PointViewType inputPoints3,
                           bool tensorPoints) const = 0;
    
    /** \brief  Evaluation of a tensor FEM basis on a <strong>reference cell</strong>; subclasses should override this.

        Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
        points in the <strong>reference cell</strong> for which the basis is defined.

        \param  outputValues      [out] - variable rank array with the basis values
        \param  inputPoints1      [in]  - rank-2 array (P1,D1) with the evaluation points for basis1
        \param  operatorType1     [in]  - the operator acting on basis1
        \param  inputPoints2      [in]  - rank-2 array (P2,D2) with the evaluation points for basis2
        \param  operatorType2     [in]  - the operator acting on basis2
        \param  inputPoints3      [in]  - rank-2 array (P3,D3) with the evaluation points for basis3
        \param  operatorType3     [in]  - the operator acting on basis3
        \param  tensorPoints      [in]  - whether the points should be interpreted as tensor components of the evaluation points, or in a one-to-one correspondence
     
       If tensorPoints is true, then the points dimension of outputValues should be (P1*P2*P3).
       If tensorPoints is false, then P1 should equal P2 and P2 should equal P3, and these should match the points dimension of outputValues.
     
     There are four variants of getValues:
     1. The three-argument version defined by Intrepid2::Basis.  TensorBasis provides an implementation of this, which calls the five-argument version (this one).
     2. The five-argument version (above), which provides partially separated point sets for the component bases.  TensorBasis3 provides an implementation of this, which calls the six-argument version.
     3. The six-argument version (this method), which fully separates the point sets for the component bases.  Subclasses should implement this; they essentially specify the decomposition of the operator.
     4. The nine-argument version (below), implemented by TensorBasis3, which provides separate point sets and operators for the component bases, as well as an optional weight.
     
     The intent is that subclasses implement this six-argument version; in that implementation, they need to do little else than call the nine-argument version below.
     
     Note that the three-argument implementation handles the OPERATOR_Dn operators directly; that is, subclasses can omit any consideration of OPERATOR_Dn operators in their implementation of the five-argument version.
    */
    void getValues( OutputViewType outputValues,
                   const PointViewType  inputPoints1, const EOperator operatorType1,
                   const PointViewType  inputPoints2, const EOperator operatorType2,
                   const PointViewType  inputPoints3, const EOperator operatorType3,
                   bool tensorPoints, double weight=1.0) const
    {
      int basisCardinality1 = basis1_.getCardinality();
      int basisCardinality2 = basis2_.getCardinality();
      int basisCardinality3 = basis3_.getCardinality();
      
      int spaceDim1 = inputPoints1.extent_int(1);
      int spaceDim2 = inputPoints2.extent_int(1);
      int spaceDim3 = inputPoints3.extent_int(1);
      
      int totalPointCount;
      int pointCount1, pointCount2, pointCount3;
      if (tensorPoints)
      {
        pointCount1 = inputPoints1.extent_int(0);
        pointCount2 = inputPoints2.extent_int(0);
        pointCount3 = inputPoints3.extent_int(0);
        totalPointCount = pointCount1 * pointCount2 * pointCount3;
      }
      else
      {
        totalPointCount = inputPoints1.extent_int(0);
        pointCount1 = totalPointCount;
        pointCount2 = totalPointCount;
        pointCount3 = totalPointCount;
        
        INTREPID2_TEST_FOR_EXCEPTION((totalPointCount != inputPoints2.extent_int(0)) || (totalPointCount != inputPoints3.extent_int(0)),
                                     std::invalid_argument, "If tensorPoints is false, the point counts must match!");
      }
      
      // structure of this implementation:
      /*
       - allocate output1, output2, output3 containers
       - either:
       1. split off the tensor functor call into its own method in TensorBasis, and
       - call it once with output1, output2, placing these in another newly allocated output12, then
       - call it again with output12, output3
       OR
       2. create a 3-argument tensor functor and call it with output1,output2,output3
       
       At the moment, the 3-argument functor seems like a better approach.  It's likely more code, but somewhat
       more efficient and easier to understand/debug.  And the code is fairly straightforward to produce.
       */
      
      // copied from the 2-argument TensorBasis implementation:
      
      OutputViewType outputValues1, outputValues2, outputValues3;
      if (outputValues.rank() == 2)
      {
        outputValues1 = getMatchingViewWithLabel(outputValues,"output values - basis 1",basisCardinality1,pointCount1);
        outputValues2 = getMatchingViewWithLabel(outputValues,"output values - basis 2",basisCardinality2,pointCount2);
        outputValues3 = getMatchingViewWithLabel(outputValues,"output values - basis 3",basisCardinality3,pointCount3);
      }
      else if (outputValues.rank() == 3)
      {
        if (spaceDim1 == 1)
        {
          // use a rank 2 container for basis1
          outputValues1 = getMatchingViewWithLabel(outputValues,"output values - basis 1",basisCardinality1,pointCount1);
        }
        else
        {
          outputValues1 = getMatchingViewWithLabel(outputValues,"output values - basis 1",basisCardinality1,pointCount1,spaceDim1);
        }
        if (spaceDim2 == 1)
        {
          // use a rank 2 container for basis2
          outputValues2 = getMatchingViewWithLabel(outputValues,"output values - basis 2",basisCardinality2,pointCount2);
        }
        else
        {
          outputValues2 = getMatchingViewWithLabel(outputValues,"output values - basis 2",basisCardinality2,pointCount2,spaceDim2);
        }
        if (spaceDim3 == 1)
        {
          // use a rank 2 container for basis2
          outputValues3 = getMatchingViewWithLabel(outputValues,"output values - basis 3",basisCardinality3,pointCount3);
        }
        else
        {
          outputValues3 = getMatchingViewWithLabel(outputValues,"output values - basis 3",basisCardinality3,pointCount3,spaceDim3);
        }
      }
      
      basis1_.getValues(outputValues1,inputPoints1,operatorType1);
      basis2_.getValues(outputValues2,inputPoints2,operatorType2);
      basis3_.getValues(outputValues3,inputPoints3,operatorType3);
      
      const int outputVectorSize = getVectorSizeForHierarchicalParallelism<OutputScalar>();
      const int pointVectorSize  = getVectorSizeForHierarchicalParallelism<PointScalar>();
      const int vectorSize = std::max(outputVectorSize,pointVectorSize);
      
      auto policy = Kokkos::TeamPolicy<ExecutionSpace>(basisCardinality1,Kokkos::AUTO(),vectorSize);
      
      using FunctorType = TensorBasis3_Functor<ExecutionSpace, OutputScalar, OutputViewType>;
      FunctorType functor(outputValues, outputValues1, outputValues2, outputValues3, tensorPoints, weight);
      Kokkos::parallel_for( policy , functor, "TensorBasis3_Functor");
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_TensorBasis_h */
