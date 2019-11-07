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

/** \file   Intrepid2_IntegratedLegendreBasis_HGRAD_LINE.hpp
    \brief  H(grad) basis on the line based on integrated Legendre polynomials.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_IntegratedLegendreBasis_HGRAD_LINE_h
#define Intrepid2_IntegratedLegendreBasis_HGRAD_LINE_h

#include <Kokkos_View.hpp>
#include <Kokkos_DynRankView.hpp>

#include <Intrepid2_config.h>

#include "Intrepid2_Polynomials.hpp"
#include "Intrepid2_Utils.hpp"

namespace Intrepid2
{
  /** \class  Intrepid2::Hierarchical_HGRAD_LINE_Functor
      \brief  Functor for computing values for the IntegratedLegendreBasis_HGRAD_LINE class.
   
   This functor is not intended for use outside of IntegratedLegendreBasis_HGRAD_LINE.
  */
  template<class ExecutionSpace, class OutputScalar, class PointScalar,
           class OutputFieldType, class InputPointsType>
  struct Hierarchical_HGRAD_LINE_Functor
  {
    using ScratchSpace       = Kokkos::DefaultExecutionSpace::scratch_memory_space;
    using OutputScratchView  = Kokkos::View<OutputScalar*,ScratchSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using PointScratchView   = Kokkos::View<PointScalar*, ScratchSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    
    using TeamPolicy = Kokkos::TeamPolicy<>;
    using TeamMember = TeamPolicy::member_type;
    
    EOperator opType_; // OPERATOR_VALUE or OPERATOR_GRAD
    
    OutputFieldType  output_;      // F,P
    InputPointsType  inputPoints_; // P,D
    
    int polyOrder_;
    bool defineVertexFunctions_;
    int numFields_, numPoints_;
    
    size_t fad_size_output_;
    
    Hierarchical_HGRAD_LINE_Functor(EOperator opType, OutputFieldType output, InputPointsType inputPoints,
                                    int polyOrder, bool defineVertexFunctions)
    : opType_(opType), output_(output), inputPoints_(inputPoints),
      polyOrder_(polyOrder), defineVertexFunctions_(defineVertexFunctions),
      fad_size_output_(getScalarDimensionForView(output))
    {
      numFields_ = output.extent_int(0);
      numPoints_ = output.extent_int(1);
      INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != inputPoints.extent_int(0), std::invalid_argument, "point counts need to match!");
      INTREPID2_TEST_FOR_EXCEPTION(numFields_ != polyOrder_+1, std::invalid_argument, "output field size does not match basis cardinality");
    }
    
    KOKKOS_INLINE_FUNCTION
    void operator()( const TeamMember & teamMember ) const
    {
      auto pointOrdinal = teamMember.league_rank();
      OutputScratchView field_values_at_point;
      if (fad_size_output_ > 0) {
        field_values_at_point = OutputScratchView(teamMember.team_shmem(), numFields_, fad_size_output_);
      }
      else {
        field_values_at_point = OutputScratchView(teamMember.team_shmem(), numFields_);
      }
      
      const auto & input_x = inputPoints_(pointOrdinal,0);
      const bool taking_derivative = (opType_ != OPERATOR_VALUE);
      const bool callingShiftedScaledLegendre = (opType_ == OPERATOR_VALUE) || (opType_ == OPERATOR_GRAD) || (opType_ == OPERATOR_D1);
      
      // shiftedScaledIntegratedLegendreValues{_dx} expects x in [0,1]
      const PointScalar x = callingShiftedScaledLegendre ? PointScalar((input_x + 1.0)/2.0) : PointScalar(input_x);
      const double legendreScaling = 1.0;
      const double outputScaling   = taking_derivative ? 0.5 : 1.0; // output scaling -- 0.5 if we take derivatives, 1.0 otherwise
      
      switch (opType_)
      {
        case OPERATOR_VALUE:
          // field values are integrated Legendre polynomials, except for the first and second field,
          // which may be 1 and x or x and 1-x, depending on whether the vertex compatibility flag is set.
          Polynomials::shiftedScaledIntegratedLegendreValues(field_values_at_point, polyOrder_, x, legendreScaling);
          
          // note that because shiftedScaledIntegratedLegendreValues determines field values recursively, there is not much
          // opportunity at that level for further parallelism
          
          if (defineVertexFunctions_)
          {
            field_values_at_point(0) = 1. - x;
            field_values_at_point(1) = x;
          }
          break;
        case OPERATOR_GRAD:
        case OPERATOR_D1:
          // field values are Legendre polynomials, except for the first and second field,
          // which may be 0 and 1 or -1 and 1, depending on whether the vertex compatibility flag is set.
          Polynomials::shiftedScaledIntegratedLegendreValues_dx(field_values_at_point, polyOrder_, x, legendreScaling);
          
          // note that because shiftedScaledIntegratedLegendreValues_dx determines field values recursively, there is not much
          // opportunity at that level for further parallelism
          
          if (defineVertexFunctions_)
          {
            field_values_at_point(0) = -1.0; // derivative of 1-x
            field_values_at_point(1) =  1.0; // derivative of x
          }
          break;
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
          auto derivativeOrder = getOperatorOrder(opType_) - 1;
          Polynomials::legendreDerivativeValues(field_values_at_point, polyOrder_, x, derivativeOrder);
          
          // L_i is defined in terms of an integral of P_(i-1), so we need to shift the values by 1
          if (numFields_ >= 3)
          {
            OutputScalar Pn_minus_one = field_values_at_point(1);
            for (int fieldOrdinal=2; fieldOrdinal<numFields_; fieldOrdinal++)
            {
              OutputScalar Pn = field_values_at_point(fieldOrdinal);
              field_values_at_point(fieldOrdinal) = Pn_minus_one;
              Pn_minus_one = Pn;
            }
          }
          if (numFields_ >= 1) field_values_at_point(0) = 0.0;
          if (numFields_ >= 2) field_values_at_point(1) = 0.0;
          // legendreDerivativeValues works on [-1,1], so no per-derivative scaling is necessary
          // however, there is a factor of 0.5 that comes from the scaling of the Legendre polynomials prior to integration
          // in the shiftedScaledIntegratedLegendreValues -- the first derivative of our integrated polynomials is 0.5 times the Legendre polynomial
          break;
        }
        default:
          // unsupported operator type
          device_assert(false);
      }
      
      // copy the values into the output container
      for (int fieldOrdinal=0; fieldOrdinal<numFields_; fieldOrdinal++)
      {
        // access() allows us to write one line that applies both to gradient (for which outputValues has rank 3, but third rank has only one entry) and to value (rank 2)
        output_.access(fieldOrdinal,pointOrdinal,0) = outputScaling * field_values_at_point(fieldOrdinal);
      }
    }
    
    // Provide the shared memory capacity.
    // This function takes the team_size as an argument,
    // which allows team_size-dependent allocations.
    size_t team_shmem_size (int team_size) const
    {
      // we want to use shared memory to create a fast buffer that we can use for basis computations
      size_t shmem_size = 0;
      if (fad_size_output_ > 0)
        shmem_size += OutputScratchView::shmem_size(numFields_, fad_size_output_);
      else
        shmem_size += OutputScratchView::shmem_size(numFields_);
      
      return shmem_size;
    }
  };
  
  /** \class  Intrepid2::IntegratedLegendreBasis_HGRAD_LINE
      \brief  Basis defining integrated Legendre basis on the line, a polynomial subspace of H(grad) on the line.

              This is used in the construction of hierarchical bases on higher-dimensional topologies.  For
              mathematical details of the construction, see:
   
               Federico Fuentes, Brendan Keith, Leszek Demkowicz, Sriram Nagaraj.
               "Orientation embedded high order shape functions for the exact sequence elements of all shapes."
               Computers & Mathematics with Applications, Volume 70, Issue 4, 2015, Pages 353-458, ISSN 0898-1221.
               https://doi.org/10.1016/j.camwa.2015.04.027.
   
               Note that the template argument defineVertexFunctions controls whether this basis is defined in a
               strictly hierarchical way.  If defineVertexFunctions is false, then the first basis function is the
               constant 1.0, and this basis will be suitable for discontinuous discretizations.  If defineVertexFunctions
               is true, then the first basis function will instead be 1.0-x, and the basis will be suitable for
               continuous discretizations.
  */
  template<typename ExecutionSpace=Kokkos::DefaultExecutionSpace,
           typename OutputScalar = double,
           typename PointScalar  = double,
           bool defineVertexFunctions = true,            // if defineVertexFunctions is true, first and second basis functions are x and 1-x.  Otherwise, they are 1 and x.
           bool useMinusOneToOneReferenceElement = true> // if useMinusOneToOneReferenceElement is true, basis is define on [-1,1].  Otherwise, [0,1].
  class IntegratedLegendreBasis_HGRAD_LINE
  : public Basis<ExecutionSpace,OutputScalar,PointScalar>
  {
  public:
    using OrdinalTypeArray1DHost = typename Basis<ExecutionSpace,OutputScalar,PointScalar>::ordinal_type_array_1d_host;
    using OrdinalTypeArray2DHost = typename Basis<ExecutionSpace,OutputScalar,PointScalar>::ordinal_type_array_2d_host;
    
    typedef typename Basis<ExecutionSpace,OutputScalar,PointScalar>::outputViewType outputViewType;
    typedef typename Basis<ExecutionSpace,OutputScalar,PointScalar>::pointViewType  pointViewType;
    typedef typename Basis<ExecutionSpace,OutputScalar,PointScalar>::scalarViewType scalarViewType;
  protected:
    int polyOrder_; // the maximum order of the polynomial
    bool defineVertexFunctions_; // if true, first and second basis functions are x and 1-x.  Otherwise, they are 1 and x.
  public:
    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order of the basis.
     
     The basis will have polyOrder + 1 members.
     
     If defineVertexFunctions is false, then all basis functions are identified with the interior of the line element, and the first two basis functions are 1 and x.
     
     If defineVertexFunctions is true, then the first two basis functions are 1-x and x, and these are identified with the left and right vertices of the cell.
     
     */
    IntegratedLegendreBasis_HGRAD_LINE(int polyOrder)
    :
    polyOrder_(polyOrder)
    {
      this->basisCardinality_  = polyOrder+1;
      this->basisDegree_       = polyOrder;
      this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >() );
      this->basisType_         = BASIS_FEM_HIERARCHICAL;
      this->basisCoordinates_  = COORDINATES_CARTESIAN;
      this->functionSpace_     = FUNCTION_SPACE_HGRAD;
      
      const int degreeLength = 1;
      this->fieldOrdinalPolynomialDegree_ = OrdinalTypeArray2DHost("Integrated Legendre H(grad) line polynomial degree lookup", this->basisCardinality_, degreeLength);
      
      for (int i=0; i<this->basisCardinality_; i++)
      {
        // for H(grad) line, if defineVertexFunctions is false, first basis member is constant, second is first-degree, etc.
        // if defineVertexFunctions is true, then the only difference is that the entry is also degree 1
        this->fieldOrdinalPolynomialDegree_(i,0) = i;
      }
      if (defineVertexFunctions)
      {
        this->fieldOrdinalPolynomialDegree_(0,0) = 1;
      }
      
      // initialize tags
      {
        const auto & cardinality = this->basisCardinality_;
        
        // Basis-dependent initializations
        const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
        const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
        const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
        const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
        
        OrdinalTypeArray1DHost tagView("tag view", cardinality*tagSize);
        
        if (defineVertexFunctions) {
          {
            const ordinal_type v0 = 0;
            tagView(v0*tagSize+0) = 0; // vertex dof
            tagView(v0*tagSize+1) = 0; // vertex id
            tagView(v0*tagSize+2) = 0; // local dof id
            tagView(v0*tagSize+3) = 1; // total number of dofs in this vertex
            
            const ordinal_type v1 = 1;
            tagView(v1*tagSize+0) = 0; // vertex dof
            tagView(v1*tagSize+1) = 1; // vertex id
            tagView(v1*tagSize+2) = 0; // local dof id
            tagView(v1*tagSize+3) = 1; // total number of dofs in this vertex
            
            const ordinal_type iend = cardinality - 2;
            for (ordinal_type i=0;i<iend;++i) {
              const auto e = i + 2;
              tagView(e*tagSize+0) = 1;    // edge dof
              tagView(e*tagSize+1) = 0;    // edge id
              tagView(e*tagSize+2) = i;    // local dof id
              tagView(e*tagSize+3) = iend; // total number of dofs in this edge
            }
          }
        } else {
          for (ordinal_type i=0;i<cardinality;++i) {
            tagView(i*tagSize+0) = 1;           // edge dof
            tagView(i*tagSize+1) = 0;           // edge id
            tagView(i*tagSize+2) = i;           // local dof id
            tagView(i*tagSize+3) = cardinality; // total number of dofs in this edge
          }
        }
        
        // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
        // tags are constructed on host
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
    using Basis<ExecutionSpace,OutputScalar,PointScalar>::getValues;
    
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
    virtual void getValues( outputViewType outputValues, const pointViewType  inputPoints,
                           const EOperator operatorType = OPERATOR_VALUE ) const override
    {
      auto numPoints = inputPoints.extent_int(0);
      
      using FunctorType = Hierarchical_HGRAD_LINE_Functor<ExecutionSpace, OutputScalar, PointScalar, outputViewType, pointViewType>;
      
      FunctorType functor(operatorType, outputValues, inputPoints, polyOrder_, defineVertexFunctions);
      
      const int outputVectorSize = getVectorSizeForHierarchicalParallelism<OutputScalar>();
      const int pointVectorSize  = getVectorSizeForHierarchicalParallelism<PointScalar>();
      const int vectorSize = std::max(outputVectorSize,pointVectorSize);
      const int teamSize = 1; // because of the way the basis functions are computed, we don't have a second level of parallelism...
      
      auto policy = Kokkos::TeamPolicy<ExecutionSpace>(numPoints,teamSize,vectorSize);
      Kokkos::parallel_for( policy , functor, "Hierarchical_HGRAD_LINE_Functor");
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_IntegratedLegendreBasis_HGRAD_LINE_h */
