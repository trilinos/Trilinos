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

/** \file   Intrepid2_LegendreBasis_HVOL_TRI.hpp
    \brief  H(vol) basis on the triangle based on integrated Legendre polynomials.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_LegendreBasis_HVOL_TRI_h
#define Intrepid2_LegendreBasis_HVOL_TRI_h

#include <Kokkos_View.hpp>
#include <Kokkos_DynRankView.hpp>

#include <Intrepid2_config.h>

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_IntegratedLegendreBasis_HGRAD_LINE.hpp"
#include "Intrepid2_Polynomials.hpp"
#include "Intrepid2_Utils.hpp"

namespace Intrepid2
{
  /** \class  Intrepid2::Hierarchical_HVOL_TRI_Functor
      \brief  Functor for computing values for the LegendreBasis_HVOL_TRI class.
   
   This functor is not intended for use outside of LegendreBasis_HVOL_TRI.
  */
  template<class DeviceType, class OutputScalar, class PointScalar,
           class OutputFieldType, class InputPointsType>
  struct Hierarchical_HVOL_TRI_Functor
  {
    using ExecutionSpace     = typename DeviceType::execution_space;
    using ScratchSpace       = typename ExecutionSpace::scratch_memory_space;
    using OutputScratchView  = Kokkos::View<OutputScalar*,ScratchSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using PointScratchView   = Kokkos::View<PointScalar*, ScratchSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    
    using TeamPolicy = Kokkos::TeamPolicy<ExecutionSpace>;
    using TeamMember = typename TeamPolicy::member_type;
    
    EOperator opType_;
    
    OutputFieldType  output_;      // F,P
    InputPointsType  inputPoints_; // P,D
    
    int polyOrder_;
    int numFields_, numPoints_;
    
    size_t fad_size_output_;
    
    static const int numVertices = 3;
    static const int numEdges    = 3;
    const int edge_start_[numEdges] = {0,1,0}; // edge i is from edge_start_[i] to edge_end_[i]
    const int edge_end_[numEdges]   = {1,2,2}; // edge i is from edge_start_[i] to edge_end_[i]
    
    Hierarchical_HVOL_TRI_Functor(EOperator opType, OutputFieldType output, InputPointsType inputPoints, int polyOrder)
    : opType_(opType), output_(output), inputPoints_(inputPoints),
      polyOrder_(polyOrder),
      fad_size_output_(getScalarDimensionForView(output))
    {
      numFields_ = output.extent_int(0);
      numPoints_ = output.extent_int(1);
      INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != inputPoints.extent_int(0), std::invalid_argument, "point counts need to match!");
      INTREPID2_TEST_FOR_EXCEPTION(numFields_ != (polyOrder_+1)*(polyOrder_+2)/2, std::invalid_argument, "output field size does not match basis cardinality");
    }
    
    KOKKOS_INLINE_FUNCTION
    void operator()( const TeamMember & teamMember ) const
    {
      auto pointOrdinal = teamMember.league_rank();
      OutputScratchView legendre_field_values_at_point, jacobi_values_at_point;
      if (fad_size_output_ > 0) {
        legendre_field_values_at_point = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1, fad_size_output_);
        jacobi_values_at_point         = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1, fad_size_output_);
      }
      else {
        legendre_field_values_at_point = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1);
        jacobi_values_at_point         = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1);
      }
      
      const auto & x = inputPoints_(pointOrdinal,0);
      const auto & y = inputPoints_(pointOrdinal,1);
      
      // write as barycentric coordinates:
      const PointScalar lambda[3]    = {1. - x - y, x, y};
      
      switch (opType_)
      {
        case OPERATOR_VALUE:
        {
          // face functions
          {
            const PointScalar tLegendre = lambda[0] + lambda[1];
            Polynomials::shiftedScaledLegendreValues(legendre_field_values_at_point, polyOrder_, lambda[1], tLegendre);

            int fieldOrdinalOffset = 0;
            const int max_ij_sum = polyOrder_;
            for (int ij_sum=0; ij_sum<=max_ij_sum; ij_sum++)
            {
              for (int i=0; i<=ij_sum; i++)
              {
                const int j = ij_sum - i;
                const auto & legendreValue = legendre_field_values_at_point(i);
                const double alpha = i*2.0+1;
                
                const PointScalar tJacobi = 1.0;// lambda[0] + lambda[1] + lambda[2];
                Polynomials::shiftedScaledJacobiValues(jacobi_values_at_point, alpha, polyOrder_, lambda[2], tJacobi);
                
                const auto & jacobiValue = jacobi_values_at_point(j);
                output_(fieldOrdinalOffset,pointOrdinal) = legendreValue * jacobiValue;
                fieldOrdinalOffset++;
              }
            }
          }
        } // end OPERATOR_VALUE
          break;
        default:
          // unsupported operator type
          device_assert(false);
      }
    }
    
    // Provide the shared memory capacity.
    // This function takes the team_size as an argument,
    // which allows team_size-dependent allocations.
    size_t team_shmem_size (int team_size) const
    {
      // TODO: edit this to match scratch that we actually need.  (What's here is copied from H^1 basis on triangles…)
      // we will use shared memory to create a fast buffer for basis computations
      size_t shmem_size = 0;
      if (fad_size_output_ > 0)
        shmem_size += 2 * OutputScratchView::shmem_size(polyOrder_ + 1, fad_size_output_);
      else
        shmem_size += 2 * OutputScratchView::shmem_size(polyOrder_ + 1);
      
      return shmem_size;
    }
  };
  
  /** \class  Intrepid2::LegendreBasis_HVOL_TRI
      \brief  Basis defining Legendre basis on the line, a polynomial subspace of H(vol) on the line: extension to triangle using Jacobi blending function.

              For mathematical details of the construction, see:
   
               Federico Fuentes, Brendan Keith, Leszek Demkowicz, Sriram Nagaraj.
               "Orientation embedded high order shape functions for the exact sequence elements of all shapes."
               Computers & Mathematics with Applications, Volume 70, Issue 4, 2015, Pages 353-458, ISSN 0898-1221.
               https://doi.org/10.1016/j.camwa.2015.04.027.
  */
  template<typename DeviceType,
           typename OutputScalar = double,
           typename PointScalar  = double>
  class LegendreBasis_HVOL_TRI
  : public Basis<DeviceType,OutputScalar,PointScalar>
  {
  public:
    using BasisBase = Basis<DeviceType,OutputScalar,PointScalar>;

    using typename BasisBase::OrdinalTypeArray1DHost;
    using typename BasisBase::OrdinalTypeArray2DHost;

    using typename BasisBase::OutputViewType;
    using typename BasisBase::PointViewType;
    using typename BasisBase::ScalarViewType;

    using typename BasisBase::ExecutionSpace;

  protected:
    int polyOrder_; // the maximum order of the polynomial
    EPointType pointType_;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order of the basis.
     
     The basis will have (polyOrder + 1)*(polyOrder + 2) / 2 members, and is in a discrete exact sequence that begins with the integrated Legendre basis of order polyOrder + 1.
     
     */
    LegendreBasis_HVOL_TRI(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    polyOrder_(polyOrder),
    pointType_(pointType)
    {
      INTREPID2_TEST_FOR_EXCEPTION(pointType!=POINTTYPE_DEFAULT,std::invalid_argument,"PointType not supported");

      this->basisCardinality_  = ((polyOrder+2) * (polyOrder+1)) / 2;
      this->basisDegree_       = polyOrder;
      this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<> >() );
      this->basisType_         = BASIS_FEM_HIERARCHICAL;
      this->basisCoordinates_  = COORDINATES_CARTESIAN;
      this->functionSpace_     = FUNCTION_SPACE_HVOL;
      
      const int degreeLength = 1;
      this->fieldOrdinalPolynomialDegree_ = OrdinalTypeArray2DHost("Integrated Legendre H(vol) triangle polynomial degree lookup", this->basisCardinality_, degreeLength);
      
      int fieldOrdinalOffset = 0;
      // **** face functions **** //
      const int max_ij_sum = polyOrder;
      for (int ij_sum=0; ij_sum<=max_ij_sum; ij_sum++)
      {
        for (int i=0; i<=ij_sum; i++)
        {
          const int j = ij_sum - i;
          this->fieldOrdinalPolynomialDegree_(fieldOrdinalOffset,0) = i+j;
          fieldOrdinalOffset++;
        }
      }
      INTREPID2_TEST_FOR_EXCEPTION(fieldOrdinalOffset != this->basisCardinality_, std::invalid_argument, "Internal error: basis enumeration is incorrect");
      
      // initialize tags
      {
        const auto & cardinality = this->basisCardinality_;
        
        // Basis-dependent initializations
        const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
        const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
        const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
        const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
        
        OrdinalTypeArray1DHost tagView("tag view", cardinality*tagSize);
        const int faceDim = 2;

        for (ordinal_type i=0;i<cardinality;++i) {
          tagView(i*tagSize+0) = faceDim;     // face dimension
          tagView(i*tagSize+1) = 0;           // face id
          tagView(i*tagSize+2) = i;           // local dof id
          tagView(i*tagSize+3) = cardinality; // total number of dofs on this face
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
    
    /** \brief  Returns basis name
     
     \return the name of the basis
     */
    const char* getName() const override {
      return "Intrepid2_LegendreBasis_HVOL_TRI";
    }
    
    /** \brief True if orientation is required
    */
    virtual bool requireOrientation() const override {
      return (this->getDegree() > 2);
    }

    // since the getValues() below only overrides the FEM variant, we specify that
    // we use the base class's getValues(), which implements the FVD variant by throwing an exception.
    // (It's an error to use the FVD variant on this basis.)
    using BasisBase::getValues;
    
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
    virtual void getValues( OutputViewType outputValues, const PointViewType inputPoints,
                           const EOperator operatorType = OPERATOR_VALUE ) const override
    {
      auto numPoints = inputPoints.extent_int(0);
      
      using FunctorType = Hierarchical_HVOL_TRI_Functor<DeviceType, OutputScalar, PointScalar, OutputViewType, PointViewType>;
      
      FunctorType functor(operatorType, outputValues, inputPoints, polyOrder_);
      
      const int outputVectorSize = getVectorSizeForHierarchicalParallelism<OutputScalar>();
      const int pointVectorSize  = getVectorSizeForHierarchicalParallelism<PointScalar>();
      const int vectorSize = std::max(outputVectorSize,pointVectorSize);
      const int teamSize = 1; // because of the way the basis functions are computed, we don't have a second level of parallelism...

      auto policy = Kokkos::TeamPolicy<ExecutionSpace>(numPoints,teamSize,vectorSize);
      Kokkos::parallel_for( policy , functor, "Hierarchical_HVOL_TRI_Functor");
    }

    /** \brief returns the basis associated to a subCell.

        The bases of the subCell are the restriction to the subCell
        of the bases of the parent cell.
        \param [in] subCellDim - dimension of subCell
        \param [in] subCellOrd - position of the subCell among of the subCells having the same dimension
        \return pointer to the subCell basis of dimension subCellDim and position subCellOrd
     */
    BasisPtr<DeviceType,OutputScalar,PointScalar>
      getSubCellRefBasis(const ordinal_type subCellDim, const ordinal_type subCellOrd) const override{
      if(subCellDim == 1) {
        return Teuchos::rcp(new
            IntegratedLegendreBasis_HGRAD_LINE<DeviceType,OutputScalar,PointScalar>
                    (this->basisDegree_));
      }
      INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Input parameters out of bounds");
    }

    /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.
     
        \return Pointer to the new Basis object.
     */
    virtual BasisPtr<typename Kokkos::HostSpace::device_type, OutputScalar, PointScalar>
    getHostBasis() const override {
      using HostDeviceType = typename Kokkos::HostSpace::device_type;
      using HostBasisType  = LegendreBasis_HVOL_TRI<HostDeviceType, OutputScalar, PointScalar>;
      return Teuchos::rcp( new HostBasisType(polyOrder_, pointType_) );
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_LegendreBasis_HVOL_TRI_h */
