// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_LegendreBasis_HVOL_LINE.hpp
    \brief  H(vol) basis on the line based on Legendre polynomials.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_LegendreBasis_HVOL_LINE_h
#define Intrepid2_LegendreBasis_HVOL_LINE_h

#include <Kokkos_DynRankView.hpp>

#include <Intrepid2_config.h>

// Sacado header that defines some fad sizing…
#ifdef HAVE_INTREPID2_SACADO
#include <KokkosExp_View_Fad.hpp>
#endif

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_DeviceAssert.hpp"
#include "Intrepid2_Polynomials.hpp"

namespace Intrepid2
{
  /** \class  Intrepid2::Hierarchical_HVOL_LINE_Functor
      \brief  Functor for computing values for the LegendreBasis_HVOL_LINE class.
   
   This functor is not intended for use outside of LegendreBasis_HVOL_LINE.
  */
  template<class DeviceType, class OutputScalar, class PointScalar,
  class OutputFieldType, class InputPointsType>
  struct Hierarchical_HVOL_LINE_Functor
  {
    using ExecutionSpace     = typename DeviceType::execution_space;
    using ScratchSpace       = typename ExecutionSpace::scratch_memory_space;
    using OutputScratchView  = Kokkos::View<OutputScalar*,ScratchSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using PointScratchView   = Kokkos::View<PointScalar*, ScratchSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    
    using TeamPolicy = Kokkos::TeamPolicy<ExecutionSpace>;
    using TeamMember = typename TeamPolicy::member_type;
    
    OutputFieldType  output_;      // F,P
    InputPointsType  inputPoints_; // P,D
    
    int polyOrder_;
    int numFields_, numPoints_;
    
    EOperator op_;
    
    size_t fad_size_output_;
    
    Hierarchical_HVOL_LINE_Functor(OutputFieldType output, InputPointsType inputPoints,
                                   int polyOrder, EOperator op)
    : output_(output), inputPoints_(inputPoints), polyOrder_(polyOrder), op_(op),
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
      
      const PointScalar x = inputPoints_(pointOrdinal,0);

      switch (op_)
      {
        case OPERATOR_VALUE:
          Polynomials::legendreValues(field_values_at_point, polyOrder_, x);
          
          // note that because legendreValues determines field values recursively, there is not much
          // opportunity at that level for further parallelism
          break;
        case OPERATOR_GRAD:
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
          auto derivativeOrder = getOperatorOrder(op_);
          Polynomials::legendreDerivativeValues(field_values_at_point, polyOrder_, x, derivativeOrder);
          break;
        }
        default:
          // unsupported operator type
          device_assert(false);
      }
      // copy the values into the output container
      for (int fieldOrdinal=0; fieldOrdinal<numFields_; fieldOrdinal++)
      {
        output_.access(fieldOrdinal,pointOrdinal,0) = field_values_at_point(fieldOrdinal);
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
  
  /** \class  Intrepid2::LegendreBasis_HVOL_LINE
      \brief  Basis defining Legendre basis on the line, a polynomial subspace of L^2 (a.k.a. H(vol)) on the line.

              This is used in the construction of hierarchical bases on higher-dimensional topologies.  For
              mathematical details of the construction, see:
   
               Federico Fuentes, Brendan Keith, Leszek Demkowicz, Sriram Nagaraj.
               "Orientation embedded high order shape functions for the exact sequence elements of all shapes."
               Computers & Mathematics with Applications, Volume 70, Issue 4, 2015, Pages 353-458, ISSN 0898-1221.
               https://doi.org/10.1016/j.camwa.2015.04.027.
  */
  template<typename DeviceType,
           typename OutputScalar = double,
           typename PointScalar  = double>
  class LegendreBasis_HVOL_LINE
  : public Basis<DeviceType,OutputScalar,PointScalar>
  {
  public:
    using BasisBase = Basis<DeviceType,OutputScalar,PointScalar>;
    using HostBasis = LegendreBasis_HVOL_LINE<typename Kokkos::HostSpace::device_type,OutputScalar,PointScalar>;

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
     
     polyOrder defines the maximum polynomial order in the basis; note that this means that to have an exact sequence
     between this and IntegratedLegendreBasis_HGRAD_LINE, the latter should be constructed with (polyOrder+1) as its
     polynomial order.
     
     */
    LegendreBasis_HVOL_LINE(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    polyOrder_(polyOrder),
    pointType_(pointType)
    {
      INTREPID2_TEST_FOR_EXCEPTION(pointType!=POINTTYPE_DEFAULT,std::invalid_argument,"PointType not supported");
      this->basisCardinality_     = polyOrder+1;
      this->basisDegree_          = polyOrder;
      this->basisCellTopologyKey_ = shards::Line<2>::key;
      this->basisType_            = BASIS_FEM_HIERARCHICAL;
      this->basisCoordinates_     = COORDINATES_CARTESIAN;
      this->functionSpace_        = FUNCTION_SPACE_HVOL;
      
      const int degreeLength = 1;
      this->fieldOrdinalPolynomialDegree_ = OrdinalTypeArray2DHost("Integrated Legendre H(grad) line polynomial degree lookup", this->basisCardinality_, degreeLength);
      this->fieldOrdinalH1PolynomialDegree_ = OrdinalTypeArray2DHost("Integrated Legendre H(grad) line polynomial degree lookup", this->basisCardinality_, degreeLength);
      
      for (int i=0; i<this->basisCardinality_; i++)
      {
        // for H(vol) line, first basis member is constant, second is first-degree, etc.
        this->fieldOrdinalPolynomialDegree_  (i,0) = i;
        this->fieldOrdinalH1PolynomialDegree_(i,0) = i+1; // H^1 degree is one greater
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
        
        for (ordinal_type i=0;i<cardinality;++i) {
          tagView(i*tagSize+0) = 1;           // edge dof
          tagView(i*tagSize+1) = 0;           // edge id
          tagView(i*tagSize+2) = i;           // local dof id
          tagView(i*tagSize+3) = cardinality; // total number of dofs in this edge
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
    using BasisBase::getValues;
    
    /** \brief  Returns basis name
     
     \return the name of the basis
     */
    virtual
    const char*
    getName() const override {
      return "Intrepid2_LegendreBasis_HVOL_LINE";
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
    virtual void getValues( OutputViewType outputValues, const PointViewType  inputPoints,
                           const EOperator operatorType = OPERATOR_VALUE ) const override
    {
      auto numPoints = inputPoints.extent_int(0);
      
      using FunctorType = Hierarchical_HVOL_LINE_Functor<DeviceType, OutputScalar, PointScalar, OutputViewType, PointViewType>;
      
      FunctorType functor(outputValues, inputPoints, polyOrder_, operatorType);
      
      const int outputVectorSize = getVectorSizeForHierarchicalParallelism<OutputScalar>();
      const int pointVectorSize  = getVectorSizeForHierarchicalParallelism<PointScalar>();
      const int vectorSize = std::max(outputVectorSize,pointVectorSize);
      const int teamSize = 1; // because of the way the basis functions are computed, we don't have a second level of parallelism...

      auto policy =
          Kokkos::TeamPolicy<ExecutionSpace>(numPoints, teamSize, vectorSize);
      Kokkos::parallel_for("Hierarchical_HVOL_LINE_Functor", policy , functor);
    }
    
    /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.
     
        \return Pointer to the new Basis object.
     */
    virtual BasisPtr<typename Kokkos::HostSpace::device_type, OutputScalar, PointScalar>
    getHostBasis() const override {
      using HostDeviceType = typename Kokkos::HostSpace::device_type;
      using HostBasisType  = LegendreBasis_HVOL_LINE<HostDeviceType, OutputScalar, PointScalar>;
      return Teuchos::rcp( new HostBasisType(polyOrder_, pointType_) );
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_LegendreBasis_HVOL_LINE_h */
