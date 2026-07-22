// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_LegendreBasis_HVOL_PYR.hpp
    \brief  H(vol) basis on the pyramid based on Legendre polynomials.
    \author Created by N.V. Roberts.
 
 Note that although this basis is derived from Legendre polynomials, it is not itself a polynomial basis, but a set of rational functions.
 
 The construction is also hierarchical, in the sense that the basis for p-1 is included in the basis for p.
 
 Intrepid2 has a pre-existing lowest-order HGRAD basis defined on the pyramid, found in Intrepid2_HGRAD_PYR_C1_FEM.hpp; this agrees precisely with this basis when p=1.
 */

#ifndef Intrepid2_LegendreBasis_HVOL_PYR_h
#define Intrepid2_LegendreBasis_HVOL_PYR_h

#include <Kokkos_DynRankView.hpp>

#include <Intrepid2_config.h>

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_Polynomials.hpp"
#include "Intrepid2_PyramidCoords.hpp"
#include "Intrepid2_Utils.hpp"

#include "Teuchos_RCP.hpp"

namespace Intrepid2
{
  /** \class  Intrepid2::Hierarchical_HVOL_PYR_Functor
      \brief  Functor for computing values for the LegendreBasis_HVOL_PYR class.
   
   This functor is not intended for use outside of LegendreBasis_HVOL_PYR.
  */
  template<class DeviceType, class OutputScalar, class PointScalar,
           class OutputFieldType, class InputPointsType>
  struct Hierarchical_HVOL_PYR_Functor
  {
    using ExecutionSpace     = typename DeviceType::execution_space;
    using ScratchSpace        = typename ExecutionSpace::scratch_memory_space;
    using OutputScratchView   = Kokkos::View<OutputScalar*,ScratchSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using OutputScratchView2D = Kokkos::View<OutputScalar**,ScratchSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using PointScratchView    = Kokkos::View<PointScalar*, ScratchSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    
    using TeamPolicy = Kokkos::TeamPolicy<ExecutionSpace>;
    using TeamMember = typename TeamPolicy::member_type;
    
    EOperator opType_;
    
    OutputFieldType  output_;      // F,P
    InputPointsType  inputPoints_; // P,D
    
    int polyOrder_;
    int numFields_, numPoints_;
    
    size_t fad_size_output_;
    
    static const int numVertices   = 5;
    static const int numMixedEdges = 4;
    static const int numTriEdges   = 4;
    static const int numEdges      = 8;
    // the following ordering of the edges matches that used by ESEAS
    // (it *looks* like this is what ESEAS uses; the basis comparison tests will clarify whether I've read their code correctly)
    // see also PyramidEdgeNodeMap in Shards_BasicTopologies.hpp
    const int edge_start_[numEdges] = {0,1,2,3,0,1,2,3}; // edge i is from edge_start_[i] to edge_end_[i]
    const int edge_end_[numEdges]   = {1,2,3,0,4,4,4,4}; // edge i is from edge_start_[i] to edge_end_[i]
    
    // quadrilateral face comes first
    static const int numQuadFaces = 1;
    static const int numTriFaces  = 4;
    
    // face ordering matches ESEAS.  (See BlendProjectPyraTF in ESEAS.)
    const int tri_face_vertex_0[numTriFaces] = {0,1,3,0}; // faces are abc where 0 ≤ a < b < c ≤ 3
    const int tri_face_vertex_1[numTriFaces] = {1,2,2,3};
    const int tri_face_vertex_2[numTriFaces] = {4,4,4,4};
    
    Hierarchical_HVOL_PYR_Functor(EOperator opType, OutputFieldType output, InputPointsType inputPoints,
                                   int polyOrder)
    : opType_(opType), output_(output), inputPoints_(inputPoints),
      polyOrder_(polyOrder),
      fad_size_output_(getScalarDimensionForView(output))
    {
      numFields_ = output.extent_int(0);
      numPoints_ = output.extent_int(1);
      const auto & p = polyOrder;
      const auto p_plus_one_cubed = (p+1) * (p+1) * (p+1);
      INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != inputPoints.extent_int(0), std::invalid_argument, "point counts need to match!");
      INTREPID2_TEST_FOR_EXCEPTION(numFields_ != p_plus_one_cubed, std::invalid_argument, "output field size does not match basis cardinality");
    }
    
    KOKKOS_INLINE_FUNCTION
    void operator()( const TeamMember & teamMember ) const
    {
      auto pointOrdinal = teamMember.league_rank();
      OutputScratchView scratch1D_1, scratch1D_2, scratch1D_3;
      OutputScratchView scratch1D_4, scratch1D_5, scratch1D_6;
      OutputScratchView scratch1D_7, scratch1D_8, scratch1D_9;
      OutputScratchView2D scratch2D_1, scratch2D_2, scratch2D_3;
      const int numAlphaValues = (polyOrder_-1 > 1) ? (polyOrder_-1) : 1; // make numAlphaValues at least 1 so we can avoid zero-extent allocations…
      if (fad_size_output_ > 0) {
        scratch1D_1 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1, fad_size_output_);
        scratch1D_2 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1, fad_size_output_);
        scratch1D_3 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1, fad_size_output_);
        scratch1D_4 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1, fad_size_output_);
        scratch1D_5 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1, fad_size_output_);
        scratch1D_6 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1, fad_size_output_);
        scratch1D_7 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1, fad_size_output_);
        scratch1D_8 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1, fad_size_output_);
        scratch1D_9 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1, fad_size_output_);
        scratch2D_1 = OutputScratchView2D(teamMember.team_shmem(), numAlphaValues, polyOrder_ + 1, fad_size_output_);
        scratch2D_2 = OutputScratchView2D(teamMember.team_shmem(), numAlphaValues, polyOrder_ + 1, fad_size_output_);
        scratch2D_3 = OutputScratchView2D(teamMember.team_shmem(), numAlphaValues, polyOrder_ + 1, fad_size_output_);
      }
      else {
        scratch1D_1 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1);
        scratch1D_2 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1);
        scratch1D_3 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1);
        scratch1D_4 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1);
        scratch1D_5 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1);
        scratch1D_6 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1);
        scratch1D_7 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1);
        scratch1D_8 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1);
        scratch1D_9 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1);
        scratch2D_1 = OutputScratchView2D(teamMember.team_shmem(), numAlphaValues, polyOrder_ + 1);
        scratch2D_2 = OutputScratchView2D(teamMember.team_shmem(), numAlphaValues, polyOrder_ + 1);
        scratch2D_3 = OutputScratchView2D(teamMember.team_shmem(), numAlphaValues, polyOrder_ + 1);
      }
      
      const auto & x = inputPoints_(pointOrdinal,0);
      const auto & y = inputPoints_(pointOrdinal,1);
      const auto & z = inputPoints_(pointOrdinal,2);
      
      // Intrepid2 uses (-1,1)^2 for x,y
      // ESEAS uses (0,1)^2
      // (Can look at what we do on the HGRAD_LINE for reference; there's a similar difference for line topology.)
      
      Kokkos::Array<PointScalar,3> coords;
      transformToESEASPyramid<>(coords[0], coords[1], coords[2], x, y, z); // map x,y coordinates from (-z,z)^2 to (0,z)^2
      
      // pyramid "affine" coordinates and gradients get stored in lambda, lambdaGrad:
      using Kokkos::Array;
      Array<PointScalar,5> lambda;
      Array<Kokkos::Array<PointScalar,3>,5> lambdaGrad;
      
      Array<Array<PointScalar,3>,2> mu; // first index is subscript; second is superscript: 0 --> (\zeta, xi_1), 1 --> (\zeta, xi_2), 2 --> (\zeta)
      Array<Array<Array<PointScalar,3>,3>,2> muGrad;
      
      Array<Array<PointScalar,2>,3> nu;
      Array<Array<Array<PointScalar,3>,2>,3> nuGrad;
      
      affinePyramid(lambda, lambdaGrad, mu, muGrad, nu, nuGrad, coords);
      
      switch (opType_)
      {
        case OPERATOR_VALUE:
        {
          // interior functions
          // rename scratch
          ordinal_type fieldOrdinalOffset = 0;
          auto & Pi = scratch1D_1;
          auto & Pj = scratch1D_2;
          auto & Pk = scratch1D_3;
          // [P_i](mu_01^{\zeta,\xi_1})
          Polynomials::shiftedScaledLegendreValues(Pi, polyOrder_, mu[1][0], mu[0][0] + mu[1][0]);
          // [P_j](mu_01^{\zeta,\xi_2})
          Polynomials::shiftedScaledLegendreValues(Pj, polyOrder_, mu[1][1], mu[0][1] + mu[1][1]);
          // [P_k](mu_01^{\zeta})
          Polynomials::shiftedScaledLegendreValues(Pk, polyOrder_, mu[1][2], mu[0][2] + mu[1][2]);
          // (grad nu_1^{\zeta,\xi_1} x grad nu_1^{\zeta,\xi_2}) \cdot grad mu_1^zeta:
          PointScalar grad_weight =
            (nuGrad[1][0][1] * nuGrad[1][1][2] - nuGrad[1][0][2] * nuGrad[1][1][1]) * muGrad[1][2][0]
          + (nuGrad[1][0][2] * nuGrad[1][1][0] - nuGrad[1][0][0] * nuGrad[1][1][2]) * muGrad[1][2][1]
          + (nuGrad[1][0][0] * nuGrad[1][1][1] - nuGrad[1][0][1] * nuGrad[1][1][0]) * muGrad[1][2][2];
          
          // following the ESEAS ordering: k increments first
          for (int k=0; k<=polyOrder_; k++)
          {
            for (int j=0; j<=polyOrder_; j++)
            {
              for (int i=0; i<=polyOrder_; i++)
              {
                output_(fieldOrdinalOffset,pointOrdinal) = Pk(k) * Pi(i) * Pj(j) * grad_weight;
                fieldOrdinalOffset++;
              }
            }
          }
        } // end OPERATOR_VALUE
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
          INTREPID2_TEST_FOR_ABORT( true,
                                   ">>> ERROR: (Intrepid2::Hierarchical_HVOL_PYR_Functor) Computing of derivatives is not supported");
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
      // we use shared memory to create a fast buffer for basis computations
      // for the (integrated) Legendre computations, we just need p+1 values stored.  For interior functions on the pyramid, we have up to 3 scratch arrays with (integrated) Legendre values stored, for each of the 3 directions (i,j,k indices): a total of 9.
      // for the (integrated) Jacobi computations, though, we want (p+1)*(# alpha values)
      // alpha is either 2i or 2(i+j), where i=2,…,p or i+j=3,…,p.  So there are at most (p-1) alpha values needed.
      // We can have up to 3 of the (integrated) Jacobi values needed at once.
      const int numAlphaValues = std::max(polyOrder_-1, 1); // make it at least 1 so we can avoid zero-extent ranks…
      size_t shmem_size = 0;
      if (fad_size_output_ > 0)
      {
        // Legendre:
        shmem_size += 9 * OutputScratchView::shmem_size(polyOrder_ + 1, fad_size_output_);
        // Jacobi:
        shmem_size += 3 * OutputScratchView2D::shmem_size(numAlphaValues, polyOrder_ + 1, fad_size_output_);
      }
      else
      {
        // Legendre:
        shmem_size += 9 * OutputScratchView::shmem_size(polyOrder_ + 1);
        // Jacobi:
        shmem_size += 3 * OutputScratchView2D::shmem_size(numAlphaValues, polyOrder_ + 1);
      }
      
      return shmem_size;
    }
  };
  
  /** \class  Intrepid2::LegendreBasis_HVOL_PYR
      \brief  Basis defining integrated Legendre basis on the line, a polynomial subspace of H(grad) on the line.

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
  class LegendreBasis_HVOL_PYR
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
     
     The basis will have polyOrder^3 + 3 * polyOrder + 1 members.
     */
    LegendreBasis_HVOL_PYR(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    polyOrder_(polyOrder),
    pointType_(pointType)
    {
      INTREPID2_TEST_FOR_EXCEPTION(pointType!=POINTTYPE_DEFAULT,std::invalid_argument,"PointType not supported");
      this->basisCardinality_     = (polyOrder + 1) * (polyOrder + 1) * (polyOrder + 1);
      this->basisDegree_          = polyOrder;
      this->basisCellTopologyKey_ = shards::Pyramid<>::key;
      this->basisType_            = BASIS_FEM_HIERARCHICAL;
      this->basisCoordinates_     = COORDINATES_CARTESIAN;
      this->functionSpace_        = FUNCTION_SPACE_HVOL;
      
      const int degreeLength = 1;
      this->fieldOrdinalPolynomialDegree_ = OrdinalTypeArray2DHost("Legendre H(vol) pyramid polynomial degree lookup", this->basisCardinality_, degreeLength);
      this->fieldOrdinalH1PolynomialDegree_ = OrdinalTypeArray2DHost("Legendre H(vol) pyramid polynomial H^1 degree lookup", this->basisCardinality_, degreeLength);
      
      int fieldOrdinalOffset = 0;
      
      // **** interior functions **** //
      const int numVolumes = 1; // interior
      for (int volumeOrdinal=0; volumeOrdinal<numVolumes; volumeOrdinal++)
      {
        // following the ESEAS ordering: k increments first
        for (int k=0; k<=polyOrder_; k++)
        {
          for (int j=0; j<=polyOrder_; j++)
          {
            for (int i=0; i<=polyOrder_; i++)
            {
              const int max_ij  = std::max(i,j);
              const int max_ijk = std::max(max_ij,k);
              this->fieldOrdinalPolynomialDegree_  (fieldOrdinalOffset,0) = max_ijk;     // L^2 degree
              this->fieldOrdinalH1PolynomialDegree_(fieldOrdinalOffset,0) = max_ijk + 1; // H^1 degree
              fieldOrdinalOffset++;
            }
          }
        }
      }
      
      INTREPID2_TEST_FOR_EXCEPTION(fieldOrdinalOffset != this->basisCardinality_, std::invalid_argument, "Internal error: basis enumeration is incorrect");
      
      // initialize tags
      {
        const auto & cardinality = this->basisCardinality_;
        
        // Basis-dependent initializations
        const ordinal_type tagSize  = 4; // size of DoF tag, i.e., number of fields in the tag
        const ordinal_type posScDim = 0; // position in the tag, counting from 0, of the subcell dim
        const ordinal_type posScOrd = 1; // position in the tag, counting from 0, of the subcell ordinal
        const ordinal_type posDfOrd = 2; // position in the tag, counting from 0, of DoF ordinal relative to the subcell
        
        OrdinalTypeArray1DHost tagView("tag view", cardinality*tagSize);
        const ordinal_type volumeDim = 3;

        for (ordinal_type i=0;i<cardinality;++i) {
          tagView(i*tagSize+0) = volumeDim;   // volume dimension
          tagView(i*tagSize+1) = 0;           // volume ordinal
          tagView(i*tagSize+2) = i;           // local dof id
          tagView(i*tagSize+3) = cardinality; // total number of dofs on this volume
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
      return "Intrepid2_LegendreBasis_HVOL_PYR";
    }
    
    /** \brief True if orientation is required
    */
    virtual bool requireOrientation() const override {
      return false;
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
    virtual void getValues( OutputViewType outputValues, const PointViewType  inputPoints,
                           const EOperator operatorType = OPERATOR_VALUE ) const override
    {
      auto numPoints = inputPoints.extent_int(0);
      
      using FunctorType = Hierarchical_HVOL_PYR_Functor<DeviceType, OutputScalar, PointScalar, OutputViewType, PointViewType>;
      
      FunctorType functor(operatorType, outputValues, inputPoints, polyOrder_);
      
      const int outputVectorSize = getVectorSizeForHierarchicalParallelism<OutputScalar>();
      const int pointVectorSize  = getVectorSizeForHierarchicalParallelism<PointScalar>();
      const int vectorSize = std::max(outputVectorSize,pointVectorSize);
      const int teamSize = 1; // because of the way the basis functions are computed, we don't have a second level of parallelism...

      auto policy = Kokkos::TeamPolicy<ExecutionSpace>(numPoints,teamSize,vectorSize);
      Kokkos::parallel_for("Hierarchical_HVOL_PYR_Functor", policy, functor);
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
      // no subcell ref basis for HVOL
      INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Input parameters out of bounds");
    }

    /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.
     
        \return Pointer to the new Basis object.
     */
    virtual BasisPtr<typename Kokkos::HostSpace::device_type, OutputScalar, PointScalar>
    getHostBasis() const override {
      using HostDeviceType = typename Kokkos::HostSpace::device_type;
      using HostBasisType  = LegendreBasis_HVOL_PYR<HostDeviceType, OutputScalar, PointScalar>;
      return Teuchos::rcp( new HostBasisType(polyOrder_, pointType_) );
    }
  };
} // end namespace Intrepid2

// do ETI with default (double) type
extern template class Intrepid2::LegendreBasis_HVOL_PYR<Kokkos::DefaultExecutionSpace::device_type,double,double>;

#endif /* Intrepid2_LegendreBasis_HVOL_PYR_h */
