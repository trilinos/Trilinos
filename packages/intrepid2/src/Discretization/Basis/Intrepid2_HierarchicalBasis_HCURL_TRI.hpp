// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HierarchicalBasis_HCURL_TRI.hpp
    \brief  H(curl) basis on the triangle using a construction involving Legendre and integrated Jacobi polynomials.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_HierarchicalBasis_HCURL_TRI_h
#define Intrepid2_HierarchicalBasis_HCURL_TRI_h

#include <Kokkos_DynRankView.hpp>

#include <Intrepid2_config.h>

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_IntegratedLegendreBasis_HGRAD_LINE.hpp"
#include "Intrepid2_Polynomials.hpp"
#include "Intrepid2_Utils.hpp"

namespace Intrepid2
{
  /** \class  Intrepid2::Hierarchical_HCURL_TRI_Functor
      \brief  Functor for computing values for the HierarchicalBasis_HCURL_TRI class.
   
   This functor is not intended for use outside of HierarchicalBasis_HCURL_TRI.
  */
  template<class DeviceType, class OutputScalar, class PointScalar,
           class OutputFieldType, class InputPointsType>
  struct Hierarchical_HCURL_TRI_Functor
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
    
    static const int numVertices     = 3;
    static const int numEdges        = 3;
    static const int numFaceFamilies = 2;
    const int edge_start_[numEdges] = {0,1,0}; // edge i is from edge_start_[i] to edge_end_[i]
    const int edge_end_[numEdges]   = {1,2,2}; // edge i is from edge_start_[i] to edge_end_[i]
    const int face_family_start_[numFaceFamilies]  = {0,1};
    const int face_family_middle_[numFaceFamilies] = {1,2};
    const int face_family_end_   [numFaceFamilies] = {2,0};
    
    Hierarchical_HCURL_TRI_Functor(EOperator opType, OutputFieldType output, InputPointsType inputPoints, int polyOrder)
    : opType_(opType), output_(output), inputPoints_(inputPoints),
      polyOrder_(polyOrder),
      fad_size_output_(getScalarDimensionForView(output))
    {
      numFields_ = output.extent_int(0);
      numPoints_ = output.extent_int(1);
      const int expectedCardinality = 3 * polyOrder_ + polyOrder_ * (polyOrder-1);
      
      INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != inputPoints.extent_int(0), std::invalid_argument, "point counts need to match!");
      INTREPID2_TEST_FOR_EXCEPTION(numFields_ != expectedCardinality, std::invalid_argument, "output field size does not match basis cardinality");
    }
    
    KOKKOS_INLINE_FUNCTION
    void operator()( const TeamMember & teamMember ) const
    {
      auto pointOrdinal = teamMember.league_rank();
      OutputScratchView edge_field_values_at_point, jacobi_values_at_point, other_values_at_point, other_values2_at_point;
      if (fad_size_output_ > 0) {
        edge_field_values_at_point = OutputScratchView(teamMember.team_shmem(), polyOrder_, fad_size_output_);
        jacobi_values_at_point     = OutputScratchView(teamMember.team_shmem(), polyOrder_, fad_size_output_);
        other_values_at_point      = OutputScratchView(teamMember.team_shmem(), polyOrder_, fad_size_output_);
        other_values2_at_point     = OutputScratchView(teamMember.team_shmem(), polyOrder_, fad_size_output_);
      }
      else {
        edge_field_values_at_point = OutputScratchView(teamMember.team_shmem(), polyOrder_);
        jacobi_values_at_point     = OutputScratchView(teamMember.team_shmem(), polyOrder_);
        other_values_at_point      = OutputScratchView(teamMember.team_shmem(), polyOrder_);
        other_values2_at_point     = OutputScratchView(teamMember.team_shmem(), polyOrder_);
      }
      
      const auto & x = inputPoints_(pointOrdinal,0);
      const auto & y = inputPoints_(pointOrdinal,1);
      
      // write as barycentric coordinates:
      const PointScalar lambda[3]    = {1. - x - y, x, y};
      const PointScalar lambda_dx[3] = {-1., 1., 0.};
      const PointScalar lambda_dy[3] = {-1., 0., 1.};
      
      const int num1DEdgeFunctions = polyOrder_; // per edge
      
      switch (opType_)
      {
        case OPERATOR_VALUE:
        {
          // edge functions
          int fieldOrdinalOffset = 0;
          for (int edgeOrdinal=0; edgeOrdinal<numEdges; edgeOrdinal++)
          {
            const auto & s0    = lambda   [edge_start_[edgeOrdinal]];
            const auto & s0_dx = lambda_dx[edge_start_[edgeOrdinal]];
            const auto & s0_dy = lambda_dy[edge_start_[edgeOrdinal]];
            
            const auto & s1    = lambda   [  edge_end_[edgeOrdinal]];
            const auto & s1_dx = lambda_dx[  edge_end_[edgeOrdinal]];
            const auto & s1_dy = lambda_dy[  edge_end_[edgeOrdinal]];
            
            Polynomials::shiftedScaledLegendreValues(edge_field_values_at_point, polyOrder_-1, PointScalar(s1), PointScalar(s0+s1));
            for (int edgeFunctionOrdinal=0; edgeFunctionOrdinal<num1DEdgeFunctions; edgeFunctionOrdinal++)
            {
              const auto & legendreValue = edge_field_values_at_point(edgeFunctionOrdinal);
              const PointScalar xWeight = s0 * s1_dx - s1 * s0_dx;
              const PointScalar yWeight = s0 * s1_dy - s1 * s0_dy;
              output_(edgeFunctionOrdinal+fieldOrdinalOffset,pointOrdinal,0) = legendreValue * xWeight;
              output_(edgeFunctionOrdinal+fieldOrdinalOffset,pointOrdinal,1) = legendreValue * yWeight;
            }
            fieldOrdinalOffset += num1DEdgeFunctions;
          }
          
          // face functions
          {
            // these functions multiply the edge functions from the 01 edge by integrated Jacobi functions, appropriately scaled
            const double jacobiScaling = 1.0; // s0 + s1 + s2
            
            const int max_ij_sum = polyOrder_ - 1;
            
            // following ESEAS, we interleave the face families.  This groups all the face dofs of a given degree together.
            const int faceFieldOrdinalOffset = fieldOrdinalOffset;
            for (int familyOrdinal=1; familyOrdinal<=2; familyOrdinal++)
            {
              int fieldOrdinal = faceFieldOrdinalOffset + familyOrdinal - 1;
              const auto &s2 = lambda[ face_family_end_[familyOrdinal-1]];
              for (int ij_sum=1; ij_sum <= max_ij_sum; ij_sum++)
              {
                for (int i=0; i<ij_sum; i++)
                {
                  const int j = ij_sum - i; // j >= 1
                  // family 1 involves edge functions from edge (0,1) (edgeOrdinal 0); family 2 involves functions from edge (1,2) (edgeOrdinal 1)
                  const int edgeBasisOrdinal = i + (familyOrdinal-1)*num1DEdgeFunctions;
                  const auto & edgeValue_x = output_(edgeBasisOrdinal,pointOrdinal,0);
                  const auto & edgeValue_y = output_(edgeBasisOrdinal,pointOrdinal,1);
                  const double alpha = i*2.0 + 1;
                  
                  Polynomials::shiftedScaledIntegratedJacobiValues(jacobi_values_at_point, alpha, polyOrder_-1, s2, jacobiScaling);
                  const auto & jacobiValue = jacobi_values_at_point(j);
                  output_(fieldOrdinal,pointOrdinal,0) = edgeValue_x * jacobiValue;
                  output_(fieldOrdinal,pointOrdinal,1) = edgeValue_y * jacobiValue;
                  
                  fieldOrdinal += 2; // 2 because there are two face families, and we interleave them.
                }
              }
            }
          }
        } // end OPERATOR_VALUE
          break;
        case OPERATOR_CURL:
        {
          // edge functions
          int fieldOrdinalOffset = 0;
          /*
           Per Fuentes et al. (see Appendix E.1, E.2), the curls of the edge functions, are
             (i+2) * [P_i](s0,s1) * (grad s0 \times grad s1)
           The P_i we have implemented in shiftedScaledLegendreValues.
           */
          // rename the scratch memory to match our usage here:
          auto & P_i      = edge_field_values_at_point;
          auto & L_2ip1_j = jacobi_values_at_point;
          for (int edgeOrdinal=0; edgeOrdinal<numEdges; edgeOrdinal++)
          {
            const auto & s0 = lambda[edge_start_[edgeOrdinal]];
            const auto & s1 = lambda[  edge_end_[edgeOrdinal]];
            
            const auto & s0_dx = lambda_dx[edge_start_[edgeOrdinal]];
            const auto & s0_dy = lambda_dy[edge_start_[edgeOrdinal]];
            const auto & s1_dx = lambda_dx[  edge_end_[edgeOrdinal]];
            const auto & s1_dy = lambda_dy[  edge_end_[edgeOrdinal]];
            
            const OutputScalar grad_s0_cross_grad_s1 = s0_dx * s1_dy - s1_dx * s0_dy;
            
            Polynomials::shiftedScaledLegendreValues(P_i, polyOrder_-1, PointScalar(s1), PointScalar(s0+s1));
            for (int i=0; i<num1DEdgeFunctions; i++)
            {
              output_(i+fieldOrdinalOffset,pointOrdinal) = (i+2) * P_i(i) * grad_s0_cross_grad_s1;
            }
            fieldOrdinalOffset += num1DEdgeFunctions;
          }
          
          /*
           Fuentes et al give the face functions as E^f_{ij}, with curl:
             [L^{2i+1}_j](s0+s1,s2) curl(E^E_i(s0,s1)) + grad[L^(2i+1)_j](s0+s1,s2) \times E^E_i(s0,s1)
           where:
           - E^E_i is the ith edge function on the edge s0 to s1
           - L^{2i+1}_j is an shifted, scaled integrated Jacobi polynomial.
           - For family 1, s0s1s2 = 012
           - For family 2, s0s1s2 = 120
           - Note that grad[L^(2i+1)_j](s0+s1,s2) is computed as [P^{2i+1}_{j-1}](s0+s1,s2) (grad s2) + [R^{2i+1}_{j-1}] grad (s0+s1+s2),
             but for triangles (s0+s1+s2) is always 1, so that the grad (s0+s1+s2) is 0.
           - Here,
               [P^{2i+1}_{j-1}](s0,s1) = P^{2i+1}_{j-1}(s1,s0+s1)
             and
               [R^{2i+1}_{j-1}(s0,s1)] = d/dt L^{2i+1}_j(s1,s0+s1)
             We have implemented P^{alpha}_{j} as shiftedScaledJacobiValues,
             and d/dt L^{alpha}_{j} as shiftedScaledIntegratedJacobiValues_dt.
           */
          // rename the scratch memory to match our usage here:
          auto & P_2ip1_j    = other_values_at_point;
          
          // following ESEAS, we interleave the face families.  This groups all the face dofs of a given degree together.
          const int faceFieldOrdinalOffset = fieldOrdinalOffset;
          for (int familyOrdinal=1; familyOrdinal<=2; familyOrdinal++)
          {
            int fieldOrdinal = faceFieldOrdinalOffset + familyOrdinal - 1;
            
            const auto &s0_index = face_family_start_ [familyOrdinal-1];
            const auto &s1_index = face_family_middle_[familyOrdinal-1];
            const auto &s2_index = face_family_end_   [familyOrdinal-1];
            const auto &s0 = lambda[s0_index];
            const auto &s1 = lambda[s1_index];
            const auto &s2 = lambda[s2_index];
            const double jacobiScaling = 1.0; // s0 + s1 + s2
            
            const auto & s0_dx = lambda_dx[s0_index];
            const auto & s0_dy = lambda_dy[s0_index];
            const auto & s1_dx = lambda_dx[s1_index];
            const auto & s1_dy = lambda_dy[s1_index];
            const auto & s2_dx = lambda_dx[s2_index];
            const auto & s2_dy = lambda_dy[s2_index];
            
            const OutputScalar grad_s0_cross_grad_s1 = s0_dx * s1_dy - s1_dx * s0_dy;
            
            Polynomials::shiftedScaledLegendreValues (P_i, polyOrder_-1, PointScalar(s1), PointScalar(s0+s1));
            // [L^{2i+1}_j](s0+s1,s2) curl(E^E_i(s0,s1)) + grad[L^(2i+1)_j](s0+s1,s2) \times E^E_i(s0,s1)
            // grad[L^(2i+1)_j](s0+s1,s2) \times E^E_i(s0,s1)
//                - Note that grad[L^(2i+1)_j](s0+s1,s2) is computed as [P^{2i+1}_{j-1}](s0+s1,s2) (grad s2) + [R^{2i+1}_{j-1}] grad (s0+s1+s2),
            const PointScalar xEdgeWeight = s0 * s1_dx - s1 * s0_dx;
            const PointScalar yEdgeWeight = s0 * s1_dy - s1 * s0_dy;
            OutputScalar grad_s2_cross_xy_edgeWeight = s2_dx * yEdgeWeight - xEdgeWeight * s2_dy;
            
            const int max_ij_sum = polyOrder_ - 1;
            for (int ij_sum=1; ij_sum <= max_ij_sum; ij_sum++)
            {
              for (int i=0; i<ij_sum; i++)
              {
                const int j = ij_sum - i; // j >= 1
                const OutputScalar edgeCurl = (i+2.) * P_i(i) * grad_s0_cross_grad_s1;
              
                const double alpha = i*2.0 + 1;
                
                Polynomials::shiftedScaledJacobiValues(P_2ip1_j, alpha, polyOrder_-1, PointScalar(s2), jacobiScaling);
                Polynomials::shiftedScaledIntegratedJacobiValues(L_2ip1_j, alpha, polyOrder_-1, s2, jacobiScaling);
              
                const PointScalar & edgeValue = P_i(i);
                output_(fieldOrdinal,pointOrdinal) = L_2ip1_j(j) * edgeCurl + P_2ip1_j(j-1) * edgeValue * grad_s2_cross_xy_edgeWeight;
                
                fieldOrdinal += 2; // 2 because there are two face families, and we interleave them.
              }
            }
          }
        }
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
                                   ">>> ERROR: (Intrepid2::Hierarchical_HCURL_TRI_Functor) Unsupported differential operator");
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
      // we will use shared memory to create a fast buffer for basis computations
      size_t shmem_size = 0;
      if (fad_size_output_ > 0)
        shmem_size += 4 * OutputScratchView::shmem_size(polyOrder_ + 1, fad_size_output_);
      else
        shmem_size += 4 * OutputScratchView::shmem_size(polyOrder_ + 1);
      
      return shmem_size;
    }
  };
  
  /** \class  Intrepid2::HierarchicalBasis_HCURL_TRI
      \brief  For mathematical details of the construction, see:
   
               Federico Fuentes, Brendan Keith, Leszek Demkowicz, Sriram Nagaraj.
               "Orientation embedded high order shape functions for the exact sequence elements of all shapes."
               Computers & Mathematics with Applications, Volume 70, Issue 4, 2015, Pages 353-458, ISSN 0898-1221.
               https://doi.org/10.1016/j.camwa.2015.04.027.
   
               Note that the template argument defineVertexFunctions controls whether this basis is defined in a
               strictly hierarchical way.  If defineVertexFunctions is false, then the first basis function is the
               constant 1.0, and this basis will be suitable for discontinuous discretizations.  If defineVertexFunctions
               is true, then the first basis function will instead be 1.0-x-y, and the basis will be suitable for
               continuous discretizations.
  */
  template<typename DeviceType,
           typename OutputScalar = double,
           typename PointScalar  = double,
           bool useCGBasis = true> // if useCGBasis is false, all basis functions will be associated with the interior
  class HierarchicalBasis_HCURL_TRI
  : public Basis<DeviceType,OutputScalar,PointScalar>
  {
  public:
    using BasisBase = Basis<DeviceType,OutputScalar,PointScalar>;
    
    using HostBasis = HierarchicalBasis_HCURL_TRI<typename Kokkos::HostSpace::device_type, OutputScalar, PointScalar, useCGBasis>;

    using typename BasisBase::OrdinalTypeArray1DHost;
    using typename BasisBase::OrdinalTypeArray2DHost;

    using typename BasisBase::OutputViewType;
    using typename BasisBase::PointViewType;
    using typename BasisBase::ScalarViewType;

    using typename BasisBase::ExecutionSpace;

  protected:
    int polyOrder_; // the maximum order of the polynomial
  public:
    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order of the basis.
        \param [in] pointType - point type for nodal basis.  Ignored here (irrelevant for hierarchical/modal basis).
     */
    HierarchicalBasis_HCURL_TRI(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    polyOrder_(polyOrder)
    {
      const int numEdgeFunctions  = polyOrder * 3;
      const int numFaceFunctions  = polyOrder * (polyOrder-1);  // two families, each with p*(p-1)/2 functions
      this->basisCardinality_     = numEdgeFunctions + numFaceFunctions;
      this->basisDegree_          = polyOrder;
      this->basisCellTopologyKey_ = shards::Triangle<>::key;
      this->basisType_            = BASIS_FEM_HIERARCHICAL;
      this->basisCoordinates_     =  COORDINATES_CARTESIAN;
      this->functionSpace_        = FUNCTION_SPACE_HCURL;
      
      const int degreeLength = 1;
      this->fieldOrdinalPolynomialDegree_ = OrdinalTypeArray2DHost("Hierarchical H(curl) triangle polynomial degree lookup", this->basisCardinality_, degreeLength);
      this->fieldOrdinalH1PolynomialDegree_ = OrdinalTypeArray2DHost("Hierarchical H(curl) triangle polynomial H^1 degree lookup", this->basisCardinality_, degreeLength);
      
      int fieldOrdinalOffset = 0;
      // **** vertex functions **** //
      // no vertex functions in H(curl)
      
      // **** edge functions **** //
      const shards::CellTopology cellTopo(shards::getCellTopologyData<shards::Triangle<> >());
      const int numFunctionsPerEdge = polyOrder; // p functions associated with each edge
      const int numEdges            = cellTopo.getEdgeCount();
      for (int edgeOrdinal=0; edgeOrdinal<numEdges; edgeOrdinal++)
      {
        for (int i=0; i<numFunctionsPerEdge; i++)
        {
          this->fieldOrdinalPolynomialDegree_(i+fieldOrdinalOffset,0) = i+1; // the multiplicands involving the gradients of the vertex functions are first degree polynomials; hence the +1 (the remaining multiplicands are order i = 0,…,p-1).
          this->fieldOrdinalH1PolynomialDegree_(i+fieldOrdinalOffset,0) = i+2;
        }
        fieldOrdinalOffset += numFunctionsPerEdge;
      }
      INTREPID2_TEST_FOR_EXCEPTION(fieldOrdinalOffset != numEdgeFunctions, std::invalid_argument, "Internal error: basis enumeration is incorrect");
      
      // **** face functions **** //
      const int max_ij_sum = polyOrder-1;
      const int faceFieldOrdinalOffset = fieldOrdinalOffset;
      for (int faceFamilyOrdinal=1; faceFamilyOrdinal<=2; faceFamilyOrdinal++)
      {
        // following ESEAS, we interleave the face families.  This groups all the face dofs of a given degree together.
        int fieldOrdinal = faceFieldOrdinalOffset + faceFamilyOrdinal - 1;
        for (int ij_sum=1; ij_sum <= max_ij_sum; ij_sum++)
        {
          for (int i=0; i<ij_sum; i++)
          {
            this->fieldOrdinalPolynomialDegree_(fieldOrdinal,0) = ij_sum+1;
            this->fieldOrdinalH1PolynomialDegree_(fieldOrdinal,0) = ij_sum+2;
            fieldOrdinal += 2; // 2 because there are two face families, and we interleave them.
          }
        }
        fieldOrdinalOffset = fieldOrdinal - 1; // due to the interleaving increment, we've gone two past the last face ordinal.  Set offset to be one past.
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
        const ordinal_type edgeDim = 1, faceDim = 2;

        if (useCGBasis) {
          {
            int tagNumber = 0;
            for (int edgeOrdinal=0; edgeOrdinal<numEdges; edgeOrdinal++)
            {
              for (int functionOrdinal=0; functionOrdinal<numFunctionsPerEdge; functionOrdinal++)
              {
                tagView(tagNumber*tagSize+0) = edgeDim;               // edge dimension
                tagView(tagNumber*tagSize+1) = edgeOrdinal;           // edge id
                tagView(tagNumber*tagSize+2) = functionOrdinal;       // local dof id
                tagView(tagNumber*tagSize+3) = numFunctionsPerEdge;   // total number of dofs on this edge
                tagNumber++;
              }
            }
            const int numFunctionsPerFace = numFaceFunctions; // just one face in the triangle
            const int numFaces = 1;
            for (int faceOrdinal=0; faceOrdinal<numFaces; faceOrdinal++)
            {
              for (int functionOrdinal=0; functionOrdinal<numFunctionsPerFace; functionOrdinal++)
              {
                tagView(tagNumber*tagSize+0) = faceDim;               // face dimension
                tagView(tagNumber*tagSize+1) = faceOrdinal;           // face id
                tagView(tagNumber*tagSize+2) = functionOrdinal;       // local dof id
                tagView(tagNumber*tagSize+3) = numFunctionsPerFace;   // total number of dofs on this face
                tagNumber++;
              }
            }
          }
        }
        else
        {
          // DG basis: all functions are associated with interior
          for (ordinal_type i=0;i<cardinality;++i) {
            tagView(i*tagSize+0) = faceDim;     // face dimension
            tagView(i*tagSize+1) = 0;           // face id
            tagView(i*tagSize+2) = i;           // local dof id
            tagView(i*tagSize+3) = cardinality; // total number of dofs on this face
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
    
    /** \brief  Returns basis name
     
     \return the name of the basis
     */
    const char* getName() const override {
      return "Intrepid2_HierarchicalBasis_HCURL_TRI";
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
      
      using FunctorType = Hierarchical_HCURL_TRI_Functor<DeviceType, OutputScalar, PointScalar, OutputViewType, PointViewType>;
      
      FunctorType functor(operatorType, outputValues, inputPoints, polyOrder_);
      
      const int outputVectorSize = getVectorSizeForHierarchicalParallelism<OutputScalar>();
      const int pointVectorSize  = getVectorSizeForHierarchicalParallelism<PointScalar>();
      const int vectorSize = std::max(outputVectorSize,pointVectorSize);
      const int teamSize = 1; // because of the way the basis functions are computed, we don't have a second level of parallelism...

      auto policy = Kokkos::TeamPolicy<ExecutionSpace>(numPoints,teamSize,vectorSize);
      Kokkos::parallel_for("Hierarchical_HCURL_TRI_Functor", policy, functor);
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
            LegendreBasis_HVOL_LINE<DeviceType,OutputScalar,PointScalar>
                    (this->basisDegree_-1));
      }
      INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Input parameters out of bounds");
    }

    /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.
     
        \return Pointer to the new Basis object.
     */
    virtual BasisPtr<typename Kokkos::HostSpace::device_type, OutputScalar, PointScalar>
    getHostBasis() const override {
      using HostDeviceType = typename Kokkos::HostSpace::device_type;
      using HostBasisType  = HierarchicalBasis_HCURL_TRI<HostDeviceType, OutputScalar, PointScalar, useCGBasis>;
      return Teuchos::rcp( new HostBasisType(polyOrder_) );
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_HierarchicalBasis_HCURL_TRI_h */
