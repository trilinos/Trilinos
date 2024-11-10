// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_IntegratedLegendreBasis_HGRAD_TRI.hpp
    \brief  H(grad) basis on the triangle based on integrated Legendre polynomials.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_IntegratedLegendreBasis_HGRAD_TRI_h
#define Intrepid2_IntegratedLegendreBasis_HGRAD_TRI_h

#include <Kokkos_DynRankView.hpp>

#include <Intrepid2_config.h>

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_IntegratedLegendreBasis_HGRAD_LINE.hpp"
#include "Intrepid2_Polynomials.hpp"
#include "Intrepid2_Utils.hpp"

namespace Intrepid2
{
  /** \class  Intrepid2::Hierarchical_HGRAD_TRI_Functor
      \brief  Functor for computing values for the IntegratedLegendreBasis_HGRAD_TRI class.
   
   This functor is not intended for use outside of IntegratedLegendreBasis_HGRAD_TRI.
  */
  template<class DeviceType, class OutputScalar, class PointScalar,
           class OutputFieldType, class InputPointsType>
  struct Hierarchical_HGRAD_TRI_Functor
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
    bool defineVertexFunctions_;
    int numFields_, numPoints_;
    
    size_t fad_size_output_;
    
    static const int numVertices = 3;
    static const int numEdges    = 3;
    const int edge_start_[numEdges] = {0,1,0}; // edge i is from edge_start_[i] to edge_end_[i]
    const int edge_end_[numEdges]   = {1,2,2}; // edge i is from edge_start_[i] to edge_end_[i]
    
    Hierarchical_HGRAD_TRI_Functor(EOperator opType, OutputFieldType output, InputPointsType inputPoints,
                                    int polyOrder, bool defineVertexFunctions)
    : opType_(opType), output_(output), inputPoints_(inputPoints),
      polyOrder_(polyOrder), defineVertexFunctions_(defineVertexFunctions),
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
      OutputScratchView edge_field_values_at_point, jacobi_values_at_point, other_values_at_point, other_values2_at_point;
      if (fad_size_output_ > 0) {
        edge_field_values_at_point = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1, fad_size_output_);
        jacobi_values_at_point     = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1, fad_size_output_);
        other_values_at_point      = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1, fad_size_output_);
        other_values2_at_point     = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1, fad_size_output_);
      }
      else {
        edge_field_values_at_point = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1);
        jacobi_values_at_point     = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1);
        other_values_at_point      = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1);
        other_values2_at_point     = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1);
      }
      
      const auto & x = inputPoints_(pointOrdinal,0);
      const auto & y = inputPoints_(pointOrdinal,1);
      
      // write as barycentric coordinates:
      const PointScalar lambda[3]    = {1. - x - y, x, y};
      const PointScalar lambda_dx[3] = {-1., 1., 0.};
      const PointScalar lambda_dy[3] = {-1., 0., 1.};
      
      const int num1DEdgeFunctions = polyOrder_ - 1;
      
      switch (opType_)
      {
        case OPERATOR_VALUE:
        {
          // vertex functions come first, according to vertex ordering: (0,0), (1,0), (0,1)
          for (int vertexOrdinal=0; vertexOrdinal<numVertices; vertexOrdinal++)
          {
            output_(vertexOrdinal,pointOrdinal) = lambda[vertexOrdinal];
          }
          if (!defineVertexFunctions_)
          {
            // "DG" basis case
            // here, we overwrite the first vertex function with 1:
            output_(0,pointOrdinal) = 1.0;
          }
          
          // edge functions
          int fieldOrdinalOffset = 3;
          for (int edgeOrdinal=0; edgeOrdinal<numEdges; edgeOrdinal++)
          {
            const auto & s0 = lambda[edge_start_[edgeOrdinal]];
            const auto & s1 = lambda[  edge_end_[edgeOrdinal]];
            
            Polynomials::shiftedScaledIntegratedLegendreValues(edge_field_values_at_point, polyOrder_, PointScalar(s1), PointScalar(s0+s1));
            for (int edgeFunctionOrdinal=0; edgeFunctionOrdinal<num1DEdgeFunctions; edgeFunctionOrdinal++)
            {
              // the first two integrated legendre functions are essentially the vertex functions; hence the +2 on on the RHS here:
              output_(edgeFunctionOrdinal+fieldOrdinalOffset,pointOrdinal) = edge_field_values_at_point(edgeFunctionOrdinal+2);
            }
            fieldOrdinalOffset += num1DEdgeFunctions;
          }
          
          // face functions
          {
            // these functions multiply the edge functions from the 01 edge by integrated Jacobi functions, appropriately scaled
            const double jacobiScaling = 1.0; // s0 + s1 + s2
            
            const int max_ij_sum = polyOrder_;
            const int min_i = 2;
            const int min_j = 1;
            const int min_ij_sum = min_i + min_j;
            for (int ij_sum = min_ij_sum; ij_sum <= max_ij_sum; ij_sum++)
            {
              for (int i=min_i; i<=ij_sum-min_j; i++)
              {
                const int j = ij_sum - i;
                const int edgeBasisOrdinal = i+numVertices-2; // i+1: where the value of the edge function is stored in output_
                const auto & edgeValue = output_(edgeBasisOrdinal,pointOrdinal);
                const double alpha = i*2.0;
                
                Polynomials::shiftedScaledIntegratedJacobiValues(jacobi_values_at_point, alpha, polyOrder_-2, lambda[2], jacobiScaling);
                const auto & jacobiValue = jacobi_values_at_point(j);
                output_(fieldOrdinalOffset,pointOrdinal) = edgeValue * jacobiValue;
                fieldOrdinalOffset++;
              }
            }
          }
        } // end OPERATOR_VALUE
          break;
        case OPERATOR_GRAD:
        case OPERATOR_D1:
        {
          // vertex functions
          if (defineVertexFunctions_)
          {
            // standard, "CG" basis case
            // first vertex function is 1-x-y
            output_(0,pointOrdinal,0) = -1.0;
            output_(0,pointOrdinal,1) = -1.0;
          }
          else
          {
            // "DG" basis case
            // here, the first "vertex" function is 1, so the derivative is 0:
            output_(0,pointOrdinal,0) = 0.0;
            output_(0,pointOrdinal,1) = 0.0;
          }
          // second vertex function is x
          output_(1,pointOrdinal,0) = 1.0;
          output_(1,pointOrdinal,1) = 0.0;
          // third vertex function is y
          output_(2,pointOrdinal,0) = 0.0;
          output_(2,pointOrdinal,1) = 1.0;
          
          // edge functions
          int fieldOrdinalOffset = 3;
          /*
           Per Fuentes et al. (see Appendix E.1, E.2), the edge functions, defined for i ≥ 2, are
             [L_i](s0,s1) = L_i(s1; s0+s1)
           and have gradients:
             grad [L_i](s0,s1) = [P_{i-1}](s0,s1) grad s1 + [R_{i-1}](s0,s1) grad (s0 + s1)
           where
             [R_{i-1}](s0,s1) = R_{i-1}(s1; s0+s1) = d/dt L_{i}(s0; s0+s1)
           The P_i we have implemented in shiftedScaledLegendreValues, while d/dt L_{i+1} is
           implemented in shiftedScaledIntegratedLegendreValues_dt.
           */
          // rename the scratch memory to match our usage here:
          auto & P_i_minus_1 = edge_field_values_at_point;
          auto & L_i_dt      = jacobi_values_at_point;
          for (int edgeOrdinal=0; edgeOrdinal<numEdges; edgeOrdinal++)
          {
            const auto & s0 = lambda[edge_start_[edgeOrdinal]];
            const auto & s1 = lambda[  edge_end_[edgeOrdinal]];
            
            const auto & s0_dx = lambda_dx[edge_start_[edgeOrdinal]];
            const auto & s0_dy = lambda_dy[edge_start_[edgeOrdinal]];
            const auto & s1_dx = lambda_dx[  edge_end_[edgeOrdinal]];
            const auto & s1_dy = lambda_dy[  edge_end_[edgeOrdinal]];
            
            Polynomials::shiftedScaledLegendreValues             (P_i_minus_1, polyOrder_-1, PointScalar(s1), PointScalar(s0+s1));
            Polynomials::shiftedScaledIntegratedLegendreValues_dt(L_i_dt,      polyOrder_,   PointScalar(s1), PointScalar(s0+s1));
            for (int edgeFunctionOrdinal=0; edgeFunctionOrdinal<num1DEdgeFunctions; edgeFunctionOrdinal++)
            {
              // the first two (integrated) Legendre functions are essentially the vertex functions; hence the +2 here:
              const int i = edgeFunctionOrdinal+2;
              output_(edgeFunctionOrdinal+fieldOrdinalOffset,pointOrdinal,0) = P_i_minus_1(i-1) * s1_dx + L_i_dt(i) * (s1_dx + s0_dx);
              output_(edgeFunctionOrdinal+fieldOrdinalOffset,pointOrdinal,1) = P_i_minus_1(i-1) * s1_dy + L_i_dt(i) * (s1_dy + s0_dy);
            }
            fieldOrdinalOffset += num1DEdgeFunctions;
          }
          
          /*
           Fuentes et al give the face functions as phi_{ij}, with gradient:
             grad phi_{ij}(s0,s1,s2) = [L^{2i}_j](s0+s1,s2) grad [L_i](s0,s1) + [L_i](s0,s1) grad [L^{2i}_j](s0+s1,s2)
           where:
           - grad [L_i](s0,s1) is the edge function gradient we computed above
           - [L_i](s0,s1) is the edge function which we have implemented above (in OPERATOR_VALUE)
           - L^{2i}_j is a Jacobi polynomial with:
               [L^{2i}_j](s0,s1) = L^{2i}_j(s1;s0+s1)
             and the gradient for j ≥ 1 is
               grad [L^{2i}_j](s0,s1) = [P^{2i}_{j-1}](s0,s1) grad s1 + [R^{2i}_{j-1}(s0,s1)] grad (s0 + s1)
           Here,
             [P^{2i}_{j-1}](s0,s1) = P^{2i}_{j-1}(s1,s0+s1)
           and
             [R^{2i}_{j-1}(s0,s1)] = d/dt L^{2i}_j(s1,s0+s1)
           We have implemented P^{alpha}_{j} as shiftedScaledJacobiValues,
           and d/dt L^{alpha}_{j} as shiftedScaledIntegratedJacobiValues_dt.
           */
          // rename the scratch memory to match our usage here:
          auto & P_2i_j_minus_1 = edge_field_values_at_point;
          auto & L_2i_j_dt      = jacobi_values_at_point;
          auto & L_i            = other_values_at_point;
          auto & L_2i_j         = other_values2_at_point;
          {
            // face functions multiply the edge functions from the 01 edge by integrated Jacobi functions, appropriately scaled
            const double jacobiScaling = 1.0; // s0 + s1 + s2

            const int max_ij_sum = polyOrder_;
            const int min_i = 2;
            const int min_j = 1;
            const int min_ij_sum = min_i + min_j;
            for (int ij_sum = min_ij_sum; ij_sum <= max_ij_sum; ij_sum++)
            {
              for (int i=min_i; i<=ij_sum-min_j; i++)
              {
                const int j = ij_sum - i;
                // the edge function here is for edge 01, in the first set of edge functions.
                const int edgeBasisOrdinal = i+numVertices-2; // i+1: where the value of the edge function is stored in output_
                const auto & grad_L_i_dx = output_(edgeBasisOrdinal,pointOrdinal,0);
                const auto & grad_L_i_dy = output_(edgeBasisOrdinal,pointOrdinal,1);
                
                const double alpha = i*2.0;

                Polynomials::shiftedScaledIntegratedLegendreValues (L_i, polyOrder_, lambda[1], lambda[0]+lambda[1]);
                Polynomials::shiftedScaledIntegratedJacobiValues_dt(L_2i_j_dt, alpha, polyOrder_, lambda[2], jacobiScaling);
                Polynomials::shiftedScaledIntegratedJacobiValues   (   L_2i_j, alpha, polyOrder_, lambda[2], jacobiScaling);
                Polynomials::shiftedScaledJacobiValues(P_2i_j_minus_1, alpha, polyOrder_-1, lambda[2], jacobiScaling);
                
                const auto & s0_dx = lambda_dx[0];
                const auto & s0_dy = lambda_dy[0];
                const auto & s1_dx = lambda_dx[1];
                const auto & s1_dy = lambda_dy[1];
                const auto & s2_dx = lambda_dx[2];
                const auto & s2_dy = lambda_dy[2];
                
                const OutputScalar basisValue_dx = L_2i_j(j) * grad_L_i_dx + L_i(i) * (P_2i_j_minus_1(j-1) * s2_dx + L_2i_j_dt(j) * (s0_dx + s1_dx + s2_dx));
                const OutputScalar basisValue_dy = L_2i_j(j) * grad_L_i_dy + L_i(i) * (P_2i_j_minus_1(j-1) * s2_dy + L_2i_j_dt(j) * (s0_dy + s1_dy + s2_dy));
                
                output_(fieldOrdinalOffset,pointOrdinal,0) = basisValue_dx;
                output_(fieldOrdinalOffset,pointOrdinal,1) = basisValue_dy;
                fieldOrdinalOffset++;
              }
            }
          }
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
          INTREPID2_TEST_FOR_ABORT( true,
                                   ">>> ERROR: (Intrepid2::Basis_HGRAD_TRI_Cn_FEM_ORTH::OrthPolynomialTri) Computing of second and higher-order derivatives is not currently supported");
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
  
  /** \class  Intrepid2::IntegratedLegendreBasis_HGRAD_TRI
      \brief  Basis defining integrated Legendre basis on the line, a polynomial subspace of H(grad) on the line: extension to triangle using Jacobi blending functions.

              For mathematical details of the construction, see:
   
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
           bool defineVertexFunctions = true>            // if defineVertexFunctions is true, first three basis functions are 1-x-y, x, and y.  Otherwise, they are 1, x, and y.
  class IntegratedLegendreBasis_HGRAD_TRI
  : public Basis<DeviceType,OutputScalar,PointScalar>
  {
  public:
    using BasisBase = Basis<DeviceType,OutputScalar,PointScalar>;
    using HostBasis = IntegratedLegendreBasis_HGRAD_TRI<typename Kokkos::HostSpace::device_type,OutputScalar,PointScalar,defineVertexFunctions>;

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
     
     The basis will have polyOrder + 1 members.
     
     If defineVertexFunctions is false, then all basis functions are identified with the interior of the line element, and the first two basis functions are 1 and x.
     
     If defineVertexFunctions is true, then the first two basis functions are 1-x and x, and these are identified with the left and right vertices of the cell.
     
     */
    IntegratedLegendreBasis_HGRAD_TRI(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    polyOrder_(polyOrder),
    pointType_(pointType)
    {
      INTREPID2_TEST_FOR_EXCEPTION(pointType!=POINTTYPE_DEFAULT,std::invalid_argument,"PointType not supported");

      this->basisCardinality_     = ((polyOrder+2) * (polyOrder+1)) / 2;
      this->basisDegree_          = polyOrder;
      this->basisCellTopologyKey_ = shards::Triangle<>::key;
      this->basisType_            = BASIS_FEM_HIERARCHICAL;
      this->basisCoordinates_     = COORDINATES_CARTESIAN;
      this->functionSpace_        = FUNCTION_SPACE_HGRAD;
      
      const int degreeLength = 1;
      this->fieldOrdinalPolynomialDegree_ = OrdinalTypeArray2DHost("Integrated Legendre H(grad) triangle polynomial degree lookup", this->basisCardinality_, degreeLength);
      this->fieldOrdinalH1PolynomialDegree_ = OrdinalTypeArray2DHost("Integrated Legendre H(grad) triangle polynomial degree lookup", this->basisCardinality_, degreeLength);
      
      int fieldOrdinalOffset = 0;
      // **** vertex functions **** //
      const int numVertices = this->getBaseCellTopology().getVertexCount();
      const int numFunctionsPerVertex = 1;
      const int numVertexFunctions = numVertices * numFunctionsPerVertex;
      for (int i=0; i<numVertexFunctions; i++)
      {
        // for H(grad) on triangle, if defineVertexFunctions is false, first three basis members are linear
        // if not, then the only difference is that the first member is constant
        this->fieldOrdinalPolynomialDegree_  (i,0) = 1;
        this->fieldOrdinalH1PolynomialDegree_(i,0) = 1;
      }
      if (!defineVertexFunctions)
      {
        this->fieldOrdinalPolynomialDegree_  (0,0) = 0;
        this->fieldOrdinalH1PolynomialDegree_(0,0) = 0;
      }
      fieldOrdinalOffset += numVertexFunctions;
      
      // **** edge functions **** //
      const int numFunctionsPerEdge = polyOrder - 1; // bubble functions: all but the vertices
      const int numEdges            = this->getBaseCellTopology().getEdgeCount();
      for (int edgeOrdinal=0; edgeOrdinal<numEdges; edgeOrdinal++)
      {
        for (int i=0; i<numFunctionsPerEdge; i++)
        {
          this->fieldOrdinalPolynomialDegree_(i+fieldOrdinalOffset,0)   = i+2; // vertex functions are 1st order; edge functions start at order 2
          this->fieldOrdinalH1PolynomialDegree_(i+fieldOrdinalOffset,0) = i+2; // vertex functions are 1st order; edge functions start at order 2
        }
        fieldOrdinalOffset += numFunctionsPerEdge;
      }
      
      // **** face functions **** //
      const int max_ij_sum = polyOrder;
      const int min_i = 2;
      const int min_j = 1;
      const int min_ij_sum = min_i + min_j;
      for (int ij_sum = min_ij_sum; ij_sum <= max_ij_sum; ij_sum++)
      {
        for (int i=min_i; i<=ij_sum-min_j; i++)
        {
          const int j = ij_sum - i;
          this->fieldOrdinalPolynomialDegree_(fieldOrdinalOffset,0)   = i+j;
          this->fieldOrdinalH1PolynomialDegree_(fieldOrdinalOffset,0) = i+j;
          fieldOrdinalOffset++;
        }
      }
      const int numFaces = 1;
      const int numFunctionsPerFace = ((polyOrder-1)*(polyOrder-2))/2;
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
        const int vertexDim = 0, edgeDim = 1, faceDim = 2;

        if (defineVertexFunctions) {
          {
            int tagNumber = 0;
            for (int vertexOrdinal=0; vertexOrdinal<numVertices; vertexOrdinal++)
            {
              for (int functionOrdinal=0; functionOrdinal<numFunctionsPerVertex; functionOrdinal++)
              {
                tagView(tagNumber*tagSize+0) = vertexDim;             // vertex dimension
                tagView(tagNumber*tagSize+1) = vertexOrdinal;         // vertex id
                tagView(tagNumber*tagSize+2) = functionOrdinal;       // local dof id
                tagView(tagNumber*tagSize+3) = numFunctionsPerVertex; // total number of dofs in this vertex
                tagNumber++;
              }
            }
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
        } else {
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
      return "Intrepid2_IntegratedLegendreBasis_HGRAD_TRI";
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
      
      using FunctorType = Hierarchical_HGRAD_TRI_Functor<DeviceType, OutputScalar, PointScalar, OutputViewType, PointViewType>;
      
      FunctorType functor(operatorType, outputValues, inputPoints, polyOrder_, defineVertexFunctions);
      
      const int outputVectorSize = getVectorSizeForHierarchicalParallelism<OutputScalar>();
      const int pointVectorSize  = getVectorSizeForHierarchicalParallelism<PointScalar>();
      const int vectorSize = std::max(outputVectorSize,pointVectorSize);
      const int teamSize = 1; // because of the way the basis functions are computed, we don't have a second level of parallelism...

      auto policy = Kokkos::TeamPolicy<ExecutionSpace>(numPoints,teamSize,vectorSize);
      Kokkos::parallel_for("Hierarchical_HGRAD_TRI_Functor", policy, functor);
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
      using HostBasisType  = IntegratedLegendreBasis_HGRAD_TRI<HostDeviceType, OutputScalar, PointScalar, defineVertexFunctions>;
      return Teuchos::rcp( new HostBasisType(polyOrder_, pointType_) );
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_IntegratedLegendreBasis_HGRAD_TRI_h */
