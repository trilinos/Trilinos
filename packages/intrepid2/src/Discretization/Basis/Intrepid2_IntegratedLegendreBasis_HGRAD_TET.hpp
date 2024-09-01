// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_IntegratedLegendreBasis_HGRAD_TET.hpp
    \brief  H(grad) basis on the tetrahedon based on integrated Legendre polynomials.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_IntegratedLegendreBasis_HGRAD_TET_h
#define Intrepid2_IntegratedLegendreBasis_HGRAD_TET_h

#include <Kokkos_DynRankView.hpp>

#include <Intrepid2_config.h>

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_IntegratedLegendreBasis_HGRAD_LINE.hpp"
#include "Intrepid2_IntegratedLegendreBasis_HGRAD_TRI.hpp"
#include "Intrepid2_Polynomials.hpp"
#include "Intrepid2_Utils.hpp"

namespace Intrepid2
{
  /** \class  Intrepid2::Hierarchical_HGRAD_TET_Functor
      \brief  Functor for computing values for the IntegratedLegendreBasis_HGRAD_TET class.
   
   This functor is not intended for use outside of IntegratedLegendreBasis_HGRAD_TET.
  */
  template<class DeviceType, class OutputScalar, class PointScalar,
           class OutputFieldType, class InputPointsType>
  struct Hierarchical_HGRAD_TET_Functor
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
    bool defineVertexFunctions_;
    int numFields_, numPoints_;
    
    size_t fad_size_output_;
    
    static const int numVertices = 4;
    static const int numEdges    = 6;
    // the following ordering of the edges matches that used by ESEAS
    const int edge_start_[numEdges] = {0,1,0,0,1,2}; // edge i is from edge_start_[i] to edge_end_[i]
    const int edge_end_[numEdges]   = {1,2,2,3,3,3}; // edge i is from edge_start_[i] to edge_end_[i]
    
    static const int numFaces    = 4;
    const int face_vertex_0[numFaces] = {0,0,1,0}; // faces are abc where 0 ≤ a < b < c ≤ 3
    const int face_vertex_1[numFaces] = {1,1,2,2};
    const int face_vertex_2[numFaces] = {2,3,3,3};
    
    // this allows us to look up the edge ordinal of the first edge of a face
    // this is useful because face functions are defined using edge basis functions of the first edge of the face
    const int face_ordinal_of_first_edge[numFaces] = {0,0,1,2};
    
    Hierarchical_HGRAD_TET_Functor(EOperator opType, OutputFieldType output, InputPointsType inputPoints,
                                    int polyOrder, bool defineVertexFunctions)
    : opType_(opType), output_(output), inputPoints_(inputPoints),
      polyOrder_(polyOrder), defineVertexFunctions_(defineVertexFunctions),
      fad_size_output_(getScalarDimensionForView(output))
    {
      numFields_ = output.extent_int(0);
      numPoints_ = output.extent_int(1);
      INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != inputPoints.extent_int(0), std::invalid_argument, "point counts need to match!");
      INTREPID2_TEST_FOR_EXCEPTION(numFields_ != (polyOrder_+1)*(polyOrder_+2)*(polyOrder_+3)/6, std::invalid_argument, "output field size does not match basis cardinality");
    }
    
    KOKKOS_INLINE_FUNCTION
    void operator()( const TeamMember & teamMember ) const
    {
      const int numFaceBasisFunctionsPerFace = (polyOrder_-2) * (polyOrder_-1) / 2;
      const int numInteriorBasisFunctions = (polyOrder_-3) * (polyOrder_-2) * (polyOrder_-1) / 6;
      
      auto pointOrdinal = teamMember.league_rank();
      OutputScratchView legendre_values1_at_point, legendre_values2_at_point;
      OutputScratchView2D jacobi_values1_at_point, jacobi_values2_at_point, jacobi_values3_at_point;
      const int numAlphaValues = (polyOrder_-1 > 1) ? (polyOrder_-1) : 1; // make numAlphaValues at least 1 so we can avoid zero-extent allocations…
      if (fad_size_output_ > 0) {
        legendre_values1_at_point = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1, fad_size_output_);
        legendre_values2_at_point = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1, fad_size_output_);
        jacobi_values1_at_point   = OutputScratchView2D(teamMember.team_shmem(), numAlphaValues, polyOrder_ + 1, fad_size_output_);
        jacobi_values2_at_point   = OutputScratchView2D(teamMember.team_shmem(), numAlphaValues, polyOrder_ + 1, fad_size_output_);
        jacobi_values3_at_point   = OutputScratchView2D(teamMember.team_shmem(), numAlphaValues, polyOrder_ + 1, fad_size_output_);
      }
      else {
        legendre_values1_at_point = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1);
        legendre_values2_at_point = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1);
        jacobi_values1_at_point   = OutputScratchView2D(teamMember.team_shmem(), numAlphaValues, polyOrder_ + 1);
        jacobi_values2_at_point   = OutputScratchView2D(teamMember.team_shmem(), numAlphaValues, polyOrder_ + 1);
        jacobi_values3_at_point   = OutputScratchView2D(teamMember.team_shmem(), numAlphaValues, polyOrder_ + 1);
      }
      
      const auto & x = inputPoints_(pointOrdinal,0);
      const auto & y = inputPoints_(pointOrdinal,1);
      const auto & z = inputPoints_(pointOrdinal,2);
      
      // write as barycentric coordinates:
      const PointScalar lambda[numVertices] = {1. - x - y - z, x, y, z};
      const PointScalar lambda_dx[numVertices] = {-1., 1., 0., 0.};
      const PointScalar lambda_dy[numVertices] = {-1., 0., 1., 0.};
      const PointScalar lambda_dz[numVertices] = {-1., 0., 0., 1.};
      
      const int num1DEdgeFunctions = polyOrder_ - 1;
      
      switch (opType_)
      {
        case OPERATOR_VALUE:
        {
          // vertex functions come first, according to vertex ordering: (0,0,0), (1,0,0), (0,1,0), (0,0,1)
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
          int fieldOrdinalOffset = numVertices;
          for (int edgeOrdinal=0; edgeOrdinal<numEdges; edgeOrdinal++)
          {
            const auto & s0 = lambda[edge_start_[edgeOrdinal]];
            const auto & s1 = lambda[  edge_end_[edgeOrdinal]];

            Polynomials::shiftedScaledIntegratedLegendreValues(legendre_values1_at_point, polyOrder_, PointScalar(s1), PointScalar(s0+s1));
            for (int edgeFunctionOrdinal=0; edgeFunctionOrdinal<num1DEdgeFunctions; edgeFunctionOrdinal++)
            {
              // the first two integrated legendre functions are essentially the vertex functions; hence the +2 on on the RHS here:
              output_(edgeFunctionOrdinal+fieldOrdinalOffset,pointOrdinal) = legendre_values1_at_point(edgeFunctionOrdinal+2);
            }
            fieldOrdinalOffset += num1DEdgeFunctions;
          }
          /*
           Face functions for face abc are the product of edge functions on their ab edge
           and a Jacobi polynomial [L^2i_j](s0+s1,s2) = L^2i_j(s2;s0+s1+s2)
           */
          for (int faceOrdinal=0; faceOrdinal<numFaces; faceOrdinal++)
          {
            const auto & s0 = lambda[face_vertex_0[faceOrdinal]];
            const auto & s1 = lambda[face_vertex_1[faceOrdinal]];
            const auto & s2 = lambda[face_vertex_2[faceOrdinal]];
            const PointScalar jacobiScaling = s0 + s1 + s2;
            
            // compute integrated Jacobi values for each desired value of alpha
            for (int n=2; n<=polyOrder_; n++)
            {
              const double alpha = n*2;
              const int alphaOrdinal = n-2;
              using Kokkos::subview;
              using Kokkos::ALL;
              auto jacobi_alpha = subview(jacobi_values1_at_point, alphaOrdinal, ALL);
              Polynomials::shiftedScaledIntegratedJacobiValues(jacobi_alpha, alpha, polyOrder_-2, s2, jacobiScaling);
            }
            
            const int edgeOrdinal = face_ordinal_of_first_edge[faceOrdinal];
            int localFaceBasisOrdinal = 0;
            for (int totalPolyOrder=3; totalPolyOrder<=polyOrder_; totalPolyOrder++)
            {
              for (int i=2; i<totalPolyOrder; i++)
              {
                const int edgeBasisOrdinal = edgeOrdinal*num1DEdgeFunctions + i-2 + numVertices;
                const auto & edgeValue = output_(edgeBasisOrdinal,pointOrdinal);
                const int alphaOrdinal = i-2;
                
                const int j = totalPolyOrder - i;
                const auto & jacobiValue = jacobi_values1_at_point(alphaOrdinal,j);
                const int fieldOrdinal = fieldOrdinalOffset + localFaceBasisOrdinal;
                output_(fieldOrdinal,pointOrdinal) = edgeValue * jacobiValue;
                
                localFaceBasisOrdinal++;
              }
            }
            fieldOrdinalOffset += numFaceBasisFunctionsPerFace;
          }
          // interior functions
          // compute integrated Jacobi values for each desired value of alpha
          for (int n=3; n<=polyOrder_; n++)
          {
            const double alpha = n*2;
            const double jacobiScaling = 1.0;
            const int alphaOrdinal = n-3;
            using Kokkos::subview;
            using Kokkos::ALL;
            auto jacobi_alpha = subview(jacobi_values1_at_point, alphaOrdinal, ALL);
            Polynomials::shiftedScaledIntegratedJacobiValues(jacobi_alpha, alpha, polyOrder_-3, lambda[3], jacobiScaling);
          }
          const int min_i  = 2;
          const int min_j  = 1;
          const int min_k  = 1;
          const int min_ij = min_i + min_j;
          const int min_ijk = min_ij + min_k;
          int localInteriorBasisOrdinal = 0;
          for (int totalPolyOrder_ijk=min_ijk; totalPolyOrder_ijk <= polyOrder_; totalPolyOrder_ijk++)
          {
            int localFaceBasisOrdinal = 0;
            for (int totalPolyOrder_ij=min_ij; totalPolyOrder_ij <= totalPolyOrder_ijk-min_j; totalPolyOrder_ij++)
            {
              for (int i=2; i <= totalPolyOrder_ij-min_j; i++)
              {
                const int j = totalPolyOrder_ij - i;
                const int k = totalPolyOrder_ijk - totalPolyOrder_ij;
                const int faceBasisOrdinal = numEdges*num1DEdgeFunctions + numVertices + localFaceBasisOrdinal;
                const auto & faceValue = output_(faceBasisOrdinal,pointOrdinal);
                const int alphaOrdinal = (i+j)-3;
                localFaceBasisOrdinal++;
              
                const int fieldOrdinal = fieldOrdinalOffset + localInteriorBasisOrdinal;
                const auto & jacobiValue = jacobi_values1_at_point(alphaOrdinal,k);
                output_(fieldOrdinal,pointOrdinal) = faceValue * jacobiValue;
                localInteriorBasisOrdinal++;
              } // end i loop
            } // end totalPolyOrder_ij loop
          } // end totalPolyOrder_ijk loop
          fieldOrdinalOffset += numInteriorBasisFunctions;
        } // end OPERATOR_VALUE
          break;
        case OPERATOR_GRAD:
        case OPERATOR_D1:
        {
          // vertex functions
          if (defineVertexFunctions_)
          {
            // standard, "CG" basis case
            // first vertex function is 1-x-y-z
            output_(0,pointOrdinal,0) = -1.0;
            output_(0,pointOrdinal,1) = -1.0;
            output_(0,pointOrdinal,2) = -1.0;
          }
          else
          {
            // "DG" basis case
            // here, the first "vertex" function is 1, so the derivative is 0:
            output_(0,pointOrdinal,0) = 0.0;
            output_(0,pointOrdinal,1) = 0.0;
            output_(0,pointOrdinal,2) = 0.0;
          }
          // second vertex function is x
          output_(1,pointOrdinal,0) = 1.0;
          output_(1,pointOrdinal,1) = 0.0;
          output_(1,pointOrdinal,2) = 0.0;
          // third vertex function is y
          output_(2,pointOrdinal,0) = 0.0;
          output_(2,pointOrdinal,1) = 1.0;
          output_(2,pointOrdinal,2) = 0.0;
          // fourth vertex function is z
          output_(3,pointOrdinal,0) = 0.0;
          output_(3,pointOrdinal,1) = 0.0;
          output_(3,pointOrdinal,2) = 1.0;

          // edge functions
          int fieldOrdinalOffset = numVertices;
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
          auto & P_i_minus_1 = legendre_values1_at_point;
          auto & L_i_dt      = legendre_values2_at_point;
          for (int edgeOrdinal=0; edgeOrdinal<numEdges; edgeOrdinal++)
          {
            const auto & s0 = lambda[edge_start_[edgeOrdinal]];
            const auto & s1 = lambda[  edge_end_[edgeOrdinal]];
            
            const auto & s0_dx = lambda_dx[edge_start_[edgeOrdinal]];
            const auto & s0_dy = lambda_dy[edge_start_[edgeOrdinal]];
            const auto & s0_dz = lambda_dz[edge_start_[edgeOrdinal]];
            const auto & s1_dx = lambda_dx[  edge_end_[edgeOrdinal]];
            const auto & s1_dy = lambda_dy[  edge_end_[edgeOrdinal]];
            const auto & s1_dz = lambda_dz[  edge_end_[edgeOrdinal]];
            
            Polynomials::shiftedScaledLegendreValues             (P_i_minus_1, polyOrder_-1, PointScalar(s1), PointScalar(s0+s1));
            Polynomials::shiftedScaledIntegratedLegendreValues_dt(L_i_dt,      polyOrder_,   PointScalar(s1), PointScalar(s0+s1));
            for (int edgeFunctionOrdinal=0; edgeFunctionOrdinal<num1DEdgeFunctions; edgeFunctionOrdinal++)
            {
              // the first two (integrated) Legendre functions are essentially the vertex functions; hence the +2 here:
              const int i = edgeFunctionOrdinal+2;
              output_(edgeFunctionOrdinal+fieldOrdinalOffset,pointOrdinal,0) = P_i_minus_1(i-1) * s1_dx + L_i_dt(i) * (s1_dx + s0_dx);
              output_(edgeFunctionOrdinal+fieldOrdinalOffset,pointOrdinal,1) = P_i_minus_1(i-1) * s1_dy + L_i_dt(i) * (s1_dy + s0_dy);
              output_(edgeFunctionOrdinal+fieldOrdinalOffset,pointOrdinal,2) = P_i_minus_1(i-1) * s1_dz + L_i_dt(i) * (s1_dz + s0_dz);
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
          auto & L_i            = legendre_values2_at_point;
          auto & L_2i_j_dt      = jacobi_values1_at_point;
          auto & L_2i_j         = jacobi_values2_at_point;
          auto & P_2i_j_minus_1 = jacobi_values3_at_point;
          
          for (int faceOrdinal=0; faceOrdinal<numFaces; faceOrdinal++)
          {
            const auto & s0 = lambda[face_vertex_0[faceOrdinal]];
            const auto & s1 = lambda[face_vertex_1[faceOrdinal]];
            const auto & s2 = lambda[face_vertex_2[faceOrdinal]];
            Polynomials::shiftedScaledIntegratedLegendreValues(L_i, polyOrder_, s1, s0+s1);
            
            const PointScalar jacobiScaling = s0 + s1 + s2;
            
            // compute integrated Jacobi values for each desired value of alpha
            for (int n=2; n<=polyOrder_; n++)
            {
              const double alpha = n*2;
              const int alphaOrdinal = n-2;
              using Kokkos::subview;
              using Kokkos::ALL;
              auto L_2i_j_dt_alpha      = subview(L_2i_j_dt,      alphaOrdinal, ALL);
              auto L_2i_j_alpha         = subview(L_2i_j,         alphaOrdinal, ALL);
              auto P_2i_j_minus_1_alpha = subview(P_2i_j_minus_1, alphaOrdinal, ALL);
              Polynomials::shiftedScaledIntegratedJacobiValues_dt(L_2i_j_dt_alpha,      alpha, polyOrder_-2, s2, jacobiScaling);
              Polynomials::shiftedScaledIntegratedJacobiValues   (L_2i_j_alpha,         alpha, polyOrder_-2, s2, jacobiScaling);
              Polynomials::shiftedScaledJacobiValues(P_2i_j_minus_1_alpha, alpha, polyOrder_-1, s2, jacobiScaling);
            }
            
            const int edgeOrdinal = face_ordinal_of_first_edge[faceOrdinal];
            int localFaceOrdinal = 0;
            for (int totalPolyOrder=3; totalPolyOrder<=polyOrder_; totalPolyOrder++)
            {
              for (int i=2; i<totalPolyOrder; i++)
              {
                const int edgeBasisOrdinal = edgeOrdinal*num1DEdgeFunctions + i-2 + numVertices;
                const auto & grad_L_i_dx = output_(edgeBasisOrdinal,pointOrdinal,0);
                const auto & grad_L_i_dy = output_(edgeBasisOrdinal,pointOrdinal,1);
                const auto & grad_L_i_dz = output_(edgeBasisOrdinal,pointOrdinal,2);
                
                const int alphaOrdinal = i-2;
                
                const auto & s0_dx = lambda_dx[face_vertex_0[faceOrdinal]];
                const auto & s0_dy = lambda_dy[face_vertex_0[faceOrdinal]];
                const auto & s0_dz = lambda_dz[face_vertex_0[faceOrdinal]];
                const auto & s1_dx = lambda_dx[face_vertex_1[faceOrdinal]];
                const auto & s1_dy = lambda_dy[face_vertex_1[faceOrdinal]];
                const auto & s1_dz = lambda_dz[face_vertex_1[faceOrdinal]];
                const auto & s2_dx = lambda_dx[face_vertex_2[faceOrdinal]];
                const auto & s2_dy = lambda_dy[face_vertex_2[faceOrdinal]];
                const auto & s2_dz = lambda_dz[face_vertex_2[faceOrdinal]];
                
                int j = totalPolyOrder - i;
 
                // put references to the entries of interest in like-named variables with lowercase first letters
                auto & l_2i_j         = L_2i_j(alphaOrdinal,j);
                auto & l_i            = L_i(i);
                auto & l_2i_j_dt      = L_2i_j_dt(alphaOrdinal,j);
                auto & p_2i_j_minus_1 = P_2i_j_minus_1(alphaOrdinal,j-1);
                
                const OutputScalar basisValue_dx = l_2i_j * grad_L_i_dx + l_i * (p_2i_j_minus_1 * s2_dx + l_2i_j_dt * (s0_dx + s1_dx + s2_dx));
                const OutputScalar basisValue_dy = l_2i_j * grad_L_i_dy + l_i * (p_2i_j_minus_1 * s2_dy + l_2i_j_dt * (s0_dy + s1_dy + s2_dy));
                const OutputScalar basisValue_dz = l_2i_j * grad_L_i_dz + l_i * (p_2i_j_minus_1 * s2_dz + l_2i_j_dt * (s0_dz + s1_dz + s2_dz));
                
                const int fieldOrdinal = fieldOrdinalOffset + localFaceOrdinal;
                
                output_(fieldOrdinal,pointOrdinal,0) = basisValue_dx;
                output_(fieldOrdinal,pointOrdinal,1) = basisValue_dy;
                output_(fieldOrdinal,pointOrdinal,2) = basisValue_dz;
                
                localFaceOrdinal++;
              }
            }
            fieldOrdinalOffset += numFaceBasisFunctionsPerFace;
          }
          // interior functions
          /*
           Per Fuentes et al. (see Appendix E.1, E.2), the interior functions, defined for i ≥ 2, are
             phi_ij(lambda_012) [L^{2(i+j)}_k](1-lambda_3,lambda_3) = phi_ij(lambda_012) L^{2(i+j)}_k (lambda_3; 1)
           and have gradients:
             L^{2(i+j)}_k (lambda_3; 1) grad (phi_ij(lambda_012)) + phi_ij(lambda_012) grad (L^{2(i+j)}_k (lambda_3; 1))
           where:
             - phi_ij(lambda_012) is the (i,j) basis function on face 012,
             - L^{alpha}_j(t0; t1) is the jth integrated Jacobi polynomial
           and the gradient of the integrated Jacobi polynomial can be computed:
           - grad L^{alpha}_j(t0; t1) = P^{alpha}_{j-1} (t0;t1) grad t0 + R^{alpha}_{j-1}(t0,t1) grad t1
           Here, t1=1, so this simplifies to:
           - grad L^{alpha}_j(t0; t1) = P^{alpha}_{j-1} (t0;t1) grad t0
           
           The P_i we have implemented in shiftedScaledJacobiValues.
           */
          // rename the scratch memory to match our usage here:
          auto & L_alpha = jacobi_values1_at_point;
          auto & P_alpha = jacobi_values2_at_point;
          
          // precompute values used in face ordinal 0:
          {
            const auto & s0 = lambda[0];
            const auto & s1 = lambda[1];
            const auto & s2 = lambda[2];
            // Legendre:
            Polynomials::shiftedScaledIntegratedLegendreValues(legendre_values1_at_point, polyOrder_, PointScalar(s1), PointScalar(s0+s1));
            
            // Jacobi for each desired alpha value:
            const PointScalar jacobiScaling = s0 + s1 + s2;
            for (int n=2; n<=polyOrder_; n++)
            {
              const double alpha = n*2;
              const int alphaOrdinal = n-2;
              using Kokkos::subview;
              using Kokkos::ALL;
              auto jacobi_alpha = subview(jacobi_values3_at_point, alphaOrdinal, ALL);
              Polynomials::shiftedScaledIntegratedJacobiValues(jacobi_alpha, alpha, polyOrder_-2, s2, jacobiScaling);
            }
          }
          
          // interior
          for (int n=3; n<=polyOrder_; n++)
          {
            const double alpha = n*2;
            const double jacobiScaling = 1.0;
            const int alphaOrdinal = n-3;
            using Kokkos::subview;
            using Kokkos::ALL;
            
            // values for interior functions:
            auto L = subview(L_alpha, alphaOrdinal, ALL);
            auto P = subview(P_alpha, alphaOrdinal, ALL);
            Polynomials::shiftedScaledIntegratedJacobiValues(L, alpha, polyOrder_-3, lambda[3], jacobiScaling);
            Polynomials::shiftedScaledJacobiValues          (P, alpha, polyOrder_-3, lambda[3], jacobiScaling);
          }
          
          const int min_i  = 2;
          const int min_j  = 1;
          const int min_k  = 1;
          const int min_ij = min_i + min_j;
          const int min_ijk = min_ij + min_k;
          int localInteriorBasisOrdinal = 0;
          for (int totalPolyOrder_ijk=min_ijk; totalPolyOrder_ijk <= polyOrder_; totalPolyOrder_ijk++)
          {
            int localFaceBasisOrdinal = 0;
            for (int totalPolyOrder_ij=min_ij; totalPolyOrder_ij <= totalPolyOrder_ijk-min_j; totalPolyOrder_ij++)
            {
              for (int i=2; i <= totalPolyOrder_ij-min_j; i++)
              {
                const int j = totalPolyOrder_ij - i;
                const int k = totalPolyOrder_ijk - totalPolyOrder_ij;
                // interior functions use basis values belonging to the first face, 012
                const int faceBasisOrdinal = numEdges*num1DEdgeFunctions + numVertices + localFaceBasisOrdinal;
                
                const auto & faceValue_dx = output_(faceBasisOrdinal,pointOrdinal,0);
                const auto & faceValue_dy = output_(faceBasisOrdinal,pointOrdinal,1);
                const auto & faceValue_dz = output_(faceBasisOrdinal,pointOrdinal,2);
                
                // determine faceValue (on face 0)
                OutputScalar faceValue;
                {
                  const auto & edgeValue = legendre_values1_at_point(i);
                  const int alphaOrdinal = i-2;
                  const auto & jacobiValue = jacobi_values3_at_point(alphaOrdinal,j);
                  faceValue = edgeValue * jacobiValue;
                }
                localFaceBasisOrdinal++;
              
                const int alphaOrdinal = (i+j)-3;
              
                const int fieldOrdinal = fieldOrdinalOffset + localInteriorBasisOrdinal;
                const auto & integratedJacobiValue = L_alpha(alphaOrdinal,k);
                const auto & jacobiValue = P_alpha(alphaOrdinal,k-1);
                output_(fieldOrdinal,pointOrdinal,0) = integratedJacobiValue * faceValue_dx + faceValue * jacobiValue * lambda_dx[3];
                output_(fieldOrdinal,pointOrdinal,1) = integratedJacobiValue * faceValue_dy + faceValue * jacobiValue * lambda_dy[3];
                output_(fieldOrdinal,pointOrdinal,2) = integratedJacobiValue * faceValue_dz + faceValue * jacobiValue * lambda_dz[3];
                
                localInteriorBasisOrdinal++;
              } // end i loop
            } // end totalPolyOrder_ij loop
          } // end totalPolyOrder_ijk loop
          fieldOrdinalOffset += numInteriorBasisFunctions;
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
                                   ">>> ERROR: (Intrepid2::Basis_HGRAD_TET_Cn_FEM_ORTH::OrthPolynomialTri) Computing of second and higher-order derivatives is not currently supported");
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
      // for the (integrated) Legendre computations, we just need p+1 values stored
      // for the (integrated) Jacobi computations, though, we want (p+1)*(# alpha values)
      // alpha is either 2i or 2(i+j), where i=2,…,p or i+j=3,…,p.  So there are at most (p-1) alpha values needed.
      // We can have up to 3 of the (integrated) Jacobi values needed at once.
      const int numAlphaValues = std::max(polyOrder_-1, 1); // make it at least 1 so we can avoid zero-extent ranks…
      size_t shmem_size = 0;
      if (fad_size_output_ > 0)
      {
        // Legendre:
        shmem_size += 2 * OutputScratchView::shmem_size(polyOrder_ + 1, fad_size_output_);
        // Jacobi:
        shmem_size += 3 * OutputScratchView2D::shmem_size(numAlphaValues, polyOrder_ + 1, fad_size_output_);
      }
      else
      {
        // Legendre:
        shmem_size += 2 * OutputScratchView::shmem_size(polyOrder_ + 1);
        // Jacobi:
        shmem_size += 3 * OutputScratchView2D::shmem_size(numAlphaValues, polyOrder_ + 1);
      }
      
      return shmem_size;
    }
  };
  
  /** \class  Intrepid2::IntegratedLegendreBasis_HGRAD_TET
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
               is true, then the first basis function will instead be 1.0-x-y, and the basis will be suitable for
               continuous discretizations.
  */
  template<typename DeviceType,
           typename OutputScalar = double,
           typename PointScalar  = double,
           bool defineVertexFunctions = true>  // if defineVertexFunctions is true, first four basis functions are 1-x-y-z, x, y, and z.  Otherwise, they are 1, x, y, and z.
  class IntegratedLegendreBasis_HGRAD_TET
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
     
     The basis will have polyOrder + 1 members.
     
     If defineVertexFunctions is false, then all basis functions are identified with the interior of the line element, and the first four basis functions are 1, x, y, and z.
     
     If defineVertexFunctions is true, then the first two basis functions are 1-x-y-z, x, y, and z, and these are identified with the left and right vertices of the cell.
     
     */
    IntegratedLegendreBasis_HGRAD_TET(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    polyOrder_(polyOrder),
    pointType_(pointType)
    {
      INTREPID2_TEST_FOR_EXCEPTION(pointType!=POINTTYPE_DEFAULT,std::invalid_argument,"PointType not supported");
      this->basisCardinality_     = ((polyOrder+1) * (polyOrder+2) * (polyOrder+3)) / 6;
      this->basisDegree_          = polyOrder;
      this->basisCellTopologyKey_ = shards::Tetrahedron<>::key;
      this->basisType_            = BASIS_FEM_HIERARCHICAL;
      this->basisCoordinates_     = COORDINATES_CARTESIAN;
      this->functionSpace_        = FUNCTION_SPACE_HGRAD;
      
      const int degreeLength = 1;
      this->fieldOrdinalPolynomialDegree_ = OrdinalTypeArray2DHost("Integrated Legendre H(grad) tetrahedron polynomial degree lookup", this->basisCardinality_, degreeLength);
      this->fieldOrdinalH1PolynomialDegree_ = OrdinalTypeArray2DHost("Integrated Legendre H(grad) tetrahedron polynomial H^1 degree lookup", this->basisCardinality_, degreeLength);
      
      int fieldOrdinalOffset = 0;
      // **** vertex functions **** //
      const int numVertices = this->getBaseCellTopology().getVertexCount();
      const int numFunctionsPerVertex = 1;
      const int numVertexFunctions = numVertices * numFunctionsPerVertex;
      for (int i=0; i<numVertexFunctions; i++)
      {
        // for H(grad) on tetrahedron, if defineVertexFunctions is false, first four basis members are linear
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
      const int numFunctionsPerFace   = ((polyOrder-1)*(polyOrder-2))/2;
      const int numFaces = 4;
      for (int faceOrdinal=0; faceOrdinal<numFaces; faceOrdinal++)
      {
        for (int totalPolyOrder=3; totalPolyOrder<=polyOrder_; totalPolyOrder++)
        {
          const int totalFaceDofs         = (totalPolyOrder-2)*(totalPolyOrder-1)/2;
          const int totalFaceDofsPrevious = (totalPolyOrder-3)*(totalPolyOrder-2)/2;
          const int faceDofsForPolyOrder  = totalFaceDofs - totalFaceDofsPrevious;
          for (int i=0; i<faceDofsForPolyOrder; i++)
          {
            this->fieldOrdinalPolynomialDegree_  (fieldOrdinalOffset,0) = totalPolyOrder;
            this->fieldOrdinalH1PolynomialDegree_(fieldOrdinalOffset,0) = totalPolyOrder;
            fieldOrdinalOffset++;
          }
        }
      }

      // **** interior functions **** //
      const int numFunctionsPerVolume = ((polyOrder-1)*(polyOrder-2)*(polyOrder-3))/6;
      const int numVolumes = 1; // interior
      for (int volumeOrdinal=0; volumeOrdinal<numVolumes; volumeOrdinal++)
      {
        for (int totalPolyOrder=4; totalPolyOrder<=polyOrder_; totalPolyOrder++)
        {
          const int totalInteriorDofs         = (totalPolyOrder-3)*(totalPolyOrder-2)*(totalPolyOrder-1)/6;
          const int totalInteriorDofsPrevious = (totalPolyOrder-4)*(totalPolyOrder-3)*(totalPolyOrder-2)/6;
          const int interiorDofsForPolyOrder  = totalInteriorDofs - totalInteriorDofsPrevious;
          
          for (int i=0; i<interiorDofsForPolyOrder; i++)
          {
            this->fieldOrdinalPolynomialDegree_  (fieldOrdinalOffset,0) = totalPolyOrder;
            this->fieldOrdinalH1PolynomialDegree_(fieldOrdinalOffset,0) = totalPolyOrder;
            fieldOrdinalOffset++;
          }
        }
      }
      
      INTREPID2_TEST_FOR_EXCEPTION(fieldOrdinalOffset != this->basisCardinality_, std::invalid_argument, "Internal error: basis enumeration is incorrect");
      
      // initialize tags
      {
        // ESEAS numbers tetrahedron faces differently from Intrepid2
        // ESEAS:     012, 013, 123, 023
        // Intrepid2: 013, 123, 032, 021
        const int intrepid2FaceOrdinals[4] {3,0,1,2}; // index is the ESEAS face ordinal; value is the intrepid2 ordinal
        
        const auto & cardinality = this->basisCardinality_;
        
        // Basis-dependent initializations
        const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
        const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
        const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
        const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
        
        OrdinalTypeArray1DHost tagView("tag view", cardinality*tagSize);
        const int vertexDim = 0, edgeDim = 1, faceDim = 2, volumeDim = 3;

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
            for (int faceOrdinalESEAS=0; faceOrdinalESEAS<numFaces; faceOrdinalESEAS++)
            {
              int faceOrdinalIntrepid2 = intrepid2FaceOrdinals[faceOrdinalESEAS];
              for (int functionOrdinal=0; functionOrdinal<numFunctionsPerFace; functionOrdinal++)
              {
                tagView(tagNumber*tagSize+0) = faceDim;               // face dimension
                tagView(tagNumber*tagSize+1) = faceOrdinalIntrepid2;  // face id
                tagView(tagNumber*tagSize+2) = functionOrdinal;       // local dof id
                tagView(tagNumber*tagSize+3) = numFunctionsPerFace;   // total number of dofs on this face
                tagNumber++;
              }
            }
            for (int volumeOrdinal=0; volumeOrdinal<numVolumes; volumeOrdinal++)
            {
              for (int functionOrdinal=0; functionOrdinal<numFunctionsPerVolume; functionOrdinal++)
              {
                tagView(tagNumber*tagSize+0) = volumeDim;               // volume dimension
                tagView(tagNumber*tagSize+1) = volumeOrdinal;           // volume id
                tagView(tagNumber*tagSize+2) = functionOrdinal;         // local dof id
                tagView(tagNumber*tagSize+3) = numFunctionsPerVolume;   // total number of dofs in this volume
                tagNumber++;
              }
            }
            INTREPID2_TEST_FOR_EXCEPTION(tagNumber != this->basisCardinality_, std::invalid_argument, "Internal error: basis tag enumeration is incorrect");
          }
        } else {
          for (ordinal_type i=0;i<cardinality;++i) {
            tagView(i*tagSize+0) = volumeDim;   // volume dimension
            tagView(i*tagSize+1) = 0;           // volume ordinal
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
      return "Intrepid2_IntegratedLegendreBasis_HGRAD_TET";
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
    virtual void getValues( OutputViewType outputValues, const PointViewType  inputPoints,
                           const EOperator operatorType = OPERATOR_VALUE ) const override
    {
      auto numPoints = inputPoints.extent_int(0);
      
      using FunctorType = Hierarchical_HGRAD_TET_Functor<DeviceType, OutputScalar, PointScalar, OutputViewType, PointViewType>;
      
      FunctorType functor(operatorType, outputValues, inputPoints, polyOrder_, defineVertexFunctions);
      
      const int outputVectorSize = getVectorSizeForHierarchicalParallelism<OutputScalar>();
      const int pointVectorSize  = getVectorSizeForHierarchicalParallelism<PointScalar>();
      const int vectorSize = std::max(outputVectorSize,pointVectorSize);
      const int teamSize = 1; // because of the way the basis functions are computed, we don't have a second level of parallelism...

      auto policy = Kokkos::TeamPolicy<ExecutionSpace>(numPoints,teamSize,vectorSize);
      Kokkos::parallel_for("Hierarchical_HGRAD_TET_Functor", policy, functor);
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
      } else if(subCellDim == 2) {
        return Teuchos::rcp(new
            IntegratedLegendreBasis_HGRAD_TRI<DeviceType,OutputScalar,PointScalar>
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
      using HostBasisType  = IntegratedLegendreBasis_HGRAD_TET<HostDeviceType, OutputScalar, PointScalar, defineVertexFunctions>;
      return Teuchos::rcp( new HostBasisType(polyOrder_, pointType_) );
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_IntegratedLegendreBasis_HGRAD_TET_h */
