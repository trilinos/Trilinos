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
// Questions? Contact Mauro Perego  (mperego@sandia.gov) or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid2_IntegratedLegendreBasis_HGRAD_PYR.hpp
    \brief  H(grad) basis on the pyramid based on integrated Legendre polynomials.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_IntegratedLegendreBasis_HGRAD_PYR_h
#define Intrepid2_IntegratedLegendreBasis_HGRAD_PYR_h

#include <Kokkos_DynRankView.hpp>

#include <Intrepid2_config.h>

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_DerivedBasis_HGRAD_QUAD.hpp"
#include "Intrepid2_IntegratedLegendreBasis_HGRAD_LINE.hpp"
#include "Intrepid2_IntegratedLegendreBasis_HGRAD_TRI.hpp"
#include "Intrepid2_Polynomials.hpp"
#include "Intrepid2_Utils.hpp"

#include "Teuchos_RCP.hpp"

namespace Intrepid2
{
  /** \class  Intrepid2::Hierarchical_HGRAD_PYR_Functor
      \brief  Functor for computing values for the IntegratedLegendreBasis_HGRAD_PYR class.
   
   This functor is not intended for use outside of IntegratedLegendreBasis_HGRAD_PYR.
  */
  template<class DeviceType, class OutputScalar, class PointScalar,
           class OutputFieldType, class InputPointsType>
  struct Hierarchical_HGRAD_PYR_Functor
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
    
    Hierarchical_HGRAD_PYR_Functor(EOperator opType, OutputFieldType output, InputPointsType inputPoints,
                                   int polyOrder, bool defineVertexFunctions)
    : opType_(opType), output_(output), inputPoints_(inputPoints),
      polyOrder_(polyOrder), defineVertexFunctions_(defineVertexFunctions),
      fad_size_output_(getScalarDimensionForView(output))
    {
      numFields_ = output.extent_int(0);
      numPoints_ = output.extent_int(1);
      const auto & p = polyOrder;
      const auto p3  = p * p * p;
      INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != inputPoints.extent_int(0), std::invalid_argument, "point counts need to match!");
      INTREPID2_TEST_FOR_EXCEPTION(numFields_ != p3 + 3 * p + 1, std::invalid_argument, "output field size does not match basis cardinality");
    }
    
    KOKKOS_INLINE_FUNCTION
    void operator()( const TeamMember & teamMember ) const
    {
      auto pointOrdinal = teamMember.league_rank();
      OutputScratchView scratch1D_1, scratch1D_2, scratch1D_3;
      OutputScratchView2D scratch2D_1, scratch2D_2, scratch2D_3;
      const int numAlphaValues = (polyOrder_-1 > 1) ? (polyOrder_-1) : 1; // make numAlphaValues at least 1 so we can avoid zero-extent allocations…
      if (fad_size_output_ > 0) {
        scratch1D_1 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1, fad_size_output_);
        scratch1D_2 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1, fad_size_output_);
        scratch1D_3 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1, fad_size_output_);
        scratch2D_1 = OutputScratchView2D(teamMember.team_shmem(), numAlphaValues, polyOrder_ + 1, fad_size_output_);
        scratch2D_2 = OutputScratchView2D(teamMember.team_shmem(), numAlphaValues, polyOrder_ + 1, fad_size_output_);
        scratch2D_3 = OutputScratchView2D(teamMember.team_shmem(), numAlphaValues, polyOrder_ + 1, fad_size_output_);
      }
      else {
        scratch1D_1 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1);
        scratch1D_2 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1);
        scratch1D_3 = OutputScratchView(teamMember.team_shmem(), polyOrder_ + 1);
        scratch2D_1 = OutputScratchView2D(teamMember.team_shmem(), numAlphaValues, polyOrder_ + 1);
        scratch2D_2 = OutputScratchView2D(teamMember.team_shmem(), numAlphaValues, polyOrder_ + 1);
        scratch2D_3 = OutputScratchView2D(teamMember.team_shmem(), numAlphaValues, polyOrder_ + 1);
      }
      
      const auto & x = inputPoints_(pointOrdinal,0);
      const auto & y = inputPoints_(pointOrdinal,1);
      const auto & z = inputPoints_(pointOrdinal,2);
      
      // Intrepid2 uses (-1,1)^2 for x,y
      // ESEAS uses (0,1)^2
      // TODO: below, scale x,y derivatives appropriately
      // (Can look at what we do on the HGRAD_LINE for reference; there's a similar difference for line topology.)
      
      Kokkos::Array<OutputScalar,3> coords {(x+1.)/2.,(y+1.)/2.,z}; // map x,y coordinates from (-1,1)^2 to (0,1)^2
      
      // pyramid "affine" coordinates and gradients get stored in lambda, lambdaGrad:
      using Kokkos::Array;
      Array<OutputScalar,5> lambda;
      Array<Kokkos::Array<OutputScalar,3>,5> lambdaGrad;
//      computeLambda    (    lambda, coords);
//      computeLambdaGrad(lambdaGrad, coords);
      
      Array<Array<OutputScalar,3>,2> mu;
      Array<Array<Array<OutputScalar,3>,3>,2> muGrad;
//      computeMu    (    mu, coords);
//      computeMuGrad(muGrad, coords);
      
      Array<Array<OutputScalar,2>,3> nu;
      Array<Array<Array<OutputScalar,3>,2>,3> nuGrad;
//      computeNu    (    nu, coords);
//      computeNuGrad(nuGrad, coords);
      
      affinePyramid(lambda, lambdaGrad, mu, muGrad, nu, nuGrad, coords);
      
      switch (opType_)
      {
        case OPERATOR_VALUE:
        {
          // vertex functions come first, according to vertex ordering: (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1)
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
          
          // rename scratch1, scratch2
          auto & Li_a1 = scratch1D_1;
          auto & Li_a2 = scratch1D_2;
          
          Polynomials::shiftedScaledIntegratedLegendreValues(Li_a1, polyOrder_, nu[0][0], nu[0][0] + nu[1][0]);
          Polynomials::shiftedScaledIntegratedLegendreValues(Li_a2, polyOrder_, nu[0][1], nu[0][1] + nu[1][1]);
          
          // edge functions
          // "mixed" edges (those shared between the quad and some tri face) first
          int fieldOrdinalOffset = numVertices;
          for (int edgeOrdinal=0; edgeOrdinal<numMixedEdges; edgeOrdinal++)
          {
            // edge 0,2 --> a=1, b=2
            // edge 1,3 --> a=2, b=1
            int a = (edgeOrdinal % 2 == 0) ? 1 : 2;
            int b = 3 - a;
            auto & Li = (a == 1) ? Li_a1 : Li_a2;
            // edge 0,3 --> c=0
            // edge 1,2 --> c=1
            int c = ((edgeOrdinal == 0) || (edgeOrdinal == 3)) ? 0 : 1;
            for (int i=2; i<=polyOrder_; i++)
            {
              output_(fieldOrdinalOffset,pointOrdinal) = mu[c][b-1] * Li(i);
              fieldOrdinalOffset++;
            }
          }
          
          // triangle edges next
          // rename scratch1
          auto & Li = scratch1D_1;
          for (int edgeOrdinal=0; edgeOrdinal<numMixedEdges; edgeOrdinal++)
          {
            Polynomials::shiftedScaledIntegratedLegendreValues(Li, polyOrder_, lambda[0], lambda[0] + lambda[4]);
            for (int i=2; i<=polyOrder_; i++)
            {
              output_(fieldOrdinalOffset,pointOrdinal) = Li(i);
              fieldOrdinalOffset++;
            }
          }
                    
          // quadrilateral face
          // mu_0 * phi_i(mu_0^{xi_1},mu_1^{xi_1}) * phi_j(mu_0^{xi_2},mu_1^{xi_2})
          
          // rename scratch
                 Li = scratch1D_1;
          auto & Lj = scratch1D_2;
          Polynomials::shiftedScaledIntegratedLegendreValues(Li, polyOrder_, mu[0][0], mu[0][0] + mu[1][0]);
          Polynomials::shiftedScaledIntegratedLegendreValues(Lj, polyOrder_, mu[0][1], mu[0][1] + mu[1][1]);
          // following the ESEAS ordering: j increments first
          for (int j=2; j<=polyOrder_; j++)
          {
            for (int i=2; i<=polyOrder_; i++)
            {
              output_(fieldOrdinalOffset,pointOrdinal) = mu[0][2] * Li(i) * Lj(j);
              fieldOrdinalOffset++;
            }
          }
          
          /*
           Face functions for triangular face (d,e,f) are the product of:
           mu_c^{zeta,xi_b},
           edge functions on their (d,e) edge,
           and a Jacobi polynomial [L^2i_j](s0+s1,s2) = L^2i_j(s2;s0+s1+s2),
           
           where s0, s1, s2 are nu[0][a-1],nu[1][a-1],nu[2][a-1],
           and (a,b) = (1,2) for faces 0,2; (a,b) = (2,1) for others;
                   c = 0     for faces 0,3;     c = 1 for others
           */
          for (int faceOrdinal=0; faceOrdinal<numTriFaces; faceOrdinal++)
          {
            // face 0,2 --> a=1, b=2
            // face 1,3 --> a=2, b=1
            int a = (faceOrdinal % 2 == 0) ? 1 : 2;
            int b = 3 - a;
            // face 0,3 --> c=0
            // face 1,2 --> c=1
            int c = ((faceOrdinal == 0) || (faceOrdinal == 3)) ? 0 : 1;
            
            const auto & s0 = nu[0][a-1];
            const auto & s1 = nu[1][a-1];
            const auto & s2 = nu[2][a-1];
            const PointScalar jacobiScaling = s0 + s1 + s2;
            
            // compute integrated Jacobi values for each desired value of alpha
            // relabel scratch2D_1
            auto & jacobi = scratch2D_1;
            for (int n=2; n<=polyOrder_; n++)
            {
              const double alpha = n*2;
              const int alphaOrdinal = n-2;
              using Kokkos::subview;
              using Kokkos::ALL;
              auto jacobi_alpha = subview(jacobi, alphaOrdinal, ALL);
              Polynomials::shiftedScaledIntegratedJacobiValues(jacobi_alpha, alpha, polyOrder_-2, s2, jacobiScaling);
            }
            // relabel scratch1D_1
            auto & edge_s0s1 = scratch1D_1;
            Polynomials::shiftedScaledIntegratedLegendreValues(edge_s0s1, polyOrder_, s0, s0 + s1);
            
            for (int totalPolyOrder=3; totalPolyOrder<=polyOrder_; totalPolyOrder++)
            {
              for (int i=2; i<totalPolyOrder; i++)
              {
                const auto & edgeValue = edge_s0s1(i);
                const int alphaOrdinal = i-2;
                
                const int j = totalPolyOrder - i;
                const auto & jacobiValue = jacobi(alphaOrdinal,j);
                output_(fieldOrdinalOffset,pointOrdinal) = mu[c][b-1] * edgeValue * jacobiValue;
                
                fieldOrdinalOffset++;
              }
            }
          }
            
          // interior functions
          // these are the product of the same quadrilateral function that we blended for the quadrilateral face and
          // [L_k](mu_{0}^zeta, mu_1^zeta)
            
          // rename scratch
                 Li = scratch1D_1;
                 Lj = scratch1D_2;
          auto & Lk = scratch1D_3;
          Polynomials::shiftedScaledIntegratedLegendreValues(Li, polyOrder_, mu[0][0], mu[0][0] + mu[1][0]);
          Polynomials::shiftedScaledIntegratedLegendreValues(Lj, polyOrder_, mu[0][1], mu[0][1] + mu[1][1]);
          Polynomials::shiftedScaledIntegratedLegendreValues(Lk, polyOrder_, mu[0][2], mu[0][2] + mu[1][2]);
          // following the ESEAS ordering: k increments first
          for (int k=2; k<=polyOrder_; k++)
          {
            for (int j=2; j<=polyOrder_; j++)
            {
              for (int i=2; i<=polyOrder_; i++)
              {
                output_(fieldOrdinalOffset,pointOrdinal) = Lk(k) * Li(i) * Lj(j);
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
          
          for (int vertexOrdinal=0; vertexOrdinal<numVertices; vertexOrdinal++)
          {
            for (int d=0; d<3; d++)
            {
              output_(vertexOrdinal,pointOrdinal,d) = lambdaGrad[vertexOrdinal][d];
            }
          }
          
          if (!defineVertexFunctions_)
          {
            // "DG" basis case
            // here, the first "vertex" function is 1, so the derivative is 0:
            output_(0,pointOrdinal,0) = 0.0;
            output_(0,pointOrdinal,1) = 0.0;
            output_(0,pointOrdinal,2) = 0.0;
          }

          // edge functions
          int fieldOrdinalOffset = numVertices;
          
          // mixed edges first
          auto & P_i_minus_1 = scratch1D_1;
          auto & L_i_dt      = scratch1D_2; // R_{i-1} = d/dt L_i
          auto & L_i         = scratch1D_3;
          
          for (int edgeOrdinal=0; edgeOrdinal<numMixedEdges; edgeOrdinal++)
          {
            // edge 0,2 --> a=1, b=2
            // edge 1,3 --> a=2, b=1
            int a = (edgeOrdinal % 2 == 0) ? 1 : 2;
            int b = 3 - a;
            
            Polynomials::shiftedScaledLegendreValues             (P_i_minus_1, polyOrder_-1, nu[0][a-1], nu[0][a-1] + nu[1][a-1]);
            Polynomials::shiftedScaledIntegratedLegendreValues_dt(L_i_dt,      polyOrder_,   nu[0][a-1], nu[0][a-1] + nu[1][a-1]);
            Polynomials::shiftedScaledIntegratedLegendreValues   (L_i,         polyOrder_,   nu[0][a-1], nu[0][a-1] + nu[1][a-1]);
            
            // edge 0,3 --> c=0
            // edge 1,2 --> c=1
            int c = ((edgeOrdinal == 0) || (edgeOrdinal == 3)) ? 0 : 1;
            for (int i=2; i<=polyOrder_; i++)
            {
              // grad (mu[c][b-1] * Li(i)) = grad (mu[c][b-1]) * Li(i) + mu[c][b-1] * grad Li(i)
              
              for (int d=0; d<3; d++)
              {
                // grad [L_i](nu_0,nu_1) = [P_{i-1}](nu_0,nu_1) * grad nu_1 + [R_{i-1}](nu_0,nu_1) * grad (nu_0 + nu_1)
                const auto & R_i_minus_1 = L_i_dt(i);
                OutputScalar grad_Li_d = P_i_minus_1(i-1) * nuGrad[1][a-1][d] + R_i_minus_1 * (nuGrad[0][a-1][d] + nuGrad[1][a-1][d]);
                output_(fieldOrdinalOffset,pointOrdinal,d) = muGrad[c][b-1][d] * L_i(i) + mu[c][b-1] * grad_Li_d;
              }
              fieldOrdinalOffset++;
            }
          }
          
          // "triangle" edges
          // TODO: implement these
          
          // TODO: quadrilateral face
          // TODO: triangle faces
          // TODO: interior functions
          
          for (int basisOrdinal=0; basisOrdinal<numFields_; basisOrdinal++)
          {
            // scale x, y derivatives by 0.5 to account for the ref space transformation: Intrepid2 uses (-1,1)^2; ESEAS uses (0,1)^2
            output_(basisOrdinal,pointOrdinal,0) *= 0.5;
            output_(basisOrdinal,pointOrdinal,1) *= 0.5;
          }
          
        } // end OPERATOR_GRAD block
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
                                   ">>> ERROR: (Intrepid2::Hierarchical_HGRAD_PYR_Functor) Computing of second and higher-order derivatives is not currently supported");
        default:
          // unsupported operator type
          device_assert(false);
      }
    }
    
    KOKKOS_INLINE_FUNCTION
    void affinePyramid(Kokkos::Array<OutputScalar,5>                                   &lambda,
                       Kokkos::Array<Kokkos::Array<OutputScalar,3>,5>                  &lambdaGrad,
                       Kokkos::Array<Kokkos::Array<OutputScalar,3>,2>                  &mu,
                       Kokkos::Array<Kokkos::Array<Kokkos::Array<OutputScalar,3>,3>,2> &muGrad,
                       Kokkos::Array<Kokkos::Array<OutputScalar,2>,3>                  &nu,
                       Kokkos::Array<Kokkos::Array<Kokkos::Array<OutputScalar,3>,2>,3> &nuGrad,
                       Kokkos::Array<PointScalar,3>  &coords) const
    {
      const auto & x = coords[0];
      const auto & y = coords[1];
      const auto & z = coords[2];
      nu[0][0]  = 1. - x - z; // nu_0^{\zeta,\xi_1}
      nu[0][1]  = 1. - y - z; // nu_0^{\zeta,\xi_2}
      nu[1][0]  =      x    ; // nu_1^{\zeta,\xi_1}
      nu[1][1]  =      y    ; // nu_1^{\zeta,\xi_2}
      nu[2][0]  =      z    ; // nu_2^{\zeta,\xi_1}
      nu[2][1]  =      z    ; // nu_2^{\zeta,\xi_2}
      
      nuGrad[0][0][0]  = -1. ; // nu_0^{\zeta,\xi_1}_dxi_1
      nuGrad[0][0][1]  =  0. ; // nu_0^{\zeta,\xi_1}_dxi_2
      nuGrad[0][0][2]  = -1. ; // nu_0^{\zeta,\xi_1}_dzeta
      
      nuGrad[0][1][0]  =  0. ; // nu_0^{\zeta,\xi_2}_dxi_1
      nuGrad[0][1][1]  = -1. ; // nu_0^{\zeta,\xi_2}_dxi_2
      nuGrad[0][1][2]  = -1. ; // nu_0^{\zeta,\xi_2}_dzeta
      
      nuGrad[1][0][0]  =  1. ; // nu_1^{\zeta,\xi_1}_dxi_1
      nuGrad[1][0][1]  =  0. ; // nu_1^{\zeta,\xi_1}_dxi_2
      nuGrad[1][0][2]  =  0. ; // nu_1^{\zeta,\xi_1}_dzeta
      
      nuGrad[1][1][0]  =  0. ; // nu_1^{\zeta,\xi_2}_dxi_1
      nuGrad[1][1][1]  =  1. ; // nu_1^{\zeta,\xi_2}_dxi_2
      nuGrad[1][1][2]  =  0. ; // nu_1^{\zeta,\xi_2}_dzeta
      
      nuGrad[2][0][0]  =  0. ; // nu_2^{\zeta,\xi_1}_dxi_1
      nuGrad[2][0][1]  =  0. ; // nu_2^{\zeta,\xi_1}_dxi_2
      nuGrad[2][0][2]  =  1. ; // nu_2^{\zeta,\xi_1}_dzeta
      
      nuGrad[2][1][0]  =  0. ; // nu_2^{\zeta,\xi_2}_dxi_1
      nuGrad[2][1][1]  =  0. ; // nu_2^{\zeta,\xi_2}_dxi_2
      nuGrad[2][1][2]  =  1. ; // nu_2^{\zeta,\xi_2}_dzeta
      
      // (1 - z) goes in denominator -- so we check for 1-z=0
      auto & muZ_0 = mu[0][2];
      auto & muZ_1 = mu[1][2];
      const double epsilon = 1e-12;
      muZ_0 = (fabs(1.-z) > epsilon) ? 1. - z :      epsilon;
      muZ_1 = (fabs(1.-z) > epsilon) ?      z : 1. - epsilon;
      PointScalar scaling = 1. / muZ_0;
      mu[0][0]  = 1. - x * scaling;
      mu[0][1]  = 1. - y * scaling;
      mu[1][0]  =      x * scaling;
      mu[1][1]  =      y * scaling;
      
      PointScalar scaling2 = scaling * scaling;
      muGrad[0][0][0]  =  -scaling     ;
      muGrad[0][0][1]  =     0.        ;
      muGrad[0][0][2]  = - x * scaling2;
      
      muGrad[0][1][0]  =     0.       ;
      muGrad[0][1][1]  =  -scaling    ;
      muGrad[0][1][2]  = -y * scaling2;
      
      muGrad[0][2][0]  =     0.   ;
      muGrad[0][2][1]  =     0.   ;
      muGrad[0][2][2]  =    -1.   ;
      
      muGrad[1][0][0]  =  scaling     ;
      muGrad[1][0][1]  =     0.       ;
      muGrad[1][0][2]  =  x * scaling2;
      
      muGrad[1][1][0]  =     0.       ;
      muGrad[1][1][1]  =   scaling    ;
      muGrad[1][1][2]  =  y * scaling2;
      
      muGrad[1][2][0]  =     0.   ;
      muGrad[1][2][1]  =     0.   ;
      muGrad[1][2][2]  =     1.   ;
      
      lambda[0] = nu[0][0] * mu[0][1];
      lambda[1] = nu[0][1] * mu[1][0];
      lambda[2] = nu[1][0] * mu[1][1];
      lambda[3] = nu[1][1] * mu[0][0];
      lambda[4] = z;
      
      for (int d=0; d<3; d++)
      {
        lambdaGrad[0][d] = nu[0][0] * muGrad[0][1][d] + nuGrad[0][0][d] * mu[0][1];
        lambdaGrad[1][d] = nu[0][1] * muGrad[1][0][d] + nuGrad[0][1][d] * mu[1][0];
        lambdaGrad[2][d] = nu[1][0] * muGrad[1][1][d] + nuGrad[1][0][d] * mu[1][1];
        lambdaGrad[3][d] = nu[1][1] * muGrad[0][0][d] + nuGrad[1][1][d] * mu[0][0];
      }
      lambdaGrad[4][0] = 0;
      lambdaGrad[4][1] = 0;
      lambdaGrad[4][2] = 1;
    }
    
    KOKKOS_INLINE_FUNCTION
    void computeLambda(Kokkos::Array<OutputScalar,5> &lambda,
                       Kokkos::Array<PointScalar,3>  &coords) const
    {
      const auto & x = coords[0];
      const auto & y = coords[1];
      const auto & z = coords[2];
      // (1 - z) goes in denominator -- so we check for 1-z=0
      const double epsilon = 1e-12;
      PointScalar one_minus_z = (fabs(1.-z) > epsilon) ? 1. - z : epsilon;
      PointScalar scaling = 1. / one_minus_z;
      PointScalar lambda0_x = 1. - z - x;
      PointScalar lambda0_y = 1. - z - y;
      lambda[0] = lambda0_x * lambda0_y * scaling;
      lambda[1] =         x * lambda0_y * scaling;
      lambda[2] =         x *         y * scaling;
      lambda[3] = lambda0_x *         y * scaling;
      lambda[4] = z;
    }
    
    KOKKOS_INLINE_FUNCTION
    void computeLambdaGrad(Kokkos::Array<Kokkos::Array<OutputScalar,3>,5> &lambdaGrad,
                           Kokkos::Array<PointScalar,3>  &coords) const
    {
      const auto & x = coords[0];
      const auto & y = coords[1];
      const auto & z = coords[2];
      // (1 - z) goes in denominator -- so we check for 1-z=0
      const double epsilon = 1e-12;
      PointScalar one_minus_z = (fabs(1.-z) > epsilon) ? 1. - z : epsilon;
      PointScalar scaling = 1. / one_minus_z;
      
      lambdaGrad[0][0] = -1. + y * scaling;
      lambdaGrad[0][1] = -1. + x * scaling;
      lambdaGrad[0][2] = x * y * scaling * scaling - 1.;
      
      lambdaGrad[1][0] =  1. - y * scaling;
      lambdaGrad[1][1] =     - x * scaling;
      lambdaGrad[1][2] = -x *  y * scaling * scaling;
      
      lambdaGrad[2][0] =       y * scaling;
      lambdaGrad[2][1] =       x * scaling;
      lambdaGrad[2][2] =  x *  y * scaling * scaling;
      
      lambdaGrad[3][0] =     - y * scaling;
      lambdaGrad[3][1] =  1. - x * scaling;
      lambdaGrad[3][2] = -x *  y * scaling * scaling;
      
      lambdaGrad[4][0] = 0;
      lambdaGrad[4][1] = 0;
      lambdaGrad[4][2] = 1;
    }
    
    //! See Fuentes et al, Appendix E.9. (Pyramid), definition of mu^{zeta,\xi_j}_i.  Here, the indices to the mu array are [i][j-1], and [i][2] corresponds to mu^{\zeta}_i.  Below, for clarity we use x,y,z as coordinate labels, rather than xi_1, xi_2, zeta.
    KOKKOS_INLINE_FUNCTION
    void computeMu(Kokkos::Array<Kokkos::Array<OutputScalar,3>,2> &mu,
                   Kokkos::Array<PointScalar,3>  &coords) const
    {
      const auto & x = coords[0];
      const auto & y = coords[1];
      // (1 - z) goes in denominator -- so we check for 1-z=0
      const double epsilon = 1e-12;
      PointScalar one_minus_z = (fabs(1.-coords[2]) > epsilon) ? 1. - coords[2] : epsilon;
      PointScalar z           = (fabs(1.-coords[2]) > epsilon) ?      coords[2] : 1. - epsilon;
      PointScalar scaling = 1. / one_minus_z;
      mu[0][0]  = 1. - x * scaling;
      mu[0][1]  = 1. - y * scaling;
      mu[0][2]  = one_minus_z;
      mu[1][0]  =      x * scaling;
      mu[1][1]  =      y * scaling;
      mu[1][2]  =      z;
    }
    
    //! See Fuentes et al, Appendix E.9. (Pyramid), definition of mu^{zeta,\xi_j}_i.  Here, the indices to the mu array are [i][j-1], and [i][2] corresponds to mu^{\zeta}_i.  Below, for clarity we use x,y,z as coordinate labels, rather than xi_1, xi_2, zeta.
    KOKKOS_INLINE_FUNCTION
    void computeMuGrad(Kokkos::Array<Kokkos::Array<Kokkos::Array<OutputScalar,3>,3>,2> &muGrad,
                       Kokkos::Array<PointScalar,3>  &coords) const
    {
      const auto & x = coords[0];
      const auto & y = coords[1];
      // (1 - z) goes in denominator -- so we check for 1-z=0
      const double epsilon = 1e-12;
      PointScalar one_minus_z = (fabs(1.-coords[2]) > epsilon) ? 1. - coords[2] : epsilon;
      PointScalar scaling  = 1. / one_minus_z;
      PointScalar scaling2 = scaling * scaling;
      muGrad[0][0][0]  =  -scaling     ;
      muGrad[0][0][1]  =     0.        ;
      muGrad[0][0][2]  = - x * scaling2;
      
      muGrad[0][1][0]  =     0.       ;
      muGrad[0][1][1]  =  -scaling    ;
      muGrad[0][1][2]  = -y * scaling2;
      
      muGrad[0][2][0]  =     0.   ;
      muGrad[0][2][1]  =     0.   ;
      muGrad[0][2][2]  =    -1.   ;
      
      muGrad[1][0][0]  =  scaling     ;
      muGrad[1][0][1]  =     0.       ;
      muGrad[1][0][2]  =  x * scaling2;
      
      muGrad[1][1][0]  =     0.       ;
      muGrad[1][1][1]  =   scaling    ;
      muGrad[1][1][2]  =  y * scaling2;
      
      muGrad[1][2][0]  =     0.   ;
      muGrad[1][2][1]  =     0.   ;
      muGrad[1][2][2]  =     1.   ;
    }
    
    //! See Fuentes et al, Appendix E.9. (Pyramid), definition of nu^{zeta,\xi_j}_j.  Here, the indices to the nu array are [i][j-1].  Below, for clarity we use x,y,z as coordinate labels, rather than xi_1, xi_2, zeta.
    KOKKOS_INLINE_FUNCTION
    void computeNu(Kokkos::Array<Kokkos::Array<OutputScalar,2>,3> &nu,
                   Kokkos::Array<PointScalar,3>  &coords) const
    {
      const auto & x = coords[0];
      const auto & y = coords[1];
      const auto & z = coords[2];
      nu[0][0]  = 1. - x - z; // nu_0^{\zeta,\xi_1}
      nu[0][1]  = 1. - y - z; // nu_0^{\zeta,\xi_2}
      nu[1][0]  =      x    ; // nu_1^{\zeta,\xi_1}
      nu[1][1]  =      y    ; // nu_1^{\zeta,\xi_2}
      nu[2][0]  =      z    ; // nu_2^{\zeta,\xi_1}
      nu[2][1]  =      z    ; // nu_2^{\zeta,\xi_2}
    }
    
    //! See Fuentes et al, Appendix E.9. (Pyramid), definition of nu^{zeta,\xi_j}_j.  Here, the indices to the nu array are [i][j-1].  Below, for clarity we use x,y,z as coordinate labels, rather than xi_1, xi_2, zeta.
    KOKKOS_INLINE_FUNCTION
    void computeNuGrad(Kokkos::Array<Kokkos::Array<Kokkos::Array<OutputScalar,3>,2>,3> &nuGrad,
                       Kokkos::Array<PointScalar,3>  &coords) const
    {
      nuGrad[0][0][0]  = -1. ; // nu_0^{\zeta,\xi_1}_dxi_1
      nuGrad[0][0][1]  =  0. ; // nu_0^{\zeta,\xi_1}_dxi_2
      nuGrad[0][0][2]  = -1. ; // nu_0^{\zeta,\xi_1}_dzeta
      
      nuGrad[0][1][0]  =  0. ; // nu_0^{\zeta,\xi_2}_dxi_1
      nuGrad[0][1][1]  = -1. ; // nu_0^{\zeta,\xi_2}_dxi_2
      nuGrad[0][1][2]  = -1. ; // nu_0^{\zeta,\xi_2}_dzeta
      
      nuGrad[1][0][0]  =  1. ; // nu_1^{\zeta,\xi_1}_dxi_1
      nuGrad[1][0][1]  =  0. ; // nu_1^{\zeta,\xi_1}_dxi_2
      nuGrad[1][0][2]  =  0. ; // nu_1^{\zeta,\xi_1}_dzeta
      
      nuGrad[1][1][0]  =  0. ; // nu_1^{\zeta,\xi_2}_dxi_1
      nuGrad[1][1][1]  =  1. ; // nu_1^{\zeta,\xi_2}_dxi_2
      nuGrad[1][1][2]  =  0. ; // nu_1^{\zeta,\xi_2}_dzeta
      
      nuGrad[2][0][0]  =  0. ; // nu_2^{\zeta,\xi_1}_dxi_1
      nuGrad[2][0][1]  =  0. ; // nu_2^{\zeta,\xi_1}_dxi_2
      nuGrad[2][0][2]  =  1. ; // nu_2^{\zeta,\xi_1}_dzeta
      
      nuGrad[2][1][0]  =  0. ; // nu_2^{\zeta,\xi_2}_dxi_1
      nuGrad[2][1][1]  =  0. ; // nu_2^{\zeta,\xi_2}_dxi_2
      nuGrad[2][1][2]  =  1. ; // nu_2^{\zeta,\xi_2}_dzeta
    }
    
    // Provide the shared memory capacity.
    // This function takes the team_size as an argument,
    // which allows team_size-dependent allocations.
    size_t team_shmem_size (int team_size) const
    {
      // we use shared memory to create a fast buffer for basis computations
      // for the (integrated) Legendre computations, we just need p+1 values stored
      // for the (integrated) Jacobi computations, though, we want (p+1)*(# alpha values)
      // alpha is either 2i or 2(i+j), where i=2,…,p or i+j=3,…,p.  So there are at most (p-1) alpha values needed.
      // We can have up to 3 of the (integrated) Jacobi values needed at once.
      const int numAlphaValues = std::max(polyOrder_-1, 1); // make it at least 1 so we can avoid zero-extent ranks…
      size_t shmem_size = 0;
      if (fad_size_output_ > 0)
      {
        // Legendre:
        shmem_size += 3 * OutputScratchView::shmem_size(polyOrder_ + 1, fad_size_output_);
        // Jacobi:
        shmem_size += 3 * OutputScratchView2D::shmem_size(numAlphaValues, polyOrder_ + 1, fad_size_output_);
      }
      else
      {
        // Legendre:
        shmem_size += 3 * OutputScratchView::shmem_size(polyOrder_ + 1);
        // Jacobi:
        shmem_size += 3 * OutputScratchView2D::shmem_size(numAlphaValues, polyOrder_ + 1);
      }
      
      return shmem_size;
    }
  };
  
  /** \class  Intrepid2::IntegratedLegendreBasis_HGRAD_PYR
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
  class IntegratedLegendreBasis_HGRAD_PYR
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
    IntegratedLegendreBasis_HGRAD_PYR(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    polyOrder_(polyOrder),
    pointType_(pointType)
    {
      INTREPID2_TEST_FOR_EXCEPTION(pointType!=POINTTYPE_DEFAULT,std::invalid_argument,"PointType not supported");
      this->basisCardinality_  = polyOrder * polyOrder * polyOrder + 3 * polyOrder + 1;
      this->basisDegree_       = polyOrder;
      this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Pyramid<> >() );
      this->basisType_         = BASIS_FEM_HIERARCHICAL;
      this->basisCoordinates_  = COORDINATES_CARTESIAN;
      this->functionSpace_     = FUNCTION_SPACE_HGRAD;
      
      const int degreeLength = 1;
      this->fieldOrdinalPolynomialDegree_ = OrdinalTypeArray2DHost("Integrated Legendre H(grad) tetrahedron polynomial degree lookup", this->basisCardinality_, degreeLength);
      this->fieldOrdinalH1PolynomialDegree_ = OrdinalTypeArray2DHost("Integrated Legendre H(grad) tetrahedron polynomial H^1 degree lookup", this->basisCardinality_, degreeLength);
      
      int fieldOrdinalOffset = 0;
      // **** vertex functions **** //
      const int numVertices = this->basisCellTopology_.getVertexCount();
      for (int i=0; i<numVertices; i++)
      {
        // for H(grad) on Pyramid, if defineVertexFunctions is false, first five basis members are linear
        // if not, then the only difference is that the first member is constant
        this->fieldOrdinalPolynomialDegree_  (i,0) = 1;
        this->fieldOrdinalH1PolynomialDegree_(i,0) = 1;
      }
      if (!defineVertexFunctions)
      {
        this->fieldOrdinalPolynomialDegree_  (0,0) = 0;
        this->fieldOrdinalH1PolynomialDegree_(0,0) = 0;
      }
      fieldOrdinalOffset += numVertices;
      
      // **** edge functions **** //
      const int numFunctionsPerEdge = polyOrder - 1; // bubble functions: all but the vertices
      const int numEdges            = this->basisCellTopology_.getEdgeCount();
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
      // one quad face
      const int numFunctionsPerQuadFace =  (polyOrder-1)*(polyOrder-1);
      
      // following the ESEAS ordering: j increments first
      for (int j=2; j<=polyOrder_; j++)
      {
        for (int i=2; i<=polyOrder_; i++)
        {
          this->fieldOrdinalPolynomialDegree_  (fieldOrdinalOffset,0) = std::max(i,j);
          this->fieldOrdinalH1PolynomialDegree_(fieldOrdinalOffset,0) = std::max(i,j);
          fieldOrdinalOffset++;
        }
      }
      
      const int numFunctionsPerTriFace = ((polyOrder-1)*(polyOrder-2))/2;
      const int numTriFaces = 4;
      for (int faceOrdinal=0; faceOrdinal<numTriFaces; faceOrdinal++)
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
      const int numFunctionsPerVolume = (polyOrder-1)*(polyOrder-1)*(polyOrder-1);
      const int numVolumes = 1; // interior
      for (int volumeOrdinal=0; volumeOrdinal<numVolumes; volumeOrdinal++)
      {
        // following the ESEAS ordering: k increments first
        for (int k=2; k<=polyOrder_; k++)
        {
          for (int j=2; j<=polyOrder_; j++)
          {
            for (int i=2; i<=polyOrder_; i++)
            {
              const int max_ij  = std::max(i,j);
              const int max_ijk = std::max(max_ij,k);
              this->fieldOrdinalPolynomialDegree_  (fieldOrdinalOffset,0) = max_ijk;
              this->fieldOrdinalH1PolynomialDegree_(fieldOrdinalOffset,0) = max_ijk;
              fieldOrdinalOffset++;
            }
          }
        }
      }
      
      INTREPID2_TEST_FOR_EXCEPTION(fieldOrdinalOffset != this->basisCardinality_, std::invalid_argument, "Internal error: basis enumeration is incorrect");
      
      // initialize tags
      {
        // Intrepid2 vertices:
        // 0: (-1,-1, 0)
        // 1: ( 1,-1, 0)
        // 2: ( 1, 1, 0)
        // 3: (-1, 1, 0)
        // 4: ( 0, 0, 1)
        
        // ESEAS vertices:
        // 0: ( 0, 0, 0)
        // 1: ( 1, 0, 0)
        // 2: ( 1, 1, 0)
        // 3: ( 0, 1, 0)
        // 4: ( 0, 0, 1)
        
        // The edge numbering appears to match between ESEAS and Intrepid2
        
        // ESEAS numbers pyramid faces differently from Intrepid2
        // See BlendProjectPyraTF in ESEAS.
        // See PyramidFaceNodeMap in Shards_BasicTopologies
        // ESEAS:     0123, 014, 124, 324, 034
        // Intrepid2: 014, 124, 234, 304, 0321
        
        const int intrepid2FaceOrdinals[5] {4,0,1,2,3}; // index is the ESEAS face ordinal; value is the intrepid2 ordinal
        
        const auto & cardinality = this->basisCardinality_;
        
        // Basis-dependent initializations
        const ordinal_type tagSize  = 4; // size of DoF tag, i.e., number of fields in the tag
        const ordinal_type posScDim = 0; // position in the tag, counting from 0, of the subcell dim
        const ordinal_type posScOrd = 1; // position in the tag, counting from 0, of the subcell ordinal
        const ordinal_type posDfOrd = 2; // position in the tag, counting from 0, of DoF ordinal relative to the subcell
        
        OrdinalTypeArray1DHost tagView("tag view", cardinality*tagSize);
        const int vertexDim = 0, edgeDim = 1, faceDim = 2, volumeDim = 3;

        if (defineVertexFunctions) {
          {
            int tagNumber = 0;
            for (int vertexOrdinal=0; vertexOrdinal<numVertices; vertexOrdinal++)
            {
              tagView(tagNumber*tagSize+0) = vertexDim;             // vertex dimension
              tagView(tagNumber*tagSize+1) = vertexOrdinal;         // vertex id
              tagView(tagNumber*tagSize+2) = 0;                     // local dof id
              tagView(tagNumber*tagSize+3) = 1;                     // total number of dofs at this vertex
              tagNumber++;
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
            {
              // quad face
              const int faceOrdinalESEAS = 0;
              const int faceOrdinalIntrepid2 = intrepid2FaceOrdinals[faceOrdinalESEAS];
              
              for (int functionOrdinal=0; functionOrdinal<numFunctionsPerQuadFace; functionOrdinal++)
              {
                tagView(tagNumber*tagSize+0) = faceDim;                 // face dimension
                tagView(tagNumber*tagSize+1) = faceOrdinalIntrepid2;    // face id
                tagView(tagNumber*tagSize+2) = functionOrdinal;         // local dof id
                tagView(tagNumber*tagSize+3) = numFunctionsPerQuadFace; // total number of dofs on this face
                tagNumber++;
              }
            }
            for (int triFaceOrdinalESEAS=0; triFaceOrdinalESEAS<numTriFaces; triFaceOrdinalESEAS++)
            {
              const int faceOrdinalESEAS     = triFaceOrdinalESEAS + 1;
              const int faceOrdinalIntrepid2 = intrepid2FaceOrdinals[faceOrdinalESEAS];
              for (int functionOrdinal=0; functionOrdinal<numFunctionsPerTriFace; functionOrdinal++)
              {
                tagView(tagNumber*tagSize+0) = faceDim;                // face dimension
                tagView(tagNumber*tagSize+1) = faceOrdinalIntrepid2;   // face id
                tagView(tagNumber*tagSize+2) = functionOrdinal;        // local dof id
                tagView(tagNumber*tagSize+3) = numFunctionsPerTriFace; // total number of dofs on this face
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
      return "Intrepid2_IntegratedLegendreBasis_HGRAD_PYR";
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
      
      using FunctorType = Hierarchical_HGRAD_PYR_Functor<DeviceType, OutputScalar, PointScalar, OutputViewType, PointViewType>;
      
      FunctorType functor(operatorType, outputValues, inputPoints, polyOrder_, defineVertexFunctions);
      
      const int outputVectorSize = getVectorSizeForHierarchicalParallelism<OutputScalar>();
      const int pointVectorSize  = getVectorSizeForHierarchicalParallelism<PointScalar>();
      const int vectorSize = std::max(outputVectorSize,pointVectorSize);
      const int teamSize = 1; // because of the way the basis functions are computed, we don't have a second level of parallelism...

      auto policy = Kokkos::TeamPolicy<ExecutionSpace>(numPoints,teamSize,vectorSize);
      Kokkos::parallel_for("Hierarchical_HGRAD_PYR_Functor", policy, functor);
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
      using HGRAD_LINE = IntegratedLegendreBasis_HGRAD_LINE<DeviceType,OutputScalar,PointScalar>;
      using HGRAD_TRI  = IntegratedLegendreBasis_HGRAD_TRI<DeviceType,OutputScalar,PointScalar>;
      using HGRAD_QUAD = Basis_Derived_HGRAD_QUAD<HGRAD_LINE>;
      const auto & p = this->basisDegree_;
      if(subCellDim == 1) // line basis
      {
        return Teuchos::rcp(new HGRAD_LINE(p));
      }
      else if (subCellDim == 2)
      {
        if (subCellOrd == 0) // quad basis
        {
          return Teuchos::rcp(new HGRAD_QUAD(p));
        }
        else // tri basis
        {
          return Teuchos::rcp(new HGRAD_TRI(p));
        }
      }
      INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Input parameters out of bounds");
    }

    /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.
     
        \return Pointer to the new Basis object.
     */
    virtual BasisPtr<typename Kokkos::HostSpace::device_type, OutputScalar, PointScalar>
    getHostBasis() const override {
      using HostDeviceType = typename Kokkos::HostSpace::device_type;
      using HostBasisType  = IntegratedLegendreBasis_HGRAD_PYR<HostDeviceType, OutputScalar, PointScalar, defineVertexFunctions>;
      return Teuchos::rcp( new HostBasisType(polyOrder_, pointType_) );
    }
  };
} // end namespace Intrepid2

// do ETI with default (double) type
extern template class Intrepid2::IntegratedLegendreBasis_HGRAD_PYR<Kokkos::DefaultExecutionSpace::device_type,double,double>;

#endif /* Intrepid2_IntegratedLegendreBasis_HGRAD_PYR_h */
