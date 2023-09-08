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

/** \file   Intrepid2_HierarchicalBasis_HDIV_PYR.hpp
    \brief  H(div) basis on the pyramid based on integrated Legendre polynomials.
    \author Created by N.V. Roberts.
 
 Note that although this basis is derived from integrated Legendre polynomials, it is not itself a polynomial basis, but a set of rational functions.
 
 The construction is also hierarchical, in the sense that the basis for p-1 is included in the basis for p.
 */

#ifndef Intrepid2_HierarchicalBasis_HDIV_PYR_h
#define Intrepid2_HierarchicalBasis_HDIV_PYR_h

#include <Kokkos_DynRankView.hpp>

#include <Intrepid2_config.h>

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_DerivedBasis_HVOL_QUAD.hpp"
#include "Intrepid2_LegendreBasis_HVOL_LINE.hpp"
#include "Intrepid2_LegendreBasis_HVOL_TRI.hpp"
#include "Intrepid2_Polynomials.hpp"
#include "Intrepid2_PyramidCoords.hpp"
#include "Intrepid2_Utils.hpp"

#include "Teuchos_RCP.hpp"

namespace Intrepid2
{
  /** \class  Intrepid2::Hierarchical_HDIV_PYR_Functor
      \brief  Functor for computing values for the HierarchicalBasis_HDIV_PYR class.
   
   This functor is not intended for use outside of HierarchicalBasis_HDIV_PYR.
  */
  template<class DeviceType, class OutputScalar, class PointScalar,
           class OutputFieldType, class InputPointsType>
  struct Hierarchical_HDIV_PYR_Functor
  {
    using ExecutionSpace      = typename DeviceType::execution_space;
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
    bool useCGBasis_;
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
    
    Hierarchical_HDIV_PYR_Functor(EOperator opType, OutputFieldType output, InputPointsType inputPoints,
                                   int polyOrder)
    : opType_(opType), output_(output), inputPoints_(inputPoints),
      polyOrder_(polyOrder),
      fad_size_output_(getScalarDimensionForView(output))
    {
      numFields_ = output.extent_int(0);
      numPoints_ = output.extent_int(1);
      const auto & p = polyOrder;
      const auto basisCardinality = p * p + 2 * p * (p+1) + 3 * p * p * (p-1);
      
      INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != inputPoints.extent_int(0), std::invalid_argument, "point counts need to match!");
      INTREPID2_TEST_FOR_EXCEPTION(numFields_ != basisCardinality, std::invalid_argument, "output field size does not match basis cardinality");
    }
    
    //! cross product: c = a x b
    KOKKOS_INLINE_FUNCTION
    void cross(Kokkos::Array<OutputScalar,3> &c,
               const Kokkos::Array<OutputScalar,3> &a,
               const Kokkos::Array<OutputScalar,3> &b) const
    {
      c[0] = a[1] * b[2] - a[2] * b[1];
      c[1] = a[2] * b[0] - a[0] * b[2];
      c[2] = a[0] * b[1] - a[1] * b[0];
    }
    
    //! dot product: c = a \cdot b
    KOKKOS_INLINE_FUNCTION
    void dot(OutputScalar &c,
             const Kokkos::Array<OutputScalar,3> &a,
             const Kokkos::Array<OutputScalar,3> &b) const
    {
      c = 0;
      for (ordinal_type d=0; d<3; d++)
      {
        c += a[d] * b[d];
      }
    }
    
    KOKKOS_INLINE_FUNCTION
    OutputScalar dot(const Kokkos::Array<OutputScalar,3> &a,
                     const Kokkos::Array<OutputScalar,3> &b) const
    {
      OutputScalar c = 0;
      for (ordinal_type d=0; d<3; d++)
      {
        c += a[d] * b[d];
      }
      return c;
    }
    
    KOKKOS_INLINE_FUNCTION
    void E_E(Kokkos::Array<OutputScalar,3> &EE,
             const ordinal_type &i,
             const OutputScratchView &PHom,
             const PointScalar &s0, const PointScalar &s1,
             const Kokkos::Array<PointScalar,3> &s0_grad,
             const Kokkos::Array<PointScalar,3> &s1_grad) const
    {
      const auto & P_i = PHom(i);
      for (ordinal_type d=0; d<3; d++)
      {
        EE[d] = P_i * (s0 * s1_grad[d] - s1 * s0_grad[d]);
      }
    }
    
    KOKKOS_INLINE_FUNCTION
    void E_E_CURL(Kokkos::Array<OutputScalar,3> &curl_EE,
                  const ordinal_type &i,
                  const OutputScratchView &PHom,
                  const PointScalar &s0, const PointScalar &s1,
                  const Kokkos::Array<PointScalar,3> &s0_grad,
                  const Kokkos::Array<PointScalar,3> &s1_grad) const
    {
      // curl (E_i^E)(s0,s1) = (i+2) [P_i](s0,s1) (grad s0 x grad s1)
      OutputScalar ip2_Pi = (i+2) * PHom(i);
      cross(curl_EE, s0_grad, s1_grad);
      for (ordinal_type d=0; d<3; d++)
      {
        curl_EE[d] *= ip2_Pi;
      }
    }
    
    //! The "quadrilateral face" H(div) functions defined by Fuentes et al., Appendix E.2., p. 433
    //! Here, PHom are the homogenized Legendre polynomials [P](s0,s1) and [P](t0,t1), given in Appendix E.1, p. 430
    KOKKOS_INLINE_FUNCTION
    void V_QUAD(Kokkos::Array<OutputScalar,3> &VQUAD,
                const ordinal_type &i, const ordinal_type &j,
                const OutputScratchView &PHom_s,
                const PointScalar &s0, const PointScalar &s1,
                const Kokkos::Array<PointScalar,3> &s0_grad,
                const Kokkos::Array<PointScalar,3> &s1_grad,
                const OutputScratchView &PHom_t,
                const PointScalar &t0, const PointScalar &t1,
                const Kokkos::Array<PointScalar,3> &t0_grad,
                const Kokkos::Array<PointScalar,3> &t1_grad) const
    {
      Kokkos::Array<OutputScalar,3> EE_i, EE_j;
      
      E_E(EE_i, i, PHom_s, s0, s1, s0_grad, s1_grad);
      E_E(EE_j, j, PHom_t, t0, t1, t0_grad, t1_grad);
      
      // VQUAD = EE_i x EE_j:
      cross(VQUAD, EE_i, EE_j);
    }
    
    //! The "quadrilateral face" H(curl) functions defined by Fuentes et al., Appendix E.2., p. 432
    //! Here, HomPi_s01, HomLi_t01 are the homogenized Legendre polynomials [P](s0,s1) and homogenized integrated Legendre polynomials [L](t0,t1), given in Appendix E.1, p. 430
    KOKKOS_INLINE_FUNCTION
    void E_QUAD(Kokkos::Array<OutputScalar,3> &EQUAD,
                const ordinal_type &i, const ordinal_type &j,
                const OutputScratchView &HomPi_s01,
                const PointScalar &s0, const PointScalar &s1,
                const Kokkos::Array<PointScalar,3> &s0_grad,
                const Kokkos::Array<PointScalar,3> &s1_grad,
                const OutputScratchView &HomLi_t01) const
    {
      const OutputScalar &phiE_j = HomLi_t01(j);
      
      Kokkos::Array<OutputScalar,3> EE_i;
      E_E(EE_i, i, HomPi_s01, s0, s1, s0_grad, s1_grad);
      
      for (ordinal_type d=0; d<3; d++)
      {
        EQUAD[d] = phiE_j * EE_i[d];
      }
    }
    
    KOKKOS_INLINE_FUNCTION
    void E_QUAD_CURL(Kokkos::Array<OutputScalar,3> &EQUAD_CURL,
                     const ordinal_type &i, const ordinal_type &j,
                     const OutputScratchView &HomPi_s01,
                     const PointScalar &s0, const PointScalar &s1,
                     const Kokkos::Array<PointScalar,3> &s0_grad,
                     const Kokkos::Array<PointScalar,3> &s1_grad,
                     const OutputScratchView &HomPj_t01,
                     const OutputScratchView &HomLj_t01,
                     const OutputScratchView &HomLj_dt_t01,
                     const Kokkos::Array<PointScalar,3> &t0_grad,
                     const Kokkos::Array<PointScalar,3> &t1_grad) const
    {
      const OutputScalar &phiE_j = HomLj_t01(j);
      
      Kokkos::Array<OutputScalar,3> curl_EE_i;
      E_E_CURL(curl_EE_i, i, HomPi_s01, s0, s1, s0_grad, s1_grad);
      
      Kokkos::Array<OutputScalar,3> EE_i;
      E_E(EE_i, i, HomPi_s01, s0, s1, s0_grad, s1_grad);
      
      Kokkos::Array<OutputScalar,3> grad_phiE_j;
      computeGradHomLi(grad_phiE_j, j, HomPj_t01, HomLj_dt_t01, t0_grad, t1_grad);
      
      cross(EQUAD_CURL, grad_phiE_j, EE_i);
      for (ordinal_type d=0; d<3; d++)
      {
        EQUAD_CURL[d] += phiE_j * curl_EE_i[d];
      }
    }
    
    // This is the "Ancillary Operator" V^{tri}_{ij} on p. 433 of Fuentes et al.
    KOKKOS_INLINE_FUNCTION
    void V_TRI(Kokkos::Array<OutputScalar,3> &VTRI,
               const ordinal_type &i, // i >= 0
               const ordinal_type &j, // j >= 0
               const OutputScratchView &P,      // container in which shiftedScaledLegendreValues have been computed for the appropriate face
               const OutputScratchView &P_2ip1, // container in which shiftedScaledJacobiValues have been computed for (2i+1) for the appropriate face
               const Kokkos::Array<PointScalar,3> &vectorWeight) const // s0 (grad s1 x grad s2) + s1 (grad s2 x grad s0) + s2 (grad s0 x grad s1) -- see computeFaceVectorWeight()
    {
      const auto &P_i      = P(i);
      const auto &P_2ip1_j = P_2ip1(j);
      
      for (ordinal_type d=0; d<3; d++)
      {
        VTRI[d] = P_i * P_2ip1_j * vectorWeight[d];
      }
    }
    
    // Divergence of the "Ancillary Operator" V^{tri}_{ij} on p. 433 of Fuentes et al.
    KOKKOS_INLINE_FUNCTION
    void V_TRI_DIV(OutputScalar &VTRI_DIV,
               const ordinal_type &i, // i >= 0
               const ordinal_type &j, // j >= 0
               const OutputScratchView &P,      // container in which shiftedScaledLegendreValues have been computed for the appropriate face
               const OutputScratchView &P_2ip1, // container in which shiftedScaledJacobiValues have been computed for (2i+1) for the appropriate face
               const OutputScalar &divWeight) const // grad s0 \dot (grad s1 x grad s2)
    {
      const auto &P_i      = P(i);
      const auto &P_2ip1_j = P_2ip1(j);
      
      VTRI_DIV = (i + j + 3.) * P_i * P_2ip1_j * divWeight;
    }
    
    KOKKOS_INLINE_FUNCTION
    void computeFaceVectorWeight(Kokkos::Array<OutputScalar,3> &vectorWeight,
                                 const ordinal_type &a,
                                 Kokkos::Array<Kokkos::Array<PointScalar,2>,3> &nu,
                                 Kokkos::Array<Kokkos::Array<Kokkos::Array<PointScalar,3>,2>,3> &nuGrad) const
    {
      // TODO: rewrite the call(s) to this method to use the other computeFaceVectorWeight()
      // compute s0 (grad s1 x grad s2) + s1 (grad s2 x grad s0) + s2 (grad s0 x grad s1)
      // where s = nu
      
      const auto & s0    = nu    [0][a-1];
      const auto & s0_dx = nuGrad[0][a-1][0];
      const auto & s0_dy = nuGrad[0][a-1][1];
      const auto & s0_dz = nuGrad[0][a-1][2];
      
      const auto & s1    = nu    [1][a-1];
      const auto & s1_dx = nuGrad[1][a-1][0];
      const auto & s1_dy = nuGrad[1][a-1][1];
      const auto & s1_dz = nuGrad[1][a-1][2];
      
      const auto & s2    = nu    [2][a-1];
      const auto & s2_dx = nuGrad[2][a-1][0];
      const auto & s2_dy = nuGrad[2][a-1][1];
      const auto & s2_dz = nuGrad[2][a-1][2];
      
      vectorWeight[0] = s0 * (s1_dy * s2_dz - s1_dz * s2_dy)
                      + s1 * (s2_dy * s0_dz - s2_dz * s0_dy)
                      + s2 * (s0_dy * s1_dz - s0_dz * s1_dy);
      
      vectorWeight[1] = s0 * (s1_dz * s2_dx - s1_dx * s2_dz)
                      + s1 * (s2_dz * s0_dx - s2_dx * s0_dz)
                      + s2 * (s0_dz * s1_dx - s0_dx * s1_dz);
      
      vectorWeight[2] = s0 * (s1_dx * s2_dy - s1_dy * s2_dx)
                      + s1 * (s2_dx * s0_dy - s2_dy * s0_dx)
                      + s2 * (s0_dx * s1_dy - s0_dy * s1_dx);
    }
    
    KOKKOS_INLINE_FUNCTION
    void computeFaceVectorWeight(Kokkos::Array<OutputScalar,3> &vectorWeight,
                                 const PointScalar &s0, const Kokkos::Array<PointScalar,3> &s0Grad,
                                 const PointScalar &s1, const Kokkos::Array<PointScalar,3> &s1Grad,
                                 const PointScalar &s2, const Kokkos::Array<PointScalar,3> &s2Grad) const
    {
      // compute s0 (grad s1 x grad s2) + s1 (grad s2 x grad s0) + s2 (grad s0 x grad s1)
      
      Kokkos::Array<Kokkos::Array<PointScalar,3>,3> cross_products;
      
      cross(cross_products[0], s1Grad, s2Grad);
      cross(cross_products[1], s2Grad, s0Grad);
      cross(cross_products[2], s0Grad, s1Grad);
      
      Kokkos::Array<PointScalar,3> s {s0,s1,s2};
      
      for (ordinal_type d=0; d<3; d++)
      {
        OutputScalar v_d = 0;
        for (ordinal_type i=0; i<3; i++)
        {
          v_d += s[i] * cross_products[i][d];
        }
        vectorWeight[d] = v_d;
      }
    }
    
    KOKKOS_INLINE_FUNCTION
    void computeFaceDivWeight(OutputScalar &divWeight,
                              const Kokkos::Array<PointScalar,3> &s0Grad,
                              const Kokkos::Array<PointScalar,3> &s1Grad,
                              const Kokkos::Array<PointScalar,3> &s2Grad) const
    {
      // grad s0 \dot (grad s1 x grad s2)
      
      Kokkos::Array<PointScalar,3> grad_s1_cross_grad_s2;
      cross(grad_s1_cross_grad_s2, s1Grad, s2Grad);
      
      dot(divWeight, s0Grad, grad_s1_cross_grad_s2);
    }
    
    KOKKOS_INLINE_FUNCTION
    void computeGradHomLi(Kokkos::Array<OutputScalar,3> &HomLi_grad, // grad [L_i](s0,s1)
                          const ordinal_type i,
                          const OutputScratchView &HomPi_s0s1,    // [P_i](s0,s1)
                          const OutputScratchView &HomLi_dt_s0s1, // [d/dt L_i](s0,s1)
                          const Kokkos::Array<PointScalar,3> &s0Grad,
                          const Kokkos::Array<PointScalar,3> &s1Grad) const
    {
//      grad [L_i](s0,s1) = [P_{i-1}](s0,s1) * grad s1 + [R_{i-1}](s0,s1) * grad (s0 + s1)
      const auto & R_i_minus_1 = HomLi_dt_s0s1(i); // d/dt L_i = R_{i-1}
      const auto & P_i_minus_1 = HomPi_s0s1(i-1);
      for (ordinal_type d=0; d<3; d++)
      {
        HomLi_grad[d] = P_i_minus_1 * s1Grad[d] + R_i_minus_1 * (s0Grad[d] + s1Grad[d]);
      }
    }
    
    KOKKOS_INLINE_FUNCTION
    void operator()( const TeamMember & teamMember ) const
    {
      // TODO: rewrite this -- copied from H(grad) implementation
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
      // DEBUGGING
//      const double x = -6./7., y = -6./7., z = 0.0;
//      const double x = 0., y = 0., z = 0.0;
      
      // Intrepid2 uses (-1,1)^2 for x,y
      // ESEAS uses (0,1)^2
      // (Can look at what we do on the HGRAD_LINE for reference; there's a similar difference for line topology.)
      
      Kokkos::Array<PointScalar,3> coords;
      transformToESEASPyramid<>(coords[0], coords[1], coords[2], x, y, z); // map x,y coordinates from (-z,z)^2 to (0,z)^2
      
      // pyramid "affine" coordinates and gradients get stored in lambda, lambdaGrad:
      using Kokkos::Array;
      Array<PointScalar,5> lambda;
      Array<Kokkos::Array<PointScalar,3>,5> lambdaGrad;
      
      Array<Array<PointScalar,3>,2> mu;
      Array<Array<Array<PointScalar,3>,3>,2> muGrad;
      
      Array<Array<PointScalar,2>,3> nu;
      Array<Array<Array<PointScalar,3>,2>,3> nuGrad;
      
      affinePyramid(lambda, lambdaGrad, mu, muGrad, nu, nuGrad, coords);
      
      switch (opType_)
      {
        case OPERATOR_VALUE:
        {
          // starting the H(div) rewrite here
          ordinal_type fieldOrdinalOffset = 0;
          // quadrilateral face
          // rename scratch1, scratch2
          auto & Pi = scratch1D_1;
          auto & Pj = scratch1D_2;
          
          auto & s0      =     mu[0][0], & s1      =     mu[1][0];
          auto & s0_grad = muGrad[0][0], & s1_grad = muGrad[1][0];
          auto & t0      =     mu[0][1], & t1      =     mu[1][1];
          auto & t0_grad = muGrad[0][1], & t1_grad = muGrad[1][1];
          
          Polynomials::shiftedScaledLegendreValues(Pi, polyOrder_-1, s1, s0 + s1);
          Polynomials::shiftedScaledLegendreValues(Pj, polyOrder_-1, t1, t0 + t1);
          
          const auto & muZ_0 = mu[0][2];
          OutputScalar mu0_cubed = muZ_0 * muZ_0 * muZ_0;
          
          // following the ESEAS ordering: j increments first
          for (int j=0; j<polyOrder_; j++)
          {
            for (int i=0; i<polyOrder_; i++)
            {
              Kokkos::Array<OutputScalar,3> VQUAD; // output from V_QUAD
              V_QUAD(VQUAD, i, j,
                     Pi, s0, s1, s0_grad, s1_grad,
                     Pj, t0, t1, t0_grad, t1_grad);
              
              for (ordinal_type d=0; d<3; d++)
              {
                output_(fieldOrdinalOffset,pointOrdinal,d) = mu0_cubed * VQUAD[d];
//                {
//                  using namespace std;
//                  cout << "VQUAD[" << d << "]: " << VQUAD[d] << endl;
//                  cout << "output_(" << fieldOrdinalOffset << "," << pointOrdinal << "," << d <<"): " << output_(fieldOrdinalOffset,pointOrdinal,d) << endl;
//                }
              }
              fieldOrdinalOffset++;
            }
          }
          
          // triangle faces
          
          /*
           Face functions for triangular face (d,e,f) given by:
           1/2 * (    mu_c^{zeta,xi_b} * V_TRI_ij(nu_012^{zeta,xi_a})
                  + 1/mu_c^{zeta,xi_b} * V_TRI_ij(nu_012^{zeta,xi_a} * mu_c^{zeta,xi_b}) )
           but the division by mu_c^{zeta,xi_b} should ideally be avoided in computations; there is
           an alternate expression, not implemented by ESEAS.  We start with the above expression
           so that we can get agreement with ESEAS.  Hopefully after that we can do the alternative,
           and confirm that we maintain agreement with ESEAS for points where ESEAS does not generate
           nans.
           
           …
           where s0, s1, s2 are nu[0][a-1],nu[1][a-1],nu[2][a-1] (nu_012 above)
           and (a,b) = (1,2) for faces 0,2; (a,b) = (2,1) for others;
                   c = 0     for faces 0,3;     c = 1 for others
           */
          const auto & P        = scratch1D_1; // for V_TRI( nu_012)
          const auto & P_2ip1   = scratch1D_2;
          const auto & Pmu      = scratch1D_3; // for V_TRI( mu * nu_012)
          const auto & Pmu_2ip1 = scratch1D_4;
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
            
            const PointScalar legendreScaling = s0 + s1;
            Polynomials::shiftedScaledLegendreValues(P, polyOrder_-1, s1, legendreScaling);
            
            const auto lambda0_index = tri_face_vertex_0[faceOrdinal];
            const auto lambda1_index = tri_face_vertex_1[faceOrdinal];
            const auto lambda2_index = tri_face_vertex_2[faceOrdinal];
            
            const auto & mu_c_b = mu[c][b-1];
//            const PointScalar mu_s0 = mu_c_b * s0;
//            const PointScalar mu_s1 = mu_c_b * s1;
//            const PointScalar mu_s2 = mu_c_b * s2;
            const auto & mu_s0 = lambda[lambda0_index];
            const auto & mu_s1 = lambda[lambda1_index];
            const auto & mu_s2 = lambda[lambda2_index];
            
            {
              // DEBUGGING
              double s0_diff = std::abs(mu_s0 - mu_c_b * s0);
              double s1_diff = std::abs(mu_s1 - mu_c_b * s1);
              const double tol = 1e-14;
              if (s0_diff > tol)
              {
                std::cout << "s0_diff: " << s0_diff << std::endl;
              }
              if (s1_diff > tol)
              {
                std::cout << "s1_diff: " << s1_diff << std::endl;
              }
            }
            
            const PointScalar muJacobiScaling = mu_s0 + mu_s1 + mu_s2;
            
            const PointScalar muLegendreScaling = mu_s0 + mu_s1;
            Polynomials::shiftedScaledLegendreValues(Pmu, polyOrder_-1, mu_s1, muLegendreScaling);
            
            Kokkos::Array<PointScalar, 3> vectorWeight;
            computeFaceVectorWeight(vectorWeight, a, nu, nuGrad);
            
            Kokkos::Array<PointScalar,3> & mu_s0_grad = lambdaGrad[lambda0_index];
            Kokkos::Array<PointScalar,3> & mu_s1_grad = lambdaGrad[lambda1_index];
            Kokkos::Array<PointScalar,3> & mu_s2_grad = lambdaGrad[lambda2_index];
            
            Kokkos::Array<PointScalar, 3> muVectorWeight;
            computeFaceVectorWeight(muVectorWeight, mu_s0, mu_s0_grad, mu_s1, mu_s1_grad, mu_s2, mu_s2_grad);
            
            for (int totalPolyOrder=0; totalPolyOrder<polyOrder_; totalPolyOrder++)
            {
              for (int i=0; i<=totalPolyOrder; i++)
              {
                const int j = totalPolyOrder - i;
                
                const double alpha = i*2.0 + 1;
                Polynomials::shiftedScaledJacobiValues(  P_2ip1, alpha, polyOrder_-1,    s2,   jacobiScaling);
                Polynomials::shiftedScaledJacobiValues(Pmu_2ip1, alpha, polyOrder_-1, mu_s2, muJacobiScaling);
                
                Kokkos::Array<OutputScalar,3> VTRI, VTRI_mu; // output from V_TRI
                
                V_TRI(VTRI,    i, j, P,   P_2ip1,     vectorWeight);
                V_TRI(VTRI_mu, i, j, Pmu, Pmu_2ip1, muVectorWeight);
                
                for (ordinal_type d=0; d<3; d++)
                {
                  output_(fieldOrdinalOffset,pointOrdinal,d) = 0.5 * (VTRI[d] * mu_c_b + VTRI_mu[d] / mu_c_b);
                }

                fieldOrdinalOffset++;
              }
            }
          }
          
          // interior functions
          {
            // label scratch
            const auto & Li_muZ01    = scratch1D_1; // used for phi_k^E values in Family I, II, IV
            const auto & Li_muX01    = scratch1D_2; // used for E_QUAD computations
            const auto & Li_muY01    = scratch1D_3; // used for E_QUAD computations
            const auto & Pi_muX01    = scratch1D_4; // used for E_QUAD computations where xi_1 comes first
            const auto & Pi_muY01    = scratch1D_5; // used for E_QUAD computations where xi_2 comes first
            const auto & Pi_muZ01    = scratch1D_6; // used for E_QUAD computations where xi_2 comes first
            const auto & Li_dt_muX01 = scratch1D_7; // used for E_QUAD computations
            const auto & Li_dt_muY01 = scratch1D_8; // used for E_QUAD computations
            const auto & Li_dt_muZ01 = scratch1D_9; // used for E_QUAD computations
            
            const auto & muX_0 = mu[0][0]; const auto & muX_0_grad = muGrad[0][0];
            const auto & muX_1 = mu[1][0]; const auto & muX_1_grad = muGrad[1][0];
            const auto & muY_0 = mu[0][1]; const auto & muY_0_grad = muGrad[0][1];
            const auto & muY_1 = mu[1][1]; const auto & muY_1_grad = muGrad[1][1];
            const auto & muZ_0 = mu[0][2]; const auto & muZ_0_grad = muGrad[0][2];
            const auto & muZ_1 = mu[1][2]; const auto & muZ_1_grad = muGrad[1][2];
            
            Polynomials::shiftedScaledIntegratedLegendreValues(Li_muX01, polyOrder_, muX_1, muX_0 + muX_1);
            Polynomials::shiftedScaledIntegratedLegendreValues(Li_muY01, polyOrder_, muY_1, muY_0 + muY_1);
            Polynomials::shiftedScaledIntegratedLegendreValues(Li_muZ01, polyOrder_, muZ_1, muZ_0 + muZ_1);
            
            Polynomials::shiftedScaledLegendreValues(Pi_muX01, polyOrder_, muX_1, muX_0 + muX_1);
            Polynomials::shiftedScaledLegendreValues(Pi_muY01, polyOrder_, muY_1, muY_0 + muY_1);
            Polynomials::shiftedScaledLegendreValues(Pi_muZ01, polyOrder_, muZ_1, muZ_0 + muZ_1);
            
            Polynomials::shiftedScaledIntegratedLegendreValues_dt(Li_dt_muX01, Pi_muX01, polyOrder_, muX_1, muX_0 + muX_1);
            Polynomials::shiftedScaledIntegratedLegendreValues_dt(Li_dt_muY01, Pi_muY01, polyOrder_, muY_1, muY_0 + muY_1);
            Polynomials::shiftedScaledIntegratedLegendreValues_dt(Li_dt_muZ01, Pi_muZ01, polyOrder_, muZ_1, muZ_0 + muZ_1);
            
            // FAMILIES I & II -- divergence-free families
            // following the ESEAS ordering: k increments first
            for (int f=0; f<2; f++)
            {
              const auto &s0 = (f==0) ? muX_0 : muY_0;  const auto & s0_grad = (f==0) ? muX_0_grad : muY_0_grad;
              const auto &s1 = (f==0) ? muX_1 : muY_1;  const auto & s1_grad = (f==0) ? muX_1_grad : muY_1_grad;
                                                        const auto & t0_grad = (f==0) ? muY_0_grad : muX_0_grad;
                                                        const auto & t1_grad = (f==0) ? muY_1_grad : muX_1_grad;
              const auto & Pi_s01    = (f==0) ? Pi_muX01    : Pi_muY01;
              const auto & Pi_t01    = (f==0) ? Pi_muY01    : Pi_muX01;
              const auto & Li_t01    = (f==0) ? Li_muY01    : Li_muX01;
              const auto & Li_dt_t01 = (f==0) ? Li_dt_muY01 : Li_dt_muX01;
              
              for (int k=2; k<=polyOrder_; k++)
              {
                const auto & phi_k = Li_muZ01(k);
                Kokkos::Array<OutputScalar,3> phi_k_grad;
                computeGradHomLi(phi_k_grad, k, Pi_muZ01, Li_dt_muZ01, muZ_0_grad, muZ_1_grad);
                
                Kokkos::Array<OutputScalar,3> muZ0_grad_phi_k_plus_phi_k_grad_muZ0;
                for (ordinal_type d=0; d<3; d++)
                {
                  muZ0_grad_phi_k_plus_phi_k_grad_muZ0[d] = muZ_0 * phi_k_grad[d] + phi_k * muZ_0_grad[d];
                }
                
                for (int j=2; j<=polyOrder_; j++)
                {
                  for (int i=0; i<polyOrder_; i++)
                  {
                    Kokkos::Array<OutputScalar,3> EQUAD_ij;
                    Kokkos::Array<OutputScalar,3> curl_EQUAD_ij;
                    
                    E_QUAD(EQUAD_ij, i, j, Pi_s01, s0, s1, s0_grad, s1_grad, Li_t01);
                    
                    E_QUAD_CURL(curl_EQUAD_ij, i, j, Pi_s01, s0, s1, s0_grad, s1_grad,
                                Pi_t01, Li_t01, Li_dt_t01, t0_grad, t1_grad);
                    
                    // first term: muZ_0 phi^E_k curl EQUAD
                    // we can reuse the memory for curl_EQUAD_ij_12; we won't need the values there again
                    Kokkos::Array<OutputScalar,3> & firstTerm = curl_EQUAD_ij;
                    for (ordinal_type d=0; d<3; d++)
                    {
                      firstTerm[d] *= muZ_0 * phi_k;
                    }
                    
                    Kokkos::Array<OutputScalar,3> secondTerm; //(muZ0 grad phi + phi grad muZ0) x EQUAD
                    
                    cross(secondTerm, muZ0_grad_phi_k_plus_phi_k_grad_muZ0, EQUAD_ij);
                    
                    for (ordinal_type d=0; d<3; d++)
                    {
                      output_(fieldOrdinalOffset,pointOrdinal,d) = firstTerm[d] + secondTerm[d];
                    }
                    
                    fieldOrdinalOffset++;
                  }
                }
              }
            } // family I, II loop
            
            // FAMILY III -- a divergence-free family
            for (int j=2; j<=polyOrder_; j++)
            {
              for (int i=2; i<=polyOrder_; i++)
              {
                fieldOrdinalOffset++;
              }
            }
            
            // FAMILY IV (non-trivial divergences)
            for (int k=2; k<=polyOrder_; k++)
            {
              for (int j=0; j<polyOrder_; j++)
              {
                for (int i=0; i<polyOrder_; i++)
                {
                  fieldOrdinalOffset++;
                }
              }
            }
            
            // FAMILY V (non-trivial divergences)
            for (int j=2; j<=polyOrder_; j++)
            {
              for (int i=2; i<=polyOrder_; i++)
              {
                fieldOrdinalOffset++;
              }
            }
            
            // FAMILY VI (non-trivial divergences)
            for (int i=2; i<=polyOrder_; i++)
            {
              fieldOrdinalOffset++;
            }
            
            // FAMILY VII (non-trivial divergences)
            for (int j=2; j<=polyOrder_; j++)
            {
              fieldOrdinalOffset++;
            }
          }
          
          
          // end rewritten portion
          
          // TODO: interior functions (below is from H^1 implementation)
          // interior functions
          // these are the product of the same quadrilateral function that we blended for the quadrilateral face and
          // [L_k](mu_{0}^zeta, mu_1^zeta)
            
          // rename scratch
//                 Li = scratch1D_1;
//                 Lj = scratch1D_2;
//          auto & Lk = scratch1D_3;
//          Polynomials::shiftedScaledIntegratedLegendreValues(Li, polyOrder_, mu[1][0], mu[0][0] + mu[1][0]);
//          Polynomials::shiftedScaledIntegratedLegendreValues(Lj, polyOrder_, mu[1][1], mu[0][1] + mu[1][1]);
//          Polynomials::shiftedScaledIntegratedLegendreValues(Lk, polyOrder_, mu[1][2], mu[0][2] + mu[1][2]);
//          // following the ESEAS ordering: k increments first
//          for (int k=2; k<=polyOrder_; k++)
//          {
//            for (int j=2; j<=polyOrder_; j++)
//            {
//              for (int i=2; i<=polyOrder_; i++)
//              {
//                output_(fieldOrdinalOffset,pointOrdinal) = Lk(k) * Li(i) * Lj(j);
//                fieldOrdinalOffset++;
//              }
//            }
//          }
        } // end OPERATOR_VALUE
          break;
        case OPERATOR_DIV:
        {
          ordinal_type fieldOrdinalOffset = 0;
          // quadrilateral face
          {
            // rename scratch1, scratch2
            auto & Pi = scratch1D_1;
            auto & Pj = scratch1D_2;
            
            auto & s0      =     mu[0][0], s1      =     mu[1][0];
            auto & s0_grad = muGrad[0][0], s1_grad = muGrad[1][0];
            auto & t0      =     mu[0][1], t1      =     mu[1][1];
            auto & t0_grad = muGrad[0][1], t1_grad = muGrad[1][1];
            
            Polynomials::shiftedScaledLegendreValues(Pi, polyOrder_-1, s1, s0 + s1);
            Polynomials::shiftedScaledLegendreValues(Pj, polyOrder_-1, t1, t0 + t1);
            
            const auto & muZ0      =     mu[0][2];
            const auto & muZ0_grad = muGrad[0][2];
            OutputScalar three_mu0_squared = 3.0 * muZ0 * muZ0;
            
            // following the ESEAS ordering: j increments first
            for (int j=0; j<polyOrder_; j++)
            {
              for (int i=0; i<polyOrder_; i++)
              {
                Kokkos::Array<OutputScalar,3> VQUAD; // output from V_QUAD
                V_QUAD(VQUAD, i, j,
                       Pi, s0, s1, s0_grad, s1_grad,
                       Pj, t0, t1, t0_grad, t1_grad);
                
                OutputScalar grad_muZ0_dot_VQUAD;
                dot(grad_muZ0_dot_VQUAD, muZ0_grad, VQUAD);
                
                output_(fieldOrdinalOffset,pointOrdinal) = three_mu0_squared * grad_muZ0_dot_VQUAD;
                fieldOrdinalOffset++;
              }
            }
          } // end quad face block
          
          // triangle faces
          {
            const auto & P        = scratch1D_1; // for V_TRI( nu_012)
            const auto & P_2ip1   = scratch1D_2;
            const auto & Pmu      = scratch1D_3; // for V_TRI( mu * nu_012)
            const auto & Pmu_2ip1 = scratch1D_4;
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
              const PointScalar jacobiScaling = s0 + s1 + s2; // we can actually assume that this is 1; see comment at bottom of p. 425 of Fuentes et al.
              
              const PointScalar legendreScaling = s0 + s1;
              Polynomials::shiftedScaledLegendreValues(P, polyOrder_-1, s1, legendreScaling);
              
              const auto lambda0_index = tri_face_vertex_0[faceOrdinal];
              const auto lambda1_index = tri_face_vertex_1[faceOrdinal];
              const auto lambda2_index = tri_face_vertex_2[faceOrdinal];
              
              const auto & mu_c_b = mu[c][b-1];
              const auto & mu_c_b_grad = muGrad[c][b-1];
              
              const auto & mu_s0 = lambda[lambda0_index];
              const auto & mu_s1 = lambda[lambda1_index];
              const auto & mu_s2 = lambda[lambda2_index]; // == s2
              
              const PointScalar muJacobiScaling = mu_s0 + mu_s1 + mu_s2;
              
              const PointScalar muLegendreScaling = mu_s0 + mu_s1;
              Polynomials::shiftedScaledLegendreValues(Pmu, polyOrder_-1, mu_s1, muLegendreScaling);
              
              Kokkos::Array<PointScalar, 3> vectorWeight;
              computeFaceVectorWeight(vectorWeight, a, nu, nuGrad);
              
              Kokkos::Array<PointScalar,3> & mu_s0_grad = lambdaGrad[lambda0_index];
              Kokkos::Array<PointScalar,3> & mu_s1_grad = lambdaGrad[lambda1_index];
              Kokkos::Array<PointScalar,3> & mu_s2_grad = lambdaGrad[lambda2_index]; // == s2_grad
              
              Kokkos::Array<PointScalar, 3> muVectorWeight;
              computeFaceVectorWeight(muVectorWeight, mu_s0, mu_s0_grad, mu_s1, mu_s1_grad, mu_s2, mu_s2_grad);
              
              OutputScalar muDivWeight;
              computeFaceDivWeight(muDivWeight, mu_s0_grad, mu_s1_grad, mu_s2_grad);
              
              for (int totalPolyOrder=0; totalPolyOrder<polyOrder_; totalPolyOrder++)
              {
                for (int i=0; i<=totalPolyOrder; i++)
                {
                  const int j = totalPolyOrder - i;
                  
                  const double alpha = i*2.0 + 1;
                  Polynomials::shiftedScaledJacobiValues(  P_2ip1, alpha, polyOrder_-1,    s2,   jacobiScaling);
                  Polynomials::shiftedScaledJacobiValues(Pmu_2ip1, alpha, polyOrder_-1, mu_s2, muJacobiScaling);
                  
                  Kokkos::Array<OutputScalar,3> VTRI, VTRI_mu; // output from V_TRI
                  
                  V_TRI(VTRI,    i, j, P,   P_2ip1,     vectorWeight);
                  V_TRI(VTRI_mu, i, j, Pmu, Pmu_2ip1, muVectorWeight);
                  
                  OutputScalar VTRI_mu_DIV;
                  V_TRI_DIV(VTRI_mu_DIV, i, j, Pmu, Pmu_2ip1, muDivWeight);
                  
                  output_(fieldOrdinalOffset,pointOrdinal) = 0.5 * (dot(mu_c_b_grad, VTRI) + VTRI_mu_DIV / mu_c_b - dot(mu_c_b_grad, VTRI_mu)/(mu_c_b * mu_c_b));
                  
                  fieldOrdinalOffset++;
                }
              }
            }
          } // end triangle face block
          // TODO: interior functions
          
          {
            // FAMILY I -- divergence free
            // following the ESEAS ordering: k increments first
            for (int k=2; k<=polyOrder_; k++)
            {
              for (int j=2; j<=polyOrder_; j++)
              {
                for (int i=0; i<polyOrder_; i++)
                {
                  output_(fieldOrdinalOffset,pointOrdinal) = 0.0;
                  fieldOrdinalOffset++;
                }
              }
            }
            
            // FAMILY II -- divergence free
            // following the ESEAS ordering: k increments first
            for (int k=2; k<=polyOrder_; k++)
            {
              for (int j=2; j<=polyOrder_; j++)
              {
                for (int i=0; i<polyOrder_; i++)
                {
                  output_(fieldOrdinalOffset,pointOrdinal) = 0.0;
                  fieldOrdinalOffset++;
                }
              }
            }
            
            // FAMILY III -- divergence free
            for (int j=2; j<=polyOrder_; j++)
            {
              for (int i=2; i<=polyOrder_; i++)
              {
                output_(fieldOrdinalOffset,pointOrdinal) = 0.0;
                fieldOrdinalOffset++;
              }
            }
            
          } // end interior function block
          
        } // end OPERATOR_DIV block
          break;
        case OPERATOR_GRAD:
        case OPERATOR_D1:
        {
          // TODO: rewrite this -- the below comes from the H^1 implementation
          // vertex functions
          
          /*for (int vertexOrdinal=0; vertexOrdinal<numVertices; vertexOrdinal++)
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
            
            Polynomials::shiftedScaledLegendreValues             (P_i_minus_1, polyOrder_-1, nu[1][a-1], nu[0][a-1] + nu[1][a-1]);
            Polynomials::shiftedScaledIntegratedLegendreValues_dt(L_i_dt,      polyOrder_,   nu[1][a-1], nu[0][a-1] + nu[1][a-1]);
            Polynomials::shiftedScaledIntegratedLegendreValues   (L_i,         polyOrder_,   nu[1][a-1], nu[0][a-1] + nu[1][a-1]);
            
            // edge 0,3 --> c=0
            // edge 1,2 --> c=1
            int c = ((edgeOrdinal == 0) || (edgeOrdinal == 3)) ? 0 : 1;
            for (int i=2; i<=polyOrder_; i++)
            {
              // grad (mu[c][b-1] * Li(i)) = grad (mu[c][b-1]) * Li(i) + mu[c][b-1] * grad Li(i)
              const auto & R_i_minus_1 = L_i_dt(i);
              
              for (int d=0; d<3; d++)
              {
                // grad [L_i](nu_0,nu_1) = [P_{i-1}](nu_0,nu_1) * grad nu_1 + [R_{i-1}](nu_0,nu_1) * grad (nu_0 + nu_1)
                
                OutputScalar grad_Li_d = P_i_minus_1(i-1) * nuGrad[1][a-1][d] + R_i_minus_1 * (nuGrad[0][a-1][d] + nuGrad[1][a-1][d]);
                output_(fieldOrdinalOffset,pointOrdinal,d) = muGrad[c][b-1][d] * L_i(i) + mu[c][b-1] * grad_Li_d;
              }
              fieldOrdinalOffset++;
            }
          }
          
          // triangle edges next
          P_i_minus_1 = scratch1D_1;
          L_i_dt      = scratch1D_2; // R_{i-1} = d/dt L_i
          L_i         = scratch1D_3;
          for (int edgeOrdinal=0; edgeOrdinal<numMixedEdges; edgeOrdinal++)
          {
            const auto & lambda_a     = lambda    [edgeOrdinal];
            const auto & lambdaGrad_a = lambdaGrad[edgeOrdinal];
            Polynomials::shiftedScaledLegendreValues             (P_i_minus_1, polyOrder_-1, lambda[4], lambda_a + lambda[4]);
            Polynomials::shiftedScaledIntegratedLegendreValues_dt(L_i_dt,      polyOrder_,   lambda[4], lambda_a + lambda[4]);
            Polynomials::shiftedScaledIntegratedLegendreValues   (L_i,         polyOrder_,   lambda[4], lambda_a + lambda[4]);
            
            for (int i=2; i<=polyOrder_; i++)
            {
              const auto & R_i_minus_1 = L_i_dt(i);
              for (int d=0; d<3; d++)
              {
                // grad [L_i](s0,s1) = [P_{i-1}](s0,s1) * grad s1 + [R_{i-1}](s0,s1) * grad (s0 + s1)
                // here, s1 = lambda[4], s0 = lambda_a
                OutputScalar grad_Li_d = P_i_minus_1(i-1) * lambdaGrad[4][d] + R_i_minus_1 * (lambdaGrad_a[d] + lambdaGrad[4][d]);
                output_(fieldOrdinalOffset,pointOrdinal,d) = grad_Li_d;
              }
              fieldOrdinalOffset++;
            }
          }
          
          // quadrilateral faces
          // rename scratch
          P_i_minus_1 = scratch1D_1;
          L_i_dt      = scratch1D_2; // R_{i-1} = d/dt L_i
          L_i         = scratch1D_3;
          auto & P_j_minus_1 = scratch1D_4;
          auto & L_j_dt      = scratch1D_5; // R_{j-1} = d/dt L_j
          auto & L_j         = scratch1D_6;
          Polynomials::shiftedScaledIntegratedLegendreValues(L_i, polyOrder_, mu[1][0], mu[0][0] + mu[1][0]);
          Polynomials::shiftedScaledIntegratedLegendreValues(L_j, polyOrder_, mu[1][1], mu[0][1] + mu[1][1]);
          
          Polynomials::shiftedScaledLegendreValues             (P_i_minus_1, polyOrder_-1, mu[1][0], mu[0][0] + mu[1][0]);
          Polynomials::shiftedScaledIntegratedLegendreValues_dt(L_i_dt,      polyOrder_,   mu[1][0], mu[0][0] + mu[1][0]);
          Polynomials::shiftedScaledIntegratedLegendreValues   (L_i,         polyOrder_,   mu[1][0], mu[0][0] + mu[1][0]);
          
          Polynomials::shiftedScaledLegendreValues             (P_j_minus_1, polyOrder_-1, mu[1][1], mu[0][1] + mu[1][1]);
          Polynomials::shiftedScaledIntegratedLegendreValues_dt(L_j_dt,      polyOrder_,   mu[1][1], mu[0][1] + mu[1][1]);
          Polynomials::shiftedScaledIntegratedLegendreValues   (L_j,         polyOrder_,   mu[1][1], mu[0][1] + mu[1][1]);
          
          // following the ESEAS ordering: j increments first
          for (int j=2; j<=polyOrder_; j++)
          {
            const auto & R_j_minus_1 = L_j_dt(j);
            
            for (int i=2; i<=polyOrder_; i++)
            {
              const auto & R_i_minus_1 = L_i_dt(i);
              
              OutputScalar phi_quad = L_i(i) * L_j(j);
              
              for (int d=0; d<3; d++)
              {
                // grad [L_j](s0,s1) = [P_{j-1}](s0,s1) * grad s1 + [R_{j-1}](s0,s1) * grad (s0 + s1)
                // here, s1 = mu[1][1], s0 = mu[0][1]
                OutputScalar grad_Lj_d = P_j_minus_1(j-1) * muGrad[1][1][d] + R_j_minus_1 * (muGrad[0][1][d] + muGrad[1][1][d]);
                // for L_i, s1 = mu[1][0], s0 = mu[0][0]
                OutputScalar grad_Li_d = P_i_minus_1(i-1) * muGrad[1][0][d] + R_i_minus_1 * (muGrad[0][0][d] + muGrad[1][0][d]);
                
                OutputScalar grad_phi_quad_d = L_i(i) * grad_Lj_d + L_j(j) * grad_Li_d;
                
                output_(fieldOrdinalOffset,pointOrdinal,d) = mu[0][2] * grad_phi_quad_d + phi_quad * muGrad[0][2][d];
              }
              
              fieldOrdinalOffset++;
            }
          }
          
          // triangle faces
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
            
            const auto & s0Grad = nuGrad[0][a-1];
            const auto & s1Grad = nuGrad[1][a-1];
            const auto & s2Grad = nuGrad[2][a-1];
            
            const PointScalar jacobiScaling = s0 + s1 + s2;
            
            // compute integrated Jacobi values for each desired value of alpha
            // relabel storage:
            // 1D containers:
            auto & P_i_minus_1 = scratch1D_1;
            auto & L_i_dt      = scratch1D_2; // R_{i-1} = d/dt L_i
            auto & L_i         = scratch1D_3;
            // 2D containers:
            auto & L_2i_j_dt      = scratch2D_1;
            auto & L_2i_j         = scratch2D_2;
            auto & P_2i_j_minus_1 = scratch2D_3;
            for (int n=2; n<=polyOrder_; n++)
            {
              const double alpha = n*2;
              const int alphaOrdinal = n-2;
              using Kokkos::subview;
              using Kokkos::ALL;
              auto L_2i_j_dt_alpha      = subview(L_2i_j_dt,      alphaOrdinal, ALL);
              auto L_2i_j_alpha         = subview(L_2i_j,         alphaOrdinal, ALL);
              auto P_2i_j_minus_1_alpha = subview(P_2i_j_minus_1, alphaOrdinal, ALL);
              Polynomials::shiftedScaledIntegratedJacobiValues_dt(L_2i_j_dt_alpha, alpha, polyOrder_-2, s2, jacobiScaling);
              Polynomials::shiftedScaledIntegratedJacobiValues   (   L_2i_j_alpha, alpha, polyOrder_-2, s2, jacobiScaling);
              Polynomials::shiftedScaledJacobiValues        (P_2i_j_minus_1_alpha, alpha, polyOrder_-1, s2, jacobiScaling);
            }
            Polynomials::shiftedScaledLegendreValues             (P_i_minus_1, polyOrder_-1, s1, s0 + s1);
            Polynomials::shiftedScaledIntegratedLegendreValues_dt(L_i_dt,      polyOrder_,   s1, s0 + s1);
            Polynomials::shiftedScaledIntegratedLegendreValues   (L_i,         polyOrder_,   s1, s0 + s1);
            
            for (int totalPolyOrder=3; totalPolyOrder<=polyOrder_; totalPolyOrder++)
            {
              for (int i=2; i<totalPolyOrder; i++)
              {
                const int alphaOrdinal = i-2;
                const int            j = totalPolyOrder - i;
                
                const auto & R_i_minus_1 = L_i_dt(i);
                OutputScalar     phi_tri = L_2i_j(alphaOrdinal,j) * L_i(i);
                
                for (int d=0; d<3; d++)
                {
                  // grad [L_i](s0,s1) = [P_{i-1}](s0,s1) * grad s1 + [R_{i-1}](s0,s1) * grad (s0 + s1)
                  OutputScalar grad_Li_d      = P_i_minus_1(i-1) * s1Grad[d] + R_i_minus_1 * (s0Grad[d] + s1Grad[d]);
                  OutputScalar grad_L2i_j_d   = P_2i_j_minus_1(alphaOrdinal,j-1) * s2Grad[d] + L_2i_j_dt(alphaOrdinal,j) * (s0Grad[d] + s1Grad[d] + s2Grad[d]);
                  OutputScalar grad_phi_tri_d = L_i(i) * grad_L2i_j_d + L_2i_j(alphaOrdinal,j) * grad_Li_d;
                  
                  output_(fieldOrdinalOffset,pointOrdinal,d) = mu[c][b-1] * grad_phi_tri_d + phi_tri * muGrad[c][b-1][d];
                }
                fieldOrdinalOffset++;
              }
            }
          }
          
          // interior functions
          P_i_minus_1 = scratch1D_1;
          L_i_dt      = scratch1D_2; // R_{i-1} = d/dt L_i
          L_i         = scratch1D_3;
          P_j_minus_1 = scratch1D_4;
          L_j_dt      = scratch1D_5; // R_{j-1} = d/dt L_j
          L_j         = scratch1D_6;
          auto & P_k_minus_1 = scratch1D_7;
          auto & L_k_dt      = scratch1D_8; // R_{k-1} = d/dt L_k
          auto & L_k         = scratch1D_9;
          
          Polynomials::shiftedScaledLegendreValues             (P_i_minus_1, polyOrder_-1, mu[1][0], mu[0][0] + mu[1][0]);
          Polynomials::shiftedScaledIntegratedLegendreValues_dt(L_i_dt,      polyOrder_,   mu[1][0], mu[0][0] + mu[1][0]);
          Polynomials::shiftedScaledIntegratedLegendreValues   (L_i,         polyOrder_,   mu[1][0], mu[0][0] + mu[1][0]);
          
          Polynomials::shiftedScaledLegendreValues             (P_j_minus_1, polyOrder_-1, mu[1][1], mu[0][1] + mu[1][1]);
          Polynomials::shiftedScaledIntegratedLegendreValues_dt(L_j_dt,      polyOrder_,   mu[1][1], mu[0][1] + mu[1][1]);
          Polynomials::shiftedScaledIntegratedLegendreValues   (L_j,         polyOrder_,   mu[1][1], mu[0][1] + mu[1][1]);
          
          Polynomials::shiftedScaledLegendreValues             (P_k_minus_1, polyOrder_-1, mu[1][2], mu[0][2] + mu[1][2]);
          Polynomials::shiftedScaledIntegratedLegendreValues_dt(L_k_dt,      polyOrder_,   mu[1][2], mu[0][2] + mu[1][2]);
          Polynomials::shiftedScaledIntegratedLegendreValues   (L_k,         polyOrder_,   mu[1][2], mu[0][2] + mu[1][2]);
          
          // following the ESEAS ordering: k increments first
          for (int k=2; k<=polyOrder_; k++)
          {
            const auto & R_k_minus_1 = L_k_dt(k);
            
            for (int j=2; j<=polyOrder_; j++)
            {
              const auto & R_j_minus_1 = L_j_dt(j);
              
              for (int i=2; i<=polyOrder_; i++)
              {
                const auto & R_i_minus_1 = L_i_dt(i);
                
                OutputScalar phi_quad = L_i(i) * L_j(j);
                
                for (int d=0; d<3; d++)
                {
                  // grad [L_i](s0,s1) = [P_{i-1}](s0,s1) * grad s1 + [R_{i-1}](s0,s1) * grad (s0 + s1)
                  // for L_i, s1 = mu[1][0], s0 = mu[0][0]
                  OutputScalar grad_Li_d = P_i_minus_1(i-1) * muGrad[1][0][d] + R_i_minus_1 * (muGrad[0][0][d] + muGrad[1][0][d]);
                  // for L_j, s1 = mu[1][1], s0 = mu[0][1]
                  OutputScalar grad_Lj_d = P_j_minus_1(j-1) * muGrad[1][1][d] + R_j_minus_1 * (muGrad[0][1][d] + muGrad[1][1][d]);
                  // for L_k, s1 = mu[1][2], s0 = mu[0][2]
                  OutputScalar grad_Lk_d = P_k_minus_1(k-1) * muGrad[1][2][d] + R_k_minus_1 * (muGrad[0][2][d] + muGrad[1][2][d]);
                  
                  OutputScalar grad_phi_quad_d = L_i(i) * grad_Lj_d + L_j(j) * grad_Li_d;
                  
                  output_(fieldOrdinalOffset,pointOrdinal,d) = L_k(k) * grad_phi_quad_d + phi_quad * grad_Lk_d;
                }
                
                fieldOrdinalOffset++;
              }
            }
          }

          for (int basisOrdinal=0; basisOrdinal<numFields_; basisOrdinal++)
          {
            // transform derivatives to account for the ref space transformation: Intrepid2 uses (-z,z)^2; ESEAS uses (0,z)^2
            const auto dx_eseas = output_(basisOrdinal,pointOrdinal,0);
            const auto dy_eseas = output_(basisOrdinal,pointOrdinal,1);
            const auto dz_eseas = output_(basisOrdinal,pointOrdinal,2);
            
            auto &dx_int2 = output_(basisOrdinal,pointOrdinal,0);
            auto &dy_int2 = output_(basisOrdinal,pointOrdinal,1);
            auto &dz_int2 = output_(basisOrdinal,pointOrdinal,2);
            
            transformFromESEASPyramidGradient(dx_int2, dy_int2, dz_int2, dx_eseas, dy_eseas, dz_eseas);
          }*/
          
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
                                   ">>> ERROR: (Intrepid2::Hierarchical_HDIV_PYR_Functor) Computing of second and higher-order derivatives is not currently supported");
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
  
  /** \class  Intrepid2::HierarchicalBasis_HDIV_PYR
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
           typename PointScalar  = double,
           bool useCGBasis = true>  // if useCGBasis is true, basis functions are associated with either faces or the interior.  If false, basis functions are all associated with interior
  class HierarchicalBasis_HDIV_PYR
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
     
     If useCGBasis is false, then all basis functions are identified with the interior of the element.
     
     If useCGBasis is true, then basis functions are are associated with either faces or the interior.
     
     */
    HierarchicalBasis_HDIV_PYR(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    polyOrder_(polyOrder),
    pointType_(pointType)
    {
      INTREPID2_TEST_FOR_EXCEPTION(pointType!=POINTTYPE_DEFAULT,std::invalid_argument,"PointType not supported");
      const auto & p           = polyOrder;
      this->basisCardinality_  = p * p + 2 * p * (p+1) + 3 * p * p * (p-1);
      this->basisDegree_       = p;
      this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Pyramid<> >() );
      this->basisType_         = BASIS_FEM_HIERARCHICAL;
      this->basisCoordinates_  = COORDINATES_CARTESIAN;
      this->functionSpace_     = FUNCTION_SPACE_HDIV;
      
      const int degreeLength = 1;
      this->fieldOrdinalPolynomialDegree_ = OrdinalTypeArray2DHost("Integrated Legendre H(grad) pyramid polynomial degree lookup", this->basisCardinality_, degreeLength);
      this->fieldOrdinalH1PolynomialDegree_ = OrdinalTypeArray2DHost("Integrated Legendre H(grad) pyramid polynomial H^1 degree lookup", this->basisCardinality_, degreeLength);
      
      int fieldOrdinalOffset = 0;
      
      // **** face functions **** //
      // one quad face
      const int numFunctionsPerQuadFace = p*p;
      
      // following the ESEAS ordering: j increments first
      for (int j=0; j<p; j++)
      {
        for (int i=0; i<p; i++)
        {
          this->fieldOrdinalPolynomialDegree_  (fieldOrdinalOffset,0) = std::max(i,j);
          this->fieldOrdinalH1PolynomialDegree_(fieldOrdinalOffset,0) = std::max(i,j)+1;
          fieldOrdinalOffset++;
        }
      }
      
      const int numFunctionsPerTriFace = 2 * p * (p+1) / 4;
      const int numTriFaces = 4;
      for (int faceOrdinal=0; faceOrdinal<numTriFaces; faceOrdinal++)
      {
        for (int totalPolyOrder=0; totalPolyOrder<polyOrder_; totalPolyOrder++)
        {
          const int totalFaceDofs         = (totalPolyOrder+1) * (totalPolyOrder+2) / 2; // number of dofs for which i+j <= totalPolyOrder
          const int totalFaceDofsPrevious =  totalPolyOrder    * (totalPolyOrder+1) / 2; // number of dofs for which i+j <= totalPolyOrder-1
          const int faceDofsForPolyOrder  = totalFaceDofs - totalFaceDofsPrevious; // number of dofs for which i+j == totalPolyOrder
          for (int i=0; i<faceDofsForPolyOrder; i++)
          {
            this->fieldOrdinalPolynomialDegree_  (fieldOrdinalOffset,0) = totalPolyOrder;
            this->fieldOrdinalH1PolynomialDegree_(fieldOrdinalOffset,0) = totalPolyOrder+1;
            fieldOrdinalOffset++;
          }
        }
      }

      // **** interior functions **** //
      const int numFunctionsPerVolume = 3 * p * p * (p-1);
      
      // FAMILY I
      // following the ESEAS ordering: k increments first
      for (int k=2; k<=polyOrder_; k++)
      {
        for (int j=2; j<=polyOrder_; j++)
        {
          for (int i=0; i<polyOrder_; i++)
          {
            const int max_jk  = std::max(j,k);
            const int max_ijk = std::max(max_jk,i);
            const int max_ip1jk = std::max(max_jk,i+1);
            this->fieldOrdinalPolynomialDegree_  (fieldOrdinalOffset,0) = max_ijk;
            this->fieldOrdinalH1PolynomialDegree_(fieldOrdinalOffset,0) = max_ip1jk;
            fieldOrdinalOffset++;
          }
        }
      }
      
      // FAMILY II
      for (int k=2; k<=polyOrder_; k++)
      {
        for (int j=2; j<=polyOrder_; j++)
        {
          for (int i=0; i<polyOrder_; i++)
          {
            const int max_jk  = std::max(j,k);
            const int max_ijk = std::max(max_jk,i);
            const int max_ip1jk = std::max(max_jk,i+1);
            this->fieldOrdinalPolynomialDegree_  (fieldOrdinalOffset,0) = max_ijk;
            this->fieldOrdinalH1PolynomialDegree_(fieldOrdinalOffset,0) = max_ip1jk;
            fieldOrdinalOffset++;
          }
        }
      }
      
      // FAMILY III
      for (int j=2; j<=polyOrder_; j++)
      {
        for (int i=2; i<=polyOrder_; i++)
        {
          const int max_ij  = std::max(i,j);
          this->fieldOrdinalPolynomialDegree_  (fieldOrdinalOffset,0) = max_ij;
          this->fieldOrdinalH1PolynomialDegree_(fieldOrdinalOffset,0) = max_ij;
          fieldOrdinalOffset++;
        }
      }
      
      // FAMILY IV
      for (int k=2; k<=polyOrder_; k++)
      {
        for (int j=0; j<polyOrder_; j++)
        {
          for (int i=0; i<polyOrder_; i++)
          {
            const int max_jk  = std::max(j,k);
            const int max_ijk = std::max(max_jk,i);
            const int max_ip1jp1k = std::max(std::max(j+1,k),i+1);
            this->fieldOrdinalPolynomialDegree_  (fieldOrdinalOffset,0) = max_ijk;
            this->fieldOrdinalH1PolynomialDegree_(fieldOrdinalOffset,0) = max_ip1jp1k;
            fieldOrdinalOffset++;
          }
        }
      }
      
      // FAMILY V
      for (int j=2; j<=polyOrder_; j++)
      {
        for (int i=2; i<=polyOrder_; i++)
        {
          const int max_ij  = std::max(i,j);
          this->fieldOrdinalPolynomialDegree_  (fieldOrdinalOffset,0) = max_ij;
          this->fieldOrdinalH1PolynomialDegree_(fieldOrdinalOffset,0) = max_ij;
          fieldOrdinalOffset++;
        }
      }
      
      // FAMILY VI
      for (int i=2; i<=polyOrder_; i++)
      {
        this->fieldOrdinalPolynomialDegree_  (fieldOrdinalOffset,0) = i;
        this->fieldOrdinalH1PolynomialDegree_(fieldOrdinalOffset,0) = i;
        fieldOrdinalOffset++;
      }
      
      // FAMILY VII
      for (int j=2; j<=polyOrder_; j++)
      {
        this->fieldOrdinalPolynomialDegree_  (fieldOrdinalOffset,0) = j;
        this->fieldOrdinalH1PolynomialDegree_(fieldOrdinalOffset,0) = j;
        fieldOrdinalOffset++;
      }
      
      if (fieldOrdinalOffset != this->basisCardinality_)
      {
        std::cout << "Internal error: basis enumeration is incorrect; fieldOrdinalOffset = " << fieldOrdinalOffset;
        std::cout << "; basisCardinality_ = " << this->basisCardinality_ << std::endl;
        
        INTREPID2_TEST_FOR_EXCEPTION(fieldOrdinalOffset != this->basisCardinality_, std::invalid_argument, "Internal error: basis enumeration is incorrect");
      }
      
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
        const int faceDim = 2, volumeDim = 3;

        if (useCGBasis) {
          {
            int tagNumber = 0;
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
            for (int functionOrdinal=0; functionOrdinal<numFunctionsPerVolume; functionOrdinal++)
            {
              tagView(tagNumber*tagSize+0) = volumeDim;               // volume dimension
              tagView(tagNumber*tagSize+1) = 0;                       // volume id
              tagView(tagNumber*tagSize+2) = functionOrdinal;         // local dof id
              tagView(tagNumber*tagSize+3) = numFunctionsPerVolume;   // total number of dofs in this volume
              tagNumber++;
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
      return "Intrepid2_HierarchicalBasis_HDIV_PYR";
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
      
      using FunctorType = Hierarchical_HDIV_PYR_Functor<DeviceType, OutputScalar, PointScalar, OutputViewType, PointViewType>;
      
      FunctorType functor(operatorType, outputValues, inputPoints, polyOrder_);
      
      const int outputVectorSize = getVectorSizeForHierarchicalParallelism<OutputScalar>();
      const int pointVectorSize  = getVectorSizeForHierarchicalParallelism<PointScalar>();
      const int vectorSize = std::max(outputVectorSize,pointVectorSize);
      const int teamSize = 1; // because of the way the basis functions are computed, we don't have a second level of parallelism...

      auto policy = Kokkos::TeamPolicy<ExecutionSpace>(numPoints,teamSize,vectorSize);
      Kokkos::parallel_for("Hierarchical_HDIV_PYR_Functor", policy, functor);
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
      const auto & p = this->basisDegree_;
      if (subCellDim == 2)
      {
        if (subCellOrd == 4) // quad basis
        {
          using HVOL_LINE = LegendreBasis_HVOL_LINE<DeviceType,OutputScalar,PointScalar>;
          using HVOL_QUAD = Basis_Derived_HVOL_QUAD<HVOL_LINE>;
          return Teuchos::rcp(new HVOL_QUAD(p-1));
        }
        else // tri basis
        {
          using HVOL_Tri = LegendreBasis_HVOL_TRI<DeviceType,OutputScalar,PointScalar>;
          return Teuchos::rcp(new HVOL_Tri(p-1));
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
      using HostBasisType  = HierarchicalBasis_HDIV_PYR<HostDeviceType, OutputScalar, PointScalar, useCGBasis>;
      return Teuchos::rcp( new HostBasisType(polyOrder_, pointType_) );
    }
  };
} // end namespace Intrepid2

// do ETI with default (double) type
extern template class Intrepid2::HierarchicalBasis_HDIV_PYR<Kokkos::DefaultExecutionSpace::device_type,double,double>;

#endif /* Intrepid2_HierarchicalBasis_HDIV_PYR_h */
