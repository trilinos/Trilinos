// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#include <math.h>

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
    
    //! See Fuentes et al. (p. 455), definition of V_{ij}^{\trianglelefteq}
    KOKKOS_INLINE_FUNCTION
    void V_LEFT_TRI(Kokkos::Array<OutputScalar,3> &VLEFTTRI,
                    const OutputScalar &phi_i, const Kokkos::Array<OutputScalar,3> &phi_i_grad,
                    const OutputScalar &phi_j, const Kokkos::Array<OutputScalar,3> &phi_j_grad,
                    const OutputScalar &t0,    const Kokkos::Array<OutputScalar,3> &t0_grad) const {
      cross(VLEFTTRI, phi_i_grad, phi_j_grad);
      const OutputScalar t0_2 = t0 * t0;
      for (ordinal_type d=0; d<3; d++)
      {
        VLEFTTRI[d] *= t0_2;
      }
      
      Kokkos::Array<OutputScalar,3> tmp_t0_grad_t0, tmp_diff, tmp_cross;
      for (ordinal_type d=0; d<3; d++)
      {
        tmp_t0_grad_t0[d] = t0 * t0_grad[d];
        tmp_diff[d] = phi_i * phi_j_grad[d] - phi_j * phi_i_grad[d];
      }
      cross(tmp_cross, tmp_t0_grad_t0, tmp_diff);
      
      for (ordinal_type d=0; d<3; d++)
      {
        VLEFTTRI[d] += tmp_cross[d];
      }
    };
    
    
    //! See Fuentes et al. (p. 455), definition of V_{ij}^{\trianglerighteq}
    KOKKOS_INLINE_FUNCTION
    void V_RIGHT_TRI(Kokkos::Array<OutputScalar,3> &VRIGHTTRI,
                     const OutputScalar &mu1,    const Kokkos::Array<OutputScalar,3> &mu1_grad,
                     const OutputScalar &phi_i,  const Kokkos::Array<OutputScalar,3> &phi_i_grad,
                     const OutputScalar &t0,     const Kokkos::Array<OutputScalar,3> &t0_grad) const {
      Kokkos::Array<OutputScalar,3> left_vector; // left vector in the cross product we take below.
      
      const OutputScalar t0_2 = t0 * t0;
      for (ordinal_type d=0; d<3; d++)
      {
        left_vector[d] = t0_2 * phi_i_grad[d] + 2. * t0 * phi_i * t0_grad[d];
      }
      cross(VRIGHTTRI, left_vector, mu1_grad);
    };
    
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
    
    //! See Equation (B.42) in Fuentes et al.  Computes 1/mu V^{tri}_ij(mu s0, mu s1, s2).
    KOKKOS_INLINE_FUNCTION
    void V_TRI_B42(Kokkos::Array<OutputScalar,3> &VTRI_mus0_mus1_s2_over_mu,
                   const Kokkos::Array<OutputScalar,3> &VTRI_00_s0_s1_s2,
                   const Kokkos::Array<OutputScalar,3> &EE_0_s0_s1,
                   const OutputScalar &s2,
                   const OutputScalar &mu, const Kokkos::Array<OutputScalar,3> &mu_grad,
                   const ordinal_type &i, // i >= 0
                   const ordinal_type &j, // j >= 0
                   const OutputScratchView &P_mus0_mus1,      // container in which shiftedScaledLegendreValues have been computed for the appropriate face, with arguments (mu s0, mu s1)
                   const OutputScratchView &P_2ip1_mus0pmus1_s2 // container in which shiftedScaledJacobiValues have been computed for (2i+1) for the appropriate face, with arguments (mu s0 + mu s1, s2)
                   ) const
    {
      const auto &Pi_mus0_mus1         = P_mus0_mus1(i);
      const auto &Pj_2ip1_mus0pmus1_s2 = P_2ip1_mus0pmus1_s2(j);
      
      // place s2 (grad mu) x E^E_0 in result container
      cross(VTRI_mus0_mus1_s2_over_mu, mu_grad, EE_0_s0_s1);
      for (ordinal_type d=0; d<3; d++)
      {
        VTRI_mus0_mus1_s2_over_mu[d] *= s2;
      }
      
      // add mu V^{tri}_00 to it:
      for (ordinal_type d=0; d<3; d++)
      {
        VTRI_mus0_mus1_s2_over_mu[d] += mu * VTRI_00_s0_s1_s2[d];
      }
      
      // multiply by [P_i, P^{2i+1}_j]
      for (ordinal_type d=0; d<3; d++)
      {
        VTRI_mus0_mus1_s2_over_mu[d] *= Pi_mus0_mus1 * Pj_2ip1_mus0pmus1_s2;
      }
    }
    
    //! See Equation (B.42) in Fuentes et al.  Computes 1/mu V^{tri}_ij(mu s0, mu s1, s2).
    KOKKOS_INLINE_FUNCTION
    void V_TRI_B42_DIV(OutputScalar &div_VTRI_mus0_mus1_s2_over_mu,
                       const Kokkos::Array<OutputScalar,3> &VTRI_00_s0_s1_s2,
                       const Kokkos::Array<OutputScalar,3> &EE_0_s0_s1,
                       const OutputScalar &s2, const Kokkos::Array<OutputScalar,3> &s2_grad,
                       const OutputScalar &mu, const Kokkos::Array<OutputScalar,3> &mu_grad,
                       const ordinal_type &i, // i >= 0
                       const ordinal_type &j, // j >= 0
                       const OutputScratchView &P_mus0_mus1,      // container in which shiftedScaledLegendreValues have been computed for the appropriate face, with arguments (mu s0, mu s1)
                       const OutputScratchView &P_2ip1_mus0pmus1_s2 // container in which shiftedScaledJacobiValues have been computed for (2i+1) for the appropriate face, with arguments (mu s0 + mu s1, s2)
                       ) const
    {
      const auto &Pi_mus0_mus1         = P_mus0_mus1(i);
      const auto &Pj_2ip1_mus0pmus1_s2 = P_2ip1_mus0pmus1_s2(j);
      
      Kokkos::Array<OutputScalar,3> vector;
      
      // place E^E_0 x grad s2 in scratch vector container
      cross(vector, EE_0_s0_s1, s2_grad);
      // multiply by (i+j+3)
      for (ordinal_type d=0; d<3; d++)
      {
        vector[d] *= i+j+3.;
      }
      // subtract V_00
      for (ordinal_type d=0; d<3; d++)
      {
        vector[d] -= VTRI_00_s0_s1_s2[d];
      }
      
      OutputScalar dot_product;
      dot(dot_product, mu_grad, vector);
      
      div_VTRI_mus0_mus1_s2_over_mu = Pi_mus0_mus1 * Pj_2ip1_mus0pmus1_s2 * dot_product;
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
      auto pointOrdinal = teamMember.league_rank();
      OutputScratchView scratch1D_1, scratch1D_2, scratch1D_3;
      OutputScratchView scratch1D_4, scratch1D_5, scratch1D_6;
      OutputScratchView scratch1D_7, scratch1D_8, scratch1D_9;
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
      }
      
      const auto & x = inputPoints_(pointOrdinal,0);
      const auto & y = inputPoints_(pointOrdinal,1);
      const auto & z = inputPoints_(pointOrdinal,2);
      
      // Intrepid2 uses (-1,1)^2 for x,y
      // ESEAS uses (0,1)^2
      
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
          ordinal_type fieldOrdinalOffset = 0;
          // quadrilateral face
          {
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
                }
                fieldOrdinalOffset++;
              }
            }
          }
          
          // triangle faces
          {
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
              
              const auto & s0 = nu[0][a-1]; const auto & s0_grad = nuGrad[0][a-1];
              const auto & s1 = nu[1][a-1]; const auto & s1_grad = nuGrad[1][a-1];
              const auto & s2 = nu[2][a-1];
              const PointScalar jacobiScaling = s0 + s1 + s2;
              
              const PointScalar legendreScaling = s0 + s1;
              Polynomials::shiftedScaledLegendreValues(P, polyOrder_-1, s1, legendreScaling);
              
              const auto lambda0_index = tri_face_vertex_0[faceOrdinal];
              const auto lambda1_index = tri_face_vertex_1[faceOrdinal];
              const auto lambda2_index = tri_face_vertex_2[faceOrdinal];
              
              const auto & mu_c_b      = mu    [c][b-1];
              const auto & mu_c_b_grad = muGrad[c][b-1];
              const auto & mu_s0 = lambda[lambda0_index];
              const auto & mu_s1 = lambda[lambda1_index];
              const auto & mu_s2 = lambda[lambda2_index];
              
              const PointScalar muJacobiScaling = mu_s0 + mu_s1 + mu_s2;
              
              const PointScalar muLegendreScaling = mu_s0 + mu_s1;
              Polynomials::shiftedScaledLegendreValues(Pmu, polyOrder_-1, mu_s1, muLegendreScaling);
              
              Kokkos::Array<PointScalar, 3> vectorWeight;
              computeFaceVectorWeight(vectorWeight, nu[0][a-1], nuGrad[0][a-1],
                                                    nu[1][a-1], nuGrad[1][a-1],
                                                    nu[2][a-1], nuGrad[2][a-1]);
              
              Kokkos::Array<OutputScalar,3> VTRI_00;
              V_TRI(VTRI_00,0,0,P,P_2ip1,vectorWeight);
              
              Kokkos::Array<OutputScalar,3> EE_0;
              E_E(EE_0, 0, P, s0, s1, s0_grad, s1_grad);
              
              for (int totalPolyOrder=0; totalPolyOrder<polyOrder_; totalPolyOrder++)
              {
                for (int i=0; i<=totalPolyOrder; i++)
                {
                  const int j = totalPolyOrder - i;
                  
                  const double alpha = i*2.0 + 1;
                  Polynomials::shiftedScaledJacobiValues(  P_2ip1, alpha, polyOrder_-1,    s2,   jacobiScaling);
                  Polynomials::shiftedScaledJacobiValues(Pmu_2ip1, alpha, polyOrder_-1, mu_s2, muJacobiScaling);
                  
                  Kokkos::Array<OutputScalar,3> VTRI; // output from V_TRI
                  V_TRI(VTRI,    i, j, P,   P_2ip1,     vectorWeight);
                  
                  Kokkos::Array<OutputScalar,3> one_over_mu_VTRI_mu; // (B.42) result
                  V_TRI_B42(one_over_mu_VTRI_mu, VTRI_00, EE_0, s2, mu_c_b, mu_c_b_grad, i, j, Pmu, Pmu_2ip1);
                  
                  for (ordinal_type d=0; d<3; d++)
                  {
                    output_(fieldOrdinalOffset,pointOrdinal,d) = 0.5 * (VTRI[d] * mu_c_b + one_over_mu_VTRI_mu[d]);
                  }

                  fieldOrdinalOffset++;
                }
              }
            }
          } // triangle faces block
          
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
                
                // for reasons that I don't entirely understand, ESEAS switches whether i or j are the fastest-moving indices depending on whether it's family I or II.  Following their code, I'm calling the outer loop variable jg, inner ig.
                // (Cross-code comparisons are considerably simpler if we number the dofs in the same way.)
                ordinal_type jg_min = (f==0) ? 2 : 0;
                ordinal_type jg_max = (f==0) ? polyOrder_ : polyOrder_-1;
                ordinal_type ig_min = (f==0) ? 0 : 2;
                ordinal_type ig_max = (f==0) ? polyOrder_ -1 : polyOrder_;
                for (ordinal_type jg=jg_min; jg<=jg_max; jg++)
                {
                  for (ordinal_type ig=ig_min; ig<=ig_max; ig++)
                  {
                    const ordinal_type &i = (f==0) ? ig : jg;
                    const ordinal_type &j = (f==0) ? jg : ig;
                    Kokkos::Array<OutputScalar,3> EQUAD_ij;
                    Kokkos::Array<OutputScalar,3> curl_EQUAD_ij;
                    
                    E_QUAD(EQUAD_ij, i, j, Pi_s01, s0, s1, s0_grad, s1_grad, Li_t01);
                    
                    E_QUAD_CURL(curl_EQUAD_ij, i, j, Pi_s01, s0, s1, s0_grad, s1_grad,
                                                     Pi_t01, Li_t01, Li_dt_t01, t0_grad, t1_grad);
                    
                    // first term: muZ_0 phi^E_k curl EQUAD
                    // we can reuse the memory for curl_EQUAD_ij; we won't need the values there again
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
              // phi_ij_QUAD: phi_i(mu_X01) * phi_j(mu_Y01) // (following the ESEAS *implementation*; Fuentes et al. (p. 454) actually have the arguments reversed, which leads to a different basis ordering)
              const auto & phi_j = Li_muY01(j);
              Kokkos::Array<OutputScalar,3> phi_j_grad;
              computeGradHomLi(phi_j_grad, j, Pi_muY01, Li_dt_muY01, muY_0_grad, muY_1_grad);
              for (int i=2; i<=polyOrder_; i++)
              {
                const auto & phi_i = Li_muX01(i);
                Kokkos::Array<OutputScalar,3> phi_i_grad;
                computeGradHomLi(phi_i_grad, i, Pi_muX01, Li_dt_muX01, muX_0_grad, muX_1_grad);
                
                Kokkos::Array<OutputScalar,3> phi_ij_grad;
                for (ordinal_type d=0; d<3; d++)
                {
                  phi_ij_grad[d] = phi_i * phi_j_grad[d] + phi_j * phi_i_grad[d];
                }
                
                Kokkos::Array<OutputScalar,3> cross_product; // phi_ij_grad x grad_muZ0
                cross(cross_product, phi_ij_grad, muZ_0_grad);
                
                ordinal_type n = max(i,j);
                OutputScalar weight = n * pow(muZ_0,n-1);
                for (ordinal_type d=0; d<3; d++)
                {
                  output_(fieldOrdinalOffset,pointOrdinal,d) = weight * cross_product[d];
                }
                fieldOrdinalOffset++;
              }
            }
            
            // FAMILY IV (non-trivial divergences)
            {
              const auto muZ_0_squared = muZ_0 * muZ_0;
              const auto &s0 = muX_0;  const auto & s0_grad = muX_0_grad;
              const auto &s1 = muX_1;  const auto & s1_grad = muX_1_grad;
              const auto &t0 = muY_0;  const auto & t0_grad = muY_0_grad;
              const auto &t1 = muY_1;  const auto & t1_grad = muY_1_grad;
              const auto &Pi = Pi_muX01;
              const auto &Pj = Pi_muY01;
              for (int k=2; k<=polyOrder_; k++)
              {
                const auto & phi_k = Li_muZ01(k);
                for (int j=0; j<polyOrder_; j++)
                {
                  for (int i=0; i<polyOrder_; i++)
                  {
                    Kokkos::Array<OutputScalar,3> VQUAD; // output from V_QUAD
                    V_QUAD(VQUAD, i, j,
                           Pi, s0, s1, s0_grad, s1_grad,
                           Pj, t0, t1, t0_grad, t1_grad);
                    
                    for (int d=0; d<3; d++)
                    {
                      output_(fieldOrdinalOffset,pointOrdinal,d) = muZ_0_squared * phi_k * VQUAD[d];
                    }
                    
                    fieldOrdinalOffset++;
                  }
                }
              }
            }
            
            // FAMILY V (non-trivial divergences)
            {
              for (int j=2; j<=polyOrder_; j++)
              {
                const auto & phi_j = Li_muY01(j);
                Kokkos::Array<OutputScalar,3> phi_j_grad;
                computeGradHomLi(phi_j_grad, j, Pi_muY01, Li_dt_muY01, muY_0_grad, muY_1_grad);
                
                for (int i=2; i<=polyOrder_; i++)
                {
                  const auto & phi_i = Li_muX01(i);
                  Kokkos::Array<OutputScalar,3> phi_i_grad;
                  computeGradHomLi(phi_i_grad, i, Pi_muX01, Li_dt_muX01, muX_0_grad, muX_1_grad);
                  
                  const int n = max(i,j);
                  const OutputScalar muZ_1_nm1 = pow(muZ_1,n-1);
                  
                  Kokkos::Array<OutputScalar,3> VLEFTTRI;
                  V_LEFT_TRI(VLEFTTRI, phi_i, phi_i_grad, phi_j, phi_j_grad, muZ_0, muZ_0_grad);
                  
                  for (int d=0; d<3; d++)
                  {
                    output_(fieldOrdinalOffset,pointOrdinal,d) = muZ_1_nm1 * VLEFTTRI[d];
                  }
                  
                  fieldOrdinalOffset++;
                }
              }
            }
            
            // FAMILY VI (non-trivial divergences)
            for (int i=2; i<=polyOrder_; i++)
            {
              const auto & phi_i = Li_muX01(i);
              Kokkos::Array<OutputScalar,3> phi_i_grad;
              computeGradHomLi(phi_i_grad, i, Pi_muX01, Li_dt_muX01, muX_0_grad, muX_1_grad);
              
              Kokkos::Array<OutputScalar,3> VRIGHTTRI;
              V_RIGHT_TRI(VRIGHTTRI, muY_1, muY_1_grad, phi_i, phi_i_grad, muZ_0, muZ_0_grad);
              
              const OutputScalar muZ_1_im1 = pow(muZ_1,i-1);
              
              for (int d=0; d<3; d++)
              {
                output_(fieldOrdinalOffset,pointOrdinal,d) = muZ_1_im1 * VRIGHTTRI[d];
              }
              
              fieldOrdinalOffset++;
            }
            
            // FAMILY VII (non-trivial divergences)
            for (int j=2; j<=polyOrder_; j++)
            {
              const auto & phi_j = Li_muY01(j);
              Kokkos::Array<OutputScalar,3> phi_j_grad;
              computeGradHomLi(phi_j_grad, j, Pi_muY01, Li_dt_muY01, muY_0_grad, muY_1_grad);
              
              Kokkos::Array<OutputScalar,3> VRIGHTTRI;
              V_RIGHT_TRI(VRIGHTTRI, muX_1, muX_1_grad, phi_j, phi_j_grad, muZ_0, muZ_0_grad);
              
              const OutputScalar muZ_1_jm1 = pow(muZ_1,j-1);
              
              for (int d=0; d<3; d++)
              {
                output_(fieldOrdinalOffset,pointOrdinal,d) = muZ_1_jm1 * VRIGHTTRI[d];
              }
              
              fieldOrdinalOffset++;
            }
          }
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
              
              const auto & s0 = nu[0][a-1]; const auto & s0_grad = nuGrad[0][a-1];
              const auto & s1 = nu[1][a-1]; const auto & s1_grad = nuGrad[1][a-1];
              const auto & s2 = nu[2][a-1]; const auto & s2_grad = nuGrad[2][a-1];
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
              computeFaceVectorWeight(vectorWeight, nu[0][a-1], nuGrad[0][a-1],
                                                    nu[1][a-1], nuGrad[1][a-1],
                                                    nu[2][a-1], nuGrad[2][a-1]);
              
              Kokkos::Array<PointScalar,3> & mu_s0_grad = lambdaGrad[lambda0_index];
              Kokkos::Array<PointScalar,3> & mu_s1_grad = lambdaGrad[lambda1_index];
              Kokkos::Array<PointScalar,3> & mu_s2_grad = lambdaGrad[lambda2_index]; // == s2_grad
              
              Kokkos::Array<PointScalar, 3> muVectorWeight;
              computeFaceVectorWeight(muVectorWeight, mu_s0, mu_s0_grad,
                                                      mu_s1, mu_s1_grad,
                                                      mu_s2, mu_s2_grad);
              
              OutputScalar muDivWeight;
              computeFaceDivWeight(muDivWeight, mu_s0_grad, mu_s1_grad, mu_s2_grad);
              
              Kokkos::Array<OutputScalar,3> VTRI_00;
              V_TRI(VTRI_00,0,0,P,P_2ip1,vectorWeight);
              
              Kokkos::Array<OutputScalar,3> EE_0;
              E_E(EE_0, 0, P, s0, s1, s0_grad, s1_grad);
              
              for (int totalPolyOrder=0; totalPolyOrder<polyOrder_; totalPolyOrder++)
              {
                for (int i=0; i<=totalPolyOrder; i++)
                {
                  const int j = totalPolyOrder - i;
                  
                  const double alpha = i*2.0 + 1;
                  Polynomials::shiftedScaledJacobiValues(  P_2ip1, alpha, polyOrder_-1,    s2,   jacobiScaling);
                  Polynomials::shiftedScaledJacobiValues(Pmu_2ip1, alpha, polyOrder_-1, mu_s2, muJacobiScaling);
                  
                  Kokkos::Array<OutputScalar,3> VTRI; // output from V_TRI
                  
                  V_TRI(VTRI, i, j, P, P_2ip1, vectorWeight);
                  
                  OutputScalar div_one_over_mu_VTRI_mu;
                  V_TRI_B42_DIV(div_one_over_mu_VTRI_mu, VTRI_00, EE_0, s2, s2_grad, mu_c_b, mu_c_b_grad, i, j, Pmu, Pmu_2ip1);
                  
                  output_(fieldOrdinalOffset,pointOrdinal) = 0.5 * (dot(mu_c_b_grad, VTRI) + div_one_over_mu_VTRI_mu);
                  
                  fieldOrdinalOffset++;
                }
              }
            }
          } // end triangle face block
          
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
            
            const auto & Li_muZ01    = scratch1D_1; // used in Family IV
            const auto & Li_muX01    = scratch1D_2; // used in Family V
            const auto & Li_muY01    = scratch1D_3; // used in Family V
            const auto & Pi_muX01    = scratch1D_4; // used in Family IV
            const auto & Pi_muY01    = scratch1D_5; // used in Family IV
            const auto & Pi_muZ01    = scratch1D_6; // used in Family IV
            const auto & Li_dt_muX01 = scratch1D_7; // used in Family V
            const auto & Li_dt_muY01 = scratch1D_8; // used in Family V
            const auto & Li_dt_muZ01 = scratch1D_9; // used in Family IV
            
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
            
            // FAMILY IV -- non-trivial divergences
            {
              const auto muZ_0_squared = muZ_0 * muZ_0;
              const auto &s0 = muX_0;  const auto & s0_grad = muX_0_grad;
              const auto &s1 = muX_1;  const auto & s1_grad = muX_1_grad;
              const auto &t0 = muY_0;  const auto & t0_grad = muY_0_grad;
              const auto &t1 = muY_1;  const auto & t1_grad = muY_1_grad;
              const auto &Pi = Pi_muX01;
              const auto &Pj = Pi_muY01;
              
              for (int k=2; k<=polyOrder_; k++)
              {
                const auto & phi_k = Li_muZ01(k);
                Kokkos::Array<OutputScalar,3> phi_k_grad;
                computeGradHomLi(phi_k_grad, k, Pi_muZ01, Li_dt_muZ01, muZ_0_grad, muZ_1_grad);
                for (int j=0; j<polyOrder_; j++)
                {
                  for (int i=0; i<polyOrder_; i++)
                  {
                    Kokkos::Array<OutputScalar,3> firstTerm{0,0,0}; // muZ_0_squared * grad phi_k + 2 muZ_0 * phi_k * grad muZ_0
                    for (int d=0; d<3; d++)
                    {
                      firstTerm[d] += muZ_0_squared * phi_k_grad[d] + 2. * muZ_0 * phi_k * muZ_0_grad[d];
                    }
                    Kokkos::Array<OutputScalar,3> VQUAD; // output from V_QUAD
                    V_QUAD(VQUAD, i, j,
                           Pi, s0, s1, s0_grad, s1_grad,
                           Pj, t0, t1, t0_grad, t1_grad);
                    
                    OutputScalar divValue;
                    dot(divValue, firstTerm, VQUAD);
                    output_(fieldOrdinalOffset,pointOrdinal) = divValue;
                    
                    fieldOrdinalOffset++;
                  }
                }
              }
            }
            
            // FAMILY V -- non-trivial divergences
            {
              for (int j=2; j<=polyOrder_; j++)
              {
                const auto & phi_j = Li_muX01(j);
                Kokkos::Array<OutputScalar,3> phi_j_grad;
                computeGradHomLi(phi_j_grad, j, Pi_muY01, Li_dt_muY01, muY_0_grad, muY_1_grad);
                
                for (int i=2; i<=polyOrder_; i++)
                {
                  const auto & phi_i = Li_muY01(i);
                  Kokkos::Array<OutputScalar,3> phi_i_grad;
                  computeGradHomLi(phi_i_grad, i, Pi_muX01, Li_dt_muX01, muX_0_grad, muX_1_grad);
                  
                  const int n = max(i,j);
                  const OutputScalar muZ_1_nm2 = pow(muZ_1,n-2);
                  
                  Kokkos::Array<OutputScalar,3> VLEFTTRI;
                  V_LEFT_TRI(VLEFTTRI, phi_i, phi_i_grad, phi_j, phi_j_grad, muZ_0, muZ_0_grad);
                  
                  OutputScalar dot_product;
                  dot(dot_product, muZ_1_grad, VLEFTTRI);
                  output_(fieldOrdinalOffset,pointOrdinal) = (n-1) * muZ_1_nm2 * dot_product;
                  
                  fieldOrdinalOffset++;
                }
              }
            }
            
            // FAMILY VI (non-trivial divergences)
            for (int i=2; i<=polyOrder_; i++)
            {
              const auto & phi_i = Li_muX01(i);
              Kokkos::Array<OutputScalar,3> phi_i_grad;
              computeGradHomLi(phi_i_grad, i, Pi_muX01, Li_dt_muX01, muX_0_grad, muX_1_grad);
              
              Kokkos::Array<OutputScalar,3> VRIGHTTRI;
              V_RIGHT_TRI(VRIGHTTRI, muY_1, muY_1_grad, phi_i, phi_i_grad, muZ_0, muZ_0_grad);
              
              OutputScalar dot_product;
              dot(dot_product, muZ_1_grad, VRIGHTTRI);
              
              const OutputScalar muZ_1_im2 = pow(muZ_1,i-2);
              output_(fieldOrdinalOffset,pointOrdinal) = (i-1) * muZ_1_im2 * dot_product;
              
              fieldOrdinalOffset++;
            }
            
            // FAMILY VII (non-trivial divergences)
            for (int j=2; j<=polyOrder_; j++)
            {
              const auto & phi_j = Li_muY01(j);
              Kokkos::Array<OutputScalar,3> phi_j_grad;
              computeGradHomLi(phi_j_grad, j, Pi_muY01, Li_dt_muY01, muY_0_grad, muY_1_grad);
              
              Kokkos::Array<OutputScalar,3> VRIGHTTRI;
              V_RIGHT_TRI(VRIGHTTRI, muX_1, muX_1_grad, phi_j, phi_j_grad, muZ_0, muZ_0_grad);
              
              OutputScalar dot_product;
              dot(dot_product, muZ_1_grad, VRIGHTTRI);
              
              const OutputScalar muZ_1_jm2 = pow(muZ_1,j-2);
              output_(fieldOrdinalOffset,pointOrdinal) = (j-1) * muZ_1_jm2 * dot_product;
              
              fieldOrdinalOffset++;
            }
            
          } // end interior function block
          
        } // end OPERATOR_DIV block
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
      size_t shmem_size = 0;
      if (fad_size_output_ > 0)
      {
        // 1D scratch views (we have up to 9 of them):
        shmem_size += 9 * OutputScratchView::shmem_size(polyOrder_ + 1, fad_size_output_);
      }
      else
      {
        // 1D scratch views (we have up to 9 of them):
        shmem_size += 9 * OutputScratchView::shmem_size(polyOrder_ + 1);
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
      const auto & p              = polyOrder;
      this->basisCardinality_     = p * p + 2 * p * (p+1) + 3 * p * p * (p-1);
      this->basisDegree_          = p;
      this->basisCellTopologyKey_ = shards::Pyramid<>::key;
      this->basisType_            = BASIS_FEM_HIERARCHICAL;
      this->basisCoordinates_     = COORDINATES_CARTESIAN;
      this->functionSpace_        = FUNCTION_SPACE_HDIV;
      
      const int degreeLength = 1;
      this->fieldOrdinalPolynomialDegree_ = OrdinalTypeArray2DHost("Integrated Legendre H(div) pyramid polynomial degree lookup", this->basisCardinality_, degreeLength);
      this->fieldOrdinalH1PolynomialDegree_ = OrdinalTypeArray2DHost("Integrated Legendre H(div) pyramid polynomial H^1 degree lookup", this->basisCardinality_, degreeLength);
      
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
        // ESEAS reverses i,j loop ordering for family II, relative to family I.
        // We follow ESEAS for convenience of cross-code comparison.
        for (int i=0; i<polyOrder_; i++)
        {
          for (int j=2; j<=polyOrder_; j++)
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
