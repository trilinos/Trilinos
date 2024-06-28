// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_PyramidCoords.hpp
    \brief  Defines several coordinates and their gradients on the pyramid; maps from Intrepid2 (shards) pyramid to the ESEAS pyramid and back.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_Intrepid2_PyramidCoords_h
#define Intrepid2_Intrepid2_PyramidCoords_h

#include <Kokkos_DynRankView.hpp>

#include <Intrepid2_config.h>

namespace Intrepid2
{
  /// Compute various affine-like coordinates on the pyramid.  See Fuentes et al, Appendix E.9 for definitions.
  template<class PointScalar>
  KOKKOS_INLINE_FUNCTION
  void affinePyramid(Kokkos::Array<PointScalar,5>                                   &lambda,
                     Kokkos::Array<Kokkos::Array<PointScalar,3>,5>                  &lambdaGrad,
                     Kokkos::Array<Kokkos::Array<PointScalar,3>,2>                  &mu,
                     Kokkos::Array<Kokkos::Array<Kokkos::Array<PointScalar,3>,3>,2> &muGrad,
                     Kokkos::Array<Kokkos::Array<PointScalar,2>,3>                  &nu,
                     Kokkos::Array<Kokkos::Array<Kokkos::Array<PointScalar,3>,2>,3> &nuGrad,
                     Kokkos::Array<PointScalar,3>  &coords)
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

  /// Transforms from the Intrepid2 pyramid, centered at the origin with base [-1,1]^2 and height 1, to ESEAS pyramid, with base [0,1]^2, height 1, with its top vertex at (0,0,1).
  template<class PointScalar>
  KOKKOS_INLINE_FUNCTION
  void transformToESEASPyramid(      PointScalar &x_eseas,       PointScalar &y_eseas,       PointScalar &z_eseas,
                               const PointScalar &x_int2,  const PointScalar &y_int2,  const PointScalar &z_int2)
  {
    x_eseas = (x_int2 + 1. - z_int2) / 2.;
    y_eseas = (y_int2 + 1. - z_int2) / 2.;
    z_eseas = z_int2;
  }

  /// Transforms gradients computed on the ESEAS pyramid to gradients on the Intrepid2 pyramid.
  template<class OutputScalar>
  KOKKOS_INLINE_FUNCTION
  void transformFromESEASPyramidGradient(      OutputScalar &dx_int2,        OutputScalar &dy_int2,        OutputScalar &dz_int2,
                                         const OutputScalar &dx_eseas, const OutputScalar &dy_eseas, const OutputScalar &dz_eseas)
  {
    dx_int2 = dx_eseas / 2.;
    dy_int2 = dy_eseas / 2.;
    dz_int2 = dz_eseas - dx_int2 - dy_int2;
  }
} // end namespace Intrepid2

#endif /* Intrepid2_Intrepid2_PyramidCoords_h */
