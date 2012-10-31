/*
//@HEADER
// ************************************************************************
//
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_HEXEXPLICITFUNCTIONS_HPP
#define KOKKOSARRAY_HEXEXPLICITFUNCTIONS_HPP

#include <math.h>

namespace Explicit {

template < class Device >
struct Hex8Functions
{
  static const unsigned SpatialDim    = 3 ;
  static const unsigned ElemNodeCount = 8 ;

  // Indices for full 3x3 tensor:

  static const unsigned K_F_XX = 0 ;
  static const unsigned K_F_YY = 1 ;
  static const unsigned K_F_ZZ = 2 ;
  static const unsigned K_F_XY = 3 ;
  static const unsigned K_F_YZ = 4 ;
  static const unsigned K_F_ZX = 5 ;
  static const unsigned K_F_YX = 6 ;
  static const unsigned K_F_ZY = 7 ;
  static const unsigned K_F_XZ = 8 ;
  static const unsigned K_F_SIZE = 9 ;

  //  Indexes into a 3 by 3 symmetric tensor stored as a length 6 vector

  static const unsigned K_S_XX = 0 ;
  static const unsigned K_S_YY = 1 ;
  static const unsigned K_S_ZZ = 2 ;
  static const unsigned K_S_XY = 3 ;
  static const unsigned K_S_YZ = 4 ;
  static const unsigned K_S_ZX = 5 ;
  static const unsigned K_S_YX = 3 ;
  static const unsigned K_S_ZY = 4 ;
  static const unsigned K_S_XZ = 5 ;
  static const unsigned K_S_SIZE = 6 ;

  //  Indexes into a 3 by 3 skew symmetric tensor stored as a length 3 vector

  static const unsigned K_V_XY = 0 ;
  static const unsigned K_V_YZ = 1 ;
  static const unsigned K_V_ZX = 2 ;
  static const unsigned K_V_SIZE = 3 ;

  //--------------------------------------------------------------------------

  template< typename Scalar >
  static KOKKOSARRAY_INLINE_FUNCTION
  Scalar dot8( const Scalar * a , const Scalar * b )
  { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3] +
           a[4] * b[4] + a[5] * b[5] + a[6] * b[6] + a[7] * b[7] ; }

  //--------------------------------------------------------------------------

  template< typename CoordScalarType , typename Scalar >
  static KOKKOSARRAY_INLINE_FUNCTION
  void grad( const CoordScalarType x[] ,
             const CoordScalarType z[] ,
                   Scalar grad_y[] )
  {
    const Scalar R42=(x[3] - x[1]);
    const Scalar R52=(x[4] - x[1]);
    const Scalar R54=(x[4] - x[3]);

    const Scalar R63=(x[5] - x[2]);
    const Scalar R83=(x[7] - x[2]);
    const Scalar R86=(x[7] - x[5]);

    const Scalar R31=(x[2] - x[0]);
    const Scalar R61=(x[5] - x[0]);
    const Scalar R74=(x[6] - x[3]);

    const Scalar R72=(x[6] - x[1]);
    const Scalar R75=(x[6] - x[4]);
    const Scalar R81=(x[7] - x[0]);

    const Scalar t1=(R63 + R54);
    const Scalar t2=(R61 + R74);
    const Scalar t3=(R72 + R81);

    const Scalar t4 =(R86 + R42);
    const Scalar t5 =(R83 + R52);
    const Scalar t6 =(R75 + R31);

    //  Calculate Y gradient from X and Z data

    grad_y[0] = (z[1] *  t1) - (z[2] * R42) - (z[3] *  t5)  + (z[4] *  t4) + (z[5] * R52) - (z[7] * R54);
    grad_y[1] = (z[2] *  t2) + (z[3] * R31) - (z[0] *  t1)  - (z[5] *  t6) + (z[6] * R63) - (z[4] * R61);
    grad_y[2] = (z[3] *  t3) + (z[0] * R42) - (z[1] *  t2)  - (z[6] *  t4) + (z[7] * R74) - (z[5] * R72);
    grad_y[3] = (z[0] *  t5) - (z[1] * R31) - (z[2] *  t3)  + (z[7] *  t6) + (z[4] * R81) - (z[6] * R83);
    grad_y[4] = (z[5] *  t3) + (z[6] * R86) - (z[7] *  t2)  - (z[0] *  t4) - (z[3] * R81) + (z[1] * R61);
    grad_y[5] = (z[6] *  t5) - (z[4] *  t3)  - (z[7] * R75) + (z[1] *  t6) - (z[0] * R52) + (z[2] * R72);
    grad_y[6] = (z[7] *  t1) - (z[5] *  t5)  - (z[4] * R86) + (z[2] *  t4) - (z[1] * R63) + (z[3] * R83);
    grad_y[7] = (z[4] *  t2) - (z[6] *  t1)  + (z[5] * R75) - (z[3] *  t6) - (z[2] * R74) + (z[0] * R54);
  }

  template< typename CoordScalarType , typename Scalar >
  static KOKKOSARRAY_INLINE_FUNCTION
  void grad( const CoordScalarType x[] ,
             const CoordScalarType y[] ,
             const CoordScalarType z[] ,
                   Scalar grad_x[] ,
                   Scalar grad_y[] ,
                   Scalar grad_z[] )
  {
    grad( x , z , grad_y );
    grad( z , y , grad_x );
    grad( y , x , grad_z );
  }

  //--------------------------------------------------------------------------

  template< typename Scalar >
  static KOKKOSARRAY_INLINE_FUNCTION
  void polar_decomp( const Scalar dt ,
                     const Scalar v_gr[] ,
                           Scalar stretch[] /* INOUT */ ,
                           Scalar str_ten[] /* OUT */ ,
                           Scalar rot[]     /* OUT */ )
  {
    const Scalar dt_half = 0.5 * dt;

    Scalar vort[ K_V_SIZE ];  // Vorticity

    //  Symmetric part
    str_ten[K_S_XX] = v_gr[K_F_XX];
    str_ten[K_S_YY] = v_gr[K_F_YY];
    str_ten[K_S_ZZ] = v_gr[K_F_ZZ];
    str_ten[K_S_XY] = 0.5 * ( v_gr[K_F_XY] + v_gr[K_F_YX] );
    str_ten[K_S_YZ] = 0.5 * ( v_gr[K_F_YZ] + v_gr[K_F_ZY] );
    str_ten[K_S_ZX] = 0.5 * ( v_gr[K_F_ZX] + v_gr[K_F_XZ] );

    //  Skew Symmetric part
    vort[K_V_XY] = 0.5 * ( v_gr[K_F_XY] - v_gr[K_F_YX] );
    vort[K_V_YZ] = 0.5 * ( v_gr[K_F_YZ] - v_gr[K_F_ZY] );
    vort[K_V_ZX] = 0.5 * ( v_gr[K_F_ZX] - v_gr[K_F_XZ] );

    //   calculate the rates of rotation via gauss elimination.

    Scalar z1 = str_ten[K_S_XY] * stretch[K_S_ZX] -
                str_ten[K_S_ZX] * stretch[K_S_XY] +
                str_ten[K_S_YY] * stretch[K_S_YZ] -
                str_ten[K_S_YZ] * stretch[K_S_YY] +
                str_ten[K_S_YZ] * stretch[K_S_ZZ] -
                str_ten[K_S_ZZ] * stretch[K_S_YZ];

    Scalar z2 = str_ten[K_S_ZX] * stretch[K_S_XX] -
                str_ten[K_S_XX] * stretch[K_S_ZX] +
                str_ten[K_S_YZ] * stretch[K_S_XY] -
                str_ten[K_S_XY] * stretch[K_S_YZ] +
                str_ten[K_S_ZZ] * stretch[K_S_ZX] -
                str_ten[K_S_ZX] * stretch[K_S_ZZ];

    Scalar z3 = str_ten[K_S_XX] * stretch[K_S_XY] -
                str_ten[K_S_XY] * stretch[K_S_XX] +
                str_ten[K_S_XY] * stretch[K_S_YY] -
                str_ten[K_S_YY] * stretch[K_S_XY] +
                str_ten[K_S_ZX] * stretch[K_S_YZ] -
                str_ten[K_S_YZ] * stretch[K_S_ZX];

    //   forward elimination

    const Scalar a1inv = 1.0 / (stretch[K_S_YY] + stretch[K_S_ZZ]);

    const Scalar a4BYa1 = -1 * stretch[K_S_XY] * a1inv;

    const Scalar a2inv = 1.0 / (stretch[K_S_ZZ] + stretch[K_S_XX] + stretch[K_S_XY] * a4BYa1);

    const Scalar a5 =  -stretch[K_S_YZ] + stretch[K_S_ZX] * a4BYa1;

    z2 -= z1 * a4BYa1;
    const Scalar a6BYa1 = -1 * stretch[K_S_ZX] * a1inv;
    const Scalar a5BYa2 = a5 * a2inv;
    z3 -= z1 * a6BYa1 - z2 * a5BYa2;

    //   backward substitution -

    z3 /= (stretch[K_S_XX] + stretch[K_S_YY] + stretch[K_S_ZX] * a6BYa1 + a5 * a5BYa2);
    z2 = (z2 - a5 * z3) * a2inv;
    z1 = (z1*a1inv - a6BYa1 * z3 -a4BYa1 * z2);

    //   calculate rotation rates - recall that spin_rate is an asymmetric tensor,
    //   so compute spin rate vector as dual of spin rate tensor,
    //   i.e   w_i = e_ijk * spin_rate_jk

    z1 += vort[K_V_YZ];
    z2 += vort[K_V_ZX];
    z3 += vort[K_V_XY];

    //   update rotation tensor:
    //  1) premultiply old rotation tensor to get right-hand side.

    Scalar r_XX = rot[K_F_XX] + dt_half*( z3 * rot[K_F_YX] - z2 * rot[K_F_ZX] );
    Scalar r_YX = rot[K_F_YX] + dt_half*( z1 * rot[K_F_ZX] - z3 * rot[K_F_XX] );
    Scalar r_ZX = rot[K_F_ZX] + dt_half*( z2 * rot[K_F_XX] - z1 * rot[K_F_YX] );
    Scalar r_XY = rot[K_F_XY] + dt_half*( z3 * rot[K_F_YY] - z2 * rot[K_F_ZY] );
    Scalar r_YY = rot[K_F_YY] + dt_half*( z1 * rot[K_F_ZY] - z3 * rot[K_F_XY] );
    Scalar r_ZY = rot[K_F_ZY] + dt_half*( z2 * rot[K_F_XY] - z1 * rot[K_F_YY] );
    Scalar r_XZ = rot[K_F_XZ] + dt_half*( z3 * rot[K_F_YZ] - z2 * rot[K_F_ZZ] );
    Scalar r_YZ = rot[K_F_YZ] + dt_half*( z1 * rot[K_F_ZZ] - z3 * rot[K_F_XZ] );
    Scalar r_ZZ = rot[K_F_ZZ] + dt_half*( z2 * rot[K_F_XZ] - z1 * rot[K_F_YZ] );


    //  2) solve for new rotation tensor via gauss elimination.
    //   forward elimination -

    Scalar a12 = - dt_half * z3;
    Scalar a13 =   dt_half * z2;
    Scalar b32 = - dt_half * z1;
    Scalar a22inv = 1.0 / (1.0 + a12 * a12);

    Scalar a13a12 = a13*a12;
    Scalar a23 = b32 + a13a12;
    r_YX += r_XX * a12;
    r_YY += r_XY * a12;
    r_YZ += r_XZ * a12;

    b32 = (b32 - a13a12) * a22inv;
    r_ZX += r_XX * a13 + r_YX * b32;
    r_ZY += r_XY * a13 + r_YY * b32;
    r_ZZ += r_XZ * a13 + r_YZ * b32;

    //   backward substitution -

    const Scalar a33inv = 1.0 / (1.0 + a13 * a13 + a23 * b32);

    rot[K_F_ZX] = r_ZX * a33inv;
    rot[K_F_ZY] = r_ZY * a33inv;
    rot[K_F_ZZ] = r_ZZ * a33inv;
    rot[K_F_YX] = ( r_YX - rot[K_F_ZX] * a23 ) * a22inv;
    rot[K_F_YY] = ( r_YY - rot[K_F_ZY] * a23 ) * a22inv;
    rot[K_F_YZ] = ( r_YZ - rot[K_F_ZZ] * a23 ) * a22inv;
    rot[K_F_XX] = r_XX - rot[K_F_ZX] * a13 - rot[K_F_YX] * a12;
    rot[K_F_XY] = r_XY - rot[K_F_ZY] * a13 - rot[K_F_YY] * a12;
    rot[K_F_XZ] = r_XZ - rot[K_F_ZZ] * a13 - rot[K_F_YZ] * a12;

    //   update stretch tensor in the new configuration -

    const Scalar a1 = str_ten[K_S_XY] + vort[K_V_XY];
    const Scalar a2 = str_ten[K_S_YZ] + vort[K_V_YZ];
    const Scalar a3 = str_ten[K_S_ZX] + vort[K_V_ZX];
    const Scalar b1 = str_ten[K_S_ZX] - vort[K_V_ZX];
    const Scalar b2 = str_ten[K_S_XY] - vort[K_V_XY];
    const Scalar b3 = str_ten[K_S_YZ] - vort[K_V_YZ];

    const Scalar s_XX = stretch[K_S_XX];
    const Scalar s_YY = stretch[K_S_YY];
    const Scalar s_ZZ = stretch[K_S_ZZ];
    const Scalar s_XY = stretch[K_S_XY];
    const Scalar s_YZ = stretch[K_S_YZ];
    const Scalar s_ZX = stretch[K_S_ZX];

    stretch[K_S_XX] += dt * (str_ten[K_S_XX] * s_XX + ( a1 + z3 ) * s_XY + ( b1 - z2 ) * s_ZX);
    stretch[K_S_YY] += dt * (str_ten[K_S_YY] * s_YY + ( a2 + z1 ) * s_YZ + ( b2 - z3 ) * s_XY);
    stretch[K_S_ZZ] += dt * (str_ten[K_S_ZZ] * s_ZZ + ( a3 + z2 ) * s_ZX + ( b3 - z1 ) * s_YZ);
    stretch[K_S_XY] += dt * (str_ten[K_S_XX] * s_XY + ( a1 )      * s_YY + ( b1      ) * s_YZ - z3 * s_XX + z1 * s_ZX);
    stretch[K_S_YZ] += dt * (str_ten[K_S_YY] * s_YZ + ( a2 )      * s_ZZ + ( b2      ) * s_ZX - z1 * s_YY + z2 * s_XY);
    stretch[K_S_ZX] += dt * (str_ten[K_S_ZZ] * s_ZX + ( a3 )      * s_XX + ( b3      ) * s_XY - z2 * s_ZZ + z3 * s_YZ);
  }

  //--------------------------------------------------------------------------

  template< typename Scalar >
  static KOKKOSARRAY_INLINE_FUNCTION
  void rotate_tensor( const Scalar str_ten[] ,
                      const Scalar rot[] ,
                            Scalar rot_str[] )
  {
    Scalar t[9];

    t[0] = str_ten[K_S_XX]*rot[K_F_XX] + str_ten[K_S_XY]*rot[K_F_YX] + str_ten[K_S_XZ]*rot[K_F_ZX];
    t[1] = str_ten[K_S_YX]*rot[K_F_XX] + str_ten[K_S_YY]*rot[K_F_YX] + str_ten[K_S_YZ]*rot[K_F_ZX];
    t[2] = str_ten[K_S_ZX]*rot[K_F_XX] + str_ten[K_S_ZY]*rot[K_F_YX] + str_ten[K_S_ZZ]*rot[K_F_ZX];

    t[3] = str_ten[K_S_XX]*rot[K_F_XY] + str_ten[K_S_XY]*rot[K_F_YY] + str_ten[K_S_XZ]*rot[K_F_ZY];
    t[4] = str_ten[K_S_YX]*rot[K_F_XY] + str_ten[K_S_YY]*rot[K_F_YY] + str_ten[K_S_YZ]*rot[K_F_ZY];
    t[5] = str_ten[K_S_ZX]*rot[K_F_XY] + str_ten[K_S_ZY]*rot[K_F_YY] + str_ten[K_S_ZZ]*rot[K_F_ZY];

    t[6] = str_ten[K_S_XX]*rot[K_F_XZ] + str_ten[K_S_XY]*rot[K_F_YZ] + str_ten[K_S_XZ]*rot[K_F_ZZ];
    t[7] = str_ten[K_S_YX]*rot[K_F_XZ] + str_ten[K_S_YY]*rot[K_F_YZ] + str_ten[K_S_YZ]*rot[K_F_ZZ];
    t[8] = str_ten[K_S_ZX]*rot[K_F_XZ] + str_ten[K_S_ZY]*rot[K_F_YZ] + str_ten[K_S_ZZ]*rot[K_F_ZZ];


    rot_str[ K_S_XX ] = rot[K_F_XX] * t[0] + rot[K_F_YX] * t[1] + rot[K_F_ZX] * t[2];
    rot_str[ K_S_YY ] = rot[K_F_XY] * t[3] + rot[K_F_YY] * t[4] + rot[K_F_ZY] * t[5];
    rot_str[ K_S_ZZ ] = rot[K_F_XZ] * t[6] + rot[K_F_YZ] * t[7] + rot[K_F_ZZ] * t[8];

    rot_str[ K_S_XY ] = rot[K_F_XX] * t[3] + rot[K_F_YX] * t[4] + rot[K_F_ZX] * t[5];
    rot_str[ K_S_YZ ] = rot[K_F_XY] * t[6] + rot[K_F_YY] * t[7] + rot[K_F_ZY] * t[8];
    rot_str[ K_S_ZX ] = rot[K_F_XZ] * t[0] + rot[K_F_YZ] * t[1] + rot[K_F_ZZ] * t[2];
  }

  //--------------------------------------------------------------------------

  template< typename Scalar >
  static KOKKOSARRAY_INLINE_FUNCTION
  void rotate_tensor_backward( const Scalar s_n[] ,
                               const Scalar r_n[] ,
                                     Scalar rot_stress[] )
  {
    //   t : temporary variables
    //   s_n : stress_new in local memory space
    //   r_n : rotation_new in local memory space
    Scalar t[9] ;

    t[0] = s_n[K_S_XX]*r_n[K_F_XX]+ s_n[K_S_XY]*r_n[K_F_XY]+ s_n[K_S_XZ]*r_n[K_F_XZ];
    t[1] = s_n[K_S_YX]*r_n[K_F_XX]+ s_n[K_S_YY]*r_n[K_F_XY]+ s_n[K_S_YZ]*r_n[K_F_XZ];
    t[2] = s_n[K_S_ZX]*r_n[K_F_XX]+ s_n[K_S_ZY]*r_n[K_F_XY]+ s_n[K_S_ZZ]*r_n[K_F_XZ];
    t[3] = s_n[K_S_XX]*r_n[K_F_YX]+ s_n[K_S_XY]*r_n[K_F_YY]+ s_n[K_S_XZ]*r_n[K_F_YZ];
    t[4] = s_n[K_S_YX]*r_n[K_F_YX]+ s_n[K_S_YY]*r_n[K_F_YY]+ s_n[K_S_YZ]*r_n[K_F_YZ];
    t[5] = s_n[K_S_ZX]*r_n[K_F_YX]+ s_n[K_S_ZY]*r_n[K_F_YY]+ s_n[K_S_ZZ]*r_n[K_F_YZ];
    t[6] = s_n[K_S_XX]*r_n[K_F_ZX]+ s_n[K_S_XY]*r_n[K_F_ZY]+ s_n[K_S_XZ]*r_n[K_F_ZZ];
    t[7] = s_n[K_S_YX]*r_n[K_F_ZX]+ s_n[K_S_YY]*r_n[K_F_ZY]+ s_n[K_S_YZ]*r_n[K_F_ZZ];
    t[8] = s_n[K_S_ZX]*r_n[K_F_ZX]+ s_n[K_S_ZY]*r_n[K_F_ZY]+ s_n[K_S_ZZ]*r_n[K_F_ZZ];

    rot_stress[ K_S_XX ] = r_n[K_F_XX]*t[0] + r_n[K_F_XY]*t[1] + r_n[K_F_XZ]*t[2];
    rot_stress[ K_S_YY ] = r_n[K_F_YX]*t[3] + r_n[K_F_YY]*t[4] + r_n[K_F_YZ]*t[5];
    rot_stress[ K_S_ZZ ] = r_n[K_F_ZX]*t[6] + r_n[K_F_ZY]*t[7] + r_n[K_F_ZZ]*t[8];

    rot_stress[ K_S_XY ] = r_n[K_F_XX]*t[3] + r_n[K_F_XY]*t[4] + r_n[K_F_XZ]*t[5];
    rot_stress[ K_S_YZ ] = r_n[K_F_YX]*t[6] + r_n[K_F_YY]*t[7] + r_n[K_F_YZ]*t[8];
    rot_stress[ K_S_ZX ] = r_n[K_F_ZX]*t[0] + r_n[K_F_ZY]*t[1] + r_n[K_F_ZZ]*t[2];
  }

  //--------------------------------------------------------------------------

  template< typename Scalar , typename ScalarStress >
  static KOKKOSARRAY_INLINE_FUNCTION
  void update_stress( const Scalar dt ,
                      const Scalar two_mu ,
                      const Scalar bulk_modulus ,
                      const Scalar rot_str[] ,
                      ScalarStress stress[] )
  {
    const Scalar e = rot_str[ K_S_XX ] + rot_str[ K_S_YY ] + rot_str[ K_S_ZZ ] ;
    const Scalar eb = e * bulk_modulus ;
    const Scalar e3 = e / 3.0 ;

    stress[K_S_XX] += dt * ( two_mu * ( rot_str[K_S_XX] - e3 ) + eb );
    stress[K_S_YY] += dt * ( two_mu * ( rot_str[K_S_YY] - e3 ) + eb );
    stress[K_S_ZZ] += dt * ( two_mu * ( rot_str[K_S_ZZ] - e3 ) + eb );

    stress[K_S_XY] += dt * two_mu * rot_str[K_S_XY];
    stress[K_S_YZ] += dt * two_mu * rot_str[K_S_YZ];
    stress[K_S_ZX] += dt * two_mu * rot_str[K_S_ZX];
  }

  //--------------------------------------------------------------------------

  template< typename Scalar >
  static KOKKOSARRAY_INLINE_FUNCTION
  void comp_force( const Scalar vx[] ,
                   const Scalar vy[] ,
                   const Scalar vz[] ,
                   const Scalar grad_x[] ,
                   const Scalar grad_y[] ,
                   const Scalar grad_z[] ,
                   const Scalar total_stress12th[] ,
                         Scalar force[][ SpatialDim ] ,
                         Scalar & energy )
  {
    Scalar internal_energy = 0 ;

    for ( unsigned inode = 0; inode < ElemNodeCount ; ++inode ) {

      force[inode][0] = total_stress12th[K_S_XX] * grad_x[inode] +
                        total_stress12th[K_S_XY] * grad_y[inode] +
                        total_stress12th[K_S_XZ] * grad_z[inode] ;

      force[inode][1] = total_stress12th[K_S_YX] * grad_x[inode] +
                        total_stress12th[K_S_YY] * grad_y[inode] +
                        total_stress12th[K_S_YZ] * grad_z[inode] ;

      force[inode][2] = total_stress12th[K_S_ZX] * grad_x[inode] +
                        total_stress12th[K_S_ZY] * grad_y[inode] +
                        total_stress12th[K_S_ZZ] * grad_z[inode] ;

      internal_energy += force[inode][0] * vx[inode] +
                         force[inode][1] * vy[inode] +
                         force[inode][2] * vz[inode] ;
    }

    energy = internal_energy ;
  }

  //--------------------------------------------------------------------------
};

} // namespace Explicit

#endif /* #ifndef KOKKOSARRAY_HEXEXPLICITFUNCTIONS_HPP */

