/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <math.h>

namespace Explicit {

//----------------------------------------------------------------------------

template<>
struct Hex8Functions< KOKKOS_MACRO_DEVICE >
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
  static inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  Scalar dot8( const Scalar * a , const Scalar * b )
  { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3] +
           a[4] * b[4] + a[5] * b[5] + a[6] * b[6] + a[7] * b[7] ; }

  //--------------------------------------------------------------------------

  template< typename CoordScalarType , typename Scalar >
  static inline
  KOKKOS_MACRO_DEVICE_FUNCTION
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
  static inline
  KOKKOS_MACRO_DEVICE_FUNCTION
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
  static inline
  KOKKOS_MACRO_DEVICE_FUNCTION
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
  static inline
  KOKKOS_MACRO_DEVICE_FUNCTION
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
  static inline
  KOKKOS_MACRO_DEVICE_FUNCTION
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
  static inline
  KOKKOS_MACRO_DEVICE_FUNCTION
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
  static inline
  KOKKOS_MACRO_DEVICE_FUNCTION
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
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename Scalar >
struct InitializeElement< Explicit::Fields< Scalar, KOKKOS_MACRO_DEVICE > >
{
  typedef KOKKOS_MACRO_DEVICE     device_type ;

  typedef Explicit::Fields< Scalar , device_type > Fields ;

  typedef Hex8Functions<device_type> ElemFunc ;

  const typename Fields::elem_node_ids_type     elem_node_connectivity ;
  const typename Fields::node_coords_type       model_coords ;

  const typename Fields::sym_tensor_array_view  stretch ;
  const typename Fields::tensor_array_view      rotation ;
  const typename Fields::tensor_array_view      rotation_new ;
  const typename Fields::property_view          elem_mass ;

  const Scalar   density ;
  const unsigned uq_count ;

  InitializeElement( const Fields & mesh_fields )
    : elem_node_connectivity( mesh_fields.elem_node_connectivity )
    , model_coords(           mesh_fields.model_coords )
    , stretch(                mesh_fields.stretch )
    , rotation(               mesh_fields.rotation )
    , rotation_new(           mesh_fields.rotation_new )
    , elem_mass(              mesh_fields.elem_mass )
    , density(                mesh_fields.density )
    , uq_count(               mesh_fields.stretch.dimension(1) )
    {
      KokkosArray::parallel_for( mesh_fields.num_elements , *this );
    }

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( unsigned ielem )const
  {
    const unsigned K_XX = 0 ;
    const unsigned K_YY = 1 ;
    const unsigned K_ZZ = 2 ;
    const Scalar ONE12TH = 1.0 / 12.0 ;

    Scalar x[ Fields::ElemNodeCount ];
    Scalar y[ Fields::ElemNodeCount ];
    Scalar z[ Fields::ElemNodeCount ];
    Scalar grad_x[ Fields::ElemNodeCount ];
    Scalar grad_y[ Fields::ElemNodeCount ];
    Scalar grad_z[ Fields::ElemNodeCount ];

    for ( unsigned i = 0 ; i < Fields::ElemNodeCount ; ++i ) {
      const unsigned n = elem_node_connectivity( ielem , i );

      x[i]  = model_coords( n , 0 );
      y[i]  = model_coords( n , 1 );
      z[i]  = model_coords( n , 2 );
    }

    ElemFunc::grad( x, y, z, grad_x, grad_y, grad_z);

    elem_mass(ielem) = ONE12TH * density * ElemFunc::dot8( x , grad_x );

    for ( unsigned jp = 0 ; jp < uq_count ; ++jp ) {

      stretch(ielem,jp,K_XX) = 1 ;
      stretch(ielem,jp,K_YY) = 1 ;
      stretch(ielem,jp,K_ZZ) = 1 ;

      rotation(ielem,jp,K_XX) = 1 ;
      rotation(ielem,jp,K_YY) = 1 ;
      rotation(ielem,jp,K_ZZ) = 1 ;

      rotation_new(ielem,jp,K_XX) = 1 ;
      rotation_new(ielem,jp,K_YY) = 1 ;
      rotation_new(ielem,jp,K_ZZ) = 1 ;
    }
  }
};


template<typename Scalar>
struct InitializeNode< Explicit::Fields< Scalar, KOKKOS_MACRO_DEVICE > >
{
  typedef KOKKOS_MACRO_DEVICE     device_type ;

  typedef Explicit::Fields< Scalar , device_type > Fields ;

  typename Fields::node_elem_ids_type      node_elem_connectivity ;
  typename Fields::property_view           nodal_mass ;
  typename Fields::property_view           elem_mass ;

  static const unsigned ElemNodeCount = Fields::ElemNodeCount ;

  InitializeNode( const Fields & mesh_fields )
    : node_elem_connectivity( mesh_fields.node_elem_connectivity )
    , nodal_mass(             mesh_fields.nodal_mass )
    , elem_mass(              mesh_fields.elem_mass )
    {
      KokkosArray::parallel_for( mesh_fields.num_nodes_owned , *this );
    }

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( unsigned inode )const
  {
    const unsigned begin = node_elem_connectivity.row_map[inode];
    const unsigned end   = node_elem_connectivity.row_map[inode+1];

    Scalar node_mass = 0;

    for( unsigned i = begin; i != end; ++i) {
      const unsigned elem_id = node_elem_connectivity.entries( i , 0 );
      node_mass += elem_mass(elem_id);
    }

    nodal_mass(inode) = node_mass / ElemNodeCount ;
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename Scalar >
struct GradFunctor< Explicit::Fields< Scalar , KOKKOS_MACRO_DEVICE > >
{
  typedef KOKKOS_MACRO_DEVICE device_type ;

  typedef Explicit::Fields< Scalar , device_type >  Fields ;
  typedef Hex8Functions<device_type> ElemFunc ;

  static const unsigned ElemNodeCount = Fields::ElemNodeCount ;
  static const unsigned K_F_SIZE      = Fields::K_F_SIZE ;
  static const unsigned K_S_SIZE      = Fields::K_S_SIZE ;

  // Global arrays used by this functor.

  const typename Fields::elem_node_ids_type  elem_node_connectivity ;
  const typename Fields::node_coords_type    model_coords ;

  const typename Fields::value_view          dt ;
  const typename Fields::geom_array_view     displacement ; 
  const typename Fields::geom_array_view     velocity ; 
  const typename Fields::tensor_array_view   vel_grad ;
  const unsigned                             uq_count ;

  // Constructor on the Host to populate this device functor.
  // All array view copies are shallow.
  GradFunctor( const Fields & fields )
    : elem_node_connectivity( fields.elem_node_connectivity )
    , model_coords(  fields.model_coords )
    , dt(            fields.dt )
    , displacement(  fields.displacement )
    , velocity(      fields.velocity )
    , vel_grad(      fields.vel_grad )
    , uq_count(      fields.displacement.dimension(1) )
    {
      KokkosArray::parallel_for( fields.num_elements , *this );
    }

  //--------------------------------------------------------------------------
  // Functor operator() which calls the three member functions.

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( unsigned ielem )const
  {
    const unsigned K_F_XX = Fields::K_F_XX ;
    const unsigned K_F_YY = Fields::K_F_YY ;
    const unsigned K_F_ZZ = Fields::K_F_ZZ ;
    const unsigned K_F_XY = Fields::K_F_XY ;
    const unsigned K_F_YZ = Fields::K_F_YZ ;
    const unsigned K_F_ZX = Fields::K_F_ZX ;
    const unsigned K_F_YX = Fields::K_F_YX ;
    const unsigned K_F_ZY = Fields::K_F_ZY ;
    const unsigned K_F_XZ = Fields::K_F_XZ ;
    const unsigned X = 0 ;
    const unsigned Y = 1 ;
    const unsigned Z = 2 ;

    const Scalar dt_scale = -0.5 * *dt;

    //  declare and reuse local data for frequently accessed data to
    //  reduce global memory reads and writes.

    unsigned elem_node[ ElemNodeCount ];

    Scalar  model_x[ ElemNodeCount ];
    Scalar  model_y[ ElemNodeCount ];
    Scalar  model_z[ ElemNodeCount ];

    Scalar  x[ ElemNodeCount ] ;
    Scalar  y[ ElemNodeCount ] ;
    Scalar  z[ ElemNodeCount ] ;

    Scalar  vx[ ElemNodeCount ] ;
    Scalar  vy[ ElemNodeCount ] ;
    Scalar  vz[ ElemNodeCount ];

    Scalar  grad_x[ ElemNodeCount ] ;
    Scalar  grad_y[ ElemNodeCount ] ;
    Scalar  grad_z[ ElemNodeCount ];

    // Read global velocity once and use many times
    // via local registers / L1 cache.
    //  store the velocity information in local memory before using,
    //  so it can be returned for other functions to use

    // Read global coordinates and velocity once and use many times
    // via local registers / L1 cache.
    // load X coordinate information and move by half time step

    for ( unsigned i = 0 ; i < ElemNodeCount ; ++i ) {
      const unsigned n = elem_node_connectivity( ielem , i );
      elem_node[i] = n ;
      model_x[i] = model_coords( n , X );
      model_y[i] = model_coords( n , Y );
      model_z[i] = model_coords( n , Z );
    }

    for ( unsigned jp = 0 ; jp < uq_count ; ++jp ) {

      for ( unsigned i = 0 ; i < ElemNodeCount ; ++i ) {

        const unsigned n = elem_node[i] ;

        vx[i] = velocity( n , jp , X );
        vy[i] = velocity( n , jp , Y );
        vz[i] = velocity( n , jp , Z );

        x[i]  = model_x[i] + displacement( n , jp , X ) + dt_scale * vx[i];
        y[i]  = model_y[i] + displacement( n , jp , Y ) + dt_scale * vy[i];
        z[i]  = model_z[i] + displacement( n , jp , Z ) + dt_scale * vz[i];
      }

      ElemFunc::grad( x, y, z, grad_x, grad_y, grad_z);

      //  Calculate hexahedral volume from x model_coords and gradient information

      const Scalar inv_vol = 1.0 / ElemFunc::dot8( x , grad_x );

      vel_grad( ielem, jp, K_F_XX ) = inv_vol * ElemFunc::dot8( vx , grad_x );
      vel_grad( ielem, jp, K_F_YX ) = inv_vol * ElemFunc::dot8( vy , grad_x );
      vel_grad( ielem, jp, K_F_ZX ) = inv_vol * ElemFunc::dot8( vz , grad_x );

      vel_grad( ielem, jp, K_F_XY ) = inv_vol * ElemFunc::dot8( vx , grad_y );
      vel_grad( ielem, jp, K_F_YY ) = inv_vol * ElemFunc::dot8( vy , grad_y );
      vel_grad( ielem, jp, K_F_ZY ) = inv_vol * ElemFunc::dot8( vz , grad_y );

      vel_grad( ielem, jp, K_F_XZ ) = inv_vol * ElemFunc::dot8( vx , grad_z );
      vel_grad( ielem, jp, K_F_YZ ) = inv_vol * ElemFunc::dot8( vy , grad_z );
      vel_grad( ielem, jp, K_F_ZZ ) = inv_vol * ElemFunc::dot8( vz , grad_z );
    }
  }
};

//----------------------------------------------------------------------------

template< typename Scalar >
struct DecompRotateFunctor< Explicit::Fields< Scalar , KOKKOS_MACRO_DEVICE > >
{
  typedef KOKKOS_MACRO_DEVICE device_type ;

  typedef Explicit::Fields< Scalar , device_type >  Fields ;

  typedef Hex8Functions< device_type > ElemFunc ;

  static const unsigned K_F_SIZE = Fields::K_F_SIZE ;
  static const unsigned K_S_SIZE = Fields::K_S_SIZE ;

  // Global arrays used by this functor.

  const typename Fields::tensor_array_view      rotation ;
  const typename Fields::tensor_array_view      rotation_new ;
  const typename Fields::tensor_array_view      vel_grad ;
  const typename Fields::sym_tensor_array_view  stretch ;
  const typename Fields::sym_tensor_array_view  rot_stretch ;
  const typename Fields::value_view             dt ;
  const unsigned uq_count ;

  DecompRotateFunctor( const Fields & mesh_fields )
    : rotation(      mesh_fields.rotation )
    , rotation_new(  mesh_fields.rotation_new )
    , vel_grad(      mesh_fields.vel_grad )
    , stretch(       mesh_fields.stretch )
    , rot_stretch(   mesh_fields.rot_stretch )
    , dt(            mesh_fields.dt )
    , uq_count(      mesh_fields.vel_grad.dimension(1) )
    {
      KokkosArray::parallel_for( mesh_fields.num_elements , *this );
    }

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( unsigned ielem ) const
  {
    const Scalar step = *dt ;

    // Local scratch space to avoid multiple accesses to global memory.

    Scalar v_gr[    K_F_SIZE ];  // Velocity gradient
    Scalar rot[     K_F_SIZE ];  // Rotation
    Scalar str[     K_S_SIZE ];  // Stretch
    Scalar str_ten[ K_S_SIZE ];  // Stretching tensor
    Scalar rot_str[ K_S_SIZE ];  // Rotated stretch

    for ( unsigned jp = 0 ; jp < uq_count ; ++jp ) {

      for ( unsigned i = 0; i < K_F_SIZE ; ++i ) v_gr[i] = vel_grad( ielem, jp, i );
      for ( unsigned i = 0; i < K_F_SIZE ; ++i ) rot[i]  = rotation( ielem, jp, i );
      for ( unsigned i = 0; i < K_S_SIZE ; ++i ) str[i]  = stretch(  ielem, jp, i );

      // Update stress, compute stretch tensor and rotated stretch
      ElemFunc::polar_decomp( step , v_gr, str, str_ten, rot );

      for ( unsigned i = 0; i < K_S_SIZE ; ++i ) stretch(      ielem, jp, i ) = str[i] ;
      for ( unsigned i = 0; i < K_F_SIZE ; ++i ) rotation_new( ielem, jp, i ) = rot[i] ;

      ElemFunc::rotate_tensor( str_ten , rot , rot_str );

      for ( unsigned i = 0 ; i < K_S_SIZE ; ++i ) rot_stretch( ielem , jp , i ) = rot_str[i] ;
    }
  }
};

//----------------------------------------------------------------------------

template< typename Scalar >
struct InternalForceFunctor< Explicit::Fields< Scalar , KOKKOS_MACRO_DEVICE > >
{
  typedef KOKKOS_MACRO_DEVICE device_type ;

  typedef Explicit::Fields< Scalar , device_type >  Fields ;
  typedef Hex8Functions< device_type > ElemFunc ;

  static const unsigned ElemNodeCount = Fields::ElemNodeCount ;
  static const unsigned SpatialDim    = Fields::SpatialDim ;
  static const unsigned K_F_SIZE = Fields::K_F_SIZE ;
  static const unsigned K_S_SIZE = Fields::K_S_SIZE ;

  //--------------------------------------------------------------------------
  // Desired next time step reduction:

  typedef Scalar value_type;

  KOKKOS_MACRO_DEVICE_FUNCTION
  static void init(value_type & next_time_step ) {
     next_time_step  = 1.0e32;
  }

  KOKKOS_MACRO_DEVICE_FUNCTION
  static void join( volatile value_type & next_time_step ,
                    const volatile value_type & source )
  {
    next_time_step = next_time_step < source ? next_time_step : source;
  }

  struct SetNextTimeStep {
    typedef KOKKOS_MACRO_DEVICE  device_type ;
    typedef Scalar               value_type;

    const typename Fields::value_view dt ;

    SetNextTimeStep( const typename Fields::value_view & arg_dt )
      : dt( arg_dt ) {}

    KOKKOS_MACRO_DEVICE_FUNCTION
    void operator()( const value_type & result ) const
    {
      *dt = result ;
    }
  };

  //--------------------------------------------------------------------------

  // Global arrays used by this functor.

  const typename Fields::elem_node_ids_type      elem_node_connectivity ;
  const typename Fields::node_coords_type        model_coords ;

  const typename Fields::value_view             dt ;
  const typename Fields::geom_array_view        displacement ;
  const typename Fields::geom_array_view        velocity ;
  const typename Fields::property_view          elem_mass ;
  const typename Fields::array_view             internal_energy ;
  const typename Fields::sym_tensor_array_view  stress ;
  const typename Fields::elem_node_geom_view    element_force ;
  const typename Fields::tensor_array_view      rotation_new ;
  const typename Fields::sym_tensor_array_view  rot_stretch ;

  const Scalar     two_mu ;
  const Scalar     bulk_modulus ;
  const Scalar     lin_bulk_visc ;
  const Scalar     quad_bulk_visc ;
  const Scalar     user_dt ;
  const unsigned   uq_count ;

  InternalForceFunctor( const Fields & mesh_fields ,
                        const Scalar arg_user_dt )
    : elem_node_connectivity( mesh_fields.elem_node_connectivity )
    , model_coords(           mesh_fields.model_coords )
    , dt(                     mesh_fields.dt )
    , displacement(           mesh_fields.displacement )
    , velocity(               mesh_fields.velocity )
    , elem_mass(              mesh_fields.elem_mass )
    , internal_energy(        mesh_fields.internal_energy )
    , stress(                 mesh_fields.stress )
    , element_force(          mesh_fields.element_force )
    , rotation_new(           mesh_fields.rotation_new )
    , rot_stretch(            mesh_fields.rot_stretch )
    , two_mu(                 mesh_fields.two_mu )
    , bulk_modulus(           mesh_fields.bulk_modulus )
    , lin_bulk_visc(          mesh_fields.lin_bulk_visc )
    , quad_bulk_visc(         mesh_fields.quad_bulk_visc )
    , user_dt(                arg_user_dt )
    , uq_count(               mesh_fields.displacement.dimension(1) )
  {
    SetNextTimeStep op_dt( mesh_fields.dt_new );

    KokkosArray::parallel_reduce( mesh_fields.num_elements, *this, op_dt );
  }

  //--------------------------------------------------------------------------

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( unsigned ielem, value_type & next_time_step )const
  {
    const Scalar ONE12TH = 1.0 / 12.0 ;

    unsigned node_id[ ElemNodeCount ];
    Scalar mx[ ElemNodeCount ] ;
    Scalar my[ ElemNodeCount ] ;
    Scalar mz[ ElemNodeCount ] ;

    Scalar x[ ElemNodeCount ] ;
    Scalar y[ ElemNodeCount ] ;
    Scalar z[ ElemNodeCount ] ;
    Scalar vx[ ElemNodeCount ] ;
    Scalar vy[ ElemNodeCount ] ;
    Scalar vz[ ElemNodeCount ];
    Scalar grad_x[ ElemNodeCount ] ;
    Scalar grad_y[ ElemNodeCount ] ;
    Scalar grad_z[ ElemNodeCount ] ;
    Scalar force[ ElemNodeCount ][ SpatialDim ] ;

    Scalar rot[ K_F_SIZE ]; // New rotation
    Scalar stress_work[ K_S_SIZE ]; // New stress
    Scalar rot_str[ K_S_SIZE ]; // rotated stretch
    Scalar total_stress12th[ K_S_SIZE ];

    const Scalar mass = elem_mass( ielem );
    const Scalar step = *dt ;

    for ( unsigned i = 0 ; i < ElemNodeCount ; ++i ) {
      const unsigned n = elem_node_connectivity(ielem,i);

      node_id[i] = n ;
      mx[i] = model_coords( n, 0 );
      my[i] = model_coords( n, 1 );
      mz[i] = model_coords( n, 2 );
    }

    for ( unsigned kp = 0 ; kp < uq_count ; ++kp ) {

      // Position and velocity:

      for ( unsigned i = 0 ; i < ElemNodeCount ; ++i ) {
        const unsigned n = node_id[i] ;

        x[i] = mx[i] + displacement(n, kp, 0 );
        y[i] = my[i] + displacement(n, kp, 1 );
        z[i] = mz[i] + displacement(n, kp, 2 );

        vx[i] = velocity(n, kp, 0 );
        vy[i] = velocity(n, kp, 1 );
        vz[i] = velocity(n, kp, 2 );
      }

      // Gradient:

      ElemFunc::grad( x , y , z , grad_x , grad_y , grad_z );

      const Scalar mid_vol = ElemFunc::dot8( x , grad_x );

      const Scalar aspect = 6.0 * mid_vol /
                            ( ElemFunc::dot8( grad_x , grad_x ) +
                              ElemFunc::dot8( grad_y , grad_y ) +
                              ElemFunc::dot8( grad_z , grad_z ) );

      const Scalar shr    = two_mu ;
      const Scalar dil    = bulk_modulus + ((2.0*shr)/3.0);
      const Scalar dtrial = sqrt( mass * aspect / dil );

      for ( unsigned i = 0 ; i < K_S_SIZE ; ++i ) {
        rot_str[i] = rot_stretch( ielem , kp , i );
      }

      const Scalar traced = rot_str[ Fields::K_S_XX ] +
                            rot_str[ Fields::K_S_YY ] +
                            rot_str[ Fields::K_S_ZZ ] ;

      const Scalar eps    = traced < 0 ? ( lin_bulk_visc - quad_bulk_visc * traced * dtrial ) : lin_bulk_visc ;
      const Scalar bulkq  = eps * dil * dtrial * traced ;

      for ( unsigned i = 0 ; i < K_F_SIZE ; ++i ) { rot[i] = rotation_new( ielem , kp , i ); }
      for ( unsigned i = 0 ; i < K_S_SIZE ; ++i ) { stress_work[i] = stress( ielem , kp , i ); }

      ElemFunc::update_stress( step , two_mu , bulk_modulus , rot_str , stress_work );

      for ( unsigned i = 0 ; i < K_S_SIZE ; ++i ) { stress( ielem , kp , i ) = stress_work[i]; }

      ElemFunc::rotate_tensor_backward( stress_work , rot , total_stress12th );

      total_stress12th[0] = ONE12TH*( total_stress12th[ 0 ] + bulkq );
      total_stress12th[1] = ONE12TH*( total_stress12th[ 1 ] + bulkq );
      total_stress12th[2] = ONE12TH*( total_stress12th[ 2 ] + bulkq );
      total_stress12th[3] = ONE12TH*( total_stress12th[ 3 ] );
      total_stress12th[4] = ONE12TH*( total_stress12th[ 4 ] );
      total_stress12th[5] = ONE12TH*( total_stress12th[ 5 ] );

      Scalar energy ;

      ElemFunc::comp_force( vx, vy, vz, grad_x, grad_y, grad_z, total_stress12th, force, energy );

      for ( unsigned i = 0 ; i < ElemNodeCount ; ++i ) {
        element_force( ielem, kp , 0 , i ) = force[i][0] ;
        element_force( ielem, kp , 1 , i ) = force[i][1] ;
        element_force( ielem, kp , 2 , i ) = force[i][2] ;
      }

      internal_energy( ielem , kp ) = energy ;

      // Next time step

      const Scalar desired_time_step =
        user_dt > 0 ? user_dt 
                    : dtrial * ( sqrt( 1.0 + eps * eps ) - eps );

      next_time_step = next_time_step < desired_time_step
                     ? next_time_step : desired_time_step ;
    }
  }
};

//----------------------------------------------------------------------------

template< typename Scalar >
struct NodalBoundary< Explicit::Fields< Scalar , KOKKOS_MACRO_DEVICE > >
{
  typedef KOKKOS_MACRO_DEVICE device_type ;
  typedef Explicit::Fields< Scalar , device_type >  Fields ;

  const typename Fields::node_coords_type  model_coords ;
  const Scalar x_bc ;

  NodalBoundary( const NodalBoundary & rhs )
    : model_coords( rhs.model_coords )
    , x_bc( rhs.x_bc )
    {}

  NodalBoundary( const Fields & mesh_fields ,
                 const Scalar arg_x_bc )
    : model_coords( mesh_fields.model_coords )
    , x_bc( arg_x_bc )
    {}

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const unsigned inode , Scalar accel[] ) const
  {
    const Scalar tol_bc = 1.0e-7;
    const bool fixed_bc = fabs( model_coords(inode,0) - x_bc ) < tol_bc ;

    if ( fixed_bc ) {
      accel[0] = 0 ;
      accel[1] = 0 ;
      accel[2] = 0 ;
    }
  }
};

template< typename Scalar , class Boundary >
struct NodalUpdateFunctor< Fields< Scalar , KOKKOS_MACRO_DEVICE > , Boundary >
{
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef device_type::size_type  size_type;

  typedef Explicit::Fields< Scalar , device_type >  Fields ;

  const Boundary                              boundary ;
  const typename Fields::node_elem_ids_type   node_elem_connectivity ;
  const typename Fields::node_coords_type     model_coords ;
  const typename Fields::property_view        nodal_mass ;

  const typename Fields::value_view           dt ;
  const typename Fields::value_view           dt_new ;
  const typename Fields::geom_array_view      displacement ;
  const typename Fields::geom_array_view      displacement_new ;
  const typename Fields::geom_array_view      velocity ;
  const typename Fields::geom_array_view      velocity_new ;
  const typename Fields::geom_array_view      acceleration ;
  const typename Fields::geom_array_view      internal_force ;
  const typename Fields::elem_node_geom_view  element_force ;
  const unsigned                              uq_count ;

  NodalUpdateFunctor( const Fields  & mesh_fields ,
                      const Boundary & arg_boundary )
   : boundary( arg_boundary )
   , node_elem_connectivity( mesh_fields.node_elem_connectivity )
   , model_coords(      mesh_fields.model_coords )
   , nodal_mass(        mesh_fields.nodal_mass )
   , dt(                mesh_fields.dt )
   , dt_new(            mesh_fields.dt_new )
   , displacement(      mesh_fields.displacement )
   , displacement_new(  mesh_fields.displacement_new )
   , velocity(          mesh_fields.velocity )
   , velocity_new(      mesh_fields.velocity_new )
   , acceleration(      mesh_fields.acceleration )
   , internal_force(    mesh_fields.internal_force )
   , element_force(     mesh_fields.element_force )
   , uq_count(          mesh_fields.displacement.dimension(1) )
   {
        //std::cout << "finish_step dt: " << dt << std::endl;
        //std::cout << "finish_step prev_dt: " << prev_dt << std::endl;

      // Only update the owned nodes:

      KokkosArray::parallel_for( mesh_fields.num_nodes_owned , *this );
   }

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( unsigned inode ) const
  {
    // Getting count as per 'CSR-like' data structure
    const unsigned begin = node_elem_connectivity.row_map[inode];
    const unsigned end   = node_elem_connectivity.row_map[inode+1];

    const Scalar m       = nodal_mass( inode );
    const Scalar dt_disp = *dt_new ;
    const Scalar dt_vel  = ( dt_disp + *dt ) / 2.0 ;

    for ( unsigned kp = 0 ; kp < uq_count ; ++kp ) {

      double local_force[] = { 0.0 , 0.0 , 0.0 };

      // Gather-sum internal force from
      // each element that a node is attached to.

      for ( unsigned i = begin; i < end ; ++i ){

        const unsigned ielem           = node_elem_connectivity.entries( i, 0);
        const unsigned elem_node_index = node_elem_connectivity.entries( i, 1);

        local_force[0] += element_force(ielem, kp, 0, elem_node_index);
        local_force[1] += element_force(ielem, kp, 1, elem_node_index);
        local_force[2] += element_force(ielem, kp, 2, elem_node_index);
      }

      internal_force(inode, kp, 0) = local_force[0];
      internal_force(inode, kp, 1) = local_force[1];
      internal_force(inode, kp, 2) = local_force[2];

      // Acceleration:

      Scalar a_current[3];

      a_current[0] = - local_force[0] / m ;
      a_current[1] = - local_force[1] / m ;
      a_current[2] = - local_force[2] / m ;

      // Boundary condition update to acceleration:

      boundary( inode , a_current );

      acceleration(inode,kp,0) = a_current[0] ;
      acceleration(inode,kp,1) = a_current[1] ;
      acceleration(inode,kp,2) = a_current[2] ;

      // Central difference time integration:

      Scalar v_new[3];

      v_new[0] = velocity(inode,kp,0) + dt_vel * a_current[0];
      v_new[1] = velocity(inode,kp,1) + dt_vel * a_current[1];
      v_new[2] = velocity(inode,kp,2) + dt_vel * a_current[2];

      velocity_new(inode,kp,0) = v_new[0] ;
      velocity_new(inode,kp,1) = v_new[1] ;
      velocity_new(inode,kp,2) = v_new[2] ;

      displacement_new(inode,kp,0) = displacement(inode,kp,0) + dt_disp * v_new[0];
      displacement_new(inode,kp,1) = displacement(inode,kp,1) + dt_disp * v_new[1];
      displacement_new(inode,kp,2) = displacement(inode,kp,2) + dt_disp * v_new[2];
    }
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename Scalar >
struct PackState< Explicit::Fields< Scalar , KOKKOS_MACRO_DEVICE > >
{
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef device_type::size_type  size_type ;

  typedef Explicit::Fields< Scalar , device_type >  Fields ;

  typedef typename Fields::geom_array_view::value_type    value_type ;
  typedef KokkosArray::View< value_type[], device_type >  buffer_type ;

  static const unsigned value_count = 6 ;

  const typename Fields::geom_array_view  displacement ;
  const typename Fields::geom_array_view  velocity ;
  const buffer_type  output ;
  const unsigned     inode_base ;
  const unsigned     uq_count ;

  PackState( const buffer_type & arg_output ,
             const Fields      & mesh_fields ,
             const unsigned      arg_begin ,
             const unsigned      arg_count )
    : displacement( mesh_fields.displacement )
    , velocity(     mesh_fields.velocity )
    , output(       arg_output )
    , inode_base(   arg_begin )
    , uq_count(     mesh_fields.displacement.dimension(1) )
    {
      KokkosArray::parallel_for( arg_count , *this );
    }

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const unsigned i ) const
  {
    const unsigned inode  = inode_base + i ;
    const unsigned jbegin = i * value_count ;

    for ( unsigned kq = 0 ; kq < uq_count ; ++kq ) {

      unsigned j = jbegin ;

      output[j++] = displacement( inode , kq , 0 );
      output[j++] = displacement( inode , kq , 1 );
      output[j++] = displacement( inode , kq , 2 );
      output[j++] = velocity( inode , kq , 0 );
      output[j++] = velocity( inode , kq , 1 );
      output[j++] = velocity( inode , kq , 2 );
    }
  }
};

template< typename Scalar >
struct UnpackState< Explicit::Fields< Scalar , KOKKOS_MACRO_DEVICE > >
{
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef device_type::size_type  size_type ;

  typedef Explicit::Fields< Scalar , device_type >  Fields ;

  typedef typename Fields::geom_array_view::value_type     value_type ;
  typedef KokkosArray::View< value_type[] , device_type >  buffer_type ;

  static const unsigned value_count = 6 ;

  const typename Fields::geom_array_view  displacement ;
  const typename Fields::geom_array_view  velocity ;
  const buffer_type  input ;
  const unsigned     inode_base ;
  const unsigned     uq_count ;

  UnpackState( const buffer_type & arg_input ,
               const Fields      & mesh_fields ,
               const unsigned      arg_begin ,
               const unsigned      arg_count )
    : displacement( mesh_fields.displacement )
    , velocity(     mesh_fields.velocity )
    , input(        arg_input )
    , inode_base(   arg_begin )
    , uq_count(     mesh_fields.displacement.dimension(1) )
    {
      KokkosArray::parallel_for( arg_count , *this );
    }

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const unsigned i ) const
  {
    const unsigned inode = inode_base + i ;
    const unsigned jbegin = i * value_count ;

    for ( unsigned kq = 0 ; kq < uq_count ; ++kq ) {
      unsigned j = jbegin ;

      displacement( inode , kq , 0 ) = input[j++] ;
      displacement( inode , kq , 1 ) = input[j++] ;
      displacement( inode , kq , 2 ) = input[j++] ;
      velocity( inode , kq , 0 ) = input[j++] ;
      velocity( inode , kq , 1 ) = input[j++] ;
      velocity( inode , kq , 2 ) = input[j++] ;
    }
  }
};

} /* namespace Explicit */



