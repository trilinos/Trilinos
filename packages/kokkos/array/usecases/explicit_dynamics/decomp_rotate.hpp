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

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#ifndef DECOMP_ROTATE
#define DECOMP_ROTATE
//
//  Indexes into a 3 by 3 symmetric tensor stored as a length 6 vector
//
#define K_S_XX 0
#define K_S_YY 1
#define K_S_ZZ 2
#define K_S_XY 3
#define K_S_YZ 4
#define K_S_ZX 5
#define K_S_YX 3
#define K_S_ZY 4
#define K_S_XZ 5
//
//  Indexes into a 3 by 3 skew symmetric tensor stored as a length 9 vector
//
#define K_V_XY 0
#define K_V_YZ 1
#define K_V_ZX 2


template< typename Scalar , class DeviceType >
struct decomp_rotate;

template<typename Scalar>
struct decomp_rotate<Scalar, KOKKOSARRAY_MACRO_DEVICE>{

  typedef KOKKOSARRAY_MACRO_DEVICE     device_type ;
  typedef typename KokkosArray::MDArray<Scalar,device_type> array_type ;

  typedef KokkosArray::Value<Scalar,device_type>     scalar;

  const array_type rotation;
  const array_type vel_grad;
  const array_type stretch;
  const array_type vorticity;
  const array_type rot_stretch;
  const scalar     dt_value;
  const int        current_state;
  const int        previous_state;

  typedef Region<Scalar,device_type> MyRegion;

  decomp_rotate(  const MyRegion & region,
          const int arg_current_state,
          const int arg_previous_state)
    : rotation( region.rotation )
    , vel_grad( region.vel_grad )
    , stretch( region.stretch )
    , vorticity( region.vorticity )
    , rot_stretch( region.rot_stretch )
    , dt_value( region.dt)
    , current_state(arg_current_state)
    , previous_state(arg_previous_state)
    {
    }



  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  void additive_decomp(int ielem, Scalar * v_gr, Scalar * str_ten)const {

  //  In addition to calculating stretching_tensor,
  //  use this as an opportunity to load global
  //  variables into a local space

    v_gr[K_F_XX] = vel_grad(ielem, K_F_XX);
    v_gr[K_F_YY] = vel_grad(ielem, K_F_YY);
    v_gr[K_F_ZZ] = vel_grad(ielem, K_F_ZZ);
    v_gr[K_F_XY] = vel_grad(ielem, K_F_XY);
    v_gr[K_F_YX] = vel_grad(ielem, K_F_YX);
    v_gr[K_F_YZ] = vel_grad(ielem, K_F_YZ);
    v_gr[K_F_ZX] = vel_grad(ielem, K_F_ZX);
    v_gr[K_F_ZY] = vel_grad(ielem, K_F_ZY);
    v_gr[K_F_XZ] = vel_grad(ielem, K_F_XZ);
  //
  //  Symmetric part
  //
    str_ten[K_S_XX] = v_gr[K_F_XX];
    str_ten[K_S_YY] = v_gr[K_F_YY];
    str_ten[K_S_ZZ] = v_gr[K_F_ZZ];
    str_ten[K_S_XY] = 0.5*(v_gr[K_F_XY] + v_gr[K_F_YX]);
    str_ten[K_S_YZ] = 0.5*(v_gr[K_F_YZ] + v_gr[K_F_ZY]);
    str_ten[K_S_ZX] = 0.5*(v_gr[K_F_ZX] + v_gr[K_F_XZ]);

  }

  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  void polar_decomp(int ielem, Scalar * v_gr, Scalar * str_ten, Scalar * str, Scalar * vort, Scalar * rot_old, Scalar * rot_new)const {

    Scalar dt = *dt_value;

    //  Skew Symmetric part
    vorticity(ielem, K_V_XY) = vort[K_V_XY] = 0.5*(v_gr[K_F_XY] - v_gr[K_F_YX]);
    vorticity(ielem, K_V_YZ) = vort[K_V_YZ] = 0.5*(v_gr[K_F_YZ] - v_gr[K_F_ZY]);
    vorticity(ielem, K_V_ZX) = vort[K_V_ZX] = 0.5*(v_gr[K_F_ZX] - v_gr[K_F_XZ]);

    //   calculate the rates of rotation via gauss elimination.
    str[K_S_XX] = stretch(ielem, K_S_XX);
    str[K_S_YY] = stretch(ielem, K_S_YY);
    str[K_S_ZZ] = stretch(ielem, K_S_ZZ);
    str[K_S_XY] = stretch(ielem, K_S_XY);
    str[K_S_YZ] = stretch(ielem, K_S_YZ);
    str[K_S_ZX] = stretch(ielem, K_S_ZX);

    Scalar z1 = str_ten[K_S_XY] * str[K_S_ZX] - str_ten[K_S_ZX] * str[K_S_XY] + str_ten[K_S_YY] * str[K_S_YZ] -
        str_ten[K_S_YZ] * str[K_S_YY] + str_ten[K_S_YZ] * str[K_S_ZZ] - str_ten[K_S_ZZ] * str[K_S_YZ];
    Scalar z2 = str_ten[K_S_ZX] * str[K_S_XX] - str_ten[K_S_XX] * str[K_S_ZX] + str_ten[K_S_YZ] * str[K_S_XY] -
        str_ten[K_S_XY] * str[K_S_YZ] + str_ten[K_S_ZZ] * str[K_S_ZX] - str_ten[K_S_ZX] * str[K_S_ZZ];
    Scalar z3 = str_ten[K_S_XX] * str[K_S_XY] - str_ten[K_S_XY] * str[K_S_XX] + str_ten[K_S_XY] * str[K_S_YY] -
        str_ten[K_S_YY] * str[K_S_XY] + str_ten[K_S_ZX] * str[K_S_YZ] - str_ten[K_S_YZ] * str[K_S_ZX];

  //   forward elimination
    const Scalar a1inv = 1.0 / (str[K_S_YY] + str[K_S_ZZ]);

    const Scalar a4BYa1 = -1 * str[K_S_XY] * a1inv;

    const Scalar a2inv = 1.0 / (str[K_S_ZZ] + str[K_S_XX] + str[K_S_XY] * a4BYa1);

    const Scalar a5 =  -str[K_S_YZ] + str[K_S_ZX] * a4BYa1;

    z2 -= z1 * a4BYa1;
    Scalar a6BYa1 = -1 * str[K_S_ZX] * a1inv;
    const Scalar a5BYa2 = a5 * a2inv;
    z3 -= z1 * a6BYa1 - z2 * a5BYa2;

  //   backward substitution -
    z3 /= (str[K_S_XX] + str[K_S_YY] + str[K_S_ZX] * a6BYa1 + a5 * a5BYa2);
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
    const Scalar dt_half = 0.5 * dt;

    rot_old[K_F_XX] = rotation(ielem, K_F_XX, previous_state);
    rot_old[K_F_YX] = rotation(ielem, K_F_YX, previous_state);
    rot_old[K_F_ZX] = rotation(ielem, K_F_ZX, previous_state);
    rot_old[K_F_XY] = rotation(ielem, K_F_XY, previous_state);
    rot_old[K_F_YY] = rotation(ielem, K_F_YY, previous_state);
    rot_old[K_F_ZY] = rotation(ielem, K_F_ZY, previous_state);
    rot_old[K_F_XZ] = rotation(ielem, K_F_XZ, previous_state);
    rot_old[K_F_YZ] = rotation(ielem, K_F_YZ, previous_state);
    rot_old[K_F_ZZ] = rotation(ielem, K_F_ZZ, previous_state);

    Scalar r_XX = rot_old[K_F_XX] + dt_half*( z3 * rot_old[K_F_YX] - z2 * rot_old[K_F_ZX] );
    Scalar r_YX = rot_old[K_F_YX] + dt_half*( z1 * rot_old[K_F_ZX] - z3 * rot_old[K_F_XX] );
    Scalar r_ZX = rot_old[K_F_ZX] + dt_half*( z2 * rot_old[K_F_XX] - z1 * rot_old[K_F_YX] );
    Scalar r_XY = rot_old[K_F_XY] + dt_half*( z3 * rot_old[K_F_YY] - z2 * rot_old[K_F_ZY] );
    Scalar r_YY = rot_old[K_F_YY] + dt_half*( z1 * rot_old[K_F_ZY] - z3 * rot_old[K_F_XY] );
    Scalar r_ZY = rot_old[K_F_ZY] + dt_half*( z2 * rot_old[K_F_XY] - z1 * rot_old[K_F_YY] );
    Scalar r_XZ = rot_old[K_F_XZ] + dt_half*( z3 * rot_old[K_F_YZ] - z2 * rot_old[K_F_ZZ] );
    Scalar r_YZ = rot_old[K_F_YZ] + dt_half*( z1 * rot_old[K_F_ZZ] - z3 * rot_old[K_F_XZ] );
    Scalar r_ZZ = rot_old[K_F_ZZ] + dt_half*( z2 * rot_old[K_F_XZ] - z1 * rot_old[K_F_YZ] );


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

    rotation(ielem, K_F_ZX, current_state)  = rot_new[K_F_ZX] = r_ZX * a33inv;
    rotation(ielem, K_F_ZY, current_state)  = rot_new[K_F_ZY] = r_ZY * a33inv;
    rotation(ielem, K_F_ZZ, current_state)  = rot_new[K_F_ZZ] = r_ZZ * a33inv;
    rotation(ielem, K_F_YX, current_state)  = rot_new[K_F_YX] = ( r_YX - rot_new[K_F_ZX] * a23 ) * a22inv;
    rotation(ielem, K_F_YY, current_state)  = rot_new[K_F_YY] = ( r_YY - rot_new[K_F_ZY] * a23 ) * a22inv;
    rotation(ielem, K_F_YZ, current_state)  = rot_new[K_F_YZ] = ( r_YZ - rot_new[K_F_ZZ] * a23 ) * a22inv;
    rotation(ielem, K_F_XX, current_state)  = rot_new[K_F_XX] = r_XX - rot_new[K_F_ZX] * a13 - rot_new[K_F_YX] * a12;
    rotation(ielem, K_F_XY, current_state)  = rot_new[K_F_XY] = r_XY - rot_new[K_F_ZY] * a13 - rot_new[K_F_YY] * a12;
    rotation(ielem, K_F_XZ, current_state)  = rot_new[K_F_XZ] = r_XZ - rot_new[K_F_ZZ] * a13 - rot_new[K_F_YZ] * a12;

  //   update stretch tensor in the new configuration -
    const Scalar a1 = str_ten[K_S_XY] + vort[K_V_XY];
    const Scalar a2 = str_ten[K_S_YZ] + vort[K_V_YZ];
    const Scalar a3 = str_ten[K_S_ZX] + vort[K_V_ZX];
    const Scalar b1 = str_ten[K_S_ZX] - vort[K_V_ZX];
    const Scalar b2 = str_ten[K_S_XY] - vort[K_V_XY];
    const Scalar b3 = str_ten[K_S_YZ] - vort[K_V_YZ];

    const Scalar s_XX = str[K_S_XX];
    const Scalar s_YY = str[K_S_YY];
    const Scalar s_ZZ = str[K_S_ZZ];
    const Scalar s_XY = str[K_S_XY];
    const Scalar s_YZ = str[K_S_YZ];
    const Scalar s_ZX = str[K_S_ZX];

    str[K_S_XX] += dt * (str_ten[K_S_XX] * s_XX + ( a1 + z3 ) * s_XY + ( b1 - z2 ) * s_ZX);
    str[K_S_YY] += dt * (str_ten[K_S_YY] * s_YY + ( a2 + z1 ) * s_YZ + ( b2 - z3 ) * s_XY);
    str[K_S_ZZ] += dt * (str_ten[K_S_ZZ] * s_ZZ + ( a3 + z2 ) * s_ZX + ( b3 - z1 ) * s_YZ);
    str[K_S_XY] += dt * (str_ten[K_S_XX] * s_XY + ( a1 )      * s_YY + ( b1      ) * s_YZ - z3 * s_XX + z1 * s_ZX);
    str[K_S_YZ] += dt * (str_ten[K_S_YY] * s_YZ + ( a2 )      * s_ZZ + ( b2      ) * s_ZX - z1 * s_YY + z2 * s_XY);
    str[K_S_ZX] += dt * (str_ten[K_S_ZZ] * s_ZX + ( a3 )      * s_XX + ( b3      ) * s_XY - z2 * s_ZZ + z3 * s_YZ);

  }


  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  void rotate_tensor(int ielem, Scalar * str_ten, Scalar * str, Scalar * rot_new)const {

    Scalar t[9];

    t[0] = str_ten[K_S_XX]*rot_new[K_F_XX] + str_ten[K_S_XY]*rot_new[K_F_YX] + str_ten[K_S_XZ]*rot_new[K_F_ZX];
    t[1] = str_ten[K_S_YX]*rot_new[K_F_XX] + str_ten[K_S_YY]*rot_new[K_F_YX] + str_ten[K_S_YZ]*rot_new[K_F_ZX];
    t[2] = str_ten[K_S_ZX]*rot_new[K_F_XX] + str_ten[K_S_ZY]*rot_new[K_F_YX] + str_ten[K_S_ZZ]*rot_new[K_F_ZX];
    t[3] = str_ten[K_S_XX]*rot_new[K_F_XY] + str_ten[K_S_XY]*rot_new[K_F_YY] + str_ten[K_S_XZ]*rot_new[K_F_ZY];
    t[4] = str_ten[K_S_YX]*rot_new[K_F_XY] + str_ten[K_S_YY]*rot_new[K_F_YY] + str_ten[K_S_YZ]*rot_new[K_F_ZY];
    t[5] = str_ten[K_S_ZX]*rot_new[K_F_XY] + str_ten[K_S_ZY]*rot_new[K_F_YY] + str_ten[K_S_ZZ]*rot_new[K_F_ZY];
    t[6] = str_ten[K_S_XX]*rot_new[K_F_XZ] + str_ten[K_S_XY]*rot_new[K_F_YZ] + str_ten[K_S_XZ]*rot_new[K_F_ZZ];
    t[7] = str_ten[K_S_YX]*rot_new[K_F_XZ] + str_ten[K_S_YY]*rot_new[K_F_YZ] + str_ten[K_S_YZ]*rot_new[K_F_ZZ];
    t[8] = str_ten[K_S_ZX]*rot_new[K_F_XZ] + str_ten[K_S_ZY]*rot_new[K_F_YZ] + str_ten[K_S_ZZ]*rot_new[K_F_ZZ];

    rot_stretch(ielem, K_S_XX) = rot_new[K_F_XX] * t[0] + rot_new[K_F_YX] * t[1] + rot_new[K_F_ZX] * t[2];
    rot_stretch(ielem, K_S_YY) = rot_new[K_F_XY] * t[3] + rot_new[K_F_YY] * t[4] + rot_new[K_F_ZY] * t[5];
    rot_stretch(ielem, K_S_ZZ) = rot_new[K_F_XZ] * t[6] + rot_new[K_F_YZ] * t[7] + rot_new[K_F_ZZ] * t[8];

    rot_stretch(ielem, K_S_XY) = rot_new[K_F_XX] * t[3] + rot_new[K_F_YX] * t[4] + rot_new[K_F_ZX] * t[5];
    rot_stretch(ielem, K_S_YZ) = rot_new[K_F_XY] * t[6] + rot_new[K_F_YY] * t[7] + rot_new[K_F_ZY] * t[8];
    rot_stretch(ielem, K_S_ZX) = rot_new[K_F_XZ] * t[0] + rot_new[K_F_YZ] * t[1] + rot_new[K_F_ZZ] * t[2];

    stretch(ielem, 0) = str[0];
    stretch(ielem, 1) = str[1];
    stretch(ielem, 2) = str[2];
    stretch(ielem, 3) = str[3];
    stretch(ielem, 4) = str[4];
    stretch(ielem, 5) = str[5];

  }

  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  void operator()( int ielem )const {

    //   Local scratch space to avoid multiple
    //   accesses to global memory.
    Scalar str_ten[6]; // Stretching tensor
    Scalar str[6];     // Stretch
    Scalar rot_old[9]; // Rotation old
    Scalar rot_new[9]; // Rotation new
    Scalar vort[3];    // Vorticity
    Scalar v_gr[9];    // Velocity gradient

    additive_decomp(ielem, v_gr, str_ten);

    polar_decomp(ielem, v_gr, str_ten, str, vort, rot_old, rot_new);

    rotate_tensor(ielem, str_ten, str, rot_new);

  }

};

#endif
