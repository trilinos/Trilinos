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

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#ifndef DIVERGENCE
#define DIVERGENCE

//#define ONE12TH 0.083333333333333333333333
#define ONE12TH (1.0/12.0)
//
//
//
#define HG_X1 0
#define HG_Y1 1
#define HG_Z1 2
#define HG_X2 3
#define HG_Y2 4
#define HG_Z2 5
#define HG_X3 6
#define HG_Y3 7
#define HG_Z3 8
#define HG_X4 9
#define HG_Y4 10
#define HG_Z4 11
//
//  Indexes into a full 3 by 3 tensor stored as a length 9 vector
//
#define K_F_XX 0
#define K_F_YY 1
#define K_F_ZZ 2
#define K_F_XY 3
#define K_F_YZ 4
#define K_F_ZX 5
#define K_F_YX 6
#define K_F_ZY 7
#define K_F_XZ 8


template< typename Scalar , class DeviceType >
struct divergence;

template<typename Scalar>
struct divergence<Scalar, KOKKOSARRAY_MACRO_DEVICE>{

  typedef KOKKOSARRAY_MACRO_DEVICE     device_type ;
  typedef typename KokkosArray::MDArray<Scalar,device_type> array_type ;
  typedef typename KokkosArray::MDArray<int,device_type>    int_array_type ;

  typedef KokkosArray::Value<Scalar,device_type>     scalar;

  typedef Scalar value_type;

  typedef Region<Scalar,device_type> MyRegion;

  const int_array_type elem_node_connectivity;

  const array_type  model_coords;
  const array_type  displacement;

  const array_type  velocity;
  const array_type  element_force;
  const array_type  vorticity;
  const array_type  rotation;
  const array_type  stress_new;
  const array_type  rot_stress;
  const array_type  rot_stretch;
  const array_type  gradop12;
  const array_type  elem_mass;
  const array_type  elem_dilmod;
  const array_type  elem_shrmod;
  //const array_type  elem_t_step;
  const array_type  internal_energy;
  const array_type  mid_vol;

  const array_type  hgop;
  const array_type  hg_resist;
  const array_type  hg_energy;

  const Scalar  two_mu;
  const Scalar  bulk_modulus;

  const Scalar     hg_stiffness;
  const Scalar     hg_viscosity;
  const Scalar     lin_bulk_visc;
  const Scalar     quad_bulk_visc;

  const Scalar     user_dt;
  const scalar     dt;

  const int        current_state;
  const int        previous_state;

  divergence(
      const MyRegion & region,
      const Scalar arg_user_dt,
      const int arg_current_state,
      const int arg_previous_state
      )
      : elem_node_connectivity(region.elem_node_connectivity)
      , model_coords(region.model_coords)
      , displacement(region.displacement)
      , velocity(region.velocity)
      , element_force(region.element_force)
      , vorticity(region.vorticity)
      , rotation(region.rotation)
      , stress_new(region.stress_new)
      , rot_stress(region.rot_stress)
      , rot_stretch(region.rot_stretch)
      , gradop12(region.gradop12)
      , elem_mass(region.elem_mass)
      , elem_dilmod(region.dilmod)
      , elem_shrmod(region.shrmod)
      //, elem_t_step(region.elem_t_step)
      , internal_energy(region.internal_energy)
      , mid_vol(region.mid_vol)
      , hgop(region.hgop)
      , hg_resist(region.hg_resist)
      , hg_energy(region.hg_energy)
      , two_mu(region.two_mu)
      , bulk_modulus(region.bulk_modulus)
      , hg_stiffness(region.hg_stiffness)
      , hg_viscosity(region.hg_viscosity)
      , lin_bulk_visc(region.lin_bulk_visc)
      , quad_bulk_visc(region.quad_bulk_visc)
      , user_dt( arg_user_dt )
      , dt( region.dt)
      , current_state(arg_current_state)
      , previous_state(arg_previous_state)
  {
  }

  KOKKOSARRAY_INLINE_FUNCTION
  static void init(value_type &update) {
    update = 1.0e32;
  }

  KOKKOSARRAY_INLINE_FUNCTION
  static void join(volatile value_type &update, const volatile value_type & source) {
    update = update < source ? update : source;
  }

  KOKKOSARRAY_INLINE_FUNCTION
    void get_nodes( int ielem, int * nodes) const
    {
      nodes[0] = elem_node_connectivity(ielem,0);
      nodes[1] = elem_node_connectivity(ielem,1);
      nodes[2] = elem_node_connectivity(ielem,2);
      nodes[3] = elem_node_connectivity(ielem,3);
      nodes[4] = elem_node_connectivity(ielem,4);
      nodes[5] = elem_node_connectivity(ielem,5);
      nodes[6] = elem_node_connectivity(ielem,6);
      nodes[7] = elem_node_connectivity(ielem,7);
    }

  KOKKOSARRAY_INLINE_FUNCTION
  void comp_grad(int ielem, int *nodes, Scalar *x, Scalar *y, Scalar *z) const {

    const int X = 0;
    const int Y = 1;
    const int Z = 2;


    // Read global coordinates once and use many times via local registers / L1 cache.
    //  load X coordinate information and move by half time step
    x[0] = model_coords(nodes[0], X) + displacement(nodes[0], X, current_state) ;
    x[1] = model_coords(nodes[1], X) + displacement(nodes[1], X, current_state) ;
    x[2] = model_coords(nodes[2], X) + displacement(nodes[2], X, current_state) ;
    x[3] = model_coords(nodes[3], X) + displacement(nodes[3], X, current_state) ;
    x[4] = model_coords(nodes[4], X) + displacement(nodes[4], X, current_state) ;
    x[5] = model_coords(nodes[5], X) + displacement(nodes[5], X, current_state) ;
    x[6] = model_coords(nodes[6], X) + displacement(nodes[6], X, current_state) ;
    x[7] = model_coords(nodes[7], X) + displacement(nodes[7], X, current_state) ;

    //   calc X difference vectors
    Scalar R42=(x[3] - x[1]);
    Scalar R52=(x[4] - x[1]);
    Scalar R54=(x[4] - x[3]);

    Scalar R63=(x[5] - x[2]);
    Scalar R83=(x[7] - x[2]);
    Scalar R86=(x[7] - x[5]);

    Scalar R31=(x[2] - x[0]);
    Scalar R61=(x[5] - x[0]);
    Scalar R74=(x[6] - x[3]);

    Scalar R72=(x[6] - x[1]);
    Scalar R75=(x[6] - x[4]);
    Scalar R81=(x[7] - x[0]);

    Scalar t1=(R63 + R54);
    Scalar t2=(R61 + R74);
    Scalar t3=(R72 + R81);

    Scalar t4 =(R86 + R42);
    Scalar t5 =(R83 + R52);
    Scalar t6 =(R75 + R31);

    //  Load Z information
    z[0] = model_coords(nodes[0], Z) + displacement(nodes[0], Z, current_state) ;
    z[1] = model_coords(nodes[1], Z) + displacement(nodes[1], Z, current_state) ;
    z[2] = model_coords(nodes[2], Z) + displacement(nodes[2], Z, current_state) ;
    z[3] = model_coords(nodes[3], Z) + displacement(nodes[3], Z, current_state) ;
    z[4] = model_coords(nodes[4], Z) + displacement(nodes[4], Z, current_state) ;
    z[5] = model_coords(nodes[5], Z) + displacement(nodes[5], Z, current_state) ;
    z[6] = model_coords(nodes[6], Z) + displacement(nodes[6], Z, current_state) ;
    z[7] = model_coords(nodes[7], Z) + displacement(nodes[7], Z, current_state) ;


    //  Calculate Y gradient from X and Z data
    gradop12(ielem, 1, 0) = (z[1] *  t1) - (z[2] * R42) - (z[3] *  t5)  + (z[4] *  t4) + (z[5] * R52) - (z[7] * R54);
    gradop12(ielem, 1, 1) = (z[2] *  t2) + (z[3] * R31) - (z[0] *  t1)  - (z[5] *  t6) + (z[6] * R63) - (z[4] * R61);
    gradop12(ielem, 1, 2) = (z[3] *  t3) + (z[0] * R42) - (z[1] *  t2)  - (z[6] *  t4) + (z[7] * R74) - (z[5] * R72);
    gradop12(ielem, 1, 3) = (z[0] *  t5) - (z[1] * R31) - (z[2] *  t3)  + (z[7] *  t6) + (z[4] * R81) - (z[6] * R83);
    gradop12(ielem, 1, 4) = (z[5] *  t3) + (z[6] * R86) - (z[7] *  t2)  - (z[0] *  t4) - (z[3] * R81) + (z[1] * R61);
    gradop12(ielem, 1, 5) = (z[6] *  t5) - (z[4] *  t3)  - (z[7] * R75) + (z[1] *  t6) - (z[0] * R52) + (z[2] * R72);
    gradop12(ielem, 1, 6) = (z[7] *  t1) - (z[5] *  t5)  - (z[4] * R86) + (z[2] *  t4) - (z[1] * R63) + (z[3] * R83);
    gradop12(ielem, 1, 7) = (z[4] *  t2) - (z[6] *  t1)  + (z[5] * R75) - (z[3] *  t6) - (z[2] * R74) + (z[0] * R54);


    //   calc Z difference vectors
    R42=(z[3] - z[1]);
    R52=(z[4] - z[1]);
    R54=(z[4] - z[3]);

    R63=(z[5] - z[2]);
    R83=(z[7] - z[2]);
    R86=(z[7] - z[5]);

     R31=(z[2] - z[0]);
    R61=(z[5] - z[0]);
    R74=(z[6] - z[3]);

    R72=(z[6] - z[1]);
    R75=(z[6] - z[4]);
    R81=(z[7] - z[0]);

    t1=(R63 + R54);
    t2=(R61 + R74);
    t3=(R72 + R81);

    t4 =(R86 + R42);
    t5 =(R83 + R52);
    t6 =(R75 + R31);

    //  Load Y information
    y[0] = model_coords(nodes[0], Y) + displacement(nodes[0], Y, current_state) ;
    y[1] = model_coords(nodes[1], Y) + displacement(nodes[1], Y, current_state) ;
    y[2] = model_coords(nodes[2], Y) + displacement(nodes[2], Y, current_state) ;
    y[3] = model_coords(nodes[3], Y) + displacement(nodes[3], Y, current_state) ;
    y[4] = model_coords(nodes[4], Y) + displacement(nodes[4], Y, current_state) ;
    y[5] = model_coords(nodes[5], Y) + displacement(nodes[5], Y, current_state) ;
    y[6] = model_coords(nodes[6], Y) + displacement(nodes[6], Y, current_state) ;
    y[7] = model_coords(nodes[7], Y) + displacement(nodes[7], Y, current_state) ;


    //  Calculate X gradient from Y and Z data
    gradop12(ielem, 0, 0) = (y[1] *  t1) - (y[2] * R42) - (y[3] *  t5) + (y[4] *  t4) + (y[5] * R52) - (y[7] * R54);
    gradop12(ielem, 0, 1) = (y[2] *  t2) + (y[3] * R31) - (y[0] *  t1) - (y[5] *  t6) + (y[6] * R63) - (y[4] * R61);
    gradop12(ielem, 0, 2) = (y[3] *  t3) + (y[0] * R42) - (y[1] *  t2) - (y[6] *  t4) + (y[7] * R74) - (y[5] * R72);
    gradop12(ielem, 0, 3) = (y[0] *  t5) - (y[1] * R31) - (y[2] *  t3) + (y[7] *  t6) + (y[4] * R81) - (y[6] * R83);
    gradop12(ielem, 0, 4) = (y[5] *  t3) + (y[6] * R86) - (y[7] *  t2) - (y[0] *  t4) - (y[3] * R81) + (y[1] * R61);
    gradop12(ielem, 0, 5) = (y[6] *  t5) - (y[4] *  t3) - (y[7] * R75) + (y[1] *  t6) - (y[0] * R52) + (y[2] * R72);
    gradop12(ielem, 0, 6) = (y[7] *  t1) - (y[5] *  t5) - (y[4] * R86) + (y[2] *  t4) - (y[1] * R63) + (y[3] * R83);
    gradop12(ielem, 0, 7) = (y[4] *  t2) - (y[6] *  t1) + (y[5] * R75) - (y[3] *  t6) - (y[2] * R74) + (y[0] * R54);


    //   calc Y difference vectors
    R42=(y[3] - y[1]);
    R52=(y[4] - y[1]);
    R54=(y[4] - y[3]);

    R63=(y[5] - y[2]);
    R83=(y[7] - y[2]);
    R86=(y[7] - y[5]);

    R31=(y[2] - y[0]);
    R61=(y[5] - y[0]);
    R74=(y[6] - y[3]);

    R72=(y[6] - y[1]);
    R75=(y[6] - y[4]);
    R81=(y[7] - y[0]);

    t1=(R63 + R54);
    t2=(R61 + R74);
    t3=(R72 + R81);

    t4 =(R86 + R42);
    t5 =(R83 + R52);
    t6 =(R75 + R31);

    //  Calculate Z gradient from X and Y data

    gradop12(ielem, 2, 0) = (x[1] *  t1) - (x[2] * R42) - (x[3] *  t5)  + (x[4] *  t4) + (x[5] * R52) - (x[7] * R54);
    gradop12(ielem, 2, 1) = (x[2] *  t2) + (x[3] * R31) - (x[0] *  t1)  - (x[5] *  t6) + (x[6] * R63) - (x[4] * R61);
    gradop12(ielem, 2, 2) = (x[3] *  t3) + (x[0] * R42) - (x[1] *  t2)  - (x[6] *  t4) + (x[7] * R74) - (x[5] * R72);
    gradop12(ielem, 2, 3) = (x[0] *  t5) - (x[1] * R31) - (x[2] *  t3)  + (x[7] *  t6) + (x[4] * R81) - (x[6] * R83);
    gradop12(ielem, 2, 4) = (x[5] *  t3) + (x[6] * R86) - (x[7] *  t2)  - (x[0] *  t4) - (x[3] * R81) + (x[1] * R61);
    gradop12(ielem, 2, 5) = (x[6] *  t5) - (x[4] *  t3)  - (x[7] * R75) + (x[1] *  t6) - (x[0] * R52) + (x[2] * R72);
    gradop12(ielem, 2, 6) = (x[7] *  t1) - (x[5] *  t5)  - (x[4] * R86) + (x[2] *  t4) - (x[1] * R63) + (x[3] * R83);
    gradop12(ielem, 2, 7) = (x[4] *  t2) - (x[6] *  t1)  + (x[5] * R75) - (x[3] *  t6) - (x[2] * R74) + (x[0] * R54);

    mid_vol(ielem) = ONE12TH * (gradop12(ielem, 0, 0) * x[0] +
                  gradop12(ielem, 0, 1) * x[1] +
                  gradop12(ielem, 0, 2) * x[2] +
                  gradop12(ielem, 0, 3) * x[3] +
                  gradop12(ielem, 0, 4) * x[4] +
                  gradop12(ielem, 0, 5) * x[5] +
                  gradop12(ielem, 0, 6) * x[6] +
                  gradop12(ielem, 0, 7) * x[7] );

  }

  KOKKOSARRAY_INLINE_FUNCTION
    void comp_hgop(  int ielem, Scalar *x, Scalar *y, Scalar *z) const {

  //   KHP: Alternatively, we could have
  //   hx0,hx1,hx2,hx3,...,hz0,hz1,hz2,hz3
    Scalar hgconst12th[12];
    Scalar inv_vol = 1.0 / (mid_vol(ielem) * 12.0);

    Scalar q0 = x[0] - x[1];
    Scalar q1 = x[2] - x[3];
    Scalar q2 = x[4] - x[5];
    Scalar q3 = x[6] - x[7];

    hgconst12th[0] = (
        (x[0] + x[1]) - (x[2]+x[3]) -
        (x[4] + x[5]) + (x[6]+x[7]) ) * inv_vol;

    hgconst12th[1] = (  q0 - q1 - q2 + q3 ) * inv_vol;
    hgconst12th[2] = (  q0 + q1 + q2 + q3 ) * inv_vol;
    hgconst12th[3] = ( -q0 - q1 + q2 + q3 ) * inv_vol;

    q0 = y[0] - y[1];
    q1 = y[2] - y[3];
    q2 = y[4] - y[5];
    q3 = y[6] - y[7];

    hgconst12th[4] = (
        (y[0] + y[1]) - (y[2]+y[3]) -
        (y[4] + y[5]) + (y[6]+y[7]) ) * inv_vol;

    hgconst12th[5] = (  q0 - q1 - q2 + q3 ) * inv_vol;
    hgconst12th[6] = (  q0 + q1 + q2 + q3 ) * inv_vol;
    hgconst12th[7] = ( -q0 - q1 + q2 + q3 ) * inv_vol;

    q0 = z[0] - z[1];
    q1 = z[2] - z[3];
    q2 = z[4] - z[5];
    q3 = z[6] - z[7];

    hgconst12th[8]  = (
        (z[0] + z[1]) - (z[2]+z[3]) -
        (z[4] + z[5]) + (z[6]+z[7]) ) * inv_vol;

    hgconst12th[9]  = (  q0 - q1 - q2 + q3 ) * inv_vol;
    hgconst12th[10] = (  q0 + q1 + q2 + q3 ) * inv_vol;
    hgconst12th[11] = ( -q0 - q1 + q2 + q3 ) * inv_vol;


    hgop(ielem,  0, 1) =  1.0 - (hgconst12th[0] * gradop12(ielem, 0, 0) + hgconst12th[4] * gradop12(ielem, 1, 0) + hgconst12th[ 8] * gradop12(ielem, 2, 0));
    hgop(ielem,  1, 1) =  1.0 - (hgconst12th[0] * gradop12(ielem, 0, 1) + hgconst12th[4] * gradop12(ielem, 1, 1) + hgconst12th[ 8] * gradop12(ielem, 2, 1));
    hgop(ielem,  2, 1) = -1.0 - (hgconst12th[0] * gradop12(ielem, 0, 2) + hgconst12th[4] * gradop12(ielem, 1, 2) + hgconst12th[ 8] * gradop12(ielem, 2, 2));
    hgop(ielem,  3, 1) = -1.0 - (hgconst12th[0] * gradop12(ielem, 0, 3) + hgconst12th[4] * gradop12(ielem, 1, 3) + hgconst12th[ 8] * gradop12(ielem, 2, 3));
    hgop(ielem,  4, 1) = -1.0 - (hgconst12th[0] * gradop12(ielem, 0, 4) + hgconst12th[4] * gradop12(ielem, 1, 4) + hgconst12th[ 8] * gradop12(ielem, 2, 4));
    hgop(ielem,  5, 1) = -1.0 - (hgconst12th[0] * gradop12(ielem, 0, 5) + hgconst12th[4] * gradop12(ielem, 1, 5) + hgconst12th[ 8] * gradop12(ielem, 2, 5));
    hgop(ielem,  6, 1) =  1.0 - (hgconst12th[0] * gradop12(ielem, 0, 6) + hgconst12th[4] * gradop12(ielem, 1, 6) + hgconst12th[ 8] * gradop12(ielem, 2, 6));
    hgop(ielem,  7, 1) =  1.0 - (hgconst12th[0] * gradop12(ielem, 0, 7) + hgconst12th[4] * gradop12(ielem, 1, 7) + hgconst12th[ 8] * gradop12(ielem, 2, 7));
    hgop(ielem,  8, 1) =  1.0 - (hgconst12th[1] * gradop12(ielem, 0, 0) + hgconst12th[5] * gradop12(ielem, 1, 0) + hgconst12th[ 9] * gradop12(ielem, 2, 0));
    hgop(ielem,  9, 1) = -1.0 - (hgconst12th[1] * gradop12(ielem, 0, 1) + hgconst12th[5] * gradop12(ielem, 1, 1) + hgconst12th[ 9] * gradop12(ielem, 2, 1));
    hgop(ielem, 10, 1) = -1.0 - (hgconst12th[1] * gradop12(ielem, 0, 2) + hgconst12th[5] * gradop12(ielem, 1, 2) + hgconst12th[ 9] * gradop12(ielem, 2, 2));
    hgop(ielem, 11, 1) =  1.0 - (hgconst12th[1] * gradop12(ielem, 0, 3) + hgconst12th[5] * gradop12(ielem, 1, 3) + hgconst12th[ 9] * gradop12(ielem, 2, 3));
    hgop(ielem, 12, 1) = -1.0 - (hgconst12th[1] * gradop12(ielem, 0, 4) + hgconst12th[5] * gradop12(ielem, 1, 4) + hgconst12th[ 9] * gradop12(ielem, 2, 4));
    hgop(ielem, 13, 1) =  1.0 - (hgconst12th[1] * gradop12(ielem, 0, 5) + hgconst12th[5] * gradop12(ielem, 1, 5) + hgconst12th[ 9] * gradop12(ielem, 2, 5));
    hgop(ielem, 14, 1) =  1.0 - (hgconst12th[1] * gradop12(ielem, 0, 6) + hgconst12th[5] * gradop12(ielem, 1, 6) + hgconst12th[ 9] * gradop12(ielem, 2, 6));
    hgop(ielem, 15, 1) = -1.0 - (hgconst12th[1] * gradop12(ielem, 0, 7) + hgconst12th[5] * gradop12(ielem, 1, 7) + hgconst12th[ 9] * gradop12(ielem, 2, 7));
    hgop(ielem, 16, 1) =  1.0 - (hgconst12th[2] * gradop12(ielem, 0, 0) + hgconst12th[6] * gradop12(ielem, 1, 0) + hgconst12th[10] * gradop12(ielem, 2, 0));
    hgop(ielem, 17, 1) = -1.0 - (hgconst12th[2] * gradop12(ielem, 0, 1) + hgconst12th[6] * gradop12(ielem, 1, 1) + hgconst12th[10] * gradop12(ielem, 2, 1));
    hgop(ielem, 18, 1) =  1.0 - (hgconst12th[2] * gradop12(ielem, 0, 2) + hgconst12th[6] * gradop12(ielem, 1, 2) + hgconst12th[10] * gradop12(ielem, 2, 2));
    hgop(ielem, 19, 1) = -1.0 - (hgconst12th[2] * gradop12(ielem, 0, 3) + hgconst12th[6] * gradop12(ielem, 1, 3) + hgconst12th[10] * gradop12(ielem, 2, 3));
    hgop(ielem, 20, 1) =  1.0 - (hgconst12th[2] * gradop12(ielem, 0, 4) + hgconst12th[6] * gradop12(ielem, 1, 4) + hgconst12th[10] * gradop12(ielem, 2, 4));
    hgop(ielem, 21, 1) = -1.0 - (hgconst12th[2] * gradop12(ielem, 0, 5) + hgconst12th[6] * gradop12(ielem, 1, 5) + hgconst12th[10] * gradop12(ielem, 2, 5));
    hgop(ielem, 22, 1) =  1.0 - (hgconst12th[2] * gradop12(ielem, 0, 6) + hgconst12th[6] * gradop12(ielem, 1, 6) + hgconst12th[10] * gradop12(ielem, 2, 6));
    hgop(ielem, 23, 1) = -1.0 - (hgconst12th[2] * gradop12(ielem, 0, 7) + hgconst12th[6] * gradop12(ielem, 1, 7) + hgconst12th[10] * gradop12(ielem, 2, 7));
    hgop(ielem, 24, 1) = -1.0 - (hgconst12th[3] * gradop12(ielem, 0, 0) + hgconst12th[7] * gradop12(ielem, 1, 0) + hgconst12th[11] * gradop12(ielem, 2, 0));
    hgop(ielem, 25, 1) =  1.0 - (hgconst12th[3] * gradop12(ielem, 0, 1) + hgconst12th[7] * gradop12(ielem, 1, 1) + hgconst12th[11] * gradop12(ielem, 2, 1));
    hgop(ielem, 26, 1) = -1.0 - (hgconst12th[3] * gradop12(ielem, 0, 2) + hgconst12th[7] * gradop12(ielem, 1, 2) + hgconst12th[11] * gradop12(ielem, 2, 2));
    hgop(ielem, 27, 1) =  1.0 - (hgconst12th[3] * gradop12(ielem, 0, 3) + hgconst12th[7] * gradop12(ielem, 1, 3) + hgconst12th[11] * gradop12(ielem, 2, 3));
    hgop(ielem, 28, 1) =  1.0 - (hgconst12th[3] * gradop12(ielem, 0, 4) + hgconst12th[7] * gradop12(ielem, 1, 4) + hgconst12th[11] * gradop12(ielem, 2, 4));
    hgop(ielem, 29, 1) = -1.0 - (hgconst12th[3] * gradop12(ielem, 0, 5) + hgconst12th[7] * gradop12(ielem, 1, 5) + hgconst12th[11] * gradop12(ielem, 2, 5));
    hgop(ielem, 30, 1) =  1.0 - (hgconst12th[3] * gradop12(ielem, 0, 6) + hgconst12th[7] * gradop12(ielem, 1, 6) + hgconst12th[11] * gradop12(ielem, 2, 6));
    hgop(ielem, 31, 1) = -1.0 - (hgconst12th[3] * gradop12(ielem, 0, 7) + hgconst12th[7] * gradop12(ielem, 1, 7) + hgconst12th[11] * gradop12(ielem, 2, 7));

    }

  KOKKOSARRAY_INLINE_FUNCTION
  void rotate_tensor_backward(int ielem)const {

    //   t : temporary variables
    //   s_n : stress_new in local memory space
    //   r_n : rotation_new in local memory space
      Scalar t[9], s_n[6], r_n[9];

    s_n[0] = stress_new(ielem, 0);
    s_n[1] = stress_new(ielem, 1);
    s_n[2] = stress_new(ielem, 2);
    s_n[3] = stress_new(ielem, 3);
    s_n[4] = stress_new(ielem, 4);
    s_n[5] = stress_new(ielem, 5);

    r_n[0] = rotation(ielem, 0, 1);
    r_n[1] = rotation(ielem, 1, 1);
    r_n[2] = rotation(ielem, 2, 1);
    r_n[3] = rotation(ielem, 3, 1);
    r_n[4] = rotation(ielem, 4, 1);
    r_n[5] = rotation(ielem, 5, 1);
    r_n[6] = rotation(ielem, 6, 1);
    r_n[7] = rotation(ielem, 7, 1);
    r_n[8] = rotation(ielem, 8, 1);

    t[0] = s_n[K_S_XX]*r_n[K_F_XX]+ s_n[K_S_XY]*r_n[K_F_XY]+ s_n[K_S_XZ]*r_n[K_F_XZ];
    t[1] = s_n[K_S_YX]*r_n[K_F_XX]+ s_n[K_S_YY]*r_n[K_F_XY]+ s_n[K_S_YZ]*r_n[K_F_XZ];
    t[2] = s_n[K_S_ZX]*r_n[K_F_XX]+ s_n[K_S_ZY]*r_n[K_F_XY]+ s_n[K_S_ZZ]*r_n[K_F_XZ];
    t[3] = s_n[K_S_XX]*r_n[K_F_YX]+ s_n[K_S_XY]*r_n[K_F_YY]+ s_n[K_S_XZ]*r_n[K_F_YZ];
    t[4] = s_n[K_S_YX]*r_n[K_F_YX]+ s_n[K_S_YY]*r_n[K_F_YY]+ s_n[K_S_YZ]*r_n[K_F_YZ];
    t[5] = s_n[K_S_ZX]*r_n[K_F_YX]+ s_n[K_S_ZY]*r_n[K_F_YY]+ s_n[K_S_ZZ]*r_n[K_F_YZ];
    t[6] = s_n[K_S_XX]*r_n[K_F_ZX]+ s_n[K_S_XY]*r_n[K_F_ZY]+ s_n[K_S_XZ]*r_n[K_F_ZZ];
    t[7] = s_n[K_S_YX]*r_n[K_F_ZX]+ s_n[K_S_YY]*r_n[K_F_ZY]+ s_n[K_S_YZ]*r_n[K_F_ZZ];
    t[8] = s_n[K_S_ZX]*r_n[K_F_ZX]+ s_n[K_S_ZY]*r_n[K_F_ZY]+ s_n[K_S_ZZ]*r_n[K_F_ZZ];

    rot_stress(ielem, K_S_XX) = r_n[K_F_XX]*t[0] + r_n[K_F_XY]*t[1] + r_n[K_F_XZ]*t[2];
    rot_stress(ielem, K_S_YY) = r_n[K_F_YX]*t[3] + r_n[K_F_YY]*t[4] + r_n[K_F_YZ]*t[5];
    rot_stress(ielem, K_S_ZZ) = r_n[K_F_ZX]*t[6] + r_n[K_F_ZY]*t[7] + r_n[K_F_ZZ]*t[8];

    rot_stress(ielem, K_S_XY) = r_n[K_F_XX]*t[3] + r_n[K_F_XY]*t[4] + r_n[K_F_XZ]*t[5];
    rot_stress(ielem, K_S_YZ) = r_n[K_F_YX]*t[6] + r_n[K_F_YY]*t[7] + r_n[K_F_YZ]*t[8];
    rot_stress(ielem, K_S_ZX) = r_n[K_F_ZX]*t[0] + r_n[K_F_ZY]*t[1] + r_n[K_F_ZZ]*t[2];

  }

  KOKKOSARRAY_INLINE_FUNCTION
  Scalar comp_aspect(int ielem)const {

    return 6.0 * 12.0 *  mid_vol(ielem) /
        (gradop12(ielem, 0, 0) * gradop12(ielem, 0, 0)+
         gradop12(ielem, 0, 1) * gradop12(ielem, 0, 1)+
         gradop12(ielem, 0, 2) * gradop12(ielem, 0, 2)+
         gradop12(ielem, 0, 3) * gradop12(ielem, 0, 3)+
         gradop12(ielem, 0, 4) * gradop12(ielem, 0, 4)+
         gradop12(ielem, 0, 5) * gradop12(ielem, 0, 5)+
         gradop12(ielem, 0, 6) * gradop12(ielem, 0, 6)+
         gradop12(ielem, 0, 7) * gradop12(ielem, 0, 7)+
         gradop12(ielem, 1, 0) * gradop12(ielem, 1, 0)+
         gradop12(ielem, 1, 1) * gradop12(ielem, 1, 1)+
         gradop12(ielem, 1, 2) * gradop12(ielem, 1, 2)+
         gradop12(ielem, 1, 3) * gradop12(ielem, 1, 3)+
         gradop12(ielem, 1, 4) * gradop12(ielem, 1, 4)+
         gradop12(ielem, 1, 5) * gradop12(ielem, 1, 5)+
         gradop12(ielem, 1, 6) * gradop12(ielem, 1, 6)+
         gradop12(ielem, 1, 7) * gradop12(ielem, 1, 7)+
         gradop12(ielem, 2, 0) * gradop12(ielem, 2, 0)+
         gradop12(ielem, 2, 1) * gradop12(ielem, 2, 1)+
         gradop12(ielem, 2, 2) * gradop12(ielem, 2, 2)+
         gradop12(ielem, 2, 3) * gradop12(ielem, 2, 3)+
         gradop12(ielem, 2, 4) * gradop12(ielem, 2, 4)+
         gradop12(ielem, 2, 5) * gradop12(ielem, 2, 5)+
         gradop12(ielem, 2, 6) * gradop12(ielem, 2, 6)+
         gradop12(ielem, 2, 7) * gradop12(ielem, 2, 7));

  }

  KOKKOSARRAY_INLINE_FUNCTION
  void comp_force(int ielem, int *nodes, Scalar fac1, Scalar fac2, Scalar * total_stress12th)const {


    //  NKC, does Presto Scalar need these spin rate terms?  Pronto appears to have dumped them....
    const Scalar dwxy = *dt * vorticity(ielem, 0);
    const Scalar dwyz = *dt * vorticity(ielem, 1);
    const Scalar dwzx = *dt * vorticity(ielem, 2);

    //  Compute new hourglass resitance by the old rotated hourglass resitance plus a hourglass rate term
    Scalar hg_resist_total[12];
    Scalar hg_temp[8];

    const int X = 0;
    const int Y = 1;
    const int Z = 2;


    Scalar vx[8], vy[8], vz[8];

    vx[0] = velocity(nodes[0], X, current_state);
    vx[1] = velocity(nodes[1], X, current_state);
    vx[2] = velocity(nodes[2], X, current_state);
    vx[3] = velocity(nodes[3], X, current_state);
    vx[4] = velocity(nodes[4], X, current_state);
    vx[5] = velocity(nodes[5], X, current_state);
    vx[6] = velocity(nodes[6], X, current_state);
    vx[7] = velocity(nodes[7], X, current_state);

    vy[0] = velocity(nodes[0], Y, current_state);
    vy[1] = velocity(nodes[1], Y, current_state);
    vy[2] = velocity(nodes[2], Y, current_state);
    vy[3] = velocity(nodes[3], Y, current_state);
    vy[4] = velocity(nodes[4], Y, current_state);
    vy[5] = velocity(nodes[5], Y, current_state);
    vy[6] = velocity(nodes[6], Y, current_state);
    vy[7] = velocity(nodes[7], Y, current_state);

    vz[0] = velocity(nodes[0], Z, current_state);
    vz[1] = velocity(nodes[1], Z, current_state);
    vz[2] = velocity(nodes[2], Z, current_state);
    vz[3] = velocity(nodes[3], Z, current_state);
    vz[4] = velocity(nodes[4], Z, current_state);
    vz[5] = velocity(nodes[5], Z, current_state);
    vz[6] = velocity(nodes[6], Z, current_state);
    vz[7] = velocity(nodes[7], Z, current_state);


    for(int i = 0; i < 4; ++i) {

      Scalar hg_rate_0 = 0, hg_rate_1 = 0, hg_rate_2 = 0;
      Scalar hg_resist_old_0 = hg_resist(ielem, (3 * i) + 0, 0);
      Scalar hg_resist_old_1 = hg_resist(ielem, (3 * i) + 1, 0);
      Scalar hg_resist_old_2 = hg_resist(ielem, (3 * i) + 2, 0);

      for(int j = 0; j < 8; j++){

        hg_temp[j] = hgop(ielem, i * 8 + j, 1);

      }

      for(int j = 0; j < 8; j++){

        hg_rate_0 += hg_temp[j] * vx[j];
        hg_rate_1 += hg_temp[j] * vy[j];
        hg_rate_2 += hg_temp[j] * vz[j];

      }

      hg_resist(ielem, (i * 3) + 0, 1)   = hg_resist_old_0 + dwxy*hg_resist_old_1 - dwzx*hg_resist_old_2 + fac1* hg_rate_0 ;
      hg_resist_total[(i * 3) + 0]      = hg_resist(ielem, (i * 3) + 0, 1) + fac2* hg_rate_0 ;

      hg_resist(ielem, (i * 3) + 1, 1)   = hg_resist_old_1 - dwxy*hg_resist_old_0 + dwyz*hg_resist_old_2 + fac1* hg_rate_1 ;
      hg_resist_total[(i * 3) + 1]       = hg_resist(ielem, (i * 3) + 1, 1) + fac2* hg_rate_1 ;

      hg_resist(ielem, (i * 3) + 2, 1)   = hg_resist_old_2 + dwzx*hg_resist_old_0 - dwyz*hg_resist_old_1 + fac1* hg_rate_2 ;
      hg_resist_total[(i * 3) + 2]       = hg_resist(ielem, (i * 3) + 2, 1) + fac2* hg_rate_2 ;

    }



    Scalar hg_force_0[8];
    Scalar hg_force_1[8];
    Scalar hg_force_2[8];

    hg_energy(ielem) = 0.0;
    internal_energy(ielem) = 0.0;

    for(int i = 0; i < 8; ++i) {

      hg_force_0[i] =
         (hg_resist_total[HG_X1] * hgop(ielem, i +  0, 0) +
        hg_resist_total[HG_X2] * hgop(ielem, i +  8, 0) +
        hg_resist_total[HG_X3] * hgop(ielem, i + 16, 0) +
        hg_resist_total[HG_X4] * hgop(ielem, i + 24, 0));

      hg_force_1[i] =
         (hg_resist_total[HG_Y1] * hgop(ielem, i +  0, 0) +
        hg_resist_total[HG_Y2] * hgop(ielem, i +  8, 0) +
        hg_resist_total[HG_Y3] * hgop(ielem, i + 16, 0) +
        hg_resist_total[HG_Y4] * hgop(ielem, i + 24, 0));

      hg_force_2[i] =
         (hg_resist_total[HG_Z1] * hgop(ielem, i +  0, 0) +
        hg_resist_total[HG_Z2] * hgop(ielem, i +  8, 0) +
        hg_resist_total[HG_Z3] * hgop(ielem, i + 16, 0) +
        hg_resist_total[HG_Z4] * hgop(ielem, i + 24, 0));

#if 0
// hg_force seems to be always zero
if ( 0 < std::fabs( hg_force_0[i] ) ||
     0 < std::fabs( hg_force_1[i] ) ||
     0 < std::fabs( hg_force_2[i] ) ) {
std::cout
    << "hg_force @ (" << ielem << "," << i << ") ="
    << " " << hg_force_0[i]
    << " " << hg_force_1[i]
    << " " << hg_force_2[i]
    << std::endl ;
}
#endif

    }

    for(int inode = 0; inode < 8; ++inode) {
      element_force(ielem, 0, inode) =
        total_stress12th[K_S_XX] * gradop12(ielem, 0, inode) +
        total_stress12th[K_S_XY] * gradop12(ielem, 1, inode) +
        total_stress12th[K_S_XZ] * gradop12(ielem, 2, inode) + hg_force_0[inode] ;

      element_force(ielem, 1, inode) =
        total_stress12th[K_S_YX] * gradop12(ielem, 0, inode) +
        total_stress12th[K_S_YY] * gradop12(ielem, 1, inode) +
        total_stress12th[K_S_YZ] * gradop12(ielem, 2, inode) + hg_force_1[inode] ;

      element_force(ielem, 2, inode) =
        total_stress12th[K_S_ZX] * gradop12(ielem, 0, inode) +
        total_stress12th[K_S_ZY] * gradop12(ielem, 1, inode) +
        total_stress12th[K_S_ZZ] * gradop12(ielem, 2, inode) + hg_force_2[inode] ;

      hg_energy(ielem)  +=  hg_force_0[inode]   * vx[inode] +
                  hg_force_1[inode]   *vy[inode] +
                  hg_force_2[inode]   *vz[inode];

      internal_energy(ielem) +=  element_force(ielem, 0, inode)*vx[inode] +
                  element_force(ielem, 1, inode)*vy[inode] +
                  element_force(ielem, 2, inode)*vz[inode];

    }
  }

  KOKKOSARRAY_INLINE_FUNCTION
    void get_stress(int ielem) const
    {
      const int kxx = 0;
      const int kyy = 1;
      const int kzz = 2;
      const int kxy = 3;
      const int kyz = 4;
      const int kzx = 5;

      const Scalar e = (rot_stretch(ielem,kxx)+rot_stretch(ielem,kyy)+rot_stretch(ielem,kzz))/3.0;

      stress_new(ielem,kxx) += *dt * (two_mu * (rot_stretch(ielem,kxx)-e)+3*bulk_modulus*e);
      stress_new(ielem,kyy) += *dt * (two_mu * (rot_stretch(ielem,kyy)-e)+3*bulk_modulus*e);
      stress_new(ielem,kzz) += *dt * (two_mu * (rot_stretch(ielem,kzz)-e)+3*bulk_modulus*e);

      stress_new(ielem,kxy) += *dt * two_mu * rot_stretch(ielem,kxy);
      stress_new(ielem,kyz) += *dt * two_mu * rot_stretch(ielem,kyz);
      stress_new(ielem,kzx) += *dt * two_mu * rot_stretch(ielem,kzx);

    }

  KOKKOSARRAY_INLINE_FUNCTION
    void operator()( int ielem, value_type & update )const {

    Scalar x[8], y[8], z[8];
    int nodes[8];

    get_nodes(ielem,nodes);

    comp_grad(ielem,nodes,x,y,z);
    comp_hgop(ielem,x,y,z);

    Scalar fac1_pre = *dt * hg_stiffness * 0.0625;
    Scalar shr = elem_shrmod(ielem) = two_mu;
    Scalar dil = elem_dilmod(ielem) =  bulk_modulus + ((2.0*shr)/3.0);


    Scalar aspect = comp_aspect(ielem);

    Scalar inv_aspect = 1.0 / aspect;

    Scalar dtrial = sqrt(elem_mass(ielem) * aspect / dil);
    Scalar traced = (rot_stretch(ielem, 0) + rot_stretch(ielem, 1) + rot_stretch(ielem, 2));

    Scalar eps = traced < 0 ? (lin_bulk_visc - quad_bulk_visc * traced * dtrial) : lin_bulk_visc ;

    Scalar bulkq = eps * dil * dtrial * traced;

    Scalar cur_time_step = dtrial * ( sqrt( 1.0 + eps * eps) - eps);

    //force fix time step
    cur_time_step = user_dt > 0 ? user_dt : cur_time_step;

    update = update < cur_time_step ? update : cur_time_step;

    //elem_t_step(ielem) = cur_time_step;

    get_stress(ielem);

    rotate_tensor_backward(ielem);

    Scalar total_stress12th[6];
    total_stress12th[0] = ONE12TH*(rot_stress(ielem, 0) + bulkq);
    total_stress12th[1] = ONE12TH*(rot_stress(ielem, 1) + bulkq);
    total_stress12th[2] = ONE12TH*(rot_stress(ielem, 2) + bulkq);
    total_stress12th[3] = ONE12TH*(rot_stress(ielem, 3));
    total_stress12th[4] = ONE12TH*(rot_stress(ielem, 4));
    total_stress12th[5] = ONE12TH*(rot_stress(ielem, 5));


    Scalar fac1 = fac1_pre * shr * inv_aspect;
    Scalar fac2 = hg_viscosity * sqrt(shr * elem_mass(ielem) * inv_aspect);

    comp_force(ielem, nodes, fac1, fac2, total_stress12th);

  }

};

template<typename Scalar , class DeviceType>
struct set_next_time_step;

template<typename Scalar>
struct set_next_time_step<Scalar ,KOKKOSARRAY_MACRO_DEVICE>{

  typedef KOKKOSARRAY_MACRO_DEVICE       device_type;
  typedef device_type::size_type    size_type;

  typedef Scalar value_type;

  typedef Region<Scalar,device_type> MyRegion;

    set_next_time_step( const MyRegion  & arg_region )
       : region(arg_region)
      {}


    KOKKOSARRAY_INLINE_FUNCTION
    void operator()(Scalar & result) const {
      *(region.prev_dt) = *(region.dt);
      *(region.dt) = result;
    }

    MyRegion   region;

}; //minimum_stable_time_step

#endif
