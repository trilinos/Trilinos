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

#ifndef GRAD_HGOP
#define GRAD_HGOP

//#define ONE12TH 0.083333333333333333333333
#define ONE12TH (1.0/12.0)
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
struct grad_hgop;

template<typename Scalar>
struct grad_hgop<Scalar, KOKKOSARRAY_MACRO_DEVICE>{

  typedef KOKKOSARRAY_MACRO_DEVICE     device_type ;
  typedef typename KokkosArray::MDArray<Scalar,device_type> array_type ;
  typedef typename KokkosArray::MDArray<int,device_type>    int_array_type ;

  typedef KokkosArray::Value<Scalar,device_type>     scalar;

  // Global arrays used by this functor.

  const int_array_type elem_node_connectivity;

  const array_type model_coords;
  const array_type displacement;
  const array_type velocity;
  const array_type vel_grad;
  const array_type hgop;

	const scalar     dt;
  const int        current_state;
  const int        previous_state;

  typedef Region<Scalar,device_type> MyRegion;

  // Constructor on the Host to populate this device functor.
  // All array view copies are shallow.
    grad_hgop( const MyRegion & region,
               const int arg_current_state,
               const int arg_previous_state)
    : elem_node_connectivity(region.elem_node_connectivity)
    , model_coords(region.model_coords)
    , displacement(region.displacement)
    , velocity(region.velocity)
    , vel_grad(region.vel_grad)
    , hgop    (region.hgop)
    , dt( region.dt)
    , current_state(arg_current_state)
    , previous_state(arg_previous_state)
    {
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
    void comp_grad( int * nodes,
                    Scalar * x,     Scalar * y,     Scalar * z,
                    Scalar * vx,     Scalar * vy,     Scalar * vz,
                    Scalar * grad_x,  Scalar * grad_y,  Scalar * grad_z) const
    {


    Scalar dt_scale = -0.5 * *dt;

    //enum { X = 0, Y = 1, Z = 2 };
    const int X = 0;
    const int Y = 1;
    const int Z = 2;

  // Read global velocity once and use many times via local registers / L1 cache.
  //  store the velocity information in local memory before using, so it can
  //  be returned for other functions to use

    vx[0] = velocity(nodes[0], X, current_state);
    vx[1] = velocity(nodes[1], X, current_state);
    vx[2] = velocity(nodes[2], X, current_state);
    vx[3] = velocity(nodes[3], X, current_state);
    vx[4] = velocity(nodes[4], X, current_state);
    vx[5] = velocity(nodes[5], X, current_state);
    vx[6] = velocity(nodes[6], X, current_state);
    vx[7] = velocity(nodes[7], X, current_state);

  // Read global coordinates once and use many times via local registers / L1 cache.
  //  load X coordinate information and move by half time step
    x[0] = model_coords(nodes[0], X) + displacement(nodes[0], X, current_state) + dt_scale * vx[0];
    x[1] = model_coords(nodes[1], X) + displacement(nodes[1], X, current_state) + dt_scale * vx[1];
    x[2] = model_coords(nodes[2], X) + displacement(nodes[2], X, current_state) + dt_scale * vx[2];
    x[3] = model_coords(nodes[3], X) + displacement(nodes[3], X, current_state) + dt_scale * vx[3];
    x[4] = model_coords(nodes[4], X) + displacement(nodes[4], X, current_state) + dt_scale * vx[4];
    x[5] = model_coords(nodes[5], X) + displacement(nodes[5], X, current_state) + dt_scale * vx[5];
    x[6] = model_coords(nodes[6], X) + displacement(nodes[6], X, current_state) + dt_scale * vx[6];
    x[7] = model_coords(nodes[7], X) + displacement(nodes[7], X, current_state) + dt_scale * vx[7];

  //  calc X difference vectors
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

    vz[0] = velocity(nodes[0], Z, current_state);
    vz[1] = velocity(nodes[1], Z, current_state);
    vz[2] = velocity(nodes[2], Z, current_state);
    vz[3] = velocity(nodes[3], Z, current_state);
    vz[4] = velocity(nodes[4], Z, current_state);
    vz[5] = velocity(nodes[5], Z, current_state);
    vz[6] = velocity(nodes[6], Z, current_state);
    vz[7] = velocity(nodes[7], Z, current_state);

  //  Load Z information
    z[0] = model_coords(nodes[0], Z) + displacement(nodes[0], Z, current_state) + dt_scale * vz[0];
    z[1] = model_coords(nodes[1], Z) + displacement(nodes[1], Z, current_state) + dt_scale * vz[1];
    z[2] = model_coords(nodes[2], Z) + displacement(nodes[2], Z, current_state) + dt_scale * vz[2];
    z[3] = model_coords(nodes[3], Z) + displacement(nodes[3], Z, current_state) + dt_scale * vz[3];
    z[4] = model_coords(nodes[4], Z) + displacement(nodes[4], Z, current_state) + dt_scale * vz[4];
    z[5] = model_coords(nodes[5], Z) + displacement(nodes[5], Z, current_state) + dt_scale * vz[5];
    z[6] = model_coords(nodes[6], Z) + displacement(nodes[6], Z, current_state) + dt_scale * vz[6];
    z[7] = model_coords(nodes[7], Z) + displacement(nodes[7], Z, current_state) + dt_scale * vz[7];

  //  Some data is used by other functors, as well as this one, so we make a
  //  local copy, but still write the values to Global memory, so that other
  //  functors can access it. By putting the local variable on the right hand
  //  side, we avoid the possibility that the values would be written to global
  //  memory, then read back from global into the local variables.


  //  Calculate Y gradient from X and Z data
    grad_y[0] = (z[1] *  t1) - (z[2] * R42) - (z[3] *  t5)  + (z[4] *  t4) + (z[5] * R52) - (z[7] * R54);
    grad_y[1] = (z[2] *  t2) + (z[3] * R31) - (z[0] *  t1)  - (z[5] *  t6) + (z[6] * R63) - (z[4] * R61);
    grad_y[2] = (z[3] *  t3) + (z[0] * R42) - (z[1] *  t2)  - (z[6] *  t4) + (z[7] * R74) - (z[5] * R72);
    grad_y[3] = (z[0] *  t5) - (z[1] * R31) - (z[2] *  t3)  + (z[7] *  t6) + (z[4] * R81) - (z[6] * R83);
    grad_y[4] = (z[5] *  t3) + (z[6] * R86) - (z[7] *  t2)  - (z[0] *  t4) - (z[3] * R81) + (z[1] * R61);
    grad_y[5] = (z[6] *  t5) - (z[4] *  t3)  - (z[7] * R75) + (z[1] *  t6) - (z[0] * R52) + (z[2] * R72);
    grad_y[6] = (z[7] *  t1) - (z[5] *  t5)  - (z[4] * R86) + (z[2] *  t4) - (z[1] * R63) + (z[3] * R83);
    grad_y[7] = (z[4] *  t2) - (z[6] *  t1)  + (z[5] * R75) - (z[3] *  t6) - (z[2] * R74) + (z[0] * R54);


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

    vy[0] = velocity(nodes[0], Y, current_state);
    vy[1] = velocity(nodes[1], Y, current_state);
    vy[2] = velocity(nodes[2], Y, current_state);
    vy[3] = velocity(nodes[3], Y, current_state);
    vy[4] = velocity(nodes[4], Y, current_state);
    vy[5] = velocity(nodes[5], Y, current_state);
    vy[6] = velocity(nodes[6], Y, current_state);
    vy[7] = velocity(nodes[7], Y, current_state);

  //  Load Y information
    y[0] = model_coords(nodes[0], Y) + displacement(nodes[0], Y, current_state) + dt_scale * vy[0];
    y[1] = model_coords(nodes[1], Y) + displacement(nodes[1], Y, current_state) + dt_scale * vy[1];
    y[2] = model_coords(nodes[2], Y) + displacement(nodes[2], Y, current_state) + dt_scale * vy[2];
    y[3] = model_coords(nodes[3], Y) + displacement(nodes[3], Y, current_state) + dt_scale * vy[3];
    y[4] = model_coords(nodes[4], Y) + displacement(nodes[4], Y, current_state) + dt_scale * vy[4];
    y[5] = model_coords(nodes[5], Y) + displacement(nodes[5], Y, current_state) + dt_scale * vy[5];
    y[6] = model_coords(nodes[6], Y) + displacement(nodes[6], Y, current_state) + dt_scale * vy[6];
    y[7] = model_coords(nodes[7], Y) + displacement(nodes[7], Y, current_state) + dt_scale * vy[7];


  //  Calculate X gradient from Y and Z data
    grad_x[0] = (y[1] *  t1) - (y[2] * R42) - (y[3] *  t5) + (y[4] *  t4) + (y[5] * R52) - (y[7] * R54);
    grad_x[1] = (y[2] *  t2) + (y[3] * R31) - (y[0] *  t1) - (y[5] *  t6) + (y[6] * R63) - (y[4] * R61);
    grad_x[2] = (y[3] *  t3) + (y[0] * R42) - (y[1] *  t2) - (y[6] *  t4) + (y[7] * R74) - (y[5] * R72);
    grad_x[3] = (y[0] *  t5) - (y[1] * R31) - (y[2] *  t3) + (y[7] *  t6) + (y[4] * R81) - (y[6] * R83);
    grad_x[4] = (y[5] *  t3) + (y[6] * R86) - (y[7] *  t2) - (y[0] *  t4) - (y[3] * R81) + (y[1] * R61);
    grad_x[5] = (y[6] *  t5) - (y[4] *  t3) - (y[7] * R75) + (y[1] *  t6) - (y[0] * R52) + (y[2] * R72);
    grad_x[6] = (y[7] *  t1) - (y[5] *  t5) - (y[4] * R86) + (y[2] *  t4) - (y[1] * R63) + (y[3] * R83);
    grad_x[7] = (y[4] *  t2) - (y[6] *  t1) + (y[5] * R75) - (y[3] *  t6) - (y[2] * R74) + (y[0] * R54);


  //  calc Y difference vectors
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
    grad_z[0] = (x[1] *  t1) - (x[2] * R42) - (x[3] *  t5)  + (x[4] *  t4) + (x[5] * R52) - (x[7] * R54);
    grad_z[1] = (x[2] *  t2) + (x[3] * R31) - (x[0] *  t1)  - (x[5] *  t6) + (x[6] * R63) - (x[4] * R61);
    grad_z[2] = (x[3] *  t3) + (x[0] * R42) - (x[1] *  t2)  - (x[6] *  t4) + (x[7] * R74) - (x[5] * R72);
    grad_z[3] = (x[0] *  t5) - (x[1] * R31) - (x[2] *  t3)  + (x[7] *  t6) + (x[4] * R81) - (x[6] * R83);
    grad_z[4] = (x[5] *  t3) + (x[6] * R86) - (x[7] *  t2)  - (x[0] *  t4) - (x[3] * R81) + (x[1] * R61);
    grad_z[5] = (x[6] *  t5) - (x[4] *  t3)  - (x[7] * R75) + (x[1] *  t6) - (x[0] * R52) + (x[2] * R72);
    grad_z[6] = (x[7] *  t1) - (x[5] *  t5)  - (x[4] * R86) + (x[2] *  t4) - (x[1] * R63) + (x[3] * R83);
    grad_z[7] = (x[4] *  t2) - (x[6] *  t1)  + (x[5] * R75) - (x[3] *  t6) - (x[2] * R74) + (x[0] * R54);

    }

    KOKKOSARRAY_INLINE_FUNCTION
    void v_grad(  int ielem,
      Scalar * vx,       Scalar * vy,       Scalar * vz,
      Scalar * grad_x,     Scalar * grad_y,     Scalar * grad_z,
      Scalar inv_vol_1_12) const {

  //   Calculate Velocity Gradients
    vel_grad(ielem, K_F_XX) = inv_vol_1_12 * (    vx[0] * grad_x[0] +
                        vx[1] * grad_x[1] +
                        vx[2] * grad_x[2] +
                        vx[3] * grad_x[3] +
                        vx[4] * grad_x[4] +
                        vx[5] * grad_x[5] +
                        vx[6] * grad_x[6] +
                        vx[7] * grad_x[7] );

    vel_grad(ielem, K_F_YX) = inv_vol_1_12 * (    vy[0] * grad_x[0] +
                        vy[1] * grad_x[1] +
                        vy[2] * grad_x[2] +
                        vy[3] * grad_x[3] +
                        vy[4] * grad_x[4] +
                        vy[5] * grad_x[5] +
                        vy[6] * grad_x[6] +
                        vy[7] * grad_x[7] );

    vel_grad(ielem, K_F_ZX) = inv_vol_1_12 * (    vz[0] * grad_x[0] +
                        vz[1] * grad_x[1] +
                        vz[2] * grad_x[2] +
                        vz[3] * grad_x[3] +
                        vz[4] * grad_x[4] +
                        vz[5] * grad_x[5] +
                        vz[6] * grad_x[6] +
                        vz[7] * grad_x[7] );



    vel_grad(ielem, K_F_XY) = inv_vol_1_12 * (    vx[0] * grad_y[0] +
                        vx[1] * grad_y[1] +
                        vx[2] * grad_y[2] +
                        vx[3] * grad_y[3] +
                        vx[4] * grad_y[4] +
                        vx[5] * grad_y[5] +
                        vx[6] * grad_y[6] +
                        vx[7] * grad_y[7] );

    vel_grad(ielem, K_F_YY) = inv_vol_1_12 * (    vy[0] * grad_y[0] +
                        vy[1] * grad_y[1] +
                        vy[2] * grad_y[2] +
                        vy[3] * grad_y[3] +
                        vy[4] * grad_y[4] +
                        vy[5] * grad_y[5] +
                        vy[6] * grad_y[6] +
                        vy[7] * grad_y[7] );

    vel_grad(ielem, K_F_ZY) = inv_vol_1_12 * (    vz[0] * grad_y[0] +
                        vz[1] * grad_y[1] +
                        vz[2] * grad_y[2] +
                        vz[3] * grad_y[3] +
                        vz[4] * grad_y[4] +
                        vz[5] * grad_y[5] +
                        vz[6] * grad_y[6] +
                        vz[7] * grad_y[7] );



    vel_grad(ielem, K_F_XZ) = inv_vol_1_12 * (    vx[0] * grad_z[0] +
                        vx[1] * grad_z[1] +
                        vx[2] * grad_z[2] +
                        vx[3] * grad_z[3] +
                        vx[4] * grad_z[4] +
                        vx[5] * grad_z[5] +
                        vx[6] * grad_z[6] +
                        vx[7] * grad_z[7] );

    vel_grad(ielem, K_F_YZ) = inv_vol_1_12 * (    vy[0] * grad_z[0] +
                        vy[1] * grad_z[1] +
                        vy[2] * grad_z[2] +
                        vy[3] * grad_z[3] +
                        vy[4] * grad_z[4] +
                        vy[5] * grad_z[5] +
                        vy[6] * grad_z[6] +
                        vy[7] * grad_z[7] );

    vel_grad(ielem, K_F_ZZ) = inv_vol_1_12 * (    vz[0] * grad_z[0] +
                        vz[1] * grad_z[1] +
                        vz[2] * grad_z[2] +
                        vz[3] * grad_z[3] +
                        vz[4] * grad_z[4] +
                        vz[5] * grad_z[5] +
                        vz[6] * grad_z[6] +
                        vz[7] * grad_z[7] );

#if 0
    std::cout << "Element " << ielem << " vx: ";
    for(int i = 0; i<8; ++i) {
      std::cout << vx[i] << ",";
    }
    std::cout << std::endl;

    std::cout << "Element " << ielem << " grad_x: ";
    for(int i = 0; i<8; ++i) {
      std::cout << grad_x[i] << ",";
    }
    std::cout << std::endl;

    std::cout << "Element " << ielem << " vel_grad(x): ";
    std::cout << vel_grad(ielem,0) << std::endl;
#endif

  }


  KOKKOSARRAY_INLINE_FUNCTION
    void comp_hgop(    int ielem,
            Scalar * x,     Scalar * y,     Scalar * z,
            Scalar * grad_x,   Scalar * grad_y,   Scalar * grad_z,
            Scalar inv_vol_1_12) const {

  //   KHP: Alternatively, we could have
  //   hx0,hx1,hx2,hx3,...,hz0,hz1,hz2,hz3
    Scalar hgconst12th[12];

    Scalar q0 = x[0] - x[1];
    Scalar q1 = x[2] - x[3];
    Scalar q2 = x[4] - x[5];
    Scalar q3 = x[6] - x[7];

    hgconst12th[0] = ( (x[0]+x[1]) - (x[2]+x[3]) - (x[4]+x[5]) + (x[6]+x[7]) ) * inv_vol_1_12;
    hgconst12th[1] = (  q0 - q1 - q2 + q3 ) * inv_vol_1_12;
    hgconst12th[2] = (  q0 + q1 + q2 + q3 ) * inv_vol_1_12;
    hgconst12th[3] = ( -q0 - q1 + q2 + q3 ) * inv_vol_1_12;

    q0 = (y[0] - y[1]);
    q1 = (y[2] - y[3]);
    q2 = (y[4] - y[5]);
    q3 = (y[6] - y[7]);

    hgconst12th[4] = ( (y[0]+y[1]) - (y[2]+y[3]) - (y[4]+y[5]) + (y[6]+y[7]) ) * inv_vol_1_12;
    hgconst12th[5] = (  q0 - q1 - q2 + q3 ) * inv_vol_1_12;
    hgconst12th[6] = (  q0 + q1 + q2 + q3 ) * inv_vol_1_12;
    hgconst12th[7] = ( -q0 - q1 + q2 + q3 ) * inv_vol_1_12;

    q0 = (z[0] - z[1]);
    q1 = (z[2] - z[3]);
    q2 = (z[4] - z[5]);
    q3 = (z[6] - z[7]);

    hgconst12th[8]  = ( (z[0]+z[1]) - (z[2]+z[3]) - (z[4]+z[5]) + (z[6]+z[7]) ) * inv_vol_1_12;
    hgconst12th[9]  = (  q0 - q1 - q2 + q3 ) * inv_vol_1_12;
    hgconst12th[10] = (  q0 + q1 + q2 + q3 ) * inv_vol_1_12;
    hgconst12th[11] = ( -q0 - q1 + q2 + q3 ) * inv_vol_1_12;


  //  In the original code, an array of 32 doubles was used to store the
  //  constant coefficients in the following calculation. By removing that
  //  array, the memory footprint for each thread is smaller, which allows
  //  higher occupancy, and frees up register space for other variables.

  //  The array of +1/-1 coefficient could be placed in constant or global memory
  //  and "re-roll" this code into a loop.

    hgop(ielem,  0, 1) =  1.0 - (hgconst12th[0] * grad_x[0] + hgconst12th[4] * grad_y[0] + hgconst12th[ 8] * grad_z[0]);
    hgop(ielem,  1, 1) =  1.0 - (hgconst12th[0] * grad_x[1] + hgconst12th[4] * grad_y[1] + hgconst12th[ 8] * grad_z[1]);
    hgop(ielem,  2, 1) = -1.0 - (hgconst12th[0] * grad_x[2] + hgconst12th[4] * grad_y[2] + hgconst12th[ 8] * grad_z[2]);
    hgop(ielem,  3, 1) = -1.0 - (hgconst12th[0] * grad_x[3] + hgconst12th[4] * grad_y[3] + hgconst12th[ 8] * grad_z[3]);
    hgop(ielem,  4, 1) = -1.0 - (hgconst12th[0] * grad_x[4] + hgconst12th[4] * grad_y[4] + hgconst12th[ 8] * grad_z[4]);
    hgop(ielem,  5, 1) = -1.0 - (hgconst12th[0] * grad_x[5] + hgconst12th[4] * grad_y[5] + hgconst12th[ 8] * grad_z[5]);
    hgop(ielem,  6, 1) =  1.0 - (hgconst12th[0] * grad_x[6] + hgconst12th[4] * grad_y[6] + hgconst12th[ 8] * grad_z[6]);
    hgop(ielem,  7, 1) =  1.0 - (hgconst12th[0] * grad_x[7] + hgconst12th[4] * grad_y[7] + hgconst12th[ 8] * grad_z[7]);
    hgop(ielem,  8, 1) =  1.0 - (hgconst12th[1] * grad_x[0] + hgconst12th[5] * grad_y[0] + hgconst12th[ 9] * grad_z[0]);
    hgop(ielem,  9, 1) = -1.0 - (hgconst12th[1] * grad_x[1] + hgconst12th[5] * grad_y[1] + hgconst12th[ 9] * grad_z[1]);
    hgop(ielem, 10, 1) = -1.0 - (hgconst12th[1] * grad_x[2] + hgconst12th[5] * grad_y[2] + hgconst12th[ 9] * grad_z[2]);
    hgop(ielem, 11, 1) =  1.0 - (hgconst12th[1] * grad_x[3] + hgconst12th[5] * grad_y[3] + hgconst12th[ 9] * grad_z[3]);
    hgop(ielem, 12, 1) = -1.0 - (hgconst12th[1] * grad_x[4] + hgconst12th[5] * grad_y[4] + hgconst12th[ 9] * grad_z[4]);
    hgop(ielem, 13, 1) =  1.0 - (hgconst12th[1] * grad_x[5] + hgconst12th[5] * grad_y[5] + hgconst12th[ 9] * grad_z[5]);
    hgop(ielem, 14, 1) =  1.0 - (hgconst12th[1] * grad_x[6] + hgconst12th[5] * grad_y[6] + hgconst12th[ 9] * grad_z[6]);
    hgop(ielem, 15, 1) = -1.0 - (hgconst12th[1] * grad_x[7] + hgconst12th[5] * grad_y[7] + hgconst12th[ 9] * grad_z[7]);
    hgop(ielem, 16, 1) =  1.0 - (hgconst12th[2] * grad_x[0] + hgconst12th[6] * grad_y[0] + hgconst12th[10] * grad_z[0]);
    hgop(ielem, 17, 1) = -1.0 - (hgconst12th[2] * grad_x[1] + hgconst12th[6] * grad_y[1] + hgconst12th[10] * grad_z[1]);
    hgop(ielem, 18, 1) =  1.0 - (hgconst12th[2] * grad_x[2] + hgconst12th[6] * grad_y[2] + hgconst12th[10] * grad_z[2]);
    hgop(ielem, 19, 1) = -1.0 - (hgconst12th[2] * grad_x[3] + hgconst12th[6] * grad_y[3] + hgconst12th[10] * grad_z[3]);
    hgop(ielem, 20, 1) =  1.0 - (hgconst12th[2] * grad_x[4] + hgconst12th[6] * grad_y[4] + hgconst12th[10] * grad_z[4]);
    hgop(ielem, 21, 1) = -1.0 - (hgconst12th[2] * grad_x[5] + hgconst12th[6] * grad_y[5] + hgconst12th[10] * grad_z[5]);
    hgop(ielem, 22, 1) =  1.0 - (hgconst12th[2] * grad_x[6] + hgconst12th[6] * grad_y[6] + hgconst12th[10] * grad_z[6]);
    hgop(ielem, 23, 1) = -1.0 - (hgconst12th[2] * grad_x[7] + hgconst12th[6] * grad_y[7] + hgconst12th[10] * grad_z[7]);
    hgop(ielem, 24, 1) = -1.0 - (hgconst12th[3] * grad_x[0] + hgconst12th[7] * grad_y[0] + hgconst12th[11] * grad_z[0]);
    hgop(ielem, 25, 1) =  1.0 - (hgconst12th[3] * grad_x[1] + hgconst12th[7] * grad_y[1] + hgconst12th[11] * grad_z[1]);
    hgop(ielem, 26, 1) = -1.0 - (hgconst12th[3] * grad_x[2] + hgconst12th[7] * grad_y[2] + hgconst12th[11] * grad_z[2]);
    hgop(ielem, 27, 1) =  1.0 - (hgconst12th[3] * grad_x[3] + hgconst12th[7] * grad_y[3] + hgconst12th[11] * grad_z[3]);
    hgop(ielem, 28, 1) =  1.0 - (hgconst12th[3] * grad_x[4] + hgconst12th[7] * grad_y[4] + hgconst12th[11] * grad_z[4]);
    hgop(ielem, 29, 1) = -1.0 - (hgconst12th[3] * grad_x[5] + hgconst12th[7] * grad_y[5] + hgconst12th[11] * grad_z[5]);
    hgop(ielem, 30, 1) =  1.0 - (hgconst12th[3] * grad_x[6] + hgconst12th[7] * grad_y[6] + hgconst12th[11] * grad_z[6]);
    hgop(ielem, 31, 1) = -1.0 - (hgconst12th[3] * grad_x[7] + hgconst12th[7] * grad_y[7] + hgconst12th[11] * grad_z[7]);

  }

  // Functor operator() which calls the three member functions.

  KOKKOSARRAY_INLINE_FUNCTION
    void operator()( int ielem )const {

    //  declare and reuse local data for frequently accessed data to
    //  reduce global memory reads and writes.

    Scalar      x[8],      y[8],      z[8];
    Scalar     vx[8],     vy[8],     vz[8];
    Scalar grad_x[8], grad_y[8], grad_z[8];

    int nodes[8];

    get_nodes(ielem,nodes);

    comp_grad(nodes, x, y, z, vx, vy, vz, grad_x, grad_y, grad_z);

  //  Calculate hexahedral volume from x model_coords and gradient information
    Scalar inv_vol_1_12;

    inv_vol_1_12 = 1.0 / ( x[0] * grad_x[0] +
                      x[1] * grad_x[1] +
                      x[2] * grad_x[2] +
                      x[3] * grad_x[3] +
                      x[4] * grad_x[4] +
                      x[5] * grad_x[5] +
                      x[6] * grad_x[6] +
                      x[7] * grad_x[7] );


    v_grad(ielem, vx, vy, vz, grad_x, grad_y, grad_z, inv_vol_1_12);

    comp_hgop(ielem, x, y, z, grad_x, grad_y, grad_z, inv_vol_1_12);

  }

};

#endif
