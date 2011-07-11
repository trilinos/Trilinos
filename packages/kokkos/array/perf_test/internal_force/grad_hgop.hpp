#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#ifndef GRAD_HGOP
#define GRAD_HGOP

#define ONE12TH 0.083333333333333333333333
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
struct grad_hgop<Scalar, KOKKOS_MACRO_DEVICE>{

  typedef KOKKOS_MACRO_DEVICE     device_type ;
    typedef typename Kokkos::MDArrayView<Scalar,device_type> array_type ;

  const array_type position;    
  const array_type velocity;
  const array_type mid_pos;
  const array_type mid_vol;
  const array_type gradop12;
  const array_type vel_grad;
  const array_type hgop;

  const Scalar     dt;

  grad_hgop( const array_type & arg_p,
             const array_type & arg_v,
             const array_type & arg_m,
             const array_type & arg_mv,
             const array_type & arg_g,
             const array_type & arg_vg,
             const array_type & arg_hg,
             const Scalar delta_t )
    : position( arg_p )
    , velocity( arg_v )
    , mid_pos ( arg_m )
    , mid_vol ( arg_mv )
    , gradop12( arg_g )
    , vel_grad( arg_vg )
    , hgop    ( arg_hg )
    , dt( delta_t )
    {}

  KOKKOS_MACRO_DEVICE_FUNCTION
    void comp_grad(   int ielem, 
          Scalar * x,     Scalar * y,     Scalar * z,
          Scalar * vx,     Scalar * vy,     Scalar * vz,
          Scalar * grad_x,  Scalar * grad_y,  Scalar * grad_z) const {


    Scalar dt_scale = -0.5 * dt;
    
    vx[0] = velocity(ielem, 0, 0);
    vx[1] = velocity(ielem, 0, 1);
    vx[2] = velocity(ielem, 0, 2);
    vx[3] = velocity(ielem, 0, 3);
    vx[4] = velocity(ielem, 0, 4);
    vx[5] = velocity(ielem, 0, 5);
    vx[6] = velocity(ielem, 0, 6);
    vx[7] = velocity(ielem, 0, 7);

    //  load X coordinate information
    mid_pos(ielem, 0, 0) = x[0] = position(ielem, 0, 0) + dt_scale * vx[0];
    mid_pos(ielem, 0, 1) = x[1] = position(ielem, 0, 1) + dt_scale * vx[1];
    mid_pos(ielem, 0, 2) = x[2] = position(ielem, 0, 2) + dt_scale * vx[2];
    mid_pos(ielem, 0, 3) = x[3] = position(ielem, 0, 3) + dt_scale * vx[3];
    mid_pos(ielem, 0, 4) = x[4] = position(ielem, 0, 4) + dt_scale * vx[4];
    mid_pos(ielem, 0, 5) = x[5] = position(ielem, 0, 5) + dt_scale * vx[5];
    mid_pos(ielem, 0, 6) = x[6] = position(ielem, 0, 6) + dt_scale * vx[6];
    mid_pos(ielem, 0, 7) = x[7] = position(ielem, 0, 7) + dt_scale * vx[7];

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


    vz[0] = velocity(ielem, 2, 0);
    vz[1] = velocity(ielem, 2, 1);
    vz[2] = velocity(ielem, 2, 2);
    vz[3] = velocity(ielem, 2, 3);
    vz[4] = velocity(ielem, 2, 4);
    vz[5] = velocity(ielem, 2, 5);
    vz[6] = velocity(ielem, 2, 6);
    vz[7] = velocity(ielem, 2, 7);
  
    //  Load Z information
    mid_pos(ielem, 2, 0) = z[0] = position(ielem, 2, 0) + dt_scale * vz[0];
    mid_pos(ielem, 2, 1) = z[1] = position(ielem, 2, 1) + dt_scale * vz[1];
    mid_pos(ielem, 2, 2) = z[2] = position(ielem, 2, 2) + dt_scale * vz[2];
    mid_pos(ielem, 2, 3) = z[3] = position(ielem, 2, 3) + dt_scale * vz[3];
    mid_pos(ielem, 2, 4) = z[4] = position(ielem, 2, 4) + dt_scale * vz[4];
    mid_pos(ielem, 2, 5) = z[5] = position(ielem, 2, 5) + dt_scale * vz[5];
    mid_pos(ielem, 2, 6) = z[6] = position(ielem, 2, 6) + dt_scale * vz[6];
    mid_pos(ielem, 2, 7) = z[7] = position(ielem, 2, 7) + dt_scale * vz[7];


    //  Calculate Y gradient from X and Z data
    grad_y[0] = gradop12(ielem, 1, 0) = (z[1] *  t1) - (z[2] * R42) - (z[3] *  t5)  + (z[4] *  t4) + (z[5] * R52) - (z[7] * R54); 
    grad_y[1] = gradop12(ielem, 1, 1) = (z[2] *  t2) + (z[3] * R31) - (z[0] *  t1)  - (z[5] *  t6) + (z[6] * R63) - (z[4] * R61); 
    grad_y[2] = gradop12(ielem, 1, 2) = (z[3] *  t3) + (z[0] * R42) - (z[1] *  t2)  - (z[6] *  t4) + (z[7] * R74) - (z[5] * R72); 
    grad_y[3] = gradop12(ielem, 1, 3) = (z[0] *  t5) - (z[1] * R31) - (z[2] *  t3)  + (z[7] *  t6) + (z[4] * R81) - (z[6] * R83); 
    grad_y[4] = gradop12(ielem, 1, 4) = (z[5] *  t3) + (z[6] * R86) - (z[7] *  t2)  - (z[0] *  t4) - (z[3] * R81) + (z[1] * R61); 
    grad_y[5] = gradop12(ielem, 1, 5) = (z[6] *  t5) - (z[4] *  t3)  - (z[7] * R75) + (z[1] *  t6) - (z[0] * R52) + (z[2] * R72);
    grad_y[6] = gradop12(ielem, 1, 6) = (z[7] *  t1) - (z[5] *  t5)  - (z[4] * R86) + (z[2] *  t4) - (z[1] * R63) + (z[3] * R83);
    grad_y[7] = gradop12(ielem, 1, 7) = (z[4] *  t2) - (z[6] *  t1)  + (z[5] * R75) - (z[3] *  t6) - (z[2] * R74) + (z[0] * R54); 


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

    vy[0] = velocity(ielem, 1, 0);
    vy[1] = velocity(ielem, 1, 1);
    vy[2] = velocity(ielem, 1, 2);
    vy[3] = velocity(ielem, 1, 3);
    vy[4] = velocity(ielem, 1, 4);
    vy[5] = velocity(ielem, 1, 5);
    vy[6] = velocity(ielem, 1, 6);
    vy[7] = velocity(ielem, 1, 7);
  
    //  Load Y information
    mid_pos(ielem, 1, 0) = y[0] = position(ielem, 1, 0) + dt_scale * vy[0];
    mid_pos(ielem, 1, 1) = y[1] = position(ielem, 1, 1) + dt_scale * vy[1];
    mid_pos(ielem, 1, 2) = y[2] = position(ielem, 1, 2) + dt_scale * vy[2];
    mid_pos(ielem, 1, 3) = y[3] = position(ielem, 1, 3) + dt_scale * vy[3];
    mid_pos(ielem, 1, 4) = y[4] = position(ielem, 1, 4) + dt_scale * vy[4];
    mid_pos(ielem, 1, 5) = y[5] = position(ielem, 1, 5) + dt_scale * vy[5];
    mid_pos(ielem, 1, 6) = y[6] = position(ielem, 1, 6) + dt_scale * vy[6];
    mid_pos(ielem, 1, 7) = y[7] = position(ielem, 1, 7) + dt_scale * vy[7];


    //  Calculate X gradient from Y and Z data
    grad_x[0] = gradop12(ielem, 0, 0) = (y[1] *  t1) - (y[2] * R42) - (y[3] *  t5) + (y[4] *  t4) + (y[5] * R52) - (y[7] * R54); 
    grad_x[1] = gradop12(ielem, 0, 1) = (y[2] *  t2) + (y[3] * R31) - (y[0] *  t1) - (y[5] *  t6) + (y[6] * R63) - (y[4] * R61); 
    grad_x[2] = gradop12(ielem, 0, 2) = (y[3] *  t3) + (y[0] * R42) - (y[1] *  t2) - (y[6] *  t4) + (y[7] * R74) - (y[5] * R72); 
    grad_x[3] = gradop12(ielem, 0, 3) = (y[0] *  t5) - (y[1] * R31) - (y[2] *  t3) + (y[7] *  t6) + (y[4] * R81) - (y[6] * R83); 
    grad_x[4] = gradop12(ielem, 0, 4) = (y[5] *  t3) + (y[6] * R86) - (y[7] *  t2) - (y[0] *  t4) - (y[3] * R81) + (y[1] * R61); 
    grad_x[5] = gradop12(ielem, 0, 5) = (y[6] *  t5) - (y[4] *  t3) - (y[7] * R75) + (y[1] *  t6) - (y[0] * R52) + (y[2] * R72);
    grad_x[6] = gradop12(ielem, 0, 6) = (y[7] *  t1) - (y[5] *  t5) - (y[4] * R86) + (y[2] *  t4) - (y[1] * R63) + (y[3] * R83);
    grad_x[7] = gradop12(ielem, 0, 7) = (y[4] *  t2) - (y[6] *  t1) + (y[5] * R75) - (y[3] *  t6) - (y[2] * R74) + (y[0] * R54); 


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

      grad_z[0] = gradop12(ielem, 2, 0) = (x[1] *  t1) - (x[2] * R42) - (x[3] *  t5)  + (x[4] *  t4) + (x[5] * R52) - (x[7] * R54); 
    grad_z[1] = gradop12(ielem, 2, 1) = (x[2] *  t2) + (x[3] * R31) - (x[0] *  t1)  - (x[5] *  t6) + (x[6] * R63) - (x[4] * R61); 
    grad_z[2] = gradop12(ielem, 2, 2) = (x[3] *  t3) + (x[0] * R42) - (x[1] *  t2)  - (x[6] *  t4) + (x[7] * R74) - (x[5] * R72); 
    grad_z[3] = gradop12(ielem, 2, 3) = (x[0] *  t5) - (x[1] * R31) - (x[2] *  t3)  + (x[7] *  t6) + (x[4] * R81) - (x[6] * R83); 
    grad_z[4] = gradop12(ielem, 2, 4) = (x[5] *  t3) + (x[6] * R86) - (x[7] *  t2)  - (x[0] *  t4) - (x[3] * R81) + (x[1] * R61); 
    grad_z[5] = gradop12(ielem, 2, 5) = (x[6] *  t5) - (x[4] *  t3)  - (x[7] * R75) + (x[1] *  t6) - (x[0] * R52) + (x[2] * R72);
    grad_z[6] = gradop12(ielem, 2, 6) = (x[7] *  t1) - (x[5] *  t5)  - (x[4] * R86) + (x[2] *  t4) - (x[1] * R63) + (x[3] * R83);
    grad_z[7] = gradop12(ielem, 2, 7) = (x[4] *  t2) - (x[6] *  t1)  + (x[5] * R75) - (x[3] *  t6) - (x[2] * R74) + (x[0] * R54); 

  }

  KOKKOS_MACRO_DEVICE_FUNCTION
    void v_grad(  int ielem, 
          Scalar * vx,     Scalar * vy,     Scalar * vz, 
          Scalar * grad_x,   Scalar * grad_y,   Scalar * grad_z, 
          Scalar inv_vol) const {

    //
    //  Calculate Velocity Gradients
    //

    vel_grad(ielem, K_F_XX) = inv_vol * (  vx[0] * grad_x[0] + \
                        vx[1] * grad_x[1] + \
                        vx[2] * grad_x[2] + \
                        vx[3] * grad_x[3] + \
                        vx[4] * grad_x[4] + \
                        vx[5] * grad_x[5] + \
                        vx[6] * grad_x[6] + \
                        vx[7] * grad_x[7] );

    vel_grad(ielem, K_F_YX) = inv_vol * (  vy[0] * grad_x[0] + \
                        vy[1] * grad_x[1] + \
                        vy[2] * grad_x[2] + \
                        vy[3] * grad_x[3] + \
                        vy[4] * grad_x[4] + \
                        vy[5] * grad_x[5] + \
                        vy[6] * grad_x[6] + \
                        vy[7] * grad_x[7] );

    vel_grad(ielem, K_F_ZX) = inv_vol * (  vz[0] * grad_x[0] + \
                        vz[1] * grad_x[1] + \
                        vz[2] * grad_x[2] + \
                        vz[3] * grad_x[3] + \
                        vz[4] * grad_x[4] + \
                        vz[5] * grad_x[5] + \
                        vz[6] * grad_x[6] + \
                        vz[7] * grad_x[7] );



    vel_grad(ielem, K_F_XY) = inv_vol * (  vx[0] * grad_y[0] + \
                        vx[1] * grad_y[1] + \
                        vx[2] * grad_y[2] + \
                        vx[3] * grad_y[3] + \
                        vx[4] * grad_y[4] + \
                        vx[5] * grad_y[5] + \
                        vx[6] * grad_y[6] + \
                        vx[7] * grad_y[7] );

    vel_grad(ielem, K_F_YY) = inv_vol * (  vy[0] * grad_y[0] + \
                        vy[1] * grad_y[1] + \
                        vy[2] * grad_y[2] + \
                        vy[3] * grad_y[3] + \
                        vy[4] * grad_y[4] + \
                        vy[5] * grad_y[5] + \
                        vy[6] * grad_y[6] + \
                        vy[7] * grad_y[7] );

    vel_grad(ielem, K_F_ZY) = inv_vol * (  vz[0] * grad_y[0] + \
                        vz[1] * grad_y[1] + \
                        vz[2] * grad_y[2] + \
                        vz[3] * grad_y[3] + \
                        vz[4] * grad_y[4] + \
                        vz[5] * grad_y[5] + \
                        vz[6] * grad_y[6] + \
                        vz[7] * grad_y[7] );



    vel_grad(ielem, K_F_XZ) = inv_vol * (  vx[0] * grad_z[0] + \
                        vx[1] * grad_z[1] + \
                        vx[2] * grad_z[2] + \
                        vx[3] * grad_z[3] + \
                        vx[4] * grad_z[4] + \
                        vx[5] * grad_z[5] + \
                        vx[6] * grad_z[6] + \
                        vx[7] * grad_z[7] );

    vel_grad(ielem, K_F_YZ) = inv_vol * (  vy[0] * grad_z[0] + \
                        vy[1] * grad_z[1] + \
                        vy[2] * grad_z[2] + \
                        vy[3] * grad_z[3] + \
                        vy[4] * grad_z[4] + \
                        vy[5] * grad_z[5] + \
                        vy[6] * grad_z[6] + \
                        vy[7] * grad_z[7] );

    vel_grad(ielem, K_F_ZZ) = inv_vol * (  vz[0] * grad_z[0] + \
                        vz[1] * grad_z[1] + \
                        vz[2] * grad_z[2] + \
                        vz[3] * grad_z[3] + \
                        vz[4] * grad_z[4] + \
                        vz[5] * grad_z[5] + \
                        vz[6] * grad_z[6] + \
                        vz[7] * grad_z[7] );

  }

  KOKKOS_MACRO_DEVICE_FUNCTION
    void comp_hgop(    int ielem, 
            Scalar * x,     Scalar * y,     Scalar * z, 
            Scalar * grad_x,   Scalar * grad_y,   Scalar * grad_z, 
            Scalar inv_vol) const {

    // KHP: Alternatively, we could have
      // hx0,hx1,hx2,hx3,...,hz0,hz1,hz2,hz3
      Scalar hgconst12th[12];

      Scalar q0 = x[0] - x[1];
      Scalar q1 = x[2] - x[3];
      Scalar q2 = x[4] - x[5];
      Scalar q3 = x[6] - x[7];

      hgconst12th[0] = ( (x[0]+x[1]) - (x[2]+x[3]) - (x[4]+x[5]) + (x[6]+x[7]) ) * inv_vol;
      hgconst12th[1] = (  q0 - q1 - q2 + q3 ) * inv_vol;
      hgconst12th[2] = (  q0 + q1 + q2 + q3 ) * inv_vol;
      hgconst12th[3] = ( -q0 - q1 + q2 + q3 ) * inv_vol;

      q0 = (y[0] - y[1]);
      q1 = (y[2] - y[3]);
      q2 = (y[4] - y[5]);
      q3 = (y[6] - y[7]);

      hgconst12th[4] = ( (y[0]+y[1]) - (y[2]+y[3]) - (y[4]+y[5]) + (y[6]+y[7]) ) * inv_vol;
      hgconst12th[5] = (  q0 - q1 - q2 + q3 ) * inv_vol;
      hgconst12th[6] = (  q0 + q1 + q2 + q3 ) * inv_vol;
      hgconst12th[7] = ( -q0 - q1 + q2 + q3 ) * inv_vol;
  
      q0 = (z[0] - z[1]);
      q1 = (z[2] - z[3]);
      q2 = (z[4] - z[5]);
      q3 = (z[6] - z[7]);
  
      hgconst12th[8]  = ( (z[0]+z[1]) - (z[2]+z[3]) - (z[4]+z[5]) + (z[6]+z[7]) ) * inv_vol;
      hgconst12th[9]  = (  q0 - q1 - q2 + q3 ) * inv_vol;
      hgconst12th[10] = (  q0 + q1 + q2 + q3 ) * inv_vol;
      hgconst12th[11] = ( -q0 - q1 + q2 + q3 ) * inv_vol;

    /*

    Alternatively:

    const Scalar hgop_arr[32] = {
     1.0, 1.0,-1.0,-1.0,-1.0,-1.0, 1.0,  1.0,
     1.0,-1.0,-1.0, 1.0,-1.0, 1.0, 1.0, -1.0,
     1.0,-1.0, 1.0,-1.0, 1.0,-1.0, 1.0, -1.0,
    -1.0, 1.0,-1.0, 1.0, 1.0,-1.0, 1.0, -1.0};

    for(int i = 0; i < 8; i++){
      hgop(ielem, i   , 1) = hgop_arr[i   ] - (hgconst12th[0] * x[i][1]  + hgconst12th[4] * y[i][1] + hgconst12th[8 ] * z[i][1]);
      hgop(ielem, i+8 , 1) = hgop_arr[i+8 ] - (hgconst12th[1] * x[i][1]  + hgconst12th[5] * y[i][1] + hgconst12th[9 ] * z[i][1]);
      hgop(ielem, i+16, 1) = hgop_arr[i+16] - (hgconst12th[2] * x[i][1]  + hgconst12th[6] * y[i][1] + hgconst12th[10] * z[i][1]);
      hgop(ielem, i+24, 1) = hgop_arr[i+24] - (hgconst12th[3] * x[i][1]  + hgconst12th[7] * y[i][1] + hgconst12th[11] * z[i][1]);
    }

    */

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

  KOKKOS_MACRO_DEVICE_FUNCTION
    void operator()( int ielem )const {

    Scalar      x[8],      y[8],      z[8];
    Scalar     vx[8],     vy[8],     vz[8];
    Scalar grad_x[8], grad_y[8], grad_z[8];

    comp_grad(ielem, x, y, z, vx, vy, vz, grad_x, grad_y, grad_z);

    //
    //  Calculate hexahedral volume from x position and gradient information
    //
  
    Scalar inv_vol;

    inv_vol = mid_vol(ielem) = ONE12TH *(  x[0] * grad_x[0] +\
                        x[1] * grad_x[1] +\
                        x[2] * grad_x[2] +\
                        x[3] * grad_x[3] +\
                        x[4] * grad_x[4] +\
                        x[5] * grad_x[5] +\
                        x[6] * grad_x[6] +\
                        x[7] * grad_x[7] );

    inv_vol = 1.0 / inv_vol;

    v_grad(ielem, vx, vy, vz, grad_x, grad_y, grad_z, inv_vol);

    inv_vol *= ONE12TH;

      comp_hgop(ielem, x, y, z, grad_x, grad_y, grad_z, inv_vol);

  }

};

#endif
