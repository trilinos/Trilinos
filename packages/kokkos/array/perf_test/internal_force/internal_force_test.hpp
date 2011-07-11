#ifndef INTERNALFORCE
#define INTERNALFORCE

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include "grad_hgop.hpp"
#include "decomp_rotate.hpp"
#include "divergence.hpp"

template<typename Scalar, class DeviceType >
struct InternalForceInitialize ;

//----------------------------------------------------------------------------
//
template<typename Scalar>
struct InternalForceInitialize<Scalar, KOKKOS_MACRO_DEVICE> {

  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef typename Kokkos::MDArrayView<Scalar,device_type> array_type ;

  array_type position;    
  array_type velocity;
  array_type rotation;    //rotation old and new

  InternalForceInitialize( const array_type & pos ,
                           const array_type & vel ,
                           const array_type & rot )
    : position( pos ), velocity( vel ), rotation( rot ) {}

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( int ielem ) const {

  //  unit cube in first octant, first node at origin

  //  X coords
    position(ielem, 0, 0) = 0.0;
    position(ielem, 0, 1) = 1.0;
    position(ielem, 0, 2) = 1.0;
    position(ielem, 0, 3) = 0.0;
    position(ielem, 0, 4) = 0.0;
    position(ielem, 0, 5) = 1.0;
    position(ielem, 0, 6) = 1.0;
    position(ielem, 0, 7) = 0.0;


  //  Y coords
    position(ielem, 1, 0) = 0.0;
    position(ielem, 1, 1) = 0.0;
    position(ielem, 1, 2) = 1.0;
    position(ielem, 1, 3) = 1.0;
    position(ielem, 1, 4) = 0.0;
    position(ielem, 1, 5) = 0.0;
    position(ielem, 1, 6) = 1.0;
    position(ielem, 1, 7) = 1.0;


  //  Z coords
    position(ielem, 2, 0) = 0.0;
    position(ielem, 2, 1) = 0.0;
    position(ielem, 2, 2) = 0.0;
    position(ielem, 2, 3) = 0.0;
    position(ielem, 2, 4) = 1.0;
    position(ielem, 2, 5) = 1.0;
    position(ielem, 2, 6) = 1.0;
    position(ielem, 2, 7) = 1.0;


  //  first node with velocity in the direction of
  //  the cube's centroid, all other nodes static

  
  //  X vel
    velocity(ielem, 0, 0) = 0.25;
    velocity(ielem, 0, 1) = 0.0;
    velocity(ielem, 0, 2) = 0.0;
    velocity(ielem, 0, 3) = 0.0;
    velocity(ielem, 0, 4) = 0.0;
    velocity(ielem, 0, 5) = 0.0;
    velocity(ielem, 0, 6) = 0.0;
    velocity(ielem, 0, 7) = 0.0;


  //  Y vel
    velocity(ielem, 1, 0) = 0.25;
    velocity(ielem, 1, 1) = 0.0;
    velocity(ielem, 1, 2) = 0.0;
    velocity(ielem, 1, 3) = 0.0;
    velocity(ielem, 1, 4) = 0.0;
    velocity(ielem, 1, 5) = 0.0;
    velocity(ielem, 1, 6) = 0.0;
    velocity(ielem, 1, 7) = 0.0;

  //  Z vel
    velocity(ielem, 2, 0) = 0.25;
    velocity(ielem, 2, 1) = 0.0;
    velocity(ielem, 2, 2) = 0.0;
    velocity(ielem, 2, 3) = 0.0;
    velocity(ielem, 2, 4) = 0.0;
    velocity(ielem, 2, 5) = 0.0;
    velocity(ielem, 2, 6) = 0.0;
    velocity(ielem, 2, 7) = 0.0;


  //  No rotations
    rotation(ielem, 0, 0) = 0.0;
    rotation(ielem, 1, 0) = 0.0;
    rotation(ielem, 2, 0) = 0.0;
    rotation(ielem, 3, 0) = 0.0;
    rotation(ielem, 4, 0) = 0.0;
    rotation(ielem, 5, 0) = 0.0;
    rotation(ielem, 6, 0) = 0.0;
    rotation(ielem, 7, 0) = 0.0;
    rotation(ielem, 8, 0) = 0.0;
  }
};

//----------------------------------------------------------------------------

template<typename Scalar>
double internal_force_test( const size_t nelem )
{
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef typename Kokkos::MDArrayView<Scalar,device_type> array_type ;

  array_type  position =      Kokkos::create_mdarray< array_type >(nelem, 3, 8);  
  array_type  velocity =      Kokkos::create_mdarray< array_type >(nelem, 3, 8);
  array_type  mid_pos =       Kokkos::create_mdarray< array_type >(nelem, 3, 8);
  array_type  hg_resist =     Kokkos::create_mdarray< array_type >(nelem, 12, 2); // old and new
  array_type  rotation =      Kokkos::create_mdarray< array_type >(nelem, 9, 2);  // rotation old and new
  array_type  gradop12 =      Kokkos::create_mdarray< array_type >(nelem, 3, 8);
  array_type  force_new =     Kokkos::create_mdarray< array_type >(nelem, 3, 8);
  array_type  hgop =          Kokkos::create_mdarray< array_type >(nelem, 32, 2); // hgop and mid_hgop
  array_type  vel_grad =      Kokkos::create_mdarray< array_type >(nelem, 9);
  array_type  stretch =       Kokkos::create_mdarray< array_type >(nelem, 6);
  array_type  s_temp =        Kokkos::create_mdarray< array_type >(nelem, 6);
  array_type  vorticity =     Kokkos::create_mdarray< array_type >(nelem, 3);
  array_type  rot_stret =     Kokkos::create_mdarray< array_type >(nelem, 6);
  array_type  stress_new =    Kokkos::create_mdarray< array_type >(nelem, 6);    
  array_type  rot_stress =    Kokkos::create_mdarray< array_type >(nelem, 6);
  array_type  mid_vol =       Kokkos::create_mdarray< array_type >(nelem);
  array_type  shrmod =        Kokkos::create_mdarray< array_type >(nelem);
  array_type  dilmod =        Kokkos::create_mdarray< array_type >(nelem);
  array_type  elem_mass =     Kokkos::create_mdarray< array_type >(nelem);
  array_type  elem_t_step =   Kokkos::create_mdarray< array_type >(nelem);
  array_type  hg_energy =     Kokkos::create_mdarray< array_type >(nelem);
  array_type  intern_energy = Kokkos::create_mdarray< array_type >(nelem);
  array_type  bulk_modulus =  Kokkos::create_mdarray< array_type >(nelem);
  array_type  two_mu =        Kokkos::create_mdarray< array_type >(nelem);

  // array_type strain_rate =   Kokkos::create_mdarray< array_type >(nelem);

  const Scalar  dt = 0.25;
  const Scalar  lin_bulk_visc = 0;
  const Scalar  quad_bulk_visc = 0;
  // const Scalar  stable_time_step = 0;
  const Scalar  hg_stiffness = 1.04;
  const Scalar  hg_viscosity = 0.93;
  // const Scalar  fac1_pre = dt * hg_stiffness * 0.0625;
  const bool    scaleHGRotation = false;

  Kokkos::parallel_for( nelem ,
    InternalForceInitialize< Scalar , device_type >( position , velocity , rotation ) );

  double total = 0.0;
  double time = 0.0;    

  Kokkos::parallel_for( nelem , grad_hgop<Scalar, device_type>
    ( position, 
      velocity, 
      mid_pos, 
      mid_vol, 
      gradop12, 
      vel_grad, 
      hgop, 
      dt )     , time );

  total += time;  

  Kokkos::parallel_for( nelem , decomp_rotate<Scalar, device_type>
    (  rotation, 
      vel_grad, 
      s_temp, 
      stretch, 
      vorticity, 
      rot_stret, 
      dt)     , time );

  total += time;

  Kokkos::parallel_for( nelem , divergence<Scalar, device_type>
    (  position, 
      velocity, 
      force_new,
      vorticity,
      rotation,
      stress_new,
      rot_stress,
      rot_stret,
      gradop12,
      elem_mass,
      dilmod,
      shrmod,
      elem_t_step,
      intern_energy,
      mid_vol,
      hgop,
      hg_resist,
      hg_energy,
      two_mu,
      bulk_modulus,
      hg_stiffness,
      hg_viscosity,
      lin_bulk_visc,
      quad_bulk_visc, 
      dt,
      scaleHGRotation)   , time );

  total += time;

  return total;
}

#endif
