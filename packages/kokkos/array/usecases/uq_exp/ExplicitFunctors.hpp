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

#ifndef KOKKOSARRAY_HEXEXPLICITFUNCTORS_HPP
#define KOKKOSARRAY_HEXEXPLICITFUNCTORS_HPP

#include <cstdlib>
#include <cmath>

#include <HexExplicitFunctions.hpp>
#include <FEMesh.hpp>

namespace Explicit {

//----------------------------------------------------------------------------

template< typename Scalar , class Device >
struct Fields {

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

  typedef Device  device_type ;

  typedef double ScalarPrecise ;
  typedef float  ScalarCompact ;

  typedef HybridFEM::FEMesh<double,ElemNodeCount,device_type>  FEMesh ;

  typedef typename FEMesh::node_coords_type      node_coords_type ;
  typedef typename FEMesh::elem_node_ids_type    elem_node_ids_type ;
  typedef typename FEMesh::node_elem_ids_type    node_elem_ids_type ;
  typedef typename KokkosArray::ParallelDataMap  parallel_data_map ;

  typedef KokkosArray::View< ScalarPrecise ** [ SpatialDim ] , device_type > spatial_precise_view ;
  typedef KokkosArray::View< ScalarCompact  ** [ SpatialDim ] , device_type > spatial_compact_view ;


  typedef KokkosArray::View< Scalar*                , device_type >  scalar_view ;
  typedef KokkosArray::View< Scalar**               , device_type >  array_view ;
  typedef KokkosArray::View< Scalar**[ K_F_SIZE ]   , device_type >  tensor_array_view ;
  typedef KokkosArray::View< Scalar**[ K_S_SIZE ]   , device_type >  sym_tensor_array_view ;

  typedef KokkosArray::View< ScalarCompact**[ SpatialDim ][ ElemNodeCount ] , device_type >
    elem_node_spatial_view ;

  typedef KokkosArray::View< Scalar ,   device_type > value_view ;
  typedef KokkosArray::View< Scalar* , device_type > property_view ;

  // Parameters:
  const unsigned num_nodes ;
  const unsigned num_nodes_owned ;
  const unsigned num_elements ;
  const unsigned uq_count ;

  const float  lin_bulk_visc;
  const float  quad_bulk_visc;
  const float  two_mu;
  const float  bulk_modulus;
  const float  density;

  // Views of mesh data:
  const elem_node_ids_type  elem_node_connectivity ;
  const node_elem_ids_type  node_elem_connectivity ;
  const node_coords_type    model_coords ;

  // Properties:
  const property_view  nodal_mass ;
  const property_view  elem_mass ;

  // Compute data:
        value_view             dt ;
        value_view             dt_new ;
        spatial_precise_view   displacement ;
        spatial_precise_view   displacement_new ;
        spatial_precise_view   velocity ;
        spatial_precise_view   velocity_new ;
        tensor_array_view      rotation ;
        tensor_array_view      rotation_new ;

  const spatial_compact_view   acceleration ;
  const spatial_compact_view   internal_force ;
  const array_view             internal_energy ;
  const sym_tensor_array_view  stress ;
  const elem_node_spatial_view    element_force ;
  const tensor_array_view      vel_grad ;
  const sym_tensor_array_view  stretch ;
  const sym_tensor_array_view  rot_stretch ;

  Fields(
      const FEMesh & mesh,
      const unsigned arg_uq_count ,
      const float    arg_lin_bulk_visc,
      const float    arg_quad_bulk_visc,
      const float    youngs_modulus,
      const float    poissons_ratio,
      const float    arg_density )
    : num_nodes(       mesh.parallel_data_map.count_owned +
                       mesh.parallel_data_map.count_receive )
    , num_nodes_owned( mesh.parallel_data_map.count_owned )
    , num_elements(    mesh.elem_node_ids.dimension_0() )
    , uq_count(        arg_uq_count )
    , lin_bulk_visc(   arg_lin_bulk_visc )
    , quad_bulk_visc(  arg_quad_bulk_visc )
    , two_mu(          youngs_modulus/(1.0+poissons_ratio) )
    , bulk_modulus(    youngs_modulus/(3*(1.0-2.0*poissons_ratio)) )
    , density(         arg_density )

    // mesh

    , elem_node_connectivity( mesh.elem_node_ids )
    , node_elem_connectivity( mesh.node_elem_ids )
    , model_coords(           mesh.node_coords )

    // Properties:
    , nodal_mass( "nodal_mass" , num_nodes_owned )
    , elem_mass(  "elem_mass" ,  num_elements )

    , dt(     "dt" )
    , dt_new( "dt" )
    , displacement(     "displacement" ,    num_nodes , uq_count )
    , displacement_new( "displacement" ,    num_nodes , uq_count )
    , velocity(         "velocity" ,        num_nodes , uq_count )
    , velocity_new(     "velocity" ,        num_nodes , uq_count )
    , rotation(         "rotation" ,        num_elements , uq_count )
    , rotation_new(     "rotation" ,        num_elements , uq_count )
    , acceleration(     "acceleration" ,    num_nodes_owned , uq_count )
    , internal_force(   "internal_force" ,  num_nodes_owned , uq_count )
    , internal_energy(  "internal_energy" , num_elements , uq_count )
    , stress(           "stress" ,          num_elements , uq_count )
    , element_force(    "element_force" ,   num_elements , uq_count )
    , vel_grad(         "vel_grad" ,        num_elements , uq_count )
    , stretch(          "stretch" ,         num_elements , uq_count )
    , rot_stretch(      "rot_stretch" ,     num_elements , uq_count )
    { }
};

template< class DataType , class LayoutType , class Device >
void swap( KokkosArray::View< DataType , LayoutType , Device > & x ,
           KokkosArray::View< DataType , LayoutType , Device > & y )
{
  const KokkosArray::View< DataType , LayoutType , Device > z = x ;
  x = y ;
  y = z ;
}

//----------------------------------------------------------------------------

template< class Fields > struct InitializeElement ;
template< class Fields > struct InitializeNode ;
template< class Fields > struct GradFunctor ;
template< class Fields > struct DecompRotateFunctor ;
template< class Fields > struct InternalForceFunctor ;
template< class Fields , class Boundary > struct NodalUpdateFunctor ;
template< class Fields > struct NodalBoundary ;
template< class Fields > struct PackState ;
template< class Fields > struct UnpackState ;


template< class Fields >
inline
void initialize_element( const Fields & arg_fields )
{ InitializeElement< Fields > op( arg_fields ); }

template< class Fields >
inline
void initialize_node( const Fields & arg_fields )
{ InitializeNode< Fields > op( arg_fields ); }

template< class Fields >
inline
void gradient( const Fields & arg_fields )
{ GradFunctor< Fields > op( arg_fields ); }

template< class Fields >
inline
void decomp_rotate( const Fields & arg_fields )
{ DecompRotateFunctor< Fields > op( arg_fields ); }

template< class Fields >
inline
void internal_force( const Fields & arg_fields , float user_dt )
{ InternalForceFunctor< Fields > op( arg_fields , user_dt ); }

template< class Fields >
inline
void nodal_update( const Fields & arg_fields , const float x_bc )
{
  typedef NodalBoundary< Fields > Boundary ;
  Boundary bc_op( arg_fields , x_bc );
  NodalUpdateFunctor< Fields , Boundary > op( arg_fields , bc_op );
}

template< typename FieldsScalar , typename ValueType , class Device >
inline
void pack_state( const Fields< FieldsScalar , Device >  & arg_fields ,
                 const KokkosArray::View< ValueType* , Device > & arg_buffer ,
                 const unsigned node_begin ,
                 const unsigned node_count )
{
  typedef Fields< FieldsScalar , Device > fields_type ;
  typedef PackState< fields_type > Pack ;
  if ( node_count ) {
    Pack( arg_buffer , arg_fields , node_begin , node_count );
  }
}

template< typename FieldsScalar , typename ValueType , class Device >
inline
void unpack_state( const Fields< FieldsScalar , Device >  & arg_fields ,
                   const KokkosArray::View< ValueType* , Device > & arg_buffer ,
                   const unsigned node_begin ,
                   const unsigned node_count )
{
  typedef Fields< FieldsScalar , Device > fields_type ;
  typedef UnpackState< fields_type > Unpack ;
  if ( node_count ) {
    Unpack( arg_buffer , arg_fields , node_begin , node_count );
  }
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class Fields >
struct InitializeElement
{
  typedef typename Fields::device_type     device_type ;

  const typename Fields::elem_node_ids_type     elem_node_connectivity ;
  const typename Fields::node_coords_type       model_coords ;

  const typename Fields::sym_tensor_array_view  stretch ;
  const typename Fields::tensor_array_view      rotation ;
  const typename Fields::tensor_array_view      rotation_new ;
  const typename Fields::property_view          elem_mass ;

  const float    density ;
  const unsigned uq_count ;
  const unsigned K_S_XX ;
  const unsigned K_S_YY ;
  const unsigned K_S_ZZ ;
  const unsigned K_F_XX ;
  const unsigned K_F_YY ;
  const unsigned K_F_ZZ ;


  InitializeElement( const Fields & mesh_fields )
    : elem_node_connectivity( mesh_fields.elem_node_connectivity )
    , model_coords(           mesh_fields.model_coords )
    , stretch(                mesh_fields.stretch )
    , rotation(               mesh_fields.rotation )
    , rotation_new(           mesh_fields.rotation_new )
    , elem_mass(              mesh_fields.elem_mass )
    , density(                mesh_fields.density )
    , uq_count(               mesh_fields.stretch.dimension_1() )
    , K_S_XX( Hex8Functions::K_S_XX )
    , K_S_YY( Hex8Functions::K_S_YY )
    , K_S_ZZ( Hex8Functions::K_S_ZZ )
    , K_F_XX( Hex8Functions::K_F_XX )
    , K_F_YY( Hex8Functions::K_F_YY )
    , K_F_ZZ( Hex8Functions::K_F_ZZ )
    {
      KokkosArray::parallel_for( mesh_fields.num_elements , *this );
    }

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( unsigned ielem )const
  {
    const float ONE12TH = 1.0 / 12.0 ;

    typename Fields::ScalarPrecise x[ Fields::ElemNodeCount ];
    typename Fields::ScalarPrecise y[ Fields::ElemNodeCount ];
    typename Fields::ScalarPrecise z[ Fields::ElemNodeCount ];
    typename Fields::ScalarCompact grad_x[ Fields::ElemNodeCount ];
    typename Fields::ScalarCompact grad_y[ Fields::ElemNodeCount ];
    typename Fields::ScalarCompact grad_z[ Fields::ElemNodeCount ];

    for ( unsigned i = 0 ; i < Fields::ElemNodeCount ; ++i ) {
      const unsigned n = elem_node_connectivity( ielem , i );

      x[i]  = model_coords( n , 0 );
      y[i]  = model_coords( n , 1 );
      z[i]  = model_coords( n , 2 );
    }

    Hex8Functions::grad( x, y, z, grad_x, grad_y, grad_z);

    elem_mass(ielem) = ONE12TH * density * Hex8Functions::dot8( x , grad_x );

    for ( unsigned jp = 0 ; jp < uq_count ; ++jp ) {

      stretch(ielem,jp,K_S_XX) = 1 ;
      stretch(ielem,jp,K_S_YY) = 1 ;
      stretch(ielem,jp,K_S_ZZ) = 1 ;

      rotation(ielem,jp,K_F_XX) = 1 ;
      rotation(ielem,jp,K_F_YY) = 1 ;
      rotation(ielem,jp,K_F_ZZ) = 1 ;

      rotation_new(ielem,jp,K_F_XX) = 1 ;
      rotation_new(ielem,jp,K_F_YY) = 1 ;
      rotation_new(ielem,jp,K_F_ZZ) = 1 ;
    }
  }
};


template<typename Scalar, class DeviceType >
struct InitializeNode< Explicit::Fields< Scalar, DeviceType > >
{
  typedef DeviceType     device_type ;

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

  KOKKOSARRAY_INLINE_FUNCTION
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

template< typename Scalar , class DeviceType >
struct GradFunctor< Explicit::Fields< Scalar , DeviceType > >
{
  typedef DeviceType device_type ;

  typedef Explicit::Fields< Scalar , device_type >  Fields ;

  static const unsigned ElemNodeCount = Fields::ElemNodeCount ;
  static const unsigned K_F_SIZE      = Fields::K_F_SIZE ;
  static const unsigned K_S_SIZE      = Fields::K_S_SIZE ;

  // Global arrays used by this functor.

  const typename Fields::elem_node_ids_type  elem_node_connectivity ;
  const typename Fields::node_coords_type    model_coords ;

  const typename Fields::value_view          dt ;
  const typename Fields::spatial_precise_view  displacement ; 
  const typename Fields::spatial_precise_view  velocity ; 
  const typename Fields::tensor_array_view     vel_grad ;
  const unsigned                               uq_count ;

  // Constructor on the Host to populate this device functor.
  // All array view copies are shallow.
  GradFunctor( const Fields & fields )
    : elem_node_connectivity( fields.elem_node_connectivity )
    , model_coords(  fields.model_coords )
    , dt(            fields.dt )
    , displacement(  fields.displacement )
    , velocity(      fields.velocity )
    , vel_grad(      fields.vel_grad )
    , uq_count(      fields.displacement.dimension_1() )
    {
      KokkosArray::parallel_for( fields.num_elements , *this );
    }

  //--------------------------------------------------------------------------
  // Functor operator() which calls the three member functions.

  KOKKOSARRAY_INLINE_FUNCTION
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

      Hex8Functions::grad( x, y, z, grad_x, grad_y, grad_z);

      //  Calculate hexahedral volume from x model_coords and gradient information

      const Scalar inv_vol = 1.0 / Hex8Functions::dot8( x , grad_x );

      vel_grad( ielem, jp, K_F_XX ) = inv_vol * Hex8Functions::dot8( vx , grad_x );
      vel_grad( ielem, jp, K_F_YX ) = inv_vol * Hex8Functions::dot8( vy , grad_x );
      vel_grad( ielem, jp, K_F_ZX ) = inv_vol * Hex8Functions::dot8( vz , grad_x );

      vel_grad( ielem, jp, K_F_XY ) = inv_vol * Hex8Functions::dot8( vx , grad_y );
      vel_grad( ielem, jp, K_F_YY ) = inv_vol * Hex8Functions::dot8( vy , grad_y );
      vel_grad( ielem, jp, K_F_ZY ) = inv_vol * Hex8Functions::dot8( vz , grad_y );

      vel_grad( ielem, jp, K_F_XZ ) = inv_vol * Hex8Functions::dot8( vx , grad_z );
      vel_grad( ielem, jp, K_F_YZ ) = inv_vol * Hex8Functions::dot8( vy , grad_z );
      vel_grad( ielem, jp, K_F_ZZ ) = inv_vol * Hex8Functions::dot8( vz , grad_z );
    }
  }
};

//----------------------------------------------------------------------------

template< typename Scalar , class DeviceType >
struct DecompRotateFunctor< Explicit::Fields< Scalar , DeviceType > >
{
  typedef DeviceType device_type ;

  typedef Explicit::Fields< Scalar , device_type >  Fields ;

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
    , uq_count(      mesh_fields.vel_grad.dimension_1() )
    {
      KokkosArray::parallel_for( mesh_fields.num_elements , *this );
    }

  KOKKOSARRAY_INLINE_FUNCTION
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
      Hex8Functions::polar_decomp( step , v_gr, str, str_ten, rot );

      for ( unsigned i = 0; i < K_S_SIZE ; ++i ) stretch(      ielem, jp, i ) = str[i] ;
      for ( unsigned i = 0; i < K_F_SIZE ; ++i ) rotation_new( ielem, jp, i ) = rot[i] ;

      Hex8Functions::rotate_tensor( str_ten , rot , rot_str );

      for ( unsigned i = 0 ; i < K_S_SIZE ; ++i ) rot_stretch( ielem , jp , i ) = rot_str[i] ;
    }
  }
};

//----------------------------------------------------------------------------

template< typename Scalar , class DeviceType >
struct InternalForceFunctor< Explicit::Fields< Scalar , DeviceType > >
{
  typedef DeviceType device_type ;

  typedef Explicit::Fields< Scalar , device_type >  Fields ;

  static const unsigned ElemNodeCount = Fields::ElemNodeCount ;
  static const unsigned SpatialDim    = Fields::SpatialDim ;
  static const unsigned K_F_SIZE = Fields::K_F_SIZE ;
  static const unsigned K_S_SIZE = Fields::K_S_SIZE ;

  //--------------------------------------------------------------------------
  // Desired next time step reduction:

  typedef Scalar value_type;

  KOKKOSARRAY_INLINE_FUNCTION
  static void init(value_type & next_time_step ) {
     next_time_step  = 1.0e32;
  }

  KOKKOSARRAY_INLINE_FUNCTION
  static void join( volatile value_type & next_time_step ,
                    const volatile value_type & source )
  {
    next_time_step = next_time_step < source ? next_time_step : source;
  }

  struct SetNextTimeStep {
    typedef DeviceType  device_type ;
    typedef Scalar      value_type;

    const typename Fields::value_view dt ;

    SetNextTimeStep( const typename Fields::value_view & arg_dt )
      : dt( arg_dt ) {}

    KOKKOSARRAY_INLINE_FUNCTION
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
  const typename Fields::spatial_precise_view   displacement ;
  const typename Fields::spatial_precise_view   velocity ;
  const typename Fields::property_view          elem_mass ;
  const typename Fields::array_view             internal_energy ;
  const typename Fields::sym_tensor_array_view  stress ;
  const typename Fields::elem_node_spatial_view    element_force ;
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
    , uq_count(               mesh_fields.displacement.dimension_1() )
  {
    SetNextTimeStep op_dt( mesh_fields.dt_new );

    KokkosArray::parallel_reduce( mesh_fields.num_elements, *this, op_dt );
  }

  //--------------------------------------------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
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

      Hex8Functions::grad( x , y , z , grad_x , grad_y , grad_z );

      const Scalar mid_vol = Hex8Functions::dot8( x , grad_x );

      const Scalar aspect = 6.0 * mid_vol /
                            ( Hex8Functions::dot8( grad_x , grad_x ) +
                              Hex8Functions::dot8( grad_y , grad_y ) +
                              Hex8Functions::dot8( grad_z , grad_z ) );

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

      Hex8Functions::update_stress( step , two_mu , bulk_modulus , rot_str , stress_work );

      for ( unsigned i = 0 ; i < K_S_SIZE ; ++i ) { stress( ielem , kp , i ) = stress_work[i]; }

      Hex8Functions::rotate_tensor_backward( stress_work , rot , total_stress12th );

      total_stress12th[0] = ONE12TH*( total_stress12th[ 0 ] + bulkq );
      total_stress12th[1] = ONE12TH*( total_stress12th[ 1 ] + bulkq );
      total_stress12th[2] = ONE12TH*( total_stress12th[ 2 ] + bulkq );
      total_stress12th[3] = ONE12TH*( total_stress12th[ 3 ] );
      total_stress12th[4] = ONE12TH*( total_stress12th[ 4 ] );
      total_stress12th[5] = ONE12TH*( total_stress12th[ 5 ] );

      Scalar energy ;

      Hex8Functions::comp_force( vx, vy, vz, grad_x, grad_y, grad_z, total_stress12th, force, energy );

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

template< typename Scalar , class DeviceType >
struct NodalBoundary< Explicit::Fields< Scalar , DeviceType > >
{
  typedef DeviceType device_type ;
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

  KOKKOSARRAY_INLINE_FUNCTION
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

template< typename Scalar , class Boundary , class DeviceType >
struct NodalUpdateFunctor< Fields< Scalar , DeviceType > , Boundary >
{
  typedef DeviceType     device_type ;
  typedef typename device_type::size_type  size_type;

  typedef Explicit::Fields< Scalar , device_type >  Fields ;

  const Boundary                              boundary ;
  const typename Fields::node_elem_ids_type   node_elem_connectivity ;
  const typename Fields::node_coords_type     model_coords ;
  const typename Fields::property_view        nodal_mass ;

  const typename Fields::value_view            dt ;
  const typename Fields::value_view            dt_new ;
  const typename Fields::spatial_precise_view  displacement ;
  const typename Fields::spatial_precise_view  displacement_new ;
  const typename Fields::spatial_precise_view  velocity ;
  const typename Fields::spatial_precise_view  velocity_new ;
  const typename Fields::spatial_compact_view  acceleration ;
  const typename Fields::spatial_compact_view  internal_force ;
  const typename Fields::elem_node_spatial_view   element_force ;
  const unsigned                               uq_count ;

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
   , uq_count(          mesh_fields.displacement.dimension_1() )
   {
        //std::cout << "finish_step dt: " << dt << std::endl;
        //std::cout << "finish_step prev_dt: " << prev_dt << std::endl;

      // Only update the owned nodes:

      KokkosArray::parallel_for( mesh_fields.num_nodes_owned , *this );
   }

  KOKKOSARRAY_INLINE_FUNCTION
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

template< typename Scalar , class DeviceType >
struct PackState< Explicit::Fields< Scalar , DeviceType > >
{
  typedef DeviceType     device_type ;
  typedef typename device_type::size_type  size_type ;

  typedef Explicit::Fields< Scalar , device_type >  Fields ;

  typedef typename Fields::spatial_precise_view::scalar_type    value_type ;
  typedef KokkosArray::View< value_type*, device_type >  buffer_type ;

  static const unsigned value_count = 6 ;

  const typename Fields::spatial_precise_view  displacement ;
  const typename Fields::spatial_precise_view  velocity ;
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
    , uq_count(     mesh_fields.displacement.dimension_1() )
    {
      KokkosArray::parallel_for( arg_count , *this );
    }

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const unsigned i ) const
  {
    const unsigned inode  = inode_base + i ;
    unsigned j = i * value_count * uq_count ;

    for ( unsigned kq = 0 ; kq < uq_count ; ++kq ) {

      output[j++] = displacement( inode , kq , 0 );
      output[j++] = displacement( inode , kq , 1 );
      output[j++] = displacement( inode , kq , 2 );
      output[j++] = velocity( inode , kq , 0 );
      output[j++] = velocity( inode , kq , 1 );
      output[j++] = velocity( inode , kq , 2 );
    }
  }
};

template< typename Scalar , class DeviceType >
struct UnpackState< Explicit::Fields< Scalar , DeviceType > >
{
  typedef DeviceType     device_type ;
  typedef typename device_type::size_type  size_type ;

  typedef Explicit::Fields< Scalar , device_type >  Fields ;

  typedef typename Fields::spatial_precise_view::scalar_type     value_type ;
  typedef KokkosArray::View< value_type* , device_type >  buffer_type ;

  static const unsigned value_count = 6 ;

  const typename Fields::spatial_precise_view  displacement ;
  const typename Fields::spatial_precise_view  velocity ;
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
    , uq_count(     mesh_fields.displacement.dimension_1() )
    {
      KokkosArray::parallel_for( arg_count , *this );
    }

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const unsigned i ) const
  {
    const unsigned inode = inode_base + i ;
    unsigned j = i * value_count * uq_count ;

    for ( unsigned kq = 0 ; kq < uq_count ; ++kq ) {

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


#endif /* #ifndef KOKKOSARRAY_HEXEXPLICITFUNCTORS_HPP */

