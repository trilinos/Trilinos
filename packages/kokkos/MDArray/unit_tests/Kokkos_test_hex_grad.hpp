/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

//each element is a unit cube
template< typename Scalar , class DeviceMap >
struct HexSimpleFill {
  typedef          DeviceMap              device_map_type ;
  typedef typename DeviceMap::device_type device_type ;
  typedef typename DeviceMap::size_type   size_type ;

  Kokkos::MDArrayView<Scalar,DeviceMap> coords ;

  HexSimpleFill(
    const Kokkos::MDArrayView<Scalar,DeviceMap> & arg_coords )
    : coords( arg_coords ) {}

  KOKKOS_DEVICE_AND_HOST_FUNCTION
  size_type work_count() const { return coords.dimension(2); }

  KOKKOS_DEVICE_FUNCTION
  void operator()( size_type ielem ) const
  {
    coords(0,0,ielem) = 0.;
    coords(1,0,ielem) = 0.;
    coords(2,0,ielem) = 0.;

    coords(0,1,ielem) = 1.;
    coords(1,1,ielem) = 0.;
    coords(2,1,ielem) = 0.;

    coords(0,2,ielem) = 1.;
    coords(1,2,ielem) = 1.;
    coords(2,2,ielem) = 0.;

    coords(0,3,ielem) = 0.;
    coords(1,3,ielem) = 1.;
    coords(2,3,ielem) = 0.;


    coords(0,4,ielem) = 0.;
    coords(1,4,ielem) = 0.;
    coords(2,4,ielem) = 1.;

    coords(0,5,ielem) = 1.;
    coords(1,5,ielem) = 0.;
    coords(2,5,ielem) = 1.;

    coords(0,6,ielem) = 1.;
    coords(1,6,ielem) = 1.;
    coords(2,6,ielem) = 1.;

    coords(0,7,ielem) = 0.;
    coords(1,7,ielem) = 1.;
    coords(2,7,ielem) = 1.;
  }
};

template< typename Scalar , class DeviceMap >
struct HexGrad {
  typedef          DeviceMap              device_map_type ;
  typedef typename DeviceMap::device_type device_type ;
  typedef typename DeviceMap::size_type   size_type ;

  Kokkos::MDArrayView<Scalar,DeviceMap> coords ; // (Space,Node,ParallelWork)
  Kokkos::MDArrayView<Scalar,DeviceMap> grad_op ;// (Space,Node,ParallelWork)

  HexGrad(
    const Kokkos::MDArrayView<Scalar,DeviceMap> & arg_coords,
    const Kokkos::MDArrayView<Scalar,DeviceMap> & arg_grad_op )
    : coords( arg_coords ) , grad_op( arg_grad_op ) {}

  KOKKOS_DEVICE_AND_HOST_FUNCTION
  size_type work_count() const { return coords.dimension(2); }

  KOKKOS_DEVICE_FUNCTION
  void operator()( size_type ielem ) const
  {
    // Repeated re-use of nodal coordinates,
    // copy them into local storage.

    Scalar x[8],y[8],z[8];

    x[0] = coords(0,0,ielem);
    x[1] = coords(0,1,ielem);
    x[2] = coords(0,2,ielem);
    x[3] = coords(0,3,ielem);
    x[4] = coords(0,4,ielem);
    x[5] = coords(0,5,ielem);
    x[6] = coords(0,6,ielem);
    x[7] = coords(0,7,ielem);

    y[0] = coords(1,0,ielem);
    y[1] = coords(1,1,ielem);
    y[2] = coords(1,2,ielem);
    y[3] = coords(1,3,ielem);
    y[4] = coords(1,4,ielem);
    y[5] = coords(1,5,ielem);
    y[6] = coords(1,6,ielem);
    y[7] = coords(1,7,ielem);


    z[0] = coords(2,0,ielem);
    z[1] = coords(2,1,ielem);
    z[2] = coords(2,2,ielem);
    z[3] = coords(2,3,ielem);
    z[4] = coords(2,4,ielem);
    z[5] = coords(2,5,ielem);
    z[6] = coords(2,6,ielem);
    z[7] = coords(2,7,ielem);

    // z difference vectors
    Scalar R42=(z[3] - z[1]);
    Scalar R52=(z[4] - z[1]);
    Scalar R54=(z[4] - z[3]);

    Scalar R63=(z[5] - z[2]);
    Scalar R83=(z[7] - z[2]);
    Scalar R86=(z[7] - z[5]);

    Scalar R31=(z[2] - z[0]);
    Scalar R61=(z[5] - z[0]);
    Scalar R74=(z[6] - z[3]);

    Scalar R72=(z[6] - z[1]);
    Scalar R75=(z[6] - z[4]);
    Scalar R81=(z[7] - z[0]);

    Scalar t1=(R63 + R54);
    Scalar t2=(R61 + R74);
    Scalar t3=(R72 + R81);

    Scalar t4 =(R86 + R42);
    Scalar t5 =(R83 + R52);
    Scalar t6 =(R75 + R31);

    grad_op(0,0,ielem) = (y[1] *  t1) - (y[2] * R42) - (y[3] *  t5)  + (y[4] *  t4) + (y[5] * R52) - (y[7] * R54);
    grad_op(0,1,ielem) = (y[2] *  t2) + (y[3] * R31) - (y[0] *  t1)  - (y[5] *  t6) + (y[6] * R63) - (y[4] * R61);
    grad_op(0,2,ielem) = (y[3] *  t3) + (y[0] * R42) - (y[1] *  t2)  - (y[6] *  t4) + (y[7] * R74) - (y[5] * R72);
    grad_op(0,3,ielem) = (y[0] *  t5) - (y[1] * R31) - (y[2] *  t3)  + (y[7] *  t6) + (y[4] * R81) - (y[6] * R83);
    grad_op(0,4,ielem) = (y[5] *  t3) + (y[6] * R86) - (y[7] *  t2)  - (y[0] *  t4) - (y[3] * R81) + (y[1] * R61);
    grad_op(0,5,ielem) = (y[6] *  t5) - (y[4] *  t3)  - (y[7] * R75) + (y[1] *  t6) - (y[0] * R52) + (y[2] * R72);
    grad_op(0,6,ielem) = (y[7] *  t1) - (y[5] *  t5)  - (y[4] * R86) + (y[2] *  t4) - (y[1] * R63) + (y[3] * R83);
    grad_op(0,7,ielem) = (y[4] *  t2) - (y[6] *  t1)  + (y[5] * R75) - (y[3] *  t6) - (y[2] * R74) + (y[0] * R54);

    R42=(x[3] - x[1]);
    R52=(x[4] - x[1]);
    R54=(x[4] - x[3]);

    R63=(x[5] - x[2]);
    R83=(x[7] - x[2]);
    R86=(x[7] - x[5]);

    R31=(x[2] - x[0]);
    R61=(x[5] - x[0]);
    R74=(x[6] - x[3]);

    R72=(x[6] - x[1]);
    R75=(x[6] - x[4]);
    R81=(x[7] - x[0]);

    t1=(R63 + R54);
    t2=(R61 + R74);
    t3=(R72 + R81);

    t4 =(R86 + R42);
    t5 =(R83 + R52);
    t6 =(R75 + R31);

    grad_op(1,0,ielem ) = (z[1] *  t1) - (z[2] * R42) - (z[3] *  t5)  + (z[4] *  t4) + (z[5] * R52) - (z[7] * R54);
    grad_op(1,1,ielem ) = (z[2] *  t2) + (z[3] * R31) - (z[0] *  t1)  - (z[5] *  t6) + (z[6] * R63) - (z[4] * R61);
    grad_op(1,2,ielem) = (z[3] *  t3) + (z[0] * R42) - (z[1] *  t2)  - (z[6] *  t4) + (z[7] * R74) - (z[5] * R72);
    grad_op(1,3,ielem) = (z[0] *  t5) - (z[1] * R31) - (z[2] *  t3)  + (z[7] *  t6) + (z[4] * R81) - (z[6] * R83);
    grad_op(1,4,ielem) = (z[5] *  t3) + (z[6] * R86) - (z[7] *  t2)  - (z[0] *  t4) - (z[3] * R81) + (z[1] * R61);
    grad_op(1,5,ielem) = (z[6] *  t5) - (z[4] *  t3)  - (z[7] * R75) + (z[1] *  t6) - (z[0] * R52) + (z[2] * R72);
    grad_op(1,6,ielem) = (z[7] *  t1) - (z[5] *  t5)  - (z[4] * R86) + (z[2] *  t4) - (z[1] * R63) + (z[3] * R83);
    grad_op(1,7,ielem) = (z[4] *  t2) - (z[6] *  t1)  + (z[5] * R75) - (z[3] *  t6) - (z[2] * R74) + (z[0] * R54);

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

    grad_op(2,0,ielem) = (x[1] *  t1) - (x[2] * R42) - (x[3] *  t5)  + (x[4] *  t4) + (x[5] * R52) - (x[7] * R54);
    grad_op(2,1,ielem) = (x[2] *  t2) + (x[3] * R31) - (x[0] *  t1)  - (x[5] *  t6) + (x[6] * R63) - (x[4] * R61);
    grad_op(2,2,ielem) = (x[3] *  t3) + (x[0] * R42) - (x[1] *  t2)  - (x[6] *  t4) + (x[7] * R74) - (x[5] * R72);
    grad_op(2,3,ielem) = (x[0] *  t5) - (x[1] * R31) - (x[2] *  t3)  + (x[7] *  t6) + (x[4] * R81) - (x[6] * R83);
    grad_op(2,4,ielem) = (x[5] *  t3) + (x[6] * R86) - (x[7] *  t2)  - (x[0] *  t4) - (x[3] * R81) + (x[1] * R61);
    grad_op(2,5,ielem) = (x[6] *  t5) - (x[4] *  t3)  - (x[7] * R75) + (x[1] *  t6) - (x[0] * R52) + (x[2] * R72);
    grad_op(2,6,ielem) = (x[7] *  t1) - (x[5] *  t5)  - (x[4] * R86) + (x[2] *  t4) - (x[1] * R63) + (x[3] * R83);
    grad_op(2,7,ielem) = (x[4] *  t2) - (x[6] *  t1)  + (x[5] * R75) - (x[3] *  t6) - (x[2] * R74) + (x[0] * R54);
  }

};



template< typename Scalar , class DeviceMap >
void test_hex_grad()
{
  const int parallel_work_length = 1000000 ;

  DeviceMap map( parallel_work_length );
  DeviceMap map_2( parallel_work_length * 2 );

  Kokkos::MDArrayView< Scalar , DeviceMap > coord , grad ;

  // Create arrays mapped onto host device
  coord = map.template create_labeled_mdarray<Scalar>( 3 , 8 , "coord" );
  grad  = map.template create_labeled_mdarray<Scalar>( 3 , 8 , "grad" );

  // Potential alternative API for allocating array:
  //
  // coord.allocate( 3 , 8 , map , "coord" );
  // grad.allocate(  3 , 8 , map , "grad" );

  // Create additional views and then destroy them for testing
  {
    Kokkos::MDArrayView< Scalar , DeviceMap > tmp1 = coord ;
    Kokkos::MDArrayView< Scalar , DeviceMap > tmp2 = tmp1 ;
  }

  // Execute the parallel kernels on the arrays:

  Kokkos::parallel_for( HexSimpleFill<Scalar,DeviceMap>( coord ) );
  Kokkos::parallel_for( HexGrad<Scalar,DeviceMap>( coord , grad ) );
}



