/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

template< typename Scalar , class DeviceType >
struct HexSimpleFill ;

//each element is a unit cube
template< typename Scalar >
struct HexSimpleFill< Scalar , KOKKOS_MACRO_DEVICE >
{
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef device_type::size_type  size_type ;

  // 3D array : ( ParallelWork , Space , Node )
  typedef typename Kokkos::MDArrayView<Scalar,device_type> array_type ;

  array_type coords ;

  HexSimpleFill( const array_type & arg_coords )
    : coords( arg_coords ) {}

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type ielem ) const
  {
    coords(ielem,0,0) = 0.;
    coords(ielem,1,0) = 0.;
    coords(ielem,2,0) = 0.;

    coords(ielem,0,1) = 1.;
    coords(ielem,1,1) = 0.;
    coords(ielem,2,1) = 0.;

    coords(ielem,0,2) = 1.;
    coords(ielem,1,2) = 1.;
    coords(ielem,2,2) = 0.;

    coords(ielem,0,3) = 0.;
    coords(ielem,1,3) = 1.;
    coords(ielem,2,3) = 0.;


    coords(ielem,0,4) = 0.;
    coords(ielem,1,4) = 0.;
    coords(ielem,2,4) = 1.;

    coords(ielem,0,5) = 1.;
    coords(ielem,1,5) = 0.;
    coords(ielem,2,5) = 1.;

    coords(ielem,0,6) = 1.;
    coords(ielem,1,6) = 1.;
    coords(ielem,2,6) = 1.;

    coords(ielem,0,7) = 0.;
    coords(ielem,1,7) = 1.;
    coords(ielem,2,7) = 1.;
  }
};

//----------------------------------------------------------------------------
template< typename Scalar , class DeviceType >
struct HexGrad ;

template< typename Scalar >
struct HexGrad< Scalar , KOKKOS_MACRO_DEVICE >
{
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef device_type::size_type  size_type ;

  // 3D array : ( ParallelWork , Space , Node )
  typedef typename Kokkos::MDArrayView<Scalar,device_type> array_type ;

  array_type coords ;
  array_type grad_op ;

  HexGrad( const array_type & arg_coords,
           const array_type & arg_grad_op )
    : coords( arg_coords )
    , grad_op( arg_grad_op ) {}

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type ielem ) const
  {
    // Repeated re-use of nodal coordinates,
    // copy them into local storage.

    Scalar x[8],y[8],z[8];

    x[0] = coords(ielem,0,0);
    x[1] = coords(ielem,0,1);
    x[2] = coords(ielem,0,2);
    x[3] = coords(ielem,0,3);
    x[4] = coords(ielem,0,4);
    x[5] = coords(ielem,0,5);
    x[6] = coords(ielem,0,6);
    x[7] = coords(ielem,0,7);

    y[0] = coords(ielem,1,0);
    y[1] = coords(ielem,1,1);
    y[2] = coords(ielem,1,2);
    y[3] = coords(ielem,1,3);
    y[4] = coords(ielem,1,4);
    y[5] = coords(ielem,1,5);
    y[6] = coords(ielem,1,6);
    y[7] = coords(ielem,1,7);


    z[0] = coords(ielem,2,0);
    z[1] = coords(ielem,2,1);
    z[2] = coords(ielem,2,2);
    z[3] = coords(ielem,2,3);
    z[4] = coords(ielem,2,4);
    z[5] = coords(ielem,2,5);
    z[6] = coords(ielem,2,6);
    z[7] = coords(ielem,2,7);

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

    grad_op(ielem,0,0) = (y[1] *  t1) - (y[2] * R42) - (y[3] *  t5)  + (y[4] *  t4) + (y[5] * R52) - (y[7] * R54);
    grad_op(ielem,0,1) = (y[2] *  t2) + (y[3] * R31) - (y[0] *  t1)  - (y[5] *  t6) + (y[6] * R63) - (y[4] * R61);
    grad_op(ielem,0,2) = (y[3] *  t3) + (y[0] * R42) - (y[1] *  t2)  - (y[6] *  t4) + (y[7] * R74) - (y[5] * R72);
    grad_op(ielem,0,3) = (y[0] *  t5) - (y[1] * R31) - (y[2] *  t3)  + (y[7] *  t6) + (y[4] * R81) - (y[6] * R83);
    grad_op(ielem,0,4) = (y[5] *  t3) + (y[6] * R86) - (y[7] *  t2)  - (y[0] *  t4) - (y[3] * R81) + (y[1] * R61);
    grad_op(ielem,0,5) = (y[6] *  t5) - (y[4] *  t3)  - (y[7] * R75) + (y[1] *  t6) - (y[0] * R52) + (y[2] * R72);
    grad_op(ielem,0,6) = (y[7] *  t1) - (y[5] *  t5)  - (y[4] * R86) + (y[2] *  t4) - (y[1] * R63) + (y[3] * R83);
    grad_op(ielem,0,7) = (y[4] *  t2) - (y[6] *  t1)  + (y[5] * R75) - (y[3] *  t6) - (y[2] * R74) + (y[0] * R54);

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

    grad_op(ielem,1,0 ) = (z[1] *  t1) - (z[2] * R42) - (z[3] *  t5)  + (z[4] *  t4) + (z[5] * R52) - (z[7] * R54);
    grad_op(ielem,1,1 ) = (z[2] *  t2) + (z[3] * R31) - (z[0] *  t1)  - (z[5] *  t6) + (z[6] * R63) - (z[4] * R61);
    grad_op(ielem,1,2) = (z[3] *  t3) + (z[0] * R42) - (z[1] *  t2)  - (z[6] *  t4) + (z[7] * R74) - (z[5] * R72);
    grad_op(ielem,1,3) = (z[0] *  t5) - (z[1] * R31) - (z[2] *  t3)  + (z[7] *  t6) + (z[4] * R81) - (z[6] * R83);
    grad_op(ielem,1,4) = (z[5] *  t3) + (z[6] * R86) - (z[7] *  t2)  - (z[0] *  t4) - (z[3] * R81) + (z[1] * R61);
    grad_op(ielem,1,5) = (z[6] *  t5) - (z[4] *  t3)  - (z[7] * R75) + (z[1] *  t6) - (z[0] * R52) + (z[2] * R72);
    grad_op(ielem,1,6) = (z[7] *  t1) - (z[5] *  t5)  - (z[4] * R86) + (z[2] *  t4) - (z[1] * R63) + (z[3] * R83);
    grad_op(ielem,1,7) = (z[4] *  t2) - (z[6] *  t1)  + (z[5] * R75) - (z[3] *  t6) - (z[2] * R74) + (z[0] * R54);

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

    grad_op(ielem,2,0) = (x[1] *  t1) - (x[2] * R42) - (x[3] *  t5)  + (x[4] *  t4) + (x[5] * R52) - (x[7] * R54);
    grad_op(ielem,2,1) = (x[2] *  t2) + (x[3] * R31) - (x[0] *  t1)  - (x[5] *  t6) + (x[6] * R63) - (x[4] * R61);
    grad_op(ielem,2,2) = (x[3] *  t3) + (x[0] * R42) - (x[1] *  t2)  - (x[6] *  t4) + (x[7] * R74) - (x[5] * R72);
    grad_op(ielem,2,3) = (x[0] *  t5) - (x[1] * R31) - (x[2] *  t3)  + (x[7] *  t6) + (x[4] * R81) - (x[6] * R83);
    grad_op(ielem,2,4) = (x[5] *  t3) + (x[6] * R86) - (x[7] *  t2)  - (x[0] *  t4) - (x[3] * R81) + (x[1] * R61);
    grad_op(ielem,2,5) = (x[6] *  t5) - (x[4] *  t3)  - (x[7] * R75) + (x[1] *  t6) - (x[0] * R52) + (x[2] * R72);
    grad_op(ielem,2,6) = (x[7] *  t1) - (x[5] *  t5)  - (x[4] * R86) + (x[2] *  t4) - (x[1] * R63) + (x[3] * R83);
    grad_op(ielem,2,7) = (x[4] *  t2) - (x[6] *  t1)  + (x[5] * R75) - (x[3] *  t6) - (x[2] * R74) + (x[0] * R54);
  }

  //--------------------------------------------------------------------------

  static double test( int count )
  {
    array_type coord = Kokkos::create_mdarray< array_type >( count , 3 , 8 );
    array_type grad  = Kokkos::create_mdarray< array_type >( count , 3 , 8 );

    // Execute the parallel kernels on the arrays:

    Kokkos::parallel_for( count , HexSimpleFill<Scalar,device_type>( coord ) );

    double seconds ;

    Kokkos::parallel_for( count , HexGrad<Scalar,device_type>( coord , grad ) , seconds );

    return seconds ;
  }
};

