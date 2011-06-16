/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

namespace Test {

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

#define TEST_HEXGRAD_NORMAL 0

template< typename Scalar >
struct HexGrad< Scalar , KOKKOS_MACRO_DEVICE >
{
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef device_type::size_type  size_type ;

  enum { N_Space = 3 , N_Node = 8 };

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
        
    Scalar a[N_Node];
    
    //Z
    a[0] = coords(ielem,2,0);
    a[1] = coords(ielem,2,1);
    a[2] = coords(ielem,2,2);
    a[3] = coords(ielem,2,3);
    a[4] = coords(ielem,2,4);
    a[5] = coords(ielem,2,5);
    a[6] = coords(ielem,2,6);
    a[7] = coords(ielem,2,7);
    
    // z difference vectors
    Scalar R42=(a[3] - a[1]);
    Scalar R52=(a[4] - a[1]);
    Scalar R54=(a[4] - a[3]);

    Scalar R63=(a[5] - a[2]);
    Scalar R83=(a[7] - a[2]);
    Scalar R86=(a[7] - a[5]);
    
    Scalar R31=(a[2] - a[0]);
    Scalar R61=(a[5] - a[0]);
    Scalar R74=(a[6] - a[3]);

    Scalar R72=(a[6] - a[1]);
    Scalar R75=(a[6] - a[4]);
    Scalar R81=(a[7] - a[0]);

    Scalar t1=(R63 + R54);
    Scalar t2=(R61 + R74);
    Scalar t3=(R72 + R81);

    Scalar t4 =(R86 + R42);
    Scalar t5 =(R83 + R52);
    Scalar t6 =(R75 + R31);
    
    //Y
    a[0] = coords(ielem,1,0);
    a[1] = coords(ielem,1,1);
    a[2] = coords(ielem,1,2);
    a[3] = coords(ielem,1,3);
    a[4] = coords(ielem,1,4);
    a[5] = coords(ielem,1,5);
    a[6] = coords(ielem,1,6);
    a[7] = coords(ielem,1,7);


    grad_op(ielem,0,0) = (a[1] *  t1) - (a[2] * R42) - (a[3] *  t5)  + (a[4] *  t4) + (a[5] * R52) - (a[7] * R54);  
    grad_op(ielem,0,1) = (a[2] *  t2) + (a[3] * R31) - (a[0] *  t1)  - (a[5] *  t6) + (a[6] * R63) - (a[4] * R61);  
    grad_op(ielem,0,2) = (a[3] *  t3) + (a[0] * R42) - (a[1] *  t2)  - (a[6] *  t4) + (a[7] * R74) - (a[5] * R72);  
    grad_op(ielem,0,3) = (a[0] *  t5) - (a[1] * R31) - (a[2] *  t3)  + (a[7] *  t6) + (a[4] * R81) - (a[6] * R83);  
    grad_op(ielem,0,4) = (a[5] *  t3) + (a[6] * R86) - (a[7] *  t2)  - (a[0] *  t4) - (a[3] * R81) + (a[1] * R61);  
    grad_op(ielem,0,5) = (a[6] *  t5) - (a[4] *  t3)  - (a[7] * R75) + (a[1] *  t6) - (a[0] * R52) + (a[2] * R72);  
    grad_op(ielem,0,6) = (a[7] *  t1) - (a[5] *  t5)  - (a[4] * R86) + (a[2] *  t4) - (a[1] * R63) + (a[3] * R83);  
    grad_op(ielem,0,7) = (a[4] *  t2) - (a[6] *  t1)  + (a[5] * R75) - (a[3] *  t6) - (a[2] * R74) + (a[0] * R54);  

    
    R42=(a[3] - a[1]);
    R52=(a[4] - a[1]);
    R54=(a[4] - a[3]);

    R63=(a[5] - a[2]);
    R83=(a[7] - a[2]);
    R86=(a[7] - a[5]);

    R31=(a[2] - a[0]);
    R61=(a[5] - a[0]);
    R74=(a[6] - a[3]);

    R72=(a[6] - a[1]);
    R75=(a[6] - a[4]);
    R81=(a[7] - a[0]);

    t1=(R63 + R54);
    t2=(R61 + R74);
    t3=(R72 + R81);

    t4 =(R86 + R42);
    t5 =(R83 + R52);
    t6 =(R75 + R31);

    //X
    a[0] = coords(ielem,0,0);
    a[1] = coords(ielem,0,1);
    a[2] = coords(ielem,0,2);
    a[3] = coords(ielem,0,3);
    a[4] = coords(ielem,0,4);
    a[5] = coords(ielem,0,5);
    a[6] = coords(ielem,0,6);
    a[7] = coords(ielem,0,7);
    
	// Z grad
    grad_op(ielem,1,7) = (a[4] *  t2) - (a[6] *  t1)  + (a[5] * R75) - (a[3] *  t6) - (a[2] * R74) + (a[0] * R54);  
    grad_op(ielem,1,6) = (a[7] *  t1) - (a[5] *  t5)  - (a[4] * R86) + (a[2] *  t4) - (a[1] * R63) + (a[3] * R83);  
    grad_op(ielem,1,5) = (a[6] *  t5) - (a[4] *  t3)  - (a[7] * R75) + (a[1] *  t6) - (a[0] * R52) + (a[2] * R72);  
    grad_op(ielem,1,4) = (a[5] *  t3) + (a[6] * R86) - (a[7] *  t2)  - (a[0] *  t4) - (a[3] * R81) + (a[1] * R61);  
    grad_op(ielem,1,3) = (a[0] *  t5) - (a[1] * R31) - (a[2] *  t3)  + (a[7] *  t6) + (a[4] * R81) - (a[6] * R83);  
    grad_op(ielem,1,2) = (a[3] *  t3) + (a[0] * R42) - (a[1] *  t2)  - (a[6] *  t4) + (a[7] * R74) - (a[5] * R72);  
    grad_op(ielem,1,1) = (a[2] *  t2) + (a[3] * R31) - (a[0] *  t1)  - (a[5] *  t6) + (a[6] * R63) - (a[4] * R61);  
    grad_op(ielem,1,0) = (a[1] *  t1) - (a[2] * R42) - (a[3] *  t5)  + (a[4] *  t4) + (a[5] * R52) - (a[7] * R54);  
    

    R42=(a[3] - a[1]);
    R52=(a[4] - a[1]);
    R54=(a[4] - a[3]);

    R63=(a[5] - a[2]);
    R83=(a[7] - a[2]);
    R86=(a[7] - a[5]);

    R31=(a[2] - a[0]);
    R61=(a[5] - a[0]);
    R74=(a[6] - a[3]);

    R72=(a[6] - a[1]);
    R75=(a[6] - a[4]);
    R81=(a[7] - a[0]);
    
    t1=(R63 + R54);
    t2=(R61 + R74);
    t3=(R72 + R81);

    t4 =(R86 + R42);
    t5 =(R83 + R52);
    t6 =(R75 + R31);
    
    //Z
    a[0] = coords(ielem,2,0);
    a[1] = coords(ielem,2,1);
    a[2] = coords(ielem,2,2);
    a[3] = coords(ielem,2,3);
    a[4] = coords(ielem,2,4);
    a[5] = coords(ielem,2,5);
    a[6] = coords(ielem,2,6);
    a[7] = coords(ielem,2,7);


    grad_op(ielem,2,0)  = (a[1] *  t1) - (a[2] * R42) - (a[3] *  t5)  + (a[4] *  t4) + (a[5] * R52) - (a[7] * R54); 
    grad_op(ielem,2,1)  = (a[2] *  t2) + (a[3] * R31) - (a[0] *  t1)  - (a[5] *  t6) + (a[6] * R63) - (a[4] * R61); 
    grad_op(ielem,2,2)  = (a[3] *  t3) + (a[0] * R42) - (a[1] *  t2)  - (a[6] *  t4) + (a[7] * R74) - (a[5] * R72); 
    grad_op(ielem,2,3)  = (a[0] *  t5) - (a[1] * R31) - (a[2] *  t3)  + (a[7] *  t6) + (a[4] * R81) - (a[6] * R83); 
    grad_op(ielem,2,4)  = (a[5] *  t3) + (a[6] * R86) - (a[7] *  t2)  - (a[0] *  t4) - (a[3] * R81) + (a[1] * R61); 
    grad_op(ielem,2,5)  = (a[6] *  t5) - (a[4] *  t3)  - (a[7] * R75) + (a[1] *  t6) - (a[0] * R52) + (a[2] * R72); 
    grad_op(ielem,2,6)  = (a[7] *  t1) - (a[5] *  t5)  - (a[4] * R86) + (a[2] *  t4) - (a[1] * R63) + (a[3] * R83); 
    grad_op(ielem,2,7)  = (a[4] *  t2) - (a[6] *  t1)  + (a[5] * R75) - (a[3] *  t6) - (a[2] * R74) + (a[0] * R54); 
    
  }

  //--------------------------------------------------------------------------

  static double test( int count )
  {
    array_type coord = Kokkos::create_mdarray< array_type >( count , 3 , 8 );
    array_type grad  = Kokkos::create_mdarray< array_type >( count , 3 , 8 );

    // Execute the parallel kernels on the arrays:

    double seconds = 0.0;
    
    Kokkos::parallel_for( count , HexSimpleFill<Scalar,device_type>( coord ) , seconds );

    Kokkos::parallel_for( count , HexGrad<Scalar,device_type>( coord , grad ) , seconds );

    return seconds ;
  }
};

}

