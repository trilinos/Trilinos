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

namespace Test {

template< class DeviceType > struct HexSimpleFill ;

//each element is a unit cube
template<>
struct HexSimpleFill< KOKKOS_MACRO_DEVICE >
{
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef device_type::size_type  size_type ;

  // 3D array : ( ParallelWork , Space , Node )
  typedef KokkosArray::MDArray<double,device_type> d_array ;

  KokkosArray::MDArray<double,device_type> coords ;

  HexSimpleFill( const KokkosArray::MDArray<double,device_type> & arg_coords )
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

template< class DeviceType > struct HexGrad ;

#define TEST_HEXGRAD_NORMAL 0

template<>
struct HexGrad< KOKKOS_MACRO_DEVICE >
{
  typedef KOKKOS_MACRO_DEVICE     device_type ;
  typedef device_type::size_type  size_type ;

  enum { N_Space = 3 , N_Node = 8 };

  // 3D array : ( ParallelWork , Space , Node )

  KokkosArray::MDArray<double,device_type> coords ;
  KokkosArray::MDArray<float,device_type>  grad_op ;

  HexGrad( const KokkosArray::MDArray<double,device_type> & arg_coords ,
           const KokkosArray::MDArray<float,device_type>  & arg_grad_op )
    : coords( arg_coords )
    , grad_op( arg_grad_op )
    {}

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type ielem ) const
  {
    // Repeated re-use of nodal coordinates,
    // copy them into local storage.
        
    double a[N_Node];
    
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
    float R42=(a[3] - a[1]);
    float R52=(a[4] - a[1]);
    float R54=(a[4] - a[3]);

    float R63=(a[5] - a[2]);
    float R83=(a[7] - a[2]);
    float R86=(a[7] - a[5]);
    
    float R31=(a[2] - a[0]);
    float R61=(a[5] - a[0]);
    float R74=(a[6] - a[3]);

    float R72=(a[6] - a[1]);
    float R75=(a[6] - a[4]);
    float R81=(a[7] - a[0]);

    float t1=(R63 + R54);
    float t2=(R61 + R74);
    float t3=(R72 + R81);

    float t4 =(R86 + R42);
    float t5 =(R83 + R52);
    float t6 =(R75 + R31);
    
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
    typedef KokkosArray::MDArray<double,device_type> darray ;
    typedef KokkosArray::MDArray<float,device_type>  farray ;
    darray coord = KokkosArray::create_mdarray< darray >( count , 3 , 8 );
    farray grad  = KokkosArray::create_mdarray< farray >( count , 3 , 8 );

    // Execute the parallel kernels on the arrays:

    double seconds = 0.0;
    
    KokkosArray::parallel_for( count , HexSimpleFill<device_type>( coord ) , seconds );

    KokkosArray::parallel_for( count , HexGrad<device_type>( coord , grad ) , seconds );

    return seconds ;
  }
};

}

