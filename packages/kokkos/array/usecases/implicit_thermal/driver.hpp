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

#include <cstddef>
#include <SparseOutput.hpp>
#include <ImplicitBoxMeshFixture.hpp>
#include <CRSGraph.hpp>

#define  PRINT_SAMPLE_OF_SOLUTION  0

namespace Test {

template< typename Scalar , class Device >
struct MiniImplTherm ;

struct PerformanceData {
  double mesh_time ;
  double elem_time ;
  size_t elem_flop ;
  double fill_time ;
  size_t fill_flop ;
  double solve_mflop_per_sec ;

  PerformanceData()
  : mesh_time(0)
  , elem_time(0)
  , elem_flop(0)
  , fill_time(0)
  , fill_flop(0)
  , solve_mflop_per_sec(0)
  {}

  void best( const PerformanceData & rhs )
  {
    if ( rhs.mesh_time < mesh_time ) mesh_time = rhs.mesh_time ;
    if ( rhs.elem_time < elem_time ) elem_time = rhs.elem_time ;
    if ( rhs.fill_time < fill_time ) fill_time = rhs.fill_time ;
    if ( solve_mflop_per_sec < rhs.solve_mflop_per_sec )
      solve_mflop_per_sec = rhs.solve_mflop_per_sec ;
  }
};

template< typename Scalar >
struct MiniImplTherm< Scalar , KOKKOS_MACRO_DEVICE > {

static void run(int x, int y, int z, PerformanceData & perf )
{
  typedef KOKKOS_MACRO_DEVICE    device_type;
  typedef device_type::size_type index_type ;

  typedef KokkosArray::MDArray<Scalar,     device_type>  scalar_array_d;
  typedef KokkosArray::MDArray<index_type, device_type>  index_array_d;

  typedef KokkosArray::MultiVector<Scalar,     device_type>  scalar_vector_d;
  typedef KokkosArray::MultiVector<index_type, device_type>  index_vector_d;

  typedef typename scalar_array_d::HostMirror  scalar_array_h ;
  typedef typename index_array_d ::HostMirror  index_array_h ;

  typedef typename scalar_vector_d::HostMirror  scalar_vector_h ;
  typedef typename index_vector_d ::HostMirror  index_vector_h ;

  // Problem coefficients

  const Scalar elem_coeff_K = 2 ;
  const Scalar elem_load_Q  = 1 ;

  //  Host Data Structures

  index_array_h   elem_graph_col_h ;
  index_vector_h  A_row_h , A_col_h ;

  //  Device Data Structures

  scalar_array_d  elem_stiffness, elem_load ;

  index_array_d  elem_graph_col_d ;
  index_vector_d A_row_d , A_col_d, dirichlet_flag_d ;
  scalar_vector_d A, b, X , dirichlet_value_d ;

  //------------------------------
  // Generate mesh and corresponding sparse matrix graph

  KokkosArray::Impl::Timer wall_clock ;

  const BoxMeshFixture< double , device_type > mesh( x , y , z );

  mesh.init_dirichlet_z( dirichlet_flag_d , dirichlet_value_d );

  elem_graph_col_d =
    KokkosArray::create_mdarray< index_array_d >(
      mesh.h_mesh.elem_node_ids.dimension(0) ,
      mesh.h_mesh.elem_node_ids.dimension(1) ,
      mesh.h_mesh.elem_node_ids.dimension(1) );

  elem_graph_col_h = create_mirror( elem_graph_col_d );

  init_crsgraph( mesh.h_mesh.node_elem_offset ,
                 mesh.h_mesh.node_elem_ids ,
                 mesh.h_mesh.elem_node_ids ,
                 elem_graph_col_h ,
                 A_row_h , A_col_h );

  deep_copy( elem_graph_col_d , elem_graph_col_h );

  // Copy sparse matrix graph to device

  A_row_d = KokkosArray::create_multivector< index_array_d >("A_row_d",A_row_h.length());
  A_col_d = KokkosArray::create_multivector< index_array_d >("A_col_d",A_col_h.length());

  KokkosArray::deep_copy(A_row_d, A_row_h);
  KokkosArray::deep_copy(A_col_d, A_col_h);

  device_type::fence();

  perf.mesh_time = wall_clock.seconds(); // Mesh and graph allocation and population.

  //------------------------------
  // Allocate device memory for linear system and element contributions.

  A = KokkosArray::create_multivector< scalar_vector_d > ("A",A_col_h.length());
  b = KokkosArray::create_multivector< scalar_vector_d > ("b",mesh.node_count);
  X = KokkosArray::create_multivector< scalar_vector_d > ("X",mesh.node_count);

  elem_stiffness =  KokkosArray::create_mdarray< scalar_array_d > (mesh.elem_count, 8, 8);
  elem_load      =  KokkosArray::create_mdarray< scalar_array_d > (mesh.elem_count, 8);

  wall_clock.reset();

  typedef ElementComp< Scalar , double , device_type > ElementFunctor ;
  typedef CRSMatrixGatherFill<Scalar, device_type> GatherFillFunctor ;

  KokkosArray::parallel_for( mesh.elem_count,
    ElementFunctor( mesh.d_mesh.elem_node_ids ,
                    mesh.d_mesh.node_coords ,
                    elem_stiffness, elem_load ,
                    elem_coeff_K , elem_load_Q ) );

  device_type::fence();

  // Element computation time and flops
  perf.elem_time = wall_clock.seconds();
  perf.elem_flop = (size_t) ElementFunctor::FLOPS_operator *
                   (size_t) mesh.elem_count ;

  wall_clock.reset();

  KokkosArray::parallel_for( mesh.node_count,
    GatherFillFunctor( A, b, A_row_d, A_col_d,
                       mesh.d_mesh.node_elem_offset ,
                       mesh.d_mesh.node_elem_ids,
                       elem_graph_col_d,
                       elem_stiffness,
                       elem_load ) );

  device_type::fence();

  // Matrix gather-fill time and flops
  perf.fill_time = wall_clock.seconds();
  perf.fill_flop = (size_t) elem_stiffness.size() + (size_t) elem_load.size();

  //------------------------------

  KokkosArray::parallel_for(mesh.node_count ,
    Dirichlet<Scalar , device_type>(A, A_row_d ,A_col_d, b,
                                    dirichlet_flag_d , dirichlet_value_d ) );


//  printSparse< scalar_vector_d , index_vector_d>("A.txt",A,A_row_d,A_col_d);

  //------------------------------
  // Solve linear sytem

  perf.solve_mflop_per_sec = CG_Solve<Scalar, device_type>::run(A , A_row_d, A_col_d , b , X );

#if  PRINT_SAMPLE_OF_SOLUTION

  scalar_vector_h X_h = KokkosArray::create_mirror( X );

  KokkosArray::deep_copy( X_h , X );

  for ( int i = 0 ; i < (int) mesh.node_count_z ; ++i ) {
    const int ix = mesh.node_count_x - 1 ;
    const int iy = mesh.node_count_y - 1 ;
    const Scalar val00 = X_h( mesh.node_id( 0 , 0 , i ) );
    const Scalar valX0 = X_h( mesh.node_id( ix , 0 , i ) );
    const Scalar valXY = X_h( mesh.node_id( ix , iy , i ) );
    const Scalar val0Y = X_h( mesh.node_id( 0 , iy , i ) );
    std::cout << "corners_00_X0_XY_Y0("<<i<<") = "
              << val00 << " " << valX0 << " "
              << valXY << " " << val0Y << std::endl ;
  }

#endif

//  printGLUT<Scalar , scalar_vector_d , scalar_array_h , index_array_h>
//      ("X.txt", X , elem_coords_h , elem_nodeIDs_h,x,y,z);
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

static void driver( const char * label , int beg , int end , int runs )
{
  std::cout << std::endl ;
  std::cout << "\"MiniImplTherm with KokkosArray " << label << "\"" << std::endl;
  std::cout << "\"Size\" ,     \"Setup\" ,    \"Element\" ,  \"Element\" , \"Fill\" ,   \"Fill\" ,  \"Solve\"" << std::endl
            << "\"elements\" , \"millisec\" , \"millisec\" , \"flops\" , \"millisec\" , \"flops\" , \"Mflop/sec\"" << std::endl ;

  for(int i = beg ; i < end; ++i )
  {
    const int ix = (int) cbrt( (double) ( 1 << i ) );
    const int iy = ix + 1 ;
    const int iz = iy + 1 ;
    const int n  = ix * iy * iz ;

    PerformanceData perf , best ;

    for(int j = 0; j < runs; j++){

     run(ix,iy,iz,perf);

     if( j == 0 ) {
       best = perf ;
     }
     else {
       best.best( perf );
     }
   }
   std::cout << std::setw(8) << n << " , "
             << std::setw(10) << best.mesh_time * 1000 << " , "
             << std::setw(10) << best.elem_time * 1000 << " , "
             << std::setw(10) << best.elem_flop << " , "
             << std::setw(10) << best.fill_time * 1000 << " , "
             << std::setw(10) << best.fill_flop << " , "
             << std::setw(10) << best.solve_mflop_per_sec
             << std::endl ;
  }
}

};

} // namespace Test

