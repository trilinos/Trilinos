/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
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

#ifndef TESTIMPLICIT_HPP
#define TESTIMPLICIT_HPP

#include <utility>
#include <iostream>
#include <iomanip>

#include <Kokkos_Core.hpp>

#include <FEMesh.hpp>
#include <BoxMeshFixture.hpp>

#include <TestGenerateSystem.hpp>
#include <TestSparseLinearSystem.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Test {

struct PerformanceData {
  double mesh_time ;
  double graph_time ;
  double elem_time ;
  double matrix_gather_fill_time ;
  double matrix_boundary_condition_time ;
  double cg_iteration_time ;

  PerformanceData()
    : mesh_time(0)
    , graph_time(0)
    , elem_time(0)
    , matrix_gather_fill_time(0)
    , matrix_boundary_condition_time(0)
    , cg_iteration_time(0)
    {}

  void best( const PerformanceData & rhs )
  {
    mesh_time = std::min( mesh_time , rhs.mesh_time );
    graph_time = std::min( graph_time , rhs.graph_time );
    elem_time = std::min( elem_time , rhs.elem_time );
    matrix_gather_fill_time = std::min( matrix_gather_fill_time , rhs.matrix_gather_fill_time );
    matrix_boundary_condition_time = std::min( matrix_boundary_condition_time , rhs.matrix_boundary_condition_time );
    cg_iteration_time = std::min( cg_iteration_time , rhs.cg_iteration_time );
  }
};

//----------------------------------------------------------------------------

template< typename Scalar , class FixtureType >
PerformanceData implicit_run( const typename FixtureType::FEMeshType & mesh , const bool verify )
{
  typedef Scalar                              scalar_type ;
  typedef FixtureType                         fixture_type ;
  typedef typename fixture_type::device_type  device_type;
  typedef typename device_type::size_type     size_type ;

  typedef typename fixture_type::FEMeshType mesh_type ;
  typedef typename fixture_type::coordinate_scalar_type coordinate_scalar_type ;

  enum { ElementNodeCount = fixture_type::element_node_count };

  const comm::Machine machine = mesh.parallel_data_map.machine ;

  const size_t iteration_limit = verify ? 2000 : 200 ;
  const double residual_tolerance = 1e-14 ;

  size_t iteration_count = 0 ;
  double residual_norm = 0 ;

  PerformanceData perf_data ;

  //------------------------------------
  // Sparse linear system types:

  typedef Kokkos::View< Scalar* , Kokkos::LayoutRight , device_type >   vector_type ;
  typedef Kokkos::CrsMatrix< Scalar , device_type >     matrix_type ;
  typedef typename matrix_type::graph_type    matrix_graph_type ;
  typedef typename matrix_type::values_type   matrix_coefficients_type ;

  //------------------------------------

  matrix_type linsys_matrix ;
  vector_type linsys_rhs ;
  vector_type linsys_x ;
  typename vector_type::HostMirror linsys_host_solution ;

  Kokkos::Impl::Timer wall_clock ;

  //------------------------------------

  linsys_matrix.graph = create_graph_from_mesh< matrix_graph_type >( mesh );

  device_type::fence();
  perf_data.graph_time = comm::max( machine , wall_clock.seconds() );

  //------------------------------------
  // Allocate linear system coefficients and rhs:

  const size_t local_owned_length =
    linsys_matrix.graph.row_map.dimension_0() - 1 ;

  linsys_matrix.values =
    matrix_coefficients_type( "coeff" , linsys_matrix.graph.entries.dimension_0() );

  linsys_rhs      = vector_type( "rhs" , local_owned_length );
  linsys_x        = vector_type( "x" ,   local_owned_length );

  linsys_host_solution = typename vector_type::HostMirror( "host_solution" , local_owned_length );

  //------------------------------------
  // Fill linear system such that: linsys_rhs = linsys_matrix * linsys_solution ;

  fill_linear_system( mesh , linsys_matrix , linsys_rhs , linsys_host_solution );

  //------------------------------------
  // Solve linear sytem

  cgsolve( mesh.parallel_data_map ,
           linsys_matrix , linsys_rhs , linsys_x ,
           iteration_count , residual_norm ,
           perf_data.cg_iteration_time ,
           iteration_limit , residual_tolerance );

  //------------------------------------

  if ( verify ) {

    const double tolerance = Kokkos::Impl::is_same<scalar_type,double>::value ? 1.0e-10 : 1.0e-4 ;

    typename vector_type::HostMirror linsys_host_x = Kokkos::create_mirror_view( linsys_x );

    Kokkos::deep_copy( linsys_host_x , linsys_x );

    const int comm_rank = comm::rank( mesh.parallel_data_map.machine );
    const int comm_size = comm::size( mesh.parallel_data_map.machine );

    size_t error_count = 0 ;

    for ( unsigned i = 0 ; i < linsys_host_solution.dimension_0() ; ++i ) {
    for ( unsigned j = 0 ; j < linsys_host_solution.dimension_1() ; ++j ) {

      const double diff = linsys_host_x(i,j) - linsys_host_solution(i,j);
      const bool large  = tolerance < fabs( linsys_host_solution(i,j) );
      const bool error  = tolerance < fabs( diff ) / ( large ? fabs( linsys_host_solution(i,j) ) : 1 );

      if ( error ) ++error_count ;

      if ( error && 0 < error_count && error_count < 10 ) {
        std::cout << "P" << comm_rank
                  << ":  error(" << i << "," << j << ") = " << diff
                  << " : x = " << linsys_host_x(i,j)
                  << " != " << linsys_host_solution(i,j) << " = solution"
                  << std::endl ;
      }
    }
    }

#if defined( KOKKOS_HAVE_MPI )
  long global[2] = { 0 , 0 };
  long local[2] = { linsys_host_solution.dimension_0() , error_count };
  MPI_Allreduce( local , global , 2 , MPI_LONG , MPI_SUM , mesh.parallel_data_map.machine.mpi_comm );
#else
  long global[2] = { linsys_host_solution.dimension_0() , error_count };
#endif

    if ( 0 == comm_rank && 0 != global[1] ) {
      std::cout << "\"Verify: Comm[" << comm_size << "]"
                << " , Vector[" << global[0] << "x" << linsys_host_solution.dimension_1() << "]"
                << " , Error count = " << global[1]
                << " , CG-residual = " << residual_norm
                << " , CG-iteration = " << iteration_count
                << "\"" << std::endl ;
    }
  }

  return perf_data ;
}

//----------------------------------------------------------------------------

template< typename Scalar , class Device >
void implicit_driver( const char * const label ,
                      comm::Machine machine ,
                      const int gang_count ,
                      const int elem_count_beg ,
                      const int elem_count_end ,
                      const int runs )
{
  typedef Scalar              scalar_type ;
  typedef Device              device_type ;
  typedef double              coordinate_scalar_type ;
  typedef FixtureElementHex8  fixture_element_type ;

  typedef BoxMeshFixture< coordinate_scalar_type ,
                          device_type ,
                          fixture_element_type > fixture_type ;

  typedef typename fixture_type::FEMeshType mesh_type ;

  const size_t proc_count = comm::size( machine );
  const size_t proc_rank  = comm::rank( machine );

  if ( elem_count_beg == 0 || elem_count_end == 0 || runs == 0 ) return ;

  //------------------------------------
  // Verification:
  {
    const int ix = 5 ;
    const int iy = ix + 1 ;
    const int iz = 2 * iy ;

    mesh_type mesh =
      fixture_type::create( proc_count , proc_rank , gang_count , ix , iy , iz );

    mesh.parallel_data_map.machine = machine ;

    implicit_run<Scalar,fixture_type>(mesh,true);
  }
  //------------------------------------

  if ( comm::rank( machine ) == 0 ) {
    std::cout << std::endl ;
    std::cout << "\"Kokkos::HybridFE::Implicit " << label << "\"" << std::endl;
    std::cout << "\"Size\" ,  \"Graphing\" , \"Element\" , \"Fill\" ,   \"Boundary\" ,  \"CG-Iter\"" << std::endl
              << "\"nodes\" , \"millisec\" , \"millisec\" , \"millisec\" , \"millisec\" , \"millisec\"" << std::endl ;
  }

  for(int i = elem_count_beg ; i < elem_count_end ; i *= 2 )
  {
    const int ix = std::max( 1 , (int) cbrt( ((double) i) / 2.0 ) );
    const int iy = ix + 1 ;
    const int iz = 2 * iy ;
    const int nN = ( ix + 1 ) * ( iy + 1 ) * ( iz + 1 );

    mesh_type mesh =
      fixture_type::create( proc_count , proc_rank , gang_count , ix , iy , iz );

    mesh.parallel_data_map.machine = machine ;

    PerformanceData perf_data , perf_best ;

    for(int j = 0; j < runs; j++){

     perf_data = implicit_run<Scalar,fixture_type>(mesh,false);

     if( j == 0 ) {
       perf_best = perf_data ;
     }
     else {
       perf_best.best( perf_data );
     }
   }

  if ( comm::rank( machine ) == 0 ) {

     std::cout << std::setw(8) << nN << " , "
               << std::setw(10) << perf_best.graph_time * 1000 << " , "
               << std::setw(10) << perf_best.elem_time * 1000 << " , "
               << std::setw(10) << perf_best.matrix_gather_fill_time * 1000 << " , "
               << std::setw(10) << perf_best.matrix_boundary_condition_time * 1000 << " , "
               << std::setw(10) << perf_best.cg_iteration_time * 1000
               << std::endl ;
    }
  }
}

//----------------------------------------------------------------------------

}


#endif /* #ifndef TESTIMPLICIT_HPP */

