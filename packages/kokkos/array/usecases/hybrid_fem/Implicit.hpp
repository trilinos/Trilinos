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

#ifndef HYBRIDFEM_IMPLICIT_HPP
#define HYBRIDFEM_IMPLICIT_HPP

#include <utility>
#include <iostream>
#include <iomanip>

#include <KokkosArray_View.hpp>
#include <SparseLinearSystem.hpp>
#include <SparseLinearSystemFill.hpp>
#include <FEMesh.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace HybridFEM {
namespace Implicit {

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

template< typename ScalarType , typename ScalarCoordType , class Device >
struct ElementComputation ;

template< typename ScalarType , typename ScalarCoordType , class Device >
struct DirichletBoundary ;

//----------------------------------------------------------------------------

template< typename Scalar , class FixtureType >
PerformanceData run( const typename FixtureType::FEMeshType & mesh ,
                     const int global_max_x ,
                     const int global_max_y ,
                     const int global_max_z ,
                     const bool print_sample )
{
  typedef Scalar                              scalar_type ;
  typedef FixtureType                         fixture_type ;
  typedef typename fixture_type::device_type  device_type;
  typedef typename device_type::size_type     size_type ;

  typedef typename fixture_type::FEMeshType mesh_type ;
  typedef typename fixture_type::coordinate_scalar_type coordinate_scalar_type ;

  enum { ElementNodeCount = fixture_type::element_node_count };

  const comm::Machine machine = mesh.parallel_data_map.machine ;

  const size_t element_count = mesh.elem_node_ids.dimension(0);

  const size_t iteration_limit = 200 ;
  const double residual_tolerance = 1e-14 ;

  size_t iteration_count = 0 ;
  double residual_norm = 0 ;

  PerformanceData perf_data ;

  //------------------------------------
  // Sparse linear system types:

  typedef KokkosArray::View< Scalar[] , device_type >   vector_type ;
  typedef KokkosArray::CrsMatrix< Scalar , device_type >     matrix_type ;
  typedef typename matrix_type::graph_type         matrix_graph_type ;
  typedef typename matrix_type::coefficients_type  matrix_coefficients_type ;

  typedef GraphFactory< matrix_graph_type , mesh_type > graph_factory ;

  //------------------------------------
  // Problem setup types:

  typedef ElementComputation< Scalar , Scalar , device_type > ElementFunctor ;
  typedef DirichletBoundary< Scalar , Scalar , device_type > BoundaryFunctor ;

  typedef typename ElementFunctor::elem_matrices_type elem_matrices_type ;
  typedef typename ElementFunctor::elem_vectors_type  elem_vectors_type ;

  typedef GatherFill< matrix_type ,
                      mesh_type ,
                      elem_matrices_type ,
                      elem_vectors_type > GatherFillFunctor ;

  //------------------------------------

  const Scalar elem_coeff_K = 2 ;
  const Scalar elem_load_Q  = 1 ;

  matrix_type linsys_matrix ;
  vector_type linsys_rhs ;
  vector_type linsys_solution ;

  typename graph_factory::element_map_type element_map ;

  KokkosArray::Impl::Timer wall_clock ;

  //------------------------------------
  // Generate sparse matrix graph and element->graph map.

  graph_factory::create( mesh , linsys_matrix.graph , element_map );

  device_type::fence();
  perf_data.graph_time = comm::max( machine , wall_clock.seconds() );

  //------------------------------------
  // Allocate linear system coefficients and rhs:

  const size_t local_owned_length =
    linsys_matrix.graph.row_map.dimension(0) - 1 ;

  linsys_matrix.coefficients =
    matrix_coefficients_type( "coeff" , linsys_matrix.graph.entries.dimension(0) );

  linsys_rhs      = vector_type( "rhs" , local_owned_length );
  linsys_solution = vector_type( "solution" , local_owned_length );

  //------------------------------------
  // Fill linear system
  {
    elem_matrices_type elem_matrices ;
    elem_vectors_type  elem_vectors ;

    if ( element_count ) {
      elem_matrices = elem_matrices_type( std::string("elem_matrices"), element_count );
      elem_vectors  = elem_vectors_type ( std::string("elem_vectors"), element_count );
    }

    //------------------------------------
    // Compute element matrices and vectors:

    wall_clock.reset();

    ElementFunctor::apply( mesh ,
                           elem_matrices , elem_vectors ,
                           elem_coeff_K , elem_load_Q );

    device_type::fence();
    perf_data.elem_time = comm::max( machine , wall_clock.seconds() );

    //------------------------------------
    // Fill linear system coefficients:

    wall_clock.reset();

    GatherFillFunctor::apply( linsys_matrix , linsys_rhs ,
               mesh , element_map , elem_matrices , elem_vectors );

    device_type::fence();
    perf_data.matrix_gather_fill_time = comm::max( machine , wall_clock.seconds() );

    // Apply boundary conditions:

    wall_clock.reset();

    BoundaryFunctor::apply( linsys_matrix , linsys_rhs , mesh , 
                            0 , global_max_z , 0 , global_max_z );

    device_type::fence();
    perf_data.matrix_boundary_condition_time = comm::max( machine , wall_clock.seconds() );
  }

  //------------------------------------
  // Solve linear sytem

  cgsolve( mesh.parallel_data_map ,
           linsys_matrix , linsys_rhs , linsys_solution ,
           iteration_count , residual_norm ,
           perf_data.cg_iteration_time ,
           iteration_limit , residual_tolerance );

  //------------------------------------

  if ( print_sample ) {

    typename mesh_type::node_coords_type::HostMirror coords_h =
      KokkosArray::create_mirror( mesh.node_coords );

    typename vector_type::HostMirror X_h =
      KokkosArray::create_mirror( linsys_solution );

    KokkosArray::deep_copy( coords_h , mesh.node_coords );
    KokkosArray::deep_copy( X_h , linsys_solution );

    for ( size_t i = 0 ; i < mesh.parallel_data_map.count_owned ; ++i ) {
      const coordinate_scalar_type x = coords_h(i,0);
      const coordinate_scalar_type y = coords_h(i,1);
      const coordinate_scalar_type z = coords_h(i,2);

      if ( x <= 0 && y <= 0 ) {
        std::cout << "  node( " << x << " " << y << " " << z << " ) = "
                  << X_h(i) << std::endl ;
      }
    }
  }

  return perf_data ;
}

//----------------------------------------------------------------------------

template< typename Scalar , class Device >
void driver( const char * label ,
             comm::Machine machine ,
             const int elem_count_beg ,
             const int elem_count_end , int runs )
{
  typedef Scalar              scalar_type ;
  typedef Device              device_type ;
  typedef double              coordinate_scalar_type ;
  typedef FixtureElementHex8  fixture_element_type ;

  typedef BoxMeshFixture< coordinate_scalar_type ,
                          device_type ,
                          fixture_element_type > fixture_type ;

  typedef typename fixture_type::FEMeshType mesh_type ;

  if ( elem_count_beg == 0 || elem_count_end == 0 || runs == 0 ) return ;

  if ( comm::rank( machine ) == 0 ) {
    std::cout << std::endl ;
    std::cout << "\"KokkosArray::HybridFE::Implicit " << label << "\"" << std::endl;
    std::cout << "\"Size\" ,  \"Graphing\" , \"Element\" , \"Fill\" ,   \"Boundary\" ,  \"CG-Iter\"" << std::endl
              << "\"elems\" , \"millisec\" , \"millisec\" , \"millisec\" , \"millisec\" , \"millisec\"" << std::endl ;
  }

  for(int i = elem_count_beg ; i < elem_count_end ; i *= 2 )
  {
    const int ix = std::max( 1 , (int) cbrt( ((double) i) / 2.0 ) );
    const int iy = ix + 1 ;
    const int iz = 2 * iy ;
    const int n  = ix * iy * iz ;

    mesh_type mesh =
      fixture_type::create( comm::size( machine ) ,
                            comm::rank( machine ) ,
                            ix , iy , iz );

    mesh.parallel_data_map.machine = machine ;

    PerformanceData perf_data , perf_best ;

    for(int j = 0; j < runs; j++){

     perf_data = run<Scalar,fixture_type>(mesh,ix,iy,iz, false );

     if( j == 0 ) {
       perf_best = perf_data ;
     }
     else {
       perf_best.best( perf_data );
     }
   }

  if ( comm::rank( machine ) == 0 ) {

     std::cout << std::setw(8) << n << " , "
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

} /* namespace Implicit */
} /* namespace HybridFEM */


#endif /* #ifndef HYBRIDFEM_IMPLICIT_HPP */

