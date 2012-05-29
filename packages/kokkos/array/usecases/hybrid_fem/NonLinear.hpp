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

#ifndef HYBRIDFEM_NONLINEAR_HPP
#define HYBRIDFEM_NONLINEAR_HPP

#include <utility>
#include <iostream>
#include <iomanip>

#include <Kokkos_Array.hpp>
#include <Kokkos_MultiVector.hpp>
#include <SparseLinearSystem.hpp>
#include <SparseLinearSystemFill.hpp>
#include <FEMesh.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace HybridFEM {
namespace NonLinear {

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
struct DirichletSolution ;

template< typename ScalarType , typename ScalarCoordType , class Device >
struct DirichletResidual ;

//----------------------------------------------------------------------------

template< typename Scalar , class Device >
PerformanceData run( comm::Machine machine ,
                     const int global_count_x ,
                     const int global_count_y ,
                     const int global_count_z ,
                     const bool print_sample )
{
  typedef Device                           device_type;
  typedef typename device_type::size_type  size_type ;
  typedef Scalar                           scalar_type ;

  //-----------------------------------
  // Convergence Criteria and perf data:

  const size_t cg_iteration_limit = 200 ;
  const double cg_tolerance = 1e-14 ;

  const size_t newton_iteration_limit = 150 ;
  const double newton_tolerance = 1e-14 ;

  size_t cg_iteration_count_total = 0 ;
  size_t cg_iteration_count = 0 ;
  double cg_iteration_time = 0 ;
  double cg_residual_norm = 0 ;

  size_t newton_iteration_count = 0 ;
  double delta_norm = 0 ;

  bool not_converged = true;

  PerformanceData perf_data ;

  //------------------------------------
  // FEMesh types:

  enum { ElementNodeCount = 8 };
  typedef double coordinate_scalar_type ;
  typedef FEMesh< coordinate_scalar_type , ElementNodeCount , device_type > mesh_type ;

  //------------------------------------
  // Sparse linear system types:

  typedef Kokkos::MultiVector< Scalar , Device >   vector_type ;
  typedef Kokkos::CrsMatrix< Scalar , Device >     matrix_type ;
  typedef typename matrix_type::graph_type         matrix_graph_type ;
  typedef typename matrix_type::coefficients_type  matrix_coefficients_type ;

  typedef Kokkos::Impl::Factory< matrix_graph_type , mesh_type > graph_factory ;

  //------------------------------------
  // Problem setup types:

  typedef ElementComputation < Scalar , Scalar , device_type >  ElementFunctor ;
  typedef DirichletSolution  < Scalar , Scalar , device_type > DirichletSolutionFunctor ;
  typedef DirichletResidual  < Scalar , Scalar , device_type > DirichletResidualFunctor ;

  //------------------------------------

  const Scalar elem_coeff_K = 1 ;
  const Scalar elem_coeff_P = 0.6 ;

  matrix_type jacobian ;
  vector_type residual ;
  vector_type delta ;
  vector_type nodal_solution ;

  typename graph_factory::element_map_type element_map ;

  //------------------------------------
  // Generate mesh and corresponding sparse matrix graph

  Kokkos::Impl::Timer wall_clock ;

  mesh_type mesh =
    box_mesh_fixture< coordinate_scalar_type , device_type >
      ( comm::size( machine ) , comm::rank( machine ) ,
        global_count_x , global_count_y , global_count_z );

  mesh.parallel_data_map.machine = machine ;

  device_type::fence();
  perf_data.mesh_time = comm::max( machine , wall_clock.seconds() );

  const size_t element_count = mesh.elem_node_ids.dimension(0);

  //------------------------------------
  // Generate sparse matrix graph and element->graph map.

  wall_clock.reset();

  graph_factory::create( mesh , jacobian.graph , element_map );

  device_type::fence();
  perf_data.graph_time = comm::max( machine , wall_clock.seconds() );

  //------------------------------------
  // Allocate linear system coefficients and rhs:

  const size_t local_owned_length = jacobian.graph.row_map.length();

  jacobian.coefficients =
    Kokkos::create_multivector< matrix_coefficients_type >( jacobian.graph.entries.dimension(0) );

  residual =
    Kokkos::create_multivector< vector_type >( local_owned_length );
  delta =
    Kokkos::create_multivector< vector_type >( local_owned_length );
  nodal_solution = 
    Kokkos::create_multivector< vector_type >( local_owned_length );

  //------------------------------------
  // Allocation of arrays to fill the linear system

  typedef Kokkos::Array< scalar_type[ElementNodeCount][ElementNodeCount] , device_type > elem_matrices_type ;
  typedef Kokkos::Array< scalar_type[ElementNodeCount]                   , device_type >  elem_vectors_type ;

  elem_matrices_type elem_matrices ; // Jacobian matrices
  elem_vectors_type  elem_vectors ;  // Residual vectors

  if ( element_count ) {
    elem_matrices =
      Kokkos::create_array< elem_matrices_type >( element_count );

    elem_vectors =
      Kokkos::create_array< elem_vectors_type >( element_count );
  }

  // For boundary condition set the correct values in the solution vector
  // One face of the rectangular mesh is held at 0 potential, and its opposite
  // face is 1.

  DirichletSolutionFunctor::apply(nodal_solution ,
                                  mesh ,
                                  0 , global_count_z - 1 ,
                                  0 , 1 );
    
  do { // Nonlinear loop

    //------------------------------------
    // Compute element matrices and vectors:

    wall_clock.reset();

    ElementFunctor::apply( mesh ,
                           elem_matrices , 
                           elem_vectors ,
                           nodal_solution ,
                           elem_coeff_K , 
                           elem_coeff_P );

    device_type::fence();
    perf_data.elem_time += comm::max( machine , wall_clock.seconds() );

    //------------------------------------
    // Fill linear system coefficients:

    wall_clock.reset();

    fill( 0 , jacobian.coefficients );
    fill( 0 , residual );

    GatherFill< matrix_type , mesh_type >::apply( jacobian , 
                                                  residual ,
                                                  mesh , 
                                                  element_map , 
                                                  elem_matrices , 
                                                  elem_vectors );

    device_type::fence();
    perf_data.matrix_gather_fill_time += comm::max( machine , wall_clock.seconds() );

    // Apply boundary conditions:

    wall_clock.reset();

    // Updates jacobian matrix to 1 on the diagonal, zero elsewhere,
    // and 0 in the residual due to the solution vector having the correct value
    DirichletResidualFunctor::apply( jacobian,
                                     residual,
                                     mesh , 
                                     0 , global_count_z - 1 );

    device_type::fence();
    perf_data.matrix_boundary_condition_time += comm::max( machine , wall_clock.seconds() );

    //------------------------------------
    // Solve linear sytem

    cgsolve( mesh.parallel_data_map ,
             jacobian , residual , delta ,
             cg_iteration_count ,
             cg_residual_norm ,
             cg_iteration_time ,
             cg_iteration_limit , cg_tolerance ) ;

    perf_data.cg_iteration_time += cg_iteration_time ;
    cg_iteration_count_total += cg_iteration_count ;

    // Update non-linear solution with delta...
    // delta is : - Dx = [Jacobian]^1 * Residual which is the negative update
    // LaTeX:
    // \vec {x}_{n+1} = \vec {x}_{n} - ( - \Delta \vec{x}_{n} )
    // text:
    // x[n+1] = x[n] + Dx

    waxpby( mesh.parallel_data_map, 1.0, nodal_solution, -1.0, delta, nodal_solution);

    // Convergence criteria:
    // if the Dx norm is computed to be less than some epsilon, 
    // or if the total number of iterations exceeds some N 
		
    delta_norm = sqrt( dot(mesh.parallel_data_map, delta) );

    ++newton_iteration_count ;

    if ( newton_iteration_count >= newton_iteration_limit ||
         delta_norm < newton_tolerance ) {
	not_converged = false;
    }
  } while( not_converged );

  if ( newton_iteration_count ) {
    perf_data.elem_time /= newton_iteration_count ;
    perf_data.matrix_gather_fill_time /= newton_iteration_count ;
    perf_data.matrix_boundary_condition_time /= newton_iteration_count ;
  }

  if ( cg_iteration_count_total ) {
    perf_data.cg_iteration_time /= cg_iteration_count_total ;
  }

  //------------------------------------

  if ( print_sample ) {

    typename mesh_type::node_coords_type::HostMirror coords_h =
      Kokkos::create_mirror( mesh.node_coords );

    typename vector_type::HostMirror X_h =
      Kokkos::create_mirror( nodal_solution );

    Kokkos::deep_copy( coords_h , mesh.node_coords );
    Kokkos::deep_copy( X_h , nodal_solution );

    for ( size_t i = 0 ; i < mesh.parallel_data_map.count_owned ; ++i ) {
      const coordinate_scalar_type x = coords_h(i,0);
      const coordinate_scalar_type y = coords_h(i,1);
      const coordinate_scalar_type z = coords_h(i,2);

      if ( x <= 0 && y <= 0 ) {
        std::cout << "  node( " << x << " " << y << " " << z << " ) =\t"
                  << X_h(i) << std::endl ;
      }
    }
  }

  return perf_data ;
}

//----------------------------------------------------------------------------

template< typename Scalar , class Device >
void driver( const char * label ,
             comm::Machine machine , int beg , int end , int runs )
{
  if ( beg == 0 || end == 0 || runs == 0 ) return ;

  if ( comm::rank( machine ) == 0 ) {
    std::cout << std::endl ;
    std::cout << "\"Kokkos::HybridFE::NonLinear " << label << "\"" << std::endl;
    std::cout << "\"Size\" ,  \"Meshing\" ,  \"Graphing\" , \"Element\" , \"Fill\" ,   \"Boundary\" ,  \"CG-Iter\"" << std::endl
              << "\"nodes\" , \"millisec\" , \"millisec\" , \"millisec\" , \"millisec\" , \"millisec\" , \"millisec\"" << std::endl ;
  }

  for(int i = beg ; i < end ; i *= 2 )
  {
    const int ix = (int) cbrt( (double) i );
    const int iy = ix + 1 ;
    const int iz = iy + 1 ;
    const int n  = ix * iy * iz ;

    PerformanceData perf_data , perf_best ;

    for(int j = 0; j < runs; j++){

     perf_data = run<Scalar,Device>(machine,ix,iy,iz, false );

     if( j == 0 ) {
       perf_best = perf_data ;
     }
     else {
       perf_best.best( perf_data );
     }
   }

   /// TODO: reduction across processors

  if ( comm::rank( machine ) == 0 ) {

     std::cout << std::setw(8) << n << " , "
               << std::setw(10) << perf_best.mesh_time * 1000 << " , "
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

} /* namespace NonLinear */
} /* namespace HybridFEM */


#endif /* #ifndef HYBRIDFEM_IMPLICIT_HPP */

