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

#ifndef KOKKOS_EXAMPLE_FENL_IMPL_HPP
#define KOKKOS_EXAMPLE_FENL_IMPL_HPP

#include <math.h>
#include <Kokkos_UnorderedMap.hpp>
#include <Kokkos_StaticCrsGraph.hpp>
#include <Kokkos_CrsMatrix.hpp>

#include <BoxElemFixture.hpp>
#include <fenl.hpp>
#include <fenlFunctors.hpp>
#include <impl/Kokkos_Timer.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {
namespace FENL {

inline
double maximum( MPI_Comm comm , double local )
{
  double global = local ;
#if defined( KOKKOS_HAVE_MPI )
  MPI_Allreduce( & local , & global , 1 , MPI_DOUBLE , MPI_MAX , comm );
#endif
  return global ;
}

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos */

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {
namespace FENL {

class ManufacturedSolution {
public:

  // Manufactured solution for one dimensional nonlinear PDE
  //
  //  -K T_zz + T^2 = 0 ; T(zmin) = T_zmin ; T(zmax) = T_zmax
  //
  //  Has an analytic solution of the form:
  //
  //    T(z) = ( a ( z - zmin ) + b )^(-2) where K = 1 / ( 6 a^2 )
  //
  //  Given T_0 and T_L compute K for this analytic solution.
  //
  //  Two analytic solutions:
  //
  //    Solution with singularity:
  //    , a( ( 1.0 / sqrt(T_zmax) + 1.0 / sqrt(T_zmin) ) / ( zmax - zmin ) )
  //    , b( -1.0 / sqrt(T_zmin) )
  //
  //    Solution without singularity:
  //    , a( ( 1.0 / sqrt(T_zmax) - 1.0 / sqrt(T_zmin) ) / ( zmax - zmin ) )
  //    , b( 1.0 / sqrt(T_zmin) )

  const double zmin ;
  const double zmax ;
  const double T_zmin ;
  const double T_zmax ;
  const double a ;
  const double b ;
  const double K ;

  ManufacturedSolution( const double arg_zmin ,
                        const double arg_zmax ,
                        const double arg_T_zmin ,
                        const double arg_T_zmax )
    : zmin( arg_zmin )
    , zmax( arg_zmax )
    , T_zmin( arg_T_zmin )
    , T_zmax( arg_T_zmax )
    , a( ( 1.0 / sqrt(T_zmax) - 1.0 / sqrt(T_zmin) ) / ( zmax - zmin ) )
    , b( 1.0 / sqrt(T_zmin) )
    , K( 1.0 / ( 6.0 * a * a ) )
    {}

  double operator()( const double z ) const
  {
    const double tmp = a * ( z - zmin ) + b ;
    return 1.0 / ( tmp * tmp );
  }
};

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos */

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {
namespace FENL {

template < class Device , BoxElemPart::ElemOrder ElemOrder >
Perf fenl(
  MPI_Comm comm ,
  const int use_print ,
  const int use_trials ,
  const int use_atomic ,
  const int use_nodes[] )
{
  typedef Kokkos::Example::BoxElemFixture< Device , ElemOrder > FixtureType ;

  typedef Kokkos::CrsMatrix< double , unsigned , Device , void , unsigned >
    SparseMatrixType ;

  typedef typename SparseMatrixType::StaticCrsGraphType
    SparseGraphType ;

  typedef Kokkos::Example::FENL::NodeNodeGraph< typename FixtureType::elem_node_type , SparseGraphType , FixtureType::ElemNode > 
     NodeNodeGraphType ;

#if 0
  typedef Kokkos::Example::FENL::ElementComputation< FixtureType , SparseMatrixType >
    ElementComputation ;

  typedef Kokkos::Example::FENL::DirichletComputation< FixtureType , SparseMatrixType >
    DirichletComputation ;
#endif

  //------------------------------------

  const int print_flag = use_print && Kokkos::Impl::is_same< Kokkos::HostSpace , typename Device::memory_space >::value ;

  int comm_rank ;
  int comm_size ;

  MPI_Comm_rank( comm , & comm_rank );
  MPI_Comm_size( comm , & comm_size );

  // Decompose by node to avoid mpi-communication for assembly

  const FixtureType fixture( BoxElemPart::DecomposeNode , comm_size , comm_rank ,
                             use_nodes[0] , use_nodes[1] , use_nodes[2] );

  if ( print_flag ) {
    std::cout << "ElemNode {" << std::endl ;
    for ( unsigned ielem = 0 ; ielem < fixture.elem_count() ; ++ielem ) {
      std::cout << "  elem[" << ielem << "]{" ;
      for ( unsigned inode = 0 ; inode < FixtureType::ElemNode ; ++inode ) {
        std::cout << " " << fixture.elem_node(ielem,inode);
      }
      std::cout << " }" << std::endl ;
    }
    std::cout << "}" << std::endl ;
  }

  //------------------------------------

  const double bc_lower_value = 1 ;
  const double bc_upper_value = 2 ;

  const Kokkos::Example::FENL::ManufacturedSolution
    manufactured_solution( 0 , 1 , bc_lower_value , bc_upper_value  );

  //------------------------------------

  std::vector< Kokkos::Example::FENL::Perf > perf_all( use_trials );

  Kokkos::Impl::Timer wall_clock ;

  for ( int i = 0 ; i < use_trials ; ++i ) {

    Kokkos::Example::FENL::Perf & perf = perf_all[i] ;

    perf.global_elem_count = fixture.elem_count_global();
    perf.global_node_count = fixture.node_count_global();

    //----------------------------------
    // Create the sparse matrix graph and element-to-graph map
    // from the element->to->node identifier array.
    // The graph only has rows for the owned nodes.

    wall_clock.reset();

    const NodeNodeGraphType
      mesh_to_graph( fixture.elem_node() , fixture.node_count_owned() );

    // Create the sparse matrix from the graph:

    SparseMatrixType jacobian( "jacobian" , mesh_to_graph.graph );

    Device::fence();

    perf.graph_time = maximum( comm , wall_clock.seconds() );

    //----------------------------------

    if ( print_flag ) {
      const unsigned nrow = mesh_to_graph.graph.row_map.dimension_0() - 1 ;
      std::cout << "JacobianGraph[ "
                << jacobian.numRows() << " x " << jacobian.numCols()
                << " ] {" << std::endl ;
      for ( unsigned irow = 0 ; irow < nrow ; ++irow ) {
        std::cout << "  row[" << irow << "]{" ;
        const unsigned entry_end = jacobian.graph.row_map(irow+1);
        for ( unsigned entry = jacobian.graph.row_map(irow) ; entry < entry_end ; ++entry ) {
          std::cout << " " << jacobian.graph.entries(entry);
        }
        std::cout << " }" << std::endl ;
      }
      std::cout << "}" << std::endl ;
    }

    //----------------------------------

#if 0

    const ElementComputation    elemcomp( fixture , use_atomic , coeff_K );
    const DirichletComputation  dirichlet( fixture );

    //----------------------------------
    // Set the solution vector

    dirichlet.apply_solution( nodal_solution );

    //----------------------------------
    // Nonlinear Newton iteration:

    for ( perf.newton_iteration = 0 ;
          perf.newton_iteration < newton_iteration_limit ;
          ++perf.newton_iteration ) {

      //--------------------------------

      comm_import.apply( nodal_solution );

      //--------------------------------
      // Element contributions to residual and jacobian

      if ( use_atomic ) {
        wall_clock.reset();
        elemcomp.apply_and_atomic_fill( jacobian , residual , nodal_solution );
        Device::fence();
        perf_data.elem_time = maximum( comm , wall_clock.seconds() );
      }
      else {
        wall_clock.reset();
        elemcomp.apply( nodal_solution );
        Device::fence();
        perf_data.elem_time = maximum( comm , wall_clock.seconds() );

        wall_clock.reset();
        gatherfill.apply( jacobian , residual , elem_jacobian , elem_residual );
        Device::fence();
        perf_data.fill_time = maximum( comm , wall_clock.seconds() );
      }

      //--------------------------------
      // Apply boundary conditions

      dirichlet.apply_residual( jacobian , residual );

      //--------------------------------
      // Evaluate convergence

      //--------------------------------
      // Solve for nonlinear update

      cgsolve( ... );

      // Update solution vector

    }

    // Per-iteration times

#endif

  }

  // Evaluate solution error

  // Performance statistics

  Perf perf_stats = perf_all[0];

  return perf_stats ;
}

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos */

#endif /* #ifndef KOKKOS_EXAMPLE_FENL_IMPL_HPP */

