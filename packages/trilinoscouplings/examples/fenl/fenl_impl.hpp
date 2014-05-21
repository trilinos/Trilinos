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

// Kokkos libraries' headers:

#include <Kokkos_UnorderedMap.hpp>
#include <Kokkos_StaticCrsGraph.hpp>
#include <Kokkos_CrsMatrix.hpp>
#include <impl/Kokkos_Timer.hpp>
#include <Kokkos_ArithTraits.hpp>

#include <Teuchos_CommHelpers.hpp>

// Examples headers:

#include <BoxElemFixture.hpp>
#include <CGSolve.hpp>
#include <BelosSolve.hpp>

#include <fenl.hpp>
#include <fenl_functors.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {
namespace FENL {

inline
double maximum( const Teuchos::RCP<const Teuchos::Comm<int> >& comm , double local )
{
  double global = 0 ;
  Teuchos::reduceAll( *comm , Teuchos::REDUCE_MAX , 1 , & local , & global );
  return global ;
}

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos */

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {
namespace FENL {

/* Builds a map from LIDs to GIDs suitable for Tpetra */
template < class Map, class Fixture >
class BuildLocalToGlobalMap {
  const Map m_lid_to_gid;
  const Fixture m_fixture;
public:
  typedef typename Map::device_type device_type;
  typedef typename device_type::size_type size_type;

  BuildLocalToGlobalMap(const Map& lid_to_gid, const Fixture& fixture) :
    m_lid_to_gid(lid_to_gid),
    m_fixture(fixture)
    {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const {
    m_lid_to_gid(i) = m_fixture.node_global_index(i);
  }
};

template <class Map, class Fixture >
void build_lid_to_gid(const Map& lid_to_gid, const Fixture& fixture) {
  typedef BuildLocalToGlobalMap<Map,Fixture> F;
  Kokkos::parallel_for(lid_to_gid.dimension_0(), F(lid_to_gid, fixture));
}

template < class Scalar, class Device , BoxElemPart::ElemOrder ElemOrder,
           class CoeffFunctionType , class ManufacturedSolutionType >
Perf fenl(
  const Teuchos::RCP<const Teuchos::Comm<int> >& comm ,
  const Teuchos::RCP<  typename ::Kokkos::Compat::KokkosDeviceWrapperNode<Device> >& node,
  const int use_print ,
  const int use_trials ,
  const int use_atomic ,
  const int use_belos ,
  const int use_muelu ,
  const int use_nodes[] ,
  const CoeffFunctionType& coeff_function ,
  const ManufacturedSolutionType& manufactured_solution ,
  const double bc_lower_value ,
  const double bc_upper_value ,
  const bool check_solution ,
  Scalar& response )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::arrayView;
  using Teuchos::ParameterList;

  typedef Kokkos::Details::ArithTraits<Scalar> KAT;
  typedef typename KAT::mag_type Magnitude;
  typedef Kokkos::Compat::KokkosDeviceWrapperNode<Device> NodeType;
  typedef Tpetra::CrsMatrix<Scalar,int,int,NodeType> GlobalMatrixType;
  typedef Tpetra::Vector<Scalar,int,int,NodeType> GlobalVectorType;
  typedef Tpetra::Map<int, int, NodeType> MapType;
  typedef RCP<const MapType> pMapType;

  typedef Kokkos::Example::BoxElemFixture< Device , ElemOrder > FixtureType ;

  //typedef Kokkos::CrsMatrix< double , unsigned , Device >
  typedef typename GlobalMatrixType::k_local_matrix_type
    LocalMatrixType ;

  typedef typename LocalMatrixType::StaticCrsGraphType
    LocalGraphType ;

  typedef Kokkos::Example::FENL::NodeNodeGraph< typename FixtureType::elem_node_type , LocalGraphType , FixtureType::ElemNode >
     NodeNodeGraphType ;

  typedef Kokkos::Example::FENL::ElementComputation< FixtureType , LocalMatrixType , CoeffFunctionType >
    ElementComputationType ;

  typedef Kokkos::Example::FENL::DirichletComputation< FixtureType , LocalMatrixType >
    DirichletComputationType ;

  typedef NodeElemGatherFill< ElementComputationType >
    NodeElemGatherFillType ;

  typedef typename ElementComputationType::vector_type LocalVectorType ;
  typedef Kokkos::DualView< Scalar** , Kokkos::LayoutLeft, Device > LocalDualVectorType;

  typedef Kokkos::Example::FENL::ResponseComputation< FixtureType , LocalVectorType >
    ResponseComputationType ;

  //------------------------------------

  const unsigned  newton_iteration_limit     = 10 ;
  const Magnitude newton_iteration_tolerance = 1e-7 ;
  const unsigned  cg_iteration_limit         = 200 ;
  const Magnitude cg_iteration_tolerance     = 1e-7 ;

  //------------------------------------

  const int print_flag = use_print && Kokkos::Impl::is_same< Kokkos::HostSpace , typename Device::memory_space >::value ;

  const int comm_rank = comm->getRank();
  const int comm_size = comm->getSize();

  // Decompose by node to avoid parallel communication in assembly

  const float bubble_x = 1.0 ;
  const float bubble_y = 1.0 ;
  const float bubble_z = 1.0 ;

  const FixtureType fixture( BoxElemPart::DecomposeNode ,
                             comm_size , comm_rank ,
                             use_nodes[0] , use_nodes[1] , use_nodes[2] ,
                             bubble_x , bubble_y , bubble_z );

  //------------------------------------

  //------------------------------------

  if ( print_flag ) {
    manufactured_solution.print( std::cout , fixture );

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

  Kokkos::Impl::Timer wall_clock ;

  Perf perf_stats = Perf() ;

  for ( int itrial = 0 ; itrial < use_trials ; ++itrial ) {

    Perf perf = Perf() ;

    perf.global_elem_count = fixture.elem_count_global();
    perf.global_node_count = fixture.node_count_global();

    //----------------------------------
    // Create the local sparse matrix graph and element-to-graph map
    // from the element->to->node identifier array.
    // The graph only has rows for the owned nodes.

    typename NodeNodeGraphType::Times graph_times;

    const NodeNodeGraphType
      mesh_to_graph( fixture.elem_node() , fixture.node_count_owned(), graph_times );

    perf.map_ratio          = maximum(comm, graph_times.ratio);
    perf.fill_node_set      = maximum(comm, graph_times.fill_node_set);
    perf.scan_node_count    = maximum(comm, graph_times.scan_node_count);
    perf.fill_graph_entries = maximum(comm, graph_times.fill_graph_entries);
    perf.sort_graph_entries = maximum(comm, graph_times.sort_graph_entries);
    perf.fill_element_graph = maximum(comm, graph_times.fill_element_graph);

    wall_clock.reset();
    // Create the local sparse matrix from the graph:

    LocalMatrixType jacobian( "jacobian" , mesh_to_graph.graph );
    // jacobian.dev_config.block_dim.x = 1;
    // jacobian.dev_config.block_dim.y = 1;
    // jacobian.dev_config.num_blocks = jacobian.numRows();

    Device::fence();

    perf.create_sparse_matrix = maximum( comm , wall_clock.seconds() );

    //----------------------------------

    if ( print_flag ) {
      const unsigned nrow = jacobian.numRows();
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

      std::cout << "ElemGraph {" << std::endl ;
      for ( unsigned ielem = 0 ; ielem < mesh_to_graph.elem_graph.dimension_0() ; ++ielem ) {
        std::cout << "  elem[" << ielem << "]{" ;
        for ( unsigned irow = 0 ; irow < mesh_to_graph.elem_graph.dimension_1() ; ++irow ) {
          std::cout << " {" ;
          for ( unsigned icol = 0 ; icol < mesh_to_graph.elem_graph.dimension_2() ; ++icol ) {
            std::cout << " " << mesh_to_graph.elem_graph(ielem,irow,icol);
          }
          std::cout << " }" ;
        }
        std::cout << " }" << std::endl ;
      }
      std::cout << "}" << std::endl ;
    }

    //----------------------------------

    // Allocate solution vector for each node in the mesh and residual vector for each owned node
    // We need dual vectors for Tpetra!!
    const LocalDualVectorType k_nodal_solution(
      "nodal_solution" , fixture.node_count(),1 );
    const LocalDualVectorType k_nodal_residual(
      "nodal_residual" , fixture.node_count_owned(),1 );
    const LocalDualVectorType k_nodal_delta(
      "nodal_delta" ,    fixture.node_count_owned(),1 );
    const LocalVectorType nodal_solution =
      Kokkos::subview<LocalVectorType>(k_nodal_solution.d_view,Kokkos::ALL(),0);
    const LocalVectorType nodal_residual =
      Kokkos::subview<LocalVectorType>(k_nodal_residual.d_view,Kokkos::ALL(),0);
    const LocalVectorType nodal_delta =
      Kokkos::subview<LocalVectorType>(k_nodal_delta.d_view,Kokkos::ALL(),0);

    // Create element computation functor
    const ElementComputationType elemcomp(
      use_atomic ? ElementComputationType( fixture , coeff_function ,
                                           nodal_solution ,
                                           mesh_to_graph.elem_graph ,
                                           jacobian , nodal_residual )
                 : ElementComputationType( fixture , coeff_function ,
                                           nodal_solution ) );

    const NodeElemGatherFillType gatherfill(
      use_atomic ? NodeElemGatherFillType()
                 : NodeElemGatherFillType( fixture.elem_node() ,
                                           mesh_to_graph.elem_graph ,
                                           nodal_residual ,
                                           jacobian ,
                                           elemcomp.elem_residuals ,
                                           elemcomp.elem_jacobians ) );

    // Create boundary condition functor
    const DirichletComputationType dirichlet(
      fixture , nodal_solution , jacobian , nodal_residual ,
      2 /* apply at 'z' ends */ ,
      bc_lower_value ,
      bc_upper_value );

    const ParameterList params();

    // Create Distributed Objects

    // Create Maps
    typedef Kokkos::View<int*,Device> lid_to_gid_type;
    lid_to_gid_type lid_to_gid_row("lig_to_gid",jacobian.numRows());
    build_lid_to_gid(lid_to_gid_row, fixture);
    typename lid_to_gid_type::HostMirror lid_to_gid_row_host =
      Kokkos::create_mirror_view(lid_to_gid_row);
    Kokkos::deep_copy(lid_to_gid_row_host, lid_to_gid_row);

    pMapType RowMap = rcp (new MapType (fixture.node_count_global(),
        arrayView(lid_to_gid_row_host.ptr_on_device(),lid_to_gid_row.dimension_0()),
        0, comm, node));

    lid_to_gid_type lid_to_gid("lig_to_gid",jacobian.numCols());
    build_lid_to_gid(lid_to_gid, fixture);
    typename lid_to_gid_type::HostMirror lid_to_gid_host =
      Kokkos::create_mirror_view(lid_to_gid);
    Kokkos::deep_copy(lid_to_gid_host, lid_to_gid);

    pMapType ColMap = rcp (new MapType (
        Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
        arrayView(lid_to_gid_host.ptr_on_device(),lid_to_gid.dimension_0()),
        0,comm, node) );

    // Create Teptra Matrix: this uses the already allocated matrix data
    GlobalMatrixType g_jacobian(RowMap,ColMap,jacobian);

    // Create Teptra Vectors: this uses the already allocated vector data
    GlobalVectorType g_nodal_solution(ColMap,k_nodal_solution);
    GlobalVectorType g_nodal_residual(RowMap,k_nodal_residual);
    GlobalVectorType g_nodal_delta(RowMap,k_nodal_delta);

    // Create a subview of just the owned data of the solution vector
    LocalDualVectorType k_nodal_solution_no_overlap =
      Kokkos::subview<LocalDualVectorType>(k_nodal_solution,std::pair<unsigned,unsigned>(0,k_nodal_delta.dimension_0()),Kokkos::ALL());
    GlobalVectorType g_nodal_solution_no_overlap(RowMap,
                                                 k_nodal_solution_no_overlap);

    typedef Tpetra::Import<typename GlobalVectorType::local_ordinal_type,
      typename GlobalVectorType::global_ordinal_type,
      typename GlobalVectorType::node_type> import_type;
    import_type import (RowMap, ColMap);

    // Create response function
    LocalVectorType nodal_solution_no_overlap =
      Kokkos::subview<LocalVectorType>(k_nodal_solution_no_overlap.d_view,
                                       Kokkos::ALL(), 0);
    const ResponseComputationType responseFunc(
      fixture, nodal_solution_no_overlap );

    //----------------------------------
    // Nonlinear Newton iteration:

    Magnitude residual_norm_init = 0 ;

    // RCP<Teuchos::FancyOStream> out =
    //   Teuchos::fancyOStream(rcp(&std::cout,false));
    // out->setShowProcRank(true);

    for ( perf.newton_iter_count = 0 ;
          perf.newton_iter_count < newton_iteration_limit ;
          ++perf.newton_iter_count ) {

      //--------------------------------

      g_nodal_solution.doImport (g_nodal_solution_no_overlap, import, Tpetra::REPLACE);

      // if (itrial == 0 && perf.newton_iter_count == 0)
      //   g_nodal_solution_no_overlap.describe(*out, Teuchos::VERB_EXTREME);

      //--------------------------------
      // Element contributions to residual and jacobian

      wall_clock.reset();

      Kokkos::deep_copy( nodal_residual , Scalar(0) );
      Kokkos::deep_copy( jacobian.values , Scalar(0) );

      elemcomp.apply();

      if ( ! use_atomic ) {
        gatherfill.apply();
      }

      Device::fence();
      perf.fill_time = maximum( comm , wall_clock.seconds() );

      //--------------------------------
      // Apply boundary conditions

      wall_clock.reset();

      dirichlet.apply();

      Device::fence();
      perf.bc_time = maximum( comm , wall_clock.seconds() );

      //--------------------------------
      // Evaluate convergence

      const Magnitude residual_norm =
          g_nodal_residual.norm2();

      perf.newton_residual = residual_norm ;

      if ( 0 == perf.newton_iter_count ) { residual_norm_init = residual_norm ; }

      if ( residual_norm < residual_norm_init * newton_iteration_tolerance ) { break ; }

      //--------------------------------
      // Solve for nonlinear update

      result_struct cgsolve;
      if (use_belos) {
        cgsolve = belos_solve(rcpFromRef(g_jacobian),
                              rcpFromRef(g_nodal_residual),
                              rcpFromRef(g_nodal_delta),
                              use_muelu,
                              cg_iteration_limit ,
                              cg_iteration_tolerance);
      }
      else {
        cgsolve = cg_solve(rcpFromRef(g_jacobian),
                           rcpFromRef(g_nodal_residual),
                           rcpFromRef(g_nodal_delta),
                           cg_iteration_limit,
                           cg_iteration_tolerance);
      }

      // Update solution vector

      g_nodal_solution_no_overlap.update(-1.0,g_nodal_delta,1.0);
      perf.cg_iter_count += cgsolve.iteration ;
      perf.cg_time       += cgsolve.iter_time ;

      //--------------------------------

      if ( print_flag ) {
        const double delta_norm =
            g_nodal_delta.norm2();

        std::cout << "Newton iteration[" << perf.newton_iter_count << "]"
                  << " residual[" << perf.newton_residual << "]"
                  << " update[" << delta_norm << "]"
                  << " cg_iteration[" << cgsolve.iteration << "]"
                  << " cg_residual[" << cgsolve.norm_res << "]"
                  << std::endl ;

        const unsigned nrow = jacobian.numRows();

        std::cout << "Residual {" ;
        for ( unsigned irow = 0 ; irow < nrow ; ++irow ) {
          std::cout << " " << nodal_residual(irow);
        }
        std::cout << " }" << std::endl ;

        std::cout << "Delta {" ;
        for ( unsigned irow = 0 ; irow < nrow ; ++irow ) {
          std::cout << " " << nodal_delta(irow);
        }
        std::cout << " }" << std::endl ;

        std::cout << "Solution {" ;
        for ( unsigned irow = 0 ; irow < nrow ; ++irow ) {
          std::cout << " " << nodal_solution(irow);
        }
        std::cout << " }" << std::endl ;

        std::cout << "Jacobian[ "
                  << jacobian.numRows() << " x " << jacobian.numCols()
                  << " ] {" << std::endl ;
        for ( unsigned irow = 0 ; irow < nrow ; ++irow ) {
          std::cout << "  {" ;
          const unsigned entry_end = jacobian.graph.row_map(irow+1);
          for ( unsigned entry = jacobian.graph.row_map(irow) ; entry < entry_end ; ++entry ) {
            std::cout << " (" << jacobian.graph.entries(entry)
                      << "," << jacobian.values(entry)
                      << ")" ;
          }
          std::cout << " }" << std::endl ;
        }
        std::cout << "}" << std::endl ;
      }

      //--------------------------------
    }

    // Evaluate response function -- currently 2-norm of solution vector

    response = responseFunc.apply();
    response = Kokkos::Example::all_reduce( response , comm );

    // Evaluate solution error

    if ( 0 == itrial ) {
      if ( check_solution ) {
        const double error_max =
          manufactured_solution.compute_error( fixture, nodal_solution );
        perf.error_max =
          std::sqrt( Kokkos::Example::all_reduce_max( error_max , comm ) );
      }

      perf_stats = perf ;
    }
    else {
      perf_stats.fill_node_set =
        std::min( perf_stats.fill_node_set , perf.fill_node_set );
      perf_stats.scan_node_count =
        std::min( perf_stats.scan_node_count , perf.scan_node_count );
      perf_stats.fill_graph_entries =
        std::min( perf_stats.fill_graph_entries , perf.fill_graph_entries );
      perf_stats.sort_graph_entries =
        std::min( perf_stats.sort_graph_entries , perf.sort_graph_entries );
      perf_stats.fill_element_graph =
        std::min( perf_stats.fill_element_graph , perf.fill_element_graph );
      perf_stats.create_sparse_matrix =
        std::min( perf_stats.create_sparse_matrix , perf.create_sparse_matrix );
      perf_stats.fill_time =
        std::min( perf_stats.fill_time , perf.fill_time );
      perf_stats.bc_time =
        std::min( perf_stats.bc_time , perf.bc_time );
      perf_stats.cg_time =
        std::min( perf_stats.cg_time , perf.cg_time );
    }
  }

  return perf_stats ;
}

#define INST_FENL( SCALAR, DEVICE, ELEMENT, COEFF, MS )              \
  template Perf                                                      \
  fenl< SCALAR, DEVICE , ELEMENT , COEFF, MS >(                      \
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm ,            \
    const Teuchos::RCP<Kokkos::Compat::KokkosDeviceWrapperNode<DEVICE> >& node,\
    const int use_print ,                                            \
    const int use_trials ,                                           \
    const int use_atomic ,                                           \
    const int use_belos ,                                            \
    const int use_muelu ,                                            \
    const int global_elems[] ,                                       \
    const COEFF& coeff_function ,                                    \
    const MS& manufactured_solution ,                                \
    const double bc_lower_value ,                                    \
    const double bc_upper_value ,                                    \
    const bool check_solution ,                                      \
    SCALAR& response);

//----------------------------------------------------------------------------

template < typename Scalar, typename MeshScalar, typename Device >
ElementComputationKLCoefficient<Scalar,MeshScalar,Device>::
ElementComputationKLCoefficient( const MeshScalar mean ,
                                 const MeshScalar variance ,
                                 const MeshScalar correlation_length ,
                                 const size_type num_rv )
  : m_mean( mean ),
    m_variance( variance ),
    m_corr_len( correlation_length ),
    m_num_rv( num_rv ),
    m_rv( "KL Random Variables", m_num_rv ),
    m_eig( "KL Eigenvalues", m_num_rv ),
    m_pi( 4.0*std::atan(1.0) )
{
  typename EigenView::HostMirror host_eig =
    Kokkos::create_mirror_view( m_eig );

  const MeshScalar a = std::sqrt( std::sqrt(m_pi)*m_corr_len );

  if (m_num_rv > 0)
    host_eig(0) = a / std::sqrt( MeshScalar(2) );

  for ( size_type i=1; i<m_num_rv; ++i ) {
    const MeshScalar b = (i+1)/2;  // floor((i+1)/2)
    const MeshScalar c = b * m_pi * m_corr_len;
    host_eig(i) = a * std::exp( -c*c / MeshScalar(8) );
  }

  Kokkos::deep_copy( m_eig , host_eig );
}

#define INST_KL( SCALAR, MESH_SCALAR, DEVICE )                        \
  template class                                                      \
  ElementComputationKLCoefficient< SCALAR, MESH_SCALAR, DEVICE >;

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos */

#endif /* #ifndef KOKKOS_EXAMPLE_FENL_IMPL_HPP */
