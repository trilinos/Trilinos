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
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

// Examples headers:

#include <BoxElemFixture.hpp>
#include <CGSolve.hpp>
#include <BelosSolve.hpp>

#include <fenl.hpp>
#include <fenl_functors.hpp>
#include <Kokkos_DefaultNode.hpp>

#include <Tpetra_Vector.hpp>
#include "Tpetra_MultiVector.hpp"


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
  typedef typename Map::execution_space execution_space;
  typedef typename execution_space::size_type size_type;

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

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos */

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {
namespace FENL {

template < class Scalar, class Device , BoxElemPart::ElemOrder ElemOrder, class CoeffFunctionType >
class Problem {
public:


  typedef BoxElemFixture< Device , ElemOrder >  FixtureType ;

  typedef typename Kokkos::Details::ArithTraits<Scalar>::mag_type  Magnitude;

  typedef Kokkos::Compat::KokkosDeviceWrapperNode< Device >  NodeType;

  typedef Kokkos::View<     Scalar * , Kokkos::LayoutLeft, Device >  LocalVectorType ;
  typedef Kokkos::DualView< Scalar** , Kokkos::LayoutLeft, Device >  LocalDualVectorType;

  typedef Tpetra::Map<int, int, NodeType>              MapType;
  typedef Tpetra::Vector<Scalar,int,int,NodeType>      GlobalVectorType;
  typedef Tpetra::CrsMatrix<Scalar,int,int,NodeType>   GlobalMatrixType;
  typedef typename GlobalMatrixType::local_matrix_type LocalMatrixType;
  typedef typename LocalMatrixType::StaticCrsGraphType LocalGraphType;


  typedef NodeNodeGraph< typename FixtureType::elem_node_type
                       , LocalGraphType
                       , FixtureType::ElemNode
                       > NodeNodeGraphType ;

  typedef typename NodeNodeGraphType::ElemGraphType  ElemGraphType ;

  typedef Teuchos::RCP<const MapType>                        rcpMapType;
  typedef Teuchos::RCP<const Teuchos::Comm<int> >            rcpCommType ;
  typedef Teuchos::RCP<Kokkos::Compat::KokkosDeviceWrapperNode<Device> > rcpNodeType ;

  typedef Tpetra::Import<
    typename GlobalVectorType::local_ordinal_type,
    typename GlobalVectorType::global_ordinal_type,
    typename GlobalVectorType::node_type
    > import_type;


private:

  typedef Kokkos::View<int*,Device> lid_to_gid_type;

  rcpMapType create_row_map()
  {
    lid_to_gid_type lid_to_gid_row( "lig_to_gid_row", fixture.node_count_owned() );

    build_lid_to_gid(lid_to_gid_row, fixture);

    typename lid_to_gid_type::HostMirror lid_to_gid_row_host =
      Kokkos::create_mirror_view(lid_to_gid_row);

    Kokkos::deep_copy(lid_to_gid_row_host, lid_to_gid_row);

    return  Teuchos::rcp (new MapType( fixture.node_count_global(),
        Teuchos::arrayView(lid_to_gid_row_host.ptr_on_device(),
                           lid_to_gid_row_host.dimension_0()),
        0, comm, node));
  }

  rcpMapType create_col_map()
  {
    lid_to_gid_type lid_to_gid_all( "lig_to_gid_all", fixture.node_count());

    build_lid_to_gid(lid_to_gid_all, fixture);

    typename lid_to_gid_type::HostMirror lid_to_gid_all_host =
      Kokkos::create_mirror_view(lid_to_gid_all);

    Kokkos::deep_copy(lid_to_gid_all_host, lid_to_gid_all);

    return Teuchos::rcp (new MapType (
        Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
        Teuchos::arrayView( lid_to_gid_all_host.ptr_on_device(),
                            lid_to_gid_all_host.dimension_0()),
        0,comm, node) );
  }

public:

  const rcpCommType  comm ;
  const rcpNodeType  node ;
  const FixtureType  fixture ;

private:

  typename NodeNodeGraphType::Times graph_times;

  const NodeNodeGraphType mesh_to_graph ;

public:

  const ElemGraphType elem_graph ;

  const rcpMapType   RowMap ;
  const rcpMapType   ColMap ;
  const import_type  import ;

  GlobalVectorType   g_nodal_solution ;
  GlobalVectorType   g_nodal_residual ;
  GlobalVectorType   g_nodal_delta ;
  GlobalVectorType   g_nodal_solution_no_overlap ;
  GlobalMatrixType   g_jacobian ;

  Scalar             response ;
  Perf               perf ;

  Problem( const rcpCommType & use_comm
         , const rcpNodeType & use_node
         , const int use_nodes[]
         , const double grid_bubble[]
         , const bool print_flag
         )
    : comm( use_comm )
    , node( use_node )
    // Decompose by node to avoid parallel communication in assembly
    , fixture( BoxElemPart::DecomposeNode
             , use_comm->getSize() , use_comm->getRank()
             , use_nodes[0] , use_nodes[1] , use_nodes[2]
             , grid_bubble[0] , grid_bubble[1] , grid_bubble[2]
             )
    , graph_times()
    , mesh_to_graph( fixture.elem_node()
                   , fixture.node_count_owned()
                   , graph_times )
    , elem_graph( mesh_to_graph.elem_graph )
    , RowMap( create_row_map() )
    , ColMap( create_col_map() )
    , import( RowMap , ColMap )
    , g_nodal_solution( ColMap, 1 )
    , g_nodal_residual( RowMap, 1 )
    , g_nodal_delta(    RowMap, 1 )
    , g_nodal_solution_no_overlap(
        RowMap ,
        Kokkos::subview( g_nodal_solution.getDualView()
                                            , std::pair<unsigned,unsigned>(0,fixture.node_count_owned())
                                            , Kokkos::ALL()
                                            ) )
    , g_jacobian( RowMap, ColMap, LocalMatrixType( "jacobian" , mesh_to_graph.graph ) )
    , response()
    , perf()
    {
      if ( maximum(comm, ( fixture.ok() ? 0 : 1 ) ) ) {
        throw std::runtime_error(std::string("Problem fixture setup failed"));
      }

      perf.global_elem_count  = fixture.elem_count_global();
      perf.global_node_count  = fixture.node_count_global();

      perf.map_ratio          = maximum(comm, graph_times.ratio);
      perf.fill_node_set      = maximum(comm, graph_times.fill_node_set);
      perf.scan_node_count    = maximum(comm, graph_times.scan_node_count);
      perf.fill_graph_entries = maximum(comm, graph_times.fill_graph_entries);
      perf.sort_graph_entries = maximum(comm, graph_times.sort_graph_entries);
      perf.fill_element_graph = maximum(comm, graph_times.fill_element_graph);

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

        //------------------------------

        LocalMatrixType jacobian = g_jacobian.getLocalMatrix();

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
        for ( unsigned ielem = 0 ; ielem < elem_graph.dimension_0() ; ++ielem ) {
          std::cout << "  elem[" << ielem << "]{" ;
          for ( unsigned irow = 0 ; irow < elem_graph.dimension_1() ; ++irow ) {
            std::cout << " {" ;
            for ( unsigned icol = 0 ; icol < elem_graph.dimension_2() ; ++icol ) {
              std::cout << " " << elem_graph(ielem,irow,icol);
            }
            std::cout << " }" ;
          }
          std::cout << " }" << std::endl ;
        }
        std::cout << "}" << std::endl ;
      }
    }

  //----------------------------------------

  void solve( const CoeffFunctionType & coeff_function
            , const double coeff_source
            , const double coeff_advection
            , const double    bc_lower_value
            , const double    bc_upper_value
            , const unsigned  newton_iteration_limit
            , const Magnitude newton_iteration_tolerance
            , const unsigned  cg_iteration_limit
            , const Magnitude cg_iteration_tolerance
            , const QuadratureData<Device>& qd
            , const bool   use_atomic
            , const bool   use_belos
            , const bool   use_muelu
            , const bool   use_mean_based
            , const bool   print_flag
            , const Teuchos::RCP<Teuchos::ParameterList>& fenlParams
            )
    {
      typedef ElementComputation< FixtureType , LocalMatrixType , CoeffFunctionType >
        ElementComputationType ;

      typedef DirichletComputation< FixtureType , LocalMatrixType >
        DirichletComputationType ;

      typedef ResponseComputation< FixtureType , LocalVectorType >
        ResponseComputationType ;

      Kokkos::Impl::Timer wall_clock ;

      LocalMatrixType jacobian = g_jacobian.getLocalMatrix();

      // Extract DualViews
      const LocalDualVectorType k_nodal_solution = g_nodal_solution.getDualView();
      const LocalDualVectorType k_nodal_residual = g_nodal_residual.getDualView();
      const LocalDualVectorType k_nodal_delta    = g_nodal_delta   .getDualView();

      const LocalVectorType nodal_solution =
        Kokkos::subview(k_nodal_solution.d_view,Kokkos::ALL(),0);
      const LocalVectorType nodal_residual =
        Kokkos::subview(k_nodal_residual.d_view,Kokkos::ALL(),0);
      const LocalVectorType nodal_delta =
        Kokkos::subview(k_nodal_delta.d_view,Kokkos::ALL(),0);

      LocalVectorType nodal_solution_no_overlap =
        Kokkos::subview(nodal_solution,std::pair<unsigned,unsigned>(0,fixture.node_count_owned()));

      // Get DeviceConfig structs used by some functors
      Kokkos::DeviceConfig dev_config_elem, dev_config_gath, dev_config_bc;

      CreateDeviceConfigs<Scalar>::eval( dev_config_elem,
                                         dev_config_gath,
                                         dev_config_bc );

      // Create element computation functor
      const ElementComputationType elemcomp( fixture , coeff_function ,
                                             coeff_source , coeff_advection ,
                                             nodal_solution ,
                                             elem_graph ,
                                             jacobian , nodal_residual ,
                                             dev_config_elem , qd );

      // Create boundary condition functor
      const DirichletComputationType dirichlet(
        fixture , nodal_solution , jacobian , nodal_residual ,
        2 /* apply at 'z' ends */ ,
        bc_lower_value ,
        bc_upper_value ,
        dev_config_bc );

      // Create response function
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

        wall_clock.reset();
        g_nodal_solution.doImport (g_nodal_solution_no_overlap, import, Tpetra::REPLACE);
        Device::fence();
        perf.import_time = maximum( comm , wall_clock.seconds() );

        // if (itrial == 0 && perf.newton_iter_count == 0)
        //   g_nodal_solution_no_overlap.describe(*out, Teuchos::VERB_EXTREME);

        //--------------------------------
        // Element contributions to residual and jacobian

        wall_clock.reset();

        Kokkos::deep_copy( nodal_residual , 0.0 );
        Kokkos::deep_copy( jacobian.values ,0.0 );

        elemcomp.apply();

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

        const Magnitude residual_norm = g_nodal_residual.norm2();

        perf.newton_residual = residual_norm ;
        perf.error_max = residual_norm ;

        if ( 0 == perf.newton_iter_count ) {
          residual_norm_init = residual_norm ;
        }

        if ( residual_norm < residual_norm_init * newton_iteration_tolerance ) {
          break ;
        }

        //--------------------------------
        // Solve for nonlinear update

        result_struct cgsolve;
        if (use_belos) {
          // Don't accumulate Belos times as the internal Teuchos timers
          // already accumulate
          cgsolve = belos_solve(rcpFromRef(g_jacobian),
                                rcpFromRef(g_nodal_residual),
                                rcpFromRef(g_nodal_delta),
                                fixture,
                                use_muelu,
                                use_mean_based ,
                                fenlParams ,
                                cg_iteration_limit ,
                                cg_iteration_tolerance);
          perf.mat_vec_time    = cgsolve.matvec_time ;
          perf.cg_iter_time    = cgsolve.iter_time ;
          perf.prec_setup_time = cgsolve.prec_setup_time ;
          perf.prec_apply_time = cgsolve.prec_apply_time ;
          perf.cg_total_time   = cgsolve.total_time ;
        }
        else {
          cgsolve = cg_solve(rcpFromRef(g_jacobian),
                             rcpFromRef(g_nodal_residual),
                             rcpFromRef(g_nodal_delta),
                             cg_iteration_limit,
                             cg_iteration_tolerance,
                             print_flag);
          perf.mat_vec_time    += cgsolve.matvec_time ;
          perf.cg_iter_time    += cgsolve.iter_time ;
          perf.prec_setup_time += cgsolve.prec_setup_time ;
          perf.prec_apply_time += cgsolve.prec_apply_time ;
          perf.cg_total_time   += cgsolve.total_time ;
        }
        perf.cg_iter_count   += cgsolve.iteration ;

        // Update solution vector

        g_nodal_solution_no_overlap.update(-1.0,g_nodal_delta,1.0);

        //--------------------------------

        if ( print_flag ) {
          const double delta_norm = g_nodal_delta.norm2();

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
      }

      // Evaluate response function -- currently 2-norm of solution vector

      response = responseFunc.apply();
      response = Kokkos::Example::all_reduce( response , comm );
    }
};


template < class Scalar, class Device , BoxElemPart::ElemOrder ElemOrder,
           class CoeffFunctionType >
Perf fenl(
  const Teuchos::RCP<const Teuchos::Comm<int> >& comm ,
  const Teuchos::RCP<  typename ::Kokkos::Compat::KokkosDeviceWrapperNode<Device> >& node,
  const std::string& fenl_xml_file,
  const int use_print ,
  const int use_trials ,
  const int use_atomic ,
  const int use_belos ,
  const int use_muelu ,
  const int use_mean_based ,
  const int use_nodes[] ,
  const CoeffFunctionType& coeff_function ,
  const double coeff_source ,
  const double coeff_advection ,
  const double bc_lower_value ,
  const double bc_upper_value ,
  Scalar& response,
  const QuadratureData<Device>& qd = QuadratureData<Device>() )
{
  typedef typename Kokkos::Details::ArithTraits<Scalar>::mag_type  Magnitude;

  typedef Problem< Scalar, Device , ElemOrder, CoeffFunctionType > ProblemType ;

  const int print_flag = use_print && Kokkos::Impl::VerifyExecutionCanAccessMemorySpace< Kokkos::HostSpace::execution_space , typename Device::memory_space >::value ;

  // Read in any params from xml file
  Teuchos::RCP<Teuchos::ParameterList> fenlParams = Teuchos::parameterList();
  Teuchos::updateParametersFromXmlFileAndBroadcast(
    fenl_xml_file, fenlParams.ptr(), *comm);

  const double geom_bubble[3] = { 1.0 , 1.0 , 1.0 };

  const unsigned  newton_iteration_limit =
    fenlParams->get("Max Nonlinear Iterations", 10) ;
  const Magnitude newton_iteration_tolerance =
    fenlParams->get("Nonlinear Solver Tolerance", 1e-7) ;
  const unsigned  cg_iteration_limit =
    fenlParams->get("Max Linear Iterations", 2000) ;
  const Magnitude cg_iteration_tolerance =
    fenlParams->get("Linear Solver Tolerance", 1e-7) ;

  //------------------------------------
  // Problem setup:

  ProblemType problem( comm , node , use_nodes , geom_bubble , print_flag );

  //------------------------------------

  Kokkos::Impl::Timer wall_clock ;

  Perf perf_stats = Perf() ;

  for ( int itrial = 0 ; itrial < use_trials ; ++itrial ) {

    problem.solve( coeff_function
                 , coeff_source
                 , coeff_advection
                 , bc_lower_value
                 , bc_upper_value
                 , newton_iteration_limit
                 , newton_iteration_tolerance
                 , cg_iteration_limit
                 , cg_iteration_tolerance
                 , qd
                 , use_atomic
                 , use_belos
                 , use_muelu
                 , use_mean_based
                 , print_flag
                 , fenlParams
                 );

    if ( 0 == itrial ) {
      response   = problem.response ;
      perf_stats = problem.perf ;
    }
    else {
      perf_stats.fill_node_set =
        std::min( perf_stats.fill_node_set , problem.perf.fill_node_set );
      perf_stats.scan_node_count =
        std::min( perf_stats.scan_node_count , problem.perf.scan_node_count );
      perf_stats.fill_graph_entries =
        std::min( perf_stats.fill_graph_entries , problem.perf.fill_graph_entries );
      perf_stats.sort_graph_entries =
        std::min( perf_stats.sort_graph_entries , problem.perf.sort_graph_entries );
      perf_stats.fill_element_graph =
        std::min( perf_stats.fill_element_graph , problem.perf.fill_element_graph );
      perf_stats.create_sparse_matrix =
        std::min( perf_stats.create_sparse_matrix , problem.perf.create_sparse_matrix );
       perf_stats.import_time =
        std::min( perf_stats.import_time , problem.perf.import_time );
      perf_stats.fill_time =
        std::min( perf_stats.fill_time , problem.perf.fill_time );
      perf_stats.bc_time =
        std::min( perf_stats.bc_time , problem.perf.bc_time );
      perf_stats.mat_vec_time =
        std::min( perf_stats.mat_vec_time , problem.perf.mat_vec_time );
      perf_stats.cg_iter_time =
        std::min( perf_stats.cg_iter_time , problem.perf.cg_iter_time );
      perf_stats.prec_setup_time =
        std::min( perf_stats.prec_setup_time , problem.perf.prec_setup_time );
      perf_stats.prec_apply_time =
        std::min( perf_stats.prec_apply_time , problem.perf.prec_apply_time );
      perf_stats.cg_total_time =
        std::min( perf_stats.cg_total_time , problem.perf.cg_total_time );
    }
  }

  return perf_stats ;
}

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos */

#endif /* #ifndef KOKKOS_EXAMPLE_FENL_IMPL_HPP */
