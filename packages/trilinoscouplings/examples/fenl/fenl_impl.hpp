//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_EXAMPLE_FENL_IMPL_HPP
#define KOKKOS_EXAMPLE_FENL_IMPL_HPP

#include <math.h>

// Kokkos libraries' headers:

#include <Kokkos_UnorderedMap.hpp>
#include <Kokkos_StaticCrsGraph.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <Kokkos_Timer.hpp>
#include <Kokkos_ArithTraits.hpp>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

// Examples headers:

#include <BoxElemFixture.hpp>
#include <fenl.hpp>
#include <fenl_functors.hpp>
#include <Tpetra_KokkosCompat_DefaultNode.hpp>

#include <Tpetra_Vector.hpp>
#include "Tpetra_MultiVector.hpp"

#include <CGSolve.hpp>
#include <BelosSolve.hpp>

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
  Kokkos::parallel_for(lid_to_gid.extent(0), F(lid_to_gid, fixture));
}

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos */

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {
namespace FENL {

template < class Scalar, class Device , BoxElemPart::ElemOrder ElemOrder >
class Problem {
public:

  typedef Scalar ScalarType;
  typedef Device DeviceType;
  typedef BoxElemFixture< Device , ElemOrder >  FixtureType ;

  typedef typename Kokkos::ArithTraits<Scalar>::mag_type  Magnitude;

  typedef Tpetra::KokkosCompat::KokkosDeviceWrapperNode< Device >  NodeType;

  typedef Kokkos::View<     Scalar * , Kokkos::LayoutLeft, Device >  LocalVectorType ;
  typedef Kokkos::View<     Scalar** , Kokkos::LayoutLeft, Device >  LocalMultiVectorType ;
  typedef Kokkos::DualView< Scalar** , Kokkos::LayoutLeft, Device >  LocalDualVectorType;

  typedef Tpetra::Map<>::local_ordinal_type LocalOrdinalType;
  typedef Tpetra::Map<>::global_ordinal_type GlobalOrdinalType;
  typedef Tpetra::Map<LocalOrdinalType,GlobalOrdinalType,NodeType> MapType;
  typedef Tpetra::Vector<Scalar,LocalOrdinalType,GlobalOrdinalType,NodeType> GlobalVectorType;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinalType,GlobalOrdinalType,NodeType> GlobalMultiVectorType;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinalType,GlobalOrdinalType,NodeType> GlobalMatrixType;
  typedef typename GlobalMatrixType::local_matrix_device_type LocalMatrixType;
  typedef typename LocalMatrixType::StaticCrsGraphType LocalGraphType;


  typedef NodeNodeGraph< typename FixtureType::elem_node_type
                       , LocalGraphType
                       , FixtureType::ElemNode
                       > NodeNodeGraphType ;

  typedef typename NodeNodeGraphType::ElemGraphType  ElemGraphType ;

  typedef Teuchos::RCP<const MapType>                        rcpMapType;
  typedef Teuchos::RCP<const Teuchos::Comm<int> >            rcpCommType ;
  typedef Teuchos::RCP<Tpetra::KokkosCompat::KokkosDeviceWrapperNode<Device> > rcpNodeType ;

  typedef Tpetra::Import<
    typename GlobalVectorType::local_ordinal_type,
    typename GlobalVectorType::global_ordinal_type,
    typename GlobalVectorType::node_type
    > import_type;

private:

  typedef Kokkos::View<GlobalOrdinalType*,Device> lid_to_gid_type;

  rcpMapType create_row_map()
  {
    lid_to_gid_type lid_to_gid_row( "lig_to_gid_row", fixture.node_count_owned() );

    build_lid_to_gid(lid_to_gid_row, fixture);

    typename lid_to_gid_type::HostMirror lid_to_gid_row_host =
      Kokkos::create_mirror_view(lid_to_gid_row);

    Kokkos::deep_copy(lid_to_gid_row_host, lid_to_gid_row);

    return  Teuchos::rcp (new MapType( fixture.node_count_global(),
        Teuchos::arrayView(lid_to_gid_row_host.data(),
                           lid_to_gid_row_host.extent(0)),
        0, comm));
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
        Teuchos::arrayView( lid_to_gid_all_host.data(),
                            lid_to_gid_all_host.extent(0)),
        0,comm) );
  }

public:

  const rcpCommType  comm ;
  const FixtureType  fixture ;

private:

  typename NodeNodeGraphType::Times graph_times;

  const NodeNodeGraphType mesh_to_graph ;

public:

  const ElemGraphType elem_graph ;

  const rcpMapType   RowMap ;
  const rcpMapType   ColMap ;
  const import_type  import ;

  GlobalVectorType      g_nodal_solution ;
  GlobalVectorType      g_nodal_residual ;
  GlobalVectorType      g_nodal_delta ;
  GlobalVectorType      g_nodal_solution_no_overlap ;
  GlobalMultiVectorType g_nodal_solution_dp ;
  GlobalMultiVectorType g_nodal_residual_dp ;
  GlobalMultiVectorType g_nodal_delta_dp ;
  GlobalMultiVectorType g_nodal_solution_no_overlap_dp ;
  GlobalMatrixType      g_jacobian ;

  Perf                   perf ;
  bool                   print_flag ;
  unsigned               num_sensitivities ;

  Problem( const rcpCommType & use_comm
         , const int use_nodes[]
         , const double grid_bubble[]
         , const bool use_print
         , const unsigned num_sens
         )
    : comm( use_comm )
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
    , g_nodal_solution_no_overlap (g_nodal_solution, RowMap)
    , g_jacobian( RowMap, ColMap, LocalMatrixType( "jacobian" , mesh_to_graph.graph , maximum_entry(mesh_to_graph.graph) + 1 ) )
    , perf()
    , num_sensitivities(num_sens)
    {
      if ( maximum(*comm, ( fixture.ok() ? 0 : 1 ) ) ) {
        throw std::runtime_error(std::string("Problem fixture setup failed"));
      }

      print_flag = use_print && Kokkos::SpaceAccessibility< Kokkos::HostSpace::execution_space , typename Device::memory_space >::accessible ;

      perf.global_elem_count  = fixture.elem_count_global();
      perf.global_node_count  = fixture.node_count_global();

      perf.map_ratio          = graph_times.ratio;
      perf.fill_node_set      = graph_times.fill_node_set;
      perf.scan_node_count    = graph_times.scan_node_count;
      perf.fill_graph_entries = graph_times.fill_graph_entries;
      perf.sort_graph_entries = graph_times.sort_graph_entries;
      perf.fill_element_graph = graph_times.fill_element_graph;

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

        LocalMatrixType jacobian = g_jacobian.getLocalMatrixDevice();

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
        for ( unsigned ielem = 0 ; ielem < elem_graph.extent(0) ; ++ielem ) {
          std::cout << "  elem[" << ielem << "]{" ;
          for ( unsigned irow = 0 ; irow < elem_graph.extent(1) ; ++irow ) {
            std::cout << " {" ;
            for ( unsigned icol = 0 ; icol < elem_graph.extent(2) ; ++icol ) {
              std::cout << " " << elem_graph(ielem,irow,icol);
            }
            std::cout << " }" ;
          }
          std::cout << " }" << std::endl ;
        }
        std::cout << "}" << std::endl ;
      }

      if (num_sensitivities > 0) {
        g_nodal_solution_dp =
          GlobalMultiVectorType( ColMap, num_sensitivities );
        g_nodal_residual_dp =
          GlobalMultiVectorType( RowMap, num_sensitivities );
        g_nodal_delta_dp =
          GlobalMultiVectorType( RowMap, num_sensitivities );
        g_nodal_solution_no_overlap_dp =
          GlobalMultiVectorType(g_nodal_solution_dp, *RowMap);
      }
    }

  //----------------------------------------

  template < class CoeffFunctionType >
  void solve( const CoeffFunctionType & coeff_function
            , const bool isotropic
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
            , const Teuchos::RCP<Teuchos::ParameterList>& fenlParams
            , Scalar& response
            , Teuchos::Array<Scalar>& response_gradient
            )
    {
      typedef ElementComputation< FixtureType , LocalMatrixType , CoeffFunctionType > ElementComputationType ;
      typedef DirichletComputation< FixtureType , LocalMatrixType > DirichletComputationType ;
      typedef ResponseComputation< FixtureType , LocalVectorType > ResponseComputationType ;

      Kokkos::Timer wall_clock ;
      Kokkos::Timer newton_clock ;
      newton_clock.reset();

      LocalMatrixType jacobian = g_jacobian.getLocalMatrixDevice();

      const auto k_nodal_solution = g_nodal_solution.getLocalViewDevice(Tpetra::Access::ReadWrite);
      const auto k_nodal_residual = g_nodal_residual.getLocalViewDevice(Tpetra::Access::ReadWrite);
      const auto k_nodal_delta    = g_nodal_delta   .getLocalViewDevice(Tpetra::Access::ReadWrite);

      const LocalVectorType nodal_solution =
        Kokkos::subview(k_nodal_solution,Kokkos::ALL(),0);
      const LocalVectorType nodal_residual =
        Kokkos::subview(k_nodal_residual,Kokkos::ALL(),0);
      const LocalVectorType nodal_delta =
        Kokkos::subview(k_nodal_delta,Kokkos::ALL(),0);

      LocalVectorType nodal_solution_no_overlap =
        Kokkos::subview(nodal_solution,std::pair<unsigned,unsigned>(0,fixture.node_count_owned()));

      // Get DeviceConfig structs used by some functors
      KokkosSparse::DeviceConfig dev_config_elem, dev_config_gath, dev_config_bc;

      CreateDeviceConfigs<Scalar>::eval( dev_config_elem,
                                         dev_config_gath,
                                         dev_config_bc );

      // Create element computation functor
      const ElementComputationType elemcomp( fixture , coeff_function , isotropic ,
                                             coeff_source , coeff_advection ,
                                             nodal_solution ,
                                             elem_graph ,
                                             jacobian , nodal_residual ,
                                             dev_config_elem , qd );

      // Create boundary condition functor
      // This also sets the boundary conditions in the solution vector
      // and zeros out non BC values
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

      // Teuchos::RCP<Teuchos::FancyOStream> out =
      //   Teuchos::fancyOStream(Teuchos::rcp(&std::cout,false));
      // out->setShowProcRank(true);

      Teuchos::RCP< Tpetra::Operator<Scalar,LocalOrdinalType,GlobalOrdinalType,NodeType> > precOp;
      for ( perf.newton_iter_count = 0 ;
            perf.newton_iter_count < newton_iteration_limit ;
            ++perf.newton_iter_count ) {

        //--------------------------------

        Device().fence();
        wall_clock.reset();
        g_nodal_solution.doImport (g_nodal_solution_no_overlap, import, Tpetra::REPLACE);

        // Take minimum import time across newton steps -- resolves strange
        // timings on titan where time after first solve is much larger
        Device().fence();
        if (perf.newton_iter_count == 0)
          perf.import_time = wall_clock.seconds();
        else
          perf.import_time = std::min( perf.import_time, wall_clock.seconds() );

        // if (itrial == 0 && perf.newton_iter_count == 0)
        //   g_nodal_solution_no_overlap.describe(*out, Teuchos::VERB_EXTREME);

        //--------------------------------
        // Element contributions to residual and jacobian

        Device().fence();
        wall_clock.reset();

        Kokkos::deep_copy( nodal_residual , 0.0 );
        Kokkos::deep_copy( jacobian.values ,0.0 );

        elemcomp.apply();

        Device().fence();
        if (perf.newton_iter_count == 0)
          perf.fill_time = wall_clock.seconds();
        else
          perf.fill_time = std::min( perf.fill_time, wall_clock.seconds() );

        //--------------------------------
        // Apply boundary conditions

        Device().fence();
        wall_clock.reset();

        dirichlet.apply();

        Device().fence();
        if (perf.newton_iter_count == 0)
          perf.bc_time = wall_clock.seconds();
        else
          perf.bc_time = std::min( perf.bc_time, wall_clock.seconds() );

        //--------------------------------
        // Evaluate convergence

        const Magnitude residual_norm = g_nodal_residual.norm2();

        perf.newton_residual = scalar_norm(residual_norm) ;
        perf.error_max = scalar_norm(residual_norm) ;

        if ( 0 == perf.newton_iter_count ) {
          residual_norm_init = residual_norm ;
        }

        if ( residual_norm < residual_norm_init * newton_iteration_tolerance ) {
          break ;
        }

        //--------------------------------
        // Solve for nonlinear update

        // Zero out newton update vector before solve
        g_nodal_delta.putScalar(0.0);

        result_struct cgsolve;
        if (use_belos) {
          // Destroy previous preconditioner (we could make reusing it optional)
          precOp = Teuchos::null;
          cgsolve = belos_solve(g_jacobian,
                                g_nodal_residual,
                                g_nodal_delta,
                                precOp,
                                fixture,
                                use_muelu,
                                use_mean_based ,
                                fenlParams ,
                                cg_iteration_limit ,
                                cg_iteration_tolerance);

          // Don't accumulate Belos times as the internal Teuchos timers
          // already accumulate
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
        const int ne = cgsolve.ensemble_its.size();
        perf.ensemble_cg_iter_count.resize(ne);
        for (int i=0; i<ne; ++i)
          perf.ensemble_cg_iter_count[i] += cgsolve.ensemble_its[i];

        // Update solution vector

        g_nodal_solution_no_overlap.update(-1.0,g_nodal_delta,1.0);

        //--------------------------------

        if ( print_flag ) {
          const Magnitude delta_norm = g_nodal_delta.norm2();

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
        //break;
      }

      // Evaluate response function -- currently 2-norm of solution vector

      response = responseFunc.apply();
      response = Kokkos::Example::all_reduce( response , comm );

      // Compute sensitivity of response function if requested
      if (response_gradient.size() > 0 && num_sensitivities > 0) {
#ifdef HAVE_TRILINOSCOUPLINGS_SACADO
        typedef Sacado::Fad::SLFad<Scalar,10> FadType;
        typedef FadCoeffFunctionTraits<CoeffFunctionType, FadType> CoeffTraits;
        typedef typename CoeffTraits::type FadCoeffFunctionType;
        typedef ParamSensitivityGatherScatterOp<LocalMultiVectorType> SensGatherScatter;
        typedef ElementComputation< FixtureType , LocalMatrixType , FadCoeffFunctionType, SensGatherScatter, FadType > FadElementComputationType ;

        // Check sizes match
        TEUCHOS_TEST_FOR_EXCEPTION(
          response_gradient.size() != num_sensitivities,
          std::logic_error,
          "Response gradient length must match number of sensitivities specified in constuctor");

        const auto nodal_solution_dp =
          g_nodal_solution_dp.getLocalViewDevice(Tpetra::Access::ReadWrite);
        const auto nodal_residual_dp =
          g_nodal_residual_dp.getLocalViewDevice(Tpetra::Access::ReadWrite);
        const auto nodal_delta_dp =
          g_nodal_delta_dp.getLocalViewDevice(Tpetra::Access::ReadWrite);

        LocalMultiVectorType nodal_solution_no_overlap_dp =
          Kokkos::subview(
            nodal_solution_dp,
            std::pair<unsigned,unsigned>(0,fixture.node_count_owned()),
            Kokkos::ALL());

        const FadCoeffFunctionType coeff_function_dp =
          CoeffTraits::eval(coeff_function);

        const SensGatherScatter gather_scatter(nodal_solution_dp,
                                               nodal_residual_dp);
        const FadElementComputationType elemcomp_dp(
          fixture , coeff_function_dp , isotropic ,
          coeff_source , coeff_advection ,
          nodal_solution ,
          elem_graph ,
          jacobian , nodal_residual ,
          dev_config_elem , qd , gather_scatter , false );

        const DirichletComputationType dirichlet_dp(
          fixture , nodal_solution , jacobian , nodal_residual ,
          2 /* apply at 'z' ends */ ,
          bc_lower_value ,
          bc_upper_value ,
          dev_config_bc ,
          nodal_solution_dp , nodal_residual_dp , false , true );

        const ResponseComputationType responseFunc_dp(
          fixture , nodal_solution_no_overlap , nodal_solution_no_overlap_dp );

        g_nodal_solution.doImport (g_nodal_solution_no_overlap, import, Tpetra::REPLACE);
        g_nodal_solution_dp.doImport (g_nodal_solution_no_overlap_dp, import, Tpetra::REPLACE);

        Device().fence();
        wall_clock.reset();

        Kokkos::deep_copy( nodal_residual , 0.0 );
        Kokkos::deep_copy( nodal_residual_dp , 0.0 );

        elemcomp_dp.apply();
        dirichlet_dp.apply();

        Device().fence();
        perf.tangent_fill_time = wall_clock.seconds();

        result_struct cgsolve;
        if (use_belos) {
          // Reuse preconditioner from forward solve
          cgsolve = belos_solve(g_jacobian,
                                g_nodal_residual_dp,
                                g_nodal_delta_dp,
                                precOp,
                                fixture,
                                use_muelu,
                                use_mean_based ,
                                fenlParams ,
                                cg_iteration_limit ,
                                cg_iteration_tolerance);

          // Don't accumulate Belos times as the internal Teuchos timers
          // already accumulate
          perf.mat_vec_time    = cgsolve.matvec_time ;
          perf.cg_iter_time    = cgsolve.iter_time ;
          perf.prec_setup_time = cgsolve.prec_setup_time ;
          perf.prec_apply_time = cgsolve.prec_apply_time ;
          perf.cg_total_time   = cgsolve.total_time ;
          perf.cg_iter_count  += cgsolve.iteration ;
        }
        else {
          for (unsigned i=0; i<num_sensitivities; ++i) {
            cgsolve = cg_solve(rcpFromRef(g_jacobian),
                               g_nodal_residual_dp.getVectorNonConst(i),
                               g_nodal_delta_dp.getVectorNonConst(i),
                               cg_iteration_limit,
                               cg_iteration_tolerance,
                               print_flag);
            perf.mat_vec_time    += cgsolve.matvec_time ;
            perf.cg_iter_time    += cgsolve.iter_time ;
            perf.prec_setup_time += cgsolve.prec_setup_time ;
            perf.prec_apply_time += cgsolve.prec_apply_time ;
            perf.cg_total_time   += cgsolve.total_time ;
            perf.cg_iter_count   += cgsolve.iteration ;
          }
        }
        const int ne = cgsolve.ensemble_its.size();
        perf.ensemble_cg_iter_count.resize(ne);
        for (int i=0; i<ne; ++i)
          perf.ensemble_cg_iter_count[i] += cgsolve.ensemble_its[i];

        g_nodal_solution_no_overlap_dp.update(-1.0,g_nodal_delta_dp,1.0);

        Teuchos::Array<Scalar> response_and_gradient =
          responseFunc_dp.apply_gradient();

        response_and_gradient =
          Kokkos::Example::all_reduce( response_and_gradient , comm );
        for (unsigned i=0; i<num_sensitivities; ++i)
          response_gradient[i] = response_and_gradient[i+1];
        response = response_and_gradient[0];
#endif
      }

      Device().fence();
      perf.newton_total_time = newton_clock.seconds();
    }
};

template < class Scalar, class Device , BoxElemPart::ElemOrder ElemOrder,
           class CoeffFunctionType >
Perf fenl(
  Problem< Scalar, Device , ElemOrder >& problem,
  const Teuchos::RCP<Teuchos::ParameterList>& fenlParams,
  const int use_print ,
  const int use_trials ,
  const int use_atomic ,
  const int use_belos ,
  const int use_muelu ,
  const int use_mean_based ,
  const CoeffFunctionType& coeff_function ,
  const bool isotropic,
  const double coeff_source ,
  const double coeff_advection ,
  const double bc_lower_value ,
  const double bc_upper_value ,
  Scalar& response,
  Teuchos::Array<Scalar>& response_gradient,
  const QuadratureData<Device>& qd = QuadratureData<Device>() )
{
  typedef typename Kokkos::ArithTraits<Scalar>::mag_type  Magnitude;

  const unsigned  newton_iteration_limit =
    fenlParams->get("Max Nonlinear Iterations", 10) ;
  const Magnitude newton_iteration_tolerance =
    fenlParams->get("Nonlinear Solver Tolerance", 1e-7) ;
  const unsigned  cg_iteration_limit =
    fenlParams->get("Max Linear Iterations", 2000) ;
  const Magnitude cg_iteration_tolerance =
    fenlParams->get("Linear Solver Tolerance", 1e-7) ;

  //------------------------------------

  Kokkos::Timer wall_clock ;

  Perf perf_stats = Perf() ;

  // Since the perf struc inside Problem is reused each time solve() is called
  // zero out some stats that we don't want to accumulate
  problem.perf.mat_vec_time = 0;
  problem.perf.cg_iter_time = 0;
  problem.perf.prec_setup_time = 0;
  problem.perf.prec_apply_time  = 0;
  problem.perf.cg_total_time = 0;
  problem.perf.newton_total_time = 0;
  problem.perf.cg_iter_count = 0;
  problem.perf.ensemble_cg_iter_count.clear();
  problem.perf.import_time = 0;
  problem.perf.fill_time = 0;
  problem.perf.tangent_fill_time = 0;
  problem.perf.bc_time = 0;

  for ( int itrial = 0 ; itrial < use_trials ; ++itrial ) {

    problem.solve( coeff_function
                 , isotropic
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
                 , fenlParams
                 , response
                 , response_gradient
                 );

    problem.perf.reduceMax(*problem.comm);

    if ( 0 == itrial ) {
      perf_stats = problem.perf ;
    }
    else {
      perf_stats.min(problem.perf);
    }
  }

  return perf_stats ;
}

template < class Scalar, class Device , BoxElemPart::ElemOrder ElemOrder,
           class CoeffFunctionType >
Perf fenl(
  const Teuchos::RCP<const Teuchos::Comm<int> >& comm ,
  const std::string& fenl_xml_file,
  const int use_print ,
  const int use_trials ,
  const int use_atomic ,
  const int use_belos ,
  const int use_muelu ,
  const int use_mean_based ,
  const int use_nodes[] ,
  const CoeffFunctionType& coeff_function ,
  const bool isotropic,
  const double coeff_source ,
  const double coeff_advection ,
  const double bc_lower_value ,
  const double bc_upper_value ,
  Scalar& response,
  Teuchos::Array<Scalar>& response_gradient,
  const QuadratureData<Device>& qd = QuadratureData<Device>() )
{

  typedef Problem< Scalar, Device , ElemOrder > ProblemType ;

  //------------------------------------
  // Read in any params from xml file
  Teuchos::RCP<Teuchos::ParameterList> fenlParams = Teuchos::parameterList();
  Teuchos::updateParametersFromXmlFileAndBroadcast(
    fenl_xml_file, fenlParams.ptr(), *comm);

  //------------------------------------
  // Problem setup:

  const double geom_bubble[3] = { 1.0 , 1.0 , 1.0 };
  ProblemType problem( comm , use_nodes , geom_bubble , use_print ,
                       response_gradient.size() );

  //------------------------------------
  // Solve

  return fenl( problem,
               fenlParams,
               use_print,
               use_trials,
               use_atomic,
               use_belos,
               use_muelu,
               use_mean_based,
               coeff_function,
               isotropic,
               coeff_source,
               coeff_advection,
               bc_lower_value,
               bc_upper_value,
               response,
               response_gradient,
               qd );
}

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos */

#endif /* #ifndef KOKKOS_EXAMPLE_FENL_IMPL_HPP */
