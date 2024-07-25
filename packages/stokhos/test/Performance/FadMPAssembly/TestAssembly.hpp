// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

// MP::Vector
#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include "Kokkos_View_MP_Vector.hpp"

// Compile-time loops
#include "Sacado_mpl_range_c.hpp"
#include "Sacado_mpl_for_each.hpp"
#include "Sacado_mpl_integral_c.hpp"

// Kokkos libraries' headers:
#include <Kokkos_UnorderedMap.hpp>
#include <Kokkos_StaticCrsGraph.hpp>
#include <Kokkos_Timer.hpp>

// Utilities
#include <Teuchos_CommHelpers.hpp>
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_VerboseObject.hpp"

// FENL
#include <BoxElemFixture.hpp>
#include <VectorImport.hpp>
#include <fenl_functors.hpp>

struct Perf {
  size_t global_elem_count ;
  size_t global_node_count ;
  double import_time ;
  double fill_time ;

  Perf() : global_elem_count(0) ,
           global_node_count(0) ,
           import_time(0) ,
           fill_time(0) {}

  void increment(const Perf& p) {
    global_elem_count = p.global_elem_count;
    global_node_count = p.global_node_count;
    import_time      += p.import_time;
    fill_time        += p.fill_time;
  }

  void scale(double s) {
    import_time *= s;
    fill_time   *= s;
  }
};

inline
double maximum( const Teuchos::RCP<const Teuchos::Comm<int> >& comm , double local )
{
  double global = 0 ;
  Teuchos::reduceAll( *comm , Teuchos::REDUCE_MAX , 1 , & local , & global );
  return global ;
}

template <typename Scalar, typename Device,
          Kokkos::Example::FENL::AssemblyMethod Method>
Perf fenl_assembly(
  const Teuchos::RCP<const Teuchos::Comm<int> >& comm ,
  const int use_print ,
  const int use_trials ,
  const int use_nodes[] ,
  Kokkos::Example::FENL::DeviceConfig dev_config ,
  Kokkos::View< Scalar* , Kokkos::LayoutLeft, Device >& nodal_residual)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::arrayView;
  using Teuchos::ParameterList;

  typedef Kokkos::Example::BoxElemFixture< Device , Kokkos::Example::BoxElemPart::ElemLinear > FixtureType ;

  typedef Kokkos::Example::FENL::CrsMatrix< Scalar , Device > LocalMatrixType ;

  typedef typename LocalMatrixType::StaticCrsGraphType LocalGraphType ;

  typedef Kokkos::Example::FENL::NodeNodeGraph< typename FixtureType::elem_node_type , LocalGraphType , FixtureType::ElemNode >
     NodeNodeGraphType ;

  //typedef Kokkos::Example::FENL::ElementComputationConstantCoefficient CoeffFunctionType;
  typedef Kokkos::Example::FENL::ExponentialKLCoefficient< Scalar, double, Device > CoeffFunctionType;
  typedef Kokkos::Example::FENL::ElementComputation< FixtureType , LocalMatrixType , Method, CoeffFunctionType >
    ElementComputationType ;

  // typedef Kokkos::Example::FENL::DirichletComputation< FixtureType , LocalMatrixType >
  //   DirichletComputationType ;

  typedef typename ElementComputationType::vector_type VectorType ;

   typedef Kokkos::Example::VectorImport<
     typename FixtureType::comm_list_type ,
     typename FixtureType::send_nodeid_type ,
     VectorType > ImportType ;

  //------------------------------------

  const int print_flag = use_print && std::is_same< Kokkos::HostSpace , typename Device::memory_space >::value ;

  const int comm_rank = comm->getRank();
  const int comm_size = comm->getSize();

  // Decompose by node to avoid parallel communication in assembly

  const double bubble_x = 1.0 ;
  const double bubble_y = 1.0 ;
  const double bubble_z = 1.0 ;

  const FixtureType fixture( Kokkos::Example::BoxElemPart::DecomposeNode ,
                             comm_size , comm_rank ,
                             use_nodes[0] , use_nodes[1] , use_nodes[2] ,
                             bubble_x , bubble_y , bubble_z );

  //------------------------------------

  const ImportType comm_nodal_import(
    comm ,
    fixture.recv_node() ,
    fixture.send_node() ,
    fixture.send_nodeid() ,
    fixture.node_count_owned() ,
    fixture.node_count() - fixture.node_count_owned() );

  //------------------------------------

  // const double bc_lower_value = 1 ;
  // const double bc_upper_value = 2 ;
  //CoeffFunctionType diffusion_coefficient( 1.0 );
  CoeffFunctionType diffusion_coefficient( 1.0, 0.1, 1.0, 5 );
  Kokkos::deep_copy( diffusion_coefficient.getRandomVariables(), 1.0 );

  //------------------------------------

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

  Kokkos::Timer wall_clock ;

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
      mesh_to_graph( fixture.elem_node() , fixture.node_count_owned(),
                     graph_times );

    // Create the local sparse matrix from the graph:
    LocalMatrixType jacobian( mesh_to_graph.graph );

    //----------------------------------

    // Allocate solution vector for each node in the mesh and residual vector for each owned node
    VectorType nodal_solution( "nodal_solution" , fixture.node_count() );
    nodal_residual = VectorType( "nodal_residual" , fixture.node_count_owned() );

    // Get DeviceConfig structs used by some functors
    Kokkos::Example::FENL::DeviceConfig dev_config_elem, dev_config_bc;
    Kokkos::Example::FENL::CreateDeviceConfigs<Scalar>::eval( dev_config_elem,
                                                              dev_config_bc );

    // Create element computation functor
    const ElementComputationType elemcomp( fixture , diffusion_coefficient ,
                                           nodal_solution ,
                                           mesh_to_graph.elem_graph ,
                                           jacobian , nodal_residual ,
                                           dev_config_elem );

    // Create boundary condition functor
    // const DirichletComputationType dirichlet(
    //   fixture , nodal_solution , jacobian , nodal_residual ,
    //   2 /* apply at 'z' ends */ ,
    //   bc_lower_value ,
    //   bc_upper_value ,
    //   dev_config_bc );

    Kokkos::deep_copy( nodal_solution , Scalar(1) );

    //--------------------------------

    wall_clock.reset();

    comm_nodal_import( nodal_solution );

    Device().fence();
    perf.import_time = maximum( comm , wall_clock.seconds() );

    //--------------------------------
    // Element contributions to residual and jacobian

    wall_clock.reset();

    Kokkos::deep_copy( nodal_residual , Scalar(0) );
    Kokkos::deep_copy( jacobian.values , Scalar(0) );

    elemcomp.apply();

    //--------------------------------
    // Apply boundary conditions

    //dirichlet.apply();

    Device().fence();
    perf.fill_time = maximum( comm , wall_clock.seconds() );

    //--------------------------------

    perf_stats.increment(perf);

  }

  return perf_stats ;
}

template <typename ScalarViewType, typename EnsembleViewType>
bool check_residuals(const ScalarViewType& scalar_residual,
                     const EnsembleViewType& ensemble_residual)
{
  const double tol = 1e-14;
  bool success = true;
  Teuchos::RCP<Teuchos::FancyOStream> out =
    Teuchos::VerboseObjectBase::getDefaultOStream();
  std::stringstream buf;
  Teuchos::FancyOStream fbuf(Teuchos::rcp(&buf,false));

  typename ScalarViewType::HostMirror host_scalar_residual =
    Kokkos::create_mirror_view(scalar_residual);
  typename EnsembleViewType::HostMirror host_ensemble_residual =
    Kokkos::create_mirror_view(ensemble_residual);
  Kokkos::deep_copy( host_scalar_residual, scalar_residual );
  Kokkos::deep_copy( host_ensemble_residual, ensemble_residual );

  TEUCHOS_TEST_EQUALITY( host_scalar_residual.extent(0),
                         host_ensemble_residual.extent(0), fbuf, success );

  const size_t num_node = host_scalar_residual.extent(0);
  const size_t num_ensemble = Kokkos::dimension_scalar(host_ensemble_residual);
  for (size_t i=0; i<num_node; ++i) {
    for (size_t j=0; j<num_ensemble; ++j) {
      TEUCHOS_TEST_FLOATING_EQUALITY(
        host_scalar_residual(i), host_ensemble_residual(i).fastAccessCoeff(j),
        tol, fbuf, success );
    }
  }

  if (!success)
    *out << buf.str();

  return success;
}

template <class Storage,
          Kokkos::Example::FENL::AssemblyMethod Method>
struct PerformanceDriverOp {
  typedef typename Storage::value_type Scalar;
  typedef typename Storage::ordinal_type Ordinal;
  typedef typename Storage::execution_space Device;
  Teuchos::RCP<const Teuchos::Comm<int> > comm ;
  const int use_print ;
  const int use_trials ;
  const int *use_nodes ;
  const bool check     ;
  Kokkos::Example::FENL::DeviceConfig dev_config;

  PerformanceDriverOp(const Teuchos::RCP<const Teuchos::Comm<int> >& comm_ ,
                      const int use_print_ ,
                      const int use_trials_ ,
                      const int use_nodes_[] ,
                      const bool check_ ,
                      Kokkos::Example::FENL::DeviceConfig dev_config_) :
    comm(comm_),
    use_print(use_print_),
    use_trials(use_trials_),
    use_nodes(use_nodes_),
    check(check_),
    dev_config(dev_config_) {}

  template <typename ArgT>
  void operator() (ArgT arg) const {
    const int ensemble = ArgT::value;
    typedef typename Storage::template apply_N<ensemble> NewStorageApply;
    typedef typename NewStorageApply::type storage_type;
    typedef Sacado::MP::Vector<storage_type> mp_vector_type;

    typedef Kokkos::View< Scalar* , Kokkos::LayoutLeft, Device > scalar_vector_type ;
    typedef Kokkos::View< mp_vector_type* , Kokkos::LayoutLeft, Device > ensemble_vector_type ;

    scalar_vector_type scalar_residual;
    Kokkos::Example::FENL::DeviceConfig scalar_dev_config(0, 1, 1);
    Perf perf_scalar =
      fenl_assembly<Scalar,Device,Method>(
        comm, use_print, use_trials*ensemble, use_nodes,
        scalar_dev_config, scalar_residual );

    ensemble_vector_type ensemble_residual;
    Kokkos::Example::FENL::DeviceConfig ensemble_dev_config = dev_config;
#if defined( KOKKOS_ENABLE_CUDA )
    const bool is_cuda = std::is_same<Device,Kokkos::Cuda>::value;
#else
    const bool is_cuda = false ;
#endif
    if (is_cuda) {
      const size_t block_size = dev_config.block_dim.x * dev_config.block_dim.y;
      ensemble_dev_config.block_dim.x = ensemble;
      ensemble_dev_config.block_dim.y = block_size/ensemble;
    }
    Perf perf_ensemble =
      fenl_assembly<mp_vector_type,Device,Method>(
        comm, use_print, use_trials, use_nodes,
        ensemble_dev_config, ensemble_residual);

    if (check)
      check_residuals( scalar_residual, ensemble_residual );

    double s =
      1000.0 / ( use_trials * ensemble * perf_scalar.global_node_count );
    perf_scalar.scale(s);
    perf_ensemble.scale(s);

    if (comm->getRank() == 0) {
      std::cout.precision(3);
      std::cout << use_nodes[0] << " , "
                << perf_scalar.global_node_count << " , "
                << std::setw(2) << ensemble << " , "
                << std::scientific
                << perf_scalar.import_time << " , "
                << perf_ensemble.import_time << " , "
                << std::fixed << std::setw(6)
                << perf_scalar.import_time / perf_ensemble.import_time << " , "
                << std::scientific
                << perf_scalar.fill_time << " , "
                << perf_ensemble.fill_time << " , "
                << std::fixed << std::setw(6)
                << perf_scalar.fill_time / perf_ensemble.fill_time << " , "
                << std::endl;
    }
  }
};

template <class Storage, int entry_min, int entry_max, int entry_step,
          Kokkos::Example::FENL::AssemblyMethod Method>
void performance_test_driver( const Teuchos::RCP<const Teuchos::Comm<int> >& comm ,
                              const int use_print ,
                              const int use_trials ,
                              const int use_nodes[] ,
                              const bool check ,
                              Kokkos::Example::FENL::DeviceConfig dev_config)
{
  if (comm->getRank() == 0) {
    std::cout.precision(8);
    std::cout << std::endl
              << "\"Grid Size\" , "
              << "\"FEM Size\" , "
              << "\"Ensemble Size\" , "
              << "\"Scalar Import Time\" , "
              << "\"Ensemble Import Time\" , "
              << "\"Ensemble Import Speedup\" , "
              << "\"Scalar Fill Time\" , "
              << "\"Ensemble Fill Time\" , "
              << "\"Ensemble Fill Speedup\" , "
              << std::endl;
  }

  // Loop over [entry_min, entry_max] vector entries per thread
  typedef Sacado::mpl::range_c< int, entry_min, entry_max+1, entry_step > Range;
  PerformanceDriverOp<Storage,Method> op(comm, use_print, use_trials,
                                         use_nodes, check, dev_config);
  Sacado::mpl::for_each_no_kokkos<Range> f(op);
}
