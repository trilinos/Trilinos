// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

// Utilities
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_VerboseObject.hpp"

// FENL
#include <BoxElemFixture.hpp>
#include <fenl_functors.hpp>

struct Perf {
  size_t global_elem_count ;
  size_t global_node_count ;
  double fill_time ;

  Perf() : global_elem_count(0) ,
           global_node_count(0) ,
           fill_time(0) {}

  void increment(const Perf& p) {
    global_elem_count = p.global_elem_count;
    global_node_count = p.global_node_count;
    fill_time        += p.fill_time;
  }

  void scale(double s) {
    fill_time   *= s;
  }
};

template <typename Scalar, typename Device,
          Kokkos::Example::BoxElemPart::ElemOrder Order,
          Kokkos::Example::FENL::AssemblyMethod Method>
Perf fenl_assembly(
  const int use_print ,
  const int use_trials ,
  const int use_nodes[] ,
  Kokkos::View< Scalar* , Device >& residual,
  Kokkos::Example::FENL::CrsMatrix<Scalar,Device>& jacobian)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::arrayView;

  typedef Kokkos::Example::BoxElemFixture< Device , Order > FixtureType ;

  typedef Kokkos::Example::FENL::CrsMatrix< Scalar , Device > LocalMatrixType ;
  typedef typename LocalMatrixType::StaticCrsGraphType LocalGraphType ;

  typedef Kokkos::Example::FENL::NodeNodeGraph< typename FixtureType::elem_node_type , LocalGraphType , FixtureType::ElemNode > NodeNodeGraphType ;

  typedef Kokkos::Example::FENL::ElementComputation< FixtureType , LocalMatrixType , Method > ElementComputationType ;

  typedef typename ElementComputationType::vector_type VectorType ;

  //------------------------------------

  // Decompose by node to avoid parallel communication in assembly

  const double bubble_x = 1.0 ;
  const double bubble_y = 1.0 ;
  const double bubble_z = 1.0 ;

  const FixtureType fixture( Kokkos::Example::BoxElemPart::DecomposeNode ,
                             1 , 0 ,
                             use_nodes[0] , use_nodes[1] , use_nodes[2] ,
                             bubble_x , bubble_y , bubble_z );

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
    jacobian = LocalMatrixType( mesh_to_graph.graph );

    //----------------------------------

    // Allocate solution vector for each node in the mesh and residual vector for each owned node
    VectorType solution( "solution" , fixture.node_count() );
    residual = VectorType( "residual" , fixture.node_count_owned() );

    // Create element computation functor
    const ElementComputationType elemcomp( fixture , solution ,
                                           mesh_to_graph.elem_graph ,
                                           jacobian , residual );

    Kokkos::deep_copy( solution , Scalar(1.2345) );

    //--------------------------------
    // Element contributions to residual and jacobian

    Kokkos::deep_copy( residual , Scalar(0) );
    Kokkos::deep_copy( jacobian.coeff , Scalar(0) );

    wall_clock.reset();

    elemcomp.apply();

    Device().fence();
    perf.fill_time = wall_clock.seconds();

    //--------------------------------

    perf_stats.increment(perf);

  }

  return perf_stats ;
}

template<class ValueType>
bool compareValues(const ValueType& a1,
                   const std::string& a1_name,
                   const ValueType&a2,
                   const std::string& a2_name,
                   const ValueType& rel_tol, const ValueType& abs_tol,
                   Teuchos::FancyOStream& out)
{
  bool success = true;

  ValueType err = std::abs(a1 - a2);
  ValueType tol = abs_tol + rel_tol*std::max(std::abs(a1),std::abs(a2));
  if (err  > tol) {
    out << "\nError, relErr(" << a1_name <<","
        << a2_name << ") = relErr(" << a1 <<"," << a2 <<") = "
        << err << " <= tol = " << tol << ": failed!\n";
    success = false;
  }

  return success;
}

template <typename VectorType, typename MatrixType>
bool check_assembly(const VectorType& analytic_residual,
                    const MatrixType& analytic_jacobian,
                    const VectorType& fad_residual,
                    const MatrixType& fad_jacobian,
                    const std::string& test_name)
{
  const double tol = 1e-14;
  bool success = true;
  Teuchos::RCP<Teuchos::FancyOStream> out =
    Teuchos::VerboseObjectBase::getDefaultOStream();
  std::stringstream buf;
  Teuchos::FancyOStream fbuf(Teuchos::rcp(&buf,false));

  typename VectorType::HostMirror host_analytic_residual =
    Kokkos::create_mirror_view(analytic_residual);
  typename VectorType::HostMirror host_fad_residual =
    Kokkos::create_mirror_view(fad_residual);
  Kokkos::deep_copy( host_analytic_residual, analytic_residual );
  Kokkos::deep_copy( host_fad_residual, fad_residual );

  fbuf << test_name << ":" << std::endl;

  if (host_analytic_residual.extent(0) != host_fad_residual.extent(0)) {
    fbuf << "Analytic residual dimension "
         << host_analytic_residual.extent(0)
         << " does not match Fad residual dimension "
         << host_fad_residual.extent(0) << std::endl;
    success = false;
  }
  else {
    const size_t num_node = host_analytic_residual.extent(0);
    for (size_t i=0; i<num_node; ++i) {
      success = success && compareValues(
        host_analytic_residual(i), "analytic residual",
        host_fad_residual(i), "Fad residual",
        tol, tol, fbuf );
    }
  }

  typename MatrixType::HostMirror host_analytic_jacobian =
    Kokkos::create_mirror_view(analytic_jacobian);
  typename MatrixType::HostMirror host_fad_jacobian =
    Kokkos::create_mirror_view(fad_jacobian);
  Kokkos::deep_copy( host_analytic_jacobian, analytic_jacobian );
  Kokkos::deep_copy( host_fad_jacobian, fad_jacobian );

  if (host_analytic_jacobian.extent(0) != host_fad_jacobian.extent(0)) {
    fbuf << "Analytic Jacobian dimension "
         << host_analytic_jacobian.extent(0)
         << " does not match Fad Jacobian dimension "
         << host_fad_jacobian.extent(0) << std::endl;
    success = false;
  }
  else {
    const size_t num_entry = host_analytic_jacobian.extent(0);
    for (size_t i=0; i<num_entry; ++i) {
      success = success && compareValues(
        host_analytic_jacobian(i), "analytic Jacobian",
        host_fad_jacobian(i), "Fad Jacobian",
        tol, tol, fbuf );
    }
  }

  if (!success)
    *out << buf.str();

  return success;
}

template <class Device>
void performance_test_driver(
  const int use_print ,
  const int use_trials ,
  const int n_begin ,
  const int n_end ,
  const int n_step ,
  const bool quadratic ,
  const bool check )
{
  using Kokkos::Example::BoxElemPart;
  using Kokkos::Example::FENL::Analytic;
  using Kokkos::Example::FENL::FadElement;
  using Kokkos::Example::FENL::FadElementOptimized;
  using Kokkos::Example::FENL::FadQuadPoint;

  std::cout.precision(8);
  std::cout << std::endl
            << "\"Grid Size\" , "
            << "\"FEM Size\" , "
            << "\"Analytic Fill Time\" , "
            << "\"Fad Element Fill Slowdown\" , "
            << "\"Fad Optimized Element Fill Slowdown\" , "
            << "\"Fad QP Fill Slowdown\" , "
            << std::endl;

  typedef Kokkos::View< double* , Device > vector_type ;
  typedef Kokkos::Example::FENL::CrsMatrix<double,Device> matrix_type;
  vector_type analytic_residual, fad_residual, fad_opt_residual,
    fad_qp_residual;
  matrix_type analytic_jacobian, fad_jacobian, fad_opt_jacobian,
    fad_qp_jacobian;

  for (int n=n_begin; n<=n_end; n+=n_step) {
    const int use_nodes[] = { n, n, n };
    Perf perf_analytic, perf_fad, perf_fad_opt, perf_fad_qp;

    if (quadratic) {
      perf_analytic =
        fenl_assembly<double,Device,BoxElemPart::ElemQuadratic,Analytic>(
          use_print, use_trials, use_nodes,
          analytic_residual, analytic_jacobian );

      perf_fad =
        fenl_assembly<double,Device,BoxElemPart::ElemQuadratic,FadElement>(
          use_print, use_trials, use_nodes,
          fad_residual, fad_jacobian);

      perf_fad_opt =
        fenl_assembly<double,Device,BoxElemPart::ElemQuadratic,FadElementOptimized>(
          use_print, use_trials, use_nodes,
          fad_opt_residual, fad_opt_jacobian);

      perf_fad_qp =
        fenl_assembly<double,Device,BoxElemPart::ElemQuadratic,FadQuadPoint>(
          use_print, use_trials, use_nodes,
          fad_qp_residual, fad_qp_jacobian);
    }
    else {
      perf_analytic =
        fenl_assembly<double,Device,BoxElemPart::ElemLinear,Analytic>(
          use_print, use_trials, use_nodes,
          analytic_residual, analytic_jacobian );

      perf_fad =
        fenl_assembly<double,Device,BoxElemPart::ElemLinear,FadElement>(
          use_print, use_trials, use_nodes,
          fad_residual, fad_jacobian);

      perf_fad_opt =
        fenl_assembly<double,Device,BoxElemPart::ElemLinear,FadElementOptimized>(
          use_print, use_trials, use_nodes,
          fad_opt_residual, fad_opt_jacobian);

      perf_fad_qp =
        fenl_assembly<double,Device,BoxElemPart::ElemLinear,FadQuadPoint>(
          use_print, use_trials, use_nodes,
          fad_qp_residual, fad_qp_jacobian);
    }
    if (check) {
      check_assembly( analytic_residual, analytic_jacobian.coeff,
                      fad_residual, fad_jacobian.coeff,
                      "Fad" );
      check_assembly( analytic_residual, analytic_jacobian.coeff,
                      fad_opt_residual, fad_opt_jacobian.coeff,
                      "Optimized Fad" );
      check_assembly( analytic_residual, analytic_jacobian.coeff,
                      fad_qp_residual, fad_qp_jacobian.coeff,
                      "QP Fad" );
    }

    double s =
      1000.0 / ( use_trials * perf_analytic.global_elem_count );
    perf_analytic.scale(s);
    perf_fad.scale(s);
    perf_fad_opt.scale(s);
    perf_fad_qp.scale(s);

    std::cout.precision(3);
    std::cout << n << " , "
              << perf_analytic.global_node_count << " , "
              << std::setw(2)
              << std::scientific
              << perf_analytic.fill_time << " , "
              << std::fixed << std::setw(6)
              << perf_fad.fill_time / perf_analytic.fill_time << " , "
              << perf_fad_opt.fill_time / perf_analytic.fill_time << " , "
              << perf_fad_qp.fill_time / perf_analytic.fill_time << " , "
              << std::endl;
  }
}
