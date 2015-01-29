// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
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
  using Teuchos::ParameterList;

  typedef Kokkos::Example::BoxElemFixture< Device , Kokkos::Example::BoxElemPart::ElemLinear > FixtureType ;

  typedef Kokkos::Example::FENL::CrsMatrix< Scalar , Device > LocalMatrixType ;
  typedef typename LocalMatrixType::StaticCrsGraphType LocalGraphType ;

  typedef Kokkos::Example::FENL::NodeNodeGraph< typename FixtureType::elem_node_type , LocalGraphType , FixtureType::ElemNode > NodeNodeGraphType ;

  typedef Kokkos::Example::FENL::ElementComputationConstantCoefficient CoeffFunctionType;
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
      mesh_to_graph( fixture.elem_node() , fixture.node_count_owned(),
                     graph_times );

    // Create the local sparse matrix from the graph:
    jacobian = LocalMatrixType( mesh_to_graph.graph );

    //----------------------------------

    // Allocate solution vector for each node in the mesh and residual vector for each owned node
    VectorType solution( "solution" , fixture.node_count() );
    residual = VectorType( "residual" , fixture.node_count_owned() );

    // Create element computation functor
    CoeffFunctionType diffusion_coefficient( 1.0 );
    const ElementComputationType elemcomp( fixture , diffusion_coefficient ,
                                           solution ,
                                           mesh_to_graph.elem_graph ,
                                           jacobian , residual );

    Kokkos::deep_copy( solution , Scalar(1.2345) );

    //--------------------------------
    // Element contributions to residual and jacobian

    Kokkos::deep_copy( residual , Scalar(0) );
    Kokkos::deep_copy( jacobian.coeff , Scalar(0) );

    wall_clock.reset();

    elemcomp.apply();

    Device::fence();
    perf.fill_time = wall_clock.seconds();

    //--------------------------------

    perf_stats.increment(perf);

  }

  return perf_stats ;
}

template <typename VectorType, typename MatrixType>
bool check_assembly(const VectorType& analytic_residual,
                    const MatrixType& analytic_jacobian,
                    const VectorType& fad_residual,
                    const MatrixType& fad_jacobian)
{
  const double tol = 1e-10;
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

  TEUCHOS_TEST_EQUALITY( host_analytic_residual.dimension_0(),
                         host_fad_residual.dimension_0(), fbuf, success );

  const size_t num_node = host_analytic_residual.dimension_0();
  for (size_t i=0; i<num_node; ++i) {
    TEUCHOS_TEST_FLOATING_EQUALITY(
      host_analytic_residual(i), host_fad_residual(i), tol, fbuf, success );
  }

  typename MatrixType::HostMirror host_analytic_jacobian =
    Kokkos::create_mirror_view(analytic_jacobian);
  typename MatrixType::HostMirror host_fad_jacobian =
    Kokkos::create_mirror_view(fad_jacobian);
  Kokkos::deep_copy( host_analytic_jacobian, analytic_jacobian );
  Kokkos::deep_copy( host_fad_jacobian, fad_jacobian );

  TEUCHOS_TEST_EQUALITY( host_analytic_jacobian.dimension_0(),
                         host_fad_jacobian.dimension_0(), fbuf, success );

  const size_t num_entry = host_analytic_jacobian.dimension_0();
  for (size_t i=0; i<num_entry; ++i) {
    TEUCHOS_TEST_FLOATING_EQUALITY(
      host_analytic_jacobian(i), host_fad_jacobian(i), tol, fbuf, success );
  }

  if (!success)
    *out << buf.str();

  return success;
}

template <class Device>
void performance_test_driver(
  const int use_print ,
  const int use_trials ,
  const int use_nodes[] ,
  const bool use_view ,
  const bool use_global ,
  const bool check )
{
  std::cout.precision(8);
  std::cout << std::endl;
  if (use_global)
    std::cout << "Use Global ";
  else
    std::cout << "Use Local ";
  if (use_view)
    std::cout << "View ";
  std::cout << "Assembly ";
  std::cout << std::endl
            << "\"Grid Size\" , "
            << "\"FEM Size\" , "
            << "\"Analytic Fill Time\" , "
            << "\"Fad Fill Time\" , "
            << "\"Fad Fill Slowdown\" , "
            << std::endl;

  typedef Kokkos::View< double* , Device > vector_type ;
  typedef Kokkos::Example::FENL::CrsMatrix<double,Device> matrix_type;
  vector_type analytic_residual, fad_residual;
  matrix_type analytic_jacobian, fad_jacobian;
  Perf perf_analytic, perf_fad;
  if (use_view) {
    if (use_global) {
      perf_analytic =
        fenl_assembly<double,Device,Kokkos::Example::FENL::AnalyticGlobalView>(
          use_print, use_trials, use_nodes,
          analytic_residual, analytic_jacobian );

      perf_fad =
        fenl_assembly<double,Device,Kokkos::Example::FENL::FadGlobalView>(
          use_print, use_trials, use_nodes,
          fad_residual, fad_jacobian);
    }
    else {
      perf_analytic =
        fenl_assembly<double,Device,Kokkos::Example::FENL::AnalyticLocalView>(
          use_print, use_trials, use_nodes,
          analytic_residual, analytic_jacobian );

      perf_fad =
        fenl_assembly<double,Device,Kokkos::Example::FENL::FadLocalView>(
          use_print, use_trials, use_nodes,
          fad_residual, fad_jacobian);
    }
  }
  else {
    perf_analytic =
      fenl_assembly<double,Device,Kokkos::Example::FENL::AnalyticLocal>(
        use_print, use_trials, use_nodes,
        analytic_residual, analytic_jacobian );

    perf_fad =
      fenl_assembly<double,Device,Kokkos::Example::FENL::FadLocal>(
        use_print, use_trials, use_nodes,
        fad_residual, fad_jacobian);
  }

  if (check)
    check_assembly( analytic_residual, analytic_jacobian.coeff,
                    fad_residual, fad_jacobian.coeff );

  double s =
    1000.0 / ( use_trials * perf_analytic.global_node_count );
  perf_analytic.scale(s);
  perf_fad.scale(s);

  std::cout.precision(3);
  std::cout << use_nodes[0] << " , "
            << perf_analytic.global_node_count << " , "
            << std::setw(2)
            << std::scientific
            << perf_analytic.fill_time << " , "
            << perf_fad.fill_time << " , "
            << std::fixed << std::setw(6)
            << perf_fad.fill_time / perf_analytic.fill_time << " , "
            << std::endl;
}
