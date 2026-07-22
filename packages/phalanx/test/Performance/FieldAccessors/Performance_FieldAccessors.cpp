// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Phalanx_KokkosViewFactory.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_MDField.hpp"

// From test/Utilities directory
#include "Traits.hpp"

/*! \brief Test to check performance of PHX::MDField vs Kokkos::View

  This test is to show that there is no performance penalty between a
  Kokkos::View and a PHX::MDField.  We wrap a Kokkos::View in the
  MDField.  The compiler should optimize away the wrapper code fro
  optimal performance.  The important comparison is the Compiletime
  MDField and the non-static Kokkos::View<double***>. At optimization
  of -O3 on gcc 4.8.1, we see no difference in runtime.

  NOTE: The actual array sizes used for timings are commented out to
  make the problem smaller.  This is so that this unit test will not
  timeout in debug builds where the executable is much slower due to
  bounds checking.  look for the declaration of the num_loops,
  num_cells, and the static Kokkos::View to change the problem size.
  These declarations are grouped together for convenience.
*/

PHX_EXTENT(P)

using size_type = PHX::exec_space::size_type;

template <typename Scalar,typename Device,typename Array>
class ComputeA {
  Array a_;
  Array b_;
  Array c_;
public:
  typedef PHX::Device execution_space;

  ComputeA(Array& a,Array& b,Array& c)
  : a_(a), b_(b), c_(c)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator () (const size_type c) const
  {
    const size_type num_ip = a_.extent(1);
    const size_type num_dim = a_.extent(2);
    for (size_type i = 0; i < num_ip; ++i) {
      for (size_type d = 0; d < num_dim; ++d) {
	a_(c,i,d) =  b_(c,i,d) * c_(c,i,d) + b_(c,i,d) + b_(c,i,d) / c_(c,i,d);
      }
    }
  }
};

TEUCHOS_UNIT_TEST(performance, ArrayAccessor)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  RCP<Time> phx_ct_time_pf = TimeMonitor::getNewTimer("MDField Compiletime Rank (device="+PHX::print<PHX::Device>()+")");
  RCP<Time> phx_rt_time_pf = TimeMonitor::getNewTimer("MDField Runtime Rank (device="+PHX::print<PHX::Device>()+")");
  RCP<Time> k_time_pf = TimeMonitor::getNewTimer("KokkosView<double***>(device="+PHX::print<PHX::Device>()+")");
  RCP<Time> k_time_pf_static = TimeMonitor::getNewTimer("KokkosView<double[][][]>(device="+PHX::print<PHX::Device>()+")");

  std::cout << std::endl << std::endl
            << "PHX::Device::size_type = "
            << PHX::print<PHX::exec_space::size_type>()
            << std::endl;

  std::cout << "PHX::index_size_type = "
            << PHX::print<PHX::index_size_type>()
            << std::endl;

  // For performance testing, build in RELEASE mode and use the
  // numbers in the comments below.  We have small numbers by
  // default for unit testing in debug mode where the code is much
  // slower.
  // const size_type num_loops = 1;
  // const size_type num_cells = 10000000;
  // using kokkos_field_static = Kokkos::View<double[10000000][25][4],PHX::MemSpace>;
  const size_type num_loops = 1;
  const size_type num_cells = 100;
  using kokkos_field_static = Kokkos::View<double[100][25][4],PHX::MemSpace>;
  const size_type num_ip = 25;
  const size_type num_dim = 4;

  // Compiletime PHX:MDField
  {
    using phx_ct_field = MDField<double,P,P,P>;
    phx_ct_field a("phx_ct_a","qp_layout",num_cells,num_ip,num_dim);
    phx_ct_field b("phx_ct_b","qp_layout",num_cells,num_ip,num_dim);
    phx_ct_field c("phx_ct_c","qp_layout",num_cells,num_ip,num_dim);
    Kokkos::deep_copy(a.get_static_view(),1.0);
    Kokkos::deep_copy(b.get_static_view(),2.0);
    Kokkos::deep_copy(c.get_static_view(),3.0);

    cout << "MDField Compiletime Rank" << endl;
    TimeMonitor tm(*phx_ct_time_pf);
    for (size_type l=0; l < num_loops; ++l) {
      Kokkos::parallel_for("MDField (Compiletime Rank)",num_cells,ComputeA<double,PHX::ExecSpace,phx_ct_field> (a,b,c));
      typename PHX::Device().fence();
    }
  }

  // Runtime PHX::MDField
  {
    using phx_rt_field = MDField<double>;
    phx_rt_field a("phx_rt_a","qp_layout",num_cells,num_ip,num_dim);
    phx_rt_field b("phx_rt_b","qp_layout",num_cells,num_ip,num_dim);
    phx_rt_field c("phx_rt_c","qp_layout",num_cells,num_ip,num_dim);
    Kokkos::deep_copy(a.get_view(),1.0);
    Kokkos::deep_copy(b.get_view(),2.0);
    Kokkos::deep_copy(c.get_view(),3.0);
    cout << "MDField Runtime Rank" << endl;
    TimeMonitor tm(*phx_rt_time_pf);
    for (size_type l=0; l < num_loops; ++l) {
      Kokkos::parallel_for("MDField (Runtime Rank)",num_cells,ComputeA<double,PHX::ExecSpace,phx_rt_field> (a,b,c));
      typename PHX::Device().fence();
    }
  }

  // Kokkos View
  {
    using kokkos_field = Kokkos::View<double***,PHX::Device>;
    kokkos_field a("kv_a",num_cells,num_ip,num_dim);
    kokkos_field b("kv_b",num_cells,num_ip,num_dim);
    kokkos_field c("kv_c",num_cells,num_ip,num_dim);
    Kokkos::deep_copy(a,1.0);
    Kokkos::deep_copy(b,2.0);
    Kokkos::deep_copy(c,3.0);
    cout << "Kokkos View" << endl;
    TimeMonitor tm(*k_time_pf);
    for (size_type l=0; l < num_loops; ++l)
      Kokkos::parallel_for("Kokkos View",num_cells,ComputeA<double,PHX::Device,kokkos_field>(a,b,c));
    typename PHX::Device().fence();
  }

  // Static Kokkos View
  {
    kokkos_field_static a("kvs_a");
    kokkos_field_static b("kvs_b");
    kokkos_field_static c("kvs_c");
    // Check that static size is consistent with runtime
    TEST_EQUALITY(num_cells,size_type(a.extent(0)));
    Kokkos::deep_copy(a,1.0);
    Kokkos::deep_copy(b,2.0);
    Kokkos::deep_copy(c,3.0);
    cout << "Static Kokkos View" << endl;
    TimeMonitor tm(*k_time_pf_static);
    for (size_type l=0; l < num_loops; ++l)
      Kokkos::parallel_for("Kokkos Static View",num_cells,ComputeA<double,PHX::Device,kokkos_field_static>(a,b,c));
    typename PHX::Device().fence();
  }

  TimeMonitor::summarize();
}
