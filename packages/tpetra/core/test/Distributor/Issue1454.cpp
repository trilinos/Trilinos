// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Distributor.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Kokkos_Core.hpp"
#include <type_traits>

// Kokkos sometimes doesn't like functors in an anonymous namespace.
namespace TpetraTest {

template<class ViewType>
class Functor1 {
public:
  Functor1 (const ViewType& proc_ids, const int comm_size) :
    proc_ids_ (proc_ids),
    comm_size_ (comm_size)
  {}

  KOKKOS_FUNCTION void operator() (const int& i) const {
    proc_ids_(i) = i % comm_size_;
  }

private:
  ViewType proc_ids_;
  int comm_size_;
};

template<class ViewType>
class Functor2 {
public:
  Functor2 (const ViewType& exports, const int comm_rank) :
    exports_ (exports),
    comm_rank_ (comm_rank)
  {}

  KOKKOS_FUNCTION void operator() (const int& i) const {
    exports_(i) = comm_rank_;
  }

private:
  ViewType exports_;
  int comm_rank_;
};

} // namespace TpetraTest

namespace { // (anonymous)

TEUCHOS_UNIT_TEST( Distributor, Issue1454 )
{
  using map_type = Tpetra::Map<>;
  using device_type = map_type::device_type;
  // mfh 01 Aug 2017: Deal with fix for #1088, by not using
  // Kokkos::CudaUVMSpace for communication buffers.
#ifdef KOKKOS_ENABLE_CUDA
  using buffer_memory_space = typename std::conditional<
  std::is_same<typename device_type::execution_space, Kokkos::Cuda>::value,
    Kokkos::CudaSpace,
    typename device_type::memory_space>::type;
#elif defined(KOKKOS_ENABLE_SYCL)
  using buffer_memory_space = typename std::conditional<
  std::is_same<typename device_type::execution_space, Kokkos::Experimental::SYCL>::value,
    Kokkos::Experimental::SYCLDeviceUSMSpace,
    typename device_type::memory_space>::type;
#else
  using buffer_memory_space = typename device_type::memory_space;
#endif
  using buffer_execution_space = typename device_type::execution_space;
  using buffer_device_type = Kokkos::Device<buffer_execution_space, buffer_memory_space>;

  auto comm = Tpetra::TestingUtilities::getDefaultComm ();
  // Create a Map just to ensure that Kokkos gets initialized and
  // finalized correctly.
  const map_type map (comm->getSize (), 1, 0, comm);

  const int comm_rank = comm->getRank();
  const int comm_size = comm->getSize();

  Tpetra::Distributor distributor( comm );
  const int n = 3 * comm_size;
  Kokkos::View<int *, device_type> proc_ids( "proc_ids", n );
  const int n_exports = proc_ids.extent( 0 );
  using ExecutionSpace = typename device_type::execution_space;
  {
    using functor_type = TpetraTest::Functor1<decltype (proc_ids) >;
    Kokkos::parallel_for ("fill_proc_ids",
                          Kokkos::RangePolicy<ExecutionSpace> (0, n),
                          functor_type (proc_ids, comm_size));
  }
  Kokkos::fence ();
  auto proc_ids_host = Kokkos::create_mirror_view (Kokkos::HostSpace (), proc_ids);
  static_assert (std::is_same<typename decltype (proc_ids_host)::memory_space,
                   Kokkos::HostSpace>::value,
                 "proc_ids_host should be a HostSpace View, but is not.");
  Kokkos::deep_copy (proc_ids_host, proc_ids);
  const int n_imports =
    distributor.createFromSends (Teuchos::ArrayView<const int> (proc_ids_host.data (), n_exports));
  Kokkos::View<int *, buffer_device_type> exports( "exports", n_exports );
  {
    typedef TpetraTest::Functor2<decltype (exports) > functor_type;
    Kokkos::parallel_for ("fill_exports",
                          Kokkos::RangePolicy<buffer_execution_space> (0, n_exports),
                          functor_type (exports, comm_rank));
  }
  Kokkos::fence();

  Kokkos::View<int *, buffer_device_type> imports( "imports", n_imports );
  auto imports_host = Kokkos::create_mirror_view (imports);
  // This assertion fails and is meaningless for HIP right now which uses HostPinnedSpace here
  #ifndef KOKKOS_ENABLE_HIP
  static_assert (std::is_same<typename decltype (imports_host)::memory_space,
                   Kokkos::HostSpace>::value,
                 "imports_host should be a HostSpace View, but is not.");
  #endif
  if (Tpetra::Details::Behavior::assumeMpiIsGPUAware ()) {
    distributor.doPostsAndWaits (exports, 1, imports);
    Kokkos::deep_copy (imports_host, imports);
  }
  else {
    Kokkos::deep_copy (imports_host, imports);
    auto exports_host = Kokkos::create_mirror_view (exports);
    // This assertion fails and is meaningless for HIP right now which uses HostPinnedSpace here
    #ifndef KOKKOS_ENABLE_HIP
    static_assert (std::is_same<typename decltype (exports_host)::memory_space,
                   Kokkos::HostSpace>::value,
      "exports_host should be a HostSpace View, but is not.");
    #endif
    Kokkos::deep_copy (exports_host, exports);
    distributor.doPostsAndWaits (exports_host, 1, imports_host);
  }

  for (int i = 0; i < n_imports; ++i) {
    TEUCHOS_ASSERT_EQUALITY( imports_host(i), i / 3 );
  }
}

} // namespace (anonymous)


