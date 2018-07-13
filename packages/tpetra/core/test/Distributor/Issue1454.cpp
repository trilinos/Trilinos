// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Chris Luchini cbluchi@sandia.gov
//
// ************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Distributor.hpp"
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
  typedef Tpetra::Map<> map_type;
  typedef map_type::device_type device_type;
  // mfh 01 Aug 2017: Deal with fix for #1088, by not using
  // Kokkos::CudaUVMSpace for communication buffers.
#ifdef KOKKOS_ENABLE_CUDA
  typedef typename std::conditional<
  std::is_same<typename device_type::execution_space, Kokkos::Cuda>::value,
    Kokkos::CudaSpace,
    typename device_type::memory_space>::type buffer_memory_space;
#else
  typedef typename device_type::memory_space buffer_memory_space;
#endif // KOKKOS_ENABLE_CUDA
  typedef typename device_type::execution_space buffer_execution_space;
  typedef Kokkos::Device<buffer_execution_space, buffer_memory_space> buffer_device_type;

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
    typedef TpetraTest::Functor1<decltype (proc_ids) > functor_type;
    Kokkos::parallel_for ("fill_proc_ids",
                          Kokkos::RangePolicy<ExecutionSpace> (0, n),
                          functor_type (proc_ids, comm_size));
  }
  Kokkos::fence ();
  auto proc_ids_host = Kokkos::create_mirror_view (proc_ids);
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

  distributor.doPostsAndWaits (exports, 1, imports);
  auto imports_host = Kokkos::create_mirror_view (imports);
  Kokkos::deep_copy (imports_host, imports);

  for (int i = 0; i < n_imports; ++i) {
    TEUCHOS_ASSERT_EQUALITY( imports_host(i), i / 3 );
  }
}

} // namespace (anonymous)


