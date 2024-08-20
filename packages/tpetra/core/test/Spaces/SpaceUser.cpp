// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Details_ExecutionSpacesUser.hpp>

/*! \file SpaceUser.cpp 
    \brief Make sure SpaceUser compiles and doesn't crash
*/

namespace { // (anonymous)

template <typename ExecSpace> struct S : public Tpetra::Details::Spaces::User {

  static constexpr size_t B = size_t(1024) * size_t(1024) * size_t(1024);

  void priority() const {
    auto e1 =
        space_instance<ExecSpace, Tpetra::Details::Spaces::Priority::low>();
    auto e2 =
        space_instance<ExecSpace, Tpetra::Details::Spaces::Priority::medium>();
    auto e3 =
        space_instance<ExecSpace, Tpetra::Details::Spaces::Priority::high>();
    auto e4 = space_instance<ExecSpace>(Tpetra::Details::Spaces::Priority::low);
    auto e5 =
        space_instance<ExecSpace>(Tpetra::Details::Spaces::Priority::medium);
    auto e6 =
        space_instance<ExecSpace>(Tpetra::Details::Spaces::Priority::high);
  }

  // request the same instance 1B times should be pretty quick
  void reuse() const {
    for (size_t i = 0; i < B; ++i) {
      auto e1 =
          space_instance<ExecSpace, Tpetra::Details::Spaces::Priority::low>();
    }
  }

}; // S

template <typename ExecutionSpace>
void test_priority(bool &success, Teuchos::FancyOStream &out) {
  S<ExecutionSpace>().priority();
  success = true;
}

template <typename ExecutionSpace>
void test_reuse(bool &success, Teuchos::FancyOStream &out) {
  S<ExecutionSpace>().reuse();
  success = true;
}

} // namespace

int main(int argc, char **argv) {
  Tpetra::ScopeGuard sg(&argc, &argv);

  bool success = false;
  auto out = Teuchos::fancyOStream (Teuchos::rcpFromRef (std::cout));
  *out << "Test SpaceUser" << std::endl;

#if defined(KOKKOS_ENABLE_SERIAL)
  test_priority<Kokkos::Serial>(success, *out);
  test_reuse<Kokkos::Serial>(success, *out);
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
  test_priority<Kokkos::OpenMP>(success, *out);
  test_reuse<Kokkos::OpenMP>(success, *out);
#endif

#if defined(KOKKOS_ENABLE_THREADS)
  test_priority<Kokkos::Threads>(success, *out);
  test_reuse<Kokkos::Threads>(success, *out);
#endif

#if defined(KOKKOS_ENABLE_CUDA)
  test_priority<Kokkos::Cuda>(success, *out);
  test_reuse<Kokkos::Cuda>(success, *out);
#endif

#if defined(KOKKOS_ENABLE_HIP)
  test_priority<Kokkos::HIP>(success, *out);
  test_reuse<Kokkos::HIP>(success, *out);
#endif

#if defined(KOKKOS_ENABLE_SYCL)
  test_priority<Kokkos::Experimental::SYCL>(success, *out);
  test_reuse<Kokkos::Experimental::SYCL>(success, *out);
#endif

  if (success) {
    std::cout << std::endl << "End Result: TEST PASSED" << std::endl;
  } else {
    std::cout << std::endl << "End Result: TEST FAILED" << std::endl;
  }
  
  Teuchos::OSTab tab1(out);
  return 0;
}