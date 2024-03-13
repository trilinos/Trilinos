/*
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

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