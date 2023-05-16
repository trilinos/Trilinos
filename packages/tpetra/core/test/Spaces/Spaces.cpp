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
#include <Tpetra_Details_ExecutionSpaces.hpp>

/*! \file Spaces.cpp 
    \brief Make sure Tpetra::Details::Spaces interfaces compile and don't crash
*/

namespace { // (anonymous)

template <typename ExecutionSpace>
void test_exec_space_wait(bool &success, Teuchos::FancyOStream &out) {
  ExecutionSpace e1, e2;
  Tpetra::Details::Spaces::exec_space_wait(e1, e2);
  Tpetra::Details::Spaces::exec_space_wait("test_exec_space_wait", e1, e2);

  success = true;
}

template <typename ExecutionSpace>
void test_make_instance(bool &success, Teuchos::FancyOStream &out) {

  using Priority = Tpetra::Details::Spaces::Priority;
  {
    ExecutionSpace e1 =
        Tpetra::Details::Spaces::make_instance<ExecutionSpace>();
    ExecutionSpace e2 =
        Tpetra::Details::Spaces::make_instance<ExecutionSpace, Priority::low>();
    ExecutionSpace e3 =
        Tpetra::Details::Spaces::make_instance<ExecutionSpace,
                                               Priority::medium>();
    ExecutionSpace e4 =
        Tpetra::Details::Spaces::make_instance<ExecutionSpace,
                                               Priority::high>();
  }
  {
    ExecutionSpace e1 =
        Tpetra::Details::Spaces::make_instance<ExecutionSpace>(Priority::low);
    ExecutionSpace e2 = Tpetra::Details::Spaces::make_instance<ExecutionSpace>(
        Priority::medium);
    ExecutionSpace e3 =
        Tpetra::Details::Spaces::make_instance<ExecutionSpace>(Priority::high);
  }

  success = true;
}

template <typename ExecutionSpace>
void test_space_instance(bool &success, Teuchos::FancyOStream &out) {
  using Priority = Tpetra::Details::Spaces::Priority;

  // unnumbered spaces are the same
  {
    Teuchos::RCP<const ExecutionSpace> e1 =
        Tpetra::Details::Spaces::space_instance<ExecutionSpace>();
    Teuchos::RCP<const ExecutionSpace> e2 =
        Tpetra::Details::Spaces::space_instance<ExecutionSpace>();
    TEUCHOS_TEST_EQUALITY(e1, e2, out, success);
  }

  // unnumbered space is the same as 0
  {
    Teuchos::RCP<const ExecutionSpace> e1 =
        Tpetra::Details::Spaces::space_instance<ExecutionSpace>();
    Teuchos::RCP<const ExecutionSpace> e2 =
        Tpetra::Details::Spaces::space_instance<ExecutionSpace>(0);
    TEUCHOS_TEST_EQUALITY(e1, e2, out, success);
  }

  // construct priorities (no relationship with each other defined yet)
  {
    Teuchos::RCP<const ExecutionSpace> e1 =
        Tpetra::Details::Spaces::space_instance<ExecutionSpace,
                                                Priority::low>();
    Teuchos::RCP<const ExecutionSpace> e2 =
        Tpetra::Details::Spaces::space_instance<ExecutionSpace,
                                                Priority::medium>();
    Teuchos::RCP<const ExecutionSpace> e3 =
        Tpetra::Details::Spaces::space_instance<ExecutionSpace,
                                                Priority::high>();
  }
  {
    Teuchos::RCP<const ExecutionSpace> e1 =
        Tpetra::Details::Spaces::space_instance<ExecutionSpace, Priority::low>(
            0);
    Teuchos::RCP<const ExecutionSpace> e2 =
        Tpetra::Details::Spaces::space_instance<ExecutionSpace,
                                                Priority::medium>(2);
    Teuchos::RCP<const ExecutionSpace> e3 =
        Tpetra::Details::Spaces::space_instance<ExecutionSpace, Priority::high>(
            10);
  }

  {
    Teuchos::RCP<const ExecutionSpace> e1 =
        Tpetra::Details::Spaces::space_instance<ExecutionSpace>(Priority::low);
    Teuchos::RCP<const ExecutionSpace> e2 =
        Tpetra::Details::Spaces::space_instance<ExecutionSpace>(
            Priority::medium);
    Teuchos::RCP<const ExecutionSpace> e3 =
        Tpetra::Details::Spaces::space_instance<ExecutionSpace>(Priority::high);
  }

  {
    Teuchos::RCP<const ExecutionSpace> e1 =
        Tpetra::Details::Spaces::space_instance<ExecutionSpace>(Priority::low,
                                                                0);
    Teuchos::RCP<const ExecutionSpace> e2 =
        Tpetra::Details::Spaces::space_instance<ExecutionSpace>(
            Priority::medium, 2);
    Teuchos::RCP<const ExecutionSpace> e3 =
        Tpetra::Details::Spaces::space_instance<ExecutionSpace>(Priority::high,
                                                                10);
  }

  success = true;
}

template <typename ExecutionSpace, bool expected>
void test_is_gpu_exec_space(bool &success, Teuchos::FancyOStream &out) {
  TEUCHOS_TEST_EQUALITY(
      expected, Tpetra::Details::Spaces::is_gpu_exec_space<ExecutionSpace>(),
      out, success);
}



} // namespace

int main(int argc, char **argv) {
  Tpetra::ScopeGuard sg(&argc, &argv);

  bool success = false;
  auto out = Teuchos::fancyOStream (Teuchos::rcpFromRef (std::cout));
  *out << "Test spaces" << std::endl;

#if defined(KOKKOS_ENABLE_SERIAL)
  test_exec_space_wait<Kokkos::Serial>(success, *out);
  test_make_instance<Kokkos::Serial>(success, *out);
  test_space_instance<Kokkos::Serial>(success, *out);
  test_is_gpu_exec_space<Kokkos::Serial, false>(success, *out);
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
  test_exec_space_wait<Kokkos::OpenMP>(success, *out);
  test_make_instance<Kokkos::OpenMP>(success, *out);
  test_space_instance<Kokkos::OpenMP>(success, *out);
  test_is_gpu_exec_space<Kokkos::OpenMP, false>(success, *out);
#endif

#if defined(KOKKOS_ENABLE_THREADS)
  test_exec_space_wait<Kokkos::Threads>(success, *out);
  test_make_instance<Kokkos::Threads>(success, *out);
  test_space_instance<Kokkos::Threads>(success, *out);
  test_is_gpu_exec_space<Kokkos::Threads, false>(success, *out);
#endif

#if defined(KOKKOS_ENABLE_CUDA)
  test_exec_space_wait<Kokkos::Cuda>(success, *out);
  test_make_instance<Kokkos::Cuda>(success, *out);
  test_space_instance<Kokkos::Cuda>(success, *out);
  test_is_gpu_exec_space<Kokkos::Cuda, true>(success, *out);
#endif

#if defined(KOKKOS_ENABLE_HIP)
  test_exec_space_wait<Kokkos::HIP>(success, *out);
  test_make_instance<Kokkos::HIP>(success, *out);
  test_space_instance<Kokkos::HIP>(success, *out);
  test_is_gpu_exec_space<Kokkos::HIP, true>(success, *out);
#endif

#if defined(KOKKOS_ENABLE_SYCL)
  test_exec_space_wait<Kokkos::Experimental::SYCL>(success, *out);
  test_make_instance<Kokkos::Experimental::SYCL>(success, *out);
  test_space_instance<Kokkos::Experimental::SYCL>(success, *out);
  test_is_gpu_exec_space<Kokkos::Experimental::SYCL, true>(success, *out);
#endif

  if (success) {
    std::cout << std::endl << "End Result: TEST PASSED" << std::endl;
  } else {
    std::cout << std::endl << "End Result: TEST FAILED" << std::endl;
  }
  
  Teuchos::OSTab tab1(out);
  return 0;
}