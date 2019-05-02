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

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_withLocalAccess.hpp"
#include "Tpetra_withLocalAccess_MultiVector.hpp"
#include "Tpetra_for_each.hpp"
#include "Tpetra_for_each_MultiVector.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Kokkos_ArithTraits.hpp"

namespace { // (anonymous)

  using Tpetra::getDefaultComm;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::endl;
  using GST = Tpetra::global_size_t;
  using map_type = Tpetra::Map<>;
  using multivec_type = Tpetra::MultiVector<>;
  using vec_type = Tpetra::Vector<>;
  using GO = map_type::global_ordinal_type;

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST( VectorHarness, GetLocalObject )
  {
    using Tpetra::Details::getMasterLocalObject;
    using Tpetra::Details::getNonowningLocalObject;
    using Tpetra::readOnly;
    using Tpetra::readWrite;
    using Tpetra::writeOnly;
    const bool debug = ::Tpetra::Details::Behavior::debug ();

    RCP<Teuchos::FancyOStream> outPtr = debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Test Tpetra::Details::{getMasterLocalObject, "
      "getNonowningLocalObject} for MultiVector and Vector" << endl;
    Teuchos::OSTab tab0 (myOut);

    myOut << "Create a Map" << endl;
    auto comm = getDefaultComm ();
    const auto INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const size_t numLocal = 13;
    const size_t numVecs  = 3;
    const GO indexBase = 0;
    auto map = rcp (new map_type (INVALID, numLocal, indexBase, comm));

    myOut << "Create a MultiVector, and make sure that it has "
      "the right number of vectors (columns)" << endl;
    multivec_type mvec (map, numVecs);
    TEST_EQUALITY( mvec.getNumVectors (), numVecs );

    myOut << "Create a Vector, and make sure that "
      "it has exactly one vector (column)" << endl;
    vec_type vec (map);
    TEST_EQUALITY_CONST(vec.getNumVectors (), size_t (1));

    // Test read-only nonowning MultiVector access.
    {
      auto X_lcl_ro_owning = getMasterLocalObject (readOnly (mvec));
      auto X_lcl_ro = getNonowningLocalObject (X_lcl_ro_owning);
      Kokkos::View<const double**,
                   Kokkos::LayoutLeft,
                   multivec_type::device_type,
                   Kokkos::MemoryUnmanaged> X_lcl_ro2 = X_lcl_ro;
      // mfh 30 Apr 2019: Commented out due to possible NVCC bug.
      //static_assert (decltype (X_lcl_ro)::Rank == 2, "Rank is not 2");
      TEST_ASSERT( size_t (X_lcl_ro.extent (0)) == numLocal );
      TEST_ASSERT( size_t (X_lcl_ro.extent (1)) == numVecs );
    }

    // Test whether read-only access works with a const MultiVector&.
    {
      using MV = multivec_type;
      auto X_lcl_ro_owning = getMasterLocalObject (readOnly (mvec));
      auto X_lcl_ro = getNonowningLocalObject (X_lcl_ro_owning);
      // Make sure X_lcl_ro can be assigned to the type we expect it to
      // be.  It doesn't have to be that type, it just has to be
      // assignable to that type.
      Kokkos::View<const double**,
                   Kokkos::LayoutLeft,
                   multivec_type::device_type,
                   Kokkos::MemoryUnmanaged> X_lcl_ro2 = X_lcl_ro;
      // mfh 30 Apr 2019: Commented out due to possible NVCC bug.
      // Errors look like this:
      //
      // error: ‘__T0’ has not been declared
#ifndef KOKKOS_ENABLE_CUDA
      static_assert (decltype (X_lcl_ro)::Rank == 2, "Rank is not 2");
#endif // KOKKOS_ENABLE_CUDA
      TEST_ASSERT( size_t (X_lcl_ro.extent (0)) == numLocal );
      TEST_ASSERT( size_t (X_lcl_ro.extent (1)) == numVecs );
    }

    // Test write-only nonowning MultiVector access.
    {
      auto X_lcl_wo_owning = getMasterLocalObject (writeOnly (mvec));
      auto X_lcl_wo = getNonowningLocalObject (X_lcl_wo_owning);
      Kokkos::View<double**,
                   Kokkos::LayoutLeft,
                   multivec_type::device_type,
                   Kokkos::MemoryUnmanaged> X_lcl_wo2 = X_lcl_wo;
      // mfh 30 Apr 2019: Commented out due to possible NVCC bug.
#ifndef KOKKOS_ENABLE_CUDA
      static_assert (decltype (X_lcl_wo)::Rank == 2, "Rank is not 2");
#endif // KOKKOS_ENABLE_CUDA
      TEST_ASSERT( size_t (X_lcl_wo.extent (0)) == numLocal );
      TEST_ASSERT( size_t (X_lcl_wo.extent (1)) == numVecs );
    }

    // Test read-write nonowning MultiVector access.
    {
      auto X_lcl_rw_owning = getMasterLocalObject (readWrite (mvec));
      auto X_lcl_rw = getNonowningLocalObject (X_lcl_rw_owning);
      Kokkos::View<double**,
                   Kokkos::LayoutLeft,
                   multivec_type::device_type,
                   Kokkos::MemoryUnmanaged> X_lcl_rw2 = X_lcl_rw;
      // mfh 30 Apr 2019: Commented out due to possible NVCC bug.
#ifndef KOKKOS_ENABLE_CUDA
      static_assert (decltype (X_lcl_rw)::Rank == 2, "Rank is not 2");
#endif // KOKKOS_ENABLE_CUDA
      TEST_ASSERT( size_t (X_lcl_rw.extent (0)) == numLocal );
      TEST_ASSERT( size_t (X_lcl_rw.extent (1)) == numVecs );
    }

    // Test read-write nonowning Vector access.
    {
      auto X_lcl_1d_ro_owning = getMasterLocalObject (readOnly (vec));
      auto X_lcl_1d_ro = getNonowningLocalObject (X_lcl_1d_ro_owning);
      Kokkos::View<const double*,
                   Kokkos::LayoutLeft,
                   multivec_type::device_type,
                   Kokkos::MemoryUnmanaged> X_lcl_1d_ro2 = X_lcl_1d_ro;
      // mfh 30 Apr 2019: Commented out due to possible NVCC bug.
#ifndef KOKKOS_ENABLE_CUDA
      static_assert (decltype (X_lcl_1d_ro)::Rank == 1, "Rank is not 1");
#endif // KOKKOS_ENABLE_CUDA
      TEST_ASSERT( size_t (X_lcl_1d_ro.extent (0)) == numLocal );
    }

    // Test write-only nonowning Vector access.
    {
      auto X_lcl_1d_wo_owning = getMasterLocalObject (writeOnly (vec));
      auto X_lcl_1d_wo = getNonowningLocalObject (X_lcl_1d_wo_owning);
      Kokkos::View<double*,
                   Kokkos::LayoutLeft,
                   multivec_type::device_type,
                   Kokkos::MemoryUnmanaged> X_lcl_1d_wo2 = X_lcl_1d_wo;
      // mfh 30 Apr 2019: Commented out due to possible NVCC bug.
#ifndef KOKKOS_ENABLE_CUDA
      static_assert (decltype (X_lcl_1d_wo)::Rank == 1, "Rank is not 1");
#endif // KOKKOS_ENABLE_CUDA
      TEST_ASSERT( size_t (X_lcl_1d_wo.extent (0)) == numLocal );
    }

    // Test read-write nonowning Vector access.
    {
      auto X_lcl_1d_wr_owning = getMasterLocalObject (readWrite (vec));
      auto X_lcl_1d_wr = getNonowningLocalObject (X_lcl_1d_wr_owning);
      Kokkos::View<double*,
                   Kokkos::LayoutLeft,
                   multivec_type::device_type,
                   Kokkos::MemoryUnmanaged> X_lcl_1d_wr2 = X_lcl_1d_wr;
      // mfh 30 Apr 2019: Commented out due to possible NVCC bug.
#ifndef KOKKOS_ENABLE_CUDA
      static_assert (decltype (X_lcl_1d_wr)::Rank == 1, "Rank is not 1");
#endif // KOKKOS_ENABLE_CUDA
      TEST_ASSERT( size_t (X_lcl_1d_wr.extent (0)) == numLocal );
    }

    //
    // Examples of using the result of getNonowningLocalObject in
    // Kokkos::parallel_for kernels.
    //

    {
      using execution_space = vec_type::device_type::execution_space;
      using memory_space = vec_type::device_type::memory_space;
      using LO = vec_type::local_ordinal_type;
      using range_type = Kokkos::RangePolicy<execution_space, LO>;

      auto X_lcl_1d_wo_owning =
        getMasterLocalObject (writeOnly (vec).on (memory_space ()));
      auto X_lcl_1d_wo = getNonowningLocalObject (X_lcl_1d_wo_owning);
      static_assert
        (std::is_same<
           decltype (X_lcl_1d_wo)::device_type::execution_space,
           vec_type::dual_view_type::t_dev::execution_space>::value,
         "Wrong execution space");
      Kokkos::parallel_for (
        "Device kernel for write-only Tpetra::Vector",
        range_type (0, LO (numLocal)),
        KOKKOS_LAMBDA (const LO lclRow) {
          X_lcl_1d_wo(lclRow) = 42.0;
        });
    }

    {
      using host_execution_space =
        vec_type::dual_view_type::t_host::execution_space;
      using LO = vec_type::local_ordinal_type;
      using range_type = Kokkos::RangePolicy<host_execution_space, LO>;

      auto X_lcl_1d_wo_owning =
        getMasterLocalObject (writeOnly (vec).on (Kokkos::HostSpace ()));
      auto X_lcl_1d_wo = getNonowningLocalObject (X_lcl_1d_wo_owning);
      static_assert
        (std::is_same<
           decltype (X_lcl_1d_wo)::device_type::execution_space,
           vec_type::dual_view_type::t_host::execution_space>::value,
         "Wrong execution space");
      // test with some not-device function
      Kokkos::parallel_for (
        "Host kernel for write-only Tpetra::Vector",
        range_type (0, LO (numLocal)),
        [=] (const LO lclRow) {
          std::pair<double, double> p {3.0, 4.0};
          X_lcl_1d_wo(lclRow) = p.first * p.second;
        });

      // Just plain modify some entries, in some sequential order.
      // Just in case LO is unsigned, start at +1.
      for (LO lclRowPlusOne = LO (numLocal);
           lclRowPlusOne > LO (0); --lclRowPlusOne) {
        const LO lclRow = lclRowPlusOne - LO (1);
        // operator[] for 1-D Kokkos::View does the same thing as
        // operator().
        X_lcl_1d_wo[lclRow] = double (lclRow) + 42.0;
      }
    }

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess,
                         outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }

  TEUCHOS_UNIT_TEST( VectorHarness, WithLocalAccess )
  {
    using LO = typename multivec_type::local_ordinal_type;
    const bool debug = ::Tpetra::Details::Behavior::debug ();

    RCP<Teuchos::FancyOStream> outPtr = debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Test Tpetra::withLocalAccess" << endl;
    Teuchos::OSTab tab0 (myOut);

    myOut << "Create a Map" << endl;
    auto comm = getDefaultComm ();
    const auto INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const size_t numLocal = 13;
    const size_t numVecs  = 3;
    const GO indexBase = 0;
    auto map = rcp (new map_type (INVALID, numLocal, indexBase, comm));

    myOut << "Create a MultiVector, and make sure that it has "
      "the right number of vectors (columns)" << endl;
    multivec_type X (map, numVecs);
    TEST_EQUALITY( X.getNumVectors (), numVecs );

    using const_lcl_mv_type =
      Kokkos::View<const double**,
                   Kokkos::LayoutLeft,
                   typename multivec_type::device_type,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using nonconst_lcl_mv_type =
      Kokkos::View<double**,
                   Kokkos::LayoutLeft,
                   typename multivec_type::device_type,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    // using lcl_vec_type =
    //   Kokkos::View<double*,
    //                Kokkos::LayoutLeft,
    //                Kokkos::MemoryUnmanaged,
    //                typename multivec_type::device_type>;
    using Tpetra::readOnly;
    using Tpetra::withLocalAccess;

    withLocalAccess
      ([&] (const const_lcl_mv_type& X_lcl) {
         TEST_EQUALITY( static_cast<size_t> (X_lcl.extent (0)),
                        static_cast<size_t> (X.getLocalLength ()) );
       },
       readOnly (X));
    withLocalAccess
      ([&] (const nonconst_lcl_mv_type& X_lcl) {
         TEST_EQUALITY( static_cast<size_t> (X_lcl.extent (0)),
                        static_cast<size_t> (X.getLocalLength ()) );
       },
       writeOnly (X));
    withLocalAccess
      ([&] (const nonconst_lcl_mv_type& X_lcl) {
         TEST_EQUALITY( static_cast<size_t> (X_lcl.extent (0)),
                        static_cast<size_t> (X.getLocalLength ()) );
       },
       readWrite (X));
  }
} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  const int errCode =
    Teuchos::UnitTestRepository::runUnitTestsFromMain (argc, argv);
  return errCode;
}
