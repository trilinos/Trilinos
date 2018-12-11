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

#include "Teuchos_UnitTestHarness.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Export.hpp"
#include "Teuchos_CommHelpers.hpp"

#include <iostream>
#include <numeric>
#include <vector>

// Test for Trilinos GitHub Issue #3968.

namespace { // (anonymous)

  using GST = Tpetra::global_size_t;
  using LO = Tpetra::Map<>::local_ordinal_type;

  // Try to get a 64-bit GlobalOrdinal type; that should be more
  // likely to manifest atomic update issues.  If we can't, though,
  // then a 32-bit type is fine.  long long is guaranteed to be at
  // least 64 bits, so prefer that over long or unsigned long.
#if defined(HAVE_TPETRA_INST_INT_LONG_LONG)
  using GO = long long; // always at least 64 bits
#elif defined(HAVE_TPETRA_INST_INT_LONG)
  using GO = long; // may be 64 or 32 bits
#elif defined(HAVE_TPETRA_INST_INT_UNSIGNED_LONG)
  using GO = unsigned long; // may be 64 or 32 bits
#else
  using GO = Tpetra::Map<>::global_ordinal_type;
#endif

  using map_type = Tpetra::Map<LO, GO>;
  using vec_type = Tpetra::Vector<GO, LO, GO>;
  using export_type = Tpetra::Export<LO, GO>;

  Teuchos::RCP<const map_type>
  createTargetMap (const Teuchos::RCP<const Teuchos::Comm<int>>& comm)
  {
    const GST gblNumInds = static_cast<GST> (comm->getSize ());
    const GO indexBase = 0;
    return Teuchos::rcp (new map_type (gblNumInds, indexBase, comm));
  }

  // As overlapping as possible, but avoid the explicitly "locally
  // replicated" path, just in case Tpetra optimizes for that (it
  // doesn't really, but it could).
  Teuchos::RCP<const map_type>
  createSourceMap (const Teuchos::RCP<const Teuchos::Comm<int>>& comm)
  {
    const GST flag = Teuchos::OrdinalTraits<GST>::invalid ();
    const LO lclNumInds = static_cast<LO> (comm->getSize ());
    std::vector<GO> gblInds (lclNumInds);
    std::iota (std::begin (gblInds), std::end (gblInds), GO (0));
    const GO indexBase = 0;
    return Teuchos::rcp (new map_type (flag, gblInds.data (), lclNumInds, indexBase, comm));
  }

  // "Base value" should be some value, not default, not negative (in
  // case GO is unsigned, which is not default but possible).  Pick
  // baseValue() > comm->getSize().
  constexpr GO baseValue () {
    return static_cast<GO> (100);
  }

  void
  fillTargetVector (vec_type& tgtVector)
  {
    tgtVector.putScalar (baseValue ());
  }

  void
  fillSourceVector (vec_type& srcVector)
  {
    srcVector.sync_host ();
    auto X_lcl_2d_h = srcVector.getLocalViewHost ();
    auto X_lcl_1d_h = Kokkos::subview (X_lcl_2d_h, Kokkos::ALL (), 0);
    const LO lclNumRows = static_cast<LO> (srcVector.getLocalLength ());

    // Combining actually goes across process ranks here, not local
    // rows.  This means that the source vector value we want depends
    // only on the rank of the calling process.  We can reuse that
    // value for all the rows.
    auto comm = srcVector.getMap ()->getComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    // "Interleave" values, to make correct implementation of atomic
    // absmax more relevant.  See GitHub Issue #3968 discussion.
    const GO incrValue = (myRank % 2 == 0) ?
      GO (numProcs - (myRank/2) - 1) :
      GO (myRank - 1) / GO (2);

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      X_lcl_1d_h(lclRow) = baseValue () + incrValue;
    }
    srcVector.sync_device ();
  }

  void
  testVectorAbsMaxExport (bool& success,
                          Teuchos::FancyOStream& out,
                          vec_type& tgtVector,
                          vec_type& srcVector,
                          const export_type& exporter,
                          const int numRepeats)
  {
    // Do some warm-up runs first that don't touch host.
    for (int repeat = 0; repeat < numRepeats; ++repeat) {
      tgtVector.putScalar (baseValue ());
      srcVector.putScalar (baseValue ());
      tgtVector.doExport (srcVector, exporter, Tpetra::ABSMAX);
    }
    // This run has the values about which we actually care.
    fillTargetVector (tgtVector);
    fillSourceVector (srcVector);
    tgtVector.doExport (srcVector, exporter, Tpetra::ABSMAX);

    tgtVector.sync_host ();
    auto X_lcl_2d_h = tgtVector.getLocalViewHost ();
    auto X_lcl_1d_h = Kokkos::subview (X_lcl_2d_h, Kokkos::ALL (), 0);

    auto comm = tgtVector.getMap ()->getComm ();
    const GO incrValue = GO (comm->getSize () - 1);
    const GO expectedValue = baseValue () + incrValue;
    TEST_EQUALITY( X_lcl_1d_h(0), expectedValue );

    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0; // output argument
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;

    reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }

  TEUCHOS_UNIT_TEST( VectorExport, Issue3968 )
  {
    auto comm = Tpetra::getDefaultComm ();
    auto srcMap = createSourceMap (comm);
    auto tgtMap = createTargetMap (comm);

    vec_type srcVector (srcMap);
    vec_type tgtVector (tgtMap);

    export_type exporter (srcMap, tgtMap);
    constexpr int numRepeats = 20;
    testVectorAbsMaxExport (success, out, tgtVector, srcVector,
                            exporter, numRepeats);
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
