// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
    auto X_lcl_2d_h = srcVector.getLocalViewHost(Tpetra::Access::OverwriteAll);
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

    auto X_lcl_2d_h = tgtVector.getLocalViewHost(Tpetra::Access::ReadOnly);
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
