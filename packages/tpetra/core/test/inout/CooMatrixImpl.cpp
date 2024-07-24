// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Details_CooMatrix.hpp"

namespace { // (anonymous)

using Tpetra::Details::Impl::CooMatrixImpl;
//using Tpetra::TestingUtilities::getDefaultComm;
using std::endl;
using std::size_t;

typedef double SC;
typedef int LO;
typedef long long GO;

// // unused:
// constexpr int expectedNumEnt = 5;
// const GO expectedGblRowInds[] = {93, 418, 31, 666, 93};
// const GO expectedGblColInds[] = {7, 4, 6, 5, 3};
// const SC expectedVals[] = {7.0, 4.0, 6.0, 5.0, 3.0};

void
testCooMatrixImpl (bool& success, Teuchos::FancyOStream& out)
{
  CooMatrixImpl<SC, GO> A;
  TEST_ASSERT( A.getLclNumEntries () == static_cast<size_t> (0) );
  {
    std::vector<GO> gblRowInds;
    A.getMyGlobalRowIndices (gblRowInds);
    TEST_ASSERT( gblRowInds.size () == static_cast<size_t> (0) );
  }
  {
    std::vector<size_t> rowOffsets;
    std::vector<GO> gblColInds;
    std::vector<SC> vals;

    A.buildCrs (rowOffsets, gblColInds.data (), vals.data ());
    TEST_ASSERT( rowOffsets.size () == static_cast<size_t> (1) );
    if (rowOffsets.size () != 0) {
      TEST_ASSERT( rowOffsets[0] == static_cast<size_t> (0) );
    }
  }
  {
    std::vector<size_t> rowOffsets;
    std::vector<LO> lclColInds;
    std::vector<SC> vals;

    A.buildLocallyIndexedCrs<size_t, LO> (rowOffsets,
                                               lclColInds.data (),
                                               vals.data (),
                                               [] (const GO gblColInd) {
                                                 return static_cast<LO> (gblColInd);
                                               });
    TEST_ASSERT( rowOffsets.size () == static_cast<size_t> (1) );
    if (rowOffsets.size () != 0) {
      TEST_ASSERT( rowOffsets[0] == static_cast<size_t> (0) );
    }
  }

  const size_t origNumEnt = 6;
  const GO origGblRowInds[] = {418, 31, 31, 31, 666, 93};
  const GO origGblColInds[] = {7, 5, 6, 5, 4, 3};
  const SC origVals[] = {7.0, 1.0, 2.0, 3.0, 10.0, 11.0};

  A.sumIntoGlobalValues (origGblRowInds, origGblColInds, origVals, origNumEnt);

  const size_t expectedNumEnt = 5;
  TEST_ASSERT( A.getLclNumEntries () == expectedNumEnt );
  if (A.getLclNumEntries () != expectedNumEnt) {
    return;
  }
  {
    const std::vector<GO> expectedGblRowInds {31, 93, 418, 666};

    std::vector<GO> gblRowInds;
    A.getMyGlobalRowIndices (gblRowInds);
    TEST_ASSERT( gblRowInds.size () == expectedGblRowInds.size () );
    if (gblRowInds.size () == expectedGblRowInds.size ()) {
      TEST_ASSERT( gblRowInds == expectedGblRowInds );
    }
  }
  {
    const std::vector<size_t> expectedRowOffsets {0, 2, 3, 4, 5};
    const std::vector<GO> expectedGblColInds {5, 6, 3, 7, 4};
    const std::vector<SC> expectedVals {4.0, 2.0, 11.0, 7.0, 10.0};

    std::vector<size_t> rowOffsets;
    std::vector<GO> gblColInds (A.getLclNumEntries ());
    std::vector<SC> vals (A.getLclNumEntries ());
    A.buildCrs (rowOffsets, gblColInds.data (), vals.data ());

    TEST_EQUALITY( rowOffsets.size (), expectedRowOffsets.size () );
    if (rowOffsets.size () == expectedRowOffsets.size ()) {
      TEST_ASSERT( rowOffsets == expectedRowOffsets );
    }
    else {
      out << "rowOffsets: [";
      for (size_t k = 0; k < rowOffsets.size (); ++k) {
        out << rowOffsets[k];
        if (k + 1 < rowOffsets.size ()) {
          out << ", ";
        }
      }
      out << "]" << endl;
    }
    TEST_EQUALITY( gblColInds.size (), expectedGblColInds.size () );
    if (gblColInds.size () == expectedGblColInds.size ()) {
      TEST_ASSERT( gblColInds == expectedGblColInds );
    }
    TEST_EQUALITY( vals.size (), expectedVals.size () );
    if (vals.size () == expectedVals.size ()) {
      TEST_ASSERT( vals == expectedVals );
    }
  }
  {
    const std::vector<size_t> expectedRowOffsets {0, 2, 3, 4, 5};
    const std::vector<LO> expectedLclColInds {5, 6, 3, 7, 4};
    const std::vector<SC> expectedVals {4.0, 2.0, 11.0, 7.0, 10.0};

    std::vector<size_t> rowOffsets;
    std::vector<LO> lclColInds (A.getLclNumEntries ());
    std::vector<SC> vals (A.getLclNumEntries ());
    A.buildLocallyIndexedCrs<size_t, LO> (rowOffsets,
                                          lclColInds.data (),
                                          vals.data (),
                                          [] (const GO gblInd) {
                                            return static_cast<LO> (gblInd); });

    TEST_EQUALITY( rowOffsets.size (), expectedRowOffsets.size () );
    if (rowOffsets.size () == expectedRowOffsets.size ()) {
      TEST_ASSERT( rowOffsets == expectedRowOffsets );
    }
    else {
      out << "rowOffsets: [";
      for (size_t k = 0; k < rowOffsets.size (); ++k) {
        out << rowOffsets[k];
        if (k + 1 < rowOffsets.size ()) {
          out << ", ";
        }
      }
      out << "]" << endl;
    }
    TEST_EQUALITY( lclColInds.size (), expectedLclColInds.size () );
    if (lclColInds.size () == expectedLclColInds.size ()) {
      TEST_ASSERT( lclColInds == expectedLclColInds );
    }
    TEST_EQUALITY( vals.size (), expectedVals.size () );
    if (vals.size () == expectedVals.size ()) {
      TEST_ASSERT( vals == expectedVals );
    }
  }
}

TEUCHOS_UNIT_TEST( CooMatrixImpl, doubleIntLongLong )
{
  out << "Testing Tpetra::Details::Impl::CooMatrixImpl" << endl;
  Teuchos::OSTab tab1 (out);

  //testThing (success, out, * (getDefaultComm ()));
  testCooMatrixImpl (success, out);
}

} // namespace (anonymous)

