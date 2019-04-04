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

