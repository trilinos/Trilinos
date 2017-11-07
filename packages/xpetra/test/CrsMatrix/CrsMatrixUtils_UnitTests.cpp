// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu        (jhu@sandia.gov)
//                    Chris Siefert      (csiefer@sandia.gov)
//                    Ray Tuminaro       (rstumin@sandia.gov)
//                    Luc Berger-Vergiat (lberge@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Xpetra_UnitTestHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>

#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
#include <KokkosSparse_CrsMatrix.hpp>
#endif

#include <Xpetra_ConfigDefs.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixUtils.hpp>

namespace {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrixUtils, sortCrsEntriesTpetra, SC, LO, GO, NO )
{
  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  std::vector<size_t> rowptr = {0, 5, 9, 13, 16};
  std::vector<LO>     colind = {5, 7, 0, 1, 2, 6, 0, 1, 3, 0, 8, 2, 3, 1, 2, 3};
  std::vector<SC>     values = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
                                13.0, 14.0, 15.0};

  Teuchos::ArrayView<size_t>  CRS_rowptr(rowptr);
  Teuchos::ArrayView<LO>      CRS_colind(colind);
  Teuchos::ArrayView<SC>      CRS_values(values);

  Xpetra::CrsMatrixUtils<SC,LO,GO,NO>::sortCrsEntries(CRS_rowptr, CRS_colind, CRS_values,
                                                      Xpetra::UseTpetra);
  bool sorted = true;
  for(decltype(CRS_rowptr.size()) row = 1; row < CRS_rowptr.size(); ++row) { // Loop over the matrix rows
    if(CRS_rowptr[row] - CRS_rowptr[row - 1] > 1) { // Check that the row has is more that one entry
      for(size_t col = CRS_rowptr[row - 1] + 1; col < CRS_rowptr[row]; ++col) {
        if(CRS_colind[col - 1] > CRS_colind[col]) {sorted = false;}
      }
    }
  }

  TEST_EQUALITY(sorted, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrixUtils, sortCrsEntriesEpetra, SC, LO, GO, NO )
{
  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  std::vector<size_t> rowptr = {0, 5, 9, 13, 16};
  std::vector<LO>     colind = {5, 7, 0, 1, 2, 6, 0, 1, 3, 0, 8, 2, 3, 1, 2, 3};
  std::vector<SC>     values = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
                                13.0, 14.0, 15.0};

  Teuchos::ArrayView<size_t>  CRS_rowptr(rowptr);
  Teuchos::ArrayView<LO>      CRS_colind(colind);
  Teuchos::ArrayView<SC>      CRS_values(values);

  Xpetra::CrsMatrixUtils<SC,LO,GO,NO>::sortCrsEntries(CRS_rowptr, CRS_colind, CRS_values,
                                                      Xpetra::UseEpetra);
  bool sorted = true;
  for(decltype(CRS_rowptr.size()) row = 1; row < CRS_rowptr.size(); ++row) { // Loop over the matrix rows
    if(CRS_rowptr[row] - CRS_rowptr[row - 1] > 1) { // Check that the row has is more that one entry
      for(size_t col = CRS_rowptr[row - 1] + 1; col < CRS_rowptr[row]; ++col) {
        if(CRS_colind[col - 1] > CRS_colind[col]) {sorted = false;}
      }
    }
  }

  TEST_EQUALITY(sorted, true);
}


TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrixUtils, sortAndMergeCrsEntriesTpetra, SC, LO, GO, NO )
{
  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  std::vector<size_t> rowptr = {0, 7, 12, 17, 21};
  std::vector<LO>     colind = {5, 7, 0, 5, 1, 1, 2,
                                6, 0, 1, 1, 3,
                                8, 0, 8, 2, 3,
                                1, 2, 3, 1};
  std::vector<SC>     values = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  Teuchos::ArrayView<size_t>  CRS_rowptr(rowptr);
  Teuchos::ArrayView<LO>      CRS_colind(colind);
  Teuchos::ArrayView<SC>      CRS_values(values);

  Xpetra::CrsMatrixUtils<SC,LO,GO,NO>::sortAndMergeCrsEntries(CRS_rowptr, CRS_colind, CRS_values,
                                                              Xpetra::UseTpetra);

  std::vector<size_t> resizedRows = {0, 5, 9, 13, 16};
  std::vector<SC> sortedAndMergedValues = {1.0, 2.0, 1.0, 2.0, 1.0,
                                           1.0, 2.0, 1.0, 1.0,
                                           1.0, 1.0, 1.0, 2.0,
                                           2.0, 1.0, 1.0};

  bool resized = true;
  for(size_t ind = 0; ind < resizedRows.size(); ++ind) {
    if(resizedRows[ind] != CRS_rowptr[ind]) {resized = false;}
  }

  bool sorted = true;
  for(decltype(CRS_rowptr.size()) row = 1; row < CRS_rowptr.size(); ++row) { // Loop over the matrix rows
    if(CRS_rowptr[row] - CRS_rowptr[row - 1] > 1) { // Check that the row has is more that one entry
      for(size_t col = CRS_rowptr[row - 1] + 1; col < CRS_rowptr[row]; ++col) {
        if(CRS_colind[col - 1] > CRS_colind[col]) {sorted = false;}
      }
    }
  }

  bool merged = true;
  for(size_t ind = 0; ind < sortedAndMergedValues.size(); ++ind) {
    if(sortedAndMergedValues[ind] != CRS_values[ind]) {merged = false;}
  }

  TEST_EQUALITY(resized, true);
  TEST_EQUALITY(sorted , true);
  TEST_EQUALITY(merged , true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrixUtils, sortAndMergeCrsEntriesEpetra, SC, LO, GO, NO )
{
  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  std::vector<size_t> rowptr = {0, 7, 12, 17, 21};
  std::vector<LO>     colind = {5, 7, 0, 5, 1, 1, 2,
                                6, 0, 1, 1, 3,
                                8, 0, 8, 2, 3,
                                1, 2, 3, 1};
  std::vector<SC>     values = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  Teuchos::ArrayView<size_t>  CRS_rowptr(rowptr);
  Teuchos::ArrayView<LO>      CRS_colind(colind);
  Teuchos::ArrayView<SC>      CRS_values(values);

  Xpetra::CrsMatrixUtils<SC,LO,GO,NO>::sortAndMergeCrsEntries(CRS_rowptr, CRS_colind, CRS_values,
                                                              Xpetra::UseEpetra);

  std::vector<size_t> resizedRows = {0, 5, 9, 13, 16};
  std::vector<SC> sortedAndMergedValues = {1.0, 2.0, 1.0, 2.0, 1.0,
                                           1.0, 2.0, 1.0, 1.0,
                                           1.0, 1.0, 1.0, 2.0,
                                           2.0, 1.0, 1.0};

  bool resized = true;
  for(size_t ind = 0; ind < resizedRows.size(); ++ind) {
    if(resizedRows[ind] != CRS_rowptr[ind]) {resized = false;}
  }

  bool sorted = true;
  for(decltype(CRS_rowptr.size()) row = 1; row < CRS_rowptr.size(); ++row) { // Loop over the matrix rows
    if(CRS_rowptr[row] - CRS_rowptr[row - 1] > 1) { // Check that the row has is more that one entry
      for(size_t col = CRS_rowptr[row - 1] + 1; col < CRS_rowptr[row]; ++col) {
        if(CRS_colind[col - 1] > CRS_colind[col]) {sorted = false;}
      }
    }
  }

  bool merged = true;
  for(size_t ind = 0; ind < sortedAndMergedValues.size(); ++ind) {
    if(sortedAndMergedValues[ind] != CRS_values[ind]) {merged = false;}
  }

  TEST_EQUALITY(resized, true);
  TEST_EQUALITY(sorted , true);
  TEST_EQUALITY(merged , true);
}



  //
  // INSTANTIATIONS
  //

  // List of tests which run only with Tpetra
#define UNIT_TEST_GROUP_ORDINAL_TPETRAONLY(SC,LO,GO,NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrixUtils, sortCrsEntriesTpetra, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrixUtils, sortAndMergeCrsEntriesTpetra, SC, LO, GO, NO)

  // List of tests which run only with Epetra
#define UNIT_TEST_GROUP_ORDINAL_EPETRAONLY(SC,LO,GO,NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrixUtils, sortCrsEntriesEpetra, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrixUtils, sortAndMergeCrsEntriesEpetra, SC, LO, GO, NO)


#if defined(HAVE_XPETRA_TPETRA)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( UNIT_TEST_GROUP_ORDINAL_TPETRAONLY )

#endif


#if defined(HAVE_XPETRA_EPETRA)

#include "Xpetra_Map.hpp" // defines EpetraNode
typedef Xpetra::EpetraNode EpetraNode;
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
UNIT_TEST_GROUP_ORDINAL_EPETRAONLY( double, int, int, EpetraNode )
#endif
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
typedef long long LongLong;
UNIT_TEST_GROUP_ORDINAL_EPETRAONLY( double, int, LongLong, EpetraNode )
#endif

#endif

} // End of namespace
