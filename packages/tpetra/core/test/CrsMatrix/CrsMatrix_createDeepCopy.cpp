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
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_RowGraph.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_createDeepCopy_CrsMatrix.hpp"
#include "Tpetra_Util.hpp"

namespace { // (anonymous)

// Test interfaces that need a RowGraph which is not just a CrsGraph.
template <class LO, class GO, class NT>
class MyRowGraph : public Tpetra::RowGraph<LO, GO, NT> {
  using base_type = Tpetra::RowGraph<LO, GO, NT>;

public:
  using local_ordinal_type = typename base_type::local_ordinal_type;
  using global_ordinal_type = typename base_type::global_ordinal_type;
  using node_type = typename base_type::node_type;

  MyRowGraph (Teuchos::RCP<const base_type> G, const bool permitRowViews) :
    G_ (G),
    supportsRowViews_ (permitRowViews)
  {}

  ~MyRowGraph () override = default;

  Teuchos::RCP<const Teuchos::Comm<int> >
  getComm () const override { return G_->getComm (); }


  Teuchos::RCP<const Tpetra::Map<LO, GO, NT>>
  getRowMap () const override {
    return G_->getRowMap ();
  }

  Teuchos::RCP<const Tpetra::Map<LO, GO, NT>>
  getColMap () const override {
    return G_->getColMap ();
  }

  Teuchos::RCP<const Tpetra::Map<LO, GO, NT>>
  getDomainMap () const override {
    return G_->getDomainMap ();
  }

  Teuchos::RCP<const Tpetra::Map<LO, GO, NT>>
  getRangeMap () const override {
    return G_->getRangeMap ();
  }

  Teuchos::RCP<const Tpetra::Import<LO, GO, NT>>
  getImporter () const override {
    return G_->getImporter ();
  }

  Teuchos::RCP<const Tpetra::Export<LO, GO, NT>>
  getExporter () const override {
    return G_->getExporter ();
  }

  Tpetra::global_size_t getGlobalNumRows () const override {
    return G_->getGlobalNumRows ();
  }

  Tpetra::global_size_t getGlobalNumCols () const override {
    return G_->getGlobalNumCols ();
  }

  size_t getNodeNumRows () const override {
    return G_->getNodeNumRows ();
  }

  size_t getNodeNumCols () const override {
    return G_->getNodeNumCols ();
  }

  GO getIndexBase () const override {
    return G_->getIndexBase ();
  }

  Tpetra::global_size_t getGlobalNumEntries () const override {
    return G_->getGlobalNumEntries ();
  }

  size_t getNodeNumEntries () const override {
    return G_->getNodeNumEntries ();
  }

  size_t getNumEntriesInGlobalRow (GO gblRow) const override {
    return G_->getNumEntriesInGlobalRow (gblRow);
  }

  size_t getNumEntriesInLocalRow (LO lclRow) const override {
    return G_->getNumEntriesInLocalRow (lclRow);
  }

  size_t getGlobalMaxNumRowEntries () const override {
    return G_->getGlobalMaxNumRowEntries ();
  }

  size_t getNodeMaxNumRowEntries () const override {
    return G_->getNodeMaxNumRowEntries ();
  }

  bool hasColMap () const override {
    return G_->hasColMap ();
  }

  bool isLocallyIndexed () const override {
    return G_->isLocallyIndexed ();
  }

  bool isGloballyIndexed () const override {
    return G_->isGloballyIndexed ();
  }

  bool isFillComplete () const override {
    return G_->isFillComplete ();
  }


  void
  getGlobalRowCopy (GO gblRow,
                    const Teuchos::ArrayView<GO>& gblColInds,
                    size_t& numColInds) const override
  {
    G_->getGlobalRowCopy (gblRow, gblColInds, numColInds);
  }

  void
  getLocalRowCopy (LO lclRow,
                   const Teuchos::ArrayView<LO>& lclColInds,
                   size_t& numColInds) const override
  {
    G_->getLocalRowCopy (lclRow, lclColInds, numColInds);
  }

  bool supportsRowViews () const override {
    return supportsRowViews_;
  }

  void
  getLocalRowView (const LO lclRow,
                   Teuchos::ArrayView<const LO>& lclColInds) const override
  {
    G_->getLocalRowView (lclRow, lclColInds);
  }

  void
  getGlobalRowView (const GO gblRow,
                    Teuchos::ArrayView<const GO>& gblColInds) const override
  {
    G_->getGlobalRowView (gblRow, gblColInds);
  }

  void
  pack (const Teuchos::ArrayView<const LO>& exportLIDs,
        Teuchos::Array<GO>& exports,
        const Teuchos::ArrayView<size_t>& numPacketsPerLID,
        size_t& constantNumPackets,
        Tpetra::Distributor& distor) const override
  {
    return G_->pack (exportLIDs, exports, numPacketsPerLID,
                     constantNumPackets, distor);
  }

private:
  Teuchos::RCP<const base_type> G_;
  bool supportsRowViews_ = false;
};

// Test interfaces that need a RowMatrix which is not just a CrsMatrix.
template <class SC, class LO, class GO, class NT>
class MyRowMatrix : public Tpetra::RowMatrix<SC, LO, GO, NT> {
  using base_type = Tpetra::RowMatrix<SC, LO, GO, NT>;

public:
  using scalar_type = typename base_type::scalar_type;
  using local_ordinal_type = typename base_type::local_ordinal_type;
  using global_ordinal_type = typename base_type::global_ordinal_type;
  using node_type = typename base_type::node_type;
  using mag_type = typename base_type::mag_type;

  MyRowMatrix (Teuchos::RCP<base_type> A,
               const bool permitRowViews) :
    A_ (A),
    G_ (new MyRowGraph<LO, GO, NT> (A->getGraph (), permitRowViews)),
    supportsRowViews_ (permitRowViews)
  {}

  ~MyRowMatrix () override = default;

  Teuchos::RCP<const Teuchos::Comm<int> >
  getComm () const override { return A_->getComm (); }


  Teuchos::RCP<const Tpetra::Map<LO, GO, NT>>
  getRowMap () const override {
    return A_->getRowMap ();
  }

  Teuchos::RCP<const Tpetra::Map<LO, GO, NT>>
  getColMap () const override {
    return A_->getColMap ();
  }

  Teuchos::RCP<const Tpetra::RowGraph<LO, GO, NT> >
  getGraph () const override {
    return G_;
  }

  Tpetra::global_size_t getGlobalNumRows () const override {
    return A_->getGlobalNumRows ();
  }

  Tpetra::global_size_t getGlobalNumCols () const override {
    return A_->getGlobalNumCols ();
  }

  size_t getNodeNumRows () const override {
    return A_->getNodeNumRows ();
  }

  size_t getNodeNumCols () const override {
    return A_->getNodeNumCols ();
  }

  GO getIndexBase () const override {
    return A_->getIndexBase ();
  }

  Tpetra::global_size_t getGlobalNumEntries () const override {
    return A_->getGlobalNumEntries ();
  }

  size_t getNodeNumEntries () const override {
    return A_->getNodeNumEntries ();
  }

  size_t getNumEntriesInGlobalRow (GO gblRow) const override {
    return A_->getNumEntriesInGlobalRow (gblRow);
  }

  size_t getNumEntriesInLocalRow (LO lclRow) const override {
    return A_->getNumEntriesInLocalRow (lclRow);
  }

  size_t getGlobalMaxNumRowEntries () const override {
    return A_->getGlobalMaxNumRowEntries ();
  }

  size_t getNodeMaxNumRowEntries () const override {
    return A_->getNodeMaxNumRowEntries ();
  }

  bool hasColMap () const override {
    return A_->hasColMap ();
  }

  bool isLocallyIndexed () const override {
    return A_->isLocallyIndexed ();
  }

  bool isGloballyIndexed () const override {
    return A_->isGloballyIndexed ();
  }

  bool isFillComplete () const override {
    return A_->isFillComplete ();
  }


  bool supportsRowViews () const override {
    return supportsRowViews_;
  }


  void
  getGlobalRowCopy (GO gblRow,
                    const Teuchos::ArrayView<GO>& gblColInds,
                    const Teuchos::ArrayView<SC>& values,
                    size_t& numColInds) const override
  {
    A_->getGlobalRowCopy (gblRow, gblColInds, values, numColInds);
  }

  void
  getLocalRowCopy (LO lclRow,
                   const Teuchos::ArrayView<LO>& lclColInds,
                   const Teuchos::ArrayView<SC>& values,
                   size_t& numColInds) const override
  {
    A_->getLocalRowCopy (lclRow, lclColInds, values, numColInds);
  }

  void
  getGlobalRowView (const GO gblRow,
                    Teuchos::ArrayView<const GO>& gblColInds,
                    Teuchos::ArrayView<const SC>& values) const override
  {
    A_->getGlobalRowView (gblRow, gblColInds, values);
  }

  void
  getLocalRowView (const LO lclRow,
                   Teuchos::ArrayView<const LO>& lclColInds,
                   Teuchos::ArrayView<const SC>& values) const override
  {
    A_->getLocalRowView (lclRow, lclColInds, values);
  }

  void
  getLocalDiagCopy (Tpetra::Vector<SC, LO, GO, NT>& diag) const override
  {
    A_->getLocalDiagCopy (diag);
  }

  void
  leftScale (const Tpetra::Vector<SC, LO, GO, NT>& x) override
  {
    A_->leftScale (x);
  }

  void
  rightScale (const Tpetra::Vector<SC, LO, GO, NT>& x) override
  {
    A_->rightScale (x);
  }

  mag_type getFrobeniusNorm () const override {
    return A_->getFrobeniusNorm ();
  }

  Teuchos::RCP<const Tpetra::Map<LO, GO, NT>>
  getDomainMap () const override {
    return A_->getDomainMap ();
  }

  Teuchos::RCP<const Tpetra::Map<LO, GO, NT>>
  getRangeMap () const override {
    return A_->getRangeMap ();
  }

  // Teuchos::RCP<const Tpetra::Import<LO, GO, NT>>
  // getImporter () const override {
  //   return A_->getImporter ();
  // }

  // Teuchos::RCP<const Tpetra::Export<LO, GO, NT>>
  // getExporter () const override {
  //   return A_->getExporter ();
  // }

  void
  pack (const Teuchos::ArrayView<const LO>& exportLIDs,
        Teuchos::Array<char>& exports,
        const Teuchos::ArrayView<size_t>& numPacketsPerLID,
        size_t& constantNumPackets,
        Tpetra::Distributor& distor) const override
  {
    return A_->pack (exportLIDs, exports, numPacketsPerLID,
                     constantNumPackets, distor);
  }

  void
  apply (const Tpetra::MultiVector<SC, LO, GO, NT>& X,
         Tpetra::MultiVector<SC, LO, GO, NT>& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         SC alpha = Teuchos::ScalarTraits<SC>::one (),
         SC beta = Teuchos::ScalarTraits<SC>::zero ()) const override
  {
    return A_->apply (X, Y, mode, alpha, beta);
  }

private:
  Teuchos::RCP<base_type> A_;
  Teuchos::RCP<MyRowGraph<LO, GO, NT>> G_;
  bool supportsRowViews_ = false;
};

template<class SC, class LO, class GO, class NT>
bool
crsMatrixInstancesEqual (const Tpetra::CrsMatrix<SC, LO, GO, NT>& A,
                         const Tpetra::CrsMatrix<SC, LO, GO, NT>& B)
{
  using Teuchos::outArg;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;

  const Teuchos::Comm<int>& comm = * (A.getMap ()->getComm ());
  int lclSuccess = 1;
  int gblSuccess = 0;

  if (A.getRowMap ().is_null () && ! B.getRowMap ().is_null ()) {
    lclSuccess = 0;
  }
  else if (! A.getRowMap ().is_null () && B.getRowMap ().is_null ()) {
    lclSuccess = 0;
  }
  else if (A.getColMap ().is_null () && ! B.getColMap ().is_null ()) {
    lclSuccess = 0;
  }
  else if (! A.getColMap ().is_null () && B.getColMap ().is_null ()) {
    lclSuccess = 0;
  }
  else if (A.isLocallyIndexed () && ! B.isLocallyIndexed ()) {
    lclSuccess = 0;
  }
  else if (A.isGloballyIndexed () && ! B.isGloballyIndexed ()) {
    lclSuccess = 0;
  }
  else if (! A.isLocallyIndexed () && ! A.isGloballyIndexed () &&
           (B.isLocallyIndexed () || B.isGloballyIndexed ())) {
    lclSuccess = 0;
  }

  reduceAll (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  if (gblSuccess != 1) {
    return false;
  }

  if (! A.getRowMap ().is_null () &&
      ! A.getRowMap ()->isSameAs (* (B.getRowMap ()))) {
    lclSuccess = 0;
  }
  if (! A.getColMap ().is_null () &&
      ! A.getColMap ()->isSameAs (* (B.getColMap ()))) {
    lclSuccess = 0;
  }
  if (! A.getDomainMap ().is_null () &&
      ! B.getDomainMap ().is_null () &&
      ! A.getDomainMap ()->isSameAs (* (B.getDomainMap ()))) {
    lclSuccess = 0;
  }
  if (! A.getRangeMap ().is_null () &&
      ! A.getRangeMap ().is_null () &&
      ! A.getRangeMap ()->isSameAs (* (B.getRangeMap ()))) {
    lclSuccess = 0;
  }

  reduceAll (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  if (gblSuccess != 1) {
    return false;
  }

  const auto& rowMap = * (A.getRowMap ());
  const LO lclNumRows = A.getNodeNumRows ();

  Teuchos::Array<SC> A_valsBuf;
  Teuchos::Array<SC> B_valsBuf;

  if (A.isLocallyIndexed ()) {
    Teuchos::Array<LO> A_lclColIndsBuf;
    Teuchos::Array<LO> B_lclColIndsBuf;

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      size_t A_numEnt = 0;
      size_t B_numEnt = 0;

      const size_t A_numEnt2 = A.getNumEntriesInLocalRow (lclRow);
      if (A_numEnt2 > size_t (A_valsBuf.size ())) {
        A_valsBuf.resize (A_numEnt2);
      }
      if (A_numEnt2 > size_t (A_lclColIndsBuf.size ())) {
        A_lclColIndsBuf.resize (A_numEnt2);
      }
      A.getLocalRowCopy (lclRow, A_lclColIndsBuf (), A_valsBuf (), A_numEnt);

      const size_t B_numEnt2 = B.getNumEntriesInLocalRow (lclRow);
      if (B_numEnt2 > size_t (B_valsBuf.size ())) {
        B_valsBuf.resize (B_numEnt2);
      }
      if (B_numEnt2 > size_t (B_lclColIndsBuf.size ())) {
        B_lclColIndsBuf.resize (B_numEnt2);
      }
      B.getLocalRowCopy (lclRow, B_lclColIndsBuf (), B_valsBuf (), B_numEnt);

      if (A_numEnt != B_numEnt) {
        lclSuccess = 0;
        break;
      }

      Teuchos::ArrayView<LO> A_lclColInds = A_lclColIndsBuf.view (0, A_numEnt);
      Teuchos::ArrayView<SC> A_vals = A_valsBuf.view (0, A_numEnt);
      Tpetra::sort2 (A_lclColInds.begin (), A_lclColInds.end (), A_vals.begin ());

      Teuchos::ArrayView<LO> B_lclColInds = B_lclColIndsBuf.view (0, B_numEnt);
      Teuchos::ArrayView<SC> B_vals = B_valsBuf.view (0, B_numEnt);
      Tpetra::sort2 (B_lclColInds.begin (), B_lclColInds.end (), B_vals.begin ());

      if (! std::equal (A_lclColInds.begin (), A_lclColInds.end (),
                        B_lclColInds.begin ())) {
        lclSuccess = 0;
        break;
      }
      if (! std::equal (A_vals.begin (), A_vals.end (),
                        B_vals.begin ())) {
        lclSuccess = 0;
        break;
      }
    }
  }
  else if (A.isGloballyIndexed ()) {
    Teuchos::Array<GO> A_gblColIndsBuf;
    Teuchos::Array<GO> B_gblColIndsBuf;

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblRow = rowMap.getGlobalElement (lclRow);
      size_t A_numEnt = 0;
      size_t B_numEnt = 0;

      const size_t A_numEnt2 = A.getNumEntriesInGlobalRow (gblRow);
      if (A_numEnt2 > size_t (A_valsBuf.size ())) {
        A_valsBuf.resize (A_numEnt2);
      }
      if (A_numEnt2 > size_t (A_gblColIndsBuf.size ())) {
        A_gblColIndsBuf.resize (A_numEnt2);
      }
      A.getGlobalRowCopy (gblRow, A_gblColIndsBuf (), A_valsBuf (), A_numEnt);

      const size_t B_numEnt2 = B.getNumEntriesInGlobalRow (gblRow);
      if (B_numEnt2 > size_t (B_valsBuf.size ())) {
        B_valsBuf.resize (B_numEnt2);
      }
      if (B_numEnt2 > size_t (B_gblColIndsBuf.size ())) {
        B_gblColIndsBuf.resize (B_numEnt2);
      }
      B.getGlobalRowCopy (gblRow, B_gblColIndsBuf (), B_valsBuf (), B_numEnt);

      if (A_numEnt != B_numEnt) {
        lclSuccess = 0;
        break;
      }

      Teuchos::ArrayView<GO> A_gblColInds = A_gblColIndsBuf.view (0, A_numEnt);
      Teuchos::ArrayView<SC> A_vals = A_valsBuf.view (0, A_numEnt);
      Tpetra::sort2 (A_gblColInds.begin (), A_gblColInds.end (), A_vals.begin ());

      Teuchos::ArrayView<GO> B_gblColInds = B_gblColIndsBuf.view (0, B_numEnt);
      Teuchos::ArrayView<SC> B_vals = B_valsBuf.view (0, B_numEnt);
      Tpetra::sort2 (B_gblColInds.begin (), B_gblColInds.end (), B_vals.begin ());

      if (! std::equal (A_gblColInds.begin (), A_gblColInds.end (),
                        B_gblColInds.begin ())) {
        gblSuccess = 0;
        break;
      }
      if (! std::equal (A_vals.begin (), A_vals.end (),
                        B_vals.begin ())) {
        gblSuccess = 0;
        break;
      }
    }
  }

  reduceAll (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  return gblSuccess == 1;
}

//
// UNIT TESTS
//

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, createDeepCopy, SC, LO, GO, NT )
{
  using Teuchos::RCP;
  using std::endl;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using map_type = Tpetra::Map<LO, GO, NT>;
  constexpr bool testCrsNotFillComplete = true;
  constexpr bool debug = true;

  RCP<Teuchos::FancyOStream> fancyOutPtr = debug ?
    Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
    Teuchos::rcpFromRef (out);
  Teuchos::FancyOStream& myOut = *fancyOutPtr;

  myOut << "Test Tpetra::createDeepCopy(RowMatrix)->CrsMatrix" << endl;
  Teuchos::OSTab tab0 (myOut);

  auto comm = Tpetra::getDefaultComm ();
  const int numProcs = comm->getSize ();
  const LO lclNumRows = 5;
  const GO gblNumRows = GO (numProcs) * GO (lclNumRows);
  const GO indexBase = 0;
  RCP<const map_type> rowMap (new map_type (gblNumRows, lclNumRows,
                                            indexBase, comm));
  {
    myOut << "Test CrsMatrix created with row Map" << endl;
    Teuchos::OSTab tab1 (myOut);
    RCP<crs_matrix_type> A
      (new crs_matrix_type (rowMap, 1));

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblRow = rowMap->getGlobalElement (lclRow);
      const GO gblCol = (gblRow + GO (1)) % GO (gblNumRows);
      const SC val = Teuchos::ScalarTraits<SC>::one ();
      A->insertGlobalValues (gblRow, LO (1), &val, &gblCol);
    }

    if (testCrsNotFillComplete) {
      myOut << "Test createDeepCopy with CrsMatrix not yet "
        "fillComplete" << endl;
      Teuchos::OSTab tab2 (myOut);
      crs_matrix_type A_copy = Tpetra::createDeepCopy (*A);
      TEST_ASSERT( crsMatrixInstancesEqual (*A, A_copy) );
    }
    for (bool hasRowViews : {false, true}) {
      myOut << "Test createDeepCopy with CrsMatrix not yet "
        "fillComplete, wrapped, with hasRowViews="
        << (hasRowViews ? "true" : "false") << endl;
      Teuchos::OSTab tab2 (myOut);
      MyRowMatrix<SC, LO, GO, NT> A_my (A, hasRowViews);
      crs_matrix_type A_copy = Tpetra::createDeepCopy (A_my);
      TEST_ASSERT( crsMatrixInstancesEqual (*A, A_copy) );
    }

    RCP<const map_type> domainMap = rowMap;
    RCP<const map_type> rangeMap = rowMap;
    A->fillComplete (domainMap, rangeMap);

    {
      myOut << "Test createDeepCopy with fillComplete CrsMatrix" << endl;
      Teuchos::OSTab tab2 (myOut);

      crs_matrix_type A_copy = Tpetra::createDeepCopy (*A);
      TEST_ASSERT( crsMatrixInstancesEqual (*A, A_copy) );
    }
    for (bool hasRowViews : {false, true}) {
      myOut << "Test createDeepCopy with fillComplete "
        "CrsMatrix, wrapped, with hasRowViews="
        << (hasRowViews ? "true" : "false") << endl;
      Teuchos::OSTab tab2 (myOut);

      MyRowMatrix<SC, LO, GO, NT> A_my (A, hasRowViews);
      crs_matrix_type A_copy = Tpetra::createDeepCopy (A_my);
      TEST_ASSERT( crsMatrixInstancesEqual (*A, A_copy) );
    }
  }

  {
    myOut << "Test CrsMatrix created with row and column Map" << endl;
    Teuchos::OSTab tab1 (myOut);

    RCP<const map_type> colMap = rowMap;
    RCP<crs_matrix_type> A
      (new crs_matrix_type (rowMap, colMap, 1));

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblRow = rowMap->getGlobalElement (lclRow);
      const GO gblCol = gblRow;
      const SC val = Teuchos::ScalarTraits<SC>::one ();
      A->insertGlobalValues (gblRow, LO (1), &val, &gblCol);
    }

    if (testCrsNotFillComplete) {
      crs_matrix_type A_copy = Tpetra::createDeepCopy (*A);
      TEST_ASSERT( crsMatrixInstancesEqual (*A, A_copy) );
    }
    for (bool hasRowViews : {false, true}) {
      MyRowMatrix<SC, LO, GO, NT> A_my (A, hasRowViews);
      crs_matrix_type A_copy = Tpetra::createDeepCopy (A_my);
      TEST_ASSERT( crsMatrixInstancesEqual (*A, A_copy) );
    }

    RCP<const map_type> domainMap = rowMap;
    RCP<const map_type> rangeMap = rowMap;
    A->fillComplete (domainMap, rangeMap);

    {
      crs_matrix_type A_copy = Tpetra::createDeepCopy (*A);
      TEST_ASSERT( crsMatrixInstancesEqual (*A, A_copy) );
    }
    for (bool hasRowViews : {false, true}) {
      MyRowMatrix<SC, LO, GO, NT> A_my (A, hasRowViews);
      crs_matrix_type A_copy = Tpetra::createDeepCopy (A_my);
      TEST_ASSERT( crsMatrixInstancesEqual (*A, A_copy) );
    }
  }
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SC, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, createDeepCopy, SC, LO, GO, NT )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP )

} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  const int errCode =
    Teuchos::UnitTestRepository::runUnitTestsFromMain (argc, argv);
  return errCode;
}
