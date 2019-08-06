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

#ifdef TPETRA_ENABLE_DEPRECATED_CODE
  TPETRA_DEPRECATED Teuchos::RCP<NT> getNode () const override {
    return G_->getNode ();
  }
#endif // TPETRA_ENABLE_DEPRECATED_CODE

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

#ifdef TPETRA_ENABLE_DEPRECATED_CODE
  global_size_t TPETRA_DEPRECATED getGlobalNumDiags () const override {
    return G_->getGlobalNumDiags ();
  }

  size_t TPETRA_DEPRECATED getNodeNumDiags () const override {
    return G_->getNodeNumDiags ();
  }

  bool TPETRA_DEPRECATED isLowerTriangular () const override {
    return G_->isLowerTriangular ();
  }

  bool TPETRA_DEPRECATED isUpperTriangular () const override {
    return G_->isUpperTriangular ();
  }
#endif // TPETRA_ENABLE_DEPRECATED_CODE

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

#ifdef TPETRA_ENABLE_DEPRECATED_CODE
  TPETRA_DEPRECATED Teuchos::RCP<NT> getNode () const override {
    return A_->getNode ();
  }
#endif // TPETRA_ENABLE_DEPRECATED_CODE

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

#ifdef TPETRA_ENABLE_DEPRECATED_CODE
  global_size_t TPETRA_DEPRECATED getGlobalNumDiags () const override {
    return A_->getGlobalNumDiags ();
  }

  size_t TPETRA_DEPRECATED getNodeNumDiags () const override {
    return A_->getNodeNumDiags ();
  }

  bool TPETRA_DEPRECATED isLowerTriangular () const override {
    return A_->isLowerTriangular ();
  }

  bool TPETRA_DEPRECATED isUpperTriangular () const override {
    return A_->isUpperTriangular ();
  }
#endif // TPETRA_ENABLE_DEPRECATED_CODE

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
#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    RCP<crs_matrix_type> A
      (new crs_matrix_type (rowMap, 1, Tpetra::StaticProfile));
#else
    RCP<crs_matrix_type> A
      (new crs_matrix_type (rowMap, 1));
#endif // TPETRA_ENABLE_DEPRECATED_CODE

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
    }
    for (bool hasRowViews : {false, true}) {
      myOut << "Test createDeepCopy with CrsMatrix not yet "
        "fillComplete, wrapped, with hasRowViews="
        << (hasRowViews ? "true" : "false") << endl;
      Teuchos::OSTab tab2 (myOut);
      MyRowMatrix<SC, LO, GO, NT> A_my (A, hasRowViews);
      crs_matrix_type A_copy = Tpetra::createDeepCopy (A_my);
    }

    RCP<const map_type> domainMap = rowMap;
    RCP<const map_type> rangeMap = rowMap;
    A->fillComplete (domainMap, rangeMap);

    {
      myOut << "Test createDeepCopy with fillComplete CrsMatrix" << endl;
      Teuchos::OSTab tab2 (myOut);

      crs_matrix_type A_copy = Tpetra::createDeepCopy (*A);
    }
    for (bool hasRowViews : {false, true}) {
      myOut << "Test createDeepCopy with fillComplete "
        "CrsMatrix, wrapped, with hasRowViews="
        << (hasRowViews ? "true" : "false") << endl;
      Teuchos::OSTab tab2 (myOut);

      MyRowMatrix<SC, LO, GO, NT> A_my (A, hasRowViews);
      crs_matrix_type A_copy = Tpetra::createDeepCopy (A_my);
    }
  }

  {
    myOut << "Test CrsMatrix created with row and column Map" << endl;
    Teuchos::OSTab tab1 (myOut);

    RCP<const map_type> colMap = rowMap;
#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    RCP<crs_matrix_type> A
      (new crs_matrix_type (rowMap, colMap, 1,
                            Tpetra::StaticProfile));
#else
    RCP<crs_matrix_type> A
      (new crs_matrix_type (rowMap, colMap, 1));
#endif // TPETRA_ENABLE_DEPRECATED_CODE

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblRow = rowMap->getGlobalElement (lclRow);
      const GO gblCol = gblRow;
      const SC val = Teuchos::ScalarTraits<SC>::one ();
      A->insertGlobalValues (gblRow, LO (1), &val, &gblCol);
    }

    if (testCrsNotFillComplete) {
      crs_matrix_type A_copy = Tpetra::createDeepCopy (*A);
    }
    for (bool hasRowViews : {false, true}) {
      MyRowMatrix<SC, LO, GO, NT> A_my (A, hasRowViews);
      crs_matrix_type A_copy = Tpetra::createDeepCopy (A_my);
    }

    RCP<const map_type> domainMap = rowMap;
    RCP<const map_type> rangeMap = rowMap;
    A->fillComplete (domainMap, rangeMap);

    {
      crs_matrix_type A_copy = Tpetra::createDeepCopy (*A);
    }
    for (bool hasRowViews : {false, true}) {
      MyRowMatrix<SC, LO, GO, NT> A_my (A, hasRowViews);
      crs_matrix_type A_copy = Tpetra::createDeepCopy (A_my);
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
