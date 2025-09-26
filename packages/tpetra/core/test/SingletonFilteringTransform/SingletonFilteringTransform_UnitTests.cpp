// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_TestingUtilities.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_CrsSingletonFilter_LinearProblem.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>
//#include <unordered_map>


namespace {

  using Tpetra::TestingUtilities::getDefaultComm;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Comm;

  // Utilities
  // --------------------------------------------------------------------------

  template <typename Scalar, typename LO, typename GO, typename Node>
  bool compareCrsMatrices(const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar, LO, GO, Node>>& A,
                          const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar, LO, GO, Node>>& B,
                          const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                          Teuchos::FancyOStream& out)
  {
      bool comparePass = true;
  
      // Compare global dimensions
      if (A->getGlobalNumRows() != B->getGlobalNumRows() ||
          A->getGlobalNumCols() != B->getGlobalNumCols()) {
          out << "Global dimensions do not match." << std::endl;
          out << "Global number of rows: A(" << A->getGlobalNumRows() << " != B(" << B->getGlobalNumRows() << ")" << std::endl;
          out << "Global number of cols: A(" << A->getGlobalNumCols() << " != B(" << B->getGlobalNumCols() << ")" << std::endl;
          comparePass = false;
      }
  
      // Compare row maps
      if (!A->getRowMap()->isSameAs(*B->getRowMap())) {
          out << "Row maps do not match." << std::endl;
          out << "A Matrix:" << std::endl;
          A->getRowMap()->describe(out, Teuchos::VERB_EXTREME);
          out << "B Matrix:" << std::endl;
          B->getRowMap()->describe(out, Teuchos::VERB_EXTREME);
          comparePass = false;
      }
  
      // Compare column maps
      if (!A->getColMap()->isSameAs(*B->getColMap())) {
          out << "Column maps do not match." << std::endl;
          out << "A Matrix:" << std::endl;
          A->getColMap()->describe(out, Teuchos::VERB_EXTREME);
          out << "B Matrix:" << std::endl;
          B->getColMap()->describe(out, Teuchos::VERB_EXTREME);
          comparePass = false;
      }
  
      // Compare entries row by row
      auto rowMap = A->getRowMap();
      for (size_t localRow = 0; localRow < rowMap->getLocalNumElements(); ++localRow) {
          // Get global row index
          auto globalRow = rowMap->getGlobalElement(localRow);
  
          // Buffers for row data
          typename Tpetra::CrsMatrix<Scalar, LO, GO, Node>::nonconst_global_inds_host_view_type indicesA, indicesB;
          typename Tpetra::CrsMatrix<Scalar, LO, GO, Node>::nonconst_values_host_view_type valuesA, valuesB;
  
          size_t numEntriesA = A->getNumEntriesInGlobalRow(globalRow);
          size_t numEntriesB = B->getNumEntriesInGlobalRow(globalRow);
  
          indicesA = typename Tpetra::CrsMatrix<Scalar, LO, GO, Node>::nonconst_global_inds_host_view_type("indicesA", numEntriesA);
          valuesA = typename Tpetra::CrsMatrix<Scalar, LO, GO, Node>::nonconst_values_host_view_type("valuesA", numEntriesA);
  
          indicesB = typename Tpetra::CrsMatrix<Scalar, LO, GO, Node>::nonconst_global_inds_host_view_type("indicesB", numEntriesB);
          valuesB = typename Tpetra::CrsMatrix<Scalar, LO, GO, Node>::nonconst_values_host_view_type("valuesB", numEntriesB);
  
          // Extract row data using getGlobalRowCopy
          A->getGlobalRowCopy(globalRow, indicesA, valuesA, numEntriesA);
          B->getGlobalRowCopy(globalRow, indicesB, valuesB, numEntriesB);
  
          // Compare indices
          if (indicesA.size() != indicesB.size()) {
              out << "Indices in row " << globalRow << " do not match." << std::endl;
              comparePass = false;
              continue;
          }
  
          for (size_t i = 0; i < indicesA.size(); ++i) {
              if (indicesA[i] != indicesB[i]) {
                  out << "Indices in row " << globalRow << " do not match." << std::endl;
                  comparePass = false;
                  break;
              }
          }
  
          // Compare values
          if (valuesA.size() != valuesB.size()) {
              out << "Values in row " << globalRow << " do not match." << std::endl;
              comparePass = false;
              continue;
          }
  
          for (size_t i = 0; i < valuesA.size(); ++i) {
              if (valuesA[i] != valuesB[i]) {
                  out << "Values in row " << globalRow << " do not match." << std::endl;
                  comparePass = false;
                  break;
              }
          }
      }
  
      // Compare communicators
      if (A->getComm() != B->getComm()) {
          out << "Communicators do not match." << std::endl;
          comparePass = false;
      }
  
      return comparePass;
  }


  template <typename Scalar, typename LO, typename GO, typename Node>
  bool compareMultiVectors(const Teuchos::RCP<const Tpetra::MultiVector<Scalar, LO, GO, Node>>& A,
                           const Teuchos::RCP<const Tpetra::MultiVector<Scalar, LO, GO, Node>>& B,
                           const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                           Teuchos::FancyOStream& out,
                           double relativeTolerance = 1e-12)
  {
      bool comparePass = true;
  
      // Compare global dimensions
      if (A->getGlobalLength() != B->getGlobalLength() ||
          A->getNumVectors() != B->getNumVectors()) {
          out << "Global dimensions do not match." << std::endl;
          out << "Global Length   : A(" << A->getGlobalLength() << " != B(" << B->getGlobalLength() << ")" << std::endl;
          out << "Global # vectors: A(" << A->getNumVectors() << " != B(" << B->getNumVectors() << ")" << std::endl;
          return false;
      }
  
      // Compare maps
      if (!A->getMap()->isSameAs(*B->getMap())) {
          out << "Maps do not match." << std::endl;
          out << "A Map:" << std::endl;
          A->getMap()->describe(out, Teuchos::VERB_EXTREME);
          out << "B Map:" << std::endl;
          B->getMap()->describe(out, Teuchos::VERB_EXTREME);
          return false;
      }
  
      // Create a temporary MultiVector to store the difference
      auto diff = Teuchos::rcp(new Tpetra::MultiVector<Scalar, LO, GO, Node>(A->getMap(), A->getNumVectors()));
  
      // Compute the difference: diff = A - B
      diff->update(1.0, *A, -1.0, *B, 0.0);
  
      // Prepare a container to store the norms
      Teuchos::Array<Scalar> norms(A->getNumVectors());
  
      // Compute the 2-norm of the difference for each vector
      diff->norm2(Teuchos::ArrayView<Scalar>(norms));
  
      // Check if the norms are below the relative tolerance
      for (size_t i = 0; i < static_cast<size_t>(norms.size()); ++i) {
          double maxNormA = A->getVector(i)->norm2();
          double maxNormB = B->getVector(i)->norm2();
          double maxNorm = std::max(maxNormA, maxNormB);
  
          if (maxNorm > 0.0) { // Avoid division by zero
              double relativeError = norms[i] / maxNorm;
              if (relativeError > relativeTolerance) {
                  out << "Vector " << i << " failed relative tolerance check." << std::endl;
                  out << "Norm of difference: " << norms[i] << std::endl;
                  out << "Relative error: " << relativeError << std::endl;
  
                  out << "A MultiVector:" << std::endl;
                  A->describe(out, Teuchos::VERB_EXTREME);
                  out << "B MultiVector:" << std::endl;
                  B->describe(out, Teuchos::VERB_EXTREME);
  
                  comparePass = false;
              }
          } else {
              // If both vectors are zero, they are considered equal
              if (norms[i] > relativeTolerance) {
                  out << "Vector " << i << " failed absolute tolerance check." << std::endl;
                  out << "Norm of difference: " << norms[i] << std::endl;
  
                  out << "A MultiVector:" << std::endl;
                  A->describe(out, Teuchos::VERB_EXTREME);
                  out << "B MultiVector:" << std::endl;
                  B->describe(out, Teuchos::VERB_EXTREME);
  
                  comparePass = false;
              }
          }
      }
  
      return comparePass;
  }

  // Unit Tests
  // --------------------------------------------------------------------------

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( SingletonFilteringTransform, SF1, LO, GO, Scalar, Node )

  {
    using CrsMatrix_t     = Tpetra::CrsMatrix    <Scalar, LO, GO, Node>;
    using MultiVector_t   = Tpetra::MultiVector  <Scalar, LO, GO, Node>;
    using Reader_t        = Tpetra::MatrixMarket::Reader<CrsMatrix_t>;
    using Map_t           = Tpetra::Map<LO, GO, Node>;
    using LinearProblem_t = Tpetra::LinearProblem<Scalar, LO, GO, Node>;
    using CrsSingletonFiltering_t = Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>;
    
    auto Comm = Tpetra::getDefaultComm();

    RCP<CrsMatrix_t> A_Original;
    RCP<MultiVector_t> LHS_Original, RHS_Original;

    A_Original   = Reader_t::readSparseFile("SF1_Matrix_Original.mm", Comm);
    RCP<const Map_t> A_Original_Map = A_Original->getRangeMap();
    LHS_Original = Reader_t::readDenseFile("SF1_LHS_Original.mm", Comm, A_Original_Map);
    RHS_Original = Reader_t::readDenseFile("SF1_RHS_Original.mm", Comm, A_Original_Map);

    RCP<MultiVector_t> x = rcp(new MultiVector_t(A_Original_Map, LHS_Original->getNumVectors() ));
    RCP<LinearProblem_t> preSingletonProblem = rcp(new LinearProblem_t( A_Original, x, RHS_Original ));
    CrsSingletonFiltering_t SingletonTrans(true);
    RCP<LinearProblem_t> postSingletonProblem = SingletonTrans( preSingletonProblem );

    SingletonTrans.fwd();

    RCP<CrsMatrix_t> A_Reduced;
    RCP<MultiVector_t> LHS_Reduced, RHS_Reduced;

    auto reducedMap = Reader_t::readMapFile("SF1_Matrix_Reduced_Map.mm", Comm);
    A_Reduced       = Reader_t::readSparseFile("SF1_Matrix_Reduced.mm", reducedMap, reducedMap, reducedMap, reducedMap);

    LHS_Reduced = Reader_t::readDenseFile("SF1_LHS_Reduced.mm", Comm, reducedMap);
    RHS_Reduced = Reader_t::readDenseFile("SF1_RHS_Reduced.mm", Comm, reducedMap);

    TEUCHOS_ASSERT(compareCrsMatrices(
      Teuchos::rcp_dynamic_cast<const CrsMatrix_t>(A_Reduced, true),
      Teuchos::rcp_dynamic_cast<const CrsMatrix_t>(postSingletonProblem->getMatrix(), true),
      Comm, out));

    TEUCHOS_ASSERT(compareMultiVectors(
      Teuchos::rcp_dynamic_cast<const MultiVector_t>(LHS_Reduced, true),
      Teuchos::rcp_dynamic_cast<const MultiVector_t>(postSingletonProblem->getLHS(), true),
      Comm, out));
    TEUCHOS_ASSERT(compareMultiVectors(
      Teuchos::rcp_dynamic_cast<const MultiVector_t>(RHS_Reduced, true),
      Teuchos::rcp_dynamic_cast<const MultiVector_t>(postSingletonProblem->getRHS(), true),
      Comm, out));

    // Uncomment to write out reduced files (e.g., to rebaseline).
    //if constexpr (std::is_same<Scalar, double>::value) {
    //  using Writer_t = Tpetra::MatrixMarket::Writer<CrsMatrix_t>;
    //  RCP<const CrsMatrix_t> reducedMatrix = Teuchos::rcp_dynamic_cast<const CrsMatrix_t>(postSingletonProblem->getMatrix(), true);
    //  Writer_t::writeSparseFile("baseline_SF1_Matrix_Reduced.mm", reducedMatrix);
    //  Writer_t::writeMapFile("baseline_SF1_Matrix_Reduced_Map.mm", *(reducedMatrix->getRowMap()));
    //  Writer_t::writeDenseFile("baseline_SF1_LHS_Reduced.mm", postSingletonProblem->getLHS());
    //  Writer_t::writeDenseFile("baseline_SF1_RHS_Reduced.mm", postSingletonProblem->getRHS());
    //}
  }

#define UNIT_TEST_GROUP(SCALAR, LO, GO, NODE) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SingletonFilteringTransform, SF1, LO, GO, SCALAR, NODE)

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN(UNIT_TEST_GROUP)

}
