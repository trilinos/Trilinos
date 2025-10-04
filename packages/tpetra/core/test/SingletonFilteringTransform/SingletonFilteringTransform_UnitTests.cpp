// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_TestingUtilities.hpp>
#include <Teuchos_UnitTestHelpers.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_CrsSingletonFilter_LinearProblem.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>
//#include <unordered_map>
#include <TpetraExt_MatrixMatrix.hpp>


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
                        const Teuchos::RCP<const Teuchos::Comm<int>>& Comm,
                        Teuchos::FancyOStream& out,
                        double relativeTolerance = 1e-12)
{
    bool comparePass = true;

    // Compare global dimensions
    if (A->getGlobalNumRows() != B->getGlobalNumRows() ||
        A->getGlobalNumCols() != B->getGlobalNumCols()) {
        out << "Global dimensions do not match." << std::endl;
        out << "Global number of rows: A(" << A->getGlobalNumRows() << ") != B(" << B->getGlobalNumRows() << ")" << std::endl;
        out << "Global number of cols: A(" << A->getGlobalNumCols() << ") != B(" << B->getGlobalNumCols() << ")" << std::endl;
        return false;
    }

    // Compare row maps
    if (!A->getRowMap()->isSameAs(*B->getRowMap())) {
        out << "Row maps do not match." << std::endl;
        out << "A Matrix:" << std::endl;
        A->getRowMap()->describe(out, Teuchos::VERB_EXTREME);
        out << "B Matrix:" << std::endl;
        B->getRowMap()->describe(out, Teuchos::VERB_EXTREME);
        return false;
    }

    // Compare column maps
    if (!A->getColMap()->isSameAs(*B->getColMap())) {
        out << "Column maps do not match." << std::endl;
        out << "A Matrix:" << std::endl;
        A->getColMap()->describe(out, Teuchos::VERB_EXTREME);
        out << "B Matrix:" << std::endl;
        B->getColMap()->describe(out, Teuchos::VERB_EXTREME);
        return false;
    }

    if constexpr (!std::is_same<Scalar, long long>::value) {
      // It appears that Reader_t::readSparseFile() reads the file for long long badly when it has real values. 
      // Skip comparing values when Scalar is long long.

      // Create a temporary matrix to store the difference
      const size_t maxNumEntriesPerRow = A->getGlobalMaxNumRowEntries();
      auto diff = Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LO, GO, Node>(A->getRowMap(), A->getColMap(), maxNumEntriesPerRow));

      // Perform the matrix-matrix addition with scalar alpha = 1.0 and beta = -1.0
      Tpetra::MatrixMatrix::Add(
          *A,    // matrix A
          false, // transpose A (don't transpose)
          1.0,  // alpha
          *B,    // matrix B
          false, // transpose B (don't transpose)
          -1.0, // beta
          diff);   // output matrix

      // Compute the norm-2 of the difference
      Scalar diffNorm = diff->getFrobeniusNorm();

      // Check if the norm is below the relative tolerance
      Scalar maxNormA = A->getFrobeniusNorm();
      Scalar maxNormB = B->getFrobeniusNorm();
      Scalar maxNorm = std::max(maxNormA, maxNormB);

      if (maxNorm > 0.0) { // Avoid division by zero
          double relativeError = diffNorm / maxNorm;
          if (relativeError > relativeTolerance) {
              out << "Matrix failed relative tolerance check." << std::endl;
              out << "Norm of difference: " << diffNorm << std::endl;
              out << "Relative error: " << relativeError << std::endl;

              out << "A Matrix:" << std::endl;
              A->describe(out, Teuchos::VERB_EXTREME);
              out << "B Matrix:" << std::endl;
              B->describe(out, Teuchos::VERB_EXTREME);

              comparePass = false;
          }
      } else {
          // If both matrices are zero, they are considered equal
          if (diffNorm > relativeTolerance) {
              out << "Matrix failed absolute tolerance check." << std::endl;
              out << "Norm of difference: " << diffNorm << std::endl;

              out << "A Matrix:" << std::endl;
              A->describe(out, Teuchos::VERB_EXTREME);
              out << "B Matrix:" << std::endl;
              B->describe(out, Teuchos::VERB_EXTREME);

              comparePass = false;
          }
      }
    }

    return comparePass;
}  


  template <typename Scalar, typename LO, typename GO, typename Node>
  bool compareMultiVectors(const Teuchos::RCP<const Tpetra::MultiVector<Scalar, LO, GO, Node>>& A,
                           const Teuchos::RCP<const Tpetra::MultiVector<Scalar, LO, GO, Node>>& B,
                           const Teuchos::RCP<const Teuchos::Comm<int>>& Comm,
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
  
      if constexpr (!std::is_same<Scalar, long long>::value) {
        // It appears that Reader_t::readSparseFile() reads the file for long long badly when it has real values. 
        // Skip comparing values when Scalar is long long.

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
      }
  
      return comparePass;
  }

  template <typename Scalar, typename LO, typename GO, typename Node>
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node>> createReducedMatrix() {
      using map_type = Tpetra::Map<LO, GO, Node>;
      using matrix_type = Tpetra::CrsMatrix<Scalar, LO, GO, Node>;
  
      auto Comm = Tpetra::getDefaultComm();
  
      // Define the row and column indices
      Teuchos::Array<GO> rowIndices = Teuchos::Array<GO>({0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 12, 13, 14, 15, 16, 17, 18});
  
      // Create row and column maps
      Teuchos::RCP<const map_type> rowMap = Teuchos::rcp(new map_type(rowIndices.size(), rowIndices(), 0, Comm));
      Teuchos::RCP<const map_type> colMap = Teuchos::rcp(new map_type(rowIndices.size(), rowIndices(), 0, Comm));
  
      // Create the CrsMatrix
      //Teuchos::RCP<matrix_type> matrix = Teuchos::rcp(new matrix_type(rowMap, colMap, 0));
      const size_t maxNumEntriesPerRow = 6;
      Teuchos::RCP<matrix_type> matrix = Teuchos::rcp(new matrix_type(rowMap, colMap, maxNumEntriesPerRow));
  
      // Insert values into the matrix
      for (size_t i = 0; i < static_cast<size_t>(rowIndices.size()); ++i) {
          GO row = rowIndices[i];
          Teuchos::Array<GO> cols;
          Teuchos::Array<Scalar> values;
  
          // Define the non-zero entries for each row
          switch (row) {
              case 0:
                  cols = Teuchos::Array<GO>({1, 2});
                  values = Teuchos::Array<Scalar>({Scalar(1), Scalar(1)});
                  break;
              case 1:
                  cols = Teuchos::Array<GO>({0, 3});
                  values = Teuchos::Array<Scalar>({Scalar(1), Scalar(-1)});
                  break;
              case 2:
                  cols = Teuchos::Array<GO>({0, 4});
                  values = Teuchos::Array<Scalar>({Scalar(1), Scalar(-1)});
                  break;
              case 3:
                  cols = Teuchos::Array<GO>({1, 5});
                  values = Teuchos::Array<Scalar>({Scalar(-1), Scalar(-1)});
                  break;
              case 4:
                  cols = Teuchos::Array<GO>({2, 6});
                  values = Teuchos::Array<Scalar>({Scalar(-1), Scalar(1)});
                  break;
              case 5:
                  cols = Teuchos::Array<GO>({3, 7});
                  values = Teuchos::Array<Scalar>({Scalar(-1), Scalar(1)});
                  break;
              case 6:
                  cols = Teuchos::Array<GO>({4});
                  values = Teuchos::Array<Scalar>({Scalar(1)});
                  break;
              case 7:
                  cols = Teuchos::Array<GO>({5, 7, 12});
                  values = Teuchos::Array<Scalar>({Scalar(1), Scalar(0.25), Scalar(-0.25)});
                  break;
              case 9:
                  cols = Teuchos::Array<GO>({11});
                  values = Teuchos::Array<Scalar>({Scalar(1)});
                  break;
              case 11:
                  cols = Teuchos::Array<GO>({9, 11, 12, 16});
                  values = Teuchos::Array<Scalar>({Scalar(1), Scalar(9.29367e-05), Scalar(-8.77169e-05), Scalar(-5.21976e-06)});
                  break;
              case 12:
                  cols = Teuchos::Array<GO>({7, 11, 12, 13, 16});
                  values = Teuchos::Array<Scalar>({Scalar(-0.25), Scalar(-5.6556e-35), Scalar(0.25), Scalar(-2.89557e-35), Scalar(-3.29088e-35)});
                  break;
              case 13:
                  cols = Teuchos::Array<GO>({12, 13, 15});
                  values = Teuchos::Array<Scalar>({Scalar(0.0012461), Scalar(0.000160258), Scalar(-1)});
                  break;
              case 14:
                  cols = Teuchos::Array<GO>({14, 17});
                  values = Teuchos::Array<Scalar>({Scalar(0.2), Scalar(-0.2)});
                  break;
              case 15:
                  cols = Teuchos::Array<GO>({13, 17});
                  values = Teuchos::Array<Scalar>({Scalar(-1), Scalar(1)});
                  break;
              case 16:
                  cols = Teuchos::Array<GO>({11, 12, 16, 18});
                  values = Teuchos::Array<Scalar>({Scalar(-9.29367e-05), Scalar(8.77169e-05), Scalar(5.21976e-06), Scalar(-1)});
                  break;
              case 17:
                  cols = Teuchos::Array<GO>({14, 15, 17, 18});
                  values = Teuchos::Array<Scalar>({Scalar(-0.2), Scalar(1), Scalar(0.2), Scalar(1)});
                  break;
              case 18:
                  cols = Teuchos::Array<GO>({16, 17});
                  values = Teuchos::Array<Scalar>({Scalar(-1), Scalar(1)});
                  break;
              default:
                  cols = Teuchos::Array<GO>();
                  values = Teuchos::Array<Scalar>();
                  break;
          }
  
          // Insert the values into the matrix
          if (!cols.empty()) {
              matrix->insertGlobalValues(row, cols(), values());
          }
      }
  
      // Finalize the matrix
      matrix->fillComplete();
  
      return matrix;
  }

template <typename Scalar, typename LO, typename GO, typename Node>
Teuchos::RCP<Tpetra::MultiVector<Scalar, LO, GO, Node>> createMultiVector() {
    using map_type = Tpetra::Map<LO, GO, Node>;
    using multivector_type = Tpetra::MultiVector<Scalar, LO, GO, Node>;

    auto Comm = Tpetra::getDefaultComm();

    // Define the global indices for the rows
    Teuchos::Array<GO> globalIndices = {0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 12, 13, 14, 15, 16, 17, 18};

    // Create the map
    Teuchos::RCP<const map_type> map = Teuchos::rcp(new map_type(globalIndices.size(), globalIndices(), 0, Comm));

    // Create the MultiVector with one vector (column)
    Teuchos::RCP<multivector_type> multiVector = Teuchos::rcp(new multivector_type(map, 1));

    // Define the values for the indices
    Teuchos::Array<Scalar> values = {
        Scalar(0), Scalar(0), Scalar(0), Scalar(0), Scalar(0), Scalar(0), Scalar(0), Scalar(0),
        Scalar(-1), Scalar(-2.31539275374700717e-05), Scalar(0), Scalar(-3.82437638430717104e-04),
        Scalar(0), Scalar(0), Scalar(2.31539275374700751e-05), Scalar(0), Scalar(0)
    };

    // Set the values in the MultiVector
    for (size_t i = 0; i < static_cast<size_t>(globalIndices.size()); ++i) {
        multiVector->replaceGlobalValue(globalIndices[i], 0, values[i]);
    }

    return multiVector;
}

  // Test Class
  // --------------------------------------------------------------------------
  // Derived test class so it is easy to access protected data and functions.
  template <typename Scalar, typename LO, typename GO, typename Node>
  class Test_CrsSingletonFilter_LinearProblem 
    : public Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node> {
    public:

      using row_matrix_type     = Tpetra::RowMatrix<Scalar, LO, GO, Node>;

      void setFullMatrix(Teuchos::RCP<row_matrix_type> matrix) { this->FullMatrix_ = matrix; }

      void test_InitFullMatrixAccess(){
        this->InitFullMatrixAccess();

        TEUCHOS_ASSERT(this->localMaxNumRowEntries_ == 6);
        TEUCHOS_ASSERT(this->FullMatrixIsCrsMatrix_ == true);
        TEUCHOS_ASSERT(this->FullCrsMatrix_ != Teuchos::null);

        for (auto i = 0; i < this->Indices_.size(); ++i) {
            TEUCHOS_ASSERT(this->Indices_[i] == 0);
        }
        for (auto i = 0; i < this->Values_.size(); ++i) {
            TEUCHOS_ASSERT(this->Values_[i] == 0);
        }
      }

      void test_GetRow3(Teuchos::FancyOStream &out, bool &success){ // 3 arguments
        using local_ordinal_type  = Tpetra::Vector<>::local_ordinal_type;
        size_t NumIndices = 1;
        Teuchos::Array<local_ordinal_type> localIndices;
        // Row 0
        this->GetRow(0, NumIndices, localIndices);
        Teuchos::Array<local_ordinal_type> ansIndices = Teuchos::Array<local_ordinal_type>({1, 2});
        TEUCHOS_ASSERT(NumIndices == static_cast<size_t>(ansIndices.size()));
        TEST_COMPARE_ARRAYS(localIndices, ansIndices);
        // Row 8
        this->GetRow(8, NumIndices, localIndices);
        ansIndices = Teuchos::Array<local_ordinal_type>({6, 8, 9, 10, 12, 13});
        TEUCHOS_ASSERT(NumIndices == static_cast<size_t>(ansIndices.size()));
        TEST_COMPARE_ARRAYS(localIndices, ansIndices);
        // Row 10 - singleton
        this->GetRow(10, NumIndices, localIndices);
        ansIndices = Teuchos::Array<local_ordinal_type>({8});
        TEUCHOS_ASSERT(NumIndices == static_cast<size_t>(ansIndices.size()));
        TEST_COMPARE_ARRAYS(localIndices, ansIndices);
      }

      void test_GetRow4(Teuchos::FancyOStream &out, bool &success){ // 4 arguments
        using local_ordinal_type  = Tpetra::Vector<>::local_ordinal_type;
        size_t NumEntries = 0;
        Teuchos::ArrayView<const Scalar> Values;
        Teuchos::ArrayView<const local_ordinal_type> localIndices;
        // Row 8
        this->GetRow(8, NumEntries, Values, localIndices);
        Teuchos::Array<local_ordinal_type> ansIndices = Teuchos::Array<local_ordinal_type>({6, 8, 9, 10, 12, 13});
        TEST_COMPARE_ARRAYS(localIndices, ansIndices);
        Scalar relTol = Scalar(1.0e-05);
        Teuchos::Array<Scalar> ansValues = Teuchos::Array<Scalar>({Scalar(-1.0), Scalar(0.00140635), Scalar(-1.0), Scalar(1.0), Scalar(-0.0012461), Scalar(-0.000160258)});
        if constexpr (std::is_same<Scalar, long long>::value) {
          // It appears that Reader_t::readSparseFile() and Scalar() round/truncate differently for long long,
          // so need to explicitly set ansValues to match Reader_t::readSparseFile().
          ansValues = Teuchos::Array<Scalar>({ Scalar(-1), Scalar(1), Scalar(-1), Scalar(1), Scalar(-1), Scalar(-1) });
        } 
        TEST_COMPARE_FLOATING_ARRAYS(Values, ansValues, relTol);
      }

      void test_GetRowGCIDs(Teuchos::FancyOStream &out, bool &success){
        using global_ordinal_type = Tpetra::Vector<>::global_ordinal_type;
        size_t NumEntries = 0;
        Teuchos::ArrayView<const Scalar> Values;
        Teuchos::Array<global_ordinal_type> Indices;
        // Row 8
        this->GetRowGCIDs(8, NumEntries, Values, Indices);
        Teuchos::Array<global_ordinal_type> ansIndices = Teuchos::Array<global_ordinal_type>({6, 8, 9, 10, 12, 13});
        TEST_COMPARE_ARRAYS(Indices, ansIndices);
        Scalar relTol = Scalar(1.0e-05);
        Teuchos::Array<Scalar> ansValues = Teuchos::Array<Scalar>({Scalar(-1.0), Scalar(0.00140635), Scalar(-1.0), Scalar(1.0), Scalar(-0.0012461), Scalar(-0.000160258)});
        if constexpr (std::is_same<Scalar, long long>::value) {
          // It appears that Reader_t::readSparseFile() and Scalar() round/truncate differently for long long,
          // so need to explicitly set ansValues to match Reader_t::readSparseFile().
          ansValues = Teuchos::Array<Scalar>({ Scalar(-1), Scalar(1), Scalar(-1), Scalar(1), Scalar(-1), Scalar(-1) });
        } 
        TEST_COMPARE_FLOATING_ARRAYS(Values, ansValues, relTol);
      }

      void test_GenerateReducedMap(Teuchos::FancyOStream &out, bool &success){
        using map_type            = Tpetra::Map<LO, GO, Node>;
        using vector_type_int     = Tpetra::Vector<int, LO, GO, Node>;

        Teuchos::RCP<vector_type_int> mapColors = rcp(new vector_type_int(this->FullMatrixRowMap()));
        Teuchos::Array<int> colors = {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0};
        for (size_t i = 0; i < static_cast<size_t>(colors.size()); ++i) {
            mapColors->replaceGlobalValue(static_cast<int>(i), colors[i]);
        }
        Teuchos::RCP<const map_type> reducedMap = this->GenerateReducedMap(this->FullMatrixRowMap(), mapColors, 0);

        auto Comm = Tpetra::getDefaultComm();
        Teuchos::Array<GO> globalIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 12, 13, 14, 15, 16, 17, 18 };
        GO numGlobalElements = 17;
        GO indexBase = 0;
        Teuchos::RCP<map_type> ansReducedMap = Teuchos::rcp(new map_type(numGlobalElements, globalIndices(), indexBase, Comm));
        
        TEST_ASSERT(reducedMap->isSameAs(*ansReducedMap));
      }

      void test_Analyze(Teuchos::FancyOStream &out, bool &success){
        using local_ordinal_type  = Tpetra::Vector<>::local_ordinal_type;
        using vector_type_int     = Tpetra::Vector<int, LO, GO, Node>;

        TEST_ASSERT(this->origObj_ != Teuchos::null);
        TEST_ASSERT(this->FullProblem() == Teuchos::null);
        TEST_ASSERT(this->AnalysisDone_ == true);
        TEST_ASSERT(this->FullMatrix()->getGlobalNumRows() == 19);
        TEST_ASSERT(this->FullMatrix()->getGlobalNumEntries() == 54);

        // Check color maps
        Teuchos::RCP<vector_type_int> mapColors = rcp(new vector_type_int(this->FullMatrixRowMap()));
        Teuchos::Array<int> colors = {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0};
        for (size_t i = 0; i < static_cast<size_t>(colors.size()); ++i) {
            mapColors->replaceGlobalValue(static_cast<int>(i), colors[i]);
        }
        auto rowMapColorsData = this->RowMapColors_->getLocalViewHost(Tpetra::Access::ReadOnly);
        auto mapColorsData = mapColors->getLocalViewHost(Tpetra::Access::ReadOnly);
        TEST_ASSERT(rowMapColorsData.extent(0) == mapColorsData.extent(0));
        TEST_ASSERT(rowMapColorsData.extent(1) == mapColorsData.extent(1));
        TEST_ASSERT(rowMapColorsData.extent(1) == 1);
        for (size_t i = 0; i < rowMapColorsData.extent(0); ++i) {
          TEST_ASSERT(rowMapColorsData(i, 0) == mapColorsData(i, 0));
        }
        auto colMapColorsData = this->ColMapColors_->getLocalViewHost(Tpetra::Access::ReadOnly);
        TEST_ASSERT(colMapColorsData.extent(0) == mapColorsData.extent(0));
        TEST_ASSERT(colMapColorsData.extent(1) == mapColorsData.extent(1));
        TEST_ASSERT(colMapColorsData.extent(1) == 1);
        for (size_t i = 0; i < colMapColorsData.extent(0); ++i) {
          TEST_ASSERT(colMapColorsData(i, 0) == mapColorsData(i, 0));
        }

        TEST_ASSERT(this->localNumSingletonRows_ == 1);
        TEST_ASSERT(this->localNumSingletonCols_ == 1);

        Teuchos::ArrayRCP<local_ordinal_type> ansColSingletonRowLIDs 
          = Teuchos::arcp(new local_ordinal_type[1]{8}, 0, 1, true);
        TEST_COMPARE_ARRAYS(this->ColSingletonRowLIDs_, ansColSingletonRowLIDs);

        Teuchos::ArrayRCP<local_ordinal_type> ansColSingletonColLIDs 
          = Teuchos::arcp(new local_ordinal_type[1]{10}, 0, 1, true);
        TEST_COMPARE_ARRAYS(this->ColSingletonColLIDs_, ansColSingletonColLIDs);
      }

      void test_Operator(Teuchos::FancyOStream &out, bool &success){
        //using local_ordinal_type  = Tpetra::Vector<>::local_ordinal_type;
        //using vector_type_int     = Tpetra::Vector<int, LO, GO, Node>;
        using map_type            = Tpetra::Map<LO, GO, Node>;

        // Tests for ConstructReducedProblem()
        TEST_ASSERT(this->HaveReducedProblem_ == true);
        TEST_ASSERT(this->FullProblem() != Teuchos::null);
        TEST_ASSERT(this->FullMatrix() != Teuchos::null);
        TEST_ASSERT(this->FullProblem()->getRHS() != Teuchos::null);
        TEST_ASSERT(this->FullProblem()->getLHS() != Teuchos::null);
        TEST_ASSERT(this->SingletonsDetected() == true);

        auto Comm = Tpetra::getDefaultComm();
        Teuchos::Array<GO> globalIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 12, 13, 14, 15, 16, 17, 18 };
        GO numGlobalElements = 17;
        GO indexBase = 0;
        Teuchos::RCP<map_type> ansReducedMap = Teuchos::rcp(new map_type(numGlobalElements, globalIndices(), indexBase, Comm));
        
        TEST_ASSERT(this->ReducedMatrixRowMap_->isSameAs(*ansReducedMap));
        TEST_ASSERT(this->ReducedMatrixColMap_->isSameAs(*ansReducedMap));

        Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node>> ansReducedMatrix = createReducedMatrix<Scalar, LO, GO, Node>();
        auto reducedMatrix = Teuchos::rcp_dynamic_cast<const Tpetra::CrsMatrix<Scalar, LO, GO, Node>>(this->ReducedMatrix(), true);
        bool matricesAreEqual = compareCrsMatrices<Scalar, LO, GO, Node>(reducedMatrix, ansReducedMatrix, Comm, out, Scalar(1.0e-05));
        TEUCHOS_ASSERT(matricesAreEqual);

        auto ansReducedRHS_ = createMultiVector<Scalar, LO, GO, Node>();
        auto reducedRHS = Teuchos::rcp_dynamic_cast<const Tpetra::MultiVector<Scalar, LO, GO, Node>>(this->ReducedRHS_, true);
        bool multiVectorsAreEqual = compareMultiVectors<Scalar, LO, GO, Node>(reducedRHS, ansReducedRHS_, Comm, out, 1.0e-5);
        TEUCHOS_ASSERT(multiVectorsAreEqual);

        // Tests for ReducedProblem()
        TEST_ASSERT(this->newObj_ != Teuchos::null);

        TEST_ASSERT(this->FullMatrix()->getGlobalNumRows() == 19);
        TEST_ASSERT(this->FullMatrix()->getGlobalNumEntries() == 54);

        TEST_ASSERT(this->ReducedMatrix()->getGlobalNumRows() == 17);
        TEST_ASSERT(this->ReducedMatrix()->getGlobalNumEntries() == 43);

      }
  };


  // Unit Tests
  // --------------------------------------------------------------------------


  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( SingletonFilteringTransform, Functions, LO, GO, Scalar, Node )
  {
    using CrsMatrix_t     = Tpetra::CrsMatrix    <Scalar, LO, GO, Node>;
    using MultiVector_t   = Tpetra::MultiVector  <Scalar, LO, GO, Node>;
    using Reader_t        = Tpetra::MatrixMarket::Reader<CrsMatrix_t>;
    
    auto Comm = Tpetra::getDefaultComm();

    RCP<CrsMatrix_t> A_Original;
    RCP<MultiVector_t> LHS_Original, RHS_Original;

    Test_CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node> test_SF;
    A_Original   = Reader_t::readSparseFile("SF1_Matrix_Original.mm", Comm);
    test_SF.setFullMatrix(A_Original);

    test_SF.test_InitFullMatrixAccess();
    test_SF.test_GetRow3(out, success);
    test_SF.test_GetRow4(out, success);
    test_SF.test_GetRowGCIDs(out, success);
    test_SF.test_GenerateReducedMap(out, success);
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( SingletonFilteringTransform, Analyze, LO, GO, Scalar, Node )
  {
    using CrsMatrix_t     = Tpetra::CrsMatrix    <Scalar, LO, GO, Node>;
    using MultiVector_t   = Tpetra::MultiVector  <Scalar, LO, GO, Node>;
    using Reader_t        = Tpetra::MatrixMarket::Reader<CrsMatrix_t>;
    using Map_t           = Tpetra::Map<LO, GO, Node>;
    using LinearProblem_t = Tpetra::LinearProblem<Scalar, LO, GO, Node>;
    
    auto Comm = Tpetra::getDefaultComm();

    RCP<CrsMatrix_t> A_Original;
    RCP<MultiVector_t> LHS_Original, RHS_Original;

    A_Original   = Reader_t::readSparseFile("SF1_Matrix_Original.mm", Comm);
    RCP<const Map_t> A_Original_Map = A_Original->getRangeMap();
    LHS_Original = Reader_t::readDenseFile("SF1_LHS_Original.mm", Comm, A_Original_Map);
    RHS_Original = Reader_t::readDenseFile("SF1_RHS_Original.mm", Comm, A_Original_Map);

    RCP<MultiVector_t> x = rcp(new MultiVector_t(A_Original_Map, LHS_Original->getNumVectors() ));
    RCP<LinearProblem_t> preSingletonProblem = rcp(new LinearProblem_t( A_Original, x, RHS_Original ));

    Test_CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node> test_SF;
    test_SF.analyze(preSingletonProblem);
    test_SF.test_Analyze(out, success);

  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( SingletonFilteringTransform, Operator, LO, GO, Scalar, Node )
  {
    // operator() just calls analyze(LinearProblem) and construct(), so this test will
    // basically cover construct() since analyze() is covered above.
    using CrsMatrix_t     = Tpetra::CrsMatrix    <Scalar, LO, GO, Node>;
    using MultiVector_t   = Tpetra::MultiVector  <Scalar, LO, GO, Node>;
    using Reader_t        = Tpetra::MatrixMarket::Reader<CrsMatrix_t>;
    using Map_t           = Tpetra::Map<LO, GO, Node>;
    using LinearProblem_t = Tpetra::LinearProblem<Scalar, LO, GO, Node>;
    
    auto Comm = Tpetra::getDefaultComm();

    RCP<CrsMatrix_t> A_Original;
    RCP<MultiVector_t> LHS_Original, RHS_Original;

    A_Original   = Reader_t::readSparseFile("SF1_Matrix_Original.mm", Comm);
    RCP<const Map_t> A_Original_Map = A_Original->getRangeMap();
    LHS_Original = Reader_t::readDenseFile("SF1_LHS_Original.mm", Comm, A_Original_Map);
    RHS_Original = Reader_t::readDenseFile("SF1_RHS_Original.mm", Comm, A_Original_Map);

    RCP<MultiVector_t> x = rcp(new MultiVector_t(A_Original_Map, LHS_Original->getNumVectors() ));
    RCP<LinearProblem_t> preSingletonProblem = rcp(new LinearProblem_t( A_Original, x, RHS_Original ));

    Test_CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node> test_SF;
    RCP<LinearProblem_t> postSingletonProblem = test_SF( preSingletonProblem );
    test_SF.test_Operator(out, success);

  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( SingletonFilteringTransform, fwd, LO, GO, Scalar, Node )
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
    CrsSingletonFiltering_t SingletonTransform(true);
    RCP<LinearProblem_t> postSingletonProblem = SingletonTransform( preSingletonProblem );

    SingletonTransform.fwd();

    RCP<CrsMatrix_t> A_Reduced;
    RCP<MultiVector_t> LHS_Reduced, RHS_Reduced;

    auto reducedMap = Reader_t::readMapFile("SF1_Matrix_Reduced_Map.mm", Comm);
    A_Reduced       = Reader_t::readSparseFile("SF1_Matrix_Reduced.mm", reducedMap, reducedMap, reducedMap, reducedMap);

    LHS_Reduced = Reader_t::readDenseFile("SF1_LHS_Reduced.mm", Comm, reducedMap);
    RHS_Reduced = Reader_t::readDenseFile("SF1_RHS_Reduced.mm", Comm, reducedMap);

    TEUCHOS_ASSERT(compareCrsMatrices(
      Teuchos::rcp_dynamic_cast<const CrsMatrix_t>(A_Reduced, true),
      Teuchos::rcp_dynamic_cast<const CrsMatrix_t>(postSingletonProblem->getMatrix(), true),
      Comm, out, Scalar(1.0e-05)));

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
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SingletonFilteringTransform, Functions, LO, GO, SCALAR, NODE) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SingletonFilteringTransform, Analyze, LO, GO, SCALAR, NODE) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SingletonFilteringTransform, Operator, LO, GO, SCALAR, NODE) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SingletonFilteringTransform, fwd, LO, GO, SCALAR, NODE) 

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN(UNIT_TEST_GROUP)

}
