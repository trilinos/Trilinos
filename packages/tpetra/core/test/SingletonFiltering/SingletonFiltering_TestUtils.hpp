// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SINGLETON_FILTERING_TEST_UTILS_HPP
#define SINGLETON_FILTERING_TEST_UTILS_HPP

#include <Teuchos_RCP.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_Comm.hpp>
#include <TpetraExt_MatrixMatrix.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_CrsSingletonFilter_LinearProblem.hpp>

// Utilities
// --------------------------------------------------------------------------

template <typename Scalar, typename LO, typename GO, typename Node>
void Display_Matrix(const Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node>>& A,
                    const Teuchos::RCP<const Teuchos::Comm<int>>& Comm,
                    Teuchos::FancyOStream& out,
                    bool printRep, bool contiguous) {
  using crs_local_inds_host_view_type = typename Tpetra::CrsMatrix<Scalar, LO, GO, Node>::local_inds_host_view_type;
  using crs_values_host_view_type     = typename Tpetra::CrsMatrix<Scalar, LO, GO, Node>::values_host_view_type;

  // Get the number of rows and columns
  Tpetra::global_size_t numRows = A->getGlobalNumRows();
  Tpetra::global_size_t numCols = A->getGlobalNumCols();

  const size_t myRank   = Comm->getRank();
  const size_t rankSize = Comm->getSize();

  if (myRank == 0)
    out << "Matrix dimensions: " << numRows << " x " << numCols << std::endl;

  if (myRank == 0)
    out << "Matrix representation (non-zero entries marked with 'X'):" << std::endl;

  if (!printRep && myRank == 0) {
    out << "        ";
    for (Tpetra::global_size_t j = 0; j < numCols; ++j) {
      int globalCol = j;
      if (!contiguous && rankSize == 1) {
        globalCol = A->getColMap()->getGlobalElement(j);
      }
      out << std::setw(12) << "Col " << std::setw(2) << globalCol;
    }
    out << " " << std::endl;
  } else if (printRep && myRank == 0) {
    out << "Col     ";
    for (Tpetra::global_size_t j = 0; j < numCols; ++j) {
      int globalCol = j;
      if (!contiguous && rankSize == 1) {
        globalCol = A->getColMap()->getGlobalElement(j);
      }
      out << std::setw(2) << globalCol;
    }
    out << " " << std::endl;
  }

  for (Tpetra::global_size_t irow = 0; irow < numRows; ++irow) {
    int globalRow = irow;
    if (!contiguous && rankSize == 1) {
      globalRow = A->getRowMap()->getGlobalElement(irow);
    }
    if (A->getRowMap()->getLocalElement(globalRow) != Teuchos::OrdinalTraits<LO>::invalid()) {
      out << "Row " << std::setw(2) << globalRow << " [ ";

      crs_local_inds_host_view_type localIndices;
      crs_values_host_view_type values;

      LO localRow = A->getRowMap()->getLocalElement(globalRow);
      A->getLocalRowView(localRow, localIndices, values);

      std::vector<bool> printed(numCols, false);
      std::vector<Scalar> printValue(numCols, false);

      LO numEntries = A->getNumEntriesInLocalRow(localRow);
      for (LO k = 0; k < numEntries; ++k) {
        int globalIndex = A->getColMap()->getGlobalElement(localIndices(k));
        if (!contiguous && rankSize == 1) {
          LO localColIndex          = A->getColMap()->getLocalElement(globalIndex);
          printed[localColIndex]    = true;
          printValue[localColIndex] = values(k);
        } else {
          printed[globalIndex]    = true;
          printValue[globalIndex] = values(k);
        }
      }

      for (Tpetra::global_size_t j = 0; j < numCols; ++j) {
        if (printed[j]) {
          if (printRep)
            out << "X ";
          else
            out << std::setw(13) << printValue[j] << " ";
        } else {
          if (printRep)
            out << "  ";
          else
            out << std::setw(13) << 0 << " ";
        }
      }

      out << "]" << std::endl;
    }
    Comm->barrier();
  }
  Comm->barrier();
}

template <typename Scalar, typename LO, typename GO, typename Node>
bool compareCrsMatrices(const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar, LO, GO, Node>>& A,
                        const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar, LO, GO, Node>>& B,
                        const Teuchos::RCP<const Teuchos::Comm<int>>& Comm,
                        Teuchos::FancyOStream& out,
                        double relativeTolerance = 1e-12) {
#if KOKKOS_VERSION >= 40799
  using impl_scalar_type = typename KokkosKernels::ArithTraits<Scalar>::val_type;
  using mag_type         = typename KokkosKernels::ArithTraits<impl_scalar_type>::mag_type;
#else
  using impl_scalar_type = typename Kokkos::ArithTraits<Scalar>::val_type;
  using mag_type         = typename Kokkos::ArithTraits<impl_scalar_type>::mag_type;
#endif

  bool comparePass = true;

  Scalar one(1.0);
  Scalar mone(-1.0);
  // Compare global dimensions
  if (A->getGlobalNumRows() != B->getGlobalNumRows() ||
      A->getGlobalNumCols() != B->getGlobalNumCols()) {
    out << "Global dimensions do not match." << std::endl;
    out << "Global number of rows: A(" << A->getGlobalNumRows() << ") != B(" << B->getGlobalNumRows() << ")" << std::endl;
    out << "Global number of cols: A(" << A->getGlobalNumCols() << ") != B(" << B->getGlobalNumCols() << ")" << std::endl;
    return false;
  }

  // Compare global entries
  if (A->getGlobalNumEntries() != B->getGlobalNumEntries()) {
    out << "Global Number of Entries do not match." << std::endl;
    out << "Global number of Entries: A(" << A->getGlobalNumEntries() << ") != B(" << B->getGlobalNumEntries() << ")" << std::endl;
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
    auto diff                        = Teuchos::rcp(new Tpetra::CrsMatrix<Scalar, LO, GO, Node>(A->getRowMap(), A->getColMap(), maxNumEntriesPerRow));

    // Perform the matrix-matrix addition with scalar alpha = 1.0 and beta = -1.0
    Tpetra::MatrixMatrix::Add(
        *A,     // matrix A
        false,  // transpose A (don't transpose)
        one,    // alpha
        *B,     // matrix B
        false,  // transpose B (don't transpose)
        mone,   // beta
        diff);  // output matrix

    // Compute the norm-2 of the difference
    mag_type diffNorm = diff->getFrobeniusNorm();

    // Check if the norm is below the relative tolerance
    mag_type maxNormA = A->getFrobeniusNorm();
    mag_type maxNormB = B->getFrobeniusNorm();
    mag_type maxNorm  = std::max(maxNormA, maxNormB);

    if (maxNorm > 0.0) {  // Avoid division by zero
      mag_type relativeError = diffNorm / maxNorm;
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
                         double relativeTolerance = 1e-12) {
#if KOKKOS_VERSION >= 40799
  using impl_scalar_type = typename KokkosKernels::ArithTraits<Scalar>::val_type;
  using mag_type         = typename KokkosKernels::ArithTraits<impl_scalar_type>::mag_type;
#else
  using impl_scalar_type = typename Kokkos::ArithTraits<Scalar>::val_type;
  using mag_type         = typename Kokkos::ArithTraits<impl_scalar_type>::mag_type;
#endif

  bool comparePass = true;

  // Compare global dimensions
  if (A->getGlobalLength() != B->getGlobalLength() ||
      A->getNumVectors() != B->getNumVectors()) {
    out << "Global dimensions do not match." << std::endl;
    out << "Global Length   : A(" << A->getGlobalLength() << ") != B(" << B->getGlobalLength() << ")" << std::endl;
    out << "Global # vectors: A(" << A->getNumVectors() << ") != B(" << B->getNumVectors() << ")" << std::endl;
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
  Teuchos::Array<mag_type> norms(A->getNumVectors());

  // Compute the 2-norm of the difference for each vector
  diff->norm2(Teuchos::ArrayView<mag_type>(norms));

  if constexpr (!std::is_same<Scalar, long long>::value) {
    // It appears that Reader_t::readSparseFile() reads the file for long long badly when it has real values.
    // Skip comparing values when Scalar is long long.

    // Check if the norms are below the relative tolerance
    for (size_t i = 0; i < static_cast<size_t>(norms.size()); ++i) {
      mag_type maxNormA = A->getVector(i)->norm2();
      mag_type maxNormB = B->getVector(i)->norm2();
      mag_type maxNorm  = std::max(maxNormA, maxNormB);

      if (maxNorm > 0.0) {  // Avoid division by zero
        mag_type relativeError = norms[i] / maxNorm;
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
void test_Singleton_fwd(std::string Matrix_Original_file, std::string LHS_Original_file, std::string RHS_Original_file,
                        std::string Matrix_Reduced_file, std::string LHS_Reduced_file, std::string RHS_Reduced_file,
                        const Teuchos::RCP<const Teuchos::Comm<int>>& Comm,
                        Teuchos::FancyOStream& out, bool& success,
                        double relativeTolerance = 1e-12,
                        bool displayMatrices     = false,
                        bool outputBaseline      = false) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Tpetra::TestingUtilities::getDefaultComm;

  using CrsMatrix_t             = Tpetra::CrsMatrix<Scalar, LO, GO, Node>;
  using MultiVector_t           = Tpetra::MultiVector<Scalar, LO, GO, Node>;
  using Reader_t                = Tpetra::MatrixMarket::Reader<CrsMatrix_t>;
  using Map_t                   = Tpetra::Map<LO, GO, Node>;
  using LinearProblem_t         = Tpetra::LinearProblem<Scalar, LO, GO, Node>;
  using CrsSingletonFiltering_t = Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>;

  RCP<CrsMatrix_t> A_Original;
  RCP<MultiVector_t> LHS_Original, RHS_Original;

  A_Original                      = Reader_t::readSparseFile(Matrix_Original_file, Comm);
  RCP<const Map_t> A_Original_Map = A_Original->getRangeMap();
  LHS_Original                    = Reader_t::readDenseFile(LHS_Original_file, Comm, A_Original_Map);
  RHS_Original                    = Reader_t::readDenseFile(RHS_Original_file, Comm, A_Original_Map);

  bool verbose                             = true;
  bool run_on_host                         = false;
  RCP<MultiVector_t> x                     = rcp(new MultiVector_t(A_Original_Map, LHS_Original->getNumVectors()));
  RCP<LinearProblem_t> preSingletonProblem = rcp(new LinearProblem_t(A_Original, x, RHS_Original));
  CrsSingletonFiltering_t SingletonTransform(run_on_host, verbose);
  RCP<LinearProblem_t> postSingletonProblem = SingletonTransform(preSingletonProblem);

  SingletonTransform.fwd();

  RCP<CrsMatrix_t> A_Reduced;
  RCP<MultiVector_t> LHS_Reduced, RHS_Reduced;

  auto reducedRowMap    = postSingletonProblem->getMatrix()->getRowMap();
  auto reducedColMap    = postSingletonProblem->getMatrix()->getColMap();
  auto reducedDomainMap = postSingletonProblem->getMatrix()->getDomainMap();
  auto reducedRangeMap  = postSingletonProblem->getMatrix()->getRangeMap();
  postSingletonProblem->getMatrix()->getRangeMap()->describe(out, Teuchos::VERB_EXTREME);
  A_Reduced = Reader_t::readSparseFile(Matrix_Reduced_file, reducedRowMap, reducedColMap, reducedDomainMap, reducedRangeMap);

  LHS_Reduced = Reader_t::readDenseFile(LHS_Reduced_file, Comm, reducedDomainMap);
  RHS_Reduced = Reader_t::readDenseFile(RHS_Reduced_file, Comm, reducedRangeMap);

  TEUCHOS_ASSERT(compareCrsMatrices(
      Teuchos::rcp_dynamic_cast<const CrsMatrix_t>(A_Reduced, true),
      Teuchos::rcp_dynamic_cast<const CrsMatrix_t>(postSingletonProblem->getMatrix(), true),
      Comm, out, 1.0e-05));

  TEUCHOS_ASSERT(compareMultiVectors(
      Teuchos::rcp_dynamic_cast<const MultiVector_t>(LHS_Reduced, true),
      Teuchos::rcp_dynamic_cast<const MultiVector_t>(postSingletonProblem->getLHS(), true),
      Comm, out));
  TEUCHOS_ASSERT(compareMultiVectors(
      Teuchos::rcp_dynamic_cast<const MultiVector_t>(RHS_Reduced, true),
      Teuchos::rcp_dynamic_cast<const MultiVector_t>(postSingletonProblem->getRHS(), true),
      Comm, out));

  if (displayMatrices) {
    Display_Matrix(Teuchos::rcp_dynamic_cast<CrsMatrix_t>(A_Original, true), Comm, out, false, true);
    Display_Matrix(Teuchos::rcp_dynamic_cast<CrsMatrix_t>(postSingletonProblem->getMatrix(), true), Comm, out, false, false);
  }

  if (outputBaseline) {
    if constexpr (std::is_same<Scalar, double>::value) {
      using Writer_t     = Tpetra::MatrixMarket::Writer<CrsMatrix_t>;
      auto reducedMatrix = Teuchos::rcp_dynamic_cast<const CrsMatrix_t>(postSingletonProblem->getMatrix(), true);
      Writer_t::writeSparseFile("baseline_" + Matrix_Reduced_file, reducedMatrix);
      Writer_t::writeMapFile("baseline_Map_" + Matrix_Reduced_file, *(reducedMatrix->getRowMap()));
      Writer_t::writeDenseFile("baseline_" + LHS_Reduced_file, postSingletonProblem->getLHS());
      Writer_t::writeDenseFile("baseline_" + RHS_Reduced_file, postSingletonProblem->getRHS());
    }
  }
}

#endif  // SINGLETON_FILTERING_TEST_UTILS_HPP
