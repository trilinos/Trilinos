// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Creates vectors with different maps; tests results of export into them
// Tests behavior of Tpetra::ADD_ASSIGN in Tpetra::MultiVector for many common
// (and a few less common) use cases.

#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Apply_Helpers.hpp"
#include "Teuchos_Array.hpp"
#include <vector>

#include "Teuchos_UnitTestHarness.hpp"
#include "TpetraCore_ETIHelperMacros.h"

namespace {

template<typename map_t, typename vector_t, typename matrix_t>
class MatrixBuilder {

public:  

  using gno_t = typename map_t::global_ordinal_type;
  using scalar_t = typename matrix_t::scalar_type;

  MatrixBuilder(const Teuchos::RCP<const Teuchos::Comm<int> > &comm_, 
                bool square_) :
    comm(comm_),
    nEntriesPerRow(3), 
    nMyRow(4), 
    nGlobalRow(nMyRow * comm_->getSize()),
    yInit(100*(comm_->getRank()+1)),
    squareMatrix(square_)
  { }
  
  // Return number of rows in generated matrix
  size_t nGlobalRows() const { return nGlobalRow; }

  // Return number of columns in generated matrix
  size_t nGlobalCols() const { 
   return (squareMatrix ? nGlobalRow : nGlobalRow + nEntriesPerRow - 1);
  }
  
  // Run tests with combinations of alpha, beta
  int runTests(
    const Teuchos::RCP<const map_t> &domainMap,
    const Teuchos::RCP<const map_t> &rangeMap,
    const char* testName)
  {
    int ierr = 0;

    // Build the matrix
    Teuchos::RCP<matrix_t> Amat = buildMatrix(domainMap, rangeMap);

    // Test with alpha = 0, beta = 0
    ierr += runAlphaBeta(Amat, scalar_t(0), scalar_t(0), testName);

    // Test with alpha != 0, beta = 0
    ierr += runAlphaBeta(Amat, scalar_t(2), scalar_t(0), testName);

    // Test with alpha = 0, beta != 0
    ierr += runAlphaBeta(Amat, scalar_t(0), scalar_t(3), testName);

    // Test with alpha != 0, beta != 0
    ierr += runAlphaBeta(Amat, scalar_t(2), scalar_t(3), testName);

    return ierr;
  }

  // Run Transpose tests with combinations of alpha, beta, transpose
  int runTestsTranspose(
    const Teuchos::RCP<const map_t> &domainMap,
    const Teuchos::RCP<const map_t> &rangeMap,
    const char* testName)
  {
    int ierr = 0;

    // Build the matrix
    Teuchos::RCP<matrix_t> Amat = buildMatrix(domainMap, rangeMap);

    // Test transpose with alpha = 0, beta = 0
    ierr += runAlphaBetaTranspose(Amat, scalar_t(0), scalar_t(0), testName);

    // Test transpose with alpha != 0, beta = 0
    ierr += runAlphaBetaTranspose(Amat, scalar_t(2), scalar_t(0), testName);

    // Test transpose with alpha = 0, beta != 0
    ierr += runAlphaBetaTranspose(Amat, scalar_t(0), scalar_t(3), testName);

    // Test transpose with alpha != 0, beta != 0
    ierr += runAlphaBetaTranspose(Amat, scalar_t(2), scalar_t(3), testName);

    return ierr;
  }

  int runTestsBatched(
    const Teuchos::RCP<const map_t> &domainMap,
    const Teuchos::RCP<const map_t> &rangeMap,
    const char* testName)
  {
    int ierr = 0;

    // Build the matrix
    Teuchos::RCP<matrix_t> Amat = buildMatrix(domainMap, rangeMap);

    // Test with alpha = 0, beta = 0
    ierr += runAlphaBetaBatched(Amat, scalar_t(0), scalar_t(0), testName);

    // Test with alpha != 0, beta = 0
    ierr += runAlphaBetaBatched(Amat, scalar_t(2), scalar_t(0), testName);

    // Test with alpha = 0, beta != 0
    ierr += runAlphaBetaBatched(Amat, scalar_t(0), scalar_t(3), testName);

    // Test with alpha != 0, beta != 0
    ierr += runAlphaBetaBatched(Amat, scalar_t(2), scalar_t(3), testName);

    return ierr;
  }
private:

  // Build non-symmetric matrix with nEntriesPerRow nonzeros (value = 1.) 
  // in each row, nMyRow rows per processor.  
  // Matrix may be square or rectangular.
  Teuchos::RCP<matrix_t>
  buildMatrix(
    const Teuchos::RCP<const map_t> &domainMap, 
    const Teuchos::RCP<const map_t> &rangeMap
  )
  {
    Teuchos::RCP<const map_t> rowMap = rcp(new map_t(nGlobalRow, 0, comm));

    Teuchos::RCP<matrix_t> Amat = rcp(new matrix_t(rowMap, nEntriesPerRow));

    Teuchos::Array<gno_t> nz(nEntriesPerRow);
    Teuchos::Array<scalar_t> val(nEntriesPerRow, scalar_t(1.));

    for (size_t i = 0; i < nMyRow; i++) {
      gno_t gid = rowMap->getGlobalElement(i);
      for (size_t j = 0; j < nEntriesPerRow; j++) 
        nz[j] = (squareMatrix ? (gid + j) % nGlobalRow : (gid + j));
      Amat->insertGlobalValues(gid, nz(), val());
    }

    Amat->fillComplete(domainMap, rangeMap);

    return Amat;
  }

  // Initialize input vector x_gid = gid
  Teuchos::RCP<vector_t> 
  getInputVector(const Teuchos::RCP<const map_t> &map)  const
  {
    Teuchos::RCP<vector_t> vec = rcp(new vector_t(map));

    auto data = vec->getLocalViewHost(Tpetra::Access::ReadWrite);
    for (size_t i = 0; i < vec->getLocalLength(); i++) {
      gno_t gid = map->getGlobalElement(i);
      data(i, 0) = scalar_t(gid);
    }
    return vec;
  }

  // Initialize output vector y_gid = yInit
  Teuchos::RCP<vector_t> 
  getOutputVector(const Teuchos::RCP<const map_t> &map)  const
  {
    Teuchos::RCP<vector_t> vec = rcp(new vector_t(map));

    vec->putScalar(yInit);

    return vec;
  }

  // Compare the result of matrix apply to analytical result
  int checkResult(
    const Teuchos::RCP<vector_t> &vec, scalar_t alpha, scalar_t beta) const
  {
    int ierr = 0;

    //vec->sync_host();
    auto data = vec->getLocalViewHost(Tpetra::Access::ReadOnly);

    for (size_t i = 0; i < vec->getLocalLength(); i++) {
      gno_t gid = vec->getMap()->getGlobalElement(i);

      scalar_t expected(0);
      for (size_t j = 0; j < nEntriesPerRow; j++) 
        expected += (squareMatrix ? (gid + j) % nGlobalRow : (gid + j));
      expected *= alpha;
      expected += beta * yInit;

      if (data(i,0) != expected) ierr++;
    }

    if (ierr > 0) 
      std::cout << comm->getRank() << " HAD " << ierr << " ERRORS" << std::endl;
    return ierr;
  }

  // Compare the result of transpose matrix apply to analytical result
  int checkResultTranspose(
    const Teuchos::RCP<vector_t> &vec, scalar_t alpha, scalar_t beta) const
  {
    int ierr = 0;

    //vec->sync_host();
    auto data = vec->getLocalViewHost(Tpetra::Access::ReadOnly);

    for (size_t i = 0; i < vec->getLocalLength(); i++) {
      ssize_t gid = ssize_t(vec->getMap()->getGlobalElement(i));

      scalar_t expected = 0;
      if (squareMatrix) {
        for (ssize_t j = 0; j < ssize_t(nEntriesPerRow); j++) {
          ssize_t idx = (nGlobalRow + gid - j) % nGlobalRow;
          expected += scalar_t(idx);
        }
      }
      else {
        for (ssize_t j = 0; j < ssize_t(nEntriesPerRow); j++) {
          ssize_t idx = gid - j;
          if (idx >= 0 && idx < ssize_t(nGlobalRow)) 
            expected += scalar_t(idx);
        }
      }
      expected *= alpha;
      expected += yInit * beta;

      if (data(i,0) != expected) ierr++;
    }

    if (ierr > 0) 
      std::cout << comm->getRank() << " HAD " << ierr << " ERRORS" << std::endl;
    return ierr;
  }

  // Compare the result of matrix batchedApply to analytical result
  int checkResultBatched(
    std::vector<vector_t*> &vec, 
    scalar_t alpha, 
    scalar_t beta) const
  {
    int ierr = 0;

    for (size_t v = 0; v < vec.size(); v++) {
      //vec[v]->sync_host();
      auto data = vec[v]->getLocalViewHost(Tpetra::Access::ReadOnly);

      for (size_t i = 0; i < vec[v]->getLocalLength(); i++) {
        gno_t gid = vec[v]->getMap()->getGlobalElement(i);

        scalar_t expected(0);
        for (size_t j = 0; j < nEntriesPerRow; j++) 
          expected += (squareMatrix ? (gid + j) % nGlobalRow : (gid + j));
        expected *= alpha;
        expected += beta * yInit;
  
        if (data(i,0) != expected) ierr++;
      }
    }

    if (ierr > 0) 
      std::cout << comm->getRank() << " HAD " << ierr << " ERRORS" << std::endl;
    return ierr;
  }

  // Apply and check result with a specific alpha, beta
  int runAlphaBeta(
    const Teuchos::RCP<matrix_t> &Amat,
    scalar_t alpha,
    scalar_t beta,
    const char *testName)
  {
    Teuchos::RCP<vector_t> xvec = getInputVector(Amat->getDomainMap());
    Teuchos::RCP<vector_t> yvec = getOutputVector(Amat->getRangeMap());
    
    Amat->apply(*xvec, *yvec, Teuchos::NO_TRANS, alpha, beta);

    return checkResult(yvec, alpha, beta);
  }

  // Apply and check result with a specific alpha, beta, transpose
  int runAlphaBetaTranspose(
    const Teuchos::RCP<matrix_t> &Amat,
    scalar_t alpha,
    scalar_t beta,
    const char *testName)
  {
    Teuchos::RCP<vector_t> xvec = getInputVector(Amat->getRangeMap());
    Teuchos::RCP<vector_t> yvec = getOutputVector(Amat->getDomainMap());

    Amat->apply(*xvec, *yvec, Teuchos::TRANS, alpha, beta);

    return checkResultTranspose(yvec, alpha, beta);
  }

  int runAlphaBetaBatched(
    const Teuchos::RCP<matrix_t> &Amat,
    scalar_t alpha,
    scalar_t beta,
    const char *testName)
  {
    Teuchos::RCP<matrix_t> A1 = Amat;
    Teuchos::RCP<matrix_t> A2 = rcp(new matrix_t(*A1));
    std::vector<matrix_t *> matrices = {A1.get(),A2.get()};

    Teuchos::RCP<vector_t> xvec = getInputVector(Amat->getDomainMap());

    Teuchos::RCP<vector_t> y1 = getOutputVector(Amat->getRangeMap());
    Teuchos::RCP<vector_t> y2 = getOutputVector(Amat->getRangeMap());
    std::vector<vector_t *> yvec = {y1.get(),y2.get()};

    Tpetra::batchedApply(matrices, *xvec, yvec, alpha, beta);

    return checkResultBatched(yvec, alpha, beta);
  }

  const Teuchos::RCP<const Teuchos::Comm<int> > comm;
  const size_t nEntriesPerRow;
  const size_t nMyRow;
  const size_t nGlobalRow;
  const scalar_t yInit;
  const bool squareMatrix;

};

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7745, RectangularDefault,
                                  Scalar, LO, GO, Node)
{
  // Use default Tpetra maps for range and domain maps
  const Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  using map_t = Tpetra::Map<LO,GO,Node>;
  using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
  using matrix_t = Tpetra::CrsMatrix<Scalar,LO,GO,Node>;

  MatrixBuilder<map_t, vector_t, matrix_t> mb(comm, false);

  // Build Trilinos-default range and domain maps
  
  Teuchos::RCP<const map_t> domainMap = 
    rcp(new map_t(mb.nGlobalCols(), 0, comm));

  Teuchos::RCP<const map_t> rangeMap = 
    rcp(new map_t(mb.nGlobalRows(), 0, comm));

  int ierr = mb.runTests(domainMap, rangeMap,
                         "RectangularDefault");

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7745, RectangularDefaultTranspose,
                                  Scalar, LO, GO, Node)
{
  // Use default Tpetra maps for range and domain maps
  const Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  using map_t = Tpetra::Map<LO,GO,Node>;
  using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
  using matrix_t = Tpetra::CrsMatrix<Scalar,LO,GO,Node>;

  MatrixBuilder<map_t, vector_t, matrix_t> mb(comm, false);

  // Build Trilinos-default range and domain maps
  
  Teuchos::RCP<const map_t> domainMap = 
    rcp(new map_t(mb.nGlobalCols(), 0, comm));

  Teuchos::RCP<const map_t> rangeMap = 
    rcp(new map_t(mb.nGlobalRows(), 0, comm));

  int ierr = mb.runTestsTranspose(domainMap, rangeMap,
                                  "RectangularDefaultTranspose");

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7745, SquareDefault,
                                  Scalar, LO, GO, Node)
{
  // Use default Tpetra maps for range and domain maps
  const Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  using map_t = Tpetra::Map<LO,GO,Node>;
  using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
  using matrix_t = Tpetra::CrsMatrix<Scalar,LO,GO,Node>;

  MatrixBuilder<map_t, vector_t, matrix_t> mb(comm, true);

  // Build Trilinos-default range and domain maps
  
  Teuchos::RCP<const map_t> domainMap = 
    rcp(new map_t(mb.nGlobalCols(), 0, comm));

  Teuchos::RCP<const map_t> rangeMap = 
    rcp(new map_t(mb.nGlobalRows(), 0, comm));

  int ierr = mb.runTests(domainMap, rangeMap, "SquareDefault");

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7745, SquareDefaultTranspose,
                                  Scalar, LO, GO, Node)
{
  // Use default Tpetra maps for range and domain maps
  const Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  using map_t = Tpetra::Map<LO,GO,Node>;
  using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
  using matrix_t = Tpetra::CrsMatrix<Scalar,LO,GO,Node>;

  MatrixBuilder<map_t, vector_t, matrix_t> mb(comm, true);

  // Build Trilinos-default range and domain maps
  
  Teuchos::RCP<const map_t> domainMap = 
    rcp(new map_t(mb.nGlobalCols(), 0, comm));

  Teuchos::RCP<const map_t> rangeMap = 
    rcp(new map_t(mb.nGlobalRows(), 0, comm));

  int ierr = mb.runTestsTranspose(domainMap, rangeMap,
                                  "SquareDefaultTranspose");

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7745, RectangularCyclic,
                                  Scalar, LO, GO, Node)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();

  using map_t = Tpetra::Map<LO,GO,Node>;
  using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
  using matrix_t = Tpetra::CrsMatrix<Scalar,LO,GO,Node>;

  MatrixBuilder<map_t, vector_t, matrix_t> mb(comm, false);
  Teuchos::Array<GO> myEntries(mb.nGlobalCols());

  // Build Trilinos-default range and domain maps
  
  // One-to-one cyclic map:  deal out entries like cards

  int nMyEntries = 0;
  for (size_t i = 0; i < mb.nGlobalCols(); i++) {
    if (int(i % np) == me) {
      myEntries[nMyEntries++] = i;
    }
  }

  Tpetra::global_size_t dummy =
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> domainMap = 
           rcp(new map_t(dummy, myEntries(0,nMyEntries), 0, comm));

  // A different one-to-one cyclic map for the rangeMap

  nMyEntries = 0;
  for (size_t i = 0; i < mb.nGlobalRows(); i++) {
    if (int((i+1) % np) == me) {
      myEntries[nMyEntries++] = i;
    }
  }

  dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> rangeMap = 
           rcp(new map_t(dummy, myEntries(0,nMyEntries), 0, comm));

  // Run the test

  int ierr = mb.runTests(domainMap, rangeMap, "RectangularCyclic");

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7745, RectangularCyclicTranspose,
                                  Scalar, LO, GO, Node)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();

  using map_t = Tpetra::Map<LO,GO,Node>;
  using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
  using matrix_t = Tpetra::CrsMatrix<Scalar,LO,GO,Node>;

  MatrixBuilder<map_t, vector_t, matrix_t> mb(comm, false);
  Teuchos::Array<GO> myEntries(mb.nGlobalCols());

  // Build Trilinos-default range and domain maps
  
  // One-to-one cyclic map:  deal out entries like cards

  int nMyEntries = 0;
  for (size_t i = 0; i < mb.nGlobalCols(); i++) {
    if (int(i % np) == me) {
      myEntries[nMyEntries++] = i;
    }
  }

  Tpetra::global_size_t dummy =
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> domainMap = 
           rcp(new map_t(dummy, myEntries(0,nMyEntries), 0, comm));

  // A different one-to-one cyclic map for the rangeMap

  nMyEntries = 0;
  for (size_t i = 0; i < mb.nGlobalRows(); i++) {
    if (int((i+1) % np) == me) {
      myEntries[nMyEntries++] = i;
    }
  }

  dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> rangeMap = 
           rcp(new map_t(dummy, myEntries(0,nMyEntries), 0, comm));

  // Run the test

  int ierr = mb.runTestsTranspose(domainMap, rangeMap,
                                  "RectangularCyclicTranspose");

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7745, SquareCyclic,
                                  Scalar, LO, GO, Node)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();

  using map_t = Tpetra::Map<LO,GO,Node>;
  using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
  using matrix_t = Tpetra::CrsMatrix<Scalar,LO,GO,Node>;

  MatrixBuilder<map_t, vector_t, matrix_t> mb(comm, true);
  Teuchos::Array<GO> myEntries(mb.nGlobalCols());

  // Build Trilinos-default range and domain maps
  
  // One-to-one cyclic map:  deal out entries like cards

  int nMyEntries = 0;
  for (size_t i = 0; i < mb.nGlobalCols(); i++) {
    if (int(i % np) == me) {
      myEntries[nMyEntries++] = i;
    }
  }

  Tpetra::global_size_t dummy =
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> domainMap = 
           rcp(new map_t(dummy, myEntries(0,nMyEntries), 0, comm));

  // A different one-to-one cyclic map for the rangeMap

  nMyEntries = 0;
  for (size_t i = 0; i < mb.nGlobalRows(); i++) {
    if (int((i+1) % np) == me) {
      myEntries[nMyEntries++] = i;
    }
  }

  dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> rangeMap = 
           rcp(new map_t(dummy, myEntries(0,nMyEntries), 0, comm));

  // Run the test

  int ierr = mb.runTests(domainMap, rangeMap, "SquareCyclic");

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7745, SquareCyclicTranspose,
                                  Scalar, LO, GO, Node)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();

  using map_t = Tpetra::Map<LO,GO,Node>;
  using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
  using matrix_t = Tpetra::CrsMatrix<Scalar,LO,GO,Node>;

  MatrixBuilder<map_t, vector_t, matrix_t> mb(comm, true);
  Teuchos::Array<GO> myEntries(mb.nGlobalCols());

  // Build Trilinos-default range and domain maps
  
  // One-to-one cyclic map:  deal out entries like cards

  int nMyEntries = 0;
  for (size_t i = 0; i < mb.nGlobalCols(); i++) {
    if (int(i % np) == me) {
      myEntries[nMyEntries++] = i;
    }
  }

  Tpetra::global_size_t dummy =
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> domainMap = 
           rcp(new map_t(dummy, myEntries(0,nMyEntries), 0, comm));

  // A different one-to-one cyclic map for the rangeMap

  nMyEntries = 0;
  for (size_t i = 0; i < mb.nGlobalRows(); i++) {
    if (int((i+1) % np) == me) {
      myEntries[nMyEntries++] = i;
    }
  }

  dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> rangeMap = 
           rcp(new map_t(dummy, myEntries(0,nMyEntries), 0, comm));

  // Run the test

  int ierr = mb.runTestsTranspose(domainMap, rangeMap,
                                  "SquareCyclicTranspose");

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}
//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7745, SquareCyclicBatched,
                                  Scalar, LO, GO, Node)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();

  using map_t = Tpetra::Map<LO,GO,Node>;
  using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
  using matrix_t = Tpetra::CrsMatrix<Scalar,LO,GO,Node>;

  MatrixBuilder<map_t, vector_t, matrix_t> mb(comm, true);
  Teuchos::Array<GO> myEntries(mb.nGlobalCols());

  // Build Trilinos-default range and domain maps
  
  // One-to-one cyclic map:  deal out entries like cards

  int nMyEntries = 0;
  for (size_t i = 0; i < mb.nGlobalCols(); i++) {
    if (int(i % np) == me) {
      myEntries[nMyEntries++] = i;
    }
  }

  Tpetra::global_size_t dummy =
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> domainMap = 
           rcp(new map_t(dummy, myEntries(0,nMyEntries), 0, comm));

  // A different one-to-one cyclic map for the rangeMap

  nMyEntries = 0;
  for (size_t i = 0; i < mb.nGlobalRows(); i++) {
    if (int((i+1) % np) == me) {
      myEntries[nMyEntries++] = i;
    }
  }

  dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> rangeMap = 
           rcp(new map_t(dummy, myEntries(0,nMyEntries), 0, comm));

  // Run the test

  int ierr = mb.runTestsBatched(domainMap, rangeMap, "SquareCyclicBatched");

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7745, RectangularDefault, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7745, RectangularDefaultTranspose, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7745, SquareDefault, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7745, SquareDefaultTranspose, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7745, RectangularCyclic, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7745, RectangularCyclicTranspose, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7745, SquareCyclic, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7745, SquareCyclicTranspose, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7745, SquareCyclicBatched, SCALAR, LO, GO, NODE) 

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )

} // namespace (anonymous)

