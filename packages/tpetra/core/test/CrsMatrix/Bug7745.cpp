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
// ************************************************************************
// @HEADER
*/

// Creates vectors with different maps; tests results of export into them
// Tests behavior of Tpetra::ADD_ASSIGN in Tpetra::MultiVector for many common
// (and a few less common) use cases.

#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Teuchos_Array.hpp"

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
    squareMatrix(square_),
    foo(Teuchos::rcp(&std::cout,false))
  { }
  
  // Return number of rows in generated matrix
  size_t nGlobalRows() const { return nGlobalRow; }

  // Return number of columns in generated matrix
  size_t nGlobalCols() const { 
   return (squareMatrix ? nGlobalRow : nGlobalRow + nEntriesPerRow - 1);
  }
  
  // Run tests with combinations of alpha, beta, transpose
  int runTests(
    const Teuchos::RCP<const map_t> &domainMap,
    const Teuchos::RCP<const map_t> &rangeMap,
    const char* testName)
  {
    int me = comm->getRank();
    int ierr = 0;

    // Build the matrix
    Teuchos::RCP<const matrix_t> Amat = buildMatrix(domainMap, rangeMap);

    // Test with alpha = 0, beta = 0
    ierr += runAlphaBeta(Amat, scalar_t(0), scalar_t(0), testName);

    // Test with alpha != 0, beta = 0
    ierr += runAlphaBeta(Amat, scalar_t(2), scalar_t(0), testName);

    // Test with alpha = 0, beta != 0
    ierr += runAlphaBeta(Amat, scalar_t(0), scalar_t(3), testName);

    // Test with alpha != 0, beta != 0
    ierr += runAlphaBeta(Amat, scalar_t(2), scalar_t(3), testName);

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

private:

  // Build non-symmetric matrix with nEntriesPerRow nonzeros (value = 1.) 
  // in each row, nMyRow rows per processor.  
  // Matrix may be square or rectangular.
  Teuchos::RCP<const matrix_t>
  buildMatrix(
    const Teuchos::RCP<const map_t> &domainMap, 
    const Teuchos::RCP<const map_t> &rangeMap
  )
  {
    int me = comm->getRank();

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
    std::cout << me << " MATRIX " << std::endl;
    Amat->describe(foo, Teuchos::VERB_EXTREME);

    return Amat;
  }

  // Initialize domain vector x_gid = gid
  Teuchos::RCP<vector_t> 
  getDomainVector(const Teuchos::RCP<const map_t> &domainMap)  const
  {
    Teuchos::RCP<vector_t> vec = rcp(new vector_t(domainMap));

    auto data = vec->getLocalViewHost();
    for (size_t i = 0; i < vec->getLocalLength(); i++) {
      gno_t gid = domainMap->getGlobalElement(i);
      data(i, 0) = scalar_t(gid);
    }
    return vec;
  }

  // Initialize range vector y_gid = yInit
  Teuchos::RCP<vector_t> 
  getRangeVector(const Teuchos::RCP<const map_t> &rangeMap)  const
  {
    Teuchos::RCP<vector_t> vec = rcp(new vector_t(rangeMap));

    vec->putScalar(yInit);

    return vec;
  }

  // Compare the result of matrix apply to analytical result
  int checkResult(
    const Teuchos::RCP<vector_t> &vec, scalar_t alpha, scalar_t beta) const
  {
    int ierr = 0;

    vec->sync_host();
    auto data = vec->getLocalViewHost();

    for (size_t i = 0; i < vec->getLocalLength(); i++) {
      gno_t gid = vec->getMap()->getGlobalElement(i);

      scalar_t expected(0);
      for (size_t j = 0; j < nEntriesPerRow; j++) 
        expected += (squareMatrix ? (gid + j) % nGlobalRow : (gid + j));
      expected *= alpha;
      expected += beta * yInit;

      if (data(i,0) != expected) ierr++;
    }

    return ierr;
  }

  // Compare the result of transpose matrix apply to analytical result
  int checkResultTranspose(
    const Teuchos::RCP<vector_t> &vec, scalar_t alpha, scalar_t beta) const
  {
    int ierr = 0;

    vec->sync_host();
    auto data = vec->getLocalViewHost();

    for (size_t i = 0; i < vec->getLocalLength(); i++) {
      gno_t gid = vec->getMap()->getGlobalElement(i);

      size_t mult;
      if (squareMatrix) {
        mult = nEntriesPerRow;
      }
      else {
        if (gid < nEntriesPerRow) 
          mult = gid+1;
        else if (gid > nGlobalCols() - nEntriesPerRow) 
          mult = nGlobalCols()-gid;
        else 
          mult = nEntriesPerRow;
      }
      scalar_t expected = scalar_t(mult) * yInit * alpha;
      expected += scalar_t(gid) * beta;

      if (data(i,0) != expected) ierr++;
    }

    return ierr;
  }

  // Apply and check result with a specific alpha, beta
  int runAlphaBeta(
    const Teuchos::RCP<const matrix_t> &Amat,
    scalar_t alpha,
    scalar_t beta,
    const char *testName)
  {
    int me = comm->getRank();

    Teuchos::RCP<vector_t> xvec = getDomainVector(Amat->getDomainMap());
    Teuchos::RCP<vector_t> yvec = getRangeVector(Amat->getRangeMap());
    
    Amat->apply(*xvec, *yvec, Teuchos::NO_TRANS, alpha, beta);

    std::cout << me << " " << testName 
              << " YVEC alpha=" << alpha 
              << " beta=" << beta << std::endl;
    yvec->describe(foo, Teuchos::VERB_EXTREME);
  
    return checkResult(yvec, alpha, beta);
  }

  // Apply and check result with a specific alpha, beta, transpose
  int runAlphaBetaTranspose(
    const Teuchos::RCP<const matrix_t> &Amat,
    scalar_t alpha,
    scalar_t beta,
    const char *testName)
  {
    int me = comm->getRank();

    Teuchos::RCP<vector_t> xvec = getRangeVector(Amat->getRangeMap());
    Teuchos::RCP<vector_t> yvec = getDomainVector(Amat->getDomainMap());

    Amat->apply(*xvec, *yvec, Teuchos::TRANS, alpha, beta);

    std::cout << me << " " << testName 
              << " YVEC Transpose alpha=" << alpha
              << " beta=" << beta << std::endl;
    yvec->describe(foo, Teuchos::VERB_EXTREME);
  
    return checkResultTranspose(yvec, alpha, beta);
  }


  const Teuchos::RCP<const Teuchos::Comm<int> > comm;
  const size_t nEntriesPerRow;
  const size_t nMyRow;
  const size_t nGlobalRow;
  const scalar_t yInit;
  const bool squareMatrix;
  Teuchos::FancyOStream foo;

};

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7745, RectangularDefaultToDefault,
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

  int ierr = mb.runTests(domainMap, rangeMap, "RectangularDefaultToDefault");

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7745, SquareDefaultToDefault,
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

  int ierr = mb.runTests(domainMap, rangeMap, "SquareDefaultToDefault");

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}
#ifdef KDD
//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7745, CyclicToDefault, Scalar,LO,GO,Node)
{
  // This case demonstrates that owned entries shared between the source and
  // target map are copied from the source vector into the target (during
  // copyAndPermute).  Owned entries that are not in the source map
  // are NOT reset; their initial values persist.  Then received shared 
  // entries are added to the owned entries.
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();
  int ierr = 0;

  Teuchos::FancyOStream foo(Teuchos::rcp(&std::cout,false));

  using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
  using map_t = Tpetra::Map<LO,GO,Node>;

  const size_t nGlobalEntries = 8 * np;
  const Scalar tgtScalar = 100. * (me+1);
  const Scalar srcScalar = 2.;
  Teuchos::Array<GO> myEntries(nGlobalEntries); 

  // Default one-to-one linear block map in Trilinos

  Teuchos::RCP<const map_t> defaultMap = 
           rcp(new map_t(nGlobalEntries, 0, comm));

  std::cout << me << " DEFAULT MAP" << std::endl;
  defaultMap->describe(foo, Teuchos::VERB_EXTREME);

  // One-to-one cyclic map:  deal out entries like cards

  int nMyEntries = 0;
  for (size_t i = 0; i < nGlobalEntries; i++) {
    if (int(i % np) == me) {
      myEntries[nMyEntries++] = i;
    }
  }

  Tpetra::global_size_t dummy =
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> cyclicMap = 
           rcp(new map_t(dummy, myEntries(0,nMyEntries), 0, comm));

  std::cout << me << " CYCLIC MAP" << std::endl;
  cyclicMap->describe(foo, Teuchos::VERB_EXTREME);

  // Create vectors

  vector_t defaultVecTgt(defaultMap);
  defaultVecTgt.putScalar(tgtScalar);

  vector_t cyclicVecSrc(cyclicMap);
  cyclicVecSrc.putScalar(srcScalar);

  // Export Cyclic-to-default

  Tpetra::Export<LO,GO,Node> cyclicToDefault(cyclicMap, defaultMap);
  defaultVecTgt.doExport(cyclicVecSrc, cyclicToDefault, Tpetra::ADD_ASSIGN);

  std::cout << me << " CYCLIC TO DEFAULT " << std::endl;
  defaultVecTgt.describe(foo, Teuchos::VERB_EXTREME);

  // Check result

  auto invalid = Teuchos::OrdinalTraits<LO>::invalid();
  auto data = defaultVecTgt.getLocalViewHost();
  for (size_t i = 0; i < defaultVecTgt.getLocalLength(); i++)
    if (data(i,0) != tgtScalar + srcScalar) ierr++;
  if (ierr > 0) 
    std::cout << "TEST FAILED:  CYCLIC-TO-DEFAULT TEST HAD " << ierr 
              << " FAILURES ON RANK " << me << std::endl;

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7745, OverlapToDefault, Scalar,LO,GO,Node)
{
  // This case demonstrates that owned entries shared between the source and
  // target map are copied from the source vector into the target (during
  // copyAndPermute).  Owned entries that are not in the source map
  // are NOT reset; their initial values persist.  Then received shared 
  // entries are added to the owned entries.
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();
  int ierr = 0;

  if (np > 1) {  // Need more than one proc to avoid duplicate entries in maps
    Teuchos::FancyOStream foo(Teuchos::rcp(&std::cout,false));

    using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
    using map_t = Tpetra::Map<LO,GO,Node>;

    const size_t nGlobalEntries = 8 * np;
    const Scalar tgtScalar = 100. * (me+1);
    const Scalar srcScalar = 2.;
    Teuchos::Array<GO> myEntries(nGlobalEntries); 

    // Default one-to-one linear block map in Trilinos

    Teuchos::RCP<const map_t> defaultMap = 
             rcp(new map_t(nGlobalEntries, 0, comm));

    std::cout << me << " DEFAULT MAP" << std::endl;
    defaultMap->describe(foo, Teuchos::VERB_EXTREME);

    // Overlap map; some entries are stored on two procs
    int nMyEntries = 0;
    for (size_t i = 0; i < defaultMap->getNodeNumElements()/2; i++) {
      myEntries[nMyEntries++] = defaultMap->getGlobalElement(i);
    }
    for (size_t i = 0; i < defaultMap->getNodeNumElements(); i++) {
      myEntries[nMyEntries++] =
        (defaultMap->getMaxGlobalIndex() + 1 + i) % nGlobalEntries;
    }
  
    Tpetra::global_size_t dummy = 
            Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    Teuchos::RCP<const map_t> overlapMap = 
             rcp(new map_t(dummy, myEntries(0,nMyEntries), 0, comm));
  
    std::cout << me << " OVERLAP MAP" << std::endl;
    overlapMap->describe(foo, Teuchos::VERB_EXTREME);

    // Create vectors

    vector_t defaultVecTgt(defaultMap);
    defaultVecTgt.putScalar(tgtScalar);

    vector_t overlapVecSrc(overlapMap);
    overlapVecSrc.putScalar(srcScalar);

    // Export Overlap-to-default

    Tpetra::Export<LO,GO,Node> overlapToDefault(overlapMap, defaultMap);
    defaultVecTgt.doExport(overlapVecSrc, overlapToDefault, Tpetra::ADD_ASSIGN);

    std::cout << me << " OVERLAP TO DEFAULT " << std::endl;
    defaultVecTgt.describe(foo, Teuchos::VERB_EXTREME);

    auto data = defaultVecTgt.getLocalViewHost();
    for (size_t i = 0; i < defaultVecTgt.getLocalLength()/2; i++) {
      // overlapped; initial target values were overwritten
      if (data(i,0) != tgtScalar + srcScalar + srcScalar) ierr++;  
    }
    for (size_t i = defaultVecTgt.getLocalLength()/2;
             i < defaultVecTgt.getLocalLength(); i++) {
      // not overlapped; initial target values were not overwritten
      if (data(i,0) != tgtScalar + srcScalar) ierr++;  
    }
    if (ierr > 0) 
      std::cout << "TEST FAILED:  OVERLAP-TO-DEFAULT TEST HAD " << ierr 
                << " FAILURES ON RANK " << me << std::endl;
  }

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7745, OddEvenToSerial, Scalar,LO,GO,Node)
{
  // Test case showing behavior when target map is all on processor zero.
  // In the source map, even numbered entries are on even numbered processors;
  // odd numbered entreis are on odd numbered processors.
  // In copyAndPermute, even numbered entries are copied from processor zero's
  // source vector to the target vector, and odd numbered entries are unchanged.
  // Then received values are added to the target vector.  The result is that
  // odd entries include the initial target values in their sum, while the 
  // even entries are not included.
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();
  int ierr = 0;

  Teuchos::FancyOStream foo(Teuchos::rcp(&std::cout,false));

  using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
  using map_t = Tpetra::Map<LO,GO,Node>;

  const size_t nGlobalEntries = 8 * np;
  const Scalar tgtScalar = 100. * (me+1);
  const Scalar srcScalar = 2.;
  Teuchos::Array<GO> myEntries(nGlobalEntries); 

  // Odd entries given to odd procs; even entries given to even procs
  int nMyEntries = 0;
  for (size_t i = 0; i < nGlobalEntries; i++) {
    if (int(i % 2) == (me % 2)) {
      myEntries[nMyEntries++] = i;
    }
  }

  Tpetra::global_size_t dummy =
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> oddEvenMap = 
           rcp(new map_t(dummy, myEntries(0,nMyEntries), 0, comm));

  std::cout << me << " ODDEVEN MAP" << std::endl;
  oddEvenMap->describe(foo, Teuchos::VERB_EXTREME);

  // Map with all entries on one processor

  dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  size_t nSerialEntries = (me == 0 ? nGlobalEntries : 0);
  Teuchos::RCP<const map_t> serialMap = 
           rcp(new map_t(dummy, nSerialEntries, 0, comm));

  std::cout << me << " SERIAL MAP" << std::endl;
  serialMap->describe(foo, Teuchos::VERB_EXTREME);

  // Create vectors

  vector_t oddEvenVecSrc(oddEvenMap);
  oddEvenVecSrc.putScalar(srcScalar);

  vector_t serialVecTgt(serialMap);
  serialVecTgt.putScalar(tgtScalar);

  // Export oddEven-to-serial

  Tpetra::Export<LO,GO,Node> oddEvenToSerial(oddEvenMap, serialMap);
  serialVecTgt.doExport(oddEvenVecSrc, oddEvenToSerial, Tpetra::ADD_ASSIGN);

  std::cout << me << " ODDEVEN TO SERIAL " << std::endl;
  serialVecTgt.describe(foo, Teuchos::VERB_EXTREME);

  // Check result

  auto data = serialVecTgt.getLocalViewHost();
  for (size_t i = 0; i < serialVecTgt.getLocalLength(); i++) {
    Scalar nCopies = Scalar(((np+1) / 2) - ((i % 2 == 1) && (np % 2 == 1)));
    if (data(i,0) != tgtScalar + srcScalar * nCopies)
      ierr++;
  }
  if (ierr > 0) 
    std::cout << "TEST FAILED:  ODDEVEN-TO-SERIAL TEST HAD " << ierr 
              << " FAILURES ON RANK " << me << std::endl;

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7745, SupersetToDefault, Scalar,LO,GO,Node)
{
  // This use case is similar to matrix assembly case in which user 
  // has a map of owned entries and a map of owned+shared entries, with the
  // owned+shared map being a superset of the owned map.  In this case, 
  // the owned values in the owned+shared vector are copied into the owned
  // vector in copyAndPermute; then received shared entries are 
  // added to the owned entries' values.
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();
  int ierr = 0;

  if (np > 1) {  // Need more than one proc to avoid duplicate entries in maps
    Teuchos::FancyOStream foo(Teuchos::rcp(&std::cout,false));

    using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
    using map_t = Tpetra::Map<LO,GO,Node>;

    const size_t nGlobalEntries = 8 * np;
    const Scalar tgtScalar = 100. * (me+1);
    const Scalar srcScalar = 2.;
    Teuchos::Array<GO> myEntries(nGlobalEntries); 

    // Default one-to-one linear block map in Trilinos

    Teuchos::RCP<const map_t> defaultMap = 
             rcp(new map_t(nGlobalEntries, 0, comm));

    std::cout << me << " DEFAULT MAP" << std::endl;
    defaultMap->describe(foo, Teuchos::VERB_EXTREME);

    // Superset map; some entries are stored on two procs
    int nMyEntries = 0;
    for (size_t i = 0; i < defaultMap->getNodeNumElements(); i++) {
      myEntries[nMyEntries++] = defaultMap->getGlobalElement(i);
    }
    for (size_t i = 0; i < defaultMap->getNodeNumElements()/2; i++) {
      myEntries[nMyEntries++] =
        (defaultMap->getMaxGlobalIndex() + 1 + i) % nGlobalEntries;
    }
  
    Tpetra::global_size_t dummy = 
            Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    Teuchos::RCP<const map_t> supersetMap =
             rcp(new map_t(dummy, myEntries(0,nMyEntries), 0, comm));
  
    std::cout << me << " SUPERSET MAP" << std::endl;
    supersetMap->describe(foo, Teuchos::VERB_EXTREME);

    // Create vectors

    vector_t defaultVecTgt(defaultMap);
    defaultVecTgt.putScalar(tgtScalar);

    vector_t supersetVecSrc(supersetMap);
    supersetVecSrc.putScalar(srcScalar);

    // Export Superset-to-default

    Tpetra::Export<LO,GO,Node> supersetToDefault(supersetMap, defaultMap);
    defaultVecTgt.doExport(supersetVecSrc, supersetToDefault,
                           Tpetra::ADD_ASSIGN);

    std::cout << me << " SUPERSET TO DEFAULT " << std::endl;
    defaultVecTgt.describe(foo, Teuchos::VERB_EXTREME);

    auto data = defaultVecTgt.getLocalViewHost();
    for (size_t i = 0; i < defaultVecTgt.getLocalLength()/2; i++)
      if (data(i,0) != tgtScalar + srcScalar + srcScalar) ierr++;  // overlapped
    for (size_t i = defaultVecTgt.getLocalLength()/2;
             i < defaultVecTgt.getLocalLength(); i++)
      if (data(i,0) != tgtScalar + srcScalar) ierr++;  // not overlapped
    if (ierr > 0) 
      std::cout << "TEST FAILED:  SUPERSET-TO-DEFAULT TEST HAD " << ierr 
                << " FAILURES ON RANK " << me << std::endl;
  }

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7745, NoSamesToDefault, Scalar,LO,GO,Node)
{
  // This use case is similar to matrix assembly case in which user 
  // has a map of owned entries and a map of shared entries, with no
  // overlap between the maps.  In this case, received shared entries are 
  // added to the owned entries' values, as copyAndPermute is never called.
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();
  int ierr = 0;

  if (np > 1) {  // Need more than one proc to avoid duplicate entries in maps
    Teuchos::FancyOStream foo(Teuchos::rcp(&std::cout,false));

    using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
    using map_t = Tpetra::Map<LO,GO,Node>;

    const size_t nGlobalEntries = 8 * np;
    const Scalar tgtScalar = 100. * (me+1);
    const Scalar srcScalar = 2.;
    Teuchos::Array<GO> myEntries(nGlobalEntries); 

    // Default one-to-one linear block map in Trilinos

    Teuchos::RCP<const map_t> defaultMap = 
             rcp(new map_t(nGlobalEntries, 0, comm));

    std::cout << me << " DEFAULT MAP" << std::endl;
    defaultMap->describe(foo, Teuchos::VERB_EXTREME);

    // Map with no sames or permutes
    int nMyEntries = 0;
    for (size_t i = 0; i < defaultMap->getNodeNumElements(); i++) {
      myEntries[nMyEntries++] =
        (defaultMap->getMaxGlobalIndex() + 1 + i) % nGlobalEntries;
    }
  
    Tpetra::global_size_t dummy = 
            Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    Teuchos::RCP<const map_t> noSamesMap = 
             rcp(new map_t(dummy, myEntries(0,nMyEntries), 0, comm));
  
    std::cout << me << " NOSAMES MAP" << std::endl;
    noSamesMap->describe(foo, Teuchos::VERB_EXTREME);

    // Create vectors

    vector_t defaultVecTgt(defaultMap);
    defaultVecTgt.putScalar(tgtScalar);

    vector_t noSamesVecSrc(noSamesMap);
    noSamesVecSrc.putScalar(srcScalar);

    // Export noSames-to-default

    Tpetra::Export<LO,GO,Node> noSamesToDefault(noSamesMap, defaultMap);
    defaultVecTgt.doExport(noSamesVecSrc, noSamesToDefault, Tpetra::ADD_ASSIGN);

    std::cout << me << " NOSAMES TO DEFAULT " << std::endl;
    defaultVecTgt.describe(foo, Teuchos::VERB_EXTREME);

    auto data = defaultVecTgt.getLocalViewHost();
    for (size_t i = 0; i < defaultVecTgt.getLocalLength(); i++)
      if (data(i,0) != tgtScalar + srcScalar) ierr++;  
    if (ierr > 0) 
      std::cout << "TEST FAILED:  NOSAMES-TO-DEFAULT TEST HAD " << ierr 
                << " FAILURES ON RANK " << me << std::endl;
  }

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}
#endif //KDD

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7745, RectangularDefaultToDefault, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7745, SquareDefaultToDefault, SCALAR, LO, GO, NODE)
#ifdef KDD 
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7745, CyclicToDefault, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7745, OverlapToDefault, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7745, OddEvenToSerial, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7745, SupersetToDefault, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7745, NoSamesToDefault, SCALAR, LO, GO, NODE)
#endif // KDD
  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )

} // namespace (anonymous)

