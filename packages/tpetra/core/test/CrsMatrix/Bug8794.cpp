// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "TpetraCore_ETIHelperMacros.h"


namespace {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug8794, InsertDenseRows,
                                  Scalar, LO, GO, Node)
{
// Test for issue #8794
// Build a matrix using insertGlobalValues
// The matrix will have some sparse rows (number of nonzeros <= 5) and
// some dense rows (number of nonzeros = 501).
// The two implementations of insert_crs_indices that differ depending 
// on the number of indices being inserted are tested.
// Multiply the matrix time a vector of global IDs and compare the result
// to expected values.

  using map_t = Tpetra::Map<>;
  using matrix_t = Tpetra::CrsMatrix<Scalar>;
  using vector_t = Tpetra::Vector<Scalar>;
  using magnitude_t = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;

  auto comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();

  int nrows = 50001;
  int divisor = 100;
  int maxNzPerRow = nrows / divisor + 1;

  // Map with rows across all processors
  Teuchos::RCP<const map_t> map = rcp(new map_t(nrows, 0, comm));

  // Vectors for SpMV and expected values
  vector_t expected(map);
  vector_t x(map);
  vector_t y(map);

  // Build matrix distributed across np-1 processors:  
  // insert nonzeros on np processors; 
  // let fillComplete migrate according to mapNpM1

  Teuchos::Array<GO> cols(maxNzPerRow);
  Teuchos::Array<Scalar> vals(maxNzPerRow, 1.);

  matrix_t Amat(map, maxNzPerRow);

  // Initialize matrix and expected value of SpMV product
  {
    expected.putScalar(0.);
    auto expectedData = Kokkos::subview(expected.getLocalViewHost(Tpetra::Access::ReadWrite), Kokkos::ALL(), 0);
    for (size_t i = 0; i < map->getLocalNumElements(); i++) {

      GO gid = map->getGlobalElement(i);
      bool denseRow = (gid % (divisor+1) == 1);

      if (!denseRow) {  // sparse row
        int nz = 0;
        cols[nz++] = gid;
        expectedData[i] += Scalar(gid);
        if (gid+1<nrows) {cols[nz++] = gid+1; expectedData[i] += Scalar(gid+1);}
        if (gid+2<nrows) {cols[nz++] = gid+2; expectedData[i] += Scalar(gid+2);}
        if (gid-1>=0)    {cols[nz++] = gid-1; expectedData[i] += Scalar(gid-1);}
        if (gid-2>=0)    {cols[nz++] = gid-2; expectedData[i] += Scalar(gid-2);}
        Amat.insertGlobalValues(gid, cols(0,nz), vals(0,nz));
      }
      else { // dense row
        if (gid % 2) {  // Insert many nonzeros all at once
          int nz = 0;
          for (int j = 0; j < nrows; j+=divisor) {
            GO tmp = (gid + j) % nrows;
            cols[nz++] = tmp;
            expectedData[i] += Scalar(tmp);
          }
          Amat.insertGlobalValues(gid, cols(0,nz), vals(0,nz));
        }
        else {  // Insert many nonzeros in batches, with some duplicates
          int nz = 0;
          int nrowsCutoff = nrows - 10*divisor;

          for (int j = 0; j < nrowsCutoff; j+=divisor) {
            GO tmp = (gid + j) % nrows;
            cols[nz++] = tmp;
            expectedData[i] += Scalar(tmp);
          }
          Amat.insertGlobalValues(gid, cols(0,nz), vals(0,nz));

          nz = 0;
          for (int j = 0; j < nrows; j+=divisor) {
            GO tmp = (gid + j) % nrows;
            cols[nz++] = tmp;
            expectedData[i] += Scalar(tmp);
          }
          Amat.insertGlobalValues(gid, cols(0,nz), vals(0,nz));
        }
      }
    }

    Amat.fillComplete();
    out << me << " of " << np << ": \n"
        << "  nrows     " << Amat.getLocalNumRows() << "\n"
        << "  nnz       " << Amat.getLocalNumEntries() << "\n"
        << "  maxPerRow " << Amat.getLocalMaxNumRowEntries() << "\n"
        << "  norm      " << Amat.getFrobeniusNorm() << "\n"
        << std::endl;
  }

  // Initialize domain vector for SpMV
  {
    auto xData = x.getLocalViewHost(Tpetra::Access::OverwriteAll);
    for (size_t i = 0; i < map->getLocalNumElements(); i++)
      xData(i, 0) = map->getGlobalElement(i);
  }

  Amat.apply(x, y);

  // Test product for correctness
  int ierr = 0;
  {
    auto expectedData = expected.getLocalViewHost(Tpetra::Access::ReadOnly);
    auto yData = y.getLocalViewHost(Tpetra::Access::ReadOnly);

    for (size_t i = 0; i < map->getLocalNumElements(); i++) {
      if (Teuchos::ScalarTraits<Scalar>::magnitude(yData(i, 0)-expectedData(i, 0))
          > 200*Teuchos::ScalarTraits<Scalar>::magnitude(expectedData(i, 0))*Teuchos::ScalarTraits<magnitude_t>::eps()) {
        out << Teuchos::ScalarTraits<Scalar>::magnitude(yData(i, 0)-expectedData(i, 0)) << " " << 100*Teuchos::ScalarTraits<Scalar>::magnitude(expectedData(i, 0))*Teuchos::ScalarTraits<magnitude_t>::eps()<< std::endl;
        out << me << " of " << np << ": y[" << map->getGlobalElement(i)
            << "] " << yData(i, 0) << " != " << expectedData(i, 0)
            << " expected" << std::endl;
        ierr++;
      }
    }
  }

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);
  if (gerr) out << "TEST FAILED with " << gerr << " errors" << std::endl;
  else out << "TEST PASSED" << std::endl;

  TEST_ASSERT(gerr == 0);
}

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug8794, InsertDenseRows, SCALAR, LO, GO, NODE) 

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )
}
