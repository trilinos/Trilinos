// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This test came about because of differences between openmpi and mpich
// implementations. Tpetra had hard coded an openmpi specific flag, causing
// crashes with mpich in some specific situations. See
// https://software.sandia.gov/bugzilla/show_bug.cgi?id=6069

#include <Tpetra_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <TpetraCore_ETIHelperMacros.h>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_TestingUtilities.hpp>

// typedef long global_ordinal_type;  //<<<<<<<<   valgrind is clean
// typedef int global_ordinal_type;   //<<<<<<<<   valgrind complains

namespace { // anonymous

template <class GO>
void GetNeighboursCartesian2d(const GO i, const GO nx, const GO ny,
                              GO& left, GO& right, GO& lower, GO& upper)
{
  GO ix, iy;
  ix = i % nx;
  iy = (i - ix) / nx;

  if (ix == 0)      left  = -1;
  else              left  = i - 1;
  if (ix == nx - 1) right = -1;
  else              right = i + 1;
  if (iy == 0)      lower = -1;
  else              lower = i - nx;
  if (iy == ny - 1) upper = -1;
  else              upper = i + nx;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, Bug6069_1, SC, LO, GO, NT)
{

  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Comm;
  using Teuchos::FancyOStream;
  using Tpetra::TestingUtilities::getDefaultComm;

  // Put these typedefs here, to avoid global shadowing warnings.
  typedef Tpetra::Map<LO, GO, NT> MapType;
  typedef Tpetra::CrsMatrix<SC, LO, GO, NT> MatrixType;

  RCP<const Comm<int> > comm = getDefaultComm();
  const size_t myRank = comm->getRank ();
  const size_t numProc = comm->getSize ();

  if (numProc != 3) {
    out << "This test must be run with exactly 3 MPI processes, but you ran "
        << "it with " << numProc << " process" << (numProc != 1 ? "es" : "")
        << "." << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                               "This test must be run with exactly "
                               "3 MPI processes, but you ran it with "
                               << numProc << " process"
                               << (numProc != 1 ? "es" : "") << ".");
  }

  GO nx=10, ny=10;
  SC a=4,b=-1,c=-1,d=-1,e=-1;

  Array<GO> myElts;

  switch (myRank) {

  case 0:
    //numMyElts=40;
    for (int i=0; i<10; ++i) {
      for (int j=0; j<4; ++j) {
        myElts.push_back(i+10*j);
      }
    }
    break;

  case 1:
    //numMyElts=30;
    for (int i=40; i<50; ++i) {
      for (int j=0; j<3; ++j) {
        myElts.push_back(i+10*j);
      }
    }
    break;

  case 2:
    //numMyElts=30;
    for (int i=70; i<80; ++i) {
      for (int j=0; j<3; ++j) {
        myElts.push_back(i+10*j);
      }
    }
    break;
  }

/*
  My global indices: [0, 10, 20, 30, 1, 11, 21, 31, 2, 12, 22, 32, 3, 13, 23, 33,
                      4, 14, 24, 34, 5, 15, 25, 35, 6, 16, 26, 36, 7, 17, 27, 37,
                      8, 18, 28, 38, 9, 19, 29, 39]

 Process 1:
  My number of entries: 30
  My minimum global index: 40
  My maximum global index: 69
  My global indices: [40, 50, 60, 41, 51, 61, 42, 52, 62, 43, 53, 63, 44, 54, 64,
                      45, 55, 65, 46, 56, 66, 47, 57, 67, 48, 58, 68, 49, 59, 69]

 Process 2:
  My number of entries: 30
  My minimum global index: 70
  My maximum global index: 99
  My global indices: [70, 80, 90, 71, 81, 91, 72, 82, 92, 73, 83, 93, 74, 84, 94,
                      75, 85, 95, 76, 86, 96, 77, 87, 97, 78, 88, 98, 79, 89, 99]

*/

  const GO indexBase = 0;
  RCP<const MapType> map = rcp (new MapType (100, myElts (), indexBase, comm));

  RCP<FancyOStream> fos = Teuchos::fancyOStream (Teuchos::rcpFromRef (std::cout));
  //fos->setOutputToRootOnly(-1);
  //map->describe(*fos,Teuchos::VERB_EXTREME);

  const size_t numMyElements = map->getLocalNumElements ();
  switch (myRank) {
  case 0:
    assert(numMyElements==40);
    break;
  case 1:
    assert(numMyElements==30);
    break;
  case 2:
    assert(numMyElements==30);
    break;
  }

  // FIXME (mfh 19 Mar 2014) Once you set this
  // ("setOutputToRootOnly"), you can't unset it.
  fos->setOutputToRootOnly(0);
  *fos << std::endl << "Creating the sparse matrix" << std::endl;

  LO nnz=5;
  RCP<MatrixType> A = rcp (new MatrixType (map, nnz));

  ArrayView<const GO> myGlobalElements = map->getLocalElementList();

  GO center, left, right, lower, upper;
  std::vector<SC> vals(nnz);
  std::vector<GO> inds(nnz);

  //    e
  //  b a c
  //    d
  for (size_t i = 0; i < numMyElements; ++i)  {
    size_t n = 0;

    center = myGlobalElements[i] - indexBase;
    GetNeighboursCartesian2d<GO>(center, nx, ny, left, right, lower, upper);

    if (left  != -1) { inds[n] = left;  vals[n++] = b; }
    if (right != -1) { inds[n] = right; vals[n++] = c; }
    if (lower != -1) { inds[n] = lower; vals[n++] = d; }
    if (upper != -1) { inds[n] = upper; vals[n++] = e; }

    // diagonal
    SC z = a;
    inds[n]   = center;
    vals[n++] = z;

    for (size_t j = 0; j < n; j++)
      inds[j] += indexBase;

    ArrayView<GO> iv(&inds[0], n);
    ArrayView<SC> av(&vals[0], n);
    A->insertGlobalValues(myGlobalElements[i], iv, av);
  }

  A->fillComplete ();

  // Make sure that all processes finished fillComplete, before
  // reporting success.
  comm->barrier ();

}

#define UNIT_TEST_GROUP_SC_LO_GO_NO( SC, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, Bug6069_1, SC, LO, GO, NT)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(UNIT_TEST_GROUP_SC_LO_GO_NO)

} // namespace (anonymous)
