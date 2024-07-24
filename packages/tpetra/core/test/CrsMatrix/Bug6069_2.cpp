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
void
GetNeighboursCartesian2d (const GO i, const GO nx, const GO ny,
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

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, Bug6069_2, SC, LO, GO, NT)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Comm;
  using Teuchos::FancyOStream;
  using Tpetra::TestingUtilities::getDefaultComm;

  typedef Tpetra::CrsMatrix<SC, LO, GO, NT> MatrixType;
  typedef Tpetra::Map<LO, GO, NT> MapType;

  // global_size_t: Tpetra defines this unsigned integer type big
  // enough to hold any global dimension or amount of data.
  typedef Tpetra::global_size_t GST;

  RCP<const Comm<int> > comm = getDefaultComm();
  const size_t myRank = comm->getRank ();
  const size_t numProc = comm->getSize ();

  if (numProc != 2) {
    out << "This test must be run with exactly 2 MPI processes, but you ran "
        << "it with " << numProc << " process" << (numProc != 1 ? "es" : "")
        << "." << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                               "This test must be run with exactly "
                               "2 MPI processes, but you ran it with "
                               << numProc << " process"
                               << (numProc != 1 ? "es" : "") << ".");
  }

  // GO nx=10, ny=10;
  // SC a=4,b=-1,c=-1,d=-1,e=-1;

  Array<GO> myElts;

  switch (myRank) {
  case 0:
    myElts.push_back(0);
    myElts.push_back(1);
    break;
  case 1:
    myElts.push_back(2);
    myElts.push_back(3);
    break;
  }

  const GST GST_INV = Teuchos::OrdinalTraits<GST>::invalid ();
  const GO indexBase = 0;
  RCP<const MapType> map = rcp (new MapType (GST_INV, myElts (), indexBase, comm));

  RCP<FancyOStream> fos = Teuchos::fancyOStream (Teuchos::rcpFromRef (std::cout));
  fos->setOutputToRootOnly (-1);
  map->describe (*fos, Teuchos::VERB_EXTREME);

  const size_t numMyElements = map->getLocalNumElements ();
  switch (myRank) {
  case 0:
    assert(numMyElements==2);
    break;
  case 1:
    assert(numMyElements==2);
    break;
  }

  // FIXME (mfh 19 Mar 2014) Once you set this
  // ("setOutputToRootOnly"), you can't unset it.
  fos->setOutputToRootOnly (0);
  *fos << std::endl << "Creating the sparse matrix" << std::endl;

  const LO nnz = 4;
  RCP<MatrixType> A = rcp (new MatrixType (map, nnz));

  ArrayView<const GO> myGlobalElements =
    map->getLocalElementList ();

  // GO center, left, right, lower, upper;

  Teuchos::Array<SC> vals (nnz);
  Teuchos::Array<GO> inds (nnz);
  for (int i = 0; i < nnz; ++i) {
    inds[i] = i;
    vals[i] = 1.0;
  }

  comm->barrier ();
  for (int k = 0; k< comm->getSize (); ++k) {
    if (comm->getRank () == k) {
      std::cout << "pid " << k << " inserting global rows" << std::endl;
      for (size_t i = 0; i < numMyElements; ++i)  {
        //size_t n = 0;

        std::cout << "   grow " << myGlobalElements[i] << " : ";
        for (int jj = 0; jj < 4; ++jj) {
          std::cout << inds[jj] << " ";
        }
        std::cout << std::endl;
        A->insertGlobalValues (myGlobalElements[i], inds (), vals ());
      }
    }
    // mfh 01 Apr 2014: sleep() is a POSIX function, not a C++
    // standard function, and thus not suitable for Tpetra tests.
    //
    //sleep (1);
    comm->barrier ();
  } //k

  A->fillComplete ();

  // Make sure that all processes finished fillComplete, before
  // reporting success.
  comm->barrier ();
}

#define UNIT_TEST_GROUP_SC_LO_GO_NO( SC, LO, GO, NT )                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, Bug6069_2, SC, LO, GO, NT)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(UNIT_TEST_GROUP_SC_LO_GO_NO)

} // anonymous
