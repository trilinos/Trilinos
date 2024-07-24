// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "TpetraCore_config.h"
#ifndef HAVE_TPETRA_EPETRA
#  error "HAVE_TPETRA_EPETRA is undefined, but Tpetra is nevertheless building this file.  This probably means that Tpetra's CMake build system has a bug.  Please report this bug to the Tpetra developers."
#else

#include "Tpetra_TestingUtilities.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_EpetraRowMatrix.hpp"
#include "Tpetra_Details_getNumDiags.hpp"

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif


namespace {
  bool testMpi = true;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
  }

  Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm()
  {
    Teuchos::RCP<const Teuchos::Comm<int> > ret;
    if (testMpi) {
      ret = Tpetra::getDefaultComm();
    }
    else {
      ret = Teuchos::rcp(new Teuchos::SerialComm<int>());
    }
    return ret;
  }

  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( EpetraRowMatrix, BasicFunctionality, LO, GO )
  {
    using Teuchos::RCP;
    using Teuchos::tuple;
    //using std::endl;
    typedef Tpetra::global_size_t GST;
    // The Epetra wrapper only works for Scalar=double.
    typedef double Scalar;

    // generate a tridiagonal matrix
    typedef Teuchos::ScalarTraits<Scalar> ST;
    typedef Tpetra::CrsMatrix<Scalar,LO,GO> MAT;
    typedef Tpetra::MultiVector<Scalar,LO,GO> MV;
    typedef typename ST::magnitudeType Mag;
    typedef Teuchos::ScalarTraits<Mag> MT;
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid();
    // get a comm
    RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const int numImages = size(*comm);
    // create a Map
    const size_t numLocal = 10;
    const size_t numVecs = 5;
    RCP<const Tpetra::Map<LO,GO> > map = Tpetra::createContigMap<LO,GO>(INVALID,numLocal,comm);
    // create a matrix, modeled closely on Chris' CrsMatrix unit-tests.
    RCP<MAT> matrix(new MAT(map, 3));
    for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
      if (r == map->getMinAllGlobalIndex()) {
        matrix->insertGlobalValues(r, tuple(r,r+1), tuple(ST::one(),ST::one()) );
      }
      else if (r == map->getMaxAllGlobalIndex()) {
        matrix->insertGlobalValues(r, tuple(r-1,r), tuple(ST::one(),ST::one()) );
      }
      else {
        matrix->insertGlobalValues(r, tuple(r-1,r,r+1), tuple(ST::one(),ST::one(),ST::one()) );
      }
    }
    for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
      // increment the diagonals
      matrix->sumIntoGlobalValues(r, tuple(r), tuple(ST::one()) );
    }
    matrix->fillComplete();
    TEST_EQUALITY( Tpetra::Details::getLocalNumDiags (*matrix), static_cast<LO> (numLocal) );
    TEST_EQUALITY( Tpetra::Details::getGlobalNumDiags (*matrix), static_cast<GO> (numImages*numLocal) );
    TEST_EQUALITY( matrix->getGlobalNumEntries(), 3*numImages*numLocal - 2 );

#ifdef HAVE_MPI
    Epetra_MpiComm ecomm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm ecomm;
#endif

    //now create the EpetraRowMatrix class, which is the real target of this
    //unit-test:
    Tpetra::EpetraRowMatrix<MAT> erowmat(matrix, ecomm);

    int myRow = numLocal/2;
    int erowmat_NumMyRowEntries = -1;
    TEST_EQUALITY_CONST( erowmat.NumMyRowEntries( myRow, erowmat_NumMyRowEntries ) == 0, true );

    size_t numEntriesLocalRow = matrix->getNumEntriesInLocalRow(myRow);

    TEST_EQUALITY_CONST( erowmat_NumMyRowEntries == (int)numEntriesLocalRow, true );

    //test the matrix-vector product by comparing the result of
    //CrsMatrix::apply with the result of Epetra_RowMatrix::Multiply:

    MV tmv1(map,numVecs,true), tmv2(map,numVecs,true);
    tmv1.randomize();
    tmv2.randomize();
    tmv1.putScalar(1.0);
    tmv2.putScalar(0.0);
    matrix->apply(tmv1,tmv2);

    Epetra_BlockMap emap(int(numImages*numLocal), 1, 0, ecomm);
    Epetra_MultiVector emv1(emap,numVecs), emv2(emap,numVecs);
    emv1.PutScalar(1.0);
    emv2.PutScalar(0.0);
    erowmat.Multiply(false, emv1, emv2);

    Teuchos::ArrayView<Scalar> vals(emv2.Values(),numLocal*numVecs);
    MV tmvres(map,vals,numLocal,numVecs);
    tmvres.update(-ST::one(),tmv2,ST::one());
    Teuchos::Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    tmvres.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  //
  // INSTANTIATIONS
  //

TPETRA_ETI_MANGLING_TYPEDEFS()

#define UNIT_TEST_GROUP( LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( EpetraRowMatrix, BasicFunctionality, LO, GO )

TPETRA_INSTANTIATE_LG( UNIT_TEST_GROUP )

}

#endif // HAVE_TPETRA_EPETRA
