// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// Ifpack: Object-Oriented Algebraic Preconditioner Package
// Copyright (2002) NTESS
// *****************************************************************************
// @HEADER

/*! \file Ifpack2_UnitTestRelaxation.cpp

\brief Ifpack2 Unit test for the Relaxation template.
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_Relaxation.hpp>
#include <Tpetra_RowMatrix.hpp>
#include <Tpetra_BlockMultiVector.hpp>

namespace { // (anonmous)

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
typedef Tpetra::global_size_t GST;
typedef tif_utest::Node Node;

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> row_matrix_type;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

  //we are now in a class method declared by the above macro, and
  //that method has these input arguments:
  //Teuchos::FancyOStream& out, bool& success

  out << "Ifpack2::Version(): " << Ifpack2::Version () << std::endl;

  GST num_rows_per_proc = 5;

  RCP<const map_type > rowmap =
    tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node> (num_rows_per_proc);
  RCP<const crs_matrix_type> crsmatrix =
    tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> (rowmap);
  Ifpack2::Relaxation<row_matrix_type> prec (crsmatrix);

  Teuchos::ParameterList params;
  params.set("relaxation: type", "Jacobi");

  TEST_NOTHROW(prec.setParameters(params));

  // Insist that the preconditioner's domain and range Maps are the
  // same as those of the matrix.

  TEST_EQUALITY( crsmatrix->getDomainMap ()->isSameAs (* (prec.getDomainMap ())), true );
  TEST_EQUALITY( crsmatrix->getRangeMap ()->isSameAs (* (prec.getRangeMap ())), true );

  prec.initialize ();
  prec.compute ();

  MV x (rowmap, 2);
  MV y (rowmap, 2);
  x.putScalar (Teuchos::ScalarTraits<Scalar>::one ());

  prec.applyMat (x, y);

  {
    Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();

    //Since crsmatrix is a diagonal matrix with 2 on the diagonal,
    //y should be full of 2's now.

    Teuchos::ArrayRCP<Scalar> twos (num_rows_per_proc*2, 2);
    TEST_COMPARE_FLOATING_ARRAYS(yview, twos(), Teuchos::ScalarTraits<Scalar>::eps());
  }

  prec.apply(x, y);
  {
    Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();

    //y should be full of 0.5's now.
    Teuchos::ArrayRCP<Scalar> halfs(num_rows_per_proc*2, 0.5);

    TEST_COMPARE_FLOATING_ARRAYS(yview, halfs(), Teuchos::ScalarTraits<Scalar>::eps());
  }
}

// Test apply() with x == y.
// When x == y, apply() need to create internally an auxiliary vector, Xcopy.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, Test1, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> row_matrix_type;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

  out << "Ifpack2::Version(): " << Ifpack2::Version () << std::endl;

  GST num_rows_per_proc = 5;

  RCP<const map_type> rowmap =
    tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node> (num_rows_per_proc);
  RCP<const crs_matrix_type> crsmatrix =
    tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> (rowmap);
  Ifpack2::Relaxation<row_matrix_type> prec (crsmatrix);

  Teuchos::ParameterList params;
  params.set ("relaxation: type", "Jacobi");
  prec.setParameters (params);
  prec.initialize ();
  prec.compute ();

  MV x (rowmap, 2);
  x.putScalar (Teuchos::ScalarTraits<Scalar>::one ());

  prec.apply (x, x);

  //y should be full of 0.5's now.
  {
    Teuchos::ArrayRCP<const Scalar> xview = x.get1dView();
    Teuchos::ArrayRCP<Scalar> halfs(num_rows_per_proc*2, 0.5);
    TEST_COMPARE_FLOATING_ARRAYS(xview, halfs(), Teuchos::ScalarTraits<Scalar>::eps());
  }
}

// Test apply() with x != y  but x and y pointing to the same memory location.
// Tpetra Multivector public constructors are always copying input data so it is harder to reach such case than with Ifpack/Epetra.
// Nevertheless, it is still possible to create two different Tpetra vectors pointing to the same memory location by using MultiVector::subView().
// I can't imagine anyone trying to do this but... in this case, apply() need also to create internally an auxiliary vector, Xcopy.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, Test2, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> row_matrix_type;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  GST num_rows_per_proc = 5;

  RCP<const map_type> rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  RCP<const crs_matrix_type> crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);
  Ifpack2::Relaxation<row_matrix_type> prec (crsmatrix);

  Teuchos::ParameterList params;
  params.set("relaxation: type", "Jacobi");
  prec.setParameters(params);

  prec.initialize();
  prec.compute();

  MV y (rowmap, 2);
  y.putScalar (Teuchos::ScalarTraits<Scalar>::one ());
  RCP<const MV> xrcp = y.subView (Teuchos::Range1D (0,1));
  const MV& x = *xrcp;

  TEST_INEQUALITY(&x, &y); // vector x and y are different
  // Vectors x and y point to the same data.
  TEST_EQUALITY(x.getLocalViewHost (Tpetra::Access::ReadOnly).data (),
                y.getLocalViewHost (Tpetra::Access::ReadOnly).data ());

  prec.apply(x, y);

  //y should be full of 0.5's now.
  {
    Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
    Teuchos::ArrayRCP<Scalar> halfs(num_rows_per_proc*2, 0.5);
    TEST_COMPARE_FLOATING_ARRAYS(yview, halfs(), Teuchos::ScalarTraits<Scalar>::eps());
  }
}

  // Test apply() with a partially "null" x and y. In parallel, it is possible that some nodes do not have any local elements.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, Test3, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> row_matrix_type;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

  out << "Ifpack2::Version(): " << Ifpack2::Version() << std::endl;

  GST num_rows_per_proc = 0;
  auto comm = Tpetra::getDefaultComm();
  if(!comm->getRank()) num_rows_per_proc=5;

  RCP<const map_type> rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  RCP<const crs_matrix_type> crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);
  Ifpack2::Relaxation<row_matrix_type> prec (crsmatrix);

  Teuchos::ParameterList params;
  params.set("relaxation: type", "Jacobi");
  prec.setParameters(params);

  prec.initialize();
  prec.compute();

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,2), y(rowmap,2);
  x.putScalar (Teuchos::ScalarTraits<Scalar>::one ());

  TEST_EQUALITY(x.getMap()->getLocalNumElements(), num_rows_per_proc);
  TEST_EQUALITY(y.getMap()->getLocalNumElements(), num_rows_per_proc);

  TEST_NOTHROW(prec.apply(x, y));
}

  // Test apply() to make sure the L1 methods work
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, Test4, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> row_matrix_type;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

  out << "Ifpack2::Version(): " << Ifpack2::Version () << std::endl;

  GST num_rows_per_proc = 5;

  RCP<const map_type > rowmap =
    tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node> (num_rows_per_proc);
  RCP<const crs_matrix_type> crsmatrix =
    tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> (rowmap);
  Ifpack2::Relaxation<row_matrix_type> prec (crsmatrix);

  Teuchos::ParameterList params;
  params.set("relaxation: type", "Jacobi");
  params.set("relaxation: use l1",true);
  prec.setParameters(params);

  prec.initialize();
  prec.compute();

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,2), y(rowmap,2);
  x.putScalar (Teuchos::ScalarTraits<Scalar>::one ());

  TEST_EQUALITY(x.getMap()->getLocalNumElements(), 5);
  TEST_EQUALITY(y.getMap()->getLocalNumElements(), 5);
  TEST_NOTHROW(prec.apply(x, y));
}

  // Test apply() to make sure the Richardson methods work
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, Richardson, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> row_matrix_type;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

  out << "Ifpack2::Version(): " << Ifpack2::Version () << std::endl;

  GST num_rows_per_proc = 5;

  RCP<const map_type > rowmap =
    tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node> (num_rows_per_proc);
  RCP<const crs_matrix_type> crsmatrix =
    tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> (rowmap);
  Ifpack2::Relaxation<row_matrix_type> prec (crsmatrix);

  Teuchos::ParameterList params;
  params.set("relaxation: type", "Richardson");
  prec.setParameters(params);

  prec.initialize();
  prec.compute();

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,2), y(rowmap,2);
  x.putScalar (Teuchos::ScalarTraits<Scalar>::one ());

  TEST_EQUALITY(x.getMap()->getLocalNumElements(), 5);
  TEST_EQUALITY(y.getMap()->getLocalNumElements(), 5);
  TEST_NOTHROW(prec.apply(x, y));
}

// Check that (symmetric) Gauss-Seidel works if there are some MPI processes with zero rows of the matrix.
// Test contributed by Jonathan Hu on 04 Jan 2013.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, SymGaussSeidelZeroRows, Scalar, LocalOrdinal, GlobalOrdinal)
{
  using Teuchos::Comm;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using std::endl;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Tpetra::Map<LO, GO, Node> map_type;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar, LO, GO, Node> row_matrix_type;
  typedef Tpetra::MultiVector<Scalar, LO, GO, Node> mv_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  out << "Ifpack2::Version(): " << Ifpack2::Version () << endl;

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
  if (comm->getSize () == 1) {
    out << "The unit test's (MPI) communicator only contains one process."
        << endl << "This test only makes sense if the communicator contains "
        << "multiple processes." << endl << "I'll let the test pass trivially."
        << endl;
    return;
  }

  // The number of rows of the matrix and vectors owned by the calling
  // process.  One process (Proc 0) owns zero rows, and the other
  // process(es) own a nonzero amount of rows.  This is the point of
  // the test -- to ensure that (symmetric) Gauss-Seidel works if some
  // (but not all) processes have zero rows.
  LO nelPerProc;
  if (comm->getRank () == 0) {
    nelPerProc = 0;
  }
  else {
    nelPerProc = 100;
  }

  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
  RCP<const map_type> rowMap (new map_type (INVALID, nelPerProc, 0, comm));
  RCP<const crs_matrix_type> crsMatrix = tif_utest::create_test_matrix<Scalar,LO,GO,Node> (rowMap);

  // We don't really need prec to be a pointer here, but it's
  // convenient for putting the constructor call in a TEST_NOTHROW
  // macro.  Otherwise the declaration of prec would be in an inner
  // scope and we couldn't use it below.
  RCP<Ifpack2::Relaxation<row_matrix_type> > prec;
  {
    TEST_NOTHROW( prec = rcp (new Ifpack2::Relaxation<row_matrix_type> (crsMatrix)) );
  }

  ParameterList params;
  params.set ("relaxation: type", "Symmetric Gauss-Seidel");
  prec->setParameters (params);
  TEST_NOTHROW( prec->initialize () );
  TEST_NOTHROW( prec->compute () );

  mv_type X (rowMap, 2);
  X.putScalar (STS::zero ());
  mv_type B (rowMap, 2);
  B.putScalar (STS::one ());

  prec->apply (B, X);
}



// Check that local (symmetric) Gauss-Seidel works if there are some MPI processes with zero rows of the matrix.
// Test contributed by Jonathan Hu on 04 Jan 2013.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, LocalSymGaussSeidelZeroRows, Scalar, LocalOrdinal, GlobalOrdinal)
{
  using Teuchos::Comm;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using std::endl;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Tpetra::Map<LO, GO, Node> map_type;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar, LO, GO, Node> row_matrix_type;
  typedef Tpetra::MultiVector<Scalar, LO, GO, Node> mv_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << endl;

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
  if (comm->getSize () == 1) {
    out << "The unit test's (MPI) communicator only contains one process."
        << endl << "This test only makes sense if the communicator contains "
        << "multiple processes." << endl << "I'll let the test pass trivially."
        << endl;
    return;
  }

  // The number of rows of the matrix and vectors owned by the calling
  // process.  One process (Proc 0) owns zero rows, and the other
  // process(es) own a nonzero amount of rows.  This is the point of
  // the test -- to ensure that (symmetric) Gauss-Seidel works if some
  // (but not all) processes have zero rows.
  LO nelPerProc;
  if (comm->getRank () == 0) {
    nelPerProc = 0;
  } else {
    nelPerProc = 100;
  }

  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
  RCP<const map_type> rowMap (new map_type (INVALID, nelPerProc, 0, comm));
  RCP<const crs_matrix_type> crsMatrix = tif_utest::create_test_matrix<Scalar,LO,GO,Node> (rowMap);

  // We don't really need prec to be a pointer here, but it's
  // convenient for putting the constructor call in a TEST_NOTHROW
  // macro.  Otherwise the declaration of prec would be in an inner
  // scope and we couldn't use it below.
  RCP<Ifpack2::Relaxation<row_matrix_type> > prec;
  {
    TEST_NOTHROW( prec = rcp (new Ifpack2::Relaxation<row_matrix_type> (crsMatrix)) );
  }

  ParameterList params;
  params.set ("relaxation: type", "Symmetric Gauss-Seidel");

  Teuchos::ArrayRCP<LO> rowIndices(crsMatrix->getLocalNumRows());
  for (size_t i = 0; i < crsMatrix->getLocalNumRows (); ++i) {
    rowIndices[i] = static_cast<LO> (i);
  }
  params.set ("relaxation: local smoothing indices", rowIndices);

  prec->setParameters (params);
  TEST_NOTHROW( prec->initialize () );
  TEST_NOTHROW( prec->compute () );

  mv_type X (rowMap, 2);
  X.putScalar (STS::zero ());
  mv_type B (rowMap, 2);
  B.putScalar (STS::one ());

  prec->apply (B, X);
}


// Check that multiple sweeps of Symmetric Gauss-Seidel (SGS) work correctly.
//
// The point is really to ensure that SGS does the Import (from the
// domain Map input X to the column Map version of X) in the right
// place, if an Import is needed.  Thus, this example uses a matrix
// with a nontrivial Import (i.e., where the domain and column Maps
// are not the same).
//md - 27 Jan 2016: The bug is fixed in which diagonals were set to 0.
//The test fails after the bug fix, which is logical, as we don't expect the
//exact same solution for distributed SGS and sequential SGS. I disabled the
//test with Siva's suggestion. (Mehmet)
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, SGS_mult_sweeps, Scalar, LocalOrdinal, GlobalOrdinal)
{
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::Comm;
  using Teuchos::ParameterList;
  using std::endl;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> crs_matrix_type;
  typedef Tpetra::Import<LO, GO, Node> import_type;
  typedef Tpetra::Map<LO, GO, Node> map_type;
  typedef Tpetra::Vector<Scalar, LO, GO, Node> vec_type;
  typedef Tpetra::RowMatrix<Scalar, LO, GO, Node> row_matrix_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef Teuchos::ScalarTraits<typename STS::magnitudeType> STM;

  out << "Test multiple Symmetric Gauss-Seidel sweeps with nontrivial Import"
      << endl;
  Teuchos::OSTab tab0 (out);

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  if (numProcs == 1) {
    out << "The unit test's (MPI) communicator only contains one process."
        << endl << "This test only makes sense if the communicator contains "
        << "multiple processes." << endl << "I'll let the test pass trivially."
        << endl;
    return;
  }

  const GST gblNumRows = static_cast<GST> (numProcs * 2);
  const GO indexBase = 0;
  RCP<const map_type> rowMap (new map_type (gblNumRows, indexBase, comm));
  RCP<const map_type> domainMap = rowMap;
  RCP<const map_type> rangeMap = rowMap;
  RCP<crs_matrix_type> A =
    rcp (new crs_matrix_type (rowMap, 6));

  {
    const size_t lclNumRows = rowMap->getLocalNumElements ();
    const Scalar ONE = STS::one ();
    const Scalar TWO = ONE + ONE;
    const Scalar FOUR = TWO + TWO;
    const Scalar EIGHT = FOUR + FOUR;
    const Scalar TWELVE = EIGHT + FOUR;
    Array<Scalar> vals (6);
    Array<GO> gblColInds (6);

    for (LO lclRow = 0; lclRow < static_cast<LO> (lclNumRows); ++lclRow) {
      const GO gblRow = rowMap->getGlobalElement (lclRow);
      const Scalar diagVal = TWELVE;
      size_t numEnt = 6;
      if (gblRow == rowMap->getMinAllGlobalIndex ()) {
        numEnt = 4;
        vals[0] = diagVal;
        vals[1] = ONE / TWO;
        vals[2] = -ONE;
        vals[3] = ONE;
        gblColInds[0] = rowMap->getMinAllGlobalIndex ();
        gblColInds[1] = gblColInds[0] + 1;
        gblColInds[2] = gblColInds[0] + 2;
        gblColInds[3] = gblColInds[0] + 3;
      }
      else if (gblRow == rowMap->getMinAllGlobalIndex () + 1) {
        numEnt = 4;
        vals[0] = -ONE / TWO;
        vals[1] = diagVal / TWO;
        vals[2] = ONE;
        vals[3] = -ONE;
        gblColInds[0] = rowMap->getMinAllGlobalIndex ();
        gblColInds[1] = gblColInds[0] + 1;
        gblColInds[2] = gblColInds[0] + 2;
        gblColInds[3] = gblColInds[0] + 3;
      }
      else if (gblRow == rowMap->getMaxAllGlobalIndex () - 1) {
        numEnt = 4;
        vals[0] = -ONE;
        vals[1] = ONE;
        vals[2] = diagVal;
        vals[3] = ONE / TWO;
        gblColInds[0] = rowMap->getMaxAllGlobalIndex () - 3;
        gblColInds[1] = rowMap->getMaxAllGlobalIndex () - 2;
        gblColInds[2] = rowMap->getMaxAllGlobalIndex () - 1;
        gblColInds[3] = rowMap->getMaxAllGlobalIndex ();
      }
      else if (gblRow == rowMap->getMaxAllGlobalIndex ()) {
        numEnt = 4;
        vals[0] = ONE;
        vals[1] = -ONE;
        vals[2] = -ONE / TWO;
        vals[3] = diagVal / TWO;
        gblColInds[0] = rowMap->getMaxAllGlobalIndex () - 3;
        gblColInds[1] = rowMap->getMaxAllGlobalIndex () - 2;
        gblColInds[2] = rowMap->getMaxAllGlobalIndex () - 1;
        gblColInds[3] = rowMap->getMaxAllGlobalIndex ();
      }
      else if (gblRow % 2 == static_cast<GO> (0)) {
        numEnt = 6;
        vals[0] = -ONE;
        vals[1] = ONE;
        vals[2] = diagVal;
        vals[3] = ONE / TWO;
        vals[4] = -ONE;
        vals[5] = ONE;
        gblColInds[0] = gblRow - 2;
        gblColInds[1] = gblRow - 1;
        gblColInds[2] = gblRow;
        gblColInds[3] = gblRow + 1;
        gblColInds[4] = gblRow + 2;
        gblColInds[5] = gblRow + 3;
      }
      else { // gblRow % 2 != 0
        numEnt = 6;
        vals[0] = ONE;
        vals[1] = -ONE;
        vals[2] = -ONE / TWO;
        vals[3] = diagVal;
        vals[4] = ONE;
        vals[5] = -ONE;
        gblColInds[0] = gblRow - 3;
        gblColInds[1] = gblRow - 2;
        gblColInds[2] = gblRow - 1;
        gblColInds[3] = gblRow;
        gblColInds[4] = gblRow + 1;
        gblColInds[5] = gblRow + 2;
      }

      ArrayView<Scalar> valsView = vals (0, numEnt);
      ArrayView<GO> gblColIndsView = gblColInds (0, numEnt);
      A->insertGlobalValues (gblRow, gblColIndsView, valsView);
    }
  }
  A->fillComplete (domainMap, rangeMap);

  RCP<const map_type> gatherRowMap;
  {
    const size_t lclNumRows =
      (myRank == 0) ? static_cast<size_t> (gblNumRows) : size_t (0);
    gatherRowMap = rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));
  }

  RCP<const map_type> gatherDomainMap = gatherRowMap;
  RCP<const map_type> gatherRangeMap = gatherRowMap;

  RCP<crs_matrix_type> A_gather = rcp (new crs_matrix_type (gatherRowMap, 6));
  import_type import (rowMap, gatherRowMap);
  A_gather->doImport (*A, import, Tpetra::INSERT);
  A_gather->fillComplete (gatherDomainMap, gatherRangeMap);

  vec_type X (domainMap);
  vec_type Y (rangeMap);
  vec_type X_gather (gatherDomainMap);
  vec_type Y_gather (gatherRangeMap);
  vec_type Y_diff (gatherRangeMap);

  // Test Symmetric Gauss-Seidel (SGS) with three sweeps.
  // Start by letting SGS set the starting solution to zero.
  ParameterList params;
  params.set ("relaxation: type", "Symmetric Gauss-Seidel");
  params.set ("relaxation: sweeps", 3);
  params.set ("relaxation: zero starting solution", true);

  Ifpack2::Relaxation<row_matrix_type> prec (A);
  prec.setParameters (params);
  TEST_NOTHROW( prec.initialize () );
  TEST_NOTHROW( prec.compute () );

  Ifpack2::Relaxation<row_matrix_type> gatherPrec (A_gather);
  gatherPrec.setParameters (params);
  TEST_NOTHROW( gatherPrec.initialize () );
  TEST_NOTHROW( gatherPrec.compute () );

  X.randomize ();
  X_gather.doImport (X, import, Tpetra::REPLACE);
  Y.randomize ();
  X_gather.doImport (X, import, Tpetra::REPLACE);

  prec.apply (X, Y);
  gatherPrec.apply (X_gather, Y_gather);
  Y_diff.doImport (Y, import, Tpetra::REPLACE);
  Y_diff.update (STS::one (), Y_gather, -STS::one ());

  typename STS::magnitudeType normInf = Y_diff.normInf ();
  TEST_EQUALITY(normInf, STM::zero ());

  out << "Repeat test without setting starting solution to zero" << endl;

  // Repeat the test without setting the starting solution to zero.
  params.set ("relaxation: zero starting solution", false);

  prec.setParameters (params);
  TEST_NOTHROW( prec.initialize () );
  TEST_NOTHROW( prec.compute () );

  gatherPrec.setParameters (params);
  TEST_NOTHROW( gatherPrec.initialize () );
  TEST_NOTHROW( gatherPrec.compute () );

  X.randomize ();
  X_gather.doImport (X, import, Tpetra::REPLACE);
  Y.randomize ();
  X_gather.doImport (X, import, Tpetra::REPLACE);

  prec.apply (X, Y);
  gatherPrec.apply (X_gather, Y_gather);

  Y_diff.doImport (Y, import, Tpetra::REPLACE);
  Y_diff.update (STS::one (), Y_gather, -STS::one ());
  normInf = Y_diff.normInf ();
  TEST_EQUALITY( normInf, STM::zero () );
}


// Macro used inside the unit test below.  It tests for global error,
// and if so, prints each process' error message and quits the test
// early.
//
// 'out' only prints on Process 0.  It's really not OK for other
// processes to print to stdout, but it usually works and we need to
// do it for debugging.
#define IFPACK2RELAXATION_REPORT_GLOBAL_ERR( WHAT_STRING ) do { \
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess)); \
  TEST_EQUALITY_CONST( gblSuccess, 1 ); \
  if (gblSuccess != 1) { \
    out << WHAT_STRING << " FAILED on one or more processes!" << endl; \
    for (int p = 0; p < numProcs; ++p) { \
      if (myRank == p && lclSuccess != 1) { \
        std::cout << errStrm.str () << std::flush; \
      } \
      comm->barrier (); \
      comm->barrier (); \
      comm->barrier (); \
    } \
    std::cerr << "TEST FAILED; RETURNING EARLY" << endl; \
    return; \
  } \
} while (false)


// Test apply() on a NotCrsMatrix, where some processes own zero rows
// of the domain and range Map vectors.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, NotCrsMatrix, Scalar, LO, GO)
{
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::endl;
  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node> row_matrix_type;
  typedef tif_utest::NotCrsMatrix<Scalar,LO,GO,Node> not_crs_matrix_type;
  typedef Tpetra::Map<LO,GO,Node> map_type;
  typedef Ifpack2::Relaxation<row_matrix_type> prec_type;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
  int lclSuccess = 1;
  int gblSuccess = 1;
  std::ostringstream errStrm; // for error collection

  out << "Ifpack2 \"NonCrsMatrix\" test" << endl;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  const GST num_rows_per_proc = (myRank == 0) ? 5 : 0;
  RCP<const map_type> rowmap;
  try {
    rowmap = tif_utest::create_tpetra_map<LO, GO, Node> (num_rows_per_proc);
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": create_tpetra_map threw exception: "
            << e.what () << endl;
  }
  IFPACK2RELAXATION_REPORT_GLOBAL_ERR( "create_tpetra_map" );

  RCP<crs_matrix_type> crsmatrix;
  try {
    crsmatrix = rcp_const_cast<crs_matrix_type> (tif_utest::create_test_matrix<Scalar, LO, GO, Node> (rowmap));
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": create_test_matrix threw exception: "
            << e.what () << endl;
  }
  IFPACK2RELAXATION_REPORT_GLOBAL_ERR( "create_test_matrix" );

  RCP<not_crs_matrix_type> notcrsmatrix;
  try {
    notcrsmatrix = rcp (new not_crs_matrix_type (crsmatrix));
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": NotCrsMatrix constructor threw exception: "
            << e.what () << endl;
  }
  IFPACK2RELAXATION_REPORT_GLOBAL_ERR( "NotCrsMatrix constructor" );

  RCP<prec_type> prec;
  try {
    prec = rcp (new prec_type (notcrsmatrix));
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": Preconditioner constructor threw exception: "
            << e.what () << endl;
  }
  IFPACK2RELAXATION_REPORT_GLOBAL_ERR( "Preconditioner constructor" );

  Teuchos::ParameterList params;
  params.set ("relaxation: type", "Jacobi");
  try {
    prec->setParameters (params);
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": prec->setParameters() threw exception: "
            << e.what () << endl;
  }
  IFPACK2RELAXATION_REPORT_GLOBAL_ERR( "prec->setParameters()" );

  try {
    prec->initialize ();
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": prec->initialize() threw exception: "
            << e.what () << endl;
  }
  IFPACK2RELAXATION_REPORT_GLOBAL_ERR( "prec->initialize()" );

  try {
    prec->compute ();
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": prec->compute() threw exception: "
            << e.what () << endl;
  }
  IFPACK2RELAXATION_REPORT_GLOBAL_ERR( "prec->compute()" );

  MV x (rowmap, 2);
  MV y (rowmap, 2);
  x.putScalar (Teuchos::ScalarTraits<Scalar>::one ());

  TEST_EQUALITY( x.getMap()->getLocalNumElements(), num_rows_per_proc );
  TEST_EQUALITY( y.getMap()->getLocalNumElements(), num_rows_per_proc );

  try {
    prec->apply (x, y);
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": prec->apply(x, y) threw exception: "
            << e.what () << endl;
  }
  IFPACK2RELAXATION_REPORT_GLOBAL_ERR( "prec->apply(x, y)" );
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, TestDiagonalBlockCrsMatrix, Scalar, LO, GO)
{
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::endl;
  typedef Tpetra::BlockCrsMatrix<Scalar,LO,GO,Node> block_crs_matrix_type;
  typedef Tpetra::BlockMultiVector<Scalar,LO,GO,Node> BMV;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node> row_matrix_type;
  typedef Tpetra::CrsGraph<LO,GO,Node> crs_graph_type;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
  typedef Ifpack2::Relaxation<row_matrix_type> prec_type;
  int lclSuccess = 1;
  int gblSuccess = 1;
  std::ostringstream errStrm; // for error collection

  out << "Ifpack2::Relaxation diagonal block matrix test" << endl;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  const int num_rows_per_proc = 5;
  const int blockSize = 3;
  RCP<crs_graph_type> crsgraph =
    tif_utest::create_diagonal_graph<LO,GO,Node> (num_rows_per_proc);

  RCP<block_crs_matrix_type> bcrsmatrix;
  bcrsmatrix = rcp_const_cast<block_crs_matrix_type> (tif_utest::create_block_diagonal_matrix<Scalar,LO,GO,Node> (crsgraph, blockSize));

  RCP<prec_type> prec;
  try {
    prec = rcp (new prec_type (bcrsmatrix));
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": Preconditioner constructor threw exception: "
            << e.what () << endl;
  }
  IFPACK2RELAXATION_REPORT_GLOBAL_ERR( "Preconditioner constructor" );

  Teuchos::ParameterList params;
  params.set ("relaxation: type", "Jacobi");
  try {
    prec->setParameters (params);
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": prec->setParameters() threw exception: "
            << e.what () << endl;
  }
  IFPACK2RELAXATION_REPORT_GLOBAL_ERR( "prec->setParameters()" );

  try {
    prec->initialize ();
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": prec->initialize() threw exception: "
            << e.what () << endl;
  }
  IFPACK2RELAXATION_REPORT_GLOBAL_ERR( "prec->initialize()" );

  try {
    prec->compute ();
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": prec->compute() threw exception: "
            << e.what () << endl;
  }
  IFPACK2RELAXATION_REPORT_GLOBAL_ERR( "prec->compute()" );

  BMV xBlock (*crsgraph->getRowMap (), blockSize, 1);
  BMV yBlock (*crsgraph->getRowMap (), blockSize, 1);
  MV x = xBlock.getMultiVectorView ();
  MV y = yBlock.getMultiVectorView ();
  x.putScalar (Teuchos::ScalarTraits<Scalar>::one ());

  TEST_EQUALITY(x.getMap()->getLocalNumElements(), blockSize*num_rows_per_proc);
  TEST_EQUALITY(y.getMap()->getLocalNumElements(), blockSize*num_rows_per_proc);

  try {
    prec->apply (x, y);
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": prec->apply(x, y) threw exception: "
            << e.what () << endl;
  }
  IFPACK2RELAXATION_REPORT_GLOBAL_ERR( "prec->apply(x, y)" );

  const Scalar exactSol = 0.2;
  using STS = Teuchos::ScalarTraits<Scalar>;
  using mag_type = typename STS::magnitudeType;
  const auto tol = mag_type(100.0) * STS::eps();

  for (int k = 0; k < num_rows_per_proc; ++k) {
    auto ylcl = yBlock.getLocalBlockHost(k, 0, Tpetra::Access::ReadOnly);
    for (int j = 0; j < blockSize; ++j) {
      TEST_FLOATING_EQUALITY(ylcl(j), exactSol, tol);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, TestLowerTriangularBlockCrsMatrix, Scalar, LO, GO)
{
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::endl;
  typedef Tpetra::BlockCrsMatrix<Scalar,LO,GO,Node> block_crs_matrix_type;
  typedef Tpetra::CrsGraph<LO,GO,Node> crs_graph_type;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node> row_matrix_type;
  typedef Tpetra::BlockMultiVector<Scalar,LO,GO,Node> BMV;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
  typedef Ifpack2::Relaxation<row_matrix_type> prec_type;
  int lclSuccess = 1;
  int gblSuccess = 1;
  std::ostringstream errStrm; // for error collection

  out << "Ifpack2::Relaxation lower triangular BlockCrsMatrix test" << endl;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  const size_t num_rows_per_proc = 3;
  RCP<crs_graph_type> crsgraph;
  try {
    crsgraph = tif_utest::create_dense_local_graph<LO, GO, Node> (num_rows_per_proc);
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": create_dense_local_graph threw exception: "
            << e.what () << endl;
  }
  IFPACK2RELAXATION_REPORT_GLOBAL_ERR( "create_dense_local_graph" );

  const int blockSize = 5;
  RCP<block_crs_matrix_type> bcrsmatrix;
  try {
    bcrsmatrix = rcp_const_cast<block_crs_matrix_type> (tif_utest::create_triangular_matrix<Scalar, LO, GO, Node, true> (crsgraph, blockSize));
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": create_triangular_matrix threw exception: "
            << e.what () << endl;
  }
  IFPACK2RELAXATION_REPORT_GLOBAL_ERR( "create_triangular_matrix" );

  RCP<prec_type> prec;
  try {
    prec = rcp (new prec_type (bcrsmatrix));
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": Preconditioner constructor threw exception: "
            << e.what () << endl;
  }
  IFPACK2RELAXATION_REPORT_GLOBAL_ERR( "Preconditioner constructor" );
  //Issue #9400: check that prec isn't holding any views to the matrix
  bcrsmatrix->getValuesHost();
  bcrsmatrix->getValuesDevice();
  bcrsmatrix->getCrsGraph().getLocalGraphDevice();
  bcrsmatrix->getCrsGraph().getLocalGraphHost();

  Teuchos::ParameterList params;
  params.set ("relaxation: type", "Gauss-Seidel");
  try {
    prec->setParameters (params);
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": prec->setParameters() threw exception: "
            << e.what () << endl;
  }
  IFPACK2RELAXATION_REPORT_GLOBAL_ERR( "prec->setParameters()" );

  try {
    prec->initialize ();
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": prec->initialize() threw exception: "
            << e.what () << endl;
  }
  IFPACK2RELAXATION_REPORT_GLOBAL_ERR( "prec->initialize()" );

  bcrsmatrix->getValuesHost();
  bcrsmatrix->getValuesDevice();
  bcrsmatrix->getCrsGraph().getLocalGraphDevice();
  bcrsmatrix->getCrsGraph().getLocalGraphHost();

  try {
    prec->compute ();
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": prec->compute() threw exception: "
            << e.what () << endl;
  }
  IFPACK2RELAXATION_REPORT_GLOBAL_ERR( "prec->compute()" );

  bcrsmatrix->getValuesHost();
  bcrsmatrix->getValuesDevice();
  bcrsmatrix->getCrsGraph().getLocalGraphDevice();
  bcrsmatrix->getCrsGraph().getLocalGraphHost();

  BMV xBlock (* (crsgraph->getRowMap ()), blockSize, 1);
  BMV yBlock (* (crsgraph->getRowMap ()), blockSize, 1);
  MV x = xBlock.getMultiVectorView ();
  MV y = yBlock.getMultiVectorView ();
  x.putScalar (Teuchos::ScalarTraits<Scalar>::one ());

  TEST_EQUALITY( x.getMap()->getLocalNumElements (), blockSize * num_rows_per_proc );
  TEST_EQUALITY( y.getMap ()->getLocalNumElements (), blockSize * num_rows_per_proc );

  try {
    prec->apply (x, y);
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": prec->apply(x, y) threw exception: "
            << e.what () << endl;
  }
  IFPACK2RELAXATION_REPORT_GLOBAL_ERR( "prec->apply(x, y)" );

  Teuchos::Array<Scalar> exactSol(num_rows_per_proc);
  exactSol[0] = 0.5;
  exactSol[1] = -0.25;
  exactSol[2] = 0.625;

  for (size_t k = 0; k < num_rows_per_proc; ++k) {
    LO lcl_row = k;
    auto ylcl = yBlock.getLocalBlockHost(lcl_row, 0, Tpetra::Access::ReadOnly);
    for (int j = 0; j < blockSize; ++j) {
      TEST_FLOATING_EQUALITY(ylcl(j), exactSol[k], 1e-14);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, TestUpperTriangularBlockCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> block_crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> row_matrix_type;
  typedef Tpetra::BlockMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> BMV;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

  out << "Ifpack2::Version(): " << Ifpack2::Version () << std::endl;

  const int num_rows_per_proc = 3;
  const int blockSize = 5;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph =
    tif_utest::create_dense_local_graph<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  RCP<block_crs_matrix_type> bcrsmatrix =
    rcp_const_cast<block_crs_matrix_type> (tif_utest::create_triangular_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,false> (crsgraph, blockSize));

  Ifpack2::Relaxation<row_matrix_type> prec (bcrsmatrix);

  Teuchos::ParameterList params;
  params.set("relaxation: type", "Symmetric Gauss-Seidel");
  prec.setParameters(params);

  prec.initialize();
  TEST_NOTHROW(prec.compute());

  BMV xBlock (*crsgraph->getRowMap (), blockSize, 1);
  BMV yBlock (*crsgraph->getRowMap (), blockSize, 1);
  MV x = xBlock.getMultiVectorView ();
  MV y = yBlock.getMultiVectorView ();
  x.putScalar (Teuchos::ScalarTraits<Scalar>::one ());

  TEST_EQUALITY(x.getMap()->getLocalNumElements(), blockSize*num_rows_per_proc);
  TEST_EQUALITY(y.getMap()->getLocalNumElements(), blockSize*num_rows_per_proc);
  TEST_NOTHROW(prec.apply(x, y));

  Teuchos::Array<Scalar> exactSol(num_rows_per_proc);
  exactSol[0] = 0.625;
  exactSol[1] = -0.25;
  exactSol[2] = 0.5;

  for (int k = 0; k < num_rows_per_proc; ++k) {
    auto ylcl = yBlock.getLocalBlockHost(k, 0, Tpetra::Access::ReadOnly);
    for (int j = 0; j < blockSize; ++j) {
      TEST_FLOATING_EQUALITY(ylcl(j), exactSol[k], 1e-14);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, GS_Crs, Scalar, LocalOrdinal, GlobalOrdinal)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using crs_matrix_type = Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  using row_matrix_type = Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  using MV = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  using map_type = Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>;
  using prec_type = Ifpack2::Relaxation<row_matrix_type>;
  using STS = Teuchos::ScalarTraits<Scalar>;
  using STM = typename STS::magnitudeType;
  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;
  //Generate banded test matrix
  RCP<const map_type> rowmap = tif_utest::create_tpetra_map<LocalOrdinal, GlobalOrdinal, Node>(100);
  RCP<const crs_matrix_type> A = tif_utest::create_banded_matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowmap, 8);
  RCP<prec_type> prec = rcp(new prec_type(A));
  //Issue #9400: make sure prec isn't holding any views for A's local matrix
  A->getLocalMatrixDevice();
  A->getLocalMatrixHost();
  ParameterList params;
  params.set("relaxation: type", "Gauss-Seidel");
  params.set("relaxation: sweeps", 3);
  prec->setParameters (params);
  prec->initialize();
  A->getLocalMatrixDevice();
  A->getLocalMatrixHost();
  prec->compute();
  A->getLocalMatrixDevice();
  A->getLocalMatrixHost();
  //Set up linear problem
  const int numVecs = 10;
  MV x(A->getDomainMap(), numVecs, true);
  MV b(rowmap, numVecs, false);
  b.randomize();
  Kokkos::View<STM*, Kokkos::HostSpace> initNorms("Initial norms", numVecs);
  //Residual norms for starting solution of zero
  b.norm2(initNorms);
  prec->apply(b, x);
  //Compute residual vector = b - Ax
  MV residual(b, Teuchos::Copy);
  A->apply(x, residual, Teuchos::NO_TRANS, -STS::one(), STS::one());
  Kokkos::View<STM*, Kokkos::HostSpace> resNorms("Residual norms", numVecs);
  residual.norm2(resNorms);
  //Make sure all residual norms are significantly smaller than initial
  for(int i = 0; i < numVecs; i++)
  {
    TEST_COMPARE(resNorms(i), <, 0.5 * initNorms(i));
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, MTSGS, Scalar, LocalOrdinal, GlobalOrdinal)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using crs_matrix_type = Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  using row_matrix_type = Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  using MV = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  using map_type = Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>;
  using prec_type = Ifpack2::Relaxation<row_matrix_type>;
  using STS = Teuchos::ScalarTraits<Scalar>;
  using STM = typename STS::magnitudeType;
  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;
  //Generate banded test matrix
  RCP<const map_type> rowmap = tif_utest::create_tpetra_map<LocalOrdinal, GlobalOrdinal, Node>(100);
  RCP<const crs_matrix_type> A = tif_utest::create_banded_matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowmap, 8);
  RCP<prec_type> prec = rcp(new prec_type(A));
  ParameterList params;
  params.set("relaxation: type", "MT Symmetric Gauss-Seidel");
  params.set("relaxation: sweeps", 3);
  prec->setParameters (params);
  prec->initialize();
  prec->compute();
  //Set up linear problem
  const int numVecs = 10;
  MV x(A->getDomainMap(), numVecs, true);
  MV b(rowmap, numVecs, false);
  b.randomize();
  Kokkos::View<STM*, Kokkos::HostSpace> initNorms("Initial norms", numVecs);
  //Residual norms for starting solution of zero
  b.norm2(initNorms);
  prec->apply(b, x);
  //Compute residual vector = b - Ax
  MV residual(b, Teuchos::Copy);
  A->apply(x, residual, Teuchos::NO_TRANS, -STS::one(), STS::one());
  Kokkos::View<STM*, Kokkos::HostSpace> resNorms("Residual norms", numVecs);
  residual.norm2(resNorms);
  //Make sure all residual norms are significantly smaller than initial
  for(int i = 0; i < numVecs; i++)
  {
    TEST_COMPARE(resNorms(i), <, 0.5 * initNorms(i));
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, MTSGS_LongRows, Scalar, LocalOrdinal, GlobalOrdinal)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using crs_matrix_type = Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  using row_matrix_type = Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  using MV = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  using map_type = Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>;
  using prec_type = Ifpack2::Relaxation<row_matrix_type>;
  using STS = Teuchos::ScalarTraits<Scalar>;
  using STM = typename STS::magnitudeType;
  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;
  //Generate banded test matrix
  RCP<const map_type> rowmap = tif_utest::create_tpetra_map<LocalOrdinal, GlobalOrdinal, Node>(100);
  RCP<const crs_matrix_type> A = tif_utest::create_banded_matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowmap, 3);
  RCP<prec_type> prec = rcp(new prec_type(A));
  ParameterList goodParams;
  goodParams.set("relaxation: type", "MT Symmetric Gauss-Seidel");
  goodParams.set("relaxation: sweeps", 3);
  goodParams.set("relaxation: long row threshold", 3);
  //Try setting up precondition with incompatible type, and make sure this throws.
  {
    ParameterList badParams = goodParams;
    badParams.set("relaxation: type", "Gauss-Seidel");
    TEST_THROW (prec->setParameters (badParams), std::invalid_argument);
  }
  //Try setting up cluster GS preconditioner with long row algorithm enabled - should also throw.
  {
    ParameterList badParams = goodParams;
    badParams.set("relaxation: mtgs cluster size", 4);
    TEST_THROW(prec->setParameters (badParams), std::invalid_argument);
  }
  prec->setParameters (goodParams);
  prec->initialize();
  prec->compute();
  //Set up linear problem
  const int numVecs = 10;
  MV x(A->getDomainMap(), numVecs, true);
  MV b(rowmap, numVecs, false);
  b.randomize();
  Kokkos::View<STM*, Kokkos::HostSpace> initNorms("Initial norms", numVecs);
  //Residual norms for starting solution of zero
  b.norm2(initNorms);
  prec->apply(b, x);
  //Compute residual vector = b - Ax
  MV residual(b, Teuchos::Copy);
  A->apply(x, residual, Teuchos::NO_TRANS, -STS::one(), STS::one());
  Kokkos::View<STM*, Kokkos::HostSpace> resNorms("Residual norms", numVecs);
  residual.norm2(resNorms);
  //Make sure all residual norms are significantly smaller than initial
  for(int i = 0; i < numVecs; i++)
  {
    TEST_COMPARE(resNorms(i), <, 0.5 * initNorms(i));
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Relaxation, ClusterMTSGS, Scalar, LocalOrdinal, GlobalOrdinal)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using crs_matrix_type = Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  using row_matrix_type = Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  using MV = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  using map_type = Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>;
  using prec_type = Ifpack2::Relaxation<row_matrix_type>;
  using STS = Teuchos::ScalarTraits<Scalar>;
  using STM = typename STS::magnitudeType;
  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;
  //Generate banded test matrix
  RCP<const map_type> rowmap = tif_utest::create_tpetra_map<LocalOrdinal, GlobalOrdinal, Node>(100);
  RCP<const crs_matrix_type> A = tif_utest::create_banded_matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowmap, 8);
  RCP<prec_type> prec = rcp(new prec_type(A));
  ParameterList params;
  params.set("relaxation: type", "MT Symmetric Gauss-Seidel");
  params.set("relaxation: sweeps", 3);
  params.set("relaxation: mtgs cluster size", 4);
  prec->setParameters (params);
  prec->initialize();
  prec->compute();
  //Set up linear problem
  const int numVecs = 10;
  MV x(A->getDomainMap(), numVecs, true);
  MV b(rowmap, numVecs, false);
  b.randomize();
  Kokkos::View<STM*, Kokkos::HostSpace> initNorms("Initial norms", numVecs);
  //Residual norms for starting solution of zero
  b.norm2(initNorms);
  prec->apply(b, x);
  //Compute residual vector = b - Ax
  MV residual(b, Teuchos::Copy);
  A->apply(x, residual, Teuchos::NO_TRANS, -STS::one(), STS::one());
  Kokkos::View<STM*, Kokkos::HostSpace> resNorms("Residual norms", numVecs);
  residual.norm2(resNorms);
  //Make sure all residual norms are smaller than initial
  for(int i = 0; i < numVecs; i++)
  {
    TEST_COMPARE(resNorms(i), <, 0.5 * initNorms(i));
  }
}

#define UNIT_TEST_GROUP_SC_LO_GO( Scalar, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, Test0, Scalar, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, Test1, Scalar, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, Test2, Scalar, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, Test3, Scalar, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, Test4, Scalar, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, Richardson, Scalar, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, SymGaussSeidelZeroRows, Scalar, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, LocalSymGaussSeidelZeroRows, Scalar, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, NotCrsMatrix, Scalar, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, TestDiagonalBlockCrsMatrix, Scalar, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, TestLowerTriangularBlockCrsMatrix, Scalar, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, TestUpperTriangularBlockCrsMatrix, Scalar, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, GS_Crs, Scalar, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, MTSGS, Scalar, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, MTSGS_LongRows, Scalar, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, ClusterMTSGS, Scalar, LO, GO )

  //TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Relaxation, SGS_mult_sweeps, Scalar, LO, GO )
#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// Test all enabled combinations of Scalar (SC), LocalOrdinal (LO),
// and GlobalOrdinal (GO) types.
//
// FIXME (21 Oct 2015) Fix test for complex Scalar types.  (The class
// itself should be fine.)

IFPACK2_INSTANTIATE_SLG_REAL( UNIT_TEST_GROUP_SC_LO_GO )

} // namespace (anonymous)


