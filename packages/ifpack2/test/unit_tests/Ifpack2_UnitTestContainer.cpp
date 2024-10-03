// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/*! \file Ifpack2_UnitTestContainer.cpp

\brief Ifpack2 Unit test for the Container classes
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_DenseContainer.hpp>
#include <Ifpack2_SparseContainer.hpp>
#include <Ifpack2_BandedContainer.hpp>
#include <Ifpack2_ILUT.hpp>

// kgd 04 Sep 2013: Commented out <float,short,short> tests
// that were causing the test build to fail.


// mfh 03 Sep 2013: Instantiate the unit tests for the very special
// case of Scalar=float, LO=short, and GO=int for the global matrix,
// and Scalar=float, LO=short, GO=short for the local matrix.  This
// case is not covered in Ifpack2_SparseContainer.cpp.  We will use
// this case below; we have to put it here because explicit
// instantiations can't be in an anonymous namespace (even if
// namespace-qualified).  We only do this if explicit instantiation is
// turned off, because we don't know how to find out if explicit
// instantiation is enabled for these types.
//#ifndef HAVE_IFPACK2_EXPLICIT_INSTANTIATION
//template class Ifpack2::SparseContainer<Tpetra::RowMatrix<float, short, int>,
//                                        Ifpack2::ILUT<Tpetra::RowMatrix<float, short, short> > >;
//#endif // HAVE_IFPACK2_EXPLICIT_INSTANTIATION


namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(SparseContainer, ILUT, Scalar, LocalOrdinal, GlobalOrdinal)
{
  using Teuchos::RCP;
  using std::endl;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> map_type;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> CRS;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> vec_type;
  typedef Ifpack2::ILUT< Tpetra::RowMatrix<Scalar,LocalOrdinal,LocalOrdinal,Node>    > ILUTlo;
  typedef Ifpack2::ILUT< Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>   > ILUTgo;
  typedef Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> ROW;

//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  out << "Ifpack2::Version(): " << Ifpack2::Version () << endl
      << "Creating test problem" << endl;

  // The simple joy of tridiagonal matrices
  global_size_t num_rows_per_proc = 5;
  RCP<const map_type> rowmap =
    tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node> (num_rows_per_proc);
  Teuchos::RCP<const CRS> crsmatrix =
    tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> (rowmap);
  vec_type x (rowmap), y (rowmap), z (rowmap), d (rowmap);
  Teuchos::ArrayRCP<Scalar> x_ptr = x.get1dViewNonConst();
  Teuchos::ArrayRCP<Scalar> y_ptr = y.get1dViewNonConst();

  out << "Filling x with pseudorandom numbers" << endl;

  // Fill x for all time
  Teuchos::ScalarTraits<double>::seedrandom(24601);
  x.randomize();

  // ====================================== //
  //        Sparse Container + ILUT         //
  // ====================================== //

  // Set IDs to grab the whole matrix
  Teuchos::Array<Teuchos::Array<typename CRS::local_ordinal_type>> blockRows(1);
  blockRows[0].resize(num_rows_per_proc);
  for (size_t i = 0; i < num_rows_per_proc; ++i) {
    blockRows[0][i] = i;
  }

  out << "SparseContainer constructor" << endl;

  Ifpack2::SparseContainer<ROW, ILUTlo> MyContainer (crsmatrix, blockRows, Teuchos::null, false);

  out << "Setting SparseContainer parameters" << endl;

  Teuchos::ParameterList params;
  params.set ("fact: ilut level-of-fill", 1.0);
  params.set ("fact: drop tolerance", 0.0);
  MyContainer.setParameters (params);

  out << "Initializing SparseContainer" << endl;
  MyContainer.initialize ();
  out << "Computing SparseContainer" << endl;
  MyContainer.compute ();

  // Reference ILUT
  out << "Setting up reference ILUT implementation" << endl;
  ILUTgo prec(crsmatrix);
  prec.setParameters(params);
  prec.initialize();
  prec.compute();

  // Apply the SparseContainer
  out << "Applying the SparseContainer" << endl;
  MyContainer.applyMV(x,y);

  // Apply raw ILUT
  out << "Applying reference ILUT" << endl;
  prec.apply(x,z);

  // Diff
  out << "Computing results" << endl;
  TEST_COMPARE_FLOATING_ARRAYS(y.get1dView(), z.get1dView(), 1e4*Teuchos::ScalarTraits<Scalar>::eps());

  // Weighted Apply the SparseContainer
  out << "Testing SparseContainer::weightedApply" << endl;
  d.putScalar(1.0);
  MyContainer.weightedApplyMV(x,y,d);

  // Diff
  out << "Computing results" << endl;
  TEST_COMPARE_FLOATING_ARRAYS(y.get1dView(), z.get1dView(), 1e4*Teuchos::ScalarTraits<Scalar>::eps());

}



// Unit test for DenseContainer.
//
// 1. Create a global test matrix A, exact solution x_exact, and
//    right-hand side b (defined as b = A*x_exact).
// 2. Define the local submatrix as the entire local matrix.
// 3. Apply DenseContainer to approximate the solution x of Ax=b.
//
// If running on only one (MPI) process, x should equal x_exact (to
// within a reasonable tolerance).
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(DenseContainer, FullMatrixSameScalar, Scalar, LocalOrdinal, GlobalOrdinal)
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::cerr;
  using std::endl;
  typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> row_matrix_type;
  typedef Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> vec_type;
  typedef Ifpack2::DenseContainer<row_matrix_type, Scalar> container_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;

  int localSuccess = 1;
  int globalSuccess = 1;

  out << "Ifpack2::Version(): " << Ifpack2::Version () << endl
      << "Creating test problem" << endl;

  global_size_t numRowsPerProc = 5;
  RCP<const map_type> rowMap =
    tif_utest::create_tpetra_map<LocalOrdinal, GlobalOrdinal, Node> (numRowsPerProc);

  out << "Creating the test matrix A" << endl;
  RCP<const crs_matrix_type> A =
    tif_utest::create_test_matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> (rowMap);

  out << "Creating an exact solution vector x" << endl;
  vec_type x_exact (rowMap);
  Teuchos::ScalarTraits<double>::seedrandom (24601);
  x_exact.randomize ();

  out << "Creating the right-hand side vector b of the linear system Ax=b" << endl;
  vec_type b (rowMap);
  A->apply (x_exact, b);

  vec_type y (rowMap);
  vec_type d (rowMap);

  // Set indices to grab the whole matrix
  Array<Array<LocalOrdinal>> blockRows(1);
  blockRows[0].resize(numRowsPerProc);
  for (size_t i = 0; i < numRowsPerProc; ++i) {
    blockRows[0][i] = i;
  }

  // For all the DenseContainer operations, we take special care to
  // ensure that all processes successfully make it through each
  // operation without throwing an exception.  This helped me a lot
  // when I was debugging DenseContainer::extract(), for example.  I
  // found that printing the exception message on each process to cerr
  // (instead of to out) actually let me read the exception message
  // before the test quit.  Otherwise, I wouldn't get to see the
  // exception message.

  out << "DenseContainer constructor" << endl;
  RCP<container_type> MyContainer;
  try {
    MyContainer = Teuchos::rcp (new container_type (A, blockRows, Teuchos::null, false));
    localSuccess = 1;
  } catch (std::exception& e) {
    localSuccess = 0;
    cerr << e.what () << endl;
  }
  reduceAll<int, int> (* (rowMap->getComm ()), REDUCE_MIN,
                       localSuccess, outArg (globalSuccess));
  TEST_EQUALITY_CONST( globalSuccess, 1 );

  out << "DenseContainer::initialize" << endl;
  try {
    MyContainer->initialize ();
    localSuccess = 1;
  } catch (std::exception& e) {
    localSuccess = 0;
    cerr << e.what () << endl;
  }
  reduceAll<int, int> (* (rowMap->getComm ()), REDUCE_MIN,
                       localSuccess, outArg (globalSuccess));
  TEST_EQUALITY_CONST( globalSuccess, 1 );

  out << "DenseContainer::compute" << endl;
  try {
    MyContainer->compute ();
    localSuccess = 1;
  } catch (std::exception& e) {
    localSuccess = 0;
    cerr << e.what () << endl;
  }
  reduceAll<int, int> (* (rowMap->getComm ()), REDUCE_MIN,
                       localSuccess, outArg (globalSuccess));
  TEST_EQUALITY_CONST( globalSuccess, 1 );

  // Apply the DenseContainer to solve the linear system Ax=b for x.
  vec_type x (rowMap);
  x.putScalar (0.0);
  out << "DenseContainer::apply" << endl;
  try {
    MyContainer->applyMV(b, x);
    localSuccess = 1;
  } catch (std::exception& e) {
    localSuccess = 0;
    cerr << e.what () << endl;
  }
  reduceAll<int, int> (* (rowMap->getComm ()), REDUCE_MIN,
                       localSuccess, outArg (globalSuccess));
  TEST_EQUALITY_CONST( globalSuccess, 1 );

  out << "Computing results:" << endl;
  magnitude_type errNorm = STM::zero ();
  {
    vec_type e (x, Teuchos::Copy);
    e.update (-1.0, x_exact, 1.0); // e = x - x_exact
    errNorm = e.norm2 ();
    out << "  ||x - x_exact||_2 = " << errNorm << endl;
  }

  // DenseContainer only solves the global system exactly
  // if there is only one MPI process in the communicator.
  if (rowMap->getComm ()->getSize () == 1) {
    localSuccess = (errNorm <= magnitude_type(1.0e2) * STS::eps ()) ? 1 : 0;
    globalSuccess = 1;
    reduceAll<int, int> (* (rowMap->getComm ()), REDUCE_MIN,
                         localSuccess, outArg (globalSuccess));
  }

  out << "DenseContainer::weightedApply" << endl;
  d.putScalar (1.0);
  MyContainer->weightedApplyMV(b, y, d);

  out << "Computing results of apply() and weightedApply() "
      << "(they should be the same in this case)" << endl;
  TEST_COMPARE_FLOATING_ARRAYS( x.get1dView(), y.get1dView(), magnitude_type(1.0e2) * STS::eps () );
}

// Unit test for BandedContainer.
//
// 1. Create a global test matrix A, exact solution x_exact, and
//    right-hand side b (defined as b = A*x_exact).
// 2. Define the local submatrix as the entire local matrix.
// 3. Apply BandedContainer to approximate the solution x of Ax=b.
//
// If running on only one (MPI) process, x should equal x_exact (to
// within a reasonable tolerance).
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(BandedContainer, FullMatrixSameScalar, Scalar, LocalOrdinal, GlobalOrdinal)
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::cerr;
  using std::endl;
  typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> row_matrix_type;
  typedef Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> vec_type;
  typedef Ifpack2::BandedContainer<row_matrix_type, Scalar> container_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;

  int localSuccess = 1;
  int globalSuccess = 1;

  out << "Ifpack2::Version(): " << Ifpack2::Version () << endl
      << "Creating test problem" << endl;

  global_size_t numRowsPerProc = 5;
  RCP<const map_type> rowMap =
    tif_utest::create_tpetra_map<LocalOrdinal, GlobalOrdinal, Node> (numRowsPerProc);

  out << "Creating the test matrix A" << endl;
  RCP<const crs_matrix_type> A =
    tif_utest::create_banded_matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> (rowMap,3);

  out << "Creating an exact solution vector x" << endl;
  vec_type x_exact (rowMap);
  Teuchos::ScalarTraits<double>::seedrandom (24601);
  x_exact.randomize ();

  out << "Creating the right-hand side vector b of the linear system Ax=b" << endl;
  vec_type b (rowMap);
  A->apply (x_exact, b);

  vec_type y (rowMap);
  vec_type d (rowMap);

  // Set indices to grab the whole matrix
  Array<Array<LocalOrdinal>> blockRows(1);
  blockRows[0].resize(numRowsPerProc);
  for (size_t i = 0; i < numRowsPerProc; ++i) {
    blockRows[0][i] = i;
  }

  // For all the BandedContainer operations, we take special care to
  // ensure that all processes successfully make it through each
  // operation without throwing an exception.  This helped me a lot
  // when I was debugging BandedContainer::extract(), for example.  I
  // found that printing the exception message on each process to cerr
  // (instead of to out) actually let me read the exception message
  // before the test quit.  Otherwise, I wouldn't get to see the
  // exception message.

  out << "BandedContainer constructor" << endl;
  RCP<container_type> MyContainer;
  try {
    MyContainer = Teuchos::rcp (new container_type (A, blockRows, Teuchos::null, false));
    localSuccess = 1;
  } catch (std::exception& e) {
    localSuccess = 0;
    cerr << e.what () << endl;
  }
  reduceAll<int, int> (* (rowMap->getComm ()), REDUCE_MIN,
                       localSuccess, outArg (globalSuccess));
  TEST_EQUALITY_CONST( globalSuccess, 1 );

  out << "BandedContainer::setParameters" << endl;
  try {
    const Teuchos::ParameterList params = Teuchos::ParameterList();
    MyContainer->setParameters(params);
    localSuccess = 1;
  } catch (std::exception& e) {
    localSuccess = 0;
    cerr << e.what () << endl;
  }
  reduceAll<int, int> (* (rowMap->getComm ()), REDUCE_MIN,
                       localSuccess, outArg (globalSuccess));
  TEST_EQUALITY_CONST( globalSuccess, 1 );

  out << "BandedContainer::initialize" << endl;
  try {
    MyContainer->initialize ();
    localSuccess = 1;
  } catch (std::exception& e) {
    localSuccess = 0;
    cerr << e.what () << endl;
  }
  reduceAll<int, int> (* (rowMap->getComm ()), REDUCE_MIN,
                       localSuccess, outArg (globalSuccess));
  TEST_EQUALITY_CONST( globalSuccess, 1 );

  out << "BandedContainer::compute" << endl;
  try {
    MyContainer->compute ();
    localSuccess = 1;
  } catch (std::exception& e) {
    localSuccess = 0;
    cerr << e.what () << endl;
  }
  reduceAll<int, int> (* (rowMap->getComm ()), REDUCE_MIN,
                       localSuccess, outArg (globalSuccess));
  TEST_EQUALITY_CONST( globalSuccess, 1 );

  // Apply the DenseContainer to solve the linear system Ax=b for x.
  vec_type x (rowMap);
  x.putScalar (0.0);
  out << "BandedContainer::apply" << endl;
  try {
    MyContainer->applyMV(b, x);
    localSuccess = 1;
  } catch (std::exception& e) {
    localSuccess = 0;
    cerr << e.what () << endl;
  }
  reduceAll<int, int> (* (rowMap->getComm ()), REDUCE_MIN,
                       localSuccess, outArg (globalSuccess));
  TEST_EQUALITY_CONST( globalSuccess, 1 );

  out << "Computing results:" << endl;
  magnitude_type errNorm = STM::zero ();
  {
    vec_type e (x, Teuchos::Copy);
    e.update (-1.0, x_exact, 1.0); // e = x - x_exact
    errNorm = e.norm2 ();
    out << "  ||x - x_exact||_2 = " << errNorm << endl;
  }

  // DenseContainer only solves the global system exactly
  // if there is only one MPI process in the communicator.
  if (rowMap->getComm ()->getSize () == 1) {
    localSuccess = (errNorm <= magnitude_type(1.0e2) * STS::eps ()) ? 1 : 0;
    globalSuccess = 1;
    reduceAll<int, int> (* (rowMap->getComm ()), REDUCE_MIN,
                         localSuccess, outArg (globalSuccess));
  }

  out << "DenseContainer::weightedApply" << endl;
  d.putScalar (1.0);
  MyContainer->weightedApplyMV(b, y, d);

  out << "Computing results of apply() and weightedApply() "
      << "(they should be the same in this case)" << endl;
  TEST_COMPARE_FLOATING_ARRAYS( x.get1dView(), y.get1dView(), magnitude_type(1.0e2) * STS::eps () );
}

// Define the set of unit tests to instantiate in this file.
#define UNIT_TEST_GROUP_SC_LO_GO(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( SparseContainer, ILUT, Scalar, LocalOrdinal, GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( BandedContainer, FullMatrixSameScalar, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( DenseContainer, FullMatrixSameScalar, Scalar, LocalOrdinal,GlobalOrdinal) \

// NOTE (mfh 21 Oct 2015) This test is special, because it wants to
// use two different GlobalOrdinal types, but Ifpack2 does not do ETI
// for SparseContainer for that case.  I think this reflects a flaw in
// the design of SparseContainer rather than a flaw in Ifpack2's ETI
// system.  It's also worrisome that the test never actually exercised
// GO != LO, because it was only ever instantiated for LO = int and GO
// = int.  Anyway, I'll protect that one instantiation for now.

#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)

// Instantiate the unit tests for Scalar=double, LO=int, and GO=int.
UNIT_TEST_GROUP_SC_LO_GO(double, int, int)

#endif // defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)

// mfh 03 Sep 2013: See the explicit instantiation at the top of this file.
//#ifndef HAVE_IFPACK2_EXPLICIT_INSTANTIATION
//UNIT_TEST_GROUP_SCALAR_ORDINAL(float, short, int)
//#endif // HAVE_IFPACK2_EXPLICIT_INSTANTIATION

} // namespace (anonymous)

