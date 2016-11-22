/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/


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
  Teuchos::Array<typename CRS::local_ordinal_type> localRows (num_rows_per_proc);
  for (size_t i = 0; i < num_rows_per_proc; ++i) {
    localRows[i] = i;
  }

  out << "SparseContainer constructor" << endl;

  Ifpack2::SparseContainer<ROW, ILUTlo> MyContainer (crsmatrix, localRows);

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
  Array<LocalOrdinal> localRows (numRowsPerProc);
  for (size_t i = 0; i < numRowsPerProc; ++i) {
    localRows[i] = i;
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
    MyContainer = Teuchos::rcp (new container_type (A, localRows));
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
    localSuccess = (errNorm <= 1.0e2 * STS::eps ()) ? 1 : 0;
    globalSuccess = 1;
    reduceAll<int, int> (* (rowMap->getComm ()), REDUCE_MIN,
                         localSuccess, outArg (globalSuccess));
  }

  out << "DenseContainer::weightedApply" << endl;
  d.putScalar (1.0);
  MyContainer->weightedApplyMV(b, y, d);

  out << "Computing results of apply() and weightedApply() "
      << "(they should be the same in this case)" << endl;
  TEST_COMPARE_FLOATING_ARRAYS( x.get1dView(), y.get1dView(), 1.0e2*STS::eps () );
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
  Array<LocalOrdinal> localRows (numRowsPerProc);
  for (size_t i = 0; i < numRowsPerProc; ++i) {
    localRows[i] = i;
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
    const Teuchos::ParameterList params = Teuchos::ParameterList();
    MyContainer = Teuchos::rcp (new container_type (A, localRows));
    MyContainer->setParameters(params);
    localSuccess = 1;
  } catch (std::exception& e) {
    localSuccess = 0;
    cerr << e.what () << endl;
  }
  reduceAll<int, int> (* (rowMap->getComm ()), REDUCE_MIN,
                       localSuccess, outArg (globalSuccess));
  TEST_EQUALITY_CONST( globalSuccess, 1 );

  out << "DenseContainer::setParameters" << endl;
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
    localSuccess = (errNorm <= 1.0e2 * STS::eps ()) ? 1 : 0;
    globalSuccess = 1;
    reduceAll<int, int> (* (rowMap->getComm ()), REDUCE_MIN,
                         localSuccess, outArg (globalSuccess));
  }

  // FIXME: why is this not working?? It seems that the gather call from X to X_local is failing (maybe due to problems with localRows?)
  // TODO second call to apply not working due to problems with local rows???
  //out << "DenseContainer::weightedApply" << endl;
  //d.putScalar (1.0);
  //MyContainer->weightedApply (b, y, d);
  //
  //out << "Computing results of apply() and weightedApply() "
  //    << "(they should be the same in this case)" << endl;
  //TEST_COMPARE_FLOATING_ARRAYS( x.get1dView(), y.get1dView(), 1.0e2*STS::eps () );
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

