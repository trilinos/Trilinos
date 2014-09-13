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

// ***********************************************************************
//
//      Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************


/*! \file Ifpack2_UnitTestContainer.cpp

\brief Ifpack2 Unit test for the Container classes
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#if defined(HAVE_IFPACK2_QD) && !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
#include <qd/dd_real.h>
#endif

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_DenseContainer.hpp>
#include <Ifpack2_SparseContainer.hpp>
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
//template class Ifpack2::SparseContainer<Tpetra::CrsMatrix<float, short, int>,
//                                        Ifpack2::ILUT<Tpetra::CrsMatrix<float, short, short> > >;
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
  typedef Ifpack2::ILUT< Tpetra::CrsMatrix<Scalar,LocalOrdinal,LocalOrdinal,Node>    > ILUTlo;
  typedef Ifpack2::ILUT< Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>   > ILUTgo;

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

  Ifpack2::SparseContainer<CRS, ILUTlo> MyContainer (crsmatrix, localRows);

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
  MyContainer.apply(x,y);

  // Apply raw ILUT
  out << "Applying reference ILUT" << endl;
  prec.apply(x,z);

  // Diff
  out << "Computing results" << endl;
  TEST_COMPARE_FLOATING_ARRAYS(y.get1dView(), z.get1dView(), 1e4*Teuchos::ScalarTraits<Scalar>::eps());

  // Weighted Apply the SparseContainer
  out << "Testing SparseContainer::weightedApply" << endl;
  d.putScalar(1.0);
  MyContainer.weightedApply(x,y,d);

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
  typedef Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> vec_type;
  typedef Ifpack2::DenseContainer<crs_matrix_type, Scalar> container_type;
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
    MyContainer->apply (b, x);
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
  MyContainer->weightedApply (b, y, d);

  out << "Computing results of apply() and weightedApply() "
      << "(they should be the same in this case)" << endl;
  TEST_COMPARE_FLOATING_ARRAYS( x.get1dView(), y.get1dView(), 1.0e2*STS::eps () );
}

// Define the set of unit tests to instantiate in this file.
#define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( SparseContainer, ILUT, Scalar, LocalOrdinal, GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( DenseContainer, FullMatrixSameScalar, Scalar, LocalOrdinal,GlobalOrdinal) \

// Instantiate the unit tests for Scalar=double, LO=int, and GO=int.
UNIT_TEST_GROUP_SCALAR_ORDINAL(double, int, int)

// mfh 03 Sep 2013: See the explicit instantiation at the top of this file.
//#ifndef HAVE_IFPACK2_EXPLICIT_INSTANTIATION
//UNIT_TEST_GROUP_SCALAR_ORDINAL(float, short, int)
//#endif // HAVE_IFPACK2_EXPLICIT_INSTANTIATION

#if defined(HAVE_IFPACK2_QD) && !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
// Instantiate the unit tests for Scalar=dd_real, LO=int, and GO=int.
// For now, we only do this if explicit instantiation is turned off,
// because we don't know how to find out if explicit instantiation is
// enabled for these types.
UNIT_TEST_GROUP_SCALAR_ORDINAL(dd_real, int, int)
#endif

}//namespace <anonymous>

