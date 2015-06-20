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


/*! \file Ifpack2_UnitTestRelaxation.cpp

\brief Ifpack2 Unit test for the OverlappingRowMatrix template.
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <iostream>

#if defined(HAVE_IFPACK2_QD) && !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
#include <qd/dd_real.h>
#endif

// Xpetra / Galeri
#ifdef HAVE_IFPACK2_XPETRA
#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_TpetraMap.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraMatrixTypes.hpp>
#endif


#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_OverlappingRowMatrix.hpp>
#include <Ifpack2_CreateOverlapGraph.hpp>

namespace {

using Teuchos::ArrayRCP;
using Teuchos::Comm;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::reduceAll;
using Teuchos::REDUCE_MIN;
using std::endl;
typedef tif_utest::Node Node;
typedef Tpetra::global_size_t GST;

// Macro used inside the unit test below.  It tests for global error,
// and if so, prints each process' error message and quits the test
// early.
//
// 'out' only prints on Process 0.  It's really not OK for other
// processes to print to stdout, but it usually works and we need to
// do it for debugging.
#define IFPACK2OVERLAPPINGROWMATRIX_REPORT_GLOBAL_ERR( WHAT_STRING ) do { \
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
    return; \
  } \
} while (false)


TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2OverlappingRowMatrix, Test0, Scalar, LO, GO)
{
  out << "Ifpack2::OverlappingRowMatrix unit test" << endl;
  Teuchos::OSTab tab0 (out);

#ifndef HAVE_IFPACK2_XPETRA
  out << "This test requires building with Xpetra enabled." << endl;
#else
  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node>       CrsType;
  typedef Xpetra::TpetraCrsMatrix<Scalar,LO,GO,Node> XCrsType;
  typedef Xpetra::Map<LO,GO,Node>                    XMapType;
  typedef Xpetra::MultiVector<Scalar,LO,GO,Node>     XMVectorType;
  typedef Tpetra::Vector<Scalar,LO,GO,Node>          VectorType;

  int lclSuccess = 1;
  int gblSuccess = 1;
  std::ostringstream errStrm; // for error collection

  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  Teuchos::CommandLineProcessor clp;
  Xpetra::Parameters xpetraParameters (clp);
  Teuchos::ParameterList GaleriList;
  int nx = 100;
  //    int nx = 6;
  size_t numElementsPerProc = nx*nx;
  GaleriList.set("nx", static_cast<GO> (nx));
  GaleriList.set("ny", static_cast<GO> (nx * numProcs));
  GaleriList.set("n", static_cast<GO> (numElementsPerProc*numProcs));

  // Short circuit --- this test should only be run in parallel.
  if (numProcs == 1) {
    out << "This test is only meaningful if run with multiple MPI processes."
        << endl;
    return;
  }

  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
  RCP<XMapType> xmap;
  RCP<Galeri::Xpetra::Problem<XMapType, XCrsType, XMVectorType> > Pr;
  RCP<XCrsType> XA;
  RCP<CrsType> A;
  try {
    xmap = Xpetra::MapFactory<LO, GO>::Build (xpetraParameters.GetLib (), INVALID,
                                              numElementsPerProc, 0, comm);
    Pr = Galeri::Xpetra::BuildProblem<Scalar, LO, GO, XMapType, XCrsType, XMVectorType> (string ("Laplace2D"), xmap, GaleriList);
    XA = Pr->BuildMatrix ();
    A = XA->getTpetra_CrsMatrixNonConst ();
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Problem setup (Galeri or Xpetra) threw exception: " << e.what () << endl;
  }
  IFPACK2OVERLAPPINGROWMATRIX_REPORT_GLOBAL_ERR( "Problem setup" );

  if (A.is_null ()) {
    lclSuccess = 0;
  }
  IFPACK2OVERLAPPINGROWMATRIX_REPORT_GLOBAL_ERR( "A is null on at least one process." );

  VectorType X (A->getRowMap ());
  VectorType Y (A->getRowMap ());
  VectorType Z (A->getRowMap ());
  ArrayRCP<ArrayRCP<Scalar> > x_ptr = X.get2dViewNonConst ();

  const int OverlapLevel = 5;

  // ======================================== //
  // Build the overlapping matrix using class //
  // Ifpack2::OverlappingRowMatrix.           //
  // ======================================== //
  RCP<Ifpack2::OverlappingRowMatrix<CrsType> > B;

  try {
    B = rcp (new Ifpack2::OverlappingRowMatrix<CrsType> (A, OverlapLevel));
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Ifpack2::OverlappingRowMatrix constructor threw exception: " << e.what () << endl;
  }
  IFPACK2OVERLAPPINGROWMATRIX_REPORT_GLOBAL_ERR( "Ifpack2::OverlappingRowMatrix constructor" );

  size_t NumGlobalRowsB = B->getGlobalNumRows ();
  size_t NumGlobalNonzerosB = B->getGlobalNumEntries ();

  for (LO i = 0 ; i < static_cast<LO> (A->getNodeNumRows ()); ++i) {
    x_ptr[0][i] = 1.0 * A->getRowMap ()->getGlobalElement (i);
  }
  Y.putScalar (0.0);

  VectorType ExtX_B (B->getRowMap ());
  VectorType ExtY_B (B->getRowMap ());
  ExtY_B.putScalar (0.0);

  try {
    B->importMultiVector (X,ExtX_B);
    B->apply (ExtX_B,ExtY_B);
    B->exportMultiVector (ExtY_B, Y, Tpetra::ADD);
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Import, apply B, and Export: " << e.what () << endl;
  }
  IFPACK2OVERLAPPINGROWMATRIX_REPORT_GLOBAL_ERR( "Import, apply B, and Export" );

  // ================================================== //
  // Build the overlapping graph using                  //
  // CreateOverlappingMatrix.                            //
  // ================================================== //
  RCP<const CrsType> C;
  try {
    C = Ifpack2::createOverlapMatrix<CrsType> (A, OverlapLevel);
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Ifpack2::createOverlapMatrix threw an exception" << e.what () << endl;
  }
  IFPACK2OVERLAPPINGROWMATRIX_REPORT_GLOBAL_ERR( "Ifpack2::createOverlapMatrix" );

  // simple checks on global quantities
  const size_t NumGlobalRowsC = C->getGlobalNumRows ();
  const size_t NumGlobalNonzerosC = C->getGlobalNumEntries ();

  TEST_EQUALITY( NumGlobalRowsB, NumGlobalRowsC );
  TEST_EQUALITY( NumGlobalNonzerosB, NumGlobalNonzerosC );

  try {
    C->apply (X, Z);
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Ifpack2::OverlappingRowMatrix::apply threw an exception" << e.what () << endl;
  }
  IFPACK2OVERLAPPINGROWMATRIX_REPORT_GLOBAL_ERR( "Ifpack2::OverlappingRowMatrix::apply" );

  TEST_COMPARE_FLOATING_ARRAYS( Y.get1dView (), Z.get1dView (), 1e4 * Teuchos::ScalarTraits<Scalar>::eps () );
#endif // HAVE_IFPACK2_XPETRA
}

#define UNIT_TEST_GROUP_SCALAR_ORDINAL( Scalar, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2OverlappingRowMatrix, Test0, Scalar, LO, GO )


UNIT_TEST_GROUP_SCALAR_ORDINAL(double, int, int)

#if defined(HAVE_IFPACK2_QD) && !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
UNIT_TEST_GROUP_SCALAR_ORDINAL(dd_real, int, int)
#endif

}//namespace <anonymous>

