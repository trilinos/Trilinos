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

/*! \file Ifpack2_UnitTestRILUK.cpp

\brief Ifpack2 parallel unit tests for the RILUK template.
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "TpetraExt_MatrixMatrix.hpp"

#if defined(HAVE_IFPACK2_QD) && !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
#include <qd/dd_real.h>
#endif

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_Experimental_RBILUK.hpp>
#include <Tpetra_RowMatrix.hpp>
#include <Tpetra_Experimental_BlockMultiVector.hpp>
#include <Tpetra_Experimental_BlockCrsMatrix.hpp>

namespace {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  typedef typename Tpetra::global_size_t GST;
  typedef typename tif_utest::Node Node;

#define IFPACK2RBILUK_REPORT_GLOBAL_ERR( WHAT_STRING ) do { \
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

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2RBILUK, TestLowerTriangularBlockCrsMatrix, Scalar, LO, GO)
{
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::endl;
  typedef Tpetra::Experimental::BlockCrsMatrix<Scalar,LO,GO,Node> block_crs_matrix_type;
  typedef Tpetra::CrsGraph<LO,GO,Node> crs_graph_type;
  typedef Tpetra::Experimental::BlockMultiVector<Scalar,LO,GO,Node> BMV;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
  typedef Ifpack2::Experimental::RBILUK<block_crs_matrix_type> prec_type;
  int lclSuccess = 1;
  int gblSuccess = 1;
  std::ostringstream errStrm; // for error collection

  out << "Ifpack2::Experimental::RBILUK lower triangular BlockCrsMatrix test" << endl;

  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
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
  IFPACK2RBILUK_REPORT_GLOBAL_ERR( "create_dense_local_graph" );

  const int blockSize = 5;
  RCP<block_crs_matrix_type> bcrsmatrix;
  try {
    bcrsmatrix = rcp_const_cast<block_crs_matrix_type> (tif_utest::create_triangular_matrix<Scalar, LO, GO, Node, true> (crsgraph, blockSize));
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": create_triangular_matrix threw exception: "
            << e.what () << endl;
  }
  IFPACK2RBILUK_REPORT_GLOBAL_ERR( "create_triangular_matrix" );

  try {
    bcrsmatrix->computeDiagonalGraph ();
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": computeDiagonalGraph() threw exception: "
            << e.what () << endl;
  }
  IFPACK2RBILUK_REPORT_GLOBAL_ERR( "computeDiagonalGraph()" );

  RCP<prec_type> prec;
  try {
    prec = rcp (new prec_type (bcrsmatrix));
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": Preconditioner constructor threw exception: "
            << e.what () << endl;
  }
  IFPACK2RBILUK_REPORT_GLOBAL_ERR( "Preconditioner constructor" );

  Teuchos::ParameterList params;
  params.set("fact: iluk level-of-fill", (LO) 0);
  params.set("fact: relax value", 0.0);

  try {
    prec->setParameters (params);
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": prec->setParameters() threw exception: "
            << e.what () << endl;
  }
  IFPACK2RBILUK_REPORT_GLOBAL_ERR( "prec->setParameters()" );

  try {
    prec->initialize ();
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": prec->initialize() threw exception: "
            << e.what () << endl;
  }
  IFPACK2RBILUK_REPORT_GLOBAL_ERR( "prec->initialize()" );

  try {
    prec->compute ();
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": prec->compute() threw exception: "
            << e.what () << endl;
  }
  IFPACK2RBILUK_REPORT_GLOBAL_ERR( "prec->compute()" );

  BMV xBlock (* (crsgraph->getRowMap ()), blockSize, 1);
  BMV yBlock (* (crsgraph->getRowMap ()), blockSize, 1);
  MV x = xBlock.getMultiVectorView ();
  MV y = yBlock.getMultiVectorView ();
  x.putScalar (Teuchos::ScalarTraits<Scalar>::one ());

  TEST_EQUALITY( x.getMap()->getNodeNumElements (), blockSize * num_rows_per_proc );
  TEST_EQUALITY( y.getMap ()->getNodeNumElements (), blockSize * num_rows_per_proc );

  try {
    prec->apply (x, y);
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": prec->apply(x, y) threw exception: "
            << e.what () << endl;
  }
  IFPACK2RBILUK_REPORT_GLOBAL_ERR( "prec->apply(x, y)" );

  Teuchos::Array<Scalar> exactSol(num_rows_per_proc);
  exactSol[0] = 0.5;
  exactSol[1] = -0.25;
  exactSol[2] = 0.625;

  for (size_t k = 0; k < num_rows_per_proc; ++k) {
    LO lcl_row = k;
    typename BMV::little_vec_type ylcl = yBlock.getLocalBlock(lcl_row,0);
    Scalar* yb = ylcl.getRawPtr();
    for (int j = 0; j < blockSize; ++j) {
      TEST_FLOATING_EQUALITY(yb[j],exactSol[k],1e-14);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2RBILUK, TestUpperTriangularBlockCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef Tpetra::Experimental::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> block_crs_matrix_type;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> crs_graph_type;
  typedef Tpetra::Experimental::BlockMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> BMV;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  typedef Ifpack2::Experimental::RBILUK<block_crs_matrix_type> prec_type;

  out << "Ifpack2::Version(): " << Ifpack2::Version () << std::endl;

  const int num_rows_per_proc = 3;
  const int blockSize = 5;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  RCP<crs_graph_type> crsgraph =
    tif_utest::create_dense_local_graph<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  RCP<block_crs_matrix_type> bcrsmatrix =
    rcp_const_cast<block_crs_matrix_type> (tif_utest::create_triangular_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,false> (crsgraph, blockSize));
  bcrsmatrix->computeDiagonalGraph();

  prec_type prec (bcrsmatrix);

  Teuchos::ParameterList params;
  params.set("fact: iluk level-of-fill", (LocalOrdinal) 0);
  params.set("fact: relax value", 0.0);
  prec.setParameters(params);

  prec.initialize();
  TEST_NOTHROW(prec.compute());

  BMV xBlock (*crsgraph->getRowMap (), blockSize, 1);
  BMV yBlock (*crsgraph->getRowMap (), blockSize, 1);
  MV x = xBlock.getMultiVectorView ();
  MV y = yBlock.getMultiVectorView ();
  x.putScalar (Teuchos::ScalarTraits<Scalar>::one ());

  TEST_EQUALITY(x.getMap()->getNodeNumElements(), blockSize*num_rows_per_proc);
  TEST_EQUALITY(y.getMap()->getNodeNumElements(), blockSize*num_rows_per_proc);
  TEST_NOTHROW(prec.apply(x, y));

  Teuchos::Array<Scalar> exactSol(num_rows_per_proc);
  exactSol[0] = 0.625;
  exactSol[1] = -0.25;
  exactSol[2] = 0.5;

  for (int k = 0; k < num_rows_per_proc; ++k) {
    typename BMV::little_vec_type ylcl = yBlock.getLocalBlock(k,0);
    Scalar* yb = ylcl.getRawPtr();
    for (int j = 0; j < blockSize; ++j) {
      TEST_FLOATING_EQUALITY(yb[j],exactSol[k],1e-14);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2RBILUK, TestFullLocalBlockCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef Tpetra::Experimental::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> block_crs_matrix_type;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> crs_graph_type;
  typedef Tpetra::Experimental::BlockMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> BMV;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  typedef Ifpack2::Experimental::RBILUK<block_crs_matrix_type> prec_type;

  out << "Ifpack2::Version(): " << Ifpack2::Version () << std::endl;

  const int num_rows_per_proc = 3;
  const int blockSize = 5;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  RCP<crs_graph_type> crsgraph =
    tif_utest::create_dense_local_graph<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  RCP<block_crs_matrix_type> bcrsmatrix =
    rcp_const_cast<block_crs_matrix_type> (tif_utest::create_full_local_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> (crsgraph, blockSize));
  bcrsmatrix->computeDiagonalGraph();

  prec_type prec (bcrsmatrix);

  Teuchos::ParameterList params;
  params.set("fact: iluk level-of-fill", (LocalOrdinal) 0);
  params.set("fact: relax value", 0.0);
  prec.setParameters(params);

  prec.initialize();
  TEST_NOTHROW(prec.compute());

  BMV xBlock (*crsgraph->getRowMap (), blockSize, 1);
  BMV yBlock (*crsgraph->getRowMap (), blockSize, 1);
  MV x = xBlock.getMultiVectorView ();
  MV y = yBlock.getMultiVectorView ();
  x.putScalar (Teuchos::ScalarTraits<Scalar>::one ());

  TEST_EQUALITY(x.getMap()->getNodeNumElements(), blockSize*num_rows_per_proc);
  TEST_EQUALITY(y.getMap()->getNodeNumElements(), blockSize*num_rows_per_proc);
  TEST_NOTHROW(prec.apply(x, y));

  Teuchos::Array<Scalar> exactSol(num_rows_per_proc);
  exactSol[0] = 3.0/14.0;
  exactSol[1] = -4.0/21.0;
  exactSol[2] = 2.0/7.0;

  for (int k = 0; k < num_rows_per_proc; ++k) {
    typename BMV::little_vec_type ylcl = yBlock.getLocalBlock(k,0);
    Scalar* yb = ylcl.getRawPtr();
    for (int j = 0; j < blockSize; ++j) {
      TEST_FLOATING_EQUALITY(yb[j],exactSol[k],1e-14);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2RBILUK, TestBandedBlockCrsMatrixWithDropping, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef Tpetra::Experimental::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> block_crs_matrix_type;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> crs_graph_type;
  typedef Tpetra::Experimental::BlockMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> BMV;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  typedef Ifpack2::Experimental::RBILUK<block_crs_matrix_type> prec_type;
  typedef Ifpack2::RILUK<crs_matrix_type> prec_crs_type;

  out << "Ifpack2::Version(): " << Ifpack2::Version () << std::endl;

  const int num_rows_per_proc = 10;
  const int blockSize = 1;

  const int lof = 0;
  const size_t rbandwidth = lof+2+2;
  RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph =
    tif_utest::create_banded_graph<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc, rbandwidth);
  RCP<block_crs_matrix_type> bcrsmatrix =
    rcp_const_cast<block_crs_matrix_type> (tif_utest::create_banded_block_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> (crsgraph, blockSize, rbandwidth));
  bcrsmatrix->computeDiagonalGraph();

  Ifpack2::Experimental::RBILUK<block_crs_matrix_type> prec (bcrsmatrix);

  Teuchos::ParameterList params;
  params.set("fact: iluk level-of-fill", (LocalOrdinal) lof);
  params.set("fact: relax value", 0.0);
  prec.setParameters(params);

  prec.initialize();
  TEST_NOTHROW(prec.compute());

  BMV xBlock (*crsgraph->getRowMap (), blockSize, 1);
  BMV yBlock (*crsgraph->getRowMap (), blockSize, 1);
  BMV zBlock (*crsgraph->getRowMap (), blockSize, 1);
  MV x = xBlock.getMultiVectorView ();
  MV y = yBlock.getMultiVectorView ();
  MV z = zBlock.getMultiVectorView ();

  x.randomize();

  TEST_EQUALITY(x.getMap()->getNodeNumElements(), blockSize*num_rows_per_proc);
  TEST_EQUALITY(y.getMap()->getNodeNumElements(), blockSize*num_rows_per_proc);
  TEST_NOTHROW(prec.apply(x, y));

  RCP<const crs_matrix_type > crsmatrix = tif_utest::create_banded_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(crsgraph->getRowMap(),rbandwidth);

  prec_crs_type prec_crs (crsmatrix);
  prec_crs.setParameters(params);

  prec_crs.initialize();
  prec_crs.compute();

  prec_crs.apply(x, z);

  Teuchos::ArrayRCP<const Scalar> zview = z.get1dView();
  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
  for (int k = 0; k < num_rows_per_proc; ++k)
  {
    Scalar yb = yview[k];
    Scalar zb = zview[k];
    TEST_FLOATING_EQUALITY(yb, zb, 1e-14);
  }


}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2RBILUK, TestBlockMatrixOps, Scalar, LocalOrdinal, GlobalOrdinal)
{

  typedef Tpetra::Experimental::LittleBlock<Scalar,LocalOrdinal> little_block_type;
  typedef Tpetra::Experimental::LittleVector<Scalar,LocalOrdinal> little_vec_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  Ifpack2::Experimental::BlockMatrixOperations<Scalar,Scalar> blockOps;

  const int blockSize = 5;
  const int blockMatSize = blockSize*blockSize;


  Teuchos::Array<Scalar> identityMatrix(blockMatSize,0.0);

  Teuchos::Array<Scalar> aMatrix(blockMatSize,0.0);
  Teuchos::Array<Scalar> bMatrix(blockMatSize,0.0);
  Teuchos::Array<Scalar> cMatrix(blockMatSize,0.0);
  Teuchos::Array<Scalar> exactMatrix(blockMatSize,0.0);

  for (int i = 0; i < blockSize; ++i)
  {
    identityMatrix[i*(blockSize+1)] = 1.0;
  }

  aMatrix[0] = 1;
  aMatrix[1] = -2;
  aMatrix[2] = 3;
  aMatrix[3] = -4;
  aMatrix[4] = 5;

  aMatrix[5] = 1;
  aMatrix[6] = -2;
  aMatrix[7] = 3;
  aMatrix[8] = -4;
  aMatrix[9] = 5;

  aMatrix[10] = 1;
  aMatrix[11] = -2;
  aMatrix[12] = 3;
  aMatrix[13] = -4;
  aMatrix[14] = 5;

  aMatrix[15] = 1;
  aMatrix[16] = -2;
  aMatrix[17] = 3;
  aMatrix[18] = -4;
  aMatrix[19] = 5;

  aMatrix[20] = 1;
  aMatrix[21] = -2;
  aMatrix[22] = 3;
  aMatrix[23] = -4;
  aMatrix[24] = 5;

  bMatrix[0] = -1;
  bMatrix[1] = 2;
  bMatrix[2] = -3;
  bMatrix[3] = 4;
  bMatrix[4] = -5;

  bMatrix[5] = -1;
  bMatrix[6] = 2;
  bMatrix[7] = -3;
  bMatrix[8] = 4;
  bMatrix[9] = -5;

  bMatrix[10] = -1;
  bMatrix[11] = 2;
  bMatrix[12] = -3;
  bMatrix[13] = 4;
  bMatrix[14] = -5;

  bMatrix[20] = -1;
  bMatrix[21] = 2;
  bMatrix[22] = -3;
  bMatrix[23] = 4;
  bMatrix[24] = -5;

  blockOps.square_matrix_matrix_multiply(aMatrix.getRawPtr(), identityMatrix.getRawPtr(), cMatrix.getRawPtr(), blockSize);

  for (int k = 0; k < blockMatSize; ++k)
  {
    TEST_FLOATING_EQUALITY(aMatrix[k], cMatrix[k], 1e-14);
  }

  blockOps.square_matrix_matrix_multiply(aMatrix.getRawPtr(), aMatrix.getRawPtr(), cMatrix.getRawPtr(), blockSize, -1.0, 1.0);
  blockOps.square_matrix_matrix_multiply(identityMatrix.getRawPtr(), identityMatrix.getRawPtr(), cMatrix.getRawPtr(), blockSize, 1.0, 1.0);
  blockOps.square_matrix_matrix_multiply(aMatrix.getRawPtr(), bMatrix.getRawPtr(), cMatrix.getRawPtr(), blockSize, 1.0, 1.0);

  exactMatrix[0] = -8.00000000000000;
  exactMatrix[1] = 18.0000000000000;
  exactMatrix[2] = -27.0000000000000;
  exactMatrix[3] = 36.0000000000000;
  exactMatrix[4] = -45.0000000000000;
  exactMatrix[5] = -9.00000000000000;
  exactMatrix[6] = 19.0000000000000;
  exactMatrix[7] = -27.0000000000000;
  exactMatrix[8] = 36.0000000000000;
  exactMatrix[9] = -45.0000000000000;
  exactMatrix[10] = -9.00000000000000;
  exactMatrix[11] = 18.0000000000000;
  exactMatrix[12] = -26.0000000000000;
  exactMatrix[13] = 36.0000000000000;
  exactMatrix[14] = -45.0000000000000;
  exactMatrix[15] = -9.00000000000000;
  exactMatrix[16] = 18.0000000000000;
  exactMatrix[17] = -27.0000000000000;
  exactMatrix[18] = 37.0000000000000;
  exactMatrix[19] = -45.0000000000000;
  exactMatrix[20] = -9.00000000000000;
  exactMatrix[21] = 18.0000000000000;
  exactMatrix[22] = -27.0000000000000;
  exactMatrix[23] = 36.0000000000000;
  exactMatrix[24] = -44.0000000000000;

  for (int k = 0; k < blockMatSize; ++k)
  {
    TEST_FLOATING_EQUALITY(exactMatrix[k], cMatrix[k], 1e-14);
  }

  typedef typename GetLapackType<Scalar>::lapack_scalar_type LST;
  typedef typename GetLapackType<Scalar>::lapack_type lapack_type;

  lapack_type lapack;

  const LocalOrdinal rowStride = blockSize;
  const LocalOrdinal colStride = 1;

  little_block_type dMat(cMatrix.getRawPtr(),blockSize,rowStride,colStride);
  LST * d_raw = dMat.getRawPtr();

  Teuchos::Array<int> ipiv(blockSize);
  Teuchos::Array<Scalar> work(1);
  int lapackInfo;
  for (int k = 0; k < blockSize; ++k) {
    ipiv[k] = 0;
  }

  lapack.GETRF(blockSize, blockSize, d_raw, blockSize, ipiv.getRawPtr(), &lapackInfo);
  TEUCHOS_TEST_FOR_EXCEPTION(
    lapackInfo != 0, std::runtime_error, "Ifpack2::Experimental::RBILUK::compute: "
    "lapackInfo = " << lapackInfo << " which indicates an error in the factorization GETRF.");

  int lwork = -1;
  lapack.GETRI(blockSize, d_raw, blockSize, ipiv.getRawPtr(), work.getRawPtr(), lwork, &lapackInfo);
  TEUCHOS_TEST_FOR_EXCEPTION(
    lapackInfo != 0, std::runtime_error, "Ifpack2::Experimental::RBILUK::compute: "
    "lapackInfo = " << lapackInfo << " which indicates an error in the matrix inverse GETRI.");

  typedef typename Kokkos::Details::ArithTraits<Scalar>::mag_type ImplMagnitudeType;
  ImplMagnitudeType worksize = Kokkos::Details::ArithTraits<Scalar>::magnitude(work[0]);
  lwork = static_cast<int>(worksize);
  work.resize(lwork);
  lapack.GETRI(blockSize, d_raw, blockSize, ipiv.getRawPtr(), work.getRawPtr(), lwork, &lapackInfo);
  TEUCHOS_TEST_FOR_EXCEPTION(
    lapackInfo != 0, std::runtime_error, "Ifpack2::Experimental::RBILUK::compute: "
    "lapackInfo = " << lapackInfo << " which indicates an error in the matrix inverse GETRI.");

  Teuchos::Array<Scalar> onevec(blockSize,STS::one());
  Teuchos::Array<Scalar> computeSolution(blockSize,STS::zero());
  Teuchos::Array<Scalar> exactSolution(blockSize,STS::zero());

  exactSolution[0] = -0.0384615384615392;
  exactSolution[1] = -0.0384615384615377;
  exactSolution[2] = -0.0384615384615388;
  exactSolution[3] = -0.0384615384615381;
  exactSolution[4] = -0.0384615384615385;

  little_vec_type bval(onevec.getRawPtr(),blockSize,1);
  little_vec_type xval(computeSolution.getRawPtr(),blockSize,1);

  xval.matvecUpdate(1.0,dMat,bval);

  for (int i = 0; i < blockSize; ++i)
    TEST_FLOATING_EQUALITY(exactSolution[i], computeSolution[i], 1e-13);


}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2RBILUK, TestDiagonalBlockCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal)
{
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::endl;
  typedef Tpetra::Experimental::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> block_crs_matrix_type;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> crs_graph_type;
  typedef Tpetra::Experimental::BlockMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> BMV;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  typedef Ifpack2::Experimental::RBILUK<block_crs_matrix_type> prec_type;
  std::ostringstream errStrm; // for error collection

  out << "Ifpack2::RBILUK diagonal block matrix test" << endl;

  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();

  const int num_rows_per_proc = 5;
  const int blockSize = 3;
  RCP<crs_graph_type> crsgraph =
    tif_utest::create_diagonal_graph<LocalOrdinal,GlobalOrdinal,Node> (num_rows_per_proc);

  RCP<block_crs_matrix_type> bcrsmatrix;
  bcrsmatrix = rcp_const_cast<block_crs_matrix_type> (tif_utest::create_block_diagonal_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> (crsgraph, blockSize));
  bcrsmatrix->computeDiagonalGraph ();

  Ifpack2::Experimental::RBILUK<block_crs_matrix_type> prec (bcrsmatrix);

  Teuchos::ParameterList params;
  params.set("fact: iluk level-of-fill", (LocalOrdinal) 0);
  params.set("fact: relax value", 0.0);
  prec.setParameters(params);

  prec.initialize();
  TEST_NOTHROW(prec.compute());

  BMV xBlock (*crsgraph->getRowMap (), blockSize, 1);
  BMV yBlock (*crsgraph->getRowMap (), blockSize, 1);
  MV x = xBlock.getMultiVectorView ();
  MV y = yBlock.getMultiVectorView ();
  x.putScalar (Teuchos::ScalarTraits<Scalar>::one ());

  prec.apply (x, y);

  const Scalar exactSol = 0.2;

  for (int k = 0; k < num_rows_per_proc; ++k) {
    typename BMV::little_vec_type ylcl = yBlock.getLocalBlock(k,0);
    Scalar* yb = ylcl.getRawPtr();
    for (int j = 0; j < blockSize; ++j) {
      TEST_FLOATING_EQUALITY(yb[j],exactSol,1e-14);
    }
  }
}

# define UNIT_TEST_GROUP_BLOCK_LGN( Scalar, LocalOrdinal, GlobalOrdinal ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2RBILUK, TestLowerTriangularBlockCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2RBILUK, TestUpperTriangularBlockCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2RBILUK, TestBandedBlockCrsMatrixWithDropping, Scalar, LocalOrdinal, GlobalOrdinal) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2RBILUK, TestBlockMatrixOps, Scalar, LocalOrdinal, GlobalOrdinal) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2RBILUK, TestDiagonalBlockCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2RBILUK, TestFullLocalBlockCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal)

typedef Tpetra::MultiVector<>::scalar_type scalar_type;
typedef Tpetra::MultiVector<>::local_ordinal_type local_ordinal_type;
typedef Tpetra::MultiVector<>::global_ordinal_type global_ordinal_type;

UNIT_TEST_GROUP_BLOCK_LGN(scalar_type, local_ordinal_type, global_ordinal_type)

}//namespace <anonymous>

