// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Ifpack2_UnitTestRILUK.cpp

\brief Ifpack2 parallel unit tests for the RILUK template.
*/


#include "Teuchos_UnitTestHarness.hpp"

#include "Ifpack2_UnitTestHelpers.hpp"
#include "Ifpack2_AdditiveSchwarz.hpp"
#include "Ifpack2_Experimental_RBILUK.hpp"
#include "Ifpack2_Version.hpp"

#include "MatrixMarket_Tpetra.hpp"
#include "TpetraExt_MatrixMatrix.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_BlockCrsMatrix.hpp"
#include "Tpetra_BlockMultiVector.hpp"
#include "Tpetra_BlockView.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "Tpetra_RowMatrix.hpp"

#include <iostream>

#if defined(HAVE_IFPACK2_QD) && !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
#include <qd/dd_real.h>
#endif

namespace { // (anonymous)

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

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(RBILUK, LowerTriangularBlockCrsMatrix, Scalar, LO, GO)
{
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::endl;
  typedef Tpetra::BlockCrsMatrix<Scalar,LO,GO,Node> block_crs_matrix_type;
  typedef Tpetra::CrsGraph<LO,GO,Node> crs_graph_type;
  typedef Tpetra::BlockMultiVector<Scalar,LO,GO,Node> BMV;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
  typedef Ifpack2::Experimental::RBILUK<block_crs_matrix_type> prec_type;
  int lclSuccess = 1;
  int gblSuccess = 1;
  std::ostringstream errStrm; // for error collection

  out << "Ifpack2::Experimental::RBILUK lower triangular BlockCrsMatrix test" << endl;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
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

  RCP<prec_type> prec;
  try {
    RCP<const block_crs_matrix_type> const_bcrsmatrix(bcrsmatrix);
    prec = rcp (new prec_type (const_bcrsmatrix));
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

  TEST_EQUALITY( x.getMap()->getLocalNumElements (), blockSize * num_rows_per_proc );
  TEST_EQUALITY( y.getMap ()->getLocalNumElements (), blockSize * num_rows_per_proc );

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
    auto ylcl = yBlock.getLocalBlockHost(lcl_row, 0, Tpetra::Access::ReadOnly);
    for (int j = 0; j < blockSize; ++j) {
      TEST_FLOATING_EQUALITY(ylcl(j),exactSol[k],1e-14);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(RBILUK, UpperTriangularBlockCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> block_crs_matrix_type;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> crs_graph_type;
  typedef Tpetra::BlockMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> BMV;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  typedef Ifpack2::Experimental::RBILUK<block_crs_matrix_type> prec_type;

  out << "Ifpack2::Version(): " << Ifpack2::Version () << std::endl;

  const int num_rows_per_proc = 3;
  const int blockSize = 5;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  RCP<crs_graph_type> crsgraph =
    tif_utest::create_dense_local_graph<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  RCP<block_crs_matrix_type> bcrsmatrix =
    rcp_const_cast<block_crs_matrix_type> (tif_utest::create_triangular_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,false> (crsgraph, blockSize));

  RCP<const block_crs_matrix_type> const_bcrsmatrix(bcrsmatrix);
  prec_type prec (const_bcrsmatrix);

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

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(RBILUK, FullLocalBlockCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> block_crs_matrix_type;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> crs_graph_type;
  typedef Tpetra::BlockMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> BMV;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  typedef Ifpack2::Experimental::RBILUK<block_crs_matrix_type> prec_type;

  out << "Ifpack2::Version(): " << Ifpack2::Version () << std::endl;

  const int num_rows_per_proc = 3;
  const int blockSize = 5;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  RCP<crs_graph_type> crsgraph =
    tif_utest::create_dense_local_graph<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  RCP<block_crs_matrix_type> bcrsmatrix =
    rcp_const_cast<block_crs_matrix_type> (tif_utest::create_full_local_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> (crsgraph, blockSize));

  RCP<const block_crs_matrix_type> const_bcrsmatrix(bcrsmatrix);
  prec_type prec (const_bcrsmatrix);

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

  TEST_EQUALITY(x.getMap()->getLocalNumElements(), blockSize*num_rows_per_proc);
  TEST_EQUALITY(y.getMap()->getLocalNumElements(), blockSize*num_rows_per_proc);
  TEST_NOTHROW(prec.apply(x, y));

  Teuchos::Array<Scalar> exactSol(num_rows_per_proc);
  exactSol[0] = 3.0/14.0;
  exactSol[1] = -4.0/21.0;
  exactSol[2] = 2.0/7.0;

  for (int k = 0; k < num_rows_per_proc; ++k) {
    auto ylcl = yBlock.getLocalBlockHost(k, 0, Tpetra::Access::ReadOnly);
    for (int j = 0; j < blockSize; ++j) {
      TEST_FLOATING_EQUALITY(ylcl(j), exactSol[k], 1e-14);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(RBILUK, BandedBlockCrsMatrixWithDropping, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> block_crs_matrix_type;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> row_matrix_type;
  //typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> crs_graph_type; // unused
  typedef Tpetra::BlockMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> BMV;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  //typedef Ifpack2::Experimental::RBILUK<block_crs_matrix_type> prec_type; // unused
  typedef Ifpack2::RILUK<row_matrix_type> prec_crs_type;

  out << "Ifpack2::Version(): " << Ifpack2::Version () << std::endl;

  const int num_rows_per_proc = 10;
  const int blockSize = 1;

  const int lof = 0;
  const size_t rbandwidth = lof+2+2;
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph =
    tif_utest::create_banded_graph<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc, rbandwidth);
  RCP<block_crs_matrix_type> bcrsmatrix =
    rcp_const_cast<block_crs_matrix_type> (tif_utest::create_banded_block_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> (crsgraph, blockSize, rbandwidth));

  RCP<const block_crs_matrix_type> const_bcrsmatrix(bcrsmatrix);
  Ifpack2::Experimental::RBILUK<block_crs_matrix_type> prec (const_bcrsmatrix);

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

  TEST_EQUALITY(x.getMap()->getLocalNumElements(), blockSize*num_rows_per_proc);
  TEST_EQUALITY(y.getMap()->getLocalNumElements(), blockSize*num_rows_per_proc);
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

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(RBILUK, BlockMatrixOps, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef Kokkos::View<Scalar**,Kokkos::LayoutRight,Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged> > little_block_type;
  typedef Kokkos::View<Scalar*,Kokkos::LayoutRight,Kokkos::HostSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged> > little_vec_type;
  typedef typename Kokkos::ArithTraits<Scalar>::val_type impl_scalar_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;

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

  little_block_type A (aMatrix.getRawPtr (), blockSize, blockSize); // row major
  little_block_type B (bMatrix.getRawPtr (), blockSize, blockSize); // row major
  little_block_type C (cMatrix.getRawPtr (), blockSize, blockSize); // row major
  little_block_type I (identityMatrix.getRawPtr (), blockSize, blockSize); // row major

  Tpetra::GEMM ("N", "N", STS::one (), A, I, STS::zero (), C);
  //blockOps.square_matrix_matrix_multiply(aMatrix.getRawPtr(), identityMatrix.getRawPtr(), cMatrix.getRawPtr(), blockSize);

  for (int k = 0; k < blockMatSize; ++k)
  {
    TEST_FLOATING_EQUALITY(aMatrix[k], cMatrix[k], 1e-14);
  }

  Tpetra::GEMM ("N", "N", -STS::one (), A, A, STS::one (), C);
  //blockOps.square_matrix_matrix_multiply(aMatrix.getRawPtr(), aMatrix.getRawPtr(), cMatrix.getRawPtr(), blockSize, -1.0, 1.0);

  Tpetra::GEMM ("N", "N", STS::one (), I, I, STS::one (), C);
  //blockOps.square_matrix_matrix_multiply(identityMatrix.getRawPtr(), identityMatrix.getRawPtr(), cMatrix.getRawPtr(), blockSize, 1.0, 1.0);

  Tpetra::GEMM ("N", "N", STS::one (), A, B, STS::one (), C);
  //blockOps.square_matrix_matrix_multiply(aMatrix.getRawPtr(), bMatrix.getRawPtr(), cMatrix.getRawPtr(), blockSize, 1.0, 1.0);

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

  const LocalOrdinal rowStride = blockSize;

  little_block_type dMat(cMatrix.getRawPtr(),blockSize,rowStride);
  Teuchos::Array<int> ipiv_teuchos(blockSize);
  Kokkos::View<int*, Kokkos::HostSpace,
    Kokkos::MemoryUnmanaged> ipiv (ipiv_teuchos.getRawPtr (), blockSize);

  Teuchos::Array<Scalar> work_teuchos (blockSize);
  Kokkos::View<impl_scalar_type*, Kokkos::HostSpace,
    Kokkos::MemoryUnmanaged> work (work_teuchos.getRawPtr (), blockSize);

  int lapackInfo;
  for (int k = 0; k < blockSize; ++k) {
    ipiv[k] = 0;
  }

  Tpetra::GETF2 (dMat, ipiv, lapackInfo);
  TEUCHOS_TEST_FOR_EXCEPTION(
    lapackInfo != 0, std::runtime_error, "Ifpack2::Experimental::RBILUK::compute: "
    "lapackInfo = " << lapackInfo << " which indicates an error in the factorization GETF2.");

  Tpetra::GETRI (dMat, ipiv, work, lapackInfo);
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

  little_vec_type bval(onevec.getRawPtr(),blockSize);
  little_vec_type xval(computeSolution.getRawPtr(),blockSize);

  //xval.matvecUpdate(1.0,dMat,bval);
  Tpetra::GEMV (1.0, dMat, bval, xval);

  for (int i = 0; i < blockSize; ++i)
    TEST_FLOATING_EQUALITY(exactSolution[i], computeSolution[i], 1e-13);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(RBILUK, DiagonalBlockCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal)
{
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::endl;
  typedef Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> block_crs_matrix_type;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> crs_graph_type;
  typedef Tpetra::BlockMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> BMV;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
  //typedef Ifpack2::Experimental::RBILUK<block_crs_matrix_type> prec_type; // unused
  std::ostringstream errStrm; // for error collection
  int lclSuccess = 1;
  int gblSuccess = 1;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();

  const bool useOut = false;
  RCP<Teuchos::oblackholestream> blackHole (new Teuchos::oblackholestream ());
  RCP<Teuchos::FancyOStream> newOutPtr;
  if (useOut) {
    newOutPtr = Teuchos::rcpFromRef (out);
  }
  else {
    newOutPtr = comm->getRank () == 0 ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::getFancyOStream (blackHole);
  }
  Teuchos::FancyOStream& newOut = *newOutPtr;

  newOut << "Ifpack2::RBILUK diagonal block matrix test" << endl;
  Teuchos::OSTab tab1 (out);

  const int num_rows_per_proc = 5;
  const int blockSize = 3;
  RCP<crs_graph_type> crsgraph =
    tif_utest::create_diagonal_graph<LocalOrdinal,GlobalOrdinal,Node> (num_rows_per_proc);

  RCP<block_crs_matrix_type> bcrsmatrix;
  try {
    bcrsmatrix = rcp_const_cast<block_crs_matrix_type> (tif_utest::create_block_diagonal_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> (crsgraph, blockSize));
  }
  catch (std::exception& e) {
    errStrm << "Proc " << comm->getRank () << ": BlockCrsMatrix constructor "
      "threw an exception: " << e.what () << endl;
    success = false;
  }
  catch (...) {
    errStrm << "Proc " << comm->getRank () << ": BlockCrsMatrix constructor "
      "threw an exception, not a subclass of std::exception" << endl;
    success = false;
  }
  lclSuccess = success ? 1 : 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  if (gblSuccess == 1) {
    newOut << "BlockCrsMatrix constructor did not throw an exception" << endl;
  }
  else {
    newOut << "BlockCrsMatrix constructor threw an exception on one or more "
      "processes!" << endl;
    Teuchos::OSTab tab2 (newOut);
    Tpetra::Details::gathervPrint (newOut, errStrm.str (), *comm);
    return; // stop the test on failure
  }

  RCP<const block_crs_matrix_type> const_bcrsmatrix(bcrsmatrix);
  RCP<Ifpack2::Experimental::RBILUK<block_crs_matrix_type> > prec;
  try {
    using Teuchos::rcp;
    prec = rcp (new Ifpack2::Experimental::RBILUK<block_crs_matrix_type> (const_bcrsmatrix));
  }
  catch (std::exception& e) {
    errStrm << "Proc " << comm->getRank () << ": RBILUK constructor threw an "
      "exception: " << e.what () << endl;
    success = false;
  }
  catch (...) {
    errStrm << "Proc " << comm->getRank () << ": RBILUK constructor threw an "
      "exception, not a subclass of std::exception" << endl;
    success = false;
  }
  lclSuccess = success ? 1 : 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  if (gblSuccess == 1) {
    newOut << "RBILUK constructor did not throw an exception" << endl;
  }
  else {
    newOut << "RBILUK constructor threw an exception on one or more processes!"
        << endl;
    Teuchos::OSTab tab2 (newOut);
    Tpetra::Details::gathervPrint (newOut, errStrm.str (), *comm);
    return; // stop the test on failure
  }

  Teuchos::ParameterList params;
  params.set("fact: iluk level-of-fill", (LocalOrdinal) 0);
  params.set("fact: relax value", 0.0);
  prec->setParameters(params);

  try {
    prec->initialize ();
  }
  catch (std::exception& e) {
    errStrm << "Proc " << comm->getRank () << ": prec->initialize() threw an "
      "exception: " << e.what () << endl;
    success = false;
  }
  catch (...) {
    errStrm << "Proc " << comm->getRank () << ": prec->initialize() threw an "
      "exception, not a subclass of std::exception" << endl;
    success = false;
  }
  lclSuccess = success ? 1 : 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  if (gblSuccess == 1) {
    newOut << "prec->initialize() did not throw an exception" << endl;
  }
  else {
    newOut << "prec->initialize() threw an exception on one or more processes!" << endl;
    Teuchos::OSTab tab2 (newOut);
    Tpetra::Details::gathervPrint (newOut, errStrm.str (), *comm);
    return; // stop the test on failure
  }

  try {
    prec->compute ();
  }
  catch (std::exception& e) {
    errStrm << "Proc " << comm->getRank () << ": prec->compute() threw an "
      "exception: " << e.what () << endl;
    success = false;
  }
  catch (...) {
    errStrm << "Proc " << comm->getRank () << ": prec->compute() threw an "
      "exception, not a subclass of std::exception" << endl;
    success = false;
  }

  lclSuccess = success ? 1 : 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  if (gblSuccess == 1) {
    newOut << "prec->compute() did not throw an exception" << endl;
  }
  else {
    newOut << "prec->compute() threw an exception on one or more processes!" << endl;
    Teuchos::OSTab tab2 (newOut);
    Tpetra::Details::gathervPrint (newOut, errStrm.str (), *comm);
    return; // stop the test on failure
  }

  BMV xBlock (*crsgraph->getRowMap (), blockSize, 1);
  BMV yBlock (*crsgraph->getRowMap (), blockSize, 1);
  MV x = xBlock.getMultiVectorView ();
  MV y = yBlock.getMultiVectorView ();
  x.putScalar (Teuchos::ScalarTraits<Scalar>::one ());
  prec->apply (x, y);

  const Scalar exactSol = 0.2;

  for (int k = 0; k < num_rows_per_proc; ++k) {
    auto ylcl = yBlock.getLocalBlockHost(k, 0, Tpetra::Access::ReadOnly);
    for (int j = 0; j < blockSize; ++j) {
      TEST_FLOATING_EQUALITY(ylcl(j), exactSol, 1e-14);
    }
  }
}


TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(RBILUK, AdditiveSchwarzSubdomainSolver, Scalar, LO, GO)
{
  using std::endl;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;

  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node> row_matrix_type;
  typedef Tpetra::BlockCrsMatrix<Scalar,LO,GO,Node> block_crs_matrix_type;

  out << "Test RBILUK as AdditiveSchwarz subdomain solver" << endl;
  Teuchos::OSTab tab1 (out);
  out << "Creating row Map and constant block CrsMatrix" << endl;

  const int num_rows_per_proc = 10;
  const int blockSize = 1;

  const int lof = 0;
  const size_t rbandwidth = lof+2+2;
  RCP<Tpetra::CrsGraph<LO,GO,Node> > crsgraph = tif_utest::create_banded_graph<LO,GO,Node>(num_rows_per_proc, rbandwidth);
  RCP<block_crs_matrix_type> bcrsmatrix =
    rcp_const_cast<block_crs_matrix_type> (tif_utest::create_banded_block_matrix<Scalar,LO,GO,Node> (crsgraph, blockSize, rbandwidth));

  RCP<const block_crs_matrix_type> const_bcrsmatrix(bcrsmatrix);

  int overlapLimit=1;

  for (int overlapLevel=0; overlapLevel<overlapLimit; ++overlapLevel) {

    out << "overlap = " << overlapLevel << std::endl;
    out << "Creating AdditiveSchwarz instance" << endl;
    Ifpack2::AdditiveSchwarz<row_matrix_type> prec(const_bcrsmatrix);

    out << "Filling in ParameterList for AdditiveSchwarz" << endl;
    Teuchos::ParameterList params;
    params.set ("inner preconditioner name", "RBILUK");
    params.set ("schwarz: overlap level", overlapLevel);
#   if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
    params.set ("schwarz: use reordering", true);
#   else
    params.set ("schwarz: use reordering", false);
#   endif

    out << "Setting AdditiveSchwarz's parameters" << endl;
    TEST_NOTHROW(prec.setParameters(params));

    out << "Testing domain and range Maps of AdditiveSchwarz" << endl;
    //trivial tests to insist that the preconditioner's domain/range maps are
    //identically those of the matrix:
    const Tpetra::Map<LO,GO,Node>* mtx_dom_map_ptr = &*bcrsmatrix->getDomainMap();
    const Tpetra::Map<LO,GO,Node>* mtx_rng_map_ptr = &*bcrsmatrix->getRangeMap();
    const Tpetra::Map<LO,GO,Node>* prec_dom_map_ptr = &*prec.getDomainMap();
    const Tpetra::Map<LO,GO,Node>* prec_rng_map_ptr = &*prec.getRangeMap();
    TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr );
    TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr );

    out << "Calling AdditiveSchwarz's initialize()" << endl;
    prec.initialize();

    out << "Calling AdditiveSchwarz's compute()" << endl;
    prec.compute();

    typedef Tpetra::BlockMultiVector<Scalar,LO,GO,Node> BMV;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;

    BMV xBlock (*crsgraph->getRowMap (), blockSize, 1);
    BMV yBlock (*crsgraph->getRowMap (), blockSize, 1);
    BMV zBlock (*crsgraph->getRowMap (), blockSize, 1);
    MV x = xBlock.getMultiVectorView ();
    MV y = yBlock.getMultiVectorView ();
    //MV z = zBlock.getMultiVectorView ();
    x.randomize();

    //Tpetra::MultiVector<Scalar,LO,GO,Node> x(rowmap,2), y(rowmap,2), z(rowmap,2);
    //x.putScalar(1);

    out << "Applying AdditiveSchwarz to a multivector" << endl;
    prec.apply (x, y);

    /*
    out << "Testing result of AdditiveSchwarz's apply" << endl;

    // The solution should now be full of 1/2s
    z.putScalar(0.5);

    Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();
    Teuchos::ArrayRCP<const Scalar> zview = z.get1dView();

    TEST_COMPARE_FLOATING_ARRAYS(yview, zview, 4*Teuchos::ScalarTraits<Scalar>::eps());
    */
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(RBILUK, AdditiveSchwarzReordering, Scalar, LocalOrdinal, GlobalOrdinal)
{
  using row_matrix_t = Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  using block_crs_matrix_type = Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;

  out << "Use AdditiveSchwarz to reorder subdomain for RBILUK" << std::endl;

  const int num_rows_per_proc = 5;
  const int blockSize = 3;
  auto crsgraph = tif_utest::create_diagonal_graph<LocalOrdinal,GlobalOrdinal,Node> (num_rows_per_proc);
  auto bcrsmatrix = tif_utest::create_block_diagonal_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> (crsgraph, blockSize);
  const auto & rowmap = *bcrsmatrix->getRowMap();

  // Apply some arbitrary reordering
  Teuchos::ArrayRCP<LocalOrdinal> reversePermuations(num_rows_per_proc);
  reversePermuations[0] = 4;
  reversePermuations[1] = 3;
  reversePermuations[2] = 1;
  reversePermuations[3] = 0;
  reversePermuations[4] = 2;

  // And build the inverse relationship
  Teuchos::ArrayRCP<LocalOrdinal> permuations(num_rows_per_proc);
  for (int i = 0; i < num_rows_per_proc; ++i)
    permuations[reversePermuations[i]] = i;
  
  // Scale each row to make sure ordering is correct
  for (int i = 0; i < num_rows_per_proc; ++i)
  {
    using inds_t = typename block_crs_matrix_type::local_inds_host_view_type;
    using vals_t = typename block_crs_matrix_type::nonconst_values_host_view_type;
    inds_t inds;
    vals_t vals;

    const Scalar scaleFact = 1.0 + i/Scalar(num_rows_per_proc);
    bcrsmatrix->getLocalRowViewNonConst(i,inds,vals);
    for (size_t j=0;j<vals.size();++j)
      vals[j] /= scaleFact;
  }

  // Now apply reordering with AdditiveSchwarz
  Ifpack2::AdditiveSchwarz<row_matrix_t> blockPrec(bcrsmatrix.getConst());
  {
    Teuchos::ParameterList params;
    params.set ("schwarz: overlap level", static_cast<int> (0));
    params.set ("schwarz: combine mode", "add");
    params.set ("inner preconditioner name", "RBILUK");
    params.set ("schwarz: zero starting solution", true);
    params.set ("schwarz: num iterations", 1);
    {
      Teuchos::ParameterList innerParams;
      innerParams.set ("fact: iluk level-of-fill", static_cast<int> (0));
      innerParams.set ("fact: iluk level-of-overlap", static_cast<int> (0));

      params.set ("inner preconditioner parameters", innerParams);
    }
#   if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
    params.set ("schwarz: use reordering", true);
    out << "Using reordering" << std::endl;
#   else
    params.set ("schwarz: use reordering", false);
    out << "Not using reordering, XPETRA and/or ZOLTAN2 missing." << std::endl;
#   endif
    {
      Teuchos::ParameterList zlist;
      zlist.set ("order_method", "user");
      zlist.set ("order_method_type", "local");
      zlist.set ("user ordering", permuations);
      zlist.set ("user reverse ordering", reversePermuations);

      params.set ("schwarz: reordering list", zlist);
    }
    blockPrec.setParameters(params);
  }

  blockPrec.initialize();
  blockPrec.compute();

  using BMV = Tpetra::BlockMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  using MV = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  BMV xB(rowmap,blockSize,2), yB(rowmap,blockSize,2);
  MV x = xB.getMultiVectorView ();
  MV y = yB.getMultiVectorView ();

  x.putScalar(1);
  blockPrec.apply(x, y);

  const Scalar exactSol = 0.2;
  for (int k = 0; k < num_rows_per_proc; ++k) {
    const Scalar scaledSol = exactSol*(1.0 + k/Scalar(num_rows_per_proc));
    auto ylcl = yB.getLocalBlockHost(k, 0, Tpetra::Access::ReadOnly);
    for (int j = 0; j < blockSize; ++j) {
      TEST_FLOATING_EQUALITY(ylcl(j), scaledSol, 1e-14);
    }
  }
}


# define UNIT_TEST_GROUP_BLOCK_LGN( Scalar, LocalOrdinal, GlobalOrdinal ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( RBILUK, LowerTriangularBlockCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( RBILUK, UpperTriangularBlockCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( RBILUK, BandedBlockCrsMatrixWithDropping, Scalar, LocalOrdinal, GlobalOrdinal) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( RBILUK, BlockMatrixOps, Scalar, LocalOrdinal, GlobalOrdinal) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( RBILUK, DiagonalBlockCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( RBILUK, FullLocalBlockCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( RBILUK, AdditiveSchwarzReordering, Scalar, LocalOrdinal,GlobalOrdinal) 
    // TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( RBILUK, AdditiveSchwarzSubdomainSolver, Scalar, LocalOrdinal, GlobalOrdinal)

typedef Tpetra::MultiVector<>::scalar_type scalar_type;
typedef Tpetra::MultiVector<>::local_ordinal_type local_ordinal_type;
typedef Tpetra::MultiVector<>::global_ordinal_type global_ordinal_type;

UNIT_TEST_GROUP_BLOCK_LGN(scalar_type, local_ordinal_type, global_ordinal_type)

} // namespace (anonymous)

