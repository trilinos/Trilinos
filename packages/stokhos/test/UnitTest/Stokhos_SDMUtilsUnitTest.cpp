// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos_SDMUtils.hpp"
#include "Stokhos_UnitTestHelpers.hpp"

// Need to test pre-pivot

namespace SDMUtilsUnitTest {

  typedef int ordinal_type;
  typedef double scalar_type;
  typedef Teuchos::SerialDenseMatrix<ordinal_type,scalar_type> SDM;
  typedef void (*qr_func_type)(ordinal_type, const SDM&, const Teuchos::Array<scalar_type>&, SDM&, SDM&);
  typedef void (*cpqr_func_type)(const SDM&, SDM&, SDM&, Teuchos::Array<ordinal_type>&);
  typedef ordinal_type (*wcpqr_func_type)(const scalar_type&, const SDM&, const Teuchos::Array<scalar_type>&, SDM&, SDM&, Teuchos::Array<ordinal_type>&);
  
  scalar_type rtol = 1.0e-11;
  scalar_type atol = 1.0e-11;

  void generateRandomMatrix(SDM& A, ordinal_type rank) {
    ordinal_type m = A.numRows();
    ordinal_type n = A.numCols();
    ordinal_type k = std::min(m,n);
    SDM B(m,m), C(n,n), S(m,n), T(m,n);
    B.random();
    C.random();
    S.putScalar(0.0);
    if (rank > k)
      rank = k;
    for (ordinal_type i=0; i<rank; i++)
      S(i,i) = 1.0;
    T.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, B, S, 0.0);
    A.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, T, C, 0.0);
  }

  bool test_QR(qr_func_type qr_func, ordinal_type m, ordinal_type n,
	       scalar_type rtol, scalar_type atol,
	       Teuchos::FancyOStream& out) {
    bool success;

    SDM A(m,n);
    A.random();
    SDM Q, R;
    Teuchos::Array<scalar_type> w(m, 1.0);
    ordinal_type k = std::min(m,n);
    qr_func(k, A, w, Q, R);

    TEUCHOS_TEST_EQUALITY(Q.numRows(), m, out, success);
    TEUCHOS_TEST_EQUALITY(Q.numCols(), k, out, success);
    TEUCHOS_TEST_EQUALITY(R.numRows(), k, out, success);
    TEUCHOS_TEST_EQUALITY(R.numCols(), k, out, success);

    // Test A = Q*R
    SDM QR(m,k);
    QR.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Q, R, 0.0);
    SDM AA(Teuchos::View, A, m, k);
    success = Stokhos::compareSDM(AA, "A", QR, "Q*R", rtol, atol, out);

    // Test Q^T*Q = I
    SDM eye(k,k), QTQ(k,k);
    for (ordinal_type i=0; i<k; i++)
      eye(i,i) = 1.0;
    QTQ.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Q, Q, 0.0);
    success = Stokhos::compareSDM(eye, "eye", QTQ, "Q^T*Q", rtol, atol, out);

    return success;
  }

  bool test_CPQR(cpqr_func_type qr_func, ordinal_type m, ordinal_type n,
		 scalar_type rtol, scalar_type atol,
		 Teuchos::FancyOStream& out) {
    bool success;

    SDM A(m,n);
    ordinal_type k = std::min(m,n);
    generateRandomMatrix(A, k);
    SDM Q, R;
    Teuchos::Array<ordinal_type> piv;
    qr_func(A, Q, R, piv);

    TEUCHOS_TEST_EQUALITY(Q.numRows(), m, out, success);
    TEUCHOS_TEST_EQUALITY(Q.numCols(), k, out, success);
    TEUCHOS_TEST_EQUALITY(R.numRows(), k, out, success);
    TEUCHOS_TEST_EQUALITY(R.numCols(), n, out, success);
    TEUCHOS_TEST_EQUALITY(piv.size(), n, out, success);

    // Test A*P = Q*R
    SDM AP(m,n), QR(m,n);
    for (ordinal_type j=0; j<n; j++)
      for (ordinal_type i=0; i<m; i++)
	AP(i,j) = A(i,piv[j]);
    QR.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Q, R, 0.0);
    success = Stokhos::compareSDM(AP, "A*P", QR, "Q*R", rtol, atol, out);

    // Test Q^T*Q = I
    SDM eye(k,k), QTQ(k,k);
    for (ordinal_type i=0; i<k; i++)
      eye(i,i) = 1.0;
    QTQ.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Q, Q, 0.0);
    success = Stokhos::compareSDM(eye, "eye", QTQ, "Q^T*Q", rtol, atol, out);

    return success;
  }

  bool test_weighted_CPQR(
    wcpqr_func_type qr_func, ordinal_type m, ordinal_type n, ordinal_type k,
    scalar_type rtol, scalar_type atol,
    Teuchos::FancyOStream& out) {
    bool success;

    SDM A(m,n);
    generateRandomMatrix(A, k);
    SDM Q, R;
    Teuchos::Array<ordinal_type> piv(n);
    for (ordinal_type i=0; i<5; i++)
      piv[i] = 1;
    Teuchos::Array<scalar_type> w(m, 1.0);
    ordinal_type rank = qr_func(1.0e-12, A, w, Q, R, piv);

    TEUCHOS_TEST_EQUALITY(rank, k, out, success);
    TEUCHOS_TEST_EQUALITY(Q.numRows(), m, out, success);
    TEUCHOS_TEST_EQUALITY(Q.numCols(), k, out, success);
    TEUCHOS_TEST_EQUALITY(R.numRows(), k, out, success);
    TEUCHOS_TEST_EQUALITY(R.numCols(), k, out, success);
    TEUCHOS_TEST_EQUALITY(piv.size(), n, out, success);

    // Test A*P = Q*R
    SDM AP(m,k), QR(m,k);
    for (ordinal_type j=0; j<k; j++)
      for (ordinal_type i=0; i<m; i++)
	AP(i,j) = A(i,piv[j]);
    QR.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Q, R, 0.0);
    success = Stokhos::compareSDM(AP, "A*P", QR, "Q*R", rtol, atol, out);

    // Test Q^T*Q = I
    SDM eye(k,k), Qt(m,k), QTQ(k,k);
    for (ordinal_type j=0; j<k; j++) {
      eye(j,j) = 1.0;
      for (ordinal_type i=0; i<m; i++)
	Qt(i,j) = w[i]*Q(i,j);
    }
    QTQ.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Qt, Q, 0.0);
    success = Stokhos::compareSDM(eye, "eye", QTQ, "Q^T*Q", rtol, atol, out);

    return success;
  }

  TEUCHOS_UNIT_TEST( Stokhos_SDMUtils, QR_CGS_TallSkinny ) {
    ordinal_type m = 100;
    ordinal_type n = 35;
    success = test_QR(Stokhos::QR_CGS, m, n, rtol, atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SDMUtils, QR_CGS_ShortFat ) {
    ordinal_type n = 100;
    ordinal_type m = 35;
    success = test_QR(Stokhos::QR_CGS, m, n, rtol, atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SDMUtils, QR_MGS_TallSkinny ) {
    ordinal_type m = 100;
    ordinal_type n = 35;
    success = test_QR(Stokhos::QR_MGS, m, n, rtol, atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SDMUtils, QR_MGS_ShortFat ) {
    ordinal_type n = 100;
    ordinal_type m = 35;
    success = test_QR(Stokhos::QR_MGS, m, n, rtol, atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SDMUtils, QR_MGS2_TallSkinny ) {
    ordinal_type m = 100;
    ordinal_type n = 35;
    success = test_QR(Stokhos::QR_MGS2, m, n, rtol, atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SDMUtils, QR_MGS2_ShortFat ) {
    ordinal_type n = 100;
    ordinal_type m = 35;
    success = test_QR(Stokhos::QR_MGS2, m, n, rtol, atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SDMUtils, QR_Householder_TallSkinny ) {
    ordinal_type m = 100;
    ordinal_type n = 35;
    success = test_QR(Stokhos::QR_Householder, m, n, rtol, atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SDMUtils, QR_Householder_ShortFat ) {
    ordinal_type n = 100;
    ordinal_type m = 35;
    success = test_QR(Stokhos::QR_Householder, m, n, rtol, atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SDMUtils, CPQR_Householder3_TallSkinny ) {
    ordinal_type m = 100;
    ordinal_type n = 35;
    success = test_CPQR(Stokhos::CPQR_Householder3, m, n, rtol, atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SDMUtils, CPQR_Householder3_ShortFat ) {
    ordinal_type n = 100;
    ordinal_type m = 35;
    success = test_CPQR(Stokhos::CPQR_Householder3, m, n, rtol, atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SDMUtils, CPQR_Householder_threshold_TallSkinny ) {
    ordinal_type m = 100;
    ordinal_type n = 35;
    ordinal_type k = 20;
    success = test_weighted_CPQR(Stokhos::CPQR_Householder_threshold, 
				 m, n, k, rtol, atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SDMUtils, CPQR_Householder_threshold_ShortFat ) {
    ordinal_type n = 100;
    ordinal_type m = 35;
    ordinal_type k = 20;
    success = test_weighted_CPQR(Stokhos::CPQR_Householder_threshold, 
				 m, n, k, rtol, atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SDMUtils, CPQR_MGS_threshold_TallSkinny ) {
    ordinal_type m = 100;
    ordinal_type n = 35;
    ordinal_type k = 20;
    success = test_weighted_CPQR(Stokhos::CPQR_MGS_threshold, 
				 m, n, k, rtol, atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SDMUtils, CPQR_MGS_threshold_ShortFat ) {
    ordinal_type n = 100;
    ordinal_type m = 35;
    ordinal_type k = 20;
    success = test_weighted_CPQR(Stokhos::CPQR_MGS_threshold, 
				 m, n, k, rtol, atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SDMUtils, CPQR_MGS_reorthog_threshold_TallSkinny ) {
    ordinal_type m = 100;
    ordinal_type n = 35;
    ordinal_type k = 20;
    success = test_weighted_CPQR(Stokhos::CPQR_MGS_reorthog_threshold, 
				 m, n, k, rtol, atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SDMUtils, CPQR_MGS_reorthog_threshold_ShortFat ) {
    ordinal_type n = 100;
    ordinal_type m = 35;
    ordinal_type k = 20;
    success = test_weighted_CPQR(Stokhos::CPQR_MGS_reorthog_threshold, 
				 m, n, k, rtol, atol, out);
  }

}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
