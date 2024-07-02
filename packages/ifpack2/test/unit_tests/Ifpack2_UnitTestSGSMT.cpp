// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Ifpack2_UnitTestSGSMT.cpp
/// \brief Unit test for multithreaded symmetric Gauss-Seidel

#include "Teuchos_UnitTestHarness.hpp"
#include "Ifpack2_UnitTestHelpers.hpp"
#include "Ifpack2_Relaxation.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include <KokkosSparse_gauss_seidel.hpp>

namespace { // (anonymous)

using Teuchos::RCP;
typedef tif_utest::Node Node;
using std::endl;

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(SGS_MT, JacobiComparison, Scalar, LO, GO)
{
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::Comm;
  using Teuchos::ParameterList;
  using std::endl;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> crs_matrix_type;
  typedef Tpetra::Map<LO, GO, Node> map_type;
  typedef Tpetra::Vector<Scalar, LO, GO, Node> vec_type;
  typedef Tpetra::RowMatrix<Scalar, LO, GO, Node> row_matrix_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  Teuchos::OSTab tab0 (out);
  out << "Test multiple threaded SGS sweeps compared w.r.t. Jacobi" << endl;
  Teuchos::OSTab tab1 (out);

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();

  Tpetra::MatrixMarket::Reader<crs_matrix_type> crs_reader;
  std::string file_name = "sherman1.mtx";

  RCP<const crs_matrix_type> A = crs_reader.readSparseFile(file_name, comm);

  RCP<const map_type> rowMap = A->getRowMap();
  RCP<const map_type> domainMap = rowMap;
  RCP<const map_type> rangeMap = rowMap;


  vec_type X_wanted (domainMap);
  vec_type Y_result (rangeMap);

  X_wanted.randomize ();
  A->apply(X_wanted, Y_result);


  // Test Gauss-Seidel (SGS) with ten sweeps.
  // Start by letting SGS set the starting solution to zero.
  ParameterList params_mt_sgs;

  params_mt_sgs.set ("relaxation: type", "MT Gauss-Seidel");
  params_mt_sgs.set ("relaxation: sweeps", 10);

  params_mt_sgs.set ("relaxation: zero starting solution", true);
  Ifpack2::Relaxation<row_matrix_type> prec_mt_sgs (A);

  prec_mt_sgs.setParameters (params_mt_sgs);
  TEST_NOTHROW( prec_mt_sgs.initialize () );
  TEST_NOTHROW( prec_mt_sgs.compute () );


  vec_type X_seek (domainMap);
  X_seek.putScalar (STS::zero ());
  vec_type initial_X_diff(X_seek, Teuchos::Copy);
  initial_X_diff.update (1, X_wanted, -1);
  double initial_normInf = initial_X_diff.normInf ();


  prec_mt_sgs.apply (Y_result, X_seek);

  vec_type X_diff_mtgs(X_seek, Teuchos::Copy);
  X_diff_mtgs.update (1, X_wanted, -1);
  double normInf_mt_gs = X_diff_mtgs.normInf ();


  // Test two-stage Gauss-Seidel (SGS) with ten sweeps, and ten inner sweeps.
  ParameterList params_sgs2;
  params_sgs2.set ("relaxation: type", "Two-stage Gauss-Seidel");
  params_sgs2.set ("relaxation: sweeps", 10);
  params_sgs2.set ("relaxation: inner sweeps", 10);
  params_sgs2.set ("relaxation: zero starting solution", true);

  Ifpack2::Relaxation<row_matrix_type> prec_sgs2 (A);

  prec_sgs2.setParameters (params_sgs2);
  TEST_NOTHROW( prec_sgs2.initialize () );
  TEST_NOTHROW( prec_sgs2.compute () );

  X_seek.putScalar (STS::zero ());
  prec_sgs2.apply (Y_result, X_seek);

  vec_type X_diff_sgs2(X_seek, Teuchos::Copy);
  X_diff_sgs2.update (1, X_wanted, -1);
  double normInf_sgs2 = X_diff_sgs2.normInf ();



  // Jacobi with ten sweeps
  ParameterList params_jacobi;
  params_jacobi.set ("relaxation: type", "Jacobi");
  params_jacobi.set ("relaxation: sweeps", 10);
  params_jacobi.set ("relaxation: zero starting solution", true);
  Ifpack2::Relaxation<row_matrix_type> prec_jacobi(A);

  prec_jacobi.setParameters (params_jacobi);
  TEST_NOTHROW( prec_jacobi.initialize () );
  TEST_NOTHROW( prec_jacobi.compute () );

  X_seek.putScalar (STS::zero ());
  prec_jacobi.apply (Y_result, X_seek);

  vec_type X_diff_jacobi(X_seek, Teuchos::Copy);
  X_diff_jacobi.update (1, X_wanted, -1);
  double normInf_jacobi = X_diff_jacobi.normInf ();

  bool sgs_improve = initial_normInf > normInf_mt_gs;
  bool sgs_better_than_jacobi = normInf_mt_gs < normInf_jacobi;

  bool sgs2_improve = initial_normInf > normInf_sgs2;
  bool sgs2_better_than_jacobi = normInf_sgs2 < normInf_jacobi;

  out << "Initial Norm:" << initial_normInf
      << " MT GS Norm:" << normInf_mt_gs
      << " SGS-2 Norm:" << normInf_sgs2
      << " Jacobi Norm:" << normInf_jacobi << std::endl;

  TEST_EQUALITY( sgs_improve, 1 );
  TEST_EQUALITY( sgs_better_than_jacobi, 1 );

  TEST_EQUALITY( sgs2_improve, 1 );
  TEST_EQUALITY( sgs2_better_than_jacobi, 1 );
}


#define UNIT_TEST_GROUP_SC_LO_GO( Scalar, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( SGS_MT, JacobiComparison, Scalar, LO, GO )
#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// TODO (mfh 24 Aug 2016) Test complex Scalar types

IFPACK2_INSTANTIATE_SLG_REAL( UNIT_TEST_GROUP_SC_LO_GO )

} // namespace (anonymous)


