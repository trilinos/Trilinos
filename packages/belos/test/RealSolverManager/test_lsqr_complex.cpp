// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <BelosConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosLSQRSolMgr.hpp>
#include <complex>

#include "MyMultiVec.hpp"
#include "MyBetterOperator.hpp"
#include "MyOperator.hpp"

//
// mfh 20 Jan 2014: This test ensures the following:
//
// 1. Belos::LSQRSolMgr can compile whether its ScalarType (first)
//    template parameter is real or complex.
// 2. Belos::LSQRSolMgr's constructor throws std::logic_error if and
//    only if its ScalarType (first) template parameter is complex.
//
// At some point, if LSQRSolMgr gets fixed so that it works with
// complex ScalarType, the second test will no longer pass.  This will
// be a good thing!  The test should still be built in that case, in
// order to demonstrate that LSQRSolMgr compiles for complex
// ScalarType.  However, in that case, the TEST_THROW macro should be
// changed to TEST_NOTHROW, and the macro's second argument should be
// removed.
//
// This test requires that Trilinos was compiled with complex
// arithmetic support enabled.
//

TEUCHOS_UNIT_TEST( LSQR, RealDoesNotThrow )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef double ST;
  typedef Belos::MultiVec<ST> MV;
  typedef Belos::Operator<ST> OP;
  // typedef Belos::MultiVecTraits<ST, MV> MVT;
  // typedef Belos::OperatorTraits<ST, MV, OP> OPT;
  typedef Belos::LSQRSolMgr<ST, MV, OP> sol_mgr_type;

  RCP<sol_mgr_type> solver;
  TEST_NOTHROW( solver = rcp (new sol_mgr_type ()) );
}

TEUCHOS_UNIT_TEST( LSQR, ComplexThrows )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef std::complex<double> ST;
  typedef Belos::MultiVec<ST> MV;
  typedef Belos::Operator<ST> OP;
  // typedef Belos::MultiVecTraits<ST, MV> MVT;
  // typedef Belos::OperatorTraits<ST, MV, OP> OPT;
  typedef Belos::LSQRSolMgr<ST, MV, OP> sol_mgr_type;

  // no longer throws due to DII system needing to make dummy constructor
  RCP<sol_mgr_type> solver = rcp (new sol_mgr_type ());
  TEST_THROW( solver->getNumIters() , std::logic_error ); // throws!
}
