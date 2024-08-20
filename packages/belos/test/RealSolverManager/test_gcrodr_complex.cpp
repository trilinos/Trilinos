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
#include <BelosGCRODRSolMgr.hpp>
#include <complex>

#include "MyMultiVec.hpp"
#include "MyBetterOperator.hpp"
#include "MyOperator.hpp"

//
// mfh 20 Jan 2014: This test ensures the following:
//
// 1. Belos::GCRODRSolMgr can compile whether its ScalarType (first)
//    template parameter is real or complex.
// 2. Belos::GCRODRSolMgr's constructor throws std::logic_error if and
//    only if its ScalarType (first) template parameter is complex.
//
// At some point, if GCRODRSolMgr gets fixed so that it works with
// complex ScalarType, the second test will no longer pass.  This will
// be a good thing!  The test should still be built in that case, in
// order to demonstrate that GCRODRSolMgr compiles for complex
// ScalarType.  However, in that case, the TEST_THROW macro should be
// changed to TEST_NOTHROW, and the macro's second argument should be
// removed.
//
// This test requires that Trilinos was compiled with complex
// arithmetic support enabled.
//
// Updates
// 4/15/2014: GCRODR is changed from RealSolverManager to SolverManager,
//            and the test is changed from THROW to NOTHROW. (phtsuji)

TEUCHOS_UNIT_TEST( GCRODR, RealDoesNotThrow )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef double ST;
  typedef Belos::MultiVec<ST> MV;
  typedef Belos::Operator<ST> OP;
  //typedef Belos::MultiVecTraits<ST, MV> MVT;
  //typedef Belos::OperatorTraits<ST, MV, OP> OPT;
  typedef Belos::GCRODRSolMgr<ST, MV, OP> sol_mgr_type;

  RCP<sol_mgr_type> solver;
  TEST_NOTHROW( solver = rcp (new sol_mgr_type ()) );
}

TEUCHOS_UNIT_TEST( GCRODR, ComplexDoesNotThrow )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef std::complex<double> ST;
  typedef Belos::MultiVec<ST> MV;
  typedef Belos::Operator<ST> OP;
  // typedef Belos::MultiVecTraits<ST, MV> MVT;
  // typedef Belos::OperatorTraits<ST, MV, OP> OPT;
  typedef Belos::GCRODRSolMgr<ST, MV, OP> sol_mgr_type;

  RCP<sol_mgr_type> solver;
  TEST_NOTHROW( solver = rcp (new sol_mgr_type ()) );
}
