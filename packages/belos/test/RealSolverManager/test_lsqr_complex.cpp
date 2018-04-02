//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
//@HEADER

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
