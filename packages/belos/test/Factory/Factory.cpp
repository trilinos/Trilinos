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

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosSolverFactory.hpp"

#include "BelosGCRODRSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosFixedPointSolMgr.hpp"
#include "BelosLSQRSolMgr.hpp"
#include "BelosPCPGSolMgr.hpp"
#include "BelosRCGSolMgr.hpp"

#include "MyMultiVec.hpp"
#include "MyBetterOperator.hpp"
#include "MyOperator.hpp"


template<class ScalarType, class FactoryType, class SolverBaseType, class SolverImplType>
void
testSolver (bool& success, Teuchos::FancyOStream& out, const std::string& solverName)
{
  using Teuchos::ParameterEntry;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::TypeNameTraits;
  using std::endl;
  typedef ScalarType ST;

  typedef SolverBaseType solver_base_type;
  typedef SolverImplType solver_impl_type;
  typedef FactoryType factory_type;
  typedef typename Teuchos::ScalarTraits<ST>::magnitudeType MT;

  Teuchos::OSTab tab0 (out);
  out << "Test Belos::SolverFactory::create for solver \"" << solverName << "\"" << endl;
  Teuchos::OSTab tab1 (out);

  out << "ScalarType: " << TypeNameTraits<ScalarType>::name () << endl
      << "FactoryType: " << TypeNameTraits<FactoryType>::name () << endl
      << "SolverBaseType: " << TypeNameTraits<SolverBaseType>::name () << endl
      << "SolverImplType: " << TypeNameTraits<SolverImplType>::name () << endl;

  factory_type factory;
  RCP<solver_base_type> solver;

  out << "Test whether factory works when input ParameterList is null" << endl;

  // It must work when the parameter list is null.
  TEST_NOTHROW( solver = factory.create (solverName, Teuchos::null) );
  TEST_ASSERT( ! solver.is_null () );
  if (! solver.is_null ()) {
    // Did we actually get the solver for which we asked?
    RCP<solver_impl_type> solverImpl = rcp_dynamic_cast<solver_impl_type> (solver);
    TEST_ASSERT( ! solverImpl.is_null () );
  }

  out << "Test whether factory works when input ParameterList is nonnull" << endl;

  // The factory must work when the parameter list is nonnull, and the
  // solver must actually read the provided parameters.
  RCP<ParameterList> plist = parameterList ("Belos");
  const MT tol = 0.99; // definitely a nondefault value
  const int maxNumIters = 42; // definitely a nondefault value

  // Both of these parameters start out as "unused."
  plist->set ("Convergence Tolerance", tol);
  plist->set ("Maximum Iterations", maxNumIters);
  TEST_ASSERT( ! plist->getEntry ("Convergence Tolerance").isUsed () );
  TEST_ASSERT( ! plist->getEntry ("Maximum Iterations").isUsed () );

  out << "Input ParameterList: " << endl;
  {
    Teuchos::OSTab tab2 (out);
    plist->print (out);
    out << endl;
  }

  solver = factory.create (solverName, plist);
  TEST_ASSERT( ! solver.is_null () );
  if (! solver.is_null ()) {
    // Did we actually get the solver for which we asked?
    RCP<solver_impl_type> solverImpl = rcp_dynamic_cast<solver_impl_type> (solver);
    TEST_ASSERT( ! solverImpl.is_null () );

    // Did the solver get the parameters that we set on input?
    RCP<const ParameterList> curParams = solver->getCurrentParameters ();
    TEST_ASSERT( ! curParams.is_null () );
    if (! curParams.is_null ()) {
      // Are the parameters' values correct?
      MT curTol = Teuchos::ScalarTraits<MT>::zero ();
      TEST_NOTHROW( curTol = curParams->get<MT> ("Convergence Tolerance") );
      TEST_EQUALITY( curTol, tol );
      int curMaxNumIters = 0;
      TEST_NOTHROW( curMaxNumIters = curParams->get<int> ("Maximum Iterations") );
      TEST_EQUALITY( curMaxNumIters, maxNumIters );

      // Did the solver or the factory actually read ("use") the
      // parameters?  Each of these getEntry calls will throw if the
      // parameter doesn't exist.
      TEST_ASSERT( curParams->getEntry ("Convergence Tolerance").isUsed () );
      TEST_ASSERT( curParams->getEntry ("Maximum Iterations").isUsed () );
    }
  }
}


// Test that Belos::SolverFactory returns a solver of the right type,
// and that the solver (or the factory) read and respected the input
// parameters.
TEUCHOS_UNIT_TEST( Factory, Bug6383 )
{
  typedef double ST;
  typedef Belos::MultiVec<ST> MV;
  typedef Belos::Operator<ST> OP;
  typedef Belos::SolverManager<ST, MV, OP> solver_base_type;
  typedef Belos::SolverFactory<ST, MV, OP> factory_type;

  {
    typedef Belos::GCRODRSolMgr<ST, MV, OP> solver_impl_type;
    const std::string solverName ("GCRODR");
    testSolver<ST, factory_type, solver_base_type, solver_impl_type> (success, out, solverName);
  }
  {
    typedef Belos::PseudoBlockGmresSolMgr<ST, MV, OP> solver_impl_type;
    const std::string solverName ("GMRES");
    testSolver<ST, factory_type, solver_base_type, solver_impl_type> (success, out, solverName);
  }
  {
    typedef Belos::PseudoBlockCGSolMgr<ST, MV, OP> solver_impl_type;
    const std::string solverName ("CG");
    testSolver<ST, factory_type, solver_base_type, solver_impl_type> (success, out, solverName);
  }
  {
    typedef Belos::BlockGmresSolMgr<ST, MV, OP> solver_impl_type;
    const std::string solverName ("Block GMRES");
    testSolver<ST, factory_type, solver_base_type, solver_impl_type> (success, out, solverName);
  }
  {
    typedef Belos::BlockCGSolMgr<ST, MV, OP> solver_impl_type;
    const std::string solverName ("Block CG");
    testSolver<ST, factory_type, solver_base_type, solver_impl_type> (success, out, solverName);
  }
  {
    typedef Belos::FixedPointSolMgr<ST, MV, OP> solver_impl_type;
    const std::string solverName ("Fixed Point");
    testSolver<ST, factory_type, solver_base_type, solver_impl_type> (success, out, solverName);
  }
  {
    typedef Belos::LSQRSolMgr<ST, MV, OP> solver_impl_type;
    const std::string solverName ("LSQR");
    testSolver<ST, factory_type, solver_base_type, solver_impl_type> (success, out, solverName);
  }
  {
    typedef Belos::PCPGSolMgr<ST, MV, OP> solver_impl_type;
    const std::string solverName ("PCPG");
    testSolver<ST, factory_type, solver_base_type, solver_impl_type> (success, out, solverName);
  }
  {
    typedef Belos::RCGSolMgr<ST, MV, OP> solver_impl_type;
    const std::string solverName ("RCG");
    testSolver<ST, factory_type, solver_base_type, solver_impl_type> (success, out, solverName);
  }
}

