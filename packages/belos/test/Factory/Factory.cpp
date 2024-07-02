// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
#include "BelosBiCGStabSolMgr.hpp"

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

  const bool testOutputFreq = false;

  Teuchos::OSTab tab0 (out);
  out << "Test Belos::SolverFactory::create for solver \"" << solverName << "\"" << endl;
  Teuchos::OSTab tab1 (out);

  out << "ScalarType: " << TypeNameTraits<ScalarType>::name () << endl
      << "FactoryType: " << TypeNameTraits<FactoryType>::name () << endl
      << "SolverBaseType: " << TypeNameTraits<SolverBaseType>::name () << endl
      << "SolverImplType: " << TypeNameTraits<SolverImplType>::name () << endl;

  RCP<factory_type> pFactory;
  RCP<solver_base_type> solver;

  out << "Test whether creating a Belos::SolverFactory works" << endl;
  try {
    pFactory = rcp (new factory_type ());
  }
  catch (std::exception& e) {
    out << "Belos::SolverFactory constructor threw an exception: "
        << e.what () << endl;
    success = false;
    return; // doesn't make sense to continue
  }
  catch (...) {
    out << "Belos::SolverFactory constructor threw an exception not a subclass "
      "of std::exception" << endl;
    success = false;
    return; // doesn't make sense to continue
  }
  factory_type& factory = *pFactory;

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
  const int verbosity = static_cast<int> (Belos::Errors) + static_cast<int> (Belos::Warnings);
  const int outputStyle = Belos::Brief; // not default (imitates AztecOO)
  const int outputFreq = 3; // definitely not default

  // Both of these parameters start out as "unused."
  plist->set ("Convergence Tolerance", tol);
  plist->set ("Maximum Iterations", maxNumIters);
  plist->set ("Verbosity", verbosity);
  plist->set ("Output Style", outputStyle);
  if (testOutputFreq) {
    plist->set ("Output Frequency", outputFreq);
  }
  TEST_ASSERT( ! plist->getEntry ("Convergence Tolerance").isUsed () );
  TEST_ASSERT( ! plist->getEntry ("Maximum Iterations").isUsed () );
  TEST_ASSERT( ! plist->getEntry ("Verbosity").isUsed () );
  TEST_ASSERT( ! plist->getEntry ("Output Style").isUsed () );
  if (testOutputFreq) {
    TEST_ASSERT( ! plist->getEntry ("Output Frequency").isUsed () );
  }

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

    // Did the solver (or the factory) actually read the parameters that we set?
    TEST_ASSERT( plist->getEntry ("Convergence Tolerance").isUsed () );
    TEST_ASSERT( plist->getEntry ("Maximum Iterations").isUsed () );
    TEST_ASSERT( plist->getEntry ("Verbosity").isUsed () );
    TEST_ASSERT( plist->getEntry ("Output Style").isUsed () );
    if (testOutputFreq) {
      TEST_ASSERT( plist->getEntry ("Output Frequency").isUsed () );
    }

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
      int curVerbosity = 0;
      TEST_NOTHROW( curVerbosity = curParams->get<int> ("Verbosity") );
      TEST_EQUALITY( curVerbosity, verbosity );
      int curOutputStyle = 0;
      TEST_NOTHROW( curOutputStyle = curParams->get<int> ("Output Style") );
      TEST_EQUALITY( curOutputStyle, outputStyle );
      if (testOutputFreq) {
        int curOutputFreq = 0;
        TEST_NOTHROW( curOutputFreq = curParams->get<int> ("Output Frequency") );
        TEST_EQUALITY( curOutputFreq, outputFreq );
      }

      // The solver (or the factory) doesn't actually need to have
      // read ("used") the parameters in curParams.  curParams could
      // be (generally is, for Belos) a deep copy of plist.
    }
  }
}



// Unfortunately, the preprocessor doesn't let me do pasting like this:
//
// typedef Belos:: ## SOLVER_CLASS ## <ST, MV, OP> solver_impl_type;
//
// That's why I still have to do the typedef outside of the macro.
#define BELOS_TEST_SOLVER( SOLVER_NAME ) \
do { \
  const std::string solverName (SOLVER_NAME); \
  try { \
    bool curSuccess = true; \
    testSolver<ST, factory_type, solver_base_type, solver_impl_type> (curSuccess, out, solverName); \
    if (! curSuccess) { \
      failedSolvers.push_back (solverName); \
    } \
    success = success && curSuccess; \
  } catch (std::exception& e) { \
    out << "*** Solver \"" << solverName << "\" threw an exception: " \
        << e.what () << std::endl; \
    success = false; \
    failedSolvers.push_back (solverName); \
  } \
} while (false)


// Test that Belos::SolverFactory returns a solver of the right type,
// and that the solver (or the factory) read and respected the input
// parameters.
TEUCHOS_UNIT_TEST( Factory, Bug6383 )
{
  using std::endl;
  typedef double ST;
  typedef Belos::MultiVec<ST> MV;
  typedef Belos::Operator<ST> OP;
  typedef Belos::SolverManager<ST, MV, OP> solver_base_type;
  typedef Belos::SolverFactory<ST, MV, OP> factory_type;

  Teuchos::OSTab tab0 (out);
  out << "Test for Bug 6383" << endl;
  Teuchos::OSTab tab1 (out);

  // List of names of solvers that failed the test.
  std::vector<std::string> failedSolvers;

  {
    typedef Belos::GCRODRSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "GCRODR" );
  }
  {
    typedef Belos::PseudoBlockGmresSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "GMRES" );
  }
  // Make sure that the factory is case insensitive (Bug 6388).
  {
    typedef Belos::PseudoBlockGmresSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "gmres" );
  }
  {
    typedef Belos::PseudoBlockCGSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "CG" );
  }
  // Make sure that the factory is case insensitive (Bug 6388).
  {
    typedef Belos::PseudoBlockCGSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "cg" );
  }
  {
    typedef Belos::BlockGmresSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "Block GMRES" );
  }
  {
    typedef Belos::BlockGmresSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "Block GMRES" );
  }
  {
    typedef Belos::BlockCGSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "Block CG" );
  }
  {
    typedef Belos::FixedPointSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "Fixed Point" );
  }

  // FIXME (mfh 05 Aug 2015) When setting "Verbosity" and/or "Output
  // Style", LSQR throws:
  //
  // .../packages/belos/src/BelosStatusTestResNormOutput.hpp:232:
  //
  // Throw test that evaluated to true: tmpComboTest == Teuchos::null
  // StatusTestResNormOutput():  test must be Belos::StatusTest[MaxIters|ResNorm|Combo].
  if (false) {
    typedef Belos::LSQRSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "LSQR" );
  }

  {
    typedef Belos::PCPGSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "PCPG" );
  }
  {
    typedef Belos::RCGSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "RCG" );
  }
  {
    typedef Belos::BiCGStabSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "BiCGStab" );
  }
  // Make sure that the factory is case insensitive (Bug 6388).
  {
    typedef Belos::BiCGStabSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "bicgstab" );
  }

#if 1
  if (success) {
    out << endl << "Test SUCCEEDED!" << endl;
  }
  else {
    out << endl << "Test FAILED!" << endl
        << "Solvers that failed: [";
    for (size_t k = 0; k < failedSolvers.size (); ++k) {
      out << "\"" << failedSolvers[k] << "\"";
      if (k + 1 < failedSolvers.size ()) {
        out << ", ";
      }
    }
    out << "]" << endl;
  }
#else
  if (! success) {
    out << endl << "*** Solvers that failed: ";
    out << "[";
    for (size_t k = 0; k < failedSolvers.size (); ++k) {
      out << "\"" << failedSolvers[k] << "\"";
      if (k + 1 < failedSolvers.size ()) {
        out << ", ";
      }
    }
    out << "]" << endl;
  }
#endif // 0
}

// Repeat the above test for the float scalar type. 
TEUCHOS_UNIT_TEST( Factory, Bug6383_Float )
{
  using std::endl;
  typedef float ST;
  typedef Belos::MultiVec<ST> MV;
  typedef Belos::Operator<ST> OP;
  typedef Belos::SolverManager<ST, MV, OP> solver_base_type;
  typedef Belos::SolverFactory<ST, MV, OP> factory_type;

  Teuchos::OSTab tab0 (out);
  out << "Test for Bug 6383" << endl;
  Teuchos::OSTab tab1 (out);

  // List of names of solvers that failed the test.
  std::vector<std::string> failedSolvers;

  {
    typedef Belos::GCRODRSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "GCRODR" );
  }
  {
    typedef Belos::PseudoBlockGmresSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "GMRES" );
  }
  // Make sure that the factory is case insensitive (Bug 6388).
  {
    typedef Belos::PseudoBlockGmresSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "gmres" );
  }
  {
    typedef Belos::PseudoBlockCGSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "CG" );
  }
  // Make sure that the factory is case insensitive (Bug 6388).
  {
    typedef Belos::PseudoBlockCGSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "cg" );
  }
  {
    typedef Belos::BlockGmresSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "Block GMRES" );
  }
  {
    typedef Belos::BlockGmresSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "Block GMRES" );
  }
  {
    typedef Belos::BlockCGSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "Block CG" );
  }
  {
    typedef Belos::FixedPointSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "Fixed Point" );
  }

  // FIXME (mfh 05 Aug 2015) When setting "Verbosity" and/or "Output
  // Style", LSQR throws:
  //
  // .../packages/belos/src/BelosStatusTestResNormOutput.hpp:232:
  //
  // Throw test that evaluated to true: tmpComboTest == Teuchos::null
  // StatusTestResNormOutput():  test must be Belos::StatusTest[MaxIters|ResNorm|Combo].
  if (false) {
    typedef Belos::LSQRSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "LSQR" );
  }

  {
    typedef Belos::PCPGSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "PCPG" );
  }
  {
    typedef Belos::RCGSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "RCG" );
  }
  {
    typedef Belos::BiCGStabSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "BiCGStab" );
  }
  // Make sure that the factory is case insensitive (Bug 6388).
  {
    typedef Belos::BiCGStabSolMgr<ST, MV, OP> solver_impl_type;
    BELOS_TEST_SOLVER( "bicgstab" );
  }

#if 1
  if (success) {
    out << endl << "Test SUCCEEDED!" << endl;
  }
  else {
    out << endl << "Test FAILED!" << endl
        << "Solvers that failed: [";
    for (size_t k = 0; k < failedSolvers.size (); ++k) {
      out << "\"" << failedSolvers[k] << "\"";
      if (k + 1 < failedSolvers.size ()) {
        out << ", ";
      }
    }
    out << "]" << endl;
  }
#else
  if (! success) {
    out << endl << "*** Solvers that failed: ";
    out << "[";
    for (size_t k = 0; k < failedSolvers.size (); ++k) {
      out << "\"" << failedSolvers[k] << "\"";
      if (k + 1 < failedSolvers.size ()) {
        out << ", ";
      }
    }
    out << "]" << endl;
  }
#endif // 0
}
