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
#include "BelosTpetraAdapter.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "BelosSolverFactory.hpp"
#include "TpetraCore_ETIHelperMacros.h"

// mfh 13 Oct 2017: Test Belos::CustomSolverFactory, for the Tpetra
// specialization of Belos::SolverFactory.
//
// 1. Write a custom Belos::SolverManager subclass.  It doesn't need
//    to do anything (we're not actually trying to solve with it; we
//    just want to test that we created an instance of it).
//
// 2. Write a custom Belos::CustomSolverFactory subclass, that can
//    create instances of the above solver.
//
// 3. Add our custom factory to Belos::SolverFactory.
//
// 4. Verify that Belos::SolverFactory can create instances of our
//    custom solver.

namespace { // (anonymous)

// Trivial Belos::SolverManager subclass for testing
// Belos::SolverFactory::addFactory.  It doesn't need to do anything,
// but it must be possible to construct and destruct instances of it.
template<class SC, class MV, class OP>
class FooSolver : public Belos::SolverManager<SC, MV, OP>
{
public:
  virtual ~FooSolver () {}

  Teuchos::RCP<Belos::SolverManager<SC, MV, OP> > clone () const override {
    return Teuchos::rcp(new FooSolver<SC, MV, OP>);
  }

  virtual const Belos::LinearProblem<SC, MV, OP>&
  getProblem () const {
    return linearProblem_;
  }

  virtual Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters () const {
    return Teuchos::null;
  }

  virtual Teuchos::RCP<const Teuchos::ParameterList>
  getCurrentParameters () const {
    return Teuchos::null;
  }

  virtual int getNumIters () const {
    return 0;
  }

  virtual bool isLOADetected () const {
    return false;
  }

  virtual void
  setProblem (const Teuchos::RCP<Belos::LinearProblem<SC, MV, OP> >& /* problem */)
  {}

  virtual void
  setParameters (const Teuchos::RCP<Teuchos::ParameterList>& /* params */)
  {}

  virtual void
  reset (const Belos::ResetType /* type */)
  {}

  virtual Belos::ReturnType solve () {
    return Belos::Unconverged;
  }

private:
  // Empty and trivial linear problem.  Exists only to avoid warnings
  // in getProblem().
  Belos::LinearProblem<SC, MV, OP> linearProblem_;
};

// Belos::CustomSolverFactory subclass for testing
// Belos::SolverFactory::addFactory.  This subclass knows how to
// create FooSolver instances (see above).
template<class SC, class MV, class OP>
class FooSolverFactory : public Belos::CustomSolverFactory<SC, MV, OP>
{
  virtual ~FooSolverFactory () {}

  virtual Teuchos::RCP<Belos::SolverManager<SC, MV, OP> >
  getSolver (const std::string& solverName,
             const Teuchos::RCP<Teuchos::ParameterList>& solverParams)
  {
    using Teuchos::rcp;
    using Teuchos::rcp_implicit_cast;
    typedef Belos::SolverManager<SC, MV, OP> base_solver_type;
    typedef FooSolver<SC, MV, OP> solver_type;

    if (solverName == "FOO") {
      return rcp_implicit_cast<base_solver_type> (rcp (new solver_type));
    }
    else {
      return Teuchos::null;
    }
  }

  virtual int numSupportedSolvers () const {
    return 1;
  }

  virtual std::vector<std::string> supportedSolverNames () const {
    return {"FOO"};
  }

  virtual bool isSupported (const std::string& solverName) const {
    return solverName == "FOO";
  }
};

//
// The actual unit test.
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CustomSolverFactory, AddFactory, SC, LO, GO, NT )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;
  typedef Tpetra::MultiVector<SC,LO,GO,NT> MV;
  typedef Tpetra::Operator<SC,LO,GO,NT> OP;
  typedef Belos::SolverManager<SC, MV, OP> solver_type;
  typedef Belos::CustomSolverFactory<SC, MV, OP> custom_factory_type;

  Belos::SolverFactory<SC, MV, OP> factory;

  // Sanity check for factory; can it create a GMRES solver?
  RCP<solver_type> solver = factory.create ("GMRES", Teuchos::null);
  TEST_ASSERT( ! solver.is_null () );
  solver = Teuchos::null;

  RCP<Teuchos::ParameterList> params (new Teuchos::ParameterList);
  solver = factory.create ("GMRES", params);
  TEST_ASSERT( ! solver.is_null () );
  solver = Teuchos::null;

  // At this point, 'factory' should NOT support a solver called "FOO".
  TEST_ASSERT( ! factory.isSupported ("FOO") );

  // At this point, 'factory' should NOT support a solver called "BAZ".
  TEST_ASSERT( ! factory.isSupported ("BAZ") );

  RCP<custom_factory_type> customFactory =
    rcp (static_cast<custom_factory_type*> (new FooSolverFactory<SC, MV, OP>));
  // Add an instance of our custom factory to the main factory.
  factory.addFactory (customFactory);

  // At this point, 'factory' should support a solver called "FOO".
  TEST_ASSERT( factory.isSupported ("FOO") );

  // At this point, 'factory' should NOT support a solver called "BAZ".
  TEST_ASSERT( ! factory.isSupported ("BAZ") );

  // Attempt to create an instance of our custom solver from 'factory'.
  solver = factory.create ("FOO", Teuchos::null);
  TEST_ASSERT( ! solver.is_null () );
  RCP<FooSolver<SC, MV, OP> > fooSolver =
    Teuchos::rcp_dynamic_cast<FooSolver<SC, MV, OP> > (solver);
  TEST_ASSERT( ! fooSolver.is_null () );
  fooSolver = Teuchos::null;

  // Create a new Belos::SolverFactory instance.  The above addFactory
  // call should have changed _all_ Belos::SolverFactory instances.
  Belos::SolverFactory<SC, MV, OP> factory2;
  TEST_ASSERT( factory2.isSupported ("FOO") );
  solver = factory2.create ("FOO", Teuchos::null);
  TEST_ASSERT( ! solver.is_null () );
  fooSolver = Teuchos::rcp_dynamic_cast<FooSolver<SC, MV, OP> > (solver);
  TEST_ASSERT( ! fooSolver.is_null () );
}

// Define typedefs that make the Tpetra macros work.
TPETRA_ETI_MANGLING_TYPEDEFS()

// Macro that instantiates the unit test
#define LCLINST( SC, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CustomSolverFactory, AddFactory, SC, LO, GO, NT )

// Tpetra's ETI will instantiate the unit test for all enabled type
// combinations.
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( LCLINST )

} // namespace (anonymous)
