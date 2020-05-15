//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_UnitTestHarness.hpp"

// NOX Objects
#include "NOX.H"
#include "NOX_Thyra.H"

// Trilinos Objects
#include "Teuchos_Comm.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Vector.hpp"

#include "BelosTypes.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_Ifpack2PreconditionerFactory.hpp"
#include "ME_Tpetra_1DFEM.hpp"

#include "NOX_Thyra_MatrixFreeJacobianOperator.hpp"
#include "NOX_MatrixFree_ModelEvaluatorDecorator.hpp"

const Tpetra::global_size_t numGlobalElements = 100;

TEUCHOS_UNIT_TEST(NOX_Tpetra_1DFEM, AnalyticJacobian_NoPrec)
{
  Teuchos::TimeMonitor::zeroOutTimers();

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  // Get default Tpetra template types
  typedef Tpetra::Vector<>::scalar_type Scalar;
  typedef Tpetra::Vector<>::local_ordinal_type LO;
  typedef Tpetra::Vector<>::global_ordinal_type GO;
  typedef Tpetra::Vector<>::node_type Node;

  // Create the model evaluator object
  Scalar x00 = 0.0;
  Scalar x01 = 1.0;
  Teuchos::RCP<EvaluatorTpetra1DFEM<Scalar,LO,GO,Node> > model =
    evaluatorTpetra1DFEM<Scalar,LO,GO,Node>(comm, numGlobalElements, x00, x01);

  Stratimikos::DefaultLinearSolverBuilder builder;

  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
  p->set("Linear Solver Type", "Belos");
  Teuchos::ParameterList& belosList = p->sublist("Linear Solver Types").sublist("Belos");
  belosList.set("Solver Type", "Pseudo Block GMRES");
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set<int>("Maximum Iterations", 200);
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set<int>("Num Blocks", 200);
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set("Verbosity", 0x7f);
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set("Output Frequency", 100);
  belosList.sublist("VerboseObject").set("Verbosity Level", "medium");
  p->set("Preconditioner Type", "None");
  builder.setParameterList(p);

  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
    lowsFactory = builder.createLinearSolveStrategy("");

  model->set_W_factory(lowsFactory);

  // Create the initial guess
  Teuchos::RCP<Thyra::VectorBase<Scalar> >
    initial_guess = model->getNominalValues().get_x()->clone_v();
  Thyra::V_S(initial_guess.ptr(),Teuchos::ScalarTraits<Scalar>::one());

  Teuchos::RCP<NOX::Thyra::Group> nox_group =
    Teuchos::rcp(new NOX::Thyra::Group(*initial_guess, model, model->create_W_op(), lowsFactory, Teuchos::null, Teuchos::null, Teuchos::null));

  nox_group->computeF();

  // Create the NOX status tests and the solver
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
    Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  converged->addStatusTest(wrms);
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // Create nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> nl_params = Teuchos::parameterList();
  nl_params->set("Nonlinear Solver", "Line Search Based");
  nl_params->sublist("Direction").sublist("Newton").sublist("Linear Solver").set("Tolerance", 1.0e-4);

  // Set output parameters
  nl_params->sublist("Printing").sublist("Output Information").set("Debug",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Warning",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Error",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Test Details",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Details",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Parameters",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Linear Solver Details",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Inner Iteration",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Outer Iteration",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Outer Iteration StatusTest",true);

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(nox_group, combo, nl_params);
  NOX::StatusTest::StatusType solvStatus = solver->solve();

  TEST_ASSERT(solvStatus == NOX::StatusTest::Converged);

  Teuchos::TimeMonitor::summarize();
}

TEUCHOS_UNIT_TEST(NOX_Tpetra_1DFEM, AnalyticJacobian_Ifpack2Prec)
{
  Teuchos::TimeMonitor::zeroOutTimers();

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  // Get default Tpetra template types
  typedef Tpetra::Vector<>::scalar_type Scalar;
  typedef Tpetra::Vector<>::local_ordinal_type LO;
  typedef Tpetra::Vector<>::global_ordinal_type GO;
  typedef Tpetra::Vector<>::node_type Node;

  // Create the model evaluator object
  Scalar x00 = 0.0;
  Scalar x01 = 1.0;
  Teuchos::RCP<EvaluatorTpetra1DFEM<Scalar,LO,GO,Node> > model =
    evaluatorTpetra1DFEM<Scalar,LO,GO,Node>(comm, numGlobalElements, x00, x01);

  Stratimikos::DefaultLinearSolverBuilder builder;
  typedef Thyra::PreconditionerFactoryBase<Scalar> Base;
  typedef Thyra::Ifpack2PreconditionerFactory<Tpetra::CrsMatrix<Scalar, LO, GO, Node> > Impl;
  builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");

  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
  p->set("Linear Solver Type", "Belos");
  Teuchos::ParameterList& belosList = p->sublist("Linear Solver Types").sublist("Belos");
  belosList.set("Solver Type", "Pseudo Block GMRES");
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set<int>("Maximum Iterations", 200);
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set<int>("Num Blocks", 200);
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set("Verbosity", 0x7f);
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set("Output Frequency", 100);
  belosList.sublist("VerboseObject").set("Verbosity Level", "medium");
  p->set("Preconditioner Type", "Ifpack2");
  Teuchos::ParameterList& ifpackList = p->sublist("Preconditioner Types").sublist("Ifpack2");
  ifpackList.set("Prec Type", "ILUT");

  builder.setParameterList(p);

  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
    lowsFactory = builder.createLinearSolveStrategy("");

  model->set_W_factory(lowsFactory);

  // Create the initial guess
  Teuchos::RCP<Thyra::VectorBase<Scalar> >
    initial_guess = model->getNominalValues().get_x()->clone_v();
  Thyra::V_S(initial_guess.ptr(),Teuchos::ScalarTraits<Scalar>::one());

  Teuchos::RCP<NOX::Thyra::Group> nox_group =
    Teuchos::rcp(new NOX::Thyra::Group(*initial_guess, model, model->create_W_op(), lowsFactory, Teuchos::null, Teuchos::null, Teuchos::null));

  nox_group->computeF();

  // Create the NOX status tests and the solver
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
    Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  converged->addStatusTest(wrms);
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // Create nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> nl_params = Teuchos::parameterList();
  nl_params->set("Nonlinear Solver", "Line Search Based");
  nl_params->sublist("Direction").sublist("Newton").sublist("Linear Solver").set("Tolerance", 1.0e-4);

  // Set output parameters
  nl_params->sublist("Printing").sublist("Output Information").set("Debug",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Warning",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Error",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Test Details",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Details",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Parameters",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Linear Solver Details",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Inner Iteration",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Outer Iteration",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Outer Iteration StatusTest",true);

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(nox_group, combo, nl_params);
  NOX::StatusTest::StatusType solvStatus = solver->solve();

  TEST_ASSERT(solvStatus == NOX::StatusTest::Converged);

  Teuchos::TimeMonitor::summarize();
}

TEUCHOS_UNIT_TEST(NOX_Tpetra_1DFEM, JFNK_NoPrec)
{
  Teuchos::TimeMonitor::zeroOutTimers();

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  // Get default Tpetra template types
  typedef Tpetra::Vector<>::scalar_type Scalar;
  typedef Tpetra::Vector<>::local_ordinal_type LO;
  typedef Tpetra::Vector<>::global_ordinal_type GO;
  typedef Tpetra::Vector<>::node_type Node;

  // Create the model evaluator object
  Scalar x00 = 0.0;
  Scalar x01 = 1.0;
  Teuchos::RCP<EvaluatorTpetra1DFEM<Scalar,LO,GO,Node> > model =
    evaluatorTpetra1DFEM<Scalar,LO,GO,Node>(comm, numGlobalElements, x00, x01);

  Stratimikos::DefaultLinearSolverBuilder builder;

  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
  p->set("Linear Solver Type", "Belos");
  Teuchos::ParameterList& belosList = p->sublist("Linear Solver Types").sublist("Belos");
  belosList.set("Solver Type", "Pseudo Block GMRES");
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set<int>("Maximum Iterations", 200);
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set<int>("Num Blocks", 200);
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set("Verbosity", 0x7f);
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set("Output Frequency", 100);
  belosList.sublist("VerboseObject").set("Verbosity Level", "medium");
  p->set("Preconditioner Type", "None");
  builder.setParameterList(p);

  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
    lowsFactory = builder.createLinearSolveStrategy("");

  model->set_W_factory(lowsFactory);

  // Create the initial guess
  Teuchos::RCP<Thyra::VectorBase<Scalar> >
    initial_guess = model->getNominalValues().get_x()->clone_v();
  Thyra::V_S(initial_guess.ptr(),Teuchos::ScalarTraits<Scalar>::one());

  // Create the JFNK operator
  Teuchos::ParameterList printParams;
  Teuchos::RCP<Teuchos::ParameterList> jfnkParams = Teuchos::parameterList();
  jfnkParams->set("Difference Type","Forward");
  jfnkParams->set("Perturbation Algorithm","KSP NOX 2001");
  jfnkParams->set("lambda",1.0e-4);
  Teuchos::RCP<NOX::Thyra::MatrixFreeJacobianOperator<Scalar> > jfnkOp =
    Teuchos::rcp(new NOX::Thyra::MatrixFreeJacobianOperator<Scalar>(printParams));
  jfnkOp->setParameterList(jfnkParams);
  jfnkParams->print(out);

  // Wrap the model evaluator in a JFNK Model Evaluator
  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > thyraModel =
    Teuchos::rcp(new NOX::MatrixFreeModelEvaluatorDecorator<Scalar>(model));

  // Create the NOX::Thyra::Group
  Teuchos::RCP<NOX::Thyra::Group> nox_group =
    Teuchos::rcp(new NOX::Thyra::Group(*initial_guess, thyraModel, jfnkOp, lowsFactory, Teuchos::null, Teuchos::null, Teuchos::null));

  nox_group->computeF();

  // VERY IMPORTANT!!!  jfnk object needs base evaluation objects.
  // This creates a circular dependency, so use a weak pointer.
  jfnkOp->setBaseEvaluationToNOXGroup(nox_group.create_weak());

  // Create a preconditioner

  // Create the NOX status tests and the solver
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
    Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  converged->addStatusTest(wrms);
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // Create nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> nl_params = Teuchos::parameterList();
  nl_params->set("Nonlinear Solver", "Line Search Based");
  nl_params->sublist("Direction").sublist("Newton").sublist("Linear Solver").set("Tolerance", 1.0e-4);

  // Set output parameters
  nl_params->sublist("Printing").sublist("Output Information").set("Debug",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Warning",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Error",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Test Details",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Details",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Parameters",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Linear Solver Details",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Inner Iteration",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Outer Iteration",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Outer Iteration StatusTest",true);

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(nox_group, combo, nl_params);
  NOX::StatusTest::StatusType solvStatus = solver->solve();

  TEST_ASSERT(solvStatus == NOX::StatusTest::Converged);

  Teuchos::TimeMonitor::summarize();
}


TEUCHOS_UNIT_TEST(NOX_Tpetra_1DFEM, JFNK_UserPrec)
{
  Teuchos::TimeMonitor::zeroOutTimers();

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  // Get default Tpetra template types
  typedef Tpetra::Vector<>::scalar_type Scalar;
  typedef Tpetra::Vector<>::local_ordinal_type LO;
  typedef Tpetra::Vector<>::global_ordinal_type GO;
  typedef Tpetra::Vector<>::node_type Node;

  // Create the model evaluator object
  Scalar x00 = 0.0;
  Scalar x01 = 1.0;
  Teuchos::RCP<EvaluatorTpetra1DFEM<Scalar,LO,GO,Node> > model =
    evaluatorTpetra1DFEM<Scalar,LO,GO,Node>(comm, numGlobalElements, x00, x01);

  Stratimikos::DefaultLinearSolverBuilder builder;

  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
  p->set("Linear Solver Type", "Belos");
  Teuchos::ParameterList& belosList = p->sublist("Linear Solver Types").sublist("Belos");
  belosList.set("Solver Type", "Pseudo Block GMRES");
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set<int>("Maximum Iterations", 200);
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set<int>("Num Blocks", 200);
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set("Verbosity", 0x7f);
  belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set("Output Frequency", 100);
  belosList.sublist("VerboseObject").set("Verbosity Level", "medium");
  p->set("Preconditioner Type", "None");
  builder.setParameterList(p);

  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
    lowsFactory = builder.createLinearSolveStrategy("");

  model->set_W_factory(lowsFactory);

  // Create the initial guess
  Teuchos::RCP<Thyra::VectorBase<Scalar> >
    initial_guess = model->getNominalValues().get_x()->clone_v();
  Thyra::V_S(initial_guess.ptr(),Teuchos::ScalarTraits<Scalar>::one());

  // Create the JFNK operator
  Teuchos::ParameterList printParams;
  Teuchos::RCP<Teuchos::ParameterList> jfnkParams = Teuchos::parameterList();
  jfnkParams->set("Difference Type","Forward");
  jfnkParams->set("Perturbation Algorithm","KSP NOX 2001");
  jfnkParams->set("lambda",1.0e-4);
  Teuchos::RCP<NOX::Thyra::MatrixFreeJacobianOperator<Scalar> > jfnkOp =
    Teuchos::rcp(new NOX::Thyra::MatrixFreeJacobianOperator<Scalar>(printParams));
  jfnkOp->setParameterList(jfnkParams);
  jfnkParams->print(out);

  // Wrap the model evaluator in a JFNK Model Evaluator
  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > thyraModel =
    Teuchos::rcp(new NOX::MatrixFreeModelEvaluatorDecorator<Scalar>(model));

  // Create the Preconditioner operator
  Teuchos::RCP<Thyra::PreconditionerBase<Scalar> > precOp =
    thyraModel->create_W_prec();

  // Create the NOX::Thyra::Group
  Teuchos::RCP<NOX::Thyra::Group> nox_group =
    Teuchos::rcp(new NOX::Thyra::Group(*initial_guess, thyraModel, jfnkOp, lowsFactory, precOp, Teuchos::null, Teuchos::null));

  nox_group->computeF();

  // VERY IMPORTANT!!!  jfnk object needs base evaluation objects.
  // This creates a circular dependency, so use a weak pointer.
  jfnkOp->setBaseEvaluationToNOXGroup(nox_group.create_weak());

  // Create the NOX status tests and the solver
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
    Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  converged->addStatusTest(wrms);
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // Create nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> nl_params = Teuchos::parameterList();
  nl_params->set("Nonlinear Solver", "Line Search Based");
  nl_params->sublist("Direction").sublist("Newton").sublist("Linear Solver").set("Tolerance", 1.0e-4);

  // Set output parameters
  nl_params->sublist("Printing").sublist("Output Information").set("Debug",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Warning",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Error",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Test Details",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Details",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Parameters",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Linear Solver Details",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Inner Iteration",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Outer Iteration",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Outer Iteration StatusTest",true);

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(nox_group, combo, nl_params);
  NOX::StatusTest::StatusType solvStatus = solver->solve();

  TEST_ASSERT(solvStatus == NOX::StatusTest::Converged);

  Teuchos::TimeMonitor::summarize();
}

