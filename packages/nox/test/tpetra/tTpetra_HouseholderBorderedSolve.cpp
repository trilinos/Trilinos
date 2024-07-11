// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_StackedTimer.hpp"

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
#include "ME_Tpetra_1DFEM.hpp"

#include "NOX_Thyra_MatrixFreeJacobianOperator.hpp"
#include "NOX_MatrixFree_ModelEvaluatorDecorator.hpp"
#include "NOX_TpetraTypedefs.hpp"
#include "LOCA_Tpetra_Factory.hpp"
#include "LOCA_Thyra_Group.H"
#include "LOCA_MultiContinuation_ConstrainedGroup.H"
#include "LOCA_Tpetra_ConstraintModelEvaluator.hpp"
#include "LOCA_Parameter_SublistParser.H"
#include "NOX_SolverStats.hpp"

// For solution io
#include "Thyra_TpetraVector.hpp"
#include <iostream>
#include <fstream>

TEUCHOS_UNIT_TEST(NOX_Tpetra_Householder, BasicSolve)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  // Get default Tpetra template types
  using Scalar = NOX::Scalar;
  using LO = NOX::LocalOrdinal;
  using GO = NOX::GlobalOrdinal;
  using Node = NOX::NodeType;

  // Create the model evaluator object
  Scalar x00 = 0.0;
  Scalar x01 = 1.0;
  const Tpetra::global_size_t numGlobalElements = 100;
  Teuchos::RCP<EvaluatorTpetra1DFEM<Scalar,LO,GO,Node> > model =
    evaluatorTpetra1DFEM<Scalar,LO,GO,Node>(comm, numGlobalElements, x00, x01);

  // Create the linear solver and register on model evaluator
  {
    Stratimikos::DefaultLinearSolverBuilder builder;

    Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
    p->set("Linear Solver Type", "Belos");
    Teuchos::ParameterList& belosList = p->sublist("Linear Solver Types").sublist("Belos");
    belosList.set("Solver Type", "Pseudo Block GMRES");
    belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set<int>("Maximum Iterations", 200);
    belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set<int>("Num Blocks", 200);
    belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set("Verbosity", Belos::Errors+Belos::IterationDetails+Belos::FinalSummary);
    belosList.sublist("Solver Types").sublist("Pseudo Block GMRES").set("Output Frequency", 5);
    belosList.sublist("VerboseObject").set("Verbosity Level", "medium");
    p->set("Preconditioner Type", "Ifpack2");
    // p->set("Preconditioner Type", "None");
    Teuchos::ParameterList& ifpackList = p->sublist("Preconditioner Types").sublist("Ifpack2");
    ifpackList.set("Prec Type", "ILUT");

    builder.setParameterList(p);

    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
      lowsFactory = builder.createLinearSolveStrategy("");

    model->set_W_factory(lowsFactory);
  }

  // Create the initial guess
  Teuchos::RCP<Thyra::VectorBase<Scalar> >
    initial_guess = model->getNominalValues().get_x()->clone_v();
  Thyra::V_S(initial_guess.ptr(),Teuchos::ScalarTraits<Scalar>::one());

  // Create top level nox/loca solver parameter list
  Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::parameterList("Top Level");

  // Create nox parameter list
  auto& nl_params = pList->sublist("NOX");
  nl_params.set("Nonlinear Solver", "Line Search Based");
  nl_params.sublist("Direction").sublist("Newton").sublist("Linear Solver").set("Tolerance", 1.0e-8);
  auto& ls_params = nl_params.sublist("Line Search");
  ls_params.set("Method","Full Step");
  auto& output_list = nl_params.sublist("Printing").sublist("Output Information");
  output_list.set("Debug",true);
  output_list.set("Warning",true);
  output_list.set("Error",true);
  output_list.set("Test Details",true);
  output_list.set("Details",true);
  output_list.set("Parameters",true);
  output_list.set("Linear Solver Details",true);
  output_list.set("Inner Iteration",true);
  output_list.set("Outer Iteration",true);
  output_list.set("Outer Iteration StatusTest",true);

  // Create the LOCA Group:
  // (NOX::Thyra::Group-->LOCA::Thyra::Group-->LOCA::Constrained::Group)
  // For Tpetra Householder, we need to actively set the
  // preconditioner and preconditioner factory so that it uses the
  // precOp separate from the Jacobian operator. Householder replaces
  // the Jacobian operator with a matrix-free version that has the
  // uv^T tacked on.
  auto explicit_jacobian = model->create_W_op();
  auto prec_matrix = Teuchos::rcp(new Thyra::DefaultPreconditioner<NOX::Scalar>(Teuchos::null,explicit_jacobian));
  TEST_ASSERT(nonnull(model->get_W_factory()->getPreconditionerFactory()));
  Teuchos::RCP<NOX::Thyra::Group> nox_group =
    Teuchos::rcp(new NOX::Thyra::Group(*initial_guess,
                                       model,
                                       explicit_jacobian,
                                       model->get_W_factory(),
                                       prec_matrix, // Reuse Jac for approx preconditioner
                                       model->get_W_factory()->getPreconditionerFactory(),
                                       Teuchos::null));

  Teuchos::RCP<LOCA::Abstract::Factory> tpetra_factory = Teuchos::rcp(new LOCA::Tpetra::Factory);

  Teuchos::RCP<LOCA::GlobalData> global_data = LOCA::createGlobalData(pList, tpetra_factory);

  Teuchos::RCP<LOCA::ParameterVector> p_vec = Teuchos::rcp(new LOCA::ParameterVector);
  p_vec->addParameter("k", 1.0); // Source term multiplier
  p_vec->addParameter("T_left", 1.2); // Source term multiplier

  std::vector<int> me_p_indices;
  me_p_indices.push_back(2);
  me_p_indices.push_back(4);
  Teuchos::RCP<LOCA::Thyra::Group> loca_group = Teuchos::rcp(new LOCA::Thyra::Group(global_data,
                                                                                    *nox_group,
                                                                                    *p_vec,
                                                                                    me_p_indices));

  auto g_names = Teuchos::rcp(new std::vector<std::string>);
  g_names->push_back("Constraint: T_right=2");
  g_names->push_back("Constraint: 2*T_left=T_right");
  auto x_thyra = ::Thyra::createMember(model->get_x_space(),"x");
  NOX::Thyra::Vector x(x_thyra);
  auto constraints = Teuchos::rcp(new LOCA::MultiContinuation::ConstraintModelEvaluator(model,*p_vec,*g_names,x));

  // Set initial parameter conditions
  constraints->setX(x);
  constraints->setParam(0,1.0);
  constraints->setParam(1,1.2);

  // Create the constraints list
  auto& locaParamsList = pList->sublist("LOCA");
  auto& constraint_list = locaParamsList.sublist("Constraints");
  constraint_list.set("Bordered Solver Method", "Householder");
  constraint_list.set<Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>>("Constraint Object", constraints);
  constraint_list.set("Constraint Parameter Names", g_names);

  auto loca_parser = Teuchos::rcp(new LOCA::Parameter::SublistParser(global_data));
  loca_parser->parseSublists(pList);

  std::vector<int> param_ids(2);
  param_ids[0] = 0;
  param_ids[1] = 1;
  auto constraint_list_ptr = Teuchos::rcpFromRef(constraint_list);
  Teuchos::RCP<LOCA::MultiContinuation::ConstrainedGroup> loca_constrained_group =
    Teuchos::rcp(new LOCA::MultiContinuation::ConstrainedGroup(global_data,
                                                               loca_parser,
                                                               constraint_list_ptr,
                                                               loca_group,
                                                               constraints,
                                                               param_ids,
                                                               false));

  loca_constrained_group->computeF();

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
    Teuchos::rcp(new NOX::StatusTest::MaxIters(10));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // Create the solver
  // auto solver = NOX::Solver::buildSolver(nox_group, combo, Teuchos::rcpFromRef(pList->sublist("NOX")));
  // auto solver = NOX::Solver::buildSolver(loca_group, combo, Teuchos::rcpFromRef(pList->sublist("NOX")));
  auto solver = NOX::Solver::buildSolver(loca_constrained_group, combo, Teuchos::rcpFromRef(pList->sublist("NOX")));

  NOX::StatusTest::StatusType solvStatus = solver->solve();

  // Output
  {
    Teuchos::TimeMonitor::getStackedTimer()->stopBaseTimer();
    Teuchos::StackedTimer::OutputOptions options;
    options.output_fraction = true;
    options.output_minmax = true;
    Teuchos::TimeMonitor::getStackedTimer()->report(out,comm,options);
  }

  // Write solution to file
  const bool printSolution = true;
  if (printSolution) {
    for (int i=0; i < comm->getSize(); ++i) {
      if (comm->getRank() == i) {
        std::ofstream file;
        if (comm->getRank() == 0)
          file.open("householder_solution.txt",std::ios::trunc);
        else
          file.open("householder_solution.txt",std::ios::app);

        const auto& final_x = solver->getSolutionGroup().getX();
        const auto& final_x_nox = *(dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(final_x).getXVec());
        const auto& final_x_thyra = dynamic_cast<const NOX::Thyra::Vector&>(final_x_nox).getThyraVector();
        const auto& final_x_tpetra_const = *(dynamic_cast<const ::Thyra::TpetraVector<NOX::Scalar,NOX::LocalOrdinal,NOX::GlobalOrdinal,NOX::NodeType>&>(final_x_thyra).getConstTpetraVector());
        auto& final_x_tpetra = const_cast<::Tpetra::Vector<NOX::Scalar,NOX::LocalOrdinal,NOX::GlobalOrdinal,NOX::NodeType>&>(final_x_tpetra_const);
        const auto& final_x_view = final_x_tpetra.getLocalViewHost(Tpetra::Access::ReadOnly);
        for (size_t j=0; j < final_x_view.extent(0); ++j)
          file << final_x_view(j,0) << std::endl;
      }
      comm->barrier();
    }
  }

  TEST_ASSERT(solvStatus == NOX::StatusTest::Converged);
  TEST_EQUALITY(solver->getSolverStatistics()->numNonlinearIterations,5);

  // Check final values
  {
    const auto& group = solver->getSolutionGroup();
    const auto& c_group = dynamic_cast<const LOCA::MultiContinuation::ConstrainedGroup&>(group);

    out << "\nFinal Parameter Value for \"k\" = " << std::setprecision(10) << c_group.getParam(0) << std::endl;
    out << "Final Parameter Value for \"T_left\" = " << std::setprecision(10) << c_group.getParam(1) << std::endl;

    const double tol = 1.0e-8;
    // Check parameter values
    TEST_FLOATING_EQUALITY(c_group.getParam(0),-0.5993277206,tol);
    TEST_FLOATING_EQUALITY(c_group.getParam(1),1.0,tol);
    // Check constraint values
    {
      // get the solution
      const auto& final_x = solver->getSolutionGroup().getX();
      const auto& final_x_nox = *(dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(final_x).getXVec());
      const auto& final_x_thyra = dynamic_cast<const NOX::Thyra::Vector&>(final_x_nox).getThyraVector();
      const auto& final_x_tpetra_const = *(dynamic_cast<const ::Thyra::TpetraVector<NOX::Scalar,NOX::LocalOrdinal,NOX::GlobalOrdinal,NOX::NodeType>&>(final_x_thyra).getConstTpetraVector());
      auto& final_x_tpetra = const_cast<::Tpetra::Vector<NOX::Scalar,NOX::LocalOrdinal,NOX::GlobalOrdinal,NOX::NodeType>&>(final_x_tpetra_const);
      const auto& final_x_view = final_x_tpetra.getLocalViewHost(Tpetra::Access::ReadOnly);

      // Left T value
      if (comm->getRank() == 0) { 
        TEST_FLOATING_EQUALITY(final_x_view(0,0),1.0,tol);
      }
      // Right T value
      if (comm->getRank() == (comm->getSize()-1)) {
        TEST_FLOATING_EQUALITY(final_x_view(final_x_view.extent(0)-1,0),2.0,tol);        
      }
    }
  }

  // Breaks RCP cyclic dependency
  LOCA::destroyGlobalData(global_data);
}
