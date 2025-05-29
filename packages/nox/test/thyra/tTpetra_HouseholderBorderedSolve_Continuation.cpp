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
#include "LOCA_Stepper.H"
#include "NOX_SolverStats.hpp"
#include "NOX_Observer.hpp"
#include "NOX_Observer_Vector.cpp"

// For solution io
#include "Thyra_TpetraVector.hpp"
#include <iostream>
#include <fstream>
#include <tuple>

// *************************************************************
// *************************************************************
// This test does a continuation in the T_left dirichlet BC. A
// constraint of 2*T_left=T_right is enforced by allowing the source
// term multiplier k to be adjusted as the independent parameter.

// A nox observer is used to check that the constraint is enforced at
// each step of the continuation.

// There are two tests. The first does a manual continuation run
// outside of LOCA. The second test uses LOCA for the continuation.
// *************************************************************
// *************************************************************

namespace test_utils {

// Returns the solution view from the nox vector. Accounts for
// potentially wrapped LOCA Vectors.
auto getConstView(const NOX::Abstract::Vector& nox_vector)
{
  const auto* nox_thyra_vector = dynamic_cast<const NOX::Thyra::Vector*>(&nox_vector);

  // If the cast above failed, pull it out of the constrained vector.
  if (nox_thyra_vector == nullptr) {
    const auto& final_x_nox = *(dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(nox_vector).getXVec());
    nox_thyra_vector = dynamic_cast<const NOX::Thyra::Vector*>(&final_x_nox);

    // If the LOCA::Stepper is building the ConstraintGroup, then there
    // is potentially some double nesting for the arc length constraint
    // depending on if you unroll the nested group.
    if (nox_thyra_vector == nullptr) {
      const auto* loca_nested_extended_vector_ptr = dynamic_cast<const LOCA::MultiContinuation::ExtendedVector*>(&final_x_nox);
      TEUCHOS_ASSERT(loca_nested_extended_vector_ptr != nullptr);
      nox_thyra_vector = Teuchos::rcp_dynamic_cast<const NOX::Thyra::Vector>(loca_nested_extended_vector_ptr->getXVec(),true).get();
    }
  }

  const auto& final_x_thyra = nox_thyra_vector->getThyraVector();
  const auto& final_x_tpetra_const = *(dynamic_cast<const ::Thyra::TpetraVector<NOX::Scalar,NOX::LocalOrdinal,NOX::GlobalOrdinal,NOX::NodeType>&>(final_x_thyra).getConstTpetraVector());
  auto& final_x_tpetra = const_cast<::Tpetra::Vector<NOX::Scalar,NOX::LocalOrdinal,NOX::GlobalOrdinal,NOX::NodeType>&>(final_x_tpetra_const);
  const auto& final_x_view = final_x_tpetra.getLocalViewHost(Tpetra::Access::ReadOnly);

  return final_x_view;
}

auto getComm(const NOX::Abstract::Vector& nox_vector)
{
  const auto* nox_thyra_vector = dynamic_cast<const NOX::Thyra::Vector*>(&nox_vector);

  // If the cast above failed, pull it out of the constrained vector.
  if (nox_thyra_vector == nullptr) {
    const auto& final_x_nox = *(dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(nox_vector).getXVec());
    nox_thyra_vector = dynamic_cast<const NOX::Thyra::Vector*>(&final_x_nox);

    // If the LOCA::Stepper is building the ConstraintGroup, then there
    // is potentially some double nesting for the arc length constraint
    // depending on if you unroll the nested group.
    if (nox_thyra_vector == nullptr) {
      const auto* loca_nested_extended_vector_ptr = dynamic_cast<const LOCA::MultiContinuation::ExtendedVector*>(&final_x_nox);
      TEUCHOS_ASSERT(loca_nested_extended_vector_ptr != nullptr);
      nox_thyra_vector = Teuchos::rcp_dynamic_cast<const NOX::Thyra::Vector>(loca_nested_extended_vector_ptr->getXVec(),true).get();
    }
  }

  const auto& final_x_thyra = nox_thyra_vector->getThyraVector();
  const auto& final_x_tpetra_const = *(dynamic_cast<const ::Thyra::TpetraVector<NOX::Scalar,NOX::LocalOrdinal,NOX::GlobalOrdinal,NOX::NodeType>&>(final_x_thyra).getConstTpetraVector());
  auto& final_x_tpetra = const_cast<::Tpetra::Vector<NOX::Scalar,NOX::LocalOrdinal,NOX::GlobalOrdinal,NOX::NodeType>&>(final_x_tpetra_const);

  return final_x_tpetra.getMap()->getComm();
}

} // namespace test_utils

class WriteParametersObserver : public NOX::Observer {
  bool first_;
  std::string file_name_;
  bool on_print_rank_;
  Teuchos::RCP<LOCA::MultiContinuation::ConstraintModelEvaluator> cme_;
  int w_;
  std::ofstream file_;

public:
  WriteParametersObserver(const std::string& file_name,
                          const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                          const Teuchos::RCP<LOCA::MultiContinuation::ConstraintModelEvaluator>& cme)
    : first_(true),file_name_(file_name),on_print_rank_(false), cme_(cme), w_(16)
  {
    // Make the last rank the print rank because that is where T_right
    // can be found (avoids mpi communication).
    if ( comm->getRank() == (comm->getSize()-1) )
      on_print_rank_ = true;
  }

  ~WriteParametersObserver() { file_.close(); }

  void runPreSolve(const NOX::Solver::Generic& solver)
  {
    if (first_ && on_print_rank_) {
      file_.open(file_name_,std::ios::trunc);
      file_ << std::left;
      file_ << std::setw(w_) << "T_left"
            << std::setw(w_) << "k"
            << std::setw(w_) << "T_right"
            << std::setw(w_) << "Error in Constraint"
            << std::endl;
      first_ = false;
    }
  }

  void runPostSolve(const NOX::Solver::Generic& solver)
  { this->print(solver); }

  void print(const NOX::Solver::Generic& solver) {
    if (on_print_rank_) {
      auto vals = test_utils::getConstView(solver.getSolutionGroup().getX());

      file_ << std::setw(w_) << cme_->getParams().getValue("T_left")
            << std::setw(w_) << cme_->getParams().getValue("k")
            << std::setw(w_) << vals(vals.extent(0)-1,0)
            << std::setw(w_) << std::fabs(vals(vals.extent(0)-1,0)-2.0*cme_->getParams().getValue("T_left"))
            << std::endl;
    }
  }
};

class CheckConstraintObserver : public NOX::Observer {
  void runPostSolve(const NOX::Solver::Generic& solver)
  {
    auto final_x_view = test_utils::getConstView(solver.getSolutionGroup().getX());
    auto comm = test_utils::getComm(solver.getSolutionGroup().getX());
    auto rank = comm->getRank();
    auto size = comm->getSize();

    double T_left{0.0};
    double T_right{0.0};
    if (rank == 0)
      T_left = final_x_view(0,0);
    if (rank == (size-1))
      T_right = final_x_view( (final_x_view.extent(0)-1) ,0);
    Teuchos::broadcast(*comm,0,1,&T_left);
    Teuchos::broadcast(*comm,size-1,1,&T_right);

    if (rank == 0) {
      std::cout << "\nCheck Constraint Observer:\n";
      std::cout << "T_left = " << T_left << "\n";
      std::cout << "T_right = " << T_right << "\n\n";
    }

    TEUCHOS_ASSERT(solver.getStatus() == NOX::StatusTest::Converged);
    TEUCHOS_ASSERT( (2.0 * T_left - T_right) < 1.0e-8);
  }
};

TEUCHOS_UNIT_TEST(NOX_Tpetra_Householder, Manual_Continuation)
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
  // Thyra::V_S(initial_guess.ptr(),Teuchos::ScalarTraits<Scalar>::one());
  Thyra::V_S(initial_guess.ptr(),2.0);

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

  // Set the observer to check the continuation values
  Teuchos::RCP<NOX::Observer> observer = Teuchos::rcp(new CheckConstraintObserver);
  nl_params.sublist("Solver Options").set("Observer",observer);

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
  me_p_indices.push_back(model->get_p_index("k").first);
  me_p_indices.push_back(model->get_p_index("T_left").first);
  Teuchos::RCP<LOCA::Thyra::Group> loca_group = Teuchos::rcp(new LOCA::Thyra::Group(global_data,
                                                                                    *nox_group,
                                                                                    *p_vec,
                                                                                    me_p_indices));

  auto g_names = Teuchos::rcp(new std::vector<std::string>);
  // g_names->push_back("Constraint: T_right=2");
  g_names->push_back("Constraint: 2*T_left=T_right");
  auto x_thyra = ::Thyra::createMember(model->get_x_space(),"x");
  NOX::Thyra::Vector x(x_thyra);
  auto constraints = Teuchos::rcp(new LOCA::MultiContinuation::ConstraintModelEvaluator(model,*p_vec,*g_names,x));

  // Set initial parameter conditions
  constraints->setX(x);
  constraints->setParam(p_vec->getIndex("k"),1.0);
  constraints->setParam(p_vec->getIndex("T_left"),1.2);

  // Create the constraints list
  auto& locaParamsList = pList->sublist("LOCA");
  auto& constraint_list = locaParamsList.sublist("Constraints");
  constraint_list.set("Bordered Solver Method", "Householder");
  constraint_list.set<Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>>("Constraint Object", constraints);
  constraint_list.set("Constraint Parameter Names", g_names);

  auto loca_parser = Teuchos::rcp(new LOCA::Parameter::SublistParser(global_data));
  loca_parser->parseSublists(pList);

  // Sets the params that are allowed to change to enforce constraints
  std::vector<int> param_ids(1);
  param_ids[0] = p_vec->getIndex("k");
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

  solver->solve();
  out << "\n\nStepping from T_left=1.2 to 2.0\n\n";
  auto& nonconst_group = const_cast<NOX::Abstract::Group&>(solver->getSolutionGroup());
  dynamic_cast<LOCA::MultiContinuation::ConstrainedGroup&>(nonconst_group).setParam("T_left",2.0);
  solver->reset();
  solver->solve();

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
          file.open("householder_continuation_manual_solution.txt",std::ios::trunc);
        else
          file.open("householder_continuation_manual_solution.txt",std::ios::app);

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

  // Second step ends in 3 newton iterations
  TEST_EQUALITY(solver->getSolverStatistics()->numNonlinearIterations,3);

  // Check final values
  {
    const auto& group = solver->getSolutionGroup();
    const auto& c_group = dynamic_cast<const LOCA::MultiContinuation::ConstrainedGroup&>(group);

    out << "\nFinal Parameter Value for \"k\" = " << std::setprecision(10) << c_group.getParam(0) << std::endl;
    out << "Final Parameter Value for \"T_left\" = " << std::setprecision(10) << c_group.getParam(p_vec->getIndex("T_left")) << std::endl;

    const double tol = 1.0e-5;
    TEST_FLOATING_EQUALITY(c_group.getParam(p_vec->getIndex("k")),-0.2996638603,tol);
    TEST_FLOATING_EQUALITY(c_group.getParam(p_vec->getIndex("T_left")),2.0,tol);
  }

  // Breaks RCP cyclic dependency
  LOCA::destroyGlobalData(global_data);
}


TEUCHOS_UNIT_TEST(NOX_Tpetra_Householder, LOCA_Continuation)
{
  Teuchos::TimeMonitor::setStackedTimer(Teuchos::rcp(new Teuchos::StackedTimer("LOCA_Continuation"))); // start a new stacked timer
  // Teuchos::TimeMonitor::getStackedTimer()->startBaseTimer();

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
  // Thyra::V_S(initial_guess.ptr(),Teuchos::ScalarTraits<Scalar>::one());
  Thyra::V_S(initial_guess.ptr(),2.0);

  // Create top level nox/loca solver parameter list
  Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::parameterList("Top Level");

  // Create the LOCA parameter list
  auto& loca_params = pList->sublist("LOCA");
  Teuchos::ParameterList& stepperList = loca_params.sublist("Stepper");
  stepperList.set("Continuation Method", "Arc Length");
  stepperList.set("Continuation Parameter", "T_left");
  stepperList.set("Initial Value", 1.0);
  stepperList.set("Max Value", 4.0);
  stepperList.set("Min Value", 1.0);
  stepperList.set("Max Steps", 50);
  stepperList.set("Max Nonlinear Iterations", 20);

  // Combine arc-length and turning point bordered rows & columns
  stepperList.set("Bordered Solver Method", "Nested");
  auto& nestedList = stepperList.sublist("Nested Bordered Solver");
  nestedList.set("Bordered Solver Method", "Householder");

  // Create step size sublist
  Teuchos::ParameterList& stepSizeList = loca_params.sublist("Step Size");
  stepSizeList.set("Initial Step Size", 0.5);
  stepSizeList.set("Min Step Size", 0.5);
  stepSizeList.set("Max Step Size", 1.0);

  // Create NOX parameter list
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
  output_list.set("Details",false);
  output_list.set("Parameters",false);
  output_list.set("Linear Solver Details",true);
  output_list.set("Inner Iteration",true);
  output_list.set("Outer Iteration",true);
  output_list.set("Outer Iteration StatusTest",true);
  output_list.set("Stepper Iteration",true);
  output_list.set("Stepper Details",true);
  output_list.set("Stepper Parameters",true);

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
  me_p_indices.push_back(model->get_p_index("k").first);
  me_p_indices.push_back(model->get_p_index("T_left").first);
  Teuchos::RCP<LOCA::Thyra::Group> loca_group = Teuchos::rcp(new LOCA::Thyra::Group(global_data,
                                                                                    *nox_group,
                                                                                    *p_vec,
                                                                                    me_p_indices));

  auto g_names = Teuchos::rcp(new std::vector<std::string>);
  // g_names->push_back("Constraint: T_right=2");
  g_names->push_back("Constraint: 2*T_left=T_right");
  auto x_thyra = ::Thyra::createMember(model->get_x_space(),"x");
  NOX::Thyra::Vector x(x_thyra);
  auto constraints = Teuchos::rcp(new LOCA::MultiContinuation::ConstraintModelEvaluator(model,*p_vec,*g_names,x));

  // Set initial parameter conditions
  constraints->setX(x);
  constraints->setParam(p_vec->getIndex("k"),1.0);
  constraints->setParam(p_vec->getIndex("T_left"),1.2);

  // Create the constraints list
  auto& locaParamsList = pList->sublist("LOCA");
  auto& constraint_list = locaParamsList.sublist("Constraints");
  constraint_list.set("Bordered Solver Method", "Householder");
  constraint_list.set<Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>>("Constraint Object", constraints);
  // Sets the params that are allowed to change to enforce constraints
  auto p_names = Teuchos::rcp(new std::vector<std::string>);
  p_names->push_back("k");
  constraint_list.set("Constraint Parameter Names", p_names);

  // Set the observer to check the continuation values
  {
    auto observer1 = Teuchos::rcp(new CheckConstraintObserver);
    auto observer2 = Teuchos::rcp(new WriteParametersObserver("householder_continuation_loca.txt",comm,constraints));
    auto observer = Teuchos::rcp(new NOX::ObserverVector);
    observer->pushBack(observer1);
    observer->pushBack(observer2);
    nl_params.sublist("Solver Options").set<Teuchos::RCP<NOX::Observer>>("Observer",observer);
  }

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

  // Create the stepper
  LOCA::Stepper stepper(global_data, loca_group, combo, pList);

  // Perform continuation run
  LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();
  TEST_EQUALITY(status, LOCA::Abstract::Iterator::Finished);

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
          file.open("householder_conintuation_loca_solution.txt",std::ios::trunc);
        else
          file.open("householder_continuation_loca_solution.txt",std::ios::app);

        auto final_x_view = test_utils::getConstView(stepper.getSolutionGroup()->getX());

        for (size_t j=0; j < final_x_view.extent(0); ++j)
          file << final_x_view(j,0) << std::endl;
      }
      comm->barrier();
    }
  }

  // Second step ends in 3 newton iterations
  TEST_EQUALITY(stepper.getSolver()->getSolverStatistics()->numNonlinearIterations,1);
  // TODO: Could LOCA be rebuilding the nox solver at every step???
  TEST_EQUALITY(stepper.getSolver()->getSolverStatistics()->numTotalNonlinearIterations,1);

  // Breaks RCP cyclic dependency
  LOCA::destroyGlobalData(global_data);
}
