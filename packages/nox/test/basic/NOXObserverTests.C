// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ParameterList.hpp>
#include "NOX_Observer_Vector.hpp"
#include "NOX_Observer_Log.hpp"
#include "NOX_Abstract_Vector.H"
#include "NOX_Solver_Generic.H"

namespace NOX_UNIT_TEST {

  class MockSolver : public NOX::Solver::Generic {
    Teuchos::RCP<NOX::Abstract::Group> nullGroup_;
    Teuchos::ParameterList dummyParamList_;
  public:
    void reset(){}
    void reset(const NOX::Abstract::Vector& initial_guess){}
    void reset(const NOX::Abstract::Vector& initial_guess,
           const Teuchos::RCP<NOX::StatusTest::Generic>& test){}
    NOX::StatusTest::StatusType step(){return NOX::StatusTest::Failed;}
    NOX::StatusTest::StatusType solve(){return NOX::StatusTest::Failed;}
    const NOX::Abstract::Group& getSolutionGroup() const
    {return *nullGroup_;}
    const NOX::Abstract::Group& getPreviousSolutionGroup() const
    {return *nullGroup_;}
    NOX::StatusTest::StatusType getStatus() const {return NOX::StatusTest::Failed;}
    int getNumIterations() const{return 0;}
    const Teuchos::ParameterList& getList() const
    {return dummyParamList_;}
    Teuchos::RCP< const NOX::Abstract::Group > getSolutionGroupPtr() const
    {return Teuchos::null;}
    Teuchos::RCP< const NOX::Abstract::Group > getPreviousSolutionGroupPtr() const
    {return Teuchos::null;}
    Teuchos::RCP< const Teuchos::ParameterList > getListPtr() const
    {return Teuchos::null;}
    Teuchos::RCP<const NOX::SolverStats> getSolverStatistics() const
    { return Teuchos::null; }
  };

  class MockVector : public NOX::Abstract::Vector {
    NOX::Abstract::Vector& init(double ) {return *this;}
    NOX::Abstract::Vector& random(bool , int ) {return *this;}
    NOX::Abstract::Vector& abs(const NOX::Abstract::Vector& ) {return *this;}
    NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& ) {return *this;}
    NOX::Abstract::Vector& reciprocal(const NOX::Abstract::Vector& ) {return *this;}
    NOX::Abstract::Vector& scale(double ) {return *this;}
    NOX::Abstract::Vector& scale(const NOX::Abstract::Vector& ) {return *this;}
    NOX::Abstract::Vector& update(double , const NOX::Abstract::Vector& , double ) {return *this;}
    NOX::Abstract::Vector& update(double , const NOX::Abstract::Vector& ,
                                  double , const NOX::Abstract::Vector& ,
                                  double ) {return *this;}
    Teuchos::RCP<NOX::Abstract::Vector>
    clone(NOX::CopyType ) const {return Teuchos::null;}
    double norm(NOX::Abstract::Vector::NormType ) const {return 0.0;}
    double norm(const NOX::Abstract::Vector& ) const {return 0.0;}
    double innerProduct(const NOX::Abstract::Vector& y) const {return 0.0;}
    NOX::size_type length() const {return 0;}
  };

  TEUCHOS_UNIT_TEST(PrePostOperatorVector, all)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    NOX::ObserverVector obs_vec;

    RCP<NOX::ObserverLog> mock_obs_1 = rcp(new NOX::ObserverLog);
    RCP<NOX::ObserverLog> mock_obs_2 = rcp(new NOX::ObserverLog);
    RCP<NOX::ObserverLog> mock_obs_3 = rcp(new NOX::ObserverLog);

    obs_vec.pushBack(mock_obs_1);
    obs_vec.pushBack(mock_obs_2);
    obs_vec.pushBack(mock_obs_3);

    obs_vec.popBack();

    MockSolver solver;
    MockVector vec;

    obs_vec.runPreIterate(solver);

    obs_vec.runPostIterate(solver);
    obs_vec.runPostIterate(solver);

    obs_vec.runPreSolve(solver);
    obs_vec.runPreSolve(solver);
    obs_vec.runPreSolve(solver);

    obs_vec.runPostSolve(solver);
    obs_vec.runPostSolve(solver);
    obs_vec.runPostSolve(solver);
    obs_vec.runPostSolve(solver);

    obs_vec.runPreSolutionUpdate(vec,solver);
    obs_vec.runPreSolutionUpdate(vec,solver);
    obs_vec.runPreSolutionUpdate(vec,solver);
    obs_vec.runPreSolutionUpdate(vec,solver);
    obs_vec.runPreSolutionUpdate(vec,solver);

    obs_vec.runPostSolutionUpdate(solver);
    obs_vec.runPostSolutionUpdate(solver);
    obs_vec.runPostSolutionUpdate(solver);
    obs_vec.runPostSolutionUpdate(solver);
    obs_vec.runPostSolutionUpdate(solver);
    obs_vec.runPostSolutionUpdate(solver);

    obs_vec.runPreLineSearch(solver);
    obs_vec.runPreLineSearch(solver);
    obs_vec.runPreLineSearch(solver);
    obs_vec.runPreLineSearch(solver);
    obs_vec.runPreLineSearch(solver);
    obs_vec.runPreLineSearch(solver);
    obs_vec.runPreLineSearch(solver);

    obs_vec.runPostLineSearch(solver);
    obs_vec.runPostLineSearch(solver);
    obs_vec.runPostLineSearch(solver);
    obs_vec.runPostLineSearch(solver);
    obs_vec.runPostLineSearch(solver);
    obs_vec.runPostLineSearch(solver);
    obs_vec.runPostLineSearch(solver);
    obs_vec.runPostLineSearch(solver);

    TEST_EQUALITY(mock_obs_1->preIterateCount(), 1);
    TEST_EQUALITY(mock_obs_1->postIterateCount(), 2);
    TEST_EQUALITY(mock_obs_1->preSolveCount(), 3);
    TEST_EQUALITY(mock_obs_1->postSolveCount(), 4);
    TEST_EQUALITY(mock_obs_1->preSolutionUpdateCount(), 5);
    TEST_EQUALITY(mock_obs_1->postSolutionUpdateCount(), 6);
    TEST_EQUALITY(mock_obs_1->preLineSearchCount(), 7);
    TEST_EQUALITY(mock_obs_1->postLineSearchCount(), 8);

    TEST_EQUALITY(mock_obs_2->preIterateCount(), 1);
    TEST_EQUALITY(mock_obs_2->postIterateCount(), 2);
    TEST_EQUALITY(mock_obs_2->preSolveCount(), 3);
    TEST_EQUALITY(mock_obs_2->postSolveCount(), 4);
    TEST_EQUALITY(mock_obs_2->preSolutionUpdateCount(), 5);
    TEST_EQUALITY(mock_obs_1->postSolutionUpdateCount(), 6);
    TEST_EQUALITY(mock_obs_2->preLineSearchCount(), 7);
    TEST_EQUALITY(mock_obs_2->postLineSearchCount(), 8);

    TEST_EQUALITY(mock_obs_3->preIterateCount(), 0);
    TEST_EQUALITY(mock_obs_3->postIterateCount(), 0);
    TEST_EQUALITY(mock_obs_3->preSolveCount(), 0);
    TEST_EQUALITY(mock_obs_3->postSolveCount(), 0);
    TEST_EQUALITY(mock_obs_3->preSolutionUpdateCount(), 0);
    TEST_EQUALITY(mock_obs_3->postSolutionUpdateCount(), 0);
    TEST_EQUALITY(mock_obs_3->preLineSearchCount(), 0);
    TEST_EQUALITY(mock_obs_3->postLineSearchCount(), 0);

    TEST_EQUALITY(mock_obs_1->getCallOrder().size(),36);
    TEST_EQUALITY(mock_obs_2->getCallOrder().size(),36);
    TEST_EQUALITY(mock_obs_3->getCallOrder().size(),0);

    const auto& call_order = mock_obs_1->getCallOrder();

    TEST_EQUALITY(call_order[0], "runPreIterate");

    TEST_EQUALITY(call_order[1], "runPostIterate");
    TEST_EQUALITY(call_order[2], "runPostIterate");

    TEST_EQUALITY(call_order[3], "runPreSolve");
    TEST_EQUALITY(call_order[4], "runPreSolve");
    TEST_EQUALITY(call_order[5], "runPreSolve");

    TEST_EQUALITY(call_order[6], "runPostSolve");
    TEST_EQUALITY(call_order[7], "runPostSolve");
    TEST_EQUALITY(call_order[8], "runPostSolve");
    TEST_EQUALITY(call_order[9], "runPostSolve");

    TEST_EQUALITY(call_order[10], "runPreSolutionUpdate");
    TEST_EQUALITY(call_order[11], "runPreSolutionUpdate");
    TEST_EQUALITY(call_order[12], "runPreSolutionUpdate");
    TEST_EQUALITY(call_order[13], "runPreSolutionUpdate");
    TEST_EQUALITY(call_order[14], "runPreSolutionUpdate");

    TEST_EQUALITY(call_order[15], "runPostSolutionUpdate");
    TEST_EQUALITY(call_order[16], "runPostSolutionUpdate");
    TEST_EQUALITY(call_order[17], "runPostSolutionUpdate");
    TEST_EQUALITY(call_order[18], "runPostSolutionUpdate");
    TEST_EQUALITY(call_order[19], "runPostSolutionUpdate");
    TEST_EQUALITY(call_order[20], "runPostSolutionUpdate");

    TEST_EQUALITY(call_order[21], "runPreLineSearch");
    TEST_EQUALITY(call_order[22], "runPreLineSearch");
    TEST_EQUALITY(call_order[23], "runPreLineSearch");
    TEST_EQUALITY(call_order[24], "runPreLineSearch");
    TEST_EQUALITY(call_order[25], "runPreLineSearch");
    TEST_EQUALITY(call_order[26], "runPreLineSearch");
    TEST_EQUALITY(call_order[27], "runPreLineSearch");

    TEST_EQUALITY(call_order[28], "runPostLineSearch");
    TEST_EQUALITY(call_order[29], "runPostLineSearch");
    TEST_EQUALITY(call_order[30], "runPostLineSearch");
    TEST_EQUALITY(call_order[31], "runPostLineSearch");
    TEST_EQUALITY(call_order[32], "runPostLineSearch");
    TEST_EQUALITY(call_order[33], "runPostLineSearch");
    TEST_EQUALITY(call_order[34], "runPostLineSearch");
    TEST_EQUALITY(call_order[35], "runPostLineSearch");

  }

}
