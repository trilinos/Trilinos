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

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "NOX.H"
#include "NOX_PrePostOperator_Vector.H"
#include "NOX_Solver_Generic.H"

namespace NOX_UNIT_TEST {

  class MockPPOp : public NOX::Abstract::PrePostOperator {

  private:

    int pre_it_count_;
    int post_it_count_;
    int pre_solve_count_;
    int post_solve_count_;

  public:

    MockPPOp() :
      pre_it_count_(0),
      post_it_count_(0),
      pre_solve_count_(0),
      post_solve_count_(0)
    {}

    void runPreIterate(const NOX::Solver::Generic& solver)
    {pre_it_count_ += 1;}

    void runPostIterate(const NOX::Solver::Generic& solver)
    {post_it_count_ += 1;}

    void runPreSolve(const NOX::Solver::Generic& solver)
    {pre_solve_count_ += 1;}

    void runPostSolve(const NOX::Solver::Generic& solver)
    {post_solve_count_ += 1;}

    int preIterateCount() const
    {return pre_it_count_;}

    int postIterateCount() const
    {return post_it_count_;}

    int preSolveCount() const
    {return pre_solve_count_;}

    int postSolveCount() const
    {return post_solve_count_;}

  };

  class MockSolver : public NOX::Solver::Generic {

  public:

    void reset(const NOX::Abstract::Vector& initial_guess){}

    void reset(const NOX::Abstract::Vector& initial_guess,
           const Teuchos::RCP<NOX::StatusTest::Generic>& test){}

    NOX::StatusTest::StatusType getStatus(){return NOX::StatusTest::Failed;}

    NOX::StatusTest::StatusType step(){return NOX::StatusTest::Failed;}

    NOX::StatusTest::StatusType solve(){return NOX::StatusTest::Failed;}

    const NOX::Abstract::Group& getSolutionGroup() const{}

    const NOX::Abstract::Group& getPreviousSolutionGroup() const{}

    int getNumIterations() const{return 0;}

    const Teuchos::ParameterList& getList() const {}

  };

  TEUCHOS_UNIT_TEST(PrePostOperatorVector, all)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    NOX::PrePostOperatorVector ppop_vec;

    RCP<MockPPOp> mock_ppop_1 = rcp(new MockPPOp);
    RCP<MockPPOp> mock_ppop_2 = rcp(new MockPPOp);
    RCP<MockPPOp> mock_ppop_3 = rcp(new MockPPOp);

    ppop_vec.pushBack(mock_ppop_1);
    ppop_vec.pushBack(mock_ppop_2);
    ppop_vec.pushBack(mock_ppop_3);

    ppop_vec.popBack();

    MockSolver solver;

    ppop_vec.runPreIterate(solver);

    ppop_vec.runPostIterate(solver);
    ppop_vec.runPostIterate(solver);

    ppop_vec.runPreSolve(solver);
    ppop_vec.runPreSolve(solver);
    ppop_vec.runPreSolve(solver);

    ppop_vec.runPostSolve(solver);
    ppop_vec.runPostSolve(solver);
    ppop_vec.runPostSolve(solver);
    ppop_vec.runPostSolve(solver);

    TEST_EQUALITY(mock_ppop_1->preIterateCount(), 1);
    TEST_EQUALITY(mock_ppop_1->postIterateCount(), 2);
    TEST_EQUALITY(mock_ppop_1->preSolveCount(), 3);
    TEST_EQUALITY(mock_ppop_1->postSolveCount(), 4);

    TEST_EQUALITY(mock_ppop_2->preIterateCount(), 1);
    TEST_EQUALITY(mock_ppop_2->postIterateCount(), 2);
    TEST_EQUALITY(mock_ppop_2->preSolveCount(), 3);
    TEST_EQUALITY(mock_ppop_2->postSolveCount(), 4);

    TEST_EQUALITY(mock_ppop_3->preIterateCount(), 0);
    TEST_EQUALITY(mock_ppop_3->postIterateCount(), 0);
    TEST_EQUALITY(mock_ppop_3->preSolveCount(), 0);
    TEST_EQUALITY(mock_ppop_3->postSolveCount(), 0);

  }

}
