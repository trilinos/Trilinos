// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "TSFNOXSolver.H"         
#include "NOX_StatusTest_SafeCombo.H"         
#include "NOX_Parameter_Teuchos2NOX.H"         
#include "TSFLinearSolverBuilder.hpp"

#ifdef HAVE_NOX_ANY

using namespace NOX;
using namespace NOX::TSF;
using namespace Teuchos;
using namespace TSFExtended;

NOXSolver::NOXSolver(const ParameterList& params,
                     const NonlinearOperator<double>& F)
  : F_(F),
    linSolver_(),
    x0_(F_.getInitialGuess()),
    soln_(),
    grp_(),
    statusTest_(),
    solver_(),
    params_(params),
    noxParams_()
{
  TEST_FOR_EXCEPTION(!params_.isSublist("NOX Solver"), runtime_error,
                     "did not find NOX Solver sublist in " << params_);
  
  ParameterList solverSublist = params_.sublist("NOX Solver");

  cerr << "solver sublist = " << solverSublist << endl;

  if (solverSublist.isSublist("Status Test"))
    {
      statusTest_ = StatusTestBuilder::makeStatusTest(solverSublist);
    }
  else
    {
      RefCountPtr<StatusTest::Generic> A = rcp(new StatusTest::NormF(1.0e-12));
      RefCountPtr<StatusTest::Generic> B = rcp(new StatusTest::MaxIters(20));
      statusTest_ = 
        rcp(new StatusTest::SafeCombo(StatusTest::SafeCombo::OR, A, B));
    }
  
  if (solverSublist.isSublist("Linear Solver"))
    {
      linSolver_ = LinearSolverBuilder::createSolver(solverSublist);
    }
  else
    {
      
    }

  grp_ = rcp(new NOX::TSF::Group(x0_, F_, linSolver_));

  NOX::Parameter::Teuchos2NOX converter;
  noxParams_ = converter.toNOX(params_);

  solver_ = rcp(new NOX::Solver::Manager(*grp_, *statusTest_, noxParams_));
}


NOX::StatusTest::StatusType NOXSolver::solve() const 
{
  NOX::StatusTest::StatusType rtn = solver_->solve();

  const NOX::TSF::Group* solnGrp 
    = dynamic_cast<const NOX::TSF::Group*>(&(solver_->getSolutionGroup()));

  TEST_FOR_EXCEPTION(solnGrp==0, runtime_error,
                     "Solution group could not be cast to NOX::TSF::Group");

  const NOX::TSF::Vector* x 
    = dynamic_cast<const NOX::TSF::Vector*>(&(solnGrp->getX()));

  TEST_FOR_EXCEPTION(x==0, runtime_error,
                     "Solution vector could not be cast to NOX::TSF::Vector");
  
  soln_ = x->getTSFVector();

  x0_ = soln_;
  F_.setEvalPt(soln_);

  grp_ = rcp(new NOX::TSF::Group(soln_, F_, linSolver_));


  NOX::Parameter::Teuchos2NOX converter;
  noxParams_ = converter.toNOX(params_);

  solver_ = rcp(new NOX::Solver::Manager(*grp_, *statusTest_, noxParams_));

  return rtn;
}

#endif
