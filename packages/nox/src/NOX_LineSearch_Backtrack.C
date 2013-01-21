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

#include "NOX_LineSearch_Backtrack.H" // class definition

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Utils.H"
#include "NOX_MeritFunction_Generic.H"
#include "NOX_StatusTest_FiniteValue.H"
#include "NOX_GlobalData.H"

NOX::LineSearch::Backtrack::
Backtrack(const Teuchos::RCP<NOX::GlobalData>& gd,
	  Teuchos::ParameterList& params)
{
  reset(gd, params);
}

NOX::LineSearch::Backtrack::~Backtrack()
{

}

bool NOX::LineSearch::Backtrack::
reset(const Teuchos::RCP<NOX::GlobalData>& gd,
      Teuchos::ParameterList& params)
{ 
  utils = gd->getUtils();
  meritFunctionPtr = gd->getMeritFunction();

  Teuchos::ParameterList& p = params.sublist("Backtrack");

  minStep = p.get("Minimum Step", 1.0e-12);
  defaultStep = p.get("Default Step", 1.0);
  recoveryStep = p.get("Recovery Step", defaultStep);
  maxIters = p.get("Max Iters", 100);

  reductionFactor = p.get("Reduction Factor", 0.5);
  if ((reductionFactor <= 0.0)  || (reductionFactor >= 1.0)) {
    utils->err() << "NOX::LineSearch::Backtrack::reset - Invalid choice \"" 
		 << reductionFactor << "\" for \"Reduction Factor\"!  " 
		 << "Value must be greater than zero and less than 1.0."
		 << std::endl;
    throw "NOX Error";
  }

  return true;
}

bool NOX::LineSearch::Backtrack::
compute(NOX::Abstract::Group& grp, double& step, 
	const NOX::Abstract::Vector& dir,
	const NOX::Solver::Generic& s)
{
  const Abstract::Group& oldGrp = s.getPreviousSolutionGroup();
  double oldF = meritFunctionPtr->computef(oldGrp);
  double newF;
  bool isFailed = false;

  step = defaultStep;
  grp.computeX(oldGrp, dir, step);

  NOX::Abstract::Group::ReturnType rtype;

  rtype = grp.computeF();    
  if (rtype != NOX::Abstract::Group::Ok)
  {
    utils->err() << "NOX::LineSearch::BackTrack::compute - Unable to compute F"
		<< std::endl;
    throw "NOX Error";
  }

  newF = meritFunctionPtr->computef(grp);
  int nIters = 1;

  if (utils->isPrintType(Utils::InnerIteration)) 
  {
   utils->out() << "\n" << Utils::fill(72) << "\n" 
	       << "-- Backtrack Line Search -- \n";
  }

  NOX::StatusTest::FiniteValue checkNAN;

  while ( ((newF >= oldF) || (checkNAN.finiteNumberTest(newF) !=0))
	 && (!isFailed)) 
  {

    if (utils->isPrintType(Utils::InnerIteration)) 
    {
      utils->out() << std::setw(3) << nIters << ":";
      utils->out() << " step = " << utils->sciformat(step);
      utils->out() << " old f = " << utils->sciformat(oldF);
      utils->out() << " new f = " << utils->sciformat(newF);
      utils->out() << std::endl;
    }

    nIters ++;
    step = step * reductionFactor;

    if ((step < minStep) || (nIters > maxIters)) 
    {
      isFailed = true;
      step = recoveryStep;
    }

    grp.computeX(oldGrp, dir, step);

    rtype = grp.computeF();    
    if (rtype != NOX::Abstract::Group::Ok)
    {
      utils->err() << "NOX::LineSearch::BackTrack::compute - Unable to compute F" << std::endl;
      throw "NOX Error";
    }

    newF = meritFunctionPtr->computef(grp);
  } 

  if (utils->isPrintType(Utils::InnerIteration)) 
  {
    utils->out() << std::setw(3) << nIters << ":";
    utils->out() << " step = " << utils->sciformat(step);
    utils->out() << " old f = " << utils->sciformat(oldF);
    utils->out() << " new f = " << utils->sciformat(newF);
    if (isFailed)
      utils->out() << " (USING RECOVERY STEP!)" << std::endl;
    else
      utils->out() << " (STEP ACCEPTED!)" << std::endl;
    utils->out() << Utils::fill(72) << "\n" << std::endl;
  }

  return (!isFailed);
}

