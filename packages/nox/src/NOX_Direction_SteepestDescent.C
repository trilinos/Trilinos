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

#include "NOX_Common.H"
#include "NOX_Direction_SteepestDescent.H" // class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Utils.H"
#include "NOX_MeritFunction_Generic.H"
#include "NOX_GlobalData.H"

NOX::Direction::SteepestDescent::
SteepestDescent(const Teuchos::RCP<NOX::GlobalData>& gd, 
		Teuchos::ParameterList& params)
{
  reset(gd, params);
}

NOX::Direction::SteepestDescent::~SteepestDescent()
{

}

bool NOX::Direction::SteepestDescent::
reset(const Teuchos::RCP<NOX::GlobalData>& gd,
      Teuchos::ParameterList& params)
{
  globalDataPtr = gd;
  utils = gd->getUtils();
  meritFuncPtr = gd->getMeritFunction();

  Teuchos::ParameterList& p = params.sublist("Steepest Descent");

  const std::string tmp = p.get("Scaling Type", "2-Norm");
  if (tmp == "2-Norm")
    scaleType = NOX::Direction::SteepestDescent::TwoNorm;
  else if (tmp == "F 2-Norm")
    scaleType = NOX::Direction::SteepestDescent::FunctionTwoNorm;
  else if (tmp == "Quadratic Model Min")
    scaleType = NOX::Direction::SteepestDescent::QuadMin;
  else if (tmp == "None")
    scaleType = NOX::Direction::SteepestDescent::None;
  else {
    utils->out() << "NOX::Direction::SteepestDescent::reset - Invalid choice "
		 << "\"" << tmp << "\" for \"Scaling Type\"" << std::endl;
    throw "NOX Error";
  }

 return true;
}

bool NOX::Direction::SteepestDescent::compute(Abstract::Vector& dir, 
				 Abstract::Group& soln, 
				 const Solver::Generic& solver) 
{
  NOX::Abstract::Group::ReturnType status;

  // Compute F at current solution
  status = soln.computeF();
  if (status != NOX::Abstract::Group::Ok) 
    throwError("compute", "Unable to compute F");

  // Compute Jacobian at current solution
  status = soln.computeJacobian();
  if (status != NOX::Abstract::Group::Ok) 
    throwError("compute", "Unable to compute Jacobian");

  // Scale
  switch (scaleType) {

  case NOX::Direction::SteepestDescent::TwoNorm:
    
    meritFuncPtr->computeGradient(soln, dir);
    dir.scale(-1.0/dir.norm());
    break;
    
  case NOX::Direction::SteepestDescent::FunctionTwoNorm:
    
    meritFuncPtr->computeGradient(soln, dir);
    dir.scale(-1.0/soln.getNormF());
    break;
    
  case NOX::Direction::SteepestDescent::QuadMin:
    
    meritFuncPtr->computeQuadraticMinimizer(soln, dir);
      
    break;
    
  case NOX::Direction::SteepestDescent::None:
    
    meritFuncPtr->computeGradient(soln, dir);
    dir.scale( -1.0 );
    break;

  default:
    
    throwError("compute", "Invalid scaleType");
    
  }

  return true;
}

bool NOX::Direction::SteepestDescent::compute(Abstract::Vector& dir, 
				 Abstract::Group& soln, 
				 const Solver::LineSearchBased& solver) 
{
  return NOX::Direction::Generic::compute( dir, soln, solver );
}

void NOX::Direction::SteepestDescent::throwError(const std::string& functionName, 
						 const std::string& errorMsg)
{
    if (utils->isPrintType(Utils::Error))
      utils->err() << "NOX::Direction::SteepestDescent::" << functionName 
	   << " - " << errorMsg << std::endl;
    throw "NOX Error";
}
