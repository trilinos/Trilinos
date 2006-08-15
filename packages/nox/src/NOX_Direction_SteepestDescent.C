// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
SteepestDescent(const Teuchos::RefCountPtr<NOX::GlobalData>& gd, 
		Teuchos::ParameterList& params)
{
  reset(gd, params);
}

NOX::Direction::SteepestDescent::~SteepestDescent()
{

}

bool NOX::Direction::SteepestDescent::
reset(const Teuchos::RefCountPtr<NOX::GlobalData>& gd,
      Teuchos::ParameterList& params)
{
  globalDataPtr = gd;
  utils = gd->getUtils();
  meritFuncPtr = gd->getMeritFunction();

  Teuchos::ParameterList& p = params.sublist("Steepest Descent");

  const string tmp = p.get("Scaling Type", "2-Norm");
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
		 << "\"" << tmp << "\" for \"Scaling Type\"" << endl;
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

void NOX::Direction::SteepestDescent::throwError(const string& functionName, 
						 const string& errorMsg)
{
    if (utils->isPrintType(Utils::Error))
      utils->err() << "NOX::Direction::SteepestDescent::" << functionName 
	   << " - " << errorMsg << endl;
    throw "NOX Error";
}
