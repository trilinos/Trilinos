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

#include "NOX_Parameter_List.H"
#include "LOCA_AnasaziOperator_Manager.H"
#include "LOCA_AnasaziOperator_JacobianInverse.H"
#include "LOCA_AnasaziOperator_ShiftInvert.H"
#include "LOCA_AnasaziOperator_Cayley.H"
#include "LOCA_ErrorCheck.H"

LOCA::AnasaziOperator::Manager::Manager(NOX::Parameter::List& eigenParams_,
					NOX::Parameter::List& solverParams_,
					NOX::Abstract::Group& grp_)
  : opName(),
    op(NULL)
{
  reset(eigenParams_, solverParams_, grp_);
}

LOCA::AnasaziOperator::Manager::~Manager()
{
  delete op;
}

NOX::Abstract::Group::ReturnType 
LOCA::AnasaziOperator::Manager::reset(NOX::Parameter::List& eigenParams,
				      NOX::Parameter::List& solverParams,
				      NOX::Abstract::Group& grp)
{
  string newOpName = eigenParams.getParameter("Operator","Jacobian Inverse");

  if (opName != newOpName) {
    opName = newOpName;
    if (op)
      delete op;

    if (opName == "Jacobian Inverse")
      op = new LOCA::AnasaziOperator::JacobianInverse(eigenParams, 
						      solverParams,
						      grp);
    else if (opName == "Shift-Invert")
      op = new LOCA::AnasaziOperator::ShiftInvert(eigenParams, 
						  solverParams,
						  grp);
    else if (opName == "Cayley")
      op = new LOCA::AnasaziOperator::Cayley(eigenParams, 
						  solverParams,
						  grp);
    else
       LOCA::ErrorCheck::throwError("LOCA::AnasaziOperator::Manager::reset()",
				    "Invalid Anasazi operator name: " + 
				    opName);
    return NOX::Abstract::Group::Ok;
  }
  else 
    return op->reset(eigenParams, solverParams, grp);
}

const string&
LOCA::AnasaziOperator::Manager::label() const
{
  return op->label();
}

NOX::Abstract::Group::ReturnType 
LOCA::AnasaziOperator::Manager::apply(const NOX::Abstract::Vector& input, 
				      NOX::Abstract::Vector& output) const
{
  return op->apply(input, output);
}

void
LOCA::AnasaziOperator::Manager::transformEigenvalue(double& ev_r, 
						    double& ev_i) const
{
  op->transformEigenvalue(ev_r, ev_i);
}

NOX::Abstract::Group::ReturnType 
LOCA::AnasaziOperator::Manager::rayleighQuotient(
				         const NOX::Abstract::Vector& evec_r,
					 const NOX::Abstract::Vector& evec_i,
					 double& rq_r, double& rq_i) const
{
  return op->rayleighQuotient(evec_r, evec_i, rq_r, rq_i);
}

LOCA::AnasaziOperator::Generic&
LOCA::AnasaziOperator::Manager::getOperator()
{
  return *op;
}

const LOCA::AnasaziOperator::Generic&
LOCA::AnasaziOperator::Manager::getOperator() const
{
  return *op;
}
