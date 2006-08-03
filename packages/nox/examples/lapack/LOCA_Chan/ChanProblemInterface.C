// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "ChanProblemInterface.H"
#include "LOCA_Parameter_Vector.H"
#include "NOX_LAPACK_Vector.H"
#include "NOX_LAPACK_Matrix.H"
#include "LOCA_GlobalData.H"
#include "NOX_Utils.H"

ChanProblemInterface::ChanProblemInterface(
		    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		    int N, double a, double b, double s)  : 
  globalData(global_data),
  initialGuess(N),
  alpha(a),
  beta(b),
  scale(s),
  n(N),
  outputFilePtr(NULL)
{
  init();
}

ChanProblemInterface::ChanProblemInterface(
		    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		    int N, double a, double b, double s, ofstream& file)  : 
  globalData(global_data),
  initialGuess(N),
  alpha(a),
  beta(b),
  scale(s),
  n(N),
  outputFilePtr(&file)
{
  init();
}

const NOX::LAPACK::Vector&
ChanProblemInterface::getInitialGuess()
{
  return initialGuess;
}

bool
ChanProblemInterface::computeF(NOX::LAPACK::Vector& f, 
			       const NOX::LAPACK::Vector &x)
{
  f(0) = x(0) - beta;
  f(n-1) = x(n-1) - beta;
  for (int i=1; i<n-1; i++)
    f(i) = (x(i-1) - 2*x(i) + x(i+1))*(n-1)*(n-1) 
      + source_param(alpha, scale)*source_term(x(i));
  
  return true;
}

bool
ChanProblemInterface::computeJacobian(NOX::LAPACK::Matrix<double>& J, 
				     const NOX::LAPACK::Vector & x)
{
  J(0,0) = 1.0;
  J(n-1,n-1) = 1.0;
  for (int i=1; i<n-1; i++) {
    J(i,i-1) = (n-1)*(n-1);
    J(i,i+1) = J(i,i-1);
    J(i,i) = -2.*J(i,i-1) + source_param(alpha, scale)*source_deriv(x(i));
  }
  return true;
}

void
ChanProblemInterface::setParams(const LOCA::ParameterVector& p) {
  alpha = p.getValue("alpha");
  beta = p.getValue("beta");
  scale = p.getValue("scale");
}

void
ChanProblemInterface::init() {

  for (int i=0; i<n; i++) 
    initialGuess(i) = 
      i*(n-1-i)*source_param(alpha, scale)/((n-1)*(n-1)) + 0.001;
}

double
ChanProblemInterface::source_term(double x) {
  return 1. + (x + 0.5*x*x)/(1. + 0.01*x*x);
}

double
ChanProblemInterface::source_deriv(double x) {
  double y = 1. + 0.01*x*x;
  return (1. + x - 0.01*x*x)/(y*y);
}

double 
ChanProblemInterface::source_param(double a, double s) {
  double as = a*s;

  //return as / (1.0 + 0.01*as*as);
  return as;
}

void
ChanProblemInterface::printSolution(const NOX::LAPACK::Vector &x,
                                    const double conParam)
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() 
      << "At parameter value: " << conParam 
      << "   the solution vector is\n";

    if (n < 8) {
      for (int i=0; i<n; i++) 
	globalData->locaUtils->out() << " " << x(i);
    }
    else {
      for (int i=0; i<6; i++) 
	globalData->locaUtils->out() << " " << x(i);
       globalData->locaUtils->out() << " ...";
      for (int i=n-2; i<n; i++)  
	globalData->locaUtils->out() << " " << x(i);
    }
    globalData->locaUtils->out() << std::endl;
  }

  if (outputFilePtr != NULL) {
    (*outputFilePtr) << conParam << " ";
    for (int i=0; i<n; i++)
      (*outputFilePtr) << x(i) << " ";
    (*outputFilePtr) << std::endl << std::endl;
  }
}
