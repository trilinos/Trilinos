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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "ChanProblemInterface.H"

ChanProblemInterface::ChanProblemInterface(int N, double a, double b, 
					   double s)  : 
  initialGuess(N),
  alpha(a),
  beta(b),
  scale(s),
  n(N)
{
  for (int i=0; i<n; i++) 
    initialGuess(i) = 
      i*(n-1-i)*source_param(alpha, scale)/((n-1)*(n-1)) + 0.001;
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
ChanProblemInterface::computeJacobian(NOX::LAPACK::Matrix& J, 
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

   cout << "At parameter value: " << conParam << "   the solution vector is\n";

   if (n < 8) {
     for (int i=0; i<n; i++)  cout << " " << x(i);
   }
   else {
     for (int i=0; i<6; i++)  cout << " " << x(i);
     cout << " ...";
     for (int i=n-2; i<n; i++)  cout << " " << x(i);
   }
   cout << endl;

}
