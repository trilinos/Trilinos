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

#include "BrusselatorProblemInterface.H"
#include "LOCA_Parameter_Vector.H"
#include "NOX_LAPACK_Vector.H"
#include "NOX_LAPACK_Matrix.H"

BrusselatorProblemInterface::BrusselatorProblemInterface(int N, double a, 
							 double b, 
							 double d1,
							 double d2,
							 ofstream& file)  : 
  initialGuess(2*N),
  alpha(a),
  beta(b),
  D1(d1),
  D2(d2),
  n(N),
  outputFile(file)
{
  for (int i=0; i<n; i++) {
    //initialGuess(i) = alpha + 0.001;
    //initialGuess(n+i) = beta/alpha + 0.001;
    initialGuess(i) = alpha;
    initialGuess(n+i) = beta/alpha;
  }
}

const NOX::LAPACK::Vector&
BrusselatorProblemInterface::getInitialGuess()
{
  return initialGuess;
}

bool
BrusselatorProblemInterface::computeF(NOX::LAPACK::Vector& f, 
				      const NOX::LAPACK::Vector &x)
{
  double h = 1.0 / double(n-1);
  double hh = h*h;

  f(0) = x(0) - alpha;
  for (int i=1; i<n-1; i++)
    f(i) = D1*(x(i-1) - 2*x(i) + x(i+1)) / hh
      + alpha - (beta + 1.0)*x(i) + x(i)*x(i)*x(i+n);
  f(n-1) = x(n-1) - alpha;

  f(n) = x(n) - beta/alpha;
  for (int i=n+1; i<2*n-1; i++)
    f(i) = D2*(x(i-1) - 2*x(i) + x(i+1)) / hh
      + beta*x(i-n) - x(i-n)*x(i-n)*x(i);
  f(2*n-1) = x(2*n-1) - beta/alpha;
  
  return true;
}

bool
BrusselatorProblemInterface::computeJacobian(NOX::LAPACK::Matrix& J, 
					     const NOX::LAPACK::Vector & x)
{
  double h = 1.0 / double(n-1);
  double hh = h*h;

  J(0,0) = 1.0;
  for (int i=1; i<n-1; i++) {
    J(i,i-1) = D1/hh;
    J(i,i) = -2.0*D1/hh - (beta + 1.0) + 2.0*x(i)*x(i+n);
    J(i,i+1) = D1/hh;
    J(i,i+n) = x(i)*x(i);
  }
  J(n-1,n-1) = 1.0;

  J(n,n) = 1.0;
  for (int i=n+1; i<2*n-1; i++) {
    J(i,i-n) = beta - 2*x(i-n)*x(i);
    J(i,i-1) = D2/hh;
    J(i,i) = -2.0*D2/hh - x(i-n)*x(i-n);
    J(i,i+1) = D2/hh;
  }
  J(2*n-1,2*n-1) = 1.0;

  return true;
}

void
BrusselatorProblemInterface::setParams(const LOCA::ParameterVector& p) {
  alpha = p.getValue("alpha");
  beta = p.getValue("beta");
  D1 = p.getValue("D1");
  D2 = p.getValue("D2");
}

void
BrusselatorProblemInterface::printSolution(const NOX::LAPACK::Vector &x,
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

   if (n < 8) {
     for (int i=n; i<2*n; i++)  cout << " " << x(i);
   }
   else {
     for (int i=n; i<n+6; i++)  cout << " " << x(i);
     cout << " ...";
     for (int i=2*n-2; i<2*n; i++)  cout << " " << x(i);
   }
   cout << endl;

   outputFile << conParam << " ";
   for (int i=0; i<2*n; i++)
     outputFile << x(i) << " ";
   outputFile << endl << endl;

}

bool
BrusselatorProblemInterface::computeMass(NOX::LAPACK::Matrix& M, 
					 const NOX::LAPACK::Vector & x)
{
  for (int i=0; i<2*n; i++)
    M(i,i) = 1.0;

  return true;
}
