// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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

#include "PitchforkProblemInterface.H"

PitchforkProblemInterface::PitchforkProblemInterface(int N, double a, 
						     double b, double l)  : 
  initialGuess(N),
  alpha(a),
  beta(b),
  lambda(l),
  n(N)
{
  for (int i=0; i<n; i++) 
    initialGuess(i) = lambda/alpha + 0.01;
}

const NOX::LAPACK::Vector&
PitchforkProblemInterface::getInitialGuess()
{
  return initialGuess;
}

bool
PitchforkProblemInterface::computeF(NOX::LAPACK::Vector& f, 
				    const NOX::LAPACK::Vector &x)
{
  double h = 2.0 / static_cast<double>(n-1);

  f(0) = 2.0*(x(1) - x(0)) + h*h*source_term(x(0));
  f(n-1) = 2.0*(x(n-2) - x(n-1)) + h*h*source_term(x(n-1));
  for (int i=1; i<n-1; i++)
    f(i) = x(i-1) - 2.0*x(i) + x(i+1) + h*h*source_term(x(i));
  
  return true;
}

bool
PitchforkProblemInterface::computeJacobian(NOX::LAPACK::Matrix<double>& J, 
					   const NOX::LAPACK::Vector & x)
{
  double h = 2.0 / static_cast<double>(n-1);
  
  J(0,0) = -2.0 + h*h*source_deriv(x(0));
  J(0,1) = 2.0;
  J(n-1,n-1) = -2.0 + h*h*source_deriv(x(n-1));
  J(n-1,n-2) = 2.0;
  for (int i=1; i<n-1; i++) {
    J(i,i-1) = 1.0;
    J(i,i+1) = 1.0;
    J(i,i) = -2.0 + h*h*source_deriv(x(i));
  }
  return true;
}

void
PitchforkProblemInterface::setParams(const LOCA::ParameterVector& p) {
  alpha = p.getValue("alpha");
  beta = p.getValue("beta");
  lambda = p.getValue("lambda");
}

double
PitchforkProblemInterface::source_term(double x) {
  return lambda*x - alpha*x*x + beta*x*x*x;
}

double
PitchforkProblemInterface::source_deriv(double x) {
  return lambda - 2.0*alpha*x + 3.0*beta*x*x;
}

void
PitchforkProblemInterface::printSolution(const NOX::LAPACK::Vector &x,
					 const double conParam)
{

   std::cout << "At parameter value: " << conParam << "   the solution vector is\n";

   if (n < 8) {
     for (int i=0; i<n; i++)  std::cout << " " << x(i);
   }
   else {
     for (int i=0; i<6; i++)  std::cout << " " << x(i);
     std::cout << " ...";
     for (int i=n-2; i<n; i++)  std::cout << " " << x(i);
   }
   std::cout << std::endl;

}
