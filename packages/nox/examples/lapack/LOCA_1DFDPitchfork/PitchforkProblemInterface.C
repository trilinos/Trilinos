// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
