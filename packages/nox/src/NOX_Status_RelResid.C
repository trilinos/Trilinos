// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NOX_Status_RelResid.H"

#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"

using namespace NOX::Status;

RelResid::RelResid(double initnorm, double tolerance)
{
  tol = tolerance;
  normrhszero = initnorm;
  reltol = tol * normrhszero;
}

RelResid::~RelResid()
{
}

StatusType RelResid::operator()(const Solver::Generic& problem) const
{
  const Abstract::Group& tmp = problem.getSolutionGroup();
  double normrhs = tmp.getNormRHS();
  if (normrhs < reltol)
    return Converged;
  else
    return Unconverged;
}

ostream& RelResid::print(ostream& stream, int indent = 0) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << "Relative Residual Norm with Tolerance = " << tol << endl;
  return stream;
}
