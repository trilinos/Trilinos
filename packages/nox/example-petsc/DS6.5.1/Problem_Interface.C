// ----------   Includes   ----------
#include <iostream>
#include "Problem_Interface.H"

// ----------   User Defined Includes   ----------
#include "DennisSchnabel.H"

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(DennisSchnabel& Problem) :
  problem(Problem)
{ }

Problem_Interface::~Problem_Interface()
{ }

bool Problem_Interface::computeF(const Vec& x, Vec& RHS)
{
  return problem.evaluate(RHS_ONLY, &x, &RHS, NULL);
}

bool Problem_Interface::computeJacobian(const Vec& x, Mat& Jac)
{
  return problem.evaluate(MATRIX_ONLY, &x, NULL, &Jac);
}

bool Problem_Interface::computePreconditioner(Mat& M)
{
  cout << "ERROR: Problem_Interface::computePreconditioner() - Use Explicit Jaciban only for this test problem!" << endl;
  throw;
}
bool Problem_Interface::preconditionVector(Vec& y)
{
  cout << "ERROR: Problem_Interface::preconditionVector() - Use Explicit Jaciban only for this test problem!" << endl;
  throw;
}
//-----------------------------------------------------------------------------

