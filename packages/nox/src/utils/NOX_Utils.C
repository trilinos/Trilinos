// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NOX_Utils.H"

int NOX::Utils::myPID = 0;
int NOX::Utils::outputLevel = 2;
int NOX::Utils::printProc = 0;
int NOX::Utils::precision = 4;
string NOX::Utils::stars = "***********************************************************************\n";

void NOX::Utils::setUtils(NOX::Parameter::List& p)
{
  outputLevel = p.getParameter("Output Level", outputLevel);
  myPID = p.getParameter("MyPID", myPID);
  printProc = p.getParameter("Output Processor", printProc);
  precision = p.getParameter("Output Precision", precision);
}

bool NOX::Utils::isPrintProc()
{
  return (printProc == myPID);
}

bool NOX::Utils::doPrint(int printLevel)
{
  return ((printProc == myPID) && (printLevel <= outputLevel));
}

bool NOX::Utils::doAllPrint(int printLevel)
{
  return (printLevel <= outputLevel);
}

int NOX::Utils::getMyPID()
{
  return myPID;
}

