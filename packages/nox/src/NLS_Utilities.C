
// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_Utilities.H"

int NLS_Utilities::myPID = 0;
int NLS_Utilities::outputLevel = 2;
int NLS_Utilities::printProc = 0;
int NLS_Utilities::precision = 4;
string NLS_Utilities::stars = "***********************************************************************\n";

void NLS_Utilities::setUtilities(NLS_ParameterList& p)
{
  outputLevel = p.getParameter("Output Level", outputLevel);
  myPID = p.getParameter("MyPID", myPID);
  printProc = p.getParameter("Output Processor", printProc);
  precision = p.getParameter("Output Precision", precision);
}

bool NLS_Utilities::isPrintProc()
{
  return (printProc == myPID);
}

bool NLS_Utilities::doPrint(int printLevel)
{
  return ((printProc == myPID) && (printLevel <= outputLevel));
}

bool NLS_Utilities::doAllPrint(int printLevel)
{
  return (printLevel <= outputLevel);
}

int NLS_Utilities::getMyPID()
{
  return myPID;
}

