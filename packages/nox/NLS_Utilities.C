
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

void NLS_Utilities::setUtilities(NLS_ParameterList& p)
{
  outputLevel = p.getParameter("Output Level", 2);
  myPID = p.getParameter("MyPID", 0);
  printProc = p.getParameter("Print Processor", 0);

  if (doPrint(4))
    cout << "Output Level " << outputLevel << "." << endl; 

  if (doPrint(5)) 
    cout << "NLS_Utilities: Processor " << myPID 
	 << " is online." << endl;  
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

