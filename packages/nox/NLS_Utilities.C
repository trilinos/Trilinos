
// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_Utilities.H"

NLS_Utilities::NLS_Utilities() 
{

}

NLS_Utilities::~NLS_Utilities() 
{
  
}

void NLS_Utilities::setUtilities(NLS_ParameterList& p)
{
  outputLevel = p.getParameter("Output Level", 2);
  myPID = p.getParameter("MyPID",0);
  printProc = p.getParameter("Print Processor", 0);

  if (isPrintProc() && isOutput(4)) 
    cout << "Output Level " << outputLevel << "." << endl; 

  if (isOutput(5)) cout << "NLS_Utilities: Processor " << myPID << 
		     " is online." << endl;  

}

bool NLS_Utilities::isPrintProc()
{
  if (printProc == myPID) return true;
  return false;
}

bool NLS_Utilities::isOutput(int printLevel)
{
  if (printLevel <= outputLevel) return true;
  return false;
}

int NLS_Utilities::getMyPID()
{
  return myPID;
}

