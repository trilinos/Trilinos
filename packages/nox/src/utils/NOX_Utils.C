// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NOX_Utils.H"
#include <iomanip>

int NOX::Utils::myPID = 0;
int NOX::Utils::outputLevel = 2;
int NOX::Utils::printProc = 0;
int NOX::Utils::precision = 4;

NOX::Utils::Fill NOX::Utils::fillobj;
NOX::Utils::Sci NOX::Utils::sciobj;

NOX::Utils::Fill& NOX::Utils::fill(int filln, char fillc) 
{
  fillobj.n = filln; 
  fillobj.c = fillc;
  return fillobj;
}

ostream& operator<<(ostream& os, const NOX::Utils::Fill& f)
{
  for (int i = 0; i < f.n; i ++)
    os << f.c;
  return os;
}

NOX::Utils::Sci& NOX::Utils::sci(double dval)
{
  sciobj.d = dval;
}

ostream& operator<<(ostream& os, const NOX::Utils::Sci& s)
{
  os.setf(ios::scientific);
  os.precision(NOX::Utils::precision);
  os << setw(NOX::Utils::precision + 6) << s.d;
  cout.unsetf(ios::scientific);
  return os;
}


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

