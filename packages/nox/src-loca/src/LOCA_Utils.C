// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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

#include "NOX_Common.H"

#include "LOCA_Utils.H"  // class definition

#include "NOX_Parameter_List.H"

int LOCA::Utils::precision = 3;
int LOCA::Utils::myPID = 0;
int LOCA::Utils::printProc = 0;
int LOCA::Utils::printTest = 0xf;
LOCA::Utils::SublistMap LOCA::Utils::sublistMap;

LOCA::Utils::Fill LOCA::Utils::fill(int filln, char fillc) 
{
  return LOCA::Utils::Fill(filln, fillc);
}

ostream& operator<<(ostream& os, const LOCA::Utils::Fill& f)
{
  for (int i = 0; i < f.n; i ++)
    os << f.c;
  return os;
}

LOCA::Utils::Sci LOCA::Utils::sci(double dval, int prec)
{
  return LOCA::Utils::Sci(dval, prec);
}

ostream& operator<<(ostream& os, const LOCA::Utils::Sci& s)
{
  os.setf(ios::scientific);
  if (s.p < 0) {
    os.precision(LOCA::Utils::precision);
    os << setw(LOCA::Utils::precision + 6) << s.d;
  }
  else {
    os.precision(s.p);
    os << setw(s.p + 6) << s.d;
  } 
  cout.unsetf(ios::scientific);
  return os;
}

void LOCA::Utils::setUtils(NOX::Parameter::List& p)
{

  initializeSublistMap(p);

  NOX::Parameter::List& utilParams = getSublist("Utilities");

  printTest = utilParams.getParameter("Output Information", printTest);
  myPID = utilParams.getParameter("MyPID", myPID);
  printProc = utilParams.getParameter("Output Processor", printProc);
  precision = utilParams.getParameter("Output Precision", precision);
}

bool LOCA::Utils::isPrintProc()
{
  return (printProc == myPID);
}

bool LOCA::Utils::doPrint(MsgType type)
{
  if (type == Error)
    return isPrintProc();

  return (isPrintProc() && ((printTest & type) != 0));
}

bool LOCA::Utils::doAllPrint(MsgType type)
{
  if (type == Error)
    return true;

  return ((printTest & type) != 0);
}

int LOCA::Utils::getMyPID()
{
  return myPID;
}

void LOCA::Utils::initializeSublistMap(NOX::Parameter::List& p) {

  // Top level sublist
  sublistMap["Top Level"] = &p;

  // LOCA sublist
  NOX::Parameter::List& locaSublist = p.sublist("LOCA");
  sublistMap["LOCA"] = &locaSublist;

  // Stepper sublist
  NOX::Parameter::List& stepperSublist = locaSublist.sublist("Stepper");
  sublistMap["Stepper"] = &stepperSublist;

  // Anasazi sublist
  NOX::Parameter::List& anasaziSublist = stepperSublist.sublist("Anasazi");
  sublistMap["Anasazi"] = &anasaziSublist;

  // Bifurcation sublist
  NOX::Parameter::List& bifurcationSublist = 
    locaSublist.sublist("Bifurcation");
  sublistMap["Bifurcation"] = &bifurcationSublist;

  // Predictor sublist
  NOX::Parameter::List& predictorSublist = locaSublist.sublist("Predictor");
  sublistMap["Predictor"] = &predictorSublist;

  // First Step Predictor sublist
  NOX::Parameter::List& fspredictorSublist = 
    predictorSublist.sublist("First Step Predictor");
  sublistMap["First Step Predictor"] = &fspredictorSublist;

  // Last Step Predictor sublist
  NOX::Parameter::List& lspredictorSublist = 
    predictorSublist.sublist("Last Step Predictor");
  sublistMap["Last Step Predictor"] = &lspredictorSublist;

  // Stepsize sublist
  NOX::Parameter::List& stepsizeSublist = locaSublist.sublist("Step Size");
  sublistMap["Step Size"] = &stepsizeSublist;

  // Utilities sublist
  NOX::Parameter::List& utilitiesSublist = locaSublist.sublist("Utilities");
  sublistMap["Utilities"] = &utilitiesSublist;

  // NOX sublist
  NOX::Parameter::List& noxSublist = p.sublist("NOX");
  sublistMap["NOX"] = &noxSublist;

  // Direction sublist
  NOX::Parameter::List& directionSublist = noxSublist.sublist("Direction");
  sublistMap["Direction"] = &directionSublist;

  // Newton sublist
  NOX::Parameter::List& newtonSublist = directionSublist.sublist("Newton");
  sublistMap["Newton"] = &newtonSublist;

  // Linear Solver sublist
  NOX::Parameter::List& lsSublist = newtonSublist.sublist("Linear Solver");
  sublistMap["Linear Solver"] = &lsSublist;

  // Line Search sublist
  NOX::Parameter::List& lineSearchSublist = noxSublist.sublist("Line Search");
  sublistMap["Line Search"] = &lineSearchSublist;

  // Printing sublist
  NOX::Parameter::List& printingSublist = noxSublist.sublist("Printing");
  sublistMap["Printing"] = &printingSublist;
}

NOX::Parameter::List&
LOCA::Utils::getSublist(const string& name) {

  // Find name in list, if it exists.
  SublistMapIterator i = sublistMap.find(name);

  // If it does exist return the sublist.
  // Otherwise, throw an error.
  if (i != sublistMap.end()) 
    return *((*i).second);
  else {
    cerr << "ERROR: Sublist " << name << " is not a valid sublist." << endl;
    throw "NOX Error";
  }
}

