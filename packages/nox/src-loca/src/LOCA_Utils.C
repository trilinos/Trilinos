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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
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
  printTest = p.getParameter("Output Information", printTest);
  myPID = p.getParameter("MyPID", myPID);
  printProc = p.getParameter("Output Processor", printProc);
  precision = p.getParameter("Output Precision", precision);
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

