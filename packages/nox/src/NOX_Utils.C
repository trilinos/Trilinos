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
#include "NOX_Utils.H"
#include "NOX_Parameter_List.H"

NOX::Utils::Utils()
{
  // Make an empty list
  NOX::Parameter::List p;
  reset(p);
}

NOX::Utils::Utils(NOX::Parameter::List& p)
{
  reset(p);
}

NOX::Utils::~Utils()
{
}

void NOX::Utils::reset(NOX::Parameter::List& p)
{
  printTest = p.getParameter("Output Information", 0xf);
  myPID = p.getParameter("MyPID", 0);
  printProc = p.getParameter("Output Processor", 0);
  precision = p.getParameter("Output Precision", 3);
}

//static
NOX::Utils::Fill NOX::Utils::fill(int filln, char fillc)
{
  return NOX::Utils::Fill(filln, fillc);
}

ostream& operator<<(ostream& os, const NOX::Utils::Fill& f)
{
  for (int i = 0; i < f.n; i ++)
    os << f.c;
  return os;
}

NOX::Utils::Sci NOX::Utils::sciformat(double dval, int p) const
{
  return NOX::Utils::Sci(dval, ((p > 0) ? p : precision) );
}

ostream& operator<<(ostream& os, const NOX::Utils::Sci& s)
{
  os.setf(ios::scientific);
  os.precision(s.p);
  os << setw(s.p + 6) << s.d;
  cout.unsetf(ios::scientific);
  return os;
}

bool NOX::Utils::isPrintProcess() const
{
  return (printProc == myPID);
}

bool NOX::Utils::isPrintProcessAndType(MsgType type) const
{
  return (isPrintProcess() && isPrintType(type));
}

bool NOX::Utils::isPrintType(MsgType type) const
{
  return ((type == NOX::Utils::Error) || ((printTest & type) != 0));
}

//private
void NOX::Utils::deprecated(const string& oldname, const string& newname)
{
#ifdef WITH_PRERELEASE
  cerr << "WARNING: NOX::Utils::" << oldname << " is deprecated!" << "\n"
       << "         Use Nox::Utils::" << newname << " instead." << endl;
#endif
}

//static & deprecated
NOX::Utils::Sci NOX::Utils::sci(double dval, int p) 
{
  deprecated("sci","sciformat");
  return  NOX::Utils::Sci(dval, ((p > 0) ? p : 3) );
}

//static & deprecated
bool NOX::Utils::isPrintProc() 
{
  deprecated("isPrintProc","isPrintProcess");
  return true;
}

//static & deprecated
bool NOX::Utils::doPrint(MsgType type) 
{
  deprecated("doPrint","isPrintProcessAndType");
  return true;
}

//static & deprecated
bool NOX::Utils::doAllPrint(MsgType type) 
{
  deprecated("doAllPrint","isPrintType");
  return true;
}

//static & deprecated
bool NOX::Utils::doPrint(int printLevel)
{
  deprecated("doPrint","isPrintProcessAndType");
  return (isPrintProc());
}

//static & deprecated
bool NOX::Utils::doAllPrint(int printLevel)
{
  deprecated("doAllPrint","isPrintType");
  return true;
}

//static & deprecated
void NOX::Utils::setUtils(NOX::Parameter::List& params)
{
  deprecated("setUtils","reset");
}
