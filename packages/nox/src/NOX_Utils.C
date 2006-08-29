// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "Teuchos_ParameterList.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

NOX::Utils::Utils(int outputInformation, int MyPID, int outputProcess, 
		  int outputPrecision, 
		  const Teuchos::RefCountPtr<std::ostream>& outputStream,
		  const Teuchos::RefCountPtr<std::ostream>& errStream) :
  precision(outputPrecision),
  myPID(MyPID),
  printTest(outputInformation),
  printProc(outputProcess),
  blackholeStream(Teuchos::rcp(new Teuchos::oblackholestream)),
  printStream(outputStream),
  myStream(),
  errorStream(errStream)
{
  if (printStream == Teuchos::null)
    printStream = Teuchos::rcp(&(std::cout), false);
  if (errorStream == Teuchos::null)
    errorStream = Teuchos::rcp(&(std::cerr), false);
  if (myPID == printProc)
    myStream = printStream;
  else
    myStream = blackholeStream;
}

NOX::Utils::Utils(Teuchos::ParameterList& p)
{
  this->reset(p);
}

NOX::Utils::Utils(const NOX::Utils& source)
{
  *this = source;
}

NOX::Utils::~Utils()
{
}

NOX::Utils& NOX::Utils::operator=(const NOX::Utils& source)
{
  printTest = source.printTest;
  myPID = source.myPID;
  printProc = source.printProc;
  precision = source.precision;
  blackholeStream = source.blackholeStream;
  printStream = source.printStream;
  myStream = source.myStream;
  errorStream = source.errorStream;

  return *this;
}

void NOX::Utils::reset(Teuchos::ParameterList& p)
{
  // Basic info

  // "Output Information" may be stored as an int or as NOX::Utils::MsgType
  if (p.isType<NOX::Utils::MsgType>("Output Information"))
    printTest = p.INVALID_TEMPLATE_QUALIFIER
      get<NOX::Utils::MsgType>("Output Information");
  else
    printTest = p.get("Output Information", 0xf);

#ifdef HAVE_MPI
  if (p.isParameter("MyPID"))
    myPID = p.get("MyPID", 0);
  else {
    MPI_Comm_rank(MPI_COMM_WORLD, &myPID);
    // Set the default in the parameter list
    p.get("MyPID", myPID);
  }
#else
  myPID = p.get("MyPID", 0);
#endif
  printProc = p.get("Output Processor", 0);
  precision = p.get("Output Precision", 3);
  
  // Output streams
  blackholeStream = Teuchos::rcp(new Teuchos::oblackholestream);

  if (p.isType< Teuchos::RefCountPtr<std::ostream> >("Output Stream"))
    printStream =
      (p).INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RefCountPtr<std::ostream> >("Output Stream");
  else 
    printStream = Teuchos::rcp(&(std::cout), false);
  myStream = (myPID == printProc) ? printStream : blackholeStream;
  
  if (p.isType< Teuchos::RefCountPtr<std::ostream> >("Error Stream"))
    errorStream =
      (p).INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RefCountPtr<std::ostream> >("Error Stream");
  else 
    errorStream = Teuchos::rcp(&(std::cerr), false);
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

NOX::Utils::Sci NOX::Utils::sciformat(double dval) const
{
  return NOX::Utils::Sci(dval, precision);
}

NOX::Utils::Sci NOX::Utils::sciformat(double dval, int p) 
{
  return NOX::Utils::Sci(dval, p);
}

ostream& operator<<(ostream& os, const NOX::Utils::Sci& s)
{
  os.setf(ios::scientific);
  os.precision(s.p);
  os << setw(s.p + 6) << s.d;
  os.unsetf(ios::scientific);
  return os;
}

bool NOX::Utils::isPrintType(MsgType type) const
{
  return ((type == NOX::Utils::Error) || ((printTest & type) != 0));
}


std::ostream& NOX::Utils::out() const
{
  return *myStream;
}

std::ostream& NOX::Utils::out(NOX::Utils::MsgType type) const
{
  if (isPrintType(type))
    return *myStream;

  return *blackholeStream;
}

std::ostream& NOX::Utils::pout() const
{
  return *printStream;
}

std::ostream& NOX::Utils::pout(NOX::Utils::MsgType type) const
{
  if (isPrintType(type))
    return *printStream;
  
  return *blackholeStream;
}

std::ostream& NOX::Utils::err() const
{
  return (myPID == printProc) ? *errorStream : *blackholeStream;
}

std::ostream& NOX::Utils::perr() const
{
  return *errorStream;
}

void NOX::Utils::print(ostream& os) const
{
  os << "NOX::Utils Printing Object" << endl;
  os << "Output Information Level = " << printTest << endl;
  os << "My PID = " << myPID << endl;
  os << "Print Processor = " << printProc << endl;
  os << "Precision = " << precision << endl;
}

ostream& operator<<(ostream& os, const NOX::Utils& utils)
{
  utils.print(os);
  return os;
}
