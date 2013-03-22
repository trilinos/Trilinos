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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
		  const Teuchos::RCP<std::ostream>& outputStream,
		  const Teuchos::RCP<std::ostream>& errStream) :
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
  printTest=0; // This must be initialized before use in reset(), or else
               // it can take on different random values on different 
               // processors leading to inconsistent program flow. 
               // KL 7 June 2009
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
  using namespace Teuchos;

  // "Output Information" may be stored as a sublist, an int, or as
  // a NOX::Utils::MsgType
  if (p.isSublist("Output Information")) {

    ParameterList& printList = p.sublist("Output Information");
    typedef std::map<std::string, NOX::Utils::MsgType> OptionMap;
    OptionMap output_options;
    output_options["Error"] = NOX::Utils::Error;
    output_options["Warning"] = NOX::Utils::Warning;
    output_options["Outer Iteration"] = NOX::Utils::OuterIteration;
    output_options["Inner Iteration"] = NOX::Utils::InnerIteration;
    output_options["Parameters"] = NOX::Utils::Parameters;
    output_options["Details"] = NOX::Utils::Details;
    output_options["Outer Iteration StatusTest"] = NOX::Utils::OuterIterationStatusTest;
    output_options["Linear Solver Details"] = NOX::Utils::LinearSolverDetails;
    output_options["Test Details"] = NOX::Utils::TestDetails;
    output_options["Stepper Iteration"] = NOX::Utils::StepperIteration;
    output_options["Stepper Details"] = NOX::Utils::StepperDetails;
    output_options["Stepper Parameters"] = NOX::Utils::StepperParameters;
    output_options["Debug"] = NOX::Utils::Debug;

    bool add_test = false;
    OptionMap::const_iterator start = output_options.begin();
    OptionMap::const_iterator stop = output_options.end();
    for (OptionMap::const_iterator i = start; i != stop; ++i) {
      add_test = printList.get(i->first, false);
      if (add_test)
	printTest += i->second;
    }
  }
  else if (isParameterType<NOX::Utils::MsgType>(p, "Output Information"))
    printTest = get<NOX::Utils::MsgType>(p, "Output Information");
  else
    printTest = p.get("Output Information", 0xf);

#ifdef HAVE_MPI
  if (p.isParameter("MyPID"))
    myPID = p.get("MyPID", 0);
  else {
    // We have to check to ensure MPI has been initialized before calling
    // MPI_Comm_rank. Relates to bug 2566. - KL, 17 Sept 2006.
    int mpiIsRunning = 0;
    MPI_Initialized(&mpiIsRunning);
    if (mpiIsRunning)
      {
        MPI_Comm_rank(MPI_COMM_WORLD, &myPID);
      }
    else
      {
        myPID=0;
      }
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

  if (p.INVALID_TEMPLATE_QUALIFIER
      isType< Teuchos::RCP<std::ostream> >("Output Stream"))
    printStream =
      (p).INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RCP<std::ostream> >("Output Stream");
  else 
    printStream = Teuchos::rcp(&(std::cout), false);
  myStream = (myPID == printProc) ? printStream : blackholeStream;
  
  if (p.INVALID_TEMPLATE_QUALIFIER 
      isType< Teuchos::RCP<std::ostream> >("Error Stream"))
    errorStream =
      (p).INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RCP<std::ostream> >("Error Stream");
  else 
    errorStream = Teuchos::rcp(&(std::cerr), false);
}

//static
NOX::Utils::Fill NOX::Utils::fill(int filln, char fillc)
{
  return NOX::Utils::Fill(filln, fillc);
}

std::ostream& NOX::operator<<(std::ostream& os, const NOX::Utils::Fill& f)
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

std::ostream& NOX::operator<<(std::ostream& os, const NOX::Utils::Sci& s)
{
  os.setf(std::ios::scientific);
  os.precision(s.p);
  os << std::setw(s.p + 6) << s.d;
  os.unsetf(std::ios::scientific);
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

void NOX::Utils::print(std::ostream& os) const
{
  os << "NOX::Utils Printing Object" << std::endl;
  os << "Output Information Level = " << printTest << std::endl;
  os << "My PID = " << myPID << std::endl;
  os << "Print Processor = " << printProc << std::endl;
  os << "Precision = " << precision << std::endl;
}

std::ostream& NOX::operator<<(std::ostream& os, const NOX::Utils& utils)
{
  utils.print(os);
  return os;
}
