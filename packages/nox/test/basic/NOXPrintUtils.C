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

// Tests the NOX's Printing object NOX::Utils 

#include "NOX.H"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

using namespace std;

int main(int argc, char *argv[])
{

  int status = 0;

  int myPID = 0;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myPID);
#endif

  // Check verbosity level
  bool verbose = false;
  if (argc > 1)
    if (argv[1][0]=='-' && argv[1][1]=='v')
      verbose = true;

  int outputInfo = NOX::Utils::Error + NOX::Utils::Warning + 
    NOX::Utils::InnerIteration;
  int printProc = 0;

  if (myPID == printProc) {
    std::cout << "Test #1: ctor 1" << std::endl;
    std::cout << "***************" << std::endl;
    std::cout << "\nBuilding utils1 using ctor #1...";
  }
  Teuchos::RCP<std::ostream> outputstream = 
    Teuchos::rcp(&(std::cout), false);
  Teuchos::RCP<std::ostream> errorstream = 
    Teuchos::rcp(&(std::cerr), false);
  NOX::Utils utils1(outputInfo, myPID, printProc, 6, outputstream, 
		    errorstream);
  if (myPID == printProc)
    std::cout << "Done!\n" << std::endl;

  std::cout << utils1 << std::endl;
  std::cout.flush();

#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (myPID == printProc) {
    std::cout << "Test #2: ctor 2 with plist" << std::endl;
    std::cout << "**************************" << std::endl;
    std::cout << "Building utils2 using ctor #2...";
  }
  Teuchos::ParameterList p;
  p.set("Output Information", outputInfo);
  p.set("MyPID", myPID);
  p.set("Output Processor", printProc);
  p.set("Output Precision", 3);
  p.set("Output Stream", outputstream);
  p.set("Error Stream", outputstream);
  NOX::Utils utils2(p);
  if (myPID == printProc)
    std::cout << "Done!" << std::endl;

  // even if not verbose, lets test the print capability
  std::cout << utils2 << std::endl;
  std::cout.flush();
  
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // Test Copy ctor
  if (myPID == printProc) {
    std::cout << "Test #3: Testing copy ctor" << std::endl;
    std::cout << "**************************" << std::endl;
  }
  NOX::Utils utils3(utils1);
  
  if (myPID == printProc)
    std::cout << "Printing utils3 - should be copy of utils1" << std::endl; 
  std::cout << utils3 << std::endl;
  std::cout.flush();

#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (myPID == printProc) {
    std::cout << "Test #4: Testing reset with plist" << std::endl;
    std::cout << "*********************************" << std::endl;
    std::cout << "Testing reset capability on utils3 - should be copy of utils2" 
	 << std::endl;
  }
  utils3.reset(p);

  std::cout << utils3 << std::endl;
  std::cout.flush();

#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // Test the std::ostreams
  if (myPID == printProc) {
    std::cout << "Test #5: Testing the out() function" << std::endl;
    std::cout << "***********************************" << std::endl;
  }
  std::string printtest = "Function out() works correctly";
  if (myPID != printProc)
    printtest = "Test Failed!";
  utils1.out() << printtest << std::endl;
  utils1.out().flush();
  
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (myPID == printProc) {
    std::cout << "\nTest #6: Testing the out(MsgType) function:" << std::endl;
    std::cout << "******************************************" << std::endl;
  }
  printtest = "Function out(MsgType) works correctly";
  if (myPID != printProc)
    printtest = "Test Failed!";
  utils1.out(NOX::Utils::OuterIteration) << printtest << std::endl;
  utils1.out(NOX::Utils::InnerIteration) << printtest << std::endl;
  utils1.out().flush();
  
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (myPID == printProc) {
    std::cout << "\nTest #7: Testing the pout() function:" << std::endl;
    std::cout << "*************************************" << std::endl;
  }
  utils1.pout() << "MyPID = " << myPID
		<< "  responded!" << std::endl;
  utils1.pout().flush();

#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (myPID == printProc) {
    std::cout << "\nTest #8: Testing the pout(MsgType) function:" << std::endl;
    std::cout << "******************************************" << std::endl;
  }
  utils1.pout(NOX::Utils::OuterIteration) << "Test Failed!" << std::endl;
  utils1.pout(NOX::Utils::InnerIteration) << "MyPID = " << myPID
					  << "  responded!" << std::endl;
  utils1.pout().flush();

#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // Test the std::ostreams
  if (myPID == printProc) {
    std::cout << "Test #9: Testing the err() function" << std::endl;
    std::cout << "***********************************" << std::endl;
  }
  printtest = "Function err() works correctly";
  if (myPID != printProc)
    printtest = "Test Failed!";
  utils1.err() << printtest << std::endl;
  
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (myPID == printProc) {
    std::cout << "\nTest #10: Testing the perr() function:" << std::endl;
    std::cout << "*************************************" << std::endl;
  }
  utils1.perr() << "MyPID = " << myPID
		<< "  responded!" << std::endl;

#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (myPID == printProc) {
    std::cout << "\nTest #11: Testing output to a file:" << std::endl;
    std::cout << "******************************************" << std::endl;
  
    Teuchos::RCP<std::ostream> fileStream = 
      Teuchos::rcp(new std::ofstream("testfile.out"));
    p.set("Output Stream", fileStream);
    NOX::Utils utils4(p);
    
    std::string testLine = "Supercalifragilisticexpialidocious";
    utils4.out() << testLine << std::endl;
    
    std::ifstream inStream("testfile.out");
    std::string read;
    inStream >> read;
    
    utils1.out() << "Text output to file: " << testLine << std::endl;
    utils1.out() << "Text read from file: " << read << std::endl;
  
    if (!testLine.compare(read)) {
      utils1.out() << "File output works correctly!" << std::endl;
    }
    else
      status = 1;
  }
  
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (myPID == printProc) {
    std::cout << "\nTest #12: ctor 2 with plist for Output Information" << std::endl;
    std::cout << "**************************" << std::endl;
    std::cout << "Building utils2 using ctor #2...";
  }
  {
    Teuchos::ParameterList pp;
    Teuchos::ParameterList& info = pp.sublist("Output Information");
    info.set("Error", true);
    info.set("Details", true);
    info.set("Outer Iteration StatusTest", true);
    info.set("Outer Iteration", true);
    pp.set("MyPID", myPID);
    pp.set("Output Processor", printProc);
    pp.set("Output Precision", 3);
    pp.set("Output Stream", outputstream);
    pp.set("Error Stream", outputstream);
    NOX::Utils utils3(pp);
    pp.print(std::cout);
    std::cout << "\n" << utils3 << std::endl;
    if (myPID == printProc)
      std::cout << "Done!" << std::endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (myPID == printProc) {
    if (status == 0) 
      std::cout << "\nTest passed!" << std::endl;
    else
      std::cout << "\nTest Failed!" << std::endl;
  }

  // Final return value (0 = successfull, non-zero = failure)
  return status;
}
