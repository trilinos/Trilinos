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
    cout << "Test #1: ctor 1" << endl;
    cout << "***************" << endl;
    cout << "\nBuilding utils1 using ctor #1...";
  }
  Teuchos::RefCountPtr<std::ostream> outputstream = 
    Teuchos::rcp(&(std::cout), false);
  Teuchos::RefCountPtr<std::ostream> errorstream = 
    Teuchos::rcp(&(std::cerr), false);
  NOX::Utils utils1(outputInfo, myPID, printProc, 6, outputstream, 
		    errorstream);
  if (myPID == printProc)
    cout << "Done!\n" << endl;

  cout << utils1 << endl;
  cout.flush();

#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (myPID == printProc) {
    cout << "Test #2: ctor 2 with plist" << endl;
    cout << "**************************" << endl;
    cout << "Building utils2 using ctor #2...";
  }
  Teuchos::ParameterList p;
  p.set("Output Information", outputInfo);
  p.set("MyPID", myPID);
  p.set("OutputProcessor", printProc);
  p.set("Output Precision", 3);
  p.set("Output Stream", outputstream);
  p.set("Error Stream", outputstream);
  NOX::Utils utils2(p);
  if (myPID == printProc)
    cout << "Done!" << endl;

  // even if not verbose, lets test the print capability
  cout << utils2 << endl;
  cout.flush();
  
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // Test Copy ctor
  if (myPID == printProc) {
    cout << "Test #3: Testing copy ctor" << endl;
    cout << "**************************" << endl;
  }
  NOX::Utils utils3(utils1);
  
  if (myPID == printProc)
    cout << "Printing utils3 - should be copy of utils1" << endl; 
  cout << utils3 << endl;
  cout.flush();

#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (myPID == printProc) {
    cout << "Test #4: Testing reset with plist" << endl;
    cout << "*********************************" << endl;
    cout << "Testing reset capability on utils3 - should be copy of utils2" 
	 << endl;
  }
  utils3.reset(p);

  cout << utils3 << endl;
  cout.flush();

#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // Test the ostreams
  if (myPID == printProc) {
    cout << "Test #5: Testing the out() function" << endl;
    cout << "***********************************" << endl;
  }
  std::string printtest = "Function out() works correctly";
  if (myPID != printProc)
    printtest = "Test Failed!";
  utils1.out() << printtest << endl;
  utils1.out().flush();
  
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (myPID == printProc) {
    cout << "\nTest #6: Testing the out(MsgType) function:" << endl;
    cout << "******************************************" << endl;
  }
  printtest = "Function out(MsgType) works correctly";
  if (myPID != printProc)
    printtest = "Test Failed!";
  utils1.out(NOX::Utils::OuterIteration) << printtest << endl;
  utils1.out(NOX::Utils::InnerIteration) << printtest << endl;
  utils1.out().flush();
  
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (myPID == printProc) {
    cout << "\nTest #7: Testing the pout() function:" << endl;
    cout << "*************************************" << endl;
  }
  utils1.pout() << "MyPID = " << myPID
		<< "  responded!" << endl;
  utils1.pout().flush();

#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (myPID == printProc) {
    cout << "\nTest #8: Testing the pout(MsgType) function:" << endl;
    cout << "******************************************" << endl;
  }
  utils1.pout(NOX::Utils::OuterIteration) << "Test Failed!" << endl;
  utils1.pout(NOX::Utils::InnerIteration) << "MyPID = " << myPID
					  << "  responded!" << endl;
  utils1.pout().flush();

#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // Test the ostreams
  if (myPID == printProc) {
    cout << "Test #9: Testing the err() function" << endl;
    cout << "***********************************" << endl;
  }
  printtest = "Function err() works correctly";
  if (myPID != printProc)
    printtest = "Test Failed!";
  utils1.err() << printtest << endl;
  
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (myPID == printProc) {
    cout << "\nTest #10: Testing the perr() function:" << endl;
    cout << "*************************************" << endl;
  }
  utils1.perr() << "MyPID = " << myPID
		<< "  responded!" << endl;

#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (myPID == printProc) {
    cout << "\nTest #11: Testing output to a file:" << endl;
    cout << "******************************************" << endl;
  
    Teuchos::RefCountPtr<std::ostream> fileStream = 
      Teuchos::rcp(new std::ofstream("testfile.out"));
    p.set("Output Stream", fileStream);
    NOX::Utils utils4(p);
    
    std::string testLine = "Supercalifragilisticexpialidocious";
    utils4.out() << testLine << endl;
    
    std::ifstream inStream("testfile.out");
    std::string read;
    inStream >> read;
    
    utils1.out() << "Text output to file: " << testLine << endl;
    utils1.out() << "Text read from file: " << read << endl;
  
    if (!testLine.compare(read)) {
      utils1.out() << "File output works correctly!" << endl;
    }
    else
      status = 1;
  }
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (myPID == printProc) {
    if (status == 0) 
      std::cout << "\nTest passed!" << endl;
    else
      std::cout << "\nTest Failed!" << endl;
  }

  // Final return value (0 = successfull, non-zero = failure)
  return status;
}
