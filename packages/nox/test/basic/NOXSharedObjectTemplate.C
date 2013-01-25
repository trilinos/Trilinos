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

// Tests the NOX's shared object template

#include "NOX.H"
#include "NOX_SharedObjectTemplate.H"
#include "NOX_MeritFunction_SumOfSquares.H"

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

  // Create a print utility
  int outputInfo = NOX::Utils::Error;
  if (verbose)
    outputInfo += NOX::Utils::Debug;
  int printProc = 0;
  Teuchos::RCP<std::ostream> outputstream = 
    Teuchos::rcp(&(std::cout), false);
  Teuchos::RCP<std::ostream> errorstream = 
    Teuchos::rcp(&(std::cerr), false);
  Teuchos::RCP<NOX::Utils> u = 
    Teuchos::rcp(new NOX::Utils(outputInfo, myPID, printProc, 6, outputstream, 
				errorstream));

  // Let's share a merit function between status tests
  Teuchos::RCP<NOX::MeritFunction::Generic> mf = 
    Teuchos::rcp(new NOX::MeritFunction::SumOfSquares(u));
  Teuchos::RCP<NOX::StatusTest::Generic> test1 = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND, 
					    u.get()));
  Teuchos::RCP<NOX::StatusTest::Generic> test2 = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND, 
					    u.get()));
  
  NOX::SharedObject<NOX::MeritFunction::Generic,NOX::StatusTest::Generic> sharedObject(mf);
  
  u->out(NOX::Utils::Debug) 
    << "Testing isOwner() via non-const get(*newOwner)...";

  sharedObject.getObject(test1.get());

  if ( !(sharedObject.isOwner(test1.get())) )
    status +=1;
  
  if ( sharedObject.isOwner(test2.get()) )
    status += 1;

  if (status == 0)
    u->out(NOX::Utils ::Debug) << "passed!" << std::endl;
  else 
    u->out(NOX::Utils ::Debug) << "failed!" << std::endl;
    
  u->out(NOX::Utils::Debug) << "Testing isOwner() via const get()...";

  sharedObject.getObject(test2.get());  // change ownership to test 2
  Teuchos::RCP<const NOX::MeritFunction::Generic> tmpPtr = 
    sharedObject.getObject();

  if ( !(sharedObject.isOwner(test2.get())) )
    status +=1;
  
  if ( sharedObject.isOwner(test1.get()) )
    status += 1;

  if (status == 0)
    u->out(NOX::Utils ::Debug) << "passed!" << std::endl;
  else 
    u->out(NOX::Utils ::Debug) << "failed!" << std::endl;
    
  

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
