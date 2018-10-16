// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
//
// Testing the DebugManager object.
//
// Verbosity levels are
//  NO_STATUS,
//  BASIC_STATUS,
//  DETAILED_STATUS,
//  VERBOSE_DETAILED_STATUS
//  NUM_STATUS_OUTPUT_LEVELS
//
//  This test can only really be verified by reading the output.
//  So we are testing that DebugManager doesn't crash.


#include <Zoltan2_DebugManager.hpp>
#include <Zoltan2_Parameters.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_DefaultComm.hpp>

#include <set>
#include <iostream>
#include <string>
#include <ostream>

using Zoltan2::DebugManager;
using Zoltan2::NO_STATUS;
using Zoltan2::BASIC_STATUS;
using Zoltan2::DETAILED_STATUS;
using Zoltan2::VERBOSE_DETAILED_STATUS;

typedef Zoltan2::MessageOutputLevel level_t;

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();
  int nprocs = comm->getSize();
  bool fail = false;

  std::set<string> basicMsgs, detailedMsgs, verboseMsgs; 
  std::set<string>::iterator next;

  std::ostringstream oss;
  oss << "Proc " << rank << ": This is a ";

  basicMsgs.insert(oss.str()+string(" basic message."));
  basicMsgs.insert(oss.str()+string("another basic message."));
  detailedMsgs.insert(oss.str()+string(" detailed message."));
  detailedMsgs.insert(oss.str()+string("another detailed message."));
  verboseMsgs.insert(oss.str()+string(" verbose message."));
  verboseMsgs.insert(oss.str()+string("another verbose message."));

  level_t numLevels = Zoltan2::NUM_STATUS_OUTPUT_LEVELS;
  DebugManager *dm = NULL;

  // all print to std::cout

  bool iPrint = (rank%2 == 0);

  comm->barrier();

  for (int i = 0; i < numLevels; i++){

    level_t level = static_cast<level_t>(i);
  
    try {
      dm = new DebugManager(rank, iPrint, std::cout, level);
    }
    catch(std::exception &e){
      fail=true;
    }

    TEST_FAIL_AND_EXIT(*comm, !fail, "constructor", 1);

    if (rank==0){
      std::cout << "\nThere are " << nprocs << " processes. ";
      std::cout << "Even ranks participate, output level is: " << level << std::endl;
    }

    comm->barrier();

    try{
      for (next=basicMsgs.begin(); next != basicMsgs.end(); ++next){
        dm->print(BASIC_STATUS, *next);
      }
      comm->barrier();
      for (next=detailedMsgs.begin(); next != detailedMsgs.end(); ++next){
        dm->print(DETAILED_STATUS, *next);
      }
      comm->barrier();
      for (next=verboseMsgs.begin(); next != verboseMsgs.end(); ++next){
        dm->print(VERBOSE_DETAILED_STATUS, *next);
      }
      comm->barrier();
    }
    catch(std::exception &e){
      fail=true;
    }

    TEST_FAIL_AND_EXIT(*comm, !fail, "print to standard output", 1);

    delete dm;
  }

  // Node zero prints to a file

  iPrint = (rank == 0);
  comm->barrier();

  for (int i = 0; i < numLevels; i++){

    level_t level = static_cast<level_t>(i);

    std::ios_base::openmode flags = std::ios_base::out & std::ios_base::trunc;

    std::ofstream outF("testFile.txt", flags);
  
    try {
      dm = new DebugManager(rank, iPrint, outF, level);
    }
    catch(std::exception &e){
      fail=true;
    }

    TEST_FAIL_AND_EXIT(*comm, !fail, "constructor", 1);

    if (rank==0){
      std::cout << "\nThere are " << nprocs << " processes. ";
      std::cout << "Rank zero only participates, output level is: ";
      std::cout << level << std::endl;
    }

    try {
      for (next=basicMsgs.begin(); next != basicMsgs.end(); ++next){
        dm->print(BASIC_STATUS, *next);
      }
      comm->barrier();

      for (next=detailedMsgs.begin(); next != detailedMsgs.end(); ++next){
        dm->print(DETAILED_STATUS, *next);
      }
      comm->barrier();

      for (next=verboseMsgs.begin(); next != verboseMsgs.end(); ++next){
        dm->print(VERBOSE_DETAILED_STATUS, *next);
      }
      comm->barrier();
    }
    catch(std::exception &e){
      fail=true;
    }

    delete dm;

    TEST_FAIL_AND_EXIT(*comm, !fail, "print to a file", 1);

    outF.close();

    comm->barrier();

    if (rank == 0){
      std::ifstream inF("testFile.txt");
      string s;
      while (getline(inF, s)){
        std::cout << s << std::endl;
      }
      inF.close();
      system("rm testFile.txt");  // \todo fix for windows
    }

    comm->barrier();
  }

  if (rank==0)
    std::cout << "PASS" << std::endl;
}
