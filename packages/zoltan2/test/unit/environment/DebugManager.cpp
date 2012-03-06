// @HEADER
// ***********************************************************************
//                Copyright message goes here.
// ***********************************************************************
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

using namespace std;
using Zoltan2::DebugManager;
using Zoltan2::NO_STATUS;
using Zoltan2::BASIC_STATUS;
using Zoltan2::DETAILED_STATUS;
using Zoltan2::VERBOSE_DETAILED_STATUS;

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();
  int nprocs = comm->getSize();
  bool fail = false;

  set<string> basicMsgs, detailedMsgs, verboseMsgs; 
  set<string>::iterator next;

  ostringstream oss;
  oss << "Proc " << rank << ": This is a ";

  basicMsgs.insert(oss.str()+string(" basic message."));
  basicMsgs.insert(oss.str()+string("another basic message."));
  detailedMsgs.insert(oss.str()+string(" detailed message."));
  detailedMsgs.insert(oss.str()+string("another detailed message."));
  verboseMsgs.insert(oss.str()+string(" verbose message."));
  verboseMsgs.insert(oss.str()+string("another verbose message."));

  int numLevels = Zoltan2::NUM_STATUS_OUTPUT_LEVELS;
  DebugManager *dm = NULL;

  // all print to cout

  bool iPrint = (rank%2 == 0);

  comm->barrier();

  for (int level = 0; level < numLevels; level++){
  
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

  for (int level = 0; level < numLevels; level++){

    ios_base::openmode flags = ios_base::out & ios_base::trunc;

    ofstream outF("testFile.txt", flags);
  
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
      ifstream inF("testFile.txt");
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
