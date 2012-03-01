// @HEADER
// ***********************************************************************
//                Copyright message goes here.
// ***********************************************************************
//
// Testing the MetricOutputManager object.
//
// Verbosity levels are
//  NO_STATUS,
//  BASIC_STATUS,
//  DETAILED_STATUS,
//  VERBOSE_DETAILED_STATUS
//  NUM_STATUS_OUTPUT_LEVELS
//
//  This test can only really be verified by reading the output.
//  So we are testing that MetricOutputManager doesn't crash.


#include <Zoltan2_MetricOutputManager.hpp>
#include <Zoltan2_Parameters.hpp>
#include <ErrorHandlingForTests.hpp>

#include <Teuchos_DefaultComm.hpp>

#include <set>
#include <iostream>
#include <string>
#include <ostream>

using namespace std;
using std::string;
using Zoltan2::MetricOutputManager;

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();
  int nprocs = comm->getSize();
  bool fail = false;

  MetricOutputManager<int> *intmom = NULL;
  MetricOutputManager<float> *floatmom = NULL;
  MetricOutputManager<double> *doublemom = NULL;

  // Even ranks print to cout

  bool iPrint = (rank%2 == 0);
  bool someOnePrints = true;

  comm->barrier();

  try {
    intmom = new MetricOutputManager<int>(
      rank, iPrint, std::cout, someOnePrints);
  }
  catch(std::exception &e){
    fail=true;
  }

  TEST_FAIL_AND_EXIT(*comm, !fail, "constructor", 1);

  if (intmom->getMetricsOn() != true)
    fail = true; 

  TEST_FAIL_AND_EXIT(*comm, !fail, "getMetricsOn", 1);

  if (rank==0){
    std::cout << "\nThere are " << nprocs << " processes. ";
    std::cout << "Even ranks only participate." << std::endl;
  }

  try{
    intmom->print(string("number of things"), string(), 10); 
    intmom->print(string("number of other things"), string(), 5); 
  }
  catch(std::exception &e){
    fail=true;
  }

  TEST_FAIL_AND_EXIT(*comm, !fail, "print to standard output", 1);

  delete intmom;

  // All print to cout

  iPrint = true;
  someOnePrints = true;
  comm->barrier();

  try {
    floatmom = new MetricOutputManager<float>(

      rank, iPrint, std::cout, someOnePrints);
  }
  catch(std::exception &e){
    fail=true;
  }

  TEST_FAIL_AND_EXIT(*comm, !fail, "constructor", 1);

  if (floatmom->getMetricsOn() != true)
    fail = true; 

  TEST_FAIL_AND_EXIT(*comm, !fail, "getMetricsOn", 1);

  if (rank==0){
    std::cout << "\nThere are " << nprocs << " processes. ";
    std::cout << "All ranks participate." << std::endl;
  }

  try{
    floatmom->print(string("Price of things"), string("dollars"), 10.10); 
    floatmom->print(string("Price of other things"), string("dollars"), 25.25); 
  }
  catch(std::exception &e){
    fail=true;
  }

  TEST_FAIL_AND_EXIT(*comm, !fail, "all print to standard output", 1);

  delete floatmom;

  // Node zero prints to a file.

  iPrint = (rank == 0);
  someOnePrints = true;
  comm->barrier();

  ios_base::openmode flags = ios_base::out & ios_base::trunc;

  ofstream outF("testMetricFile.txt", flags);

  try {
    doublemom = new MetricOutputManager<double>( rank, iPrint, outF, someOnePrints);
  }
  catch(std::exception &e){
    fail=true;
  }

  TEST_FAIL_AND_EXIT(*comm, !fail, "constructor", 1);

  if (rank==0){
    std::cout << "\nThere are " << nprocs << " processes. ";
    std::cout << "Rank zero only participates" << std::endl;
  }

  try{
    doublemom->print(string("Time to do something"), string("microseconds"), 
      10.101012345); 
    doublemom->print(string("Time to do something else"), string("microseconds"), 
      25.2500024); 
  }
  catch(std::exception &e){
    fail=true;
  }

  TEST_FAIL_AND_EXIT(*comm, !fail, "printing to a file", 1);

  outF.close();

  comm->barrier();

  if (rank == 0){
    ifstream inF("testMetricFile.txt");
    string s;
    while (getline(inF, s)){
      std::cout << s << std::endl;
    }
    inF.close();
    system("rm testMetricFile.txt");  // \todo fix for windows
  }

  comm->barrier();

  delete doublemom;

  if (rank==0)
   std::cout << "PASS" << std::endl;
}
