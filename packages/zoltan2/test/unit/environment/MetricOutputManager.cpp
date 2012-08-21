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
#include <Zoltan2_TestHelpers.hpp>

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
      rank, iPrint, std::cout, someOnePrints, string("units"), 10);
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
    intmom->print(string("number of things"), 10); 
    intmom->print(string("number of other things"), 5); 
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

      rank, iPrint, std::cout, someOnePrints, string("dollars"), 10);
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
    floatmom->print(string("Price of things"), 10.10); 
    floatmom->print(string("Price of other things"), 25.25); 
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
    doublemom = new MetricOutputManager<double>( rank, iPrint, outF, someOnePrints, string("microseconds"), 10);
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
    doublemom->print(string("Time to do something"), 
      10.101012345); 
    doublemom->print(string("Time to do something else"), 
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
