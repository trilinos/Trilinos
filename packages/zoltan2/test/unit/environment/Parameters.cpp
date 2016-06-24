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
// Testing integer range list parameters.  Serial test.

#include <Zoltan2_config.h>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_IntegerRangeList.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterEntryValidator.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_VerySimpleVectorAdapter.hpp>
#include <Zoltan2_MappingProblem.hpp>
#include <Zoltan2_BasicIdentifierAdapter.hpp>

typedef Teuchos::Array<int> rangeList_t;

Teuchos::ParameterList runTest(Teuchos::ParameterList inputParams)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // a test problem pulled out of block.cpp - later we will need to make all types and test them here
  typedef Tpetra::Map<> Map_t;
  typedef Map_t::local_ordinal_type localId_t;
  typedef Map_t::global_ordinal_type globalId_t;
  typedef double scalar_t;
  size_t localCount = 40;
  globalId_t *globalIds = new globalId_t [localCount];
  globalId_t offset = 0;
  for (size_t i=0; i < localCount; i++)
    globalIds[i] = offset++;
  typedef Zoltan2::BasicUserTypes<scalar_t, localId_t, globalId_t> myTypes;
  typedef Zoltan2::BasicIdentifierAdapter<myTypes> inputAdapter_t;

  std::vector<const scalar_t *> noWeights;
  std::vector<int> noStrides;

  inputAdapter_t ia(localCount, globalIds, noWeights, noStrides);
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &inputParams);
  return problem.getEnvironment()->getParameters();
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();

  if (rank > 0)
    return 0;

  // Set a few parameters, and then validate them.
  Teuchos::ParameterList myParams("testParameterList");
  myParams.set("debug_level", "detailed_status");        
  myParams.set("debug_procs", "all");   
  myParams.set("debug_output_stream", "std::cout");
  myParams.set("timer_output_file", "appPerformance.txt");

  // Normally an application would not call this.  The
  // Environment object will validate the entered parameters.
  // Since debug_procs is an IntegerRangeList,
  // this call will convert it to a Teuchos::Array that uses
  // a special flag to indicate "all" or "none".
  Teuchos::ParameterList outputParams;
  try{
    outputParams = runTest(myParams);
  }
  catch(std::exception &e){
    std::cerr << "Validate parameters generated an error:" << std::endl;
    std::cerr << e.what() << std::endl;
    std::cerr << "FAIL" << std::endl;
    return 1;
  }

  std::cout << std::endl;
  std::cout << "A few parameters after validation: " << std::endl;
  std::cout << outputParams << std::endl;

  // MDM - Note that this was a reference and it will crash for me when a1 is set again later
  // Didn't work out the pipeline yet but switching to a non reference seems to resolve it
  rangeList_t a1 = outputParams.get<rangeList_t>("debug_procs");
  std::cout << "debug_procs translation: ";
  Zoltan2::printIntegralRangeList(std::cout, *a1);
  std::cout << std::endl;

  // Now let's enter a bad value for a parameter and make sure
  // we get an error.

  Teuchos::ParameterList faultyParams("badParameterList");
  faultyParams.set("debug_procs", "not-even-remotely-an-integer-range");
  bool failed = false;
  try{
    outputParams = runTest(faultyParams);
  }
  catch(std::exception &e){
    std::cout << std::endl;
    std::cout << "Invalid parameter correctly generated an error:" << std::endl;
    std::cout << e.what() << std::endl;
    failed = true;
  }

  if (!failed){
    std::cerr << "Bad parameter was not detected in parameter list." << std::endl;
    return 1;
  }

  // Now set every parameter to a reasonable value

  Teuchos::ParameterList all("setAllParametersList");


  all.set("debug_level", "basic_status");

  all.set("debug_procs", "1,2,5-10,2");
  all.set("memory_procs", "1,2,3,4,all");

  all.set("debug_output_stream", "std::cerr");
  all.set("timer_output_stream", "std::cout");
  all.set("memory_output_stream", "/dev/null");


  all.set("debug_output_file", "/home/me/debug.txt");
  all.set("timer_output_file", "/home/me/performance.txt");
  all.set("memory_output_file", "/home/me/memoryUsed.txt");

  all.set("speed_versus_quality", "speed");
  all.set("memory_versus_speed", "memory");

  all.set("error_check_level", "basic_assertions");

  all.set("random_seed", .12121212);

  all.set("topology", "2,6,6");

  all.set("randomize_input", "true");

  all.set("partitioning_objective", "minimize_cut_edge_weight");

  all.set("imbalance_tolerance", 1.2);

  all.set("num_global_parts", 12);
  all.set("num_local_parts", 2);

  all.set("partitioning_approach", "partition");

  all.set("objects_to_partition", "graph_vertices");

  all.set("model", "hypergraph");

  all.set("algorithm", "phg");

  all.set("symmetrize_input", "no");
  all.set("subset_graph", "false");

  try {
    outputParams = runTest(all);
  }
  catch(std::exception &e){
    std::cerr << "Validate parameters generated an error:" << std::endl;
    std::cerr << e.what() << std::endl;
    std::cerr << "FAIL" << std::endl;
    return 1;
  }

  std::cout << std::endl;
  std::cout << "All parameters validated and modified: ";
  std::cout << outputParams << std::endl;

  a1 = all.getPtr<rangeList_t>("debug_procs");
  std::cout << "debug_procs translation: ";
  Zoltan2::printIntegralRangeList(std::cout, *a1);
  std::cout << std::endl;

  a1 = all.getPtr<rangeList_t>("memory_procs");
  std::cout << "memory_procs translation: ";
  Zoltan2::printIntegralRangeList(std::cout, *a1);
  std::cout << std::endl;

  // Print out all the documentation

  std::cout << std::endl;
  std::cout << "Parameter documentation:" << std::endl;

  // MDM - will refactor this - should probably pull the validParameters from the Problem environment
  // Perhaps we want to do this for each test we run - eventually this needs to cover all the Problem types
  // Zoltan2::printListDocumentation(validParameters, std::cout, std::string());
  std::cout << "PASS"  << std::endl;
  return 0;
}
