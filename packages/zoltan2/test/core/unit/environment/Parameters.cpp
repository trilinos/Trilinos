// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

typedef Teuchos::Array<int> rangeList_t;

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();

  if (rank > 0)
    return 0;

  // Set a few parameters, and then validate them.

  Teuchos::ParameterList validParameters;

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

  try{
    Zoltan2::createValidatorList(myParams, validParameters);
    myParams.validateParametersAndSetDefaults(validParameters);
  }
  catch(std::exception &e){
    std::cerr << "Validate parameters generated an error:" << std::endl;
    std::cerr << e.what() << std::endl;
    std::cerr << "FAIL" << std::endl;
    return 1;
  }

  validParameters = Teuchos::ParameterList();

  std::cout << std::endl;
  std::cout << "A few parameters after validation: " << std::endl;
  std::cout << myParams << std::endl;

  rangeList_t *a1 = myParams.getPtr<rangeList_t>("debug_procs");
  std::cout << "debug_procs translation: ";
  Zoltan2::printIntegralRangeList(std::cout, *a1);
  std::cout << std::endl;

  // Now let's enter a bad value for a parameter and make sure
  // we get an error.

  Teuchos::ParameterList faultyParams("badParameterList");
  faultyParams.set("debug_procs", "not-even-remotely-an-integer-range");
  bool failed = false;
  try{
    Zoltan2::createValidatorList(faultyParams, validParameters);
    faultyParams.validateParametersAndSetDefaults(validParameters);
  }
  catch(std::exception &e){
    std::cout << std::endl;
    std::cout << "Invalid parameter correctly generated an error:" << std::endl;
    std::cout << e.what() << std::endl;
    failed = true;
  }

  validParameters = Teuchos::ParameterList();

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

  all.set("error_check_level", "basic_assertions");

  all.set("partitioning_objective", "minimize_cut_edge_weight");

  all.set("imbalance_tolerance", 1.2);

  all.set("num_global_parts", 12);
  all.set("num_local_parts", 2);

  all.set("partitioning_approach", "partition");

  all.set("objects_to_partition", "graph_vertices");

  all.set("model", "hypergraph");

  all.set("algorithm", "phg");

  all.set("symmetrize_input", "no");
  all.set("subset_graph", false); // bool parameter

  try{
    Zoltan2::createValidatorList(all, validParameters);
    all.validateParametersAndSetDefaults(validParameters);
  }
  catch(std::exception &e){
    std::cerr << "Validate parameters generated an error:" << std::endl;
    std::cerr << e.what() << std::endl;
    std::cerr << "FAIL" << std::endl;
    return 1;
  }

  std::cout << std::endl;
  std::cout << "All parameters validated and modified: ";
  std::cout << all << std::endl;

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
  Zoltan2::printListDocumentation(validParameters, std::cout, std::string()); 

  std::cout << "PASS"  << std::endl;
  return 0;
}
