// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
// ***********************************************************************
//
// Testing initialization of parameters.  Serial test.

#include <Zoltan2_config.h>
#include <Zoltan2_Parameters.hpp>
#include <Zoltan2_Environment.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterEntryValidator.hpp>

typedef Teuchos::Array<int> rangeList_t;

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();

  if (rank > 0)
    return 0;

  // Set a few parameters, and then validate them.

  Teuchos::ParameterList validParameters;

  Teuchos::ParameterList myParams("testParameterList");

  myParams.set("debug_level", "detailed_status");        
  myParams.set("debug_procs", "all");   
  myParams.set("debug_output_stream", "std::cout");

  myParams.set("timing_level", "basic_status");  
  myParams.set("timing_procs", "0"); 
  myParams.set("timing_output_file", "appPerformance.txt");

  // Normally an application would not call this.  The
  // Environment object will validate the entered parameters.
  // Since debug_procs and timing_procs are IntegerRangeLists,
  // this call will convert them to Teuchos::Arrays that use
  // a special flag to indicate "all" or "none".

  try{
    Zoltan2::createValidatorList(myParams, validParameters);
    myParams.validateParametersAndSetDefaults(validParameters);
    Zoltan2::Environment::convertStringToInt(myParams);
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

  rangeList_t &a1 = myParams.get<rangeList_t>("debug_procs");
  std::cout << "debug_procs translation: ";
  Zoltan2::printIntegralRangeList(std::cout, a1);
  std::cout << std::endl;

  a1 = myParams.get<rangeList_t>("timing_procs");
  std::cout << "timing_procs translation: ";
  Zoltan2::printIntegralRangeList(std::cout, a1); 
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
  all.set("timing_level", "no_status");
  all.set("memory_profiling_level", "detailed_status");

  all.set("debug_procs", "1,2,5-10,2");
  all.set("timing_procs", "10-2");
  all.set("memory_profiling_procs", "1,2,3,4,all");

  all.set("debug_output_stream", "std::cerr");
  all.set("timing_output_stream", "std::cout");
  all.set("memory_profiling_output_stream", "/dev/null");


  all.set("debug_output_file", "/home/me/debug.txt");
  all.set("timing_output_file", "/home/me/performance.txt");
  all.set("memory_profiling_output_file", "/home/me/memoryUsed.txt");

  all.set("speed_versus_quality", "speed");
  all.set("memory_versus_speed", "memory");

  all.set("error_check_level", "basic_assertions");

  all.set("random_seed", .12121212);

  Teuchos::ParameterList &parParams = all.sublist("partitioning");

  parParams.set("topology", "2,6,6");

  parParams.set("randomize_input", "true");

  parParams.set("objective", "minimize_cut_edge_weight");

  parParams.set("imbalance_tolerance", 1.2);

  parParams.set("num_global_parts", 12);
  parParams.set("num_local_parts", 2);

  parParams.set("approach", "partition");

  parParams.set("objects", "graph_vertices");

  parParams.set("model", "hypergraph");

  parParams.set("algorithm", "phg");

  Teuchos::ParameterList &graphParams = parParams.sublist("graph");

  graphParams.set("symmetrize_input", "no");
  graphParams.set("subset_graph", "false");

  try{
    Zoltan2::createValidatorList(all, validParameters);
    all.validateParametersAndSetDefaults(validParameters);
    Zoltan2::Environment::convertStringToInt(all);
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

  a1 = all.get<rangeList_t>("debug_procs");
  std::cout << "debug_procs translation: ";
  Zoltan2::printIntegralRangeList(std::cout, a1);
  std::cout << std::endl;

  a1 = all.get<rangeList_t>("timing_procs");
  std::cout << "timing_procs translation: ";
  Zoltan2::printIntegralRangeList(std::cout, a1);
  std::cout << std::endl;

  a1 = all.get<rangeList_t>("memory_profiling_procs");
  std::cout << "memory_profiling_procs translation: ";
  Zoltan2::printIntegralRangeList(std::cout, a1);
  std::cout << std::endl;

  // Print out all the documentation

  std::cout << std::endl;
  std::cout << "Parameter documentation:" << std::endl;
  Zoltan2::printListDocumentation(validParameters, std::cout, std::string()); 

  std::cout << "PASS"  << std::endl;
  return 0;
}
