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
#include <Teuchos_ParameterList.hpp>

int main(int argc, char *argv[])
{
  Teuchos::ParameterList validParameters;

  Zoltan2::createValidParameterList(validParameters);
  std::cout << "Default, validating parameter list created: ";
  std::cout << validParameters.name() << std::endl;
  std::cout << validParameters << std::endl;

  // Set a few parameters, and then validate them.

  Teuchos::ParameterList myParams("testParameterList");

  myParams.set("debug_level", 2);        
  myParams.set("debug_procs", "all");   
  myParams.set("debug_output_stream", "std::cout");

  myParams.set("profiling_level", 1);  
  myParams.set("profiling_procs", "0"); 
  myParams.set("profiling_output_file", "appPerformance.txt");

  // Normally an application would not call this.  The
  // Environment object will validate the entered parameters.
  // Since debug_procs and profiling_procs are IntegerRangeLists,
  // this call will convert them to Teuchos::Arrays that use
  // a special flag to indicate "all" or "none".

  try{
    myParams.validateParametersAndSetDefaults(validParameters);
  }
  catch(std::exception &e){
    std::cerr << "Valid parameters generated an error:" << std::endl;
    std::cerr << e.what() << std::endl;
    std::cerr << "FAIL" << std::endl;
    return 1;
  }

  // Now let's see the validated parameters that Zoltan2 will see.

  int &debugLevel = myParams.get<int>(std::string("debug_level"));
  int &profilingLevel = myParams.get<int>(std::string("profiling_level"));

  Teuchos::Array<int> &debugProcs = 
    myParams.get<Teuchos::Array<int> >(std::string("debug_procs"));

  Teuchos::Array<int> &profilerProcs = 
    myParams.get<Teuchos::Array<int> >(std::string("profiling_procs"));

  std::string &debugOstream = myParams.get<std::string>(std::string("debug_output_stream"));
  std::string &debugFile = myParams.get<std::string>(std::string("debug_output_file"));

  std::string &profilingOstream = myParams.get<std::string>(std::string("profiling_output_stream"));
  std::string &profilingFile = myParams.get<std::string>(std::string("profiling_output_file"));

  std::cout << "Debug level " << debugLevel;
  std::cout << ", debug procs ";
  Zoltan2::printIntegralRangeList(std::cout, debugProcs);
  std::cout << std::endl;
  if (debugFile.size() > 0)
    std::cout << "Write debug output to " << debugFile << std::endl;
  else
    std::cout << "Write debug output to " << debugOstream << std::endl;

  std::cout << "Profiling level " << profilingLevel;
  std::cout << ", profiling procs ";
  Zoltan2::printIntegralRangeList(std::cout, profilerProcs);
  std::cout << std::endl;
  if (profilingFile.size() > 0)
    std::cout << "Write profiling output to " << profilingFile << std::endl;
  else
    std::cout << "Write profiling output to " << profilingOstream << std::endl;

  // Now let's leave some parameters unset and verify they get
  // set to the default.

  Teuchos::ParameterList lazyParams("incompleteParameterList");
  lazyParams.set("debug_procs", "1-10");
  try{
    lazyParams.validateParametersAndSetDefaults(validParameters);
  }
  catch(std::exception &e){
    std::cerr << "Valid single parameter generated an error:" << std::endl;
    std::cerr << e.what() << std::endl;
    std::cerr << "FAIL" << std::endl;
    return 1;
  }

  debugLevel = lazyParams.get<int>(std::string("debug_level"));
  profilingLevel = lazyParams.get<int>(std::string("profiling_level"));
  debugProcs = 
    lazyParams.get<Teuchos::Array<int> >(std::string("debug_procs"));
  profilerProcs = 
    lazyParams.get<Teuchos::Array<int> >(std::string("profiling_procs"));
  debugOstream = lazyParams.get<std::string>(std::string("debug_output_stream"));
  debugFile = lazyParams.get<std::string>(std::string("debug_output_file"));
  profilingOstream = lazyParams.get<std::string>(std::string("profiling_output_stream"));
  profilingFile = lazyParams.get<std::string>(std::string("profiling_output_file"));

  std::cout << "Debug level " << debugLevel;
  std::cout << ", debug procs ";
  Zoltan2::printIntegralRangeList(std::cout, debugProcs);
  std::cout << std::endl;
  if (debugFile.size() > 0)
    std::cout << "Write debug output to " << debugFile << std::endl;
  else
    std::cout << "Write debug output to " << debugOstream << std::endl;

  std::cout << "Profiling level " << profilingLevel;
  std::cout << ", profiling procs ";
  Zoltan2::printIntegralRangeList(std::cout, profilerProcs);
  std::cout << std::endl;
  if (profilingFile.size() > 0)
    std::cout << "Write profiling output to " << profilingFile << std::endl;
  else
    std::cout << "Write profiling output to " << profilingOstream << std::endl;

  // Now let's entry a bad value for a parameter and make sure
  // we get an error.

  Teuchos::ParameterList faultyParams("badParameterList");
  faultyParams.set("debug_procs", "not-even-remotely-an-integer-range");
  bool failed = false;
  try{
    faultyParams.validateParametersAndSetDefaults(validParameters);
  }
  catch(std::exception &e){
    std::cout << "Invalid parameter correctly generated an error:" << std::endl;
    std::cout << e.what() << std::endl;
    failed = true;
  }

  if (!failed){
    std::cerr << "Bad parameter was not detected in parameter list." << std::endl;
    return 1;
  }

  // Let's look at the entire parameter list that was set for us.

  std::cout << lazyParams << std::endl;

  std::cout << "PASS"  << std::endl;
}
