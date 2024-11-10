// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
// Testing parameters.  Serial test.

#include <Zoltan2_config.h>
#include <Zoltan2_Environment.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <string>

using std::string;

//
// For all parameters:
//   parameter name, a valid value, an invalid value

// FileNameValidator - a number is invalid

#define NUMFN 3
static string fnParams[NUMFN][3]={
  {"debug_output_file", "temp.txt", "5"},
  {"timer_output_file", "timerInfo.txt", "10.3"},
  {"memory_output_file", "memory.txt", "3.33"}
};

// Value is a particular string
#ifdef HAVE_ZOLTAN2_PULP
#define NUMSTR 38
#else
#define NUMSTR 35
#endif

static string strParams[NUMSTR][3]={
  {"error_check_level", "basic_assertions", "invalid_assertion_request"},
  {"debug_level", "basic_status", "invalid_status"},
  {"timer_type", "no_timers", "invalid_timers"},
  {"debug_output_stream", "std::cout", "invalid_stream"},
  {"timer_output_stream", "/dev/null", "invalid_stream"},
  {"memory_output_stream", "std::cerr", "invalid_stream"},
  {"debug_procs", "all", "not_a_valid_list_of_any_type"},
  {"mj_parts", "2,3,4", "not_a_valid_list_of_any_type"},
  {"memory_procs", "2-10", "not_a_valid_list_of_any_type"},
  {"order_method", "rcm", "invalid_method"},
  {"order_method_type", "local", "invalid_method_type"},
  {"order_package", "amd", "invalid_package"},
  {"partitioning_objective", "balance_object_weight", "invalid_objective"},
  {"partitioning_approach", "repartition", "invalid_approach"},
  {"objects_to_partition", "graph_vertices", "invalid_objects"},
  {"model", "graph", "invalid_model"},
  {"algorithm", "rcb", "invalid_algorithm"},
  {"symmetrize_input", "transpose", "invalid_option"},
  {"symmetrize_input", "transpose", "invalid_option"},
  {"mj_concurrent_part_count", "0", "invalid_value"},          // AnyNumberParameterEntryValidator
  {"mj_recursion_depth", "0", "invalid_value"},                // AnyNumberParameterEntryValidator
  {"mapping_type", "0", "invalid_value"},                      // AnyNumberParameterEntryValidator
  {"imbalance_tolerance", "1.1", "invalid_option"},            // AnyNumberParameterEntryValidator
  {"mj_minimum_migration_imbalance", "1.1", "invalid_option"}, // AnyNumberParameterEntryValidator
#ifdef HAVE_ZOLTAN2_PULP
  {"pulp_vert_imbalance", "1.1", "invalid_option"},            // AnyNumberParameterEntryValidator
  {"pulp_edge_imbalance", "1.1", "invalid_option"},            // AnyNumberParameterEntryValidator
  {"pulp_imbalance", "1.1", "invalid_option"},            // AnyNumberParameterEntryValidator
#endif // HAVE_ZOLTAN2_PULP
  {"scotch_imbalance_ratio", "1.1", "invalid_option"},         // AnyNumberParameterEntryValidator
  {"compute_metrics", "false", "invalid_bool_setting"},        // BoolParameterEntryValidator - accepts true/false/"true"/"false"
  {"rectilinear", "false", "invalid_bool_setting"},            // BoolParameterEntryValidator - accepts true/false/"true"/"false"
  {"subset_graph", "false", "invalid_bool_setting"},           // BoolParameterEntryValidator - accepts true/false/"true"/"false"
  {"mj_enable_rcb", "true", "invalid_bool_setting"},           // BoolParameterEntryValidator - accepts true/false/"true"/"false"
  {"mj_keep_part_boxes", "true", "invalid_bool_setting"},      // BoolParameterEntryValidator - accepts true/false/"true"/"false"
  {"num_global_parts", "1", "invalid_value"},                  // EnhancedNumberValidator
  {"num_local_parts", "0", "invalid_value"},                   // EnhancedNumberValidator
  {"mj_premigration_option", "1", "invalid_value"},               // EnhancedNumberValidator
  {"mj_migration_option", "2", "invalid_value"},               // EnhancedNumberValidator
  {"mj_num_teams", "60", "invalid_value"},                      // EnhancedNumberValidator
};

template <typename T>
int testInvalidValue( Teuchos::ParameterList &pl, 
  string paramName, T badValue)
{
  Teuchos::ParameterList validParameters;
  pl.set(paramName, badValue);
  std::cout << std::endl;
  std::cout << paramName << " = " << badValue << std::endl;

  bool failed = false;
  try{
    Zoltan2::createValidatorList(pl, validParameters);
    pl.validateParametersAndSetDefaults(validParameters);
  }
  catch(std::exception &e){
    std::cout << "Correctly generated an error:" << std::endl;
    std::cout << e.what() << std::endl;
    failed = true;
  }

  if (!failed){
    std::cerr << "Bad parameter value was not detected in parameter list." << std::endl;
    return 1;
  }
  return 0;
}

// this we can remove later
// kept here temporarily for reference
int testForIssue612()
{
  // Testing AnyNumberParameterEntryValidator

  // Create a parameter list
  Teuchos::ParameterList valid("valid parameter list");

  // Create parameter using validator
  typedef Teuchos::AnyNumberParameterEntryValidator validator_t;
  Teuchos::RCP<const validator_t> anyNumVal = Teuchos::rcp(new validator_t);

  ///////////////////////////////////////////////////////////////////
  // Initial test:  set parameter to 0.5 in valid parameter list.
  // Need to use the *expected* data type here.
  std::cout << "set good default value" << std::endl;

  std::string parameterName("parameterName");
  try {
    valid.set(parameterName, 5.0, "parameterDoc", anyNumVal);
  }
  catch (std::exception &e) {
    std::cout << "FAIL  error setting good default value "
              << e.what() << std::endl;
    return -1;
  }

  double dd = valid.getEntry(parameterName).getValue<double>(&dd);
  std::cout << "good default value <double> = " << dd << std::endl;

  // User creates his own parameter list and passes things of various types.
  // The user list must be validated against the valid list.
  Teuchos::ParameterList user("user");

  // This one should work
  std::cout << "test good user value" << std::endl;
  user.set(parameterName, "0.123");
  try {
    user.validateParametersAndSetDefaults(valid);
  }
  catch (std::exception &e) {
    std::cout << "FAIL  " << e.what() << std::endl;
    return -1;
  }

  dd = user.getEntry(parameterName).getValue<double>(&dd);
  std::cout << "good user value <double> = " << dd << std::endl;

  // This one should not work; the user's string is not a number
  std::cout << "test bogus user value" << std::endl;
  bool aok = false;

  // MDM note - the fail point will now be on user.set
  // std::stod will throw on this since we have added the validator
  try {
    user.set(parameterName, "bogus");
  }
  catch(std::exception &e) {
    // correct behavior
    std::cout << "Parameter list correctly rejected bogus user value."
              << std::endl;
    aok = true;
  }

  if (!aok) {
    std::cout << "FAIL  parameter list accepted a bogus user value"
              << std::endl;
    return -1;
  }

  // Test the valid with a bogus default value.  This operation should also
  // not work.  The validator should catch the bogus input.
  std::cout << "set bogus default value" << std::endl;

  std::string parameterNameToo("parameterNameToo");
  aok = false;
  try {
    valid.set(parameterNameToo, "bogus", "parameterDoc", anyNumVal);
  }
  catch (std::exception &e) {
    // correct behavior
    std::cout << "Parameter list correctly rejected bogus default value."
              << std::endl;
    aok = true;
  }

  if (!aok) {
    std::cout << "FAIL  parameter list accepted a bogus default value"
              << std::endl;
    return -1;
  }

  std::cout << "PASS" << std::endl;
  return 0;
}

  // Print out all the documentation

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();

  if (rank > 0)
    return 0;

  // short term reference - to delete later
  // this was an example proposed for issue #612
  // just keeping it here as a reference point to be deleted in the future
  int tempTest = testForIssue612();
  if( tempTest != 0 ) {
    return tempTest;
  }

  // Create a valid parameter list.

  Teuchos::ParameterList validParameters;
  Teuchos::ParameterList myParams("testParameterList");

  for (int i=0; i < NUMSTR; i++){
    myParams.set(strParams[i][0], strParams[i][1]);
  }

  for (int i=0; i < NUMFN; i++){
    myParams.set(fnParams[i][0], fnParams[i][1]);
  }

  Teuchos::ParameterList origParams(myParams);

  // Normally an application would not call this.  The
  // Environment object will validate the entered parameters.

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

  std::cout << std::endl;
  std::cout << "Parameters after validation: " << std::endl;
  std::cout << myParams << std::endl;

  // Try invalid parameter values
  for (int i=0; i < NUMSTR; i++){
    Teuchos::ParameterList badParams(origParams);
    int fail = 
      testInvalidValue<string>(badParams, strParams[i][0], strParams[i][2]);
    if (fail){
      std::cout << "FAIL" << std::endl;
      return 1;
    }
  }

  for (int i=0; i < NUMFN; i++){
    Teuchos::ParameterList badParams(origParams);
    std::istringstream iss(fnParams[i][2]);
    double badVal;
    iss >> badVal;
    int fail = 
       testInvalidValue<double>(badParams, fnParams[i][0], badVal);
    if (fail){
      std::cout << "FAIL" << std::endl;
      return 1;
    }
  }


  // Print out all the documentation

  std::cout << std::endl;
  std::cout << "Parameter documentation:" << std::endl;
  Zoltan2::printListDocumentation(validParameters, std::cout, std::string()); 

  std::cout << "PASS"  << std::endl;
  return 0;
}
