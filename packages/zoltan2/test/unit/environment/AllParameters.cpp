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
// Testing parameters.  Serial test.

#include <Zoltan2_config.h>
#include <Zoltan2_Environment.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_DefaultComm.hpp>

// For all parameters:
//   parameter name, a valid value, an invalid value

#define NUMPARAMS 36

static char *(paramList[NUMPARAMS][3])={
  {"error_check_level", "basic_assertions", "invalid_assertion_request"},
  {"debug_level", "basic_status", "bogus_status"},
  {"timer_type", "no_timers", "invalid_timers"},
  {"debug_output_stream", "cout", "invalid_stream"},
  {"timer_output_stream", "/dev/null", "invalid_stream"},
  {"memory_output_stream", "cerr", "invalid_stream"},
  {"debug_output_file", "temp.txt", "----"},
  {"timer_output_file", "temp.txt", "----"},
  {"memory_output_file", "temp.txt", "----"},
  {"debug_procs", "all", "not_a_valid_list_of_any_type"},
  {"pqParts", "2,3,4", "not_a_valid_list_of_any_type"},
  {"memory_procs", "2-10", "not_a_valid_list_of_any_type"},
  {"speed_versus_quality", "balance", "amazing_performance"},
  {"memory_versus_speed", "memory", "impossible_performance"},
  {"random_seed", "9.999", "nan"},
  {"order_method", "rcm", "rcmNew"},
  {"order_package", "amd", "amdNew"},
  {"compute_metrics", "yes", "maybe"},
  {"topology", "2,3,6", "I_don't_know"},
  {"randomize_input", "1", "22"},
  {"partitioning_objective", "balance_object_weight", "get_curry"},
  {"imbalance_tolerance", "1.1", "intolerant"},
  {"num_global_parts", "12", "xxx"},
  {"num_local_parts", "1", "no_idea"},
  {"partitioning_approach", "repartition", "cut_it_up"},
  {"objects_to_partition", "vertices", "nothing"},
  {"model", "graph", "manifold"},
  {"algorithm", "rcb", "xyz"},
  {"rectilinear_blocks", "yes", "dont_know"},
  {"average_cuts", "false", "dont_know"},
  {"bisection_num_test_cuts", "3", "dont_know"},
  {"symmetrize_input", "true", "dont_know"},
  {"subset_graph", "1", "dont_know"},
  {"force_binary_search", "true", "dont_know"},
  {"force_linear_search", "yes", "dont_know"},
  {"parallel_part_calculation_count", "2", "dont_know"}
};

typedef Teuchos::Array<int> rangeList_t;

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();

  if (rank > 0)
    return 0;

  // Create a valid parameter list.

  Teuchos::ParameterList validParameters;
  Teuchos::ParameterList myParams("testParameterList");

  for (int i=0; i < NUMPARAMS; i++)
    myParams.set(paramList[i][0], paramList[i][1]);

  Teuchos::ParameterList origParams(myParams);

  // Normally an application would not call this.  The
  // Environment object will validate the entered parameters.

  try{
    Zoltan2::createValidatorList(myParams, validParameters);
    myParams.validateParametersAndSetDefaults(validParameters);
    Zoltan2::Environment::convertStringToInt(myParams);
  }
  catch(std::exception &e){
    std::cerr << "Validate parameters generated an error:" << endl;
    std::cerr << e.what() << endl;
    std::cerr << "FAIL" << endl;
    return 1;
  }

  cout << endl;
  cout << "Parameters after validation: " << endl;
  cout << myParams << endl;

  // Try invalid parameter values

  for (int i=0; i < NUMPARAMS; i++){
    validParameters = Teuchos::ParameterList();
    Teuchos::ParameterList badParams(origParams);
    badParams.set(paramList[i][0], paramList[i][2]);
    cout << endl;
    cout << paramList[i][0] << " = " << paramList[i][2] << endl;

    bool failed = false;
    try{
      Zoltan2::createValidatorList(badParams, validParameters);
      badParams.validateParametersAndSetDefaults(validParameters);
    }
    catch(std::exception &e){
      cout << "Correctly generated an error:" << endl;
      cout << e.what() << endl;
      failed = true;
    }

    if (!failed){
      cerr << "Bad parameter was not detected in parameter list." << endl;
      return 1;
    }
  }

  // Print out all the documentation

  cout << endl;
  cout << "Parameter documentation:" << endl;
  Zoltan2::printListDocumentation(validParameters, cout, std::string()); 

  cout << "PASS"  << endl;
  return 0;
}
