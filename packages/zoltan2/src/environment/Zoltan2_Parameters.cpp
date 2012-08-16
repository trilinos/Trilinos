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

/*! \file Zoltan2_Parameters.cpp
    \brief Methods that support the Zoltan2 ParameterList
*/

#include <Zoltan2_Parameters.hpp>
#include <Zoltan2_IntegerRangeList.hpp>
#include <Zoltan2_Standards.hpp>

#include <ostream>
#include <sstream>
#include <Teuchos_StandardParameterEntryValidators.hpp>

using namespace std;

namespace Zoltan2 {

/*! \brief  Create a list of all parameters.
 *
 *  \param pList  on return all Zoltan2 parameters with validators
 *                  have been added to pList
 *
 *  Parameters should be validating and documented.
 *
 *  New parameters should be added to createAllParameters.  Changes
 *  to existing parameters or their documentation should also be
 *  made in createAllParameters.
 *
 *  The comments in the top of createAllParameters describe the
 *  validators we use and what type of parameter each is for.  Find
 *  a type that matches the type of your parameter, and copy that
 *  definition of that parameter.
 *
 *  Include the comment above the parameterName definition.  Then you can
 *  grep on "topLevel" in this source file to see the parameter structure.
 *
 * To add a parameter to a parameter list "pl":
 *
 *    \code
 *    pl.set<data type of second argument>(
 *        nameOfParameter,
 *        anyValidValue,
 *        documentationString,
 *        validator)
 *    \endcode
 *
 * "anyValidValue" is not a default value.  Defaults are set in
 * library code if the parameter is not set by the user. The value
 * needs to be valid, but it is ignored.
 *
 * "documentationString" explains everything the user needs to no in
 * order to correctly use the parameter.  Validators can provide much
 * of the documentationString.
 *
 * The validators are described in the source code.  They are run when
 * the Environment::commitParameters() method is called.  They throw
 * an error with an explanation if the parameter is invalid.
 *
 *  \todo imbalance_tolerance should be a list of tolerances, one for
 *              each weight
 *  \todo Should all of our parameter values be strings, in order
 *          to support reading in parameters from a file?
 *  \todo Many parameters remain to be added here.
 */

void createAllParameters(Teuchos::ParameterList &pList)
{
  using Teuchos::tuple;
  using std::string;

  ////////////////////////////////////////////////////////////
  // Validators:

  // AnyNumber allows entry of a string, real or integer value.
  // The value can be retrieved with getInt, getDouble, getString.
  // NOTE: Errors seem to occur unless you retrieve the value as a double.
  // TODO: Make sure we are not using defaults because we don't find
  //    the int value.

  using Teuchos::AnyNumberParameterEntryValidator;
  RCP<const AnyNumberParameterEntryValidator> anyNumValidatorP;

  // EnhancedNumber specifies one data type and can specify
  // a valid minimum, maximum and step size.

  using Teuchos::EnhancedNumberValidator;
  RCP<const EnhancedNumberValidator<int> > intValidatorP;
  RCP<const EnhancedNumberValidator<size_t> > sizetValidatorP;
  RCP<const EnhancedNumberValidator<double> > doubleValidatorP;

  // StringValidator has string input and output. You can specify the 
  // list of valid strings.

  using Teuchos::StringValidator;
  RCP<const StringValidator> strValidatorP;

  // FileNameValidator expects a file name.  You can
  // specify whether the file must exist or not.

  using Teuchos::FileNameValidator;
  RCP<const FileNameValidator > fnameValidatorP;

  // StringToIntegral allows you to specify a list of strings,
  // a corresponding list of documentation strings, and a list
  // of integers that the strings map to.  The value entered is
  // a string. Retrieval methods are getIntegralValue, getStringValue,
  // and getStringDocs.
  //
  // We use it for on/off boolean values, for levels of verbosity,
  // and for output streams.

  using Teuchos::StringToIntegralParameterEntryValidator;
  typedef StringToIntegralParameterEntryValidator<int> str2intValidator;
  RCP<const str2intValidator> str2intValidatorP;

  // allowed values for boolean-type parameters

  Tuple<string,8> yesNoStrings = 
    tuple<string>( "true", "yes", "1", "on", "false", "no", "0", "off");

  Tuple<int,8> yesNoIntegrals = 
    tuple<int>( 1, 1, 1, 1, 0, 0, 0, 0);

  // allowed values for output streams

  Tuple<string,8> ostreamStrings = 
    tuple<string>( "std::cout", "cout", "stdout",
                   "std::cerr", "cerr", "stderr",
                   "/dev/null", "null");
                  
  Tuple<int,8> ostreamIntegrals = 
    tuple<int>( 0, 0, 0, 1, 1, 1, 2, 2);

  // IntegerRangeList is a Zoltan2 validator defined in this file's
  // header file.  An integer range list is entered as a string of 
  // comma separted tokens.  Tokens can be integers, integer ranges 
  // delimited by a dash ("-"), or the word "all".
  //
  // A specific data type is required, and valid minimum and maximum
  // integers may be specified.
  //
  // Examples:
  //    "1, 5, 10, 2, 12-15"
  //    "all"
  //    "0"
  //    "2-5"
  //    ""
  //
  // The integer range list is encoded as an Array<int>.  The last
  // integer in the array indicates whether the user wants "all", "none"
  // or whether the values are listed.  The list can be interpreted with
  // helper methods allValuesAreInRangeList, noValuesAreInRangeList,
  // IsInRangeList, and printIntegralRangeList.
  //
  // We use an IntegerRangeList for lists of MPI ranks that participate
  // in output of profiling information, and for lists of fixed vertices.
  //
  // If desired, strides can be added to integer range lists, but at this 
  // point in time we don't see a need for them.

  RCP<const IntegerRangeListValidator<int> > intRangeValidatorP;

  ///////////////////////////////////////////////////////////
  // Define every parameter.
  ///////////////////////////////////////////////////////////

  string parameterName;
  std::ostringstream docString;

  ///////////////////////////////////////////////////////////
  // LEVEL: top level, library configuration parameters
  ///////////////////////////////////////////////////////////

  ////////// topLevel/error_check_level
  parameterName = string("error_check_level");

  str2intValidatorP = rcp(new str2intValidator(

    tuple<string>("no_assertions",
                  "basic_assertions",
                  "complex_assertions",
                  "debug_mode_assertions"),

    tuple<string>(
      "no assertions will be performed",
      "typical checks of argument validity (fast, default)",
      "additional checks, i.e. is input graph a valid graph)",
      "check for everything including logic errors (slowest)"),

    tuple<int>(NO_ASSERTIONS,
               BASIC_ASSERTION, 
               COMPLEX_ASSERTION, 
               DEBUG_MODE_ASSERTION),

    parameterName));

  string omitInfo("the amount of error checking performed\n");
  omitInfo.append("(If the compile flag Z2_OMIT_ALL_ERROR_CHECKING was set,\n");
  omitInfo.append("then error checking code is not executed at runtime.)\n");
  docString.str("");
  str2intValidatorP->printDoc(omitInfo, docString); 
      
  pList.set<string>(parameterName, "basic_assertions", docString.str(),
    str2intValidatorP);

  ////////// topLevel/debug_level
  parameterName = string("debug_level");

  str2intValidatorP = rcp(new str2intValidator(
     tuple<string>("no_status",
                  "basic_status",
                  "detailed_status",
                  "verbose_detailed_status"),

     tuple<string>(
      "library outputs no status information",
      "library outputs basic status information (default)",
      "library outputs detailed information",
      "library outputs very detailed information"),

     tuple<int>(
        NO_STATUS,
        BASIC_STATUS,
        DETAILED_STATUS,
        VERBOSE_DETAILED_STATUS),
   
     parameterName));

  omitInfo = string("the amount of status/warning/debugging info printed\n");
  omitInfo.append("(If the compile flag Z2_OMIT_ALL_STATUS_MESSAGES was set,\n");
  omitInfo.append("then message output code is not executed at runtime.)\n");
  docString.str("");
  str2intValidatorP->printDoc(
    "the amount of status/debugging/warning information to print\n", docString);

  pList.set<string>(parameterName, "basic_status", docString.str(), 
    str2intValidatorP);

  ////////// topLevel/timer_type

  parameterName = string("timer_type");

  str2intValidatorP = rcp(new str2intValidator(
    tuple<string>(
     "no_timers", "macro_timers", "micro_timers", "both_timers", "test_timers"),

    tuple<string>(
      "No timing data will be collected (the default).",
      "Time an algorithm (or other entity) as a whole.",
      "Time the substeps of an entity.",
      "Run both MACRO and MICRO timers.",
      "Run timers added to code for testing, removed later"),

    tuple<int>(NO_TIMERS, MACRO_TIMERS, MICRO_TIMERS, BOTH_TIMERS, TEST_TIMERS),
   
     parameterName));

  omitInfo = string("the type of timing information to collect\n");
  omitInfo.append("(If the compile flag Z2_OMIT_ALL_PROFILING was set,\n");
  omitInfo.append("then the timing code is not executed at runtime.)\n");
  docString.str("");
  str2intValidatorP->printDoc(omitInfo, docString);

  pList.set<string>(parameterName, "no_timers", docString.str(), 
    str2intValidatorP);

  //////////
  // Debug, timing, and memory profiling output streams
  //////////

  ////////// topLevel/debug_output_stream
  parameterName = string("debug_output_stream");

  str2intValidatorP =
    rcp(new str2intValidator(ostreamStrings, ostreamIntegrals, parameterName));

  docString.str("");
  str2intValidatorP->printDoc(
    "output stream for debug/status/warning messages (default cout)\n",
    docString);

  pList.set<string>(parameterName, "cout", docString.str(),
    str2intValidatorP);

  ////////// topLevel/timer_output_stream
  parameterName = string("timer_output_stream");

  str2intValidatorP =
    rcp(new str2intValidator(ostreamStrings, ostreamIntegrals, parameterName));

  docString.str("");
  str2intValidatorP->printDoc(
    "output stream for timing report (default cout)\n",
    docString);

  pList.set<string>(parameterName, "cout", docString.str(),
    str2intValidatorP);

  ////////// topLevel/memory_output_stream
  parameterName = string("memory_output_stream");

  str2intValidatorP =
    rcp(new str2intValidator(ostreamStrings, ostreamIntegrals, parameterName));

  docString.str("");
  str2intValidatorP->printDoc(
    "output stream for memory usage messages (default cout)\n",
    docString);

  pList.set<string>(parameterName, "cout", docString.str(),
    str2intValidatorP);

  //////////
  // Debug, timing and memory profiling file name parameters
  //////////

  fnameValidatorP = rcp(new FileNameValidator(false));

  ////////// topLevel/debug_output_file
  parameterName = string("debug_output_file");

  docString.str("");
  fnameValidatorP->printDoc(
    "name of file to which debug/status messages should be written\n"
    "(process rank will be included in file name)\n",
     docString);

  pList.set<string>(parameterName, "/dev/null", docString.str(), fnameValidatorP);

  ////////// topLevel/timer_output_file
  parameterName = string("timer_output_file");

  docString.str("");
  fnameValidatorP->printDoc(
    "name of file to which timing information should be written\n"
    "(process rank will be included in file name)\n",
     docString);

  pList.set<string>(parameterName, "/dev/null", docString.str(), fnameValidatorP);

  ////////// topLevel/memory_output_file
  parameterName = string("memory_output_file");

  docString.str("");
  fnameValidatorP->printDoc(
    "name of file to which memory profiling information should be written\n"
    "(process rank will be included in file name)\n",
     docString);

  pList.set<string>(parameterName, "/dev/null", docString.str(), fnameValidatorP);

  //////////
  // Debug, timing and memory profiling: list of processes that print
  //////////


  RCP<const IntegerRangeListValidator<int> > intRangeValidatorUnsertedP = rcp(new IntegerRangeListValidator<int>(true));

  ////////// topLevel/debug_procs
  parameterName = string("debug_procs");

  docString.str("");
  intRangeValidatorUnsertedP->printDoc(
    "list of ranks that output debugging/status messages (default \"0\")\n",
     docString);

  pList.set<string>(parameterName, "0", docString.str(), intRangeValidatorUnsertedP);



  ////////// topLevel/debug_procs



  RCP<const IntegerRangeListValidator<partId_t> > pqRangeRangeValidatorP = rcp(new IntegerRangeListValidator<partId_t>(true));
  parameterName = string("pqParts");

  docString.str("");
  pqRangeRangeValidatorP->printDoc(
    "list of parts for pqJagged partitioning algorithm. As many as the dimension count.\n",
     docString);

  pList.set<string>(parameterName, "0", docString.str(), pqRangeRangeValidatorP);

  ////////// topLevel/memory_procs
  intRangeValidatorP = rcp(new IntegerRangeListValidator<int>(true));
  parameterName = string("memory_procs");

  docString.str("");
  intRangeValidatorP->printDoc(
    "list of ranks that memory profiling information (default \"0\")\n",
     docString);

  pList.set<string>(parameterName, "0", docString.str(), intRangeValidatorP);



  parameterName = string("parallel_part_calculation_count");

  anyNumValidatorP = rcp(new AnyNumberParameterEntryValidator);

  docString.str("");
  anyNumValidatorP->printDoc(
      "The number of parts whose cut coordinates will be calculated concurently.",
      docString);

  pList.set<string>(parameterName, "1", docString.str(),
    anyNumValidatorP);
  ///////////////////////////////////////////////////////////
  // LEVEL: top level, general problem parameters
  ///////////////////////////////////////////////////////////

  ////////// topLevel/speed_versus_quality
  parameterName = string("speed_versus_quality");

  strValidatorP = rcp(new StringValidator(
    tuple<string>("speed", "balance", "quality")));

  docString.str("");
  strValidatorP->printDoc(
    "When algorithm choices exist, opt for speed or solution quality?\n"
    "(Default is a balance of speed and quality)\n",
     docString);

  pList.set<string>(parameterName, "balance", docString.str(), strValidatorP);

  ////////// topLevel/memory_versus_speed
  parameterName = string("memory_versus_speed");

  strValidatorP = rcp(new StringValidator(
    tuple<string>("memory", "balance", "speed")));

  docString.str("");
  strValidatorP->printDoc(
    "When algorithm choices exist, opt for the use of less memory\n" 
    "at the expense of runtime\n" 
    "(Default is a balance of memory conservation and speed)\n",
     docString);

  pList.set<string>(parameterName, "balance", docString.str(), strValidatorP);

  ////////// topLevel/random_seed
  parameterName = string("random_seed");

  anyNumValidatorP = rcp(new AnyNumberParameterEntryValidator);

  docString.str("");
  anyNumValidatorP->printDoc("random seed\n", docString);

  pList.set<string>(parameterName, "0.5", docString.str(), anyNumValidatorP);

  ///////////////////////////////////////////////////////////
  // LEVEL: Sub list, ordering problem parameters
  ///////////////////////////////////////////////////////////

  ParameterList &ordering = pList.sublist("ordering", false,
    string("Ordering problem parameters"));

  ////////// topLevel/ordering/order_method
  parameterName = string("order_method");
  strValidatorP = rcp(new StringValidator(
    tuple<string>("rcm", "minimum_degree", "method3")));

  docString.str("");
  strValidatorP->printDoc(
    "Document the order_method parameter here\n"
    "(Default is ?)\n",
     docString);

  ordering.set<string>(parameterName, "rcm", docString.str(), strValidatorP);

  ////////// topLevel/ordering/order_package
  parameterName = string("order_package");
  strValidatorP = rcp(new StringValidator(
    tuple<string>("amd", "package2", "package3")));

  docString.str("");
  strValidatorP->printDoc(
    "Document the order_package parameter here\n"
    "(Default is ?)\n",
     docString);

  ordering.set<string>(parameterName, "amd", docString.str(), strValidatorP);


  ///////////////////////////////////////////////////////////
  // LEVEL: Sub list, partitioning problem parameters
  ///////////////////////////////////////////////////////////

  ParameterList &partitioning = pList.sublist("partitioning", false, 
    string("Partitioning problem parameters"));

 ////////// topLevel/partitioning/compute_metrics
  parameterName = string("compute_metrics");

  str2intValidatorP =
    rcp(new str2intValidator(yesNoStrings, yesNoIntegrals, parameterName));

  docString.str("");
  str2intValidatorP->printDoc("Compute metrics after computing solution\n",
    docString);

  partitioning.set<string>(parameterName, "no", docString.str(),
    str2intValidatorP);

  ////////// topLevel/partitioning/topology
  parameterName = string("topology");

  strValidatorP = rcp(new StringValidator);

  docString.str("");
  docString << "Topology of node to be used in hierarchical partitioning\n";
  docString << "  \"2,4\"  for dual-socket quad-core\n";
  docString << "  \"2,2,6\"  for dual-socket, dual-die, six-core\n";
  docString << "  \"2,2,3\"  for dual-socket, dual-die, six-core but\n";
  docString << "             with only three partitions per die\n";

  partitioning.set<string>(parameterName, "", docString.str(), strValidatorP);
  
  ////////// topLevel/partitioning/randomize_input
  parameterName = string("randomize_input");
    
  str2intValidatorP = 
    rcp(new str2intValidator(yesNoStrings, yesNoIntegrals, parameterName));

  docString.str("");
  str2intValidatorP->printDoc("randomize input prior to partitioning\n",
    docString);
  
  partitioning.set<string>(parameterName, "no", docString.str(),
    str2intValidatorP);

  ////////// topLevel/partitioning/objective
  parameterName = string("objective");

  strValidatorP = rcp(new StringValidator(
    tuple<string>( 
        "balance_object_count",
        "balance_object_weight",
        "multicriteria_minimize_total_weight",
        "multicriteria_minimize_maximum_weight",
        "multicriteria_balance_total_maximum",
        "minimize_cut_edge_count",
        "minimize_cut_edge_weight",
        "minimize_neighboring_parts",
        "minimize_boundary_vertices")));

  docString.str("");
  strValidatorP->printDoc(
    "objective of partitioning (default depends on algorithm)\n",
     docString);

  partitioning.set<string>(parameterName, "balance_object_weight", 
    docString.str(), strValidatorP);

  ////////// topLevel/partitioning/imbalance_tolerance
  parameterName = string("imbalance_tolerance");

  anyNumValidatorP = rcp(new AnyNumberParameterEntryValidator);

  docString.str("");
  anyNumValidatorP->printDoc(
   "imbalance tolerance, ratio of maximum load over average load"
   " (default 1.1)\n", 
      docString);

  partitioning.set<string>(parameterName, "1.1", docString.str(), 
    anyNumValidatorP);
  





  ////////// topLevel/partitioning/num_global_parts
  parameterName = string("num_global_parts");

  anyNumValidatorP = rcp(new AnyNumberParameterEntryValidator);

  docString.str("");
  anyNumValidatorP->printDoc(
   "global number of parts to compute (default is number of processes)\n",
      docString);

  partitioning.set<string>(parameterName, "0", docString.str(), 
    anyNumValidatorP);
  
  ////////// topLevel/partitioning/num_local_parts
  parameterName = string("num_local_parts");

  anyNumValidatorP = rcp(new AnyNumberParameterEntryValidator);

  docString.str("");
  anyNumValidatorP->printDoc(
   "number of parts to compute for this process(default is one)\n",
      docString);

  partitioning.set<string>(parameterName, "0", docString.str(), 
    anyNumValidatorP);

  ////////// topLevel/partitioning/approach
  parameterName = string("approach");

  strValidatorP = rcp(new StringValidator(
    tuple<string>("partition", "repartition", "maximize_overlap")));

  docString.str("");
  strValidatorP->printDoc(
    "Partition from scratch, partition incrementally from current\n"
    "partition, of partition from scratch but maximize overlap\n"
    "with the current partition (default is \"partition\" from scratch)\n",
     docString);

  partitioning.set<string>(parameterName, "partition", docString.str(), 
    strValidatorP);
  
  ////////// topLevel/partitioning/objects
  parameterName = string("objects");

  strValidatorP = rcp(new StringValidator(
    tuple<string>(
      "matrix_rows",
      "matrix_columns",
      "matrix_nonzeros",
      "mesh_elements",
      "mesh_nodes",
      "graph_edges",
      "graph_vertices",
      "coordinates",
      "identifiers")));

  docString.str("");
  strValidatorP->printDoc(
    "Objects to be partitioned (defaults are \"matrix_rows\" for\n"
    "matrix input, \"mesh_nodes\" for mesh input, and \"graph_vertices\"\n"
    "for graph input)\n",
     docString);

  partitioning.set<string>(parameterName, "graph_vertices", docString.str(), 
    strValidatorP);
  
  ////////// topLevel/partitioning/model
  parameterName = string("model");

  strValidatorP = rcp(new StringValidator(
    tuple<string>("hypergraph", "graph", "geometry", "ids")));

  docString.str("");
  strValidatorP->printDoc(
    "This is a low level parameter.  Normally the library will choose\n"
    "a computational model based on the algorithm or objective specified\n"
    "by the user.\n",
     docString);
  
  partitioning.set<string>(parameterName, "graph", docString.str(), 
    strValidatorP);
  
  ////////// topLevel/partitioning/algorithm
  parameterName = string("algorithm");

  strValidatorP = rcp(new StringValidator(
    tuple<string>(
      "rcb",
      "PQJagged",
      "rib",
      "hsfc",
      "patoh",
      "phg",
      "metis",
      "parmetis",
      "scotch",
      "ptscotch",
      "block",
      "cyclic",
      "random")));

  docString.str("");
  strValidatorP->printDoc("partitioning algorithm\n", docString);
  
  partitioning.set<string>(parameterName, "random", docString.str(), 
    strValidatorP);

  ///////////////////////////////////////////////////////////
  // LEVEL: Sub sub list, geometric partitioning parameters
  ///////////////////////////////////////////////////////////

  ParameterList &geom = partitioning.sublist("geometric", false, 
    string("geometric partitioning problem parameters"));

  ////////// topLevel/partitioning/geometric/rectilinear_blocks
  parameterName = string("rectilinear_blocks");

  str2intValidatorP =
    rcp(new str2intValidator(yesNoStrings, yesNoIntegrals, parameterName));

  docString.str("");
  str2intValidatorP->printDoc(
    "If true, then when a cut is made, all of the dots located on the cut\n"
    "are moved to the same side of the cut. The resulting regions are then\n"
    "rectilinear.  The resulting load balance may not be as good as when\n"
    "the group of dots is split by the cut. Default is false.\n",
    docString);

  geom.set<string>(parameterName, "no", docString.str(),
    str2intValidatorP);

  ////////// topLevel/partitioning/geometric/average_cuts
  parameterName = string("average_cuts");

  str2intValidatorP =
    rcp(new str2intValidator(yesNoStrings, yesNoIntegrals, parameterName));

  docString.str("");
  str2intValidatorP->printDoc(
    "When true, coordinates of RCB cutting planes are computed to be \n"
    "the average of the coordinates of the closest object on each side \n"
    "of the cut. Otherwise, coordinates of cutting planes may equal \n"
    "those of one of the closest objects. Default is false.\n",
    docString);

  geom.set<string>(parameterName, "no", docString.str(),
    str2intValidatorP);

  ////////// topLevel/partitioning/geometric/bisection_num_test_cuts
  parameterName = string("bisection_num_test_cuts");

  intValidatorP = rcp(new EnhancedNumberValidator<int>(1,250,1));

  docString.str("");
  intValidatorP->printDoc(
   "Experimental: number of test cuts to do simultaneously (default is 1)\n",
      docString);

  geom.set<int>(parameterName, 1, docString.str(), intValidatorP);

  ///////////////////////////////////////////////////////////
  // LEVEL: Sub sub list, hypergraph partitioning parameters
  ///////////////////////////////////////////////////////////

#if 0
  ParameterList &hg = partitioning.sublist("hypergraph", false, 
    string("hypergraph partitioning problem parameters"));
#endif

  ///////////////////////////////////////////////////////////
  // LEVEL: Sub sub list, graph partitioning parameters
  ///////////////////////////////////////////////////////////

  ParameterList &graph = partitioning.sublist("graph", false, 
    string("graph partitioning problem parameters"));
    
  ////////// topLevel/partitioning/graph/symmetrize_input
  parameterName = string("symmetrize_input");

  strValidatorP = rcp(new StringValidator(
    tuple<string>( "no", "transpose", "bipartite")));

  docString.str("");
  strValidatorP->printDoc(
    "Symmetrize input prior to partitioning.  If \"transpose\",\n"
    "symmetrize A by computing A plus ATranspose.  If \"bipartite\",\n"
    "A becomes [[0 A][ATranspose 0]].  \n",
    docString);

  graph.set<string>(parameterName, "no", docString.str(), 
    strValidatorP);

  ////////// topLevel/partitioning/graph/subset_graph
  parameterName = string("subset_graph");
    
  str2intValidatorP = 
    rcp(new str2intValidator(yesNoStrings, yesNoIntegrals, parameterName));

  docString.str("");
  str2intValidatorP->printDoc(
    "If \"yes\", the graph input is to be subsetted.  If a vertex neighbor\n"
    "is not a valid vertex, it will be omitted from the graph.  Otherwise,\n"
    "an invalid neighbor identifier is considered an error.\n",
    docString);
  
  graph.set<string>(parameterName, "no", docString.str(),
    str2intValidatorP);
  
  // TODO: Parmetis and Scotch parameters.

  ///////////////////////////////////////////////////////////
  // LEVEL: Sub list, ordering problem parameters  TODO
  ///////////////////////////////////////////////////////////

#if 0
  ParameterList &ordering = pList.sublist("ordering", false, 
    string("Ordering problem parameters"));
#endif

  ///////////////////////////////////////////////////////////
  // LEVEL: Sub list, coloring problem parameters  TODO
  ///////////////////////////////////////////////////////////

#if 0
  ParameterList &coloring = pList.sublist("coloring", false, 
    string("Coloring problem parameters"));
#endif

  ///////////////////////////////////////////////////////////
  // LEVEL: Sub list, matching problem parameters  TODO
  ///////////////////////////////////////////////////////////

#if 0
  ParameterList &matching = pList.sublist("matching", false, 
    string("Matching problem parameters"));
#endif
}

/*! \brief  Create a parameter list that can validate a 
 *          specific list of parameters.
 *
 *   \param plSome   the user's parameter list
 *   \param plAll    an empty parameter list that on return will contain
 *                    validating parameter entries for each
 *                    parameter in the user's parameter list.
 *   \param plVal     on return is a parameter list with all of the
 *                       user's parameters, but with validators.
 *
 *  Environment::commitParameters() calls 
 *  validateParametersAndSetDefaults() on the user's parameter list
 *  rather than validateParameters() because we want the validators' 
 *  validateAndModify() methods to be called rather than the validate()
 *  methods.  But unfortunately, validateAndModify() in addition to
 *  modifying the parameter if necessary also sets it to a default if
 *  the parameter does not appear in the user's parameter list.
 *
 *  We want the user's parameters to be modified, but we do not want unset 
 *  parameters to appear in the validated user's parameter list. To
 *  achieve this, we create a validating list that contains only the
 *  parameters that appear in the user's parameter list.
 *  
 */
static void setValidatorsInList(
  const Teuchos::ParameterList &plSome,   // in: user's parameters
  const Teuchos::ParameterList &plAll,    // in: validators for all params
  Teuchos::ParameterList &plVal)          // out: validators for user's params
{
  ParameterList::ConstIterator next = plSome.begin();

  while (next != plSome.end()){

    const std::string &name = next->first;
    const ParameterEntry &entrySome = plSome.getEntry(name);
    const ParameterEntry &entryAll = plAll.getEntry(name);

    if (entrySome.isList()){
      const ParameterList &sublistSome = plSome.sublist(name);   // get
      const ParameterList &sublistAll = plAll.sublist(name);     // get
      ParameterList &sublistVal = plVal.sublist(name);     // create & get
      setValidatorsInList(sublistSome, sublistAll, sublistVal);
    }
    else{
      plVal.setEntry(name, entryAll);
    }

    ++next;
  }
}

void createValidatorList(
   const Teuchos::ParameterList &plIn,
   Teuchos::ParameterList &plOut)
{
  ParameterList allParameters;

  createAllParameters(allParameters);

  setValidatorsInList(plIn, allParameters, plOut);
}

// Why isn't there a Teuchos method that does this?

void printListDocumentation(
  const Teuchos::ParameterList &pl,
  std::ostream &os,
  std::string listNames)
{
  using std::string;

  if (listNames.size() == 0)
    listNames = string("top");

  Array<string> subLists;
  ParameterList::ConstIterator next = pl.begin();

  while (next != pl.end()){
    const string &name = next->first;
    const ParameterEntry &entry = pl.getEntry(name);

    if (entry.isList()){
      subLists.append(name);
    }
    else{
      string doc = entry.docString();
      os << "List: "<< listNames << ", parameter: " << name << "\n";
      if (doc.size())
        os << doc << "\n";
    }

    ++next;
  }

  for (int i=0; i < subLists.size(); i++){
    string newListName = listNames + string("/") + subLists[i];
    const ParameterList &sublist = pl.sublist(subLists[i]);
    printListDocumentation(sublist, os, newListName);
  }
}


}  //namespace Zoltan2
