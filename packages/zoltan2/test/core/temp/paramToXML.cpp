// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \brief utility that creates an initial XML file of parameters.
 *   The goal is that after initial creation, we will edit the
 *   XML file when adding new parameters.
 */

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Teuchos_XMLParameterListWriter.hpp>
#include <Teuchos_ParameterXMLFileReader.hpp>
#include <Teuchos_ValidatorXMLConverterDB.hpp>

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Parameters.hpp>
#include <Zoltan2_IntegerRangeList.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using std::string;
using std::ostringstream;
using std::endl;
using std::cout;
using std::cout;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Tuple;



void createAllParameters(Teuchos::ParameterList &pList);

int main(int argc, char *argv[])
{
  string xmlFile("initParams.xml");
  Teuchos::XMLParameterListWriter plw;

  // XML converter for Zoltan2::IntegerRangeListValidator must be added
  // once to the validator database.

  typedef Zoltan2::IntegerRangeListValidatorXMLConverter<int> irlConverter_t;
  RCP<irlConverter_t > converter = rcp(new irlConverter_t);
  typedef Zoltan2::IntegerRangeListValidator<int> irl_t;
  RCP<const irl_t> intRangeValidatorP = rcp(new irl_t);

  Teuchos::ValidatorXMLConverterDB::addConverter(
       intRangeValidatorP,    // can be a dummy of this type
       converter);

  Teuchos::ParameterList pl("zoltan2ValidatingParameters");
  createAllParameters(pl);
    
  // Write out to XML
  Teuchos::XMLObject obj = plw.toXML(pl);

  cout << "Parameter list: " << endl;
  cout << obj << endl;
  
  std::ofstream of;
  of.open(xmlFile.c_str());
  of << obj << endl;
  of.close();
}

void createAllParameters(Teuchos::ParameterList &pList)
{
  using Teuchos::tuple;
  using std::string;

  using Teuchos::AnyNumberParameterEntryValidator;
  RCP<const AnyNumberParameterEntryValidator> anyNumValidatorP;

  using Teuchos::EnhancedNumberValidator;
  RCP<const EnhancedNumberValidator<int> > intValidatorP;

  using Teuchos::StringValidator;
  RCP<const StringValidator> strValidatorP;

  using Teuchos::FileNameValidator;
  RCP<const FileNameValidator > fnameValidatorP;

  using Teuchos::StringToIntegralParameterEntryValidator;
  typedef StringToIntegralParameterEntryValidator<int> str2intValidator;
  RCP<const str2intValidator> str2intValidatorP;

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

  RCP<const Zoltan2::IntegerRangeListValidator<int> > intRangeValidatorP;

  string parameterName;
  std::ostringstream docString;

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

    tuple<int>(Zoltan2::NO_ASSERTIONS,
               Zoltan2::BASIC_ASSERTION, 
               Zoltan2::COMPLEX_ASSERTION, 
               Zoltan2::DEBUG_MODE_ASSERTION),

    parameterName));

  string omitInfo("the amount of error checking performed\n");
  omitInfo.append("(If the compile flag Z2_OMIT_ALL_ERROR_CHECKING was set,\n");
  omitInfo.append("then error checking code is not executed at runtime.)\n");
  docString.str("");
  str2intValidatorP->printDoc(omitInfo, docString); 
      
  pList.set<string>(parameterName, "basic_assertions", docString.str(),
    str2intValidatorP);

  ////////// debug_level
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
        Zoltan2::NO_STATUS,
        Zoltan2::BASIC_STATUS,
        Zoltan2::DETAILED_STATUS,
        Zoltan2::VERBOSE_DETAILED_STATUS),
   
     parameterName));

  omitInfo = string("the amount of status/warning/debugging info printed\n");
  omitInfo.append("(If the compile flag Z2_OMIT_ALL_STATUS_MESSAGES was set,\n");
  omitInfo.append("then message output code is not executed at runtime.)\n");
  docString.str("");
  str2intValidatorP->printDoc(
    "the amount of status/debugging/warning information to print\n", docString);

  pList.set<string>(parameterName, "basic_status", docString.str(), 
    str2intValidatorP);

  ////////// timer_type

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

    tuple<int>(Zoltan2::NO_TIMERS, Zoltan2::MACRO_TIMERS, Zoltan2::MICRO_TIMERS, Zoltan2::BOTH_TIMERS, Zoltan2::TEST_TIMERS),
   
     parameterName));

  omitInfo = string("the type of timing information to collect\n");
  omitInfo.append("(If the compile flag Z2_OMIT_ALL_PROFILING was set,\n");
  omitInfo.append("then the timing code is not executed at runtime.)\n");
  docString.str("");
  str2intValidatorP->printDoc(omitInfo, docString);

  pList.set<string>(parameterName, "no_timers", docString.str(), 
    str2intValidatorP);

  ////////// debug_output_stream
  parameterName = string("debug_output_stream");

  str2intValidatorP =
    rcp(new str2intValidator(ostreamStrings, ostreamIntegrals, parameterName));

  docString.str("");
  str2intValidatorP->printDoc(
    "output stream for debug/status/warning messages (default cout)\n",
    docString);

  pList.set<string>(parameterName, "cout", docString.str(),
    str2intValidatorP);

  ////////// timer_output_stream
  parameterName = string("timer_output_stream");

  str2intValidatorP =
    rcp(new str2intValidator(ostreamStrings, ostreamIntegrals, parameterName));

  docString.str("");
  str2intValidatorP->printDoc(
    "output stream for timing report (default cout)\n",
    docString);

  pList.set<string>(parameterName, "cout", docString.str(),
    str2intValidatorP);

  ////////// memory_output_stream
  parameterName = string("memory_output_stream");

  str2intValidatorP =
    rcp(new str2intValidator(ostreamStrings, ostreamIntegrals, parameterName));

  docString.str("");
  str2intValidatorP->printDoc(
    "output stream for memory usage messages (default cout)\n",
    docString);

  pList.set<string>(parameterName, "cout", docString.str(),
    str2intValidatorP);


  ////////// debug_output_file
  parameterName = string("debug_output_file");

  fnameValidatorP = rcp(new FileNameValidator(false));
  docString.str("");
  fnameValidatorP->printDoc(
    "name of file to which debug/status messages should be written\n"
    "(process rank will be included in file name)\n",
     docString);

  pList.set<string>(parameterName, "/dev/null", docString.str(), fnameValidatorP);

  ////////// timer_output_file
  parameterName = string("timer_output_file");

  docString.str("");
  fnameValidatorP->printDoc(
    "name of file to which timing information should be written\n"
    "(process rank will be included in file name)\n",
     docString);

  pList.set<string>(parameterName, "/dev/null", docString.str(), fnameValidatorP);

  ////////// memory_output_file
  parameterName = string("memory_output_file");

  docString.str("");
  fnameValidatorP->printDoc(
    "name of file to which memory profiling information should be written\n"
    "(process rank will be included in file name)\n",
     docString);

  pList.set<string>(parameterName, "/dev/null", docString.str(), fnameValidatorP);

  ////////// debug_procs
  parameterName = string("debug_procs");

  RCP<const Zoltan2::IntegerRangeListValidator<int> > intRangeValidatorUnsertedP = rcp(new Zoltan2::IntegerRangeListValidator<int>(true));

  docString.str("");
  intRangeValidatorUnsertedP->printDoc(
    "list of ranks that output debugging/status messages (default \"0\")\n",
     docString);

  pList.set<string>(parameterName, "0", docString.str(), intRangeValidatorUnsertedP);

  ////////// pqParts
  parameterName = string("pqParts");

  RCP<const Zoltan2::IntegerRangeListValidator<int> > pqRangeRangeValidatorP = rcp(new Zoltan2::IntegerRangeListValidator<int>(true));

  docString.str("");
  pqRangeRangeValidatorP->printDoc(
    "list of parts for pqJagged partitioning algorithm. As many as the dimension count.\n",
     docString);

  pList.set<string>(parameterName, "0", docString.str(), pqRangeRangeValidatorP);

  ////////// memory_procs
  intRangeValidatorP = rcp(new Zoltan2::IntegerRangeListValidator<int>(true));
  parameterName = string("memory_procs");

  docString.str("");
  intRangeValidatorP->printDoc(
    "list of ranks that memory profiling information (default \"0\")\n",
     docString);

  pList.set<string>(parameterName, "0", docString.str(), intRangeValidatorP);

  ////////// speed_versus_quality
  parameterName = string("speed_versus_quality");

  strValidatorP = rcp(new StringValidator(
    tuple<string>("speed", "balance", "quality")));

  docString.str("");
  strValidatorP->printDoc(
    "When algorithm choices exist, opt for speed or solution quality?\n"
    "(Default is a balance of speed and quality)\n",
     docString);

  pList.set<string>(parameterName, "balance", docString.str(), strValidatorP);

  ////////// memory_versus_speed
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

  ////////// random_seed
  parameterName = string("random_seed");

  anyNumValidatorP = rcp(new AnyNumberParameterEntryValidator);

  docString.str("");
  anyNumValidatorP->printDoc("random seed\n", docString);

  pList.set<string>(parameterName, "0.5", docString.str(), anyNumValidatorP);

  ////////// order_method
  parameterName = string("order_method");
  strValidatorP = rcp(new StringValidator(
    tuple<string>("rcm", "minimum_degree", "method3")));

  docString.str("");
  strValidatorP->printDoc(
    "Document the order_method parameter here\n"
    "(Default is ?)\n",
     docString);

  pList.set<string>(parameterName, "rcm", docString.str(), strValidatorP);

  ////////// order_package
  parameterName = string("order_package");
  strValidatorP = rcp(new StringValidator(
    tuple<string>("amd", "package2", "package3")));

  docString.str("");
  strValidatorP->printDoc(
    "Document the order_package parameter here\n"
    "(Default is ?)\n",
     docString);

  pList.set<string>(parameterName, "amd", docString.str(), strValidatorP);

 ////////// compute_metrics
  parameterName = string("compute_metrics");

  str2intValidatorP =
    rcp(new str2intValidator(yesNoStrings, yesNoIntegrals, parameterName));

  docString.str("");
  str2intValidatorP->printDoc("Compute metrics after computing solution\n",
    docString);

  pList.set<string>(parameterName, "no", docString.str(),
    str2intValidatorP);

  ////////// topology
  parameterName = string("topology");

  strValidatorP = rcp(new StringValidator);

  docString.str("");
  docString << "Topology of node to be used in hierarchical partitioning\n";
  docString << "  \"2,4\"  for dual-socket quad-core\n";
  docString << "  \"2,2,6\"  for dual-socket, dual-die, six-core\n";
  docString << "  \"2,2,3\"  for dual-socket, dual-die, six-core but\n";
  docString << "             with only three partitions per die\n";

  pList.set<string>(parameterName, "", docString.str(), strValidatorP);
  
  ////////// randomize_input
  parameterName = string("randomize_input");
    
  str2intValidatorP = 
    rcp(new str2intValidator(yesNoStrings, yesNoIntegrals, parameterName));

  docString.str("");
  str2intValidatorP->printDoc("randomize input prior to partitioning\n",
    docString);
  
  pList.set<string>(parameterName, "no", docString.str(),
    str2intValidatorP);

  ////////// objective
  parameterName = string("partitioning_objective");   // TODO

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

  pList.set<string>(parameterName, "balance_object_weight", 
    docString.str(), strValidatorP);

  ////////// imbalance_tolerance
  parameterName = string("imbalance_tolerance");

  anyNumValidatorP = rcp(new AnyNumberParameterEntryValidator);

  docString.str("");
  anyNumValidatorP->printDoc(
   "imbalance tolerance, ratio of maximum load over average load"
   " (default 1.1)\n", 
      docString);

  pList.set<string>(parameterName, "1.1", docString.str(), 
    anyNumValidatorP);
  
  ////////// num_global_parts
  parameterName = string("num_global_parts");

  anyNumValidatorP = rcp(new AnyNumberParameterEntryValidator);

  docString.str("");
  anyNumValidatorP->printDoc(
   "global number of parts to compute (default is number of processes)\n",
      docString);

  pList.set<string>(parameterName, "0", docString.str(), 
    anyNumValidatorP);
  
  ////////// num_local_parts
  parameterName = string("num_local_parts");

  anyNumValidatorP = rcp(new AnyNumberParameterEntryValidator);

  docString.str("");
  anyNumValidatorP->printDoc(
   "number of parts to compute for this process(default is one)\n",
      docString);

  pList.set<string>(parameterName, "0", docString.str(), 
    anyNumValidatorP);

  ////////// approach
  parameterName = string("partitioning_approach");   // TODO

  strValidatorP = rcp(new StringValidator(
    tuple<string>("partition", "repartition", "maximize_overlap")));

  docString.str("");
  strValidatorP->printDoc(
    "Partition from scratch, partition incrementally from current\n"
    "partition, of partition from scratch but maximize overlap\n"
    "with the current partition (default is \"partition\" from scratch)\n",
     docString);

  pList.set<string>(parameterName, "partition", docString.str(), 
    strValidatorP);
  
  ////////// objects
  parameterName = string("objects_to_partition");   // TODO

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

  pList.set<string>(parameterName, "graph_vertices", docString.str(), 
    strValidatorP);
  
  ////////// model
  parameterName = string("model");

  strValidatorP = rcp(new StringValidator(
    tuple<string>("hypergraph", "graph", "geometry", "ids")));

  docString.str("");
  strValidatorP->printDoc(
    "This is a low level parameter.  Normally the library will choose\n"
    "a computational model based on the algorithm or objective specified\n"
    "by the user.\n",
     docString);
  
  pList.set<string>(parameterName, "graph", docString.str(), 
    strValidatorP);
  
  ////////// algorithm
  parameterName = string("algorithm");

  strValidatorP = rcp(new StringValidator(
    tuple<string>(
      "rcb",
      "multijagged",
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
  
  pList.set<string>(parameterName, "random", docString.str(), 
    strValidatorP);

  ////////// rectilinear
  parameterName = string("rectilinear");

  str2intValidatorP =
    rcp(new str2intValidator(yesNoStrings, yesNoIntegrals, parameterName));

  docString.str("");
  str2intValidatorP->printDoc(
    "If true, then when a cut is made, all of the dots located on the cut\n"
    "are moved to the same side of the cut. The resulting regions are then\n"
    "rectilinear.  The resulting load balance may not be as good as when\n"
    "the group of dots is split by the cut. Default is false.\n",
    docString);

  pList.set<string>(parameterName, "no", docString.str(),
    str2intValidatorP);

  ////////// average_cuts
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

  pList.set<string>(parameterName, "no", docString.str(),
    str2intValidatorP);

  ////////// symmetrize_input
  parameterName = string("symmetrize_input");

  strValidatorP = rcp(new StringValidator(
    tuple<string>( "no", "transpose", "bipartite")));

  docString.str("");
  strValidatorP->printDoc(
    "Symmetrize input prior to pList.  If \"transpose\",\n"
    "symmetrize A by computing A plus ATranspose.  If \"bipartite\",\n"
    "A becomes [[0 A][ATranspose 0]].  \n",
    docString);

  pList.set<string>(parameterName, "no", docString.str(), 
    strValidatorP);

  ////////// subset_graph
  parameterName = string("subset_graph");
    
  str2intValidatorP = 
    rcp(new str2intValidator(yesNoStrings, yesNoIntegrals, parameterName));

  docString.str("");
  str2intValidatorP->printDoc(
    "If \"yes\", the graph input is to be subsetted.  If a vertex neighbor\n"
    "is not a valid vertex, it will be omitted from the pList.  Otherwise,\n"
    "an invalid neighbor identifier is considered an error.\n",
    docString);
  
  pList.set<string>(parameterName, "no", docString.str(),
    str2intValidatorP);
}
