// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Parameters.cpp
  
  \brief Set up the valid parameters, their documentation, 
        their defaults, and their validators.

*/
#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Parameters.hpp>

#include <ostream>
#include <sstream>
#include <Teuchos_StandardParameterEntryValidators.hpp>

using namespace std;

// TODO the ordering, coloring and non-graph partitioning parameters
//    and many documentation strings
//
// This is in a .cpp file because it's not templated.  If the function
// needs to be templated on, say, AppGID for fixed vertex id lists, then
// it will need to be in Zoltan2_Parameters_def.hpp.

namespace Zoltan2 {

void createValidParameterList(Teuchos::ParameterList &pl, const Comm<int> &comm)
{
  using Teuchos::EnhancedNumberValidator;
  using Teuchos::StringValidator;
  using Teuchos::StringToIntegralParameterEntryValidator;
  using Teuchos::AnyNumberParameterEntryValidator;
  using Teuchos::FileNameValidator;

  using std::string;

  bool isDefault=true, isNotDefault=false;
  bool isList=true, isNotList=false;
  string parameterName;
  std::ostringstream docString;
  ParameterEntry entry;
  int intDefault;
  double doubleDefault;
  string strDefault;
  Array<int> intArrayDefault;
  Array<string> validChoices;

  typedef StringToIntegralParameterEntryValidator<int> str2intValidator;

  RCP<const EnhancedNumberValidator<int> > intValidatorP;
  RCP<const EnhancedNumberValidator<size_t> > sizetValidatorP;
  RCP<const EnhancedNumberValidator<double> > doubleValidatorP;
  RCP<const IntegerRangeListValidator<int> > intRangeValidatorP;
  RCP<const StringValidator> strValidatorP;
  RCP<const AnyNumberParameterEntryValidator> anyNumValidatorP;
  RCP<const str2intValidator> str2intValidatorP;
  RCP<const FileNameValidator > fnameValidatorP;

  Array<string> yesNoWords;   // for string to int transformations
  Array<int> yesNoNumbers;
  yesNoWords.push_back("true");
  yesNoWords.push_back("yes");
  yesNoWords.push_back("1");
  yesNoNumbers.push_back(1);
  yesNoNumbers.push_back(1);
  yesNoNumbers.push_back(1);
  yesNoWords.push_back("false");
  yesNoWords.push_back("no");
  yesNoWords.push_back("0");
  yesNoNumbers.push_back(0);
  yesNoNumbers.push_back(0);
  yesNoNumbers.push_back(0);

  // levels for error checking at run-time

  std::ostringstream levels;
  levels << BASIC_ASSERTION << " - basic checking of function arguments\n";
  levels << COMPLEX_ASSERTION << 
    " - more in-depth checks (i.e. input graph is valid)\n";
  levels << DEBUG_MODE_ASSERTION << 
    " - check everything regardless of the time required\n";
  levels << "(To turn off all run-time checking,";
  levels << " compile with the Z2_OMIT_ALL_ERROR_CHECKING flag.)\n";
  string assertionLevelDoc(levels.str());

  // levels of verbosity for debugging, timing and memory usage messages

  levels.str("");
  levels << NO_STATUS << " - no status output\n";
  levels << BASIC_STATUS << " - basic status information\n";
  levels << DETAILED_STATUS << 
    " - detailed status output (i.e. methods entered and exited for debugging,";
  levels << "   or sub-steps in an algorithm for timing)\n";
  levels << VERBOSE_DETAILED_STATUS << 
    " - verbose debugging output (i.e. entire graph, entire solution)\n";
  string statusLevelDoc(levels.str());

  // does user want local status only or global reductions?

  levels.str("");
  levels << LOCAL_SUMMARY << " - display local status only\n";
  levels << GLOBAL_SUMMARY << " - gather global results for status reporting\n";
  string statusSummaryDoc(levels.str());

  pl.setName("validatingDefaultingZoltan2Parameters");

  ///////////////////////////////////////////////////////////
  // Library configuration parameters
  ///////////////////////////////////////////////////////////

  //   The user can request:
  //        debug messages
  //        timing messages
  //        memory profiling messages
  //
  //   Each type of message has: 
  //     a level of verbosity ("*_level")
  //     one or more processes that do the output ("*_procs")
  //     an output stream like std::cout and friends  ("*_output_stream")
  //     or an output file name ("*_output_file")

  string debugDoc("For debugging messages");
  string timingDoc("For timing messages");
  string memoryDoc("For memory usage messages");

    /*-----------------------------------------*/
    /* level of verbosity                      */
    /*-----------------------------------------*/

  intDefault = NO_STATUS;
  intValidatorP = Teuchos::rcp(
    new EnhancedNumberValidator<int>(0,NUM_STATUS_OUTPUT_LEVELS-1));
  docString.str("amount of information to display at run-time\n");
  docString << statusLevelDoc;

  parameterName = string("debug_level");  
  entry = ParameterEntry( intDefault, isDefault, isNotList,
    debugDoc + docString.str(), intValidatorP); 
  pl.setEntry(parameterName, entry);

  parameterName = string("timing_level");  
  entry = ParameterEntry( intDefault, isDefault, isNotList,
    timingDoc + docString.str(), intValidatorP);
  pl.setEntry(parameterName, entry);

  parameterName = string("memory_profiling_level");  
  entry = ParameterEntry( intDefault, isDefault, isNotList,
    memoryDoc + docString.str(), intValidatorP);
  pl.setEntry(parameterName, entry);

    /*-----------------------------------------*/
    /* which processes do the output?          */
    /*-----------------------------------------*/

  intArrayDefault = Array<int>(2);  
  intArrayDefault[0] = 0;                 // node zero only
  intArrayDefault[1] = RANGE_IS_LISTED;   // flag
  intRangeValidatorP = Teuchos::rcp(new IntegerRangeListValidator<int>());
  docString.str("");
  intRangeValidatorP->printDoc(
    string("the processes that will collect and output information\n"), 
    docString); 

  parameterName = string("debug_procs");
  entry = ParameterEntry(intArrayDefault, isDefault, isList,
    debugDoc + docString.str(), intRangeValidatorP);
  pl.setEntry(parameterName, entry);

  parameterName = string("timing_procs");
  entry = ParameterEntry( intArrayDefault, isDefault, isList,
    timingDoc + docString.str(), intRangeValidatorP);
  pl.setEntry(parameterName, entry);

  parameterName = string("memory_profiling_procs");
  entry = ParameterEntry( intArrayDefault, isDefault, isList,
    memoryDoc + docString.str(), intRangeValidatorP);

    /*-----------------------------------------*/
    /* output to stream or a file?             */
    /*-----------------------------------------*/

  strDefault = string("std::cout");
  validChoices.clear();
  validChoices.push_back("std::cout");
  validChoices.push_back("std::cerr");
  validChoices.push_back("/dev/null");
  strValidatorP = Teuchos::rcp(new StringValidator(validChoices));
  docString.str("");
  strValidatorP->printDoc(string("output stream for messages"),
    docString); 

  parameterName = string("debug_output_stream");
  entry = ParameterEntry( strDefault, isDefault, isNotList,
        debugDoc + docString.str(), strValidatorP);
  pl.setEntry(parameterName, entry);

  parameterName = string("timing_output_stream");
  entry = ParameterEntry( strDefault, isDefault, isNotList,
        timingDoc + docString.str(), strValidatorP);
  pl.setEntry(parameterName, entry);

  parameterName = string("memory_profiling_output_stream");
  entry = ParameterEntry( strDefault, isDefault, isNotList,
        memoryDoc+docString.str(), strValidatorP);
  pl.setEntry(parameterName, entry);

  strDefault = string();
  docString.str("name of file to which messages should be written");
  fnameValidatorP = Teuchos::rcp(new FileNameValidator(false));

  parameterName = string("debug_output_file");
  entry = ParameterEntry( strDefault, isDefault, isNotList,
     debugDoc + docString.str(), fnameValidatorP);
  pl.setEntry(parameterName, entry);

  parameterName = string("timing_output_file");
  entry = ParameterEntry( strDefault, isDefault, isNotList,
     timingDoc + docString.str(), fnameValidatorP);
  pl.setEntry(parameterName, entry);

  parameterName = string("memory_profiling_output_file");
  entry = ParameterEntry( strDefault, isDefault, isNotList,
     memoryDoc + docString.str(), fnameValidatorP);
  pl.setEntry(parameterName, entry);

  ///////////////////////////////////////////////////////////
  // General problem parameters
  ///////////////////////////////////////////////////////////

  parameterName = string("speed_versus_quality");
  validChoices.clear();
  validChoices.push_back("speed");
  validChoices.push_back("balance");
  validChoices.push_back("quality");
  strDefault = string("balance");
  strValidatorP = Teuchos::rcp(new StringValidator(validChoices));
  docString.str("todo");
  entry = ParameterEntry(
        strDefault, isDefault, isNotList,
        docString.str(), 
        strValidatorP);

  pl.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = string("memory_footprint_versus_runtime");
  validChoices.clear();
  validChoices.push_back("memory_footprint");
  validChoices.push_back("memory");
  validChoices.push_back("balance");
  validChoices.push_back("runtime");
  strDefault = string("balance");
  strValidatorP = Teuchos::rcp(new StringValidator(validChoices));
  docString.str("todo");
  entry = ParameterEntry(
        strDefault, isDefault, isNotList,
        docString.str(), 
        strValidatorP);

  pl.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = string("error_check_level");
  intDefault = BASIC_ASSERTION;
  intValidatorP = 
    Teuchos::rcp(new EnhancedNumberValidator<int>(0,NUM_ASSERTION_LEVELS-1));
  entry = ParameterEntry(intDefault, isDefault, isNotList,
        assertionLevelDoc, intValidatorP);
  pl.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = string("random_seed");
  doubleDefault = .5;
  // false: don't allow all types
  Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes types(false);
  types.allowInt();
  types.allowDouble();
  AnyNumberParameterEntryValidator::EPreferredType preference =
    AnyNumberParameterEntryValidator::PREFER_DOUBLE;
  anyNumValidatorP = 
    Teuchos::rcp(new AnyNumberParameterEntryValidator(preference, types));
  docString.str("todo");
  entry = ParameterEntry(
        doubleDefault, isDefault, isNotList,
        docString.str(), 
        anyNumValidatorP);

  pl.setEntry(parameterName, entry);

  ///////////////////////////////////////////////////////////
  // Partitioning problem parameters
  ///////////////////////////////////////////////////////////

  Teuchos::ParameterList &partitioning = 
    pl.sublist("partitioning", false, string("Partitioning problem parameters"));
    /*-----------------------------------------*/

  parameterName = string("topology");
  strDefault = string("all");
  intRangeValidatorP = rcp(new IntegerRangeListValidator<int>(1, 10000));
  docString.str("topology of compute node for hierarchical partitioning TODO");
  entry = ParameterEntry(
        strDefault, !isDefault, isNotList,
        docString.str(), 
        intRangeValidatorP);

  partitioning.setEntry(parameterName, entry);
  
    /*-----------------------------------------*/

  parameterName = string("randomize_input");
  intDefault = 0;
  str2intValidatorP = 
    Teuchos::rcp(new StringToIntegralParameterEntryValidator<int>(
      yesNoWords, yesNoNumbers, string("randomize_input")));
  docString.str("todo");
  entry = ParameterEntry(
        intDefault, isDefault, isNotList,
        docString.str(), 
        str2intValidatorP);

  partitioning.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = string("objective");
  strDefault = string("balance_object_weight");
  validChoices.clear();
  validChoices.push_back("balance_object_count");
  validChoices.push_back("balance_object_weight");
  validChoices.push_back("redistribution_commvolume_tradeoff");  // doesn't this need its own param
  validChoices.push_back("multicriteria_minimize_total_weight");
  validChoices.push_back("multicriteria_minimize_maximum_weight");
  validChoices.push_back("minimize_cut_edge_count");
  validChoices.push_back("minimize_cut_edge_weight");
  validChoices.push_back("minimize_neighboring_parts");
  validChoices.push_back("minimize_boundary_vertices");

  docString.str("todo");
  strValidatorP = Teuchos::rcp(new StringValidator(validChoices));
  entry = ParameterEntry(
        strDefault, isDefault, isList,
        docString.str(), 
        strValidatorP);

  partitioning.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = string("imbalance_tolerance");  
  doubleDefault = 1.1;
  doubleValidatorP = 
    Teuchos::rcp(new EnhancedNumberValidator<double>(0,2.0));
  docString.str("todo");
  entry = ParameterEntry(
        doubleDefault, isDefault, isNotList,
        docString.str(), 
        doubleValidatorP);

  partitioning.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  // num_global_parts, if set must be one or greater.  It is
  // processed in PartitioningSolution::setPartDistribution().
  // If set, num_global_parts must be set on all processes to
  // the same value.

  size_t veryLarge = size_t(1) << (sizeof(size_t)*8 - 2);
  size_t minVal(1);
  sizetValidatorP = 
    Teuchos::rcp(new EnhancedNumberValidator<size_t>(minVal,veryLarge));

  parameterName = string("num_global_parts");  
  docString.str("");
  entry = ParameterEntry(minVal, !isDefault, isNotList, docString.str(), 
    sizetValidatorP) ;

  partitioning.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  // num_local_parts, if set must be zero or greater.  It is
  // processed in PartitioningSolution::setPartDistribution().
  // If set on any process, it must be set on all processes.
  // If num_global_parts is also set, num_global_parts must be 
  // the sum of the value of num_local_parts on each process.

  minVal = 0;
  sizetValidatorP = 
    Teuchos::rcp(new EnhancedNumberValidator<size_t>(minVal,veryLarge));

  parameterName = string("num_local_parts");  
  docString.str("");
  entry = ParameterEntry(minVal, !isDefault, isNotList, docString.str(), 
    sizetValidatorP);

  partitioning.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = string("approach");
  strDefault = string("partition");
  validChoices.clear();
  validChoices.push_back("partition");
  validChoices.push_back("repartition");
  validChoices.push_back("maximize_overlap");
  docString.str("todo");
  strValidatorP = Teuchos::rcp(new StringValidator(validChoices));
  entry = ParameterEntry(
        strDefault, isDefault, isList,
        docString.str(), 
        strValidatorP);

  partitioning.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = string("objects");
  strDefault = string("no default");
  validChoices.clear();
  validChoices.push_back("matrix_rows");
  validChoices.push_back("matrix_columns");
  validChoices.push_back("matrix_nonzeros");
  validChoices.push_back("mesh_elements");
  validChoices.push_back("mesh_nodes");
  validChoices.push_back("graph_vertices");
  validChoices.push_back("coordinates");
  validChoices.push_back("identifiers");
  docString.str("todo");
  strValidatorP = Teuchos::rcp(new StringValidator(validChoices));
  entry = ParameterEntry(
        strDefault, isNotDefault, isList,
        docString.str(), 
        strValidatorP);

  partitioning.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = string("model");
  strDefault = string("no default");
  validChoices.clear();
  validChoices.push_back("hypergraph");
  validChoices.push_back("graph");
  validChoices.push_back("geometry");
  validChoices.push_back("ids");
  docString.str("todo");
  strValidatorP = Teuchos::rcp(new StringValidator(validChoices));
  entry = ParameterEntry(
        strDefault, isNotDefault, isList,
        docString.str(), 
        strValidatorP);

  partitioning.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = string("algorithm");
  strDefault = string("no default");
  validChoices.clear();
  validChoices.push_back("rcb");
  validChoices.push_back("rib");
  validChoices.push_back("hsfc");
  validChoices.push_back("patoh");
  validChoices.push_back("phg");
  validChoices.push_back("metis");
  validChoices.push_back("parmetis");
  validChoices.push_back("scotch");
  validChoices.push_back("ptscotch");
  validChoices.push_back("block");
  validChoices.push_back("cyclic");
  validChoices.push_back("random");
  docString.str("todo");
  strValidatorP = Teuchos::rcp(new StringValidator(validChoices));
  entry = ParameterEntry(
        strDefault, isNotDefault, isList,
        docString.str(), 
        strValidatorP);

  partitioning.setEntry(parameterName, entry);

  ///////////////////////////////////////////////////////////
  // Geometric partitioning problem parameters   TODO
  ///////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////
  // Hypergraph partitioning problem parameters    TODO
  ///////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////
  // Graph partitioning problem parameters     TODO
  ///////////////////////////////////////////////////////////

  Teuchos::ParameterList &graph= 
    partitioning.sublist("graph", false, 
      string("graph partitioning problem parameters"));

    /*-----------------------------------------*/

  parameterName = string("symmetrize_input");
  strDefault = string("no");
  validChoices.clear();
  validChoices.push_back("no");
  validChoices.push_back("transpose");
  validChoices.push_back("bipartite");
  docString.str("todo");
  strValidatorP = Teuchos::rcp(new StringValidator(validChoices));
  entry = ParameterEntry(
        strDefault, isNotDefault, isList,
        docString.str(), 
        strValidatorP);

  graph.setEntry(parameterName, entry);

  // TODO: Parmetis and Scotch parameters.
}


}  //namespace Zoltan2
