// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Parameters.cpp
  
  \brief Set up the valid parameters, their defaults, their validators.

*/

#include <Zoltan2_Parameters.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

// TODO the docString for each parameter.  
// TODO the ordering, coloring and non-graph partitioning parameters
//
// This is in a .cpp file because it's not templated.  If the function
// needs to be templated on, say, AppGID for fixed vertex id lists, then
// it will need to be in Zoltan2_Parameters_def.hpp.

namespace Zoltan2 {

void createValidParameterList(Teuchos::ParameterList &pl)
{
  using Teuchos::EnhancedNumberValidator;
  using Teuchos::StringValidator;
  using Teuchos::StringToIntegralParameterEntryValidator;
  using Teuchos::AnyNumberParameterEntryValidator;
  using Teuchos::FileNameValidator;
  using Teuchos::ParameterEntry;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::Array;

  bool isDefault=true, isNotDefault=false;
  bool isList=true, isNotList=false;
  std::string parameterName;
  std::ostringstream docString;
  ParameterEntry entry;
  int intDefault;
  double doubleDefault;
  std::string strDefault;
  Array<int> intArrayDefault;
  Array<std::string> validChoices;
  RCP<const EnhancedNumberValidator<int> > intValidatorP;
  RCP<const EnhancedNumberValidator<double> > doubleValidatorP;
  RCP<const IntegerRangeListValidator<int> > intRangeValidatorP;
  RCP<const StringValidator> strValidatorP;
  RCP<const AnyNumberParameterEntryValidator> anyNumValidatorP;
  RCP<const StringToIntegralParameterEntryValidator<int> > str2intValidatorP;
  RCP<const FileNameValidator > fnameValidatorP;

  Array<std::string> yesNoWords;
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

  ///////////////////////////////////////////////////////////
  // Library configuration parameters
  ///////////////////////////////////////////////////////////

  pl.setName("validatingDefaultingZoltan2Parameters");

    /*-----------------------------------------*/

  parameterName = std::string("debug_level");  
  intDefault = 0;
  intValidatorP = 
    Teuchos::rcp(new EnhancedNumberValidator<int>(0,NUM_DEBUG_LEVELS-1));
  docString.str("");
  docString << "level of debugging output, an integer ranging from 0 to ";
  docString << NUM_DEBUG_LEVELS-1 << ".";
  entry = ParameterEntry(
        intDefault, isDefault, isNotList,
        docString.str(), 
        intValidatorP);

  pl.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = std::string("debug_procs");
  intArrayDefault = Array<int>(2);  
  intArrayDefault[0] = 0;
  intArrayDefault[1] = RANGE_IS_LISTED;   // flag
  intRangeValidatorP = Teuchos::rcp(new IntegerRangeListValidator<int>());
  docString.str("");
  intRangeValidatorP->printDoc(
    std::string("the processes that will display debugging output"), 
    docString); 
  entry = ParameterEntry(
        intArrayDefault, isDefault, isList,
        docString.str(), 
        intRangeValidatorP);

  pl.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = std::string("debug_ostream");
  strDefault = std::string("std::cout");
  docString.str("todo");
  validChoices.clear();
  validChoices.push_back("std::cout");
  validChoices.push_back("std::cerr");
  validChoices.push_back("/dev/null");
  strValidatorP = Teuchos::rcp(new StringValidator(validChoices));
  entry = ParameterEntry(
        strDefault, isDefault, isList,
        docString.str(), 
        strValidatorP);

  pl.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = std::string("debug_file");
  strDefault = std::string();
  docString.str("todo");
  fnameValidatorP = Teuchos::rcp(new FileNameValidator(false));
  entry = ParameterEntry(
        strDefault, isDefault, isList,
        docString.str(), 
        fnameValidatorP);

  pl.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  // TODO - I think we really want timing specifiers.
  parameterName = std::string("timing_level");
  intDefault = 0;
  intValidatorP = 
    Teuchos::rcp(new EnhancedNumberValidator<int>(0,NUM_TIMING_LEVELS-1));
  docString.str("");
  docString << "level of profiling output, an integer ranging from 0 to ";
  docString << NUM_TIMING_LEVELS-1 << ".";
  entry = ParameterEntry(
        intDefault, isDefault, isNotList,
        docString.str(), 
        intValidatorP);

  pl.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = std::string("timing_procs");
  intArrayDefault = Array<int>(1,RANGE_INCLUDES_ALL);  
  intRangeValidatorP = Teuchos::rcp(new IntegerRangeListValidator<int>());
  docString.str("");
  intRangeValidatorP->printDoc(
    std::string("the processes that will display profiling output"), 
    docString); 
  entry = ParameterEntry(
        intArrayDefault, isDefault, isList,
        docString.str(), 
        intRangeValidatorP);

  pl.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = std::string("timing_ostream");
  strDefault = std::string("std::cout");
  docString.str("todo");
  validChoices.clear();
  validChoices.push_back("std::cout");
  validChoices.push_back("std::cerr");
  validChoices.push_back("/dev/null");
  strValidatorP = Teuchos::rcp(new StringValidator(validChoices));
  entry = ParameterEntry(
        strDefault, isDefault, isList,
        docString.str(), 
        strValidatorP);

  pl.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = std::string("timing_file");
  strDefault = std::string();
  docString.str("todo");
  fnameValidatorP = Teuchos::rcp(new FileNameValidator(false));
  entry = ParameterEntry(
        strDefault, isDefault, isList,
        docString.str(), 
        fnameValidatorP);

  pl.setEntry(parameterName, entry);


  ///////////////////////////////////////////////////////////
  // General problem parameters
  ///////////////////////////////////////////////////////////

  parameterName = std::string("speed_versus_quality");
  validChoices.clear();
  validChoices.push_back("speed");
  validChoices.push_back("balance");
  validChoices.push_back("quality");
  strDefault = std::string("balance");
  strValidatorP = Teuchos::rcp(new StringValidator(validChoices));
  docString.str("todo");
  entry = ParameterEntry(
        strDefault, isDefault, isNotList,
        docString.str(), 
        strValidatorP);

  pl.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = std::string("memory_footprint_versus_runtime");
  validChoices.clear();
  validChoices.push_back("memory_footprint");
  validChoices.push_back("memory");
  validChoices.push_back("balance");
  validChoices.push_back("runtime");
  strDefault = std::string("balance");
  strValidatorP = Teuchos::rcp(new StringValidator(validChoices));
  docString.str("todo");
  entry = ParameterEntry(
        strDefault, isDefault, isNotList,
        docString.str(), 
        strValidatorP);

  pl.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = std::string("error_check_level");
  intDefault = 1;
  intValidatorP = 
    Teuchos::rcp(new EnhancedNumberValidator<int>(0,NUM_ASSERTION_LEVELS-1));
  docString.str("");
  docString << "The level of error checking to perform at runtime.\n";
  docString << "valid range: 0 - " << NUM_ASSERTION_LEVELS-1 << "\n";
  docString << "default level: " << intDefault << "\n";
  entry = ParameterEntry(
        intDefault, isDefault, isNotList,
        docString.str(), 
        intValidatorP);

  pl.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = std::string("random_seed");
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

  ParameterList &partitioning = 
    pl.sublist("partitioning", false, std::string("Partitioning problem parameters"));

    /*-----------------------------------------*/

  parameterName = std::string("randomize_input");
  intDefault = 0;
  str2intValidatorP = 
    Teuchos::rcp(new StringToIntegralParameterEntryValidator<int>(
      yesNoWords, yesNoNumbers, std::string("randomize_input")));
  docString.str("todo");
  entry = ParameterEntry(
        intDefault, isDefault, isNotList,
        docString.str(), 
        str2intValidatorP);

  partitioning.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = std::string("objective");
  strDefault = std::string("balance_object_weight");
  validChoices.clear();
  validChoices.push_back("balance_object_count");
  validChoices.push_back("balance_object_weight");
  validChoices.push_back("redistribution_commvolume_tradeoff");  // doesn't this need it's own param
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

  parameterName = std::string("use_node_topology");
  intDefault = 0;
  str2intValidatorP = 
    Teuchos::rcp(new StringToIntegralParameterEntryValidator<int>(
      yesNoWords, yesNoNumbers, std::string("randomize_input")));
  docString.str("todo");
  entry = ParameterEntry(
        intDefault, isDefault, isNotList,
        docString.str(), 
        str2intValidatorP);

  partitioning.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = std::string("imbalance_tolerance");  
  doubleDefault = 15.0;
  doubleValidatorP = 
    Teuchos::rcp(new EnhancedNumberValidator<double>(0,100));
  docString.str("todo");
  entry = ParameterEntry(
        doubleDefault, isDefault, isNotList,
        docString.str(), 
        doubleValidatorP);

  partitioning.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = std::string("num_global_parts");  
  intValidatorP = 
    Teuchos::rcp(new EnhancedNumberValidator<int>(-1,INT_MAX));
  docString.str("");
  entry = ParameterEntry(
        0, isNotDefault, isNotList,
        docString.str(), 
        intValidatorP);

  partitioning.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = std::string("num_local_parts");  
  intDefault = 1;
  intValidatorP = 
    Teuchos::rcp(new EnhancedNumberValidator<int>(-1,INT_MAX));
  docString.str("");
  entry = ParameterEntry(
        intDefault, isDefault, isNotList,
        docString.str(), 
        intValidatorP);

  partitioning.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = std::string("approach");
  strDefault = std::string("partition");
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

  parameterName = std::string("objects");
  strDefault = std::string("no default");
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

  parameterName = std::string("model");
  strDefault = std::string("no default");
  validChoices.clear();
  validChoices.push_back("hypergraph");
  validChoices.push_back("graph");
  validChoices.push_back("geometry");
  docString.str("todo");
  strValidatorP = Teuchos::rcp(new StringValidator(validChoices));
  entry = ParameterEntry(
        strDefault, isNotDefault, isList,
        docString.str(), 
        strValidatorP);

  partitioning.setEntry(parameterName, entry);

    /*-----------------------------------------*/

  parameterName = std::string("algorithm");
  strDefault = std::string("no default");
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

  ParameterList &graph= 
    partitioning.sublist("graph", false, 
      std::string("graph partitioning problem parameters"));

    /*-----------------------------------------*/

  parameterName = std::string("symmetrize_input");
  strDefault = std::string("no");
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
