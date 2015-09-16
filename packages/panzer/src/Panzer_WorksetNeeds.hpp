#ifndef __Panzer_WorsetNeeds_hpp__
#define __Panzer_WorsetNeeds_hpp__

#include "Teuchos_RCP.hpp"

#include "Panzer_ConfigDefs.hpp"
#include "Panzer_CellData.hpp"

namespace panzer {

class PureBasis;
class IntegrationRule;

/** This class provides a simplified interface to the objects
  * required to specify a Workset. In paritcular this is all
  * "meta" data that describes which basis functions are need,
  * which integration rules are needed and the shape of the
  * cell.
  * 
  * This is intended to be specified for each element block
  * and side set based on the integration rules and basis functions
  * that are needed.
  */ 
struct WorksetNeeds {
  CellData cellData;
  std::vector<Teuchos::RCP<const IntegrationRule> > int_rules;
  std::vector<Teuchos::RCP<const PureBasis> > bases;
  std::vector<std::string> rep_field_name; // representative field name
};

} // end namespace panzer

#endif
