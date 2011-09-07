#ifndef __Panzer_ResponseUtilities_hpp__
#define __Panzer_ResponseUtilities_hpp__

#include <map>
#include <set>
#include <vector>
#include <string>
#include <iostream>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_ParameterEntryValidator.hpp"
#include "Teuchos_RCP.hpp"

#include "Panzer_Response.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

namespace panzer {

/** This class validates a response type. Essentially
  * it is used to make sure the parameter value is correctly
  * formatted. 
  */
class ResponseEntryValidator : public Teuchos::ParameterEntryValidator {
public:
  ResponseEntryValidator() {}

  ValidStringsList validStringValues() const
  { return Teuchos::null; }

  void validate(const Teuchos::ParameterEntry & entry, 
                const std::string & paramName,
                const std::string & sublistName) const;

  const std::string getXMLTypeName() const
  { return "response"; }

  void printDoc(const std::string & docString, std::ostream &out) const;

  //! Utility function for tokenizing
  static void split(const std::string & str,
                    const std::string & delim,
                    std::vector<std::string> & tokens);
};

/** Given a parameter list, loop over all the entries
  * and parse their value field and build a response ID field and a 
  * set of types to evaluate.
  */
void buildResponseMap(const Teuchos::ParameterList & p,
                      std::map<std::string,std::pair<ResponseId,std::set<std::string> > > & responses);


}

#endif
