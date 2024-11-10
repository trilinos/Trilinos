// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_CommaSeparatedEntryValidator.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Teuchos_StrUtils.hpp"

namespace panzer {

void CommaSeparatedEntryValidator::
split(const std::string & str,
      const std::string & delim,
      std::vector<std::string> & output)
{
   output.clear();

   // typedef boost::tokenizer<boost::char_separator<char> > 
   //         tokenizer;

   // boost::char_separator<char> sep(delim.c_str());
   // tokenizer tokens(str, sep);
   // for(tokenizer::iterator tok_iter = tokens.begin();
   //     tok_iter != tokens.end(); ++tok_iter) {
   //    // extract token, remove spaces
   //    std::string s = *tok_iter;
   //    boost::trim(s);
   //    if(s.length()!=0)
   //       output.push_back(s);
   // }

   panzer::StringTokenizer(output, str, delim, true);
}

void 
CommaSeparatedEntryValidator::
validate(const Teuchos::ParameterEntry & entry,  
	 const std::string & paramName,
	 const std::string & sublistName) const
{
  const std::string &entryName = entry.getAny(false).typeName();
  Teuchos::any anyValue = entry.getAny(true);
  
  // type passed, validate value
  TEUCHOS_TEST_FOR_EXCEPTION(!(anyValue.type() == typeid(std::string) ),
    Teuchos::Exceptions::InvalidParameterType,
    "Sorry but it looks like the \"" << paramName << "\"" <<
    " parameter in the \"" << sublistName <<
    "\" sublist does not exist." << std::endl << std::endl <<
    "Error: The value that you entered was the wrong type." << std::endl <<
    "Parameter: " << paramName << std::endl <<
    "Type specified: " << entryName << std::endl <<
    "Type accepted: " << typeid(std::string).name() <<
    std::endl << std::endl);

  const std::string & value = Teuchos::any_cast<std::string>(anyValue);

  std::vector<std::string> tokens;
  split(value,",",tokens);  

  if(!allowEmpty_) {
     const std::string errorStr = "The value for \"string-list\" type parameter in sublist \""+sublistName+"\" named \""+paramName+"\" "
                                  "is incorrectly formatted. The expected format is\n"
                                  "   \"<string>[, <string>]*\" "
                                  "your value is \""+value+"\"";

     // verify that their is a response type and an evaluation type
     TEUCHOS_TEST_FOR_EXCEPTION(tokens.size()==0,
        Teuchos::Exceptions::InvalidParameterValue,errorStr);
  }
}


void CommaSeparatedEntryValidator::printDoc(
  std::string const &docString, std::ostream &out) const
{
  Teuchos::StrUtils::printLines(out,"# ",docString);
  out << "#  Validator Used: " << std::endl;
  out << "#  CommaSeparatedEntry Validator" << std::endl;
}

}
