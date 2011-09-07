#include "Panzer_ResponseUtilities.hpp"

#include "Teuchos_StrUtils.hpp"

namespace panzer {

void ResponseEntryValidator::split(const std::string & str,const std::string & delim,
                                   std::vector<std::string> & output)
{
   output.clear();

   typedef boost::tokenizer<boost::char_separator<char> > 
           tokenizer;

   boost::char_separator<char> sep(delim.c_str());
   tokenizer tokens(str, sep);
   for(tokenizer::iterator tok_iter = tokens.begin();
       tok_iter != tokens.end(); ++tok_iter) {
      // extract token, remove spaces
      std::string s = *tok_iter;
      boost::trim(s);
      output.push_back(s);
   }
}

void ResponseEntryValidator::validate(const Teuchos::ParameterEntry & entry,  
                                      const std::string & paramName,
                                      const std::string & sublistName) const
{
  const std::string &entryName = entry.getAny(false).typeName();
  Teuchos::any anyValue = entry.getAny(true);
  
  // type passed, validate value
  TEST_FOR_EXCEPTION(!(anyValue.type() == typeid(std::string) ),
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
  split(value,":",tokens);  

  const std::string errorStr = "The value for \"response\" type parameter named \""+paramName+"\" "
                               "is incorrectly formatted. The expected format is\n"
                               "   \"<Response type>(<Field name>): <Evaluation type> [, <Evaluation type>]*\""
                               "your value is \""+value+"\"";

  // verify that their is a response type and a evaluation type
  TEST_FOR_EXCEPTION(tokens.size()!=2,
     Teuchos::Exceptions::InvalidParameterValue,errorStr);

  std::string respType = tokens[0];
  std::string evalTypes = tokens[1];

  TEST_FOR_EXCEPTION(respType.length()==0 || evalTypes.length()==0,
     Teuchos::Exceptions::InvalidParameterValue,errorStr);

  // parse response: <Response Type>(<Field name>)
  split(respType,"()",tokens);
  TEST_FOR_EXCEPTION(tokens.size()!=2,
     Teuchos::Exceptions::InvalidParameterValue,errorStr);

  respType = tokens[0];
  std::string respName = tokens[1];

  TEST_FOR_EXCEPTION(respType.length()==0 || respName.length()==0,
     Teuchos::Exceptions::InvalidParameterValue,errorStr);

  // parse evaluation fields: <Evaluation type> [,<Evaluation type>]*
  split(evalTypes,",",tokens);  
  TEST_FOR_EXCEPTION(tokens.size()==0,
     Teuchos::Exceptions::InvalidParameterValue,errorStr);
}


void ResponseEntryValidator::printDoc(
  std::string const &docString, std::ostream &out) const
{
  Teuchos::StrUtils::printLines(out,"# ",docString);
  out << "#  Validator Used: " << std::endl;
  out << "#  ResponseEntry Validator" << std::endl;
}

void buildResponseMap(const Teuchos::ParameterList & p,std::map<std::string,std::pair<ResponseId,std::set<std::string> > > & responses)
{
   ResponseEntryValidator validator;
   const std::string & sublistName = p.name();
   std::vector<std::string> tokens;

   responses.clear();

   // loop over entries of parameter list, must satisfy response formatting
   for(Teuchos::ParameterList::ConstIterator itr=p.begin();
       itr!=p.end();++itr) {
 
      const std::string & paramName = itr->first;
      const Teuchos::ParameterEntry & pe = itr->second;

      // validate value is formatted correctly: if so (this doesn't throw)
      // then we can grind on w/o a care in the world!
      validator.validate(pe,paramName,sublistName);

      // extract value string for parsing
      Teuchos::any anyValue = pe.getAny(true);
      const std::string & value = Teuchos::any_cast<std::string>(anyValue);

      // parse the value string
      ResponseEntryValidator::split(value,":",tokens);  
      std::string respTypeAndName = tokens[0];
      std::string evalTypes = tokens[1];

      ResponseEntryValidator::split(respTypeAndName,"()",tokens);  
      std::string respType = tokens[0];
      std::string respName = tokens[1];
      
      ResponseEntryValidator::split(evalTypes,",",tokens);  

      // build response id, and evaluation set pair
      std::pair<ResponseId,std::set<std::string> > & respPair = responses[paramName];
    
      respPair.first.name = respName;
      respPair.first.type = respType; 
      respPair.second.insert(tokens.begin(),tokens.end()); 
   }
}

}
