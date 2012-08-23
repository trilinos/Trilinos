// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include "Panzer_ResponseUtilities.hpp"

#include "Teuchos_StrUtils.hpp"

namespace panzer {

void CommaSeperatedEntryValidator::split(const std::string & str,const std::string & delim,
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
      if(s.length()!=0)
         output.push_back(s);
   }
}

void CommaSeperatedEntryValidator::validate(const Teuchos::ParameterEntry & entry,  
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


void CommaSeperatedEntryValidator::printDoc(
  std::string const &docString, std::ostream &out) const
{
  Teuchos::StrUtils::printLines(out,"# ",docString);
  out << "#  Validator Used: " << std::endl;
  out << "#  CommaSeperatedEntry Validator" << std::endl;
}

void buildResponseMap(const Teuchos::ParameterList & p,std::map<std::string,std::pair<ResponseId,std::set<std::string> > > & responses)
{
   TEUCHOS_ASSERT(false);
}

void buildResponseMap(const Teuchos::ParameterList & p,
                      std::map<std::string,std::pair<ResponseId,std::pair<std::list<std::string>,std::list<std::string> > > > & responses)
{
   static Teuchos::RCP<const Teuchos::ParameterList> validList;

   // build valid parameter list
   if(validList==Teuchos::null) {
      Teuchos::RCP<Teuchos::ParameterList> tmpList = Teuchos::rcp(new Teuchos::ParameterList);
      tmpList->set<std::string>("Type","");
      tmpList->set<std::string>("Field Name","");
      tmpList->set<std::string>("Element Blocks","empty","Element blocks for this response",Teuchos::rcp(new CommaSeperatedEntryValidator));
      tmpList->set<std::string>("Evaluation Types","empty","Evaluation types for this response",Teuchos::rcp(new CommaSeperatedEntryValidator));

      validList = tmpList;
   }
 
   CommaSeperatedEntryValidator validator;
   const std::string & sublistName = p.name();
   std::vector<std::string> tokens;

   responses.clear();

   // loop over entries of parameter list, must satisfy response formatting
   for(Teuchos::ParameterList::ConstIterator itr=p.begin(); itr!=p.end();++itr) {
 
      const std::string & paramName = itr->first;
      const Teuchos::ParameterEntry & pe = itr->second;

      // make sure this is a parameter list
      TEUCHOS_TEST_FOR_EXCEPTION(!pe.isList(),Teuchos::Exceptions::InvalidParameterValue,
                         "In list \""+sublistName+"\", the parameter \""+paramName+"\" is expected "
                         "to be a sublist. Response map cannot be built!");

      // extract parameter list and validate
      const Teuchos::ParameterList & respList = Teuchos::getValue<Teuchos::ParameterList>(pe); 
      respList.validateParameters(*validList);

      const std::string & respLabel = paramName;
      ResponseId & rid = responses[respLabel].first;
      std::list<std::string> & eBlocks = responses[respLabel].second.first; // element blocks
      std::list<std::string> & eTypes = responses[respLabel].second.second;  // evaluation types

      rid.type = respList.get<std::string>("Type");
      rid.name = respList.get<std::string>("Field Name");

      CommaSeperatedEntryValidator::split(respList.get<std::string>("Element Blocks"),",",tokens);
      eBlocks.assign(tokens.begin(),tokens.end()); // this should automatically wipe out old values
      
      CommaSeperatedEntryValidator::split(respList.get<std::string>("Evaluation Types"),",",tokens);
      eTypes.assign(tokens.begin(),tokens.end()); // this should automatically wipe out old values
   }
}

}
