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

#include "Panzer_STK_PeriodicBC_Parser.hpp"

#include "Panzer_STK_PeriodicBC_MatchConditions.hpp"

#include "Teuchos_ParameterListExceptions.hpp"

namespace panzer_stk {

PeriodicBC_Parser::PeriodicBC_Parser()
   : countStr_("Count")
   , condPrefix_("Periodic Condition ")
{
}

const std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > &
PeriodicBC_Parser::getMatchers() const
{
   return matchers_;
}

void PeriodicBC_Parser::setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pl)
{
   if(not pl->isParameter(countStr_)) {
      bool validEntry = false;
      TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(
        !validEntry, Teuchos::Exceptions::InvalidParameterName,
        "Error, the parameter {name=\"" << countStr_ << "\","
        "type=\"int\""
        "\nis required in parameter (sub)list \""<< pl->name() <<"\"."
        "\n\nThe valid parameters and types are:\n"
        << getValidParameters()->currentParametersString() << "\n\n"
        << "Passed parameter list: \n" << pl->currentParametersString()
      );
   }

   int numBCs = pl->get<int>(countStr_);
   pl->validateParameters(*getValidParameters(numBCs));

   // loop over boundary conditions
   for(int i=1;i<=numBCs;i++) {
      std::stringstream ss;

      ss << condPrefix_ << i;
      std::string cond = pl->get<std::string>(ss.str());
      matchers_.push_back(buildMatcher(cond));
   }

   storedPL_ = pl;
}

Teuchos::RCP<Teuchos::ParameterList>  PeriodicBC_Parser::getNonconstParameterList()
{
   return storedPL_;
}

Teuchos::RCP<Teuchos::ParameterList>  PeriodicBC_Parser::unsetParameterList()
{
   Teuchos::RCP<Teuchos::ParameterList> pl = storedPL_;
   storedPL_ = Teuchos::null;
   return pl;
}

Teuchos::RCP<const Teuchos::ParameterList> PeriodicBC_Parser::getValidParameters() const
{
   static Teuchos::RCP<Teuchos::ParameterList> pl;

   // build a sample parameter list with a single preconditioner
   if(pl==Teuchos::null) {
      std::stringstream ss;
      ss << condPrefix_ << 1 << std::endl;

      pl = Teuchos::rcp(new Teuchos::ParameterList);
      pl->set<int>(countStr_,1,
                   "Number of set periodic boundary conditions");
      pl->set<std::string>(ss.str(),"MatchCondition bndry1;bndry2",
                   "Boundary condition fairs formatted: <MatchCondition> <bndry1>;<bndry2>");
   }

   return pl.getConst();
}

Teuchos::RCP<Teuchos::ParameterList> PeriodicBC_Parser::getValidParameters(int count) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(count<0,std::logic_error,
                      "PeriodicBC requires a positive number (or none) of periodic boundary conditions.");

   Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
   pl->set(countStr_,count);

   for(int k=1;k<=count;k++) {
      std::stringstream ss;
      ss << condPrefix_ << k;

      pl->set<std::string>(ss.str(),"MatchCondition bndry1;bndry2");
   }

   return pl;
}

// basic string utilities to help wit parsing (only local)
/////////////////////////////////////////////////////////////

static std::string trim_left(const std::string & s)
{
   std::string::size_type beg = s.find_first_not_of(' ');

   return s.substr(beg,s.length()-beg);
}

static std::string trim_right(const std::string & s)
{
   std::string::size_type end = s.find_last_not_of(' ');

   return s.substr(0,end+1);
}

static std::string trim(const std::string & s)
{
   return trim_right(trim_left(s));
}

/////////////////////////////////////////////////////////////

void PeriodicBC_Parser::buildMatcher_Tokenize(const std::string & buildStr,
                                             std::string & matcher,
                                             std::string & bndry1,
                                             std::string & bndry2) const
{
   std::string::size_type endMatch = buildStr.find_first_of(' ');
   std::string::size_type begBndry = buildStr.find_first_of(';');

   matcher = trim(buildStr.substr(0,endMatch));
   bndry1 = trim(buildStr.substr(endMatch,begBndry-endMatch));
   bndry2 = trim(buildStr.substr(begBndry+1,buildStr.length()));
}

bool PeriodicBC_Parser::buildMatcher_Tokenize_withParams(const std::string & buildStr,
                                                         std::string & matcher,
                                                         std::vector<std::string> & params,
                                                         std::string & bndry1,
                                                         std::string & bndry2) const
{
   std::string::size_type endMatchAndParams = buildStr.find_first_of(':');
   std::string::size_type begBndry = buildStr.find_first_of(';');

   // no parameters: default to old style input
   if(endMatchAndParams==std::string::npos) {
      buildMatcher_Tokenize(buildStr,matcher,bndry1,bndry2); 
      return false;
   }

   bndry1 = trim(buildStr.substr(endMatchAndParams+1,begBndry-(endMatchAndParams+1)));
   bndry2 = trim(buildStr.substr(begBndry+1,buildStr.length()));

   std::string matchAndParams = trim(buildStr.substr(0,endMatchAndParams));
   std::string::size_type endMatch = matchAndParams.find_first_of(' ');

   // no parameters included
   if(endMatch==std::string::npos) {
      matcher = matchAndParams;
      return true;
   }

   // find parameters
   /////////////////////////////////////////////////////////

   // check matching conditions
   matcher = trim(matchAndParams.substr(0,endMatch));
   matchAndParams = matchAndParams.substr(endMatch+1);

   std::string::size_type comma = matchAndParams.find_first_of(',');
   while(comma!=std::string::npos) {
      std::string p = trim(matchAndParams.substr(0,comma));

      TEUCHOS_TEST_FOR_EXCEPTION(p.length()<1,std::logic_error,
                          "Error parsing periodic boundary condition \"" + buildStr + "\"");

      params.push_back(p);
      matchAndParams = matchAndParams.substr(comma+1);
      comma = matchAndParams.find_first_of(',');
   }

   std::string finalParam = trim(matchAndParams);
   if(finalParam.length()>0)
      params.push_back(finalParam);

   return true;
}

Teuchos::RCP<const PeriodicBC_MatcherBase> 
PeriodicBC_Parser::buildMatcher(const std::string & buildStr) const
{
   std::vector<std::string> params;
   std::string matcher, bndry1, bndry2;

   buildMatcher_Tokenize_withParams(buildStr,matcher,params,bndry1,bndry2);

   if(matcher=="x-coord") {
     panzer_stk::CoordMatcher matcher(0,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher);
   }

   if(matcher=="y-coord") {
     panzer_stk::CoordMatcher matcher(1,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher);
   }

   if(matcher=="z-coord") {
     panzer_stk::CoordMatcher matcher(2,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher);
   }

   if(matcher=="xy-coord" || matcher=="yx-coord") {
     panzer_stk::PlaneMatcher matcher(0,1,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher);
   }

   if(matcher=="xz-coord" || matcher=="zx-coord") {
     panzer_stk::PlaneMatcher matcher(0,2,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher);
   }

   if(matcher=="yz-coord" || matcher=="zy-coord") {
     panzer_stk::PlaneMatcher matcher(1,2,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher);
   }

   TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
       "Failed parsing parameter list: could not find periodic boundary "
       "condition matcher \"" << matcher << "\" "
       "in string \"" << buildStr << "\"");

   return Teuchos::null;
}

}
