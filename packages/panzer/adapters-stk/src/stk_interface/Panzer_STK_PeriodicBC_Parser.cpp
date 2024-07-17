// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_STK_PeriodicBC_Parser.hpp"

#include "Panzer_STK_PeriodicBC_MatchConditions.hpp"

#include "Teuchos_ParameterListExceptions.hpp"

namespace panzer_stk {

PeriodicBC_Parser::PeriodicBC_Parser()
   : countStr_("Count")
   , condPrefix_("Periodic Condition ")
   , searchStr_("Use BBox Search")
   , useBBoxSearch_(false) // TODO swap this to change default search (see also STK_interface.hpp)
{
}

const std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > &
PeriodicBC_Parser::getMatchers() const
{
   return matchers_;
}

const bool & PeriodicBC_Parser::useBoundingBoxSearch() const
{
   return useBBoxSearch_;
}

void PeriodicBC_Parser::setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pl)
{
   if(!pl->isParameter(countStr_)) {
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
   useBBoxSearch_ = pl->get<bool>(searchStr_,useBBoxSearch_);
   pl->validateParameters(*getValidParameters(numBCs));

   // loop over boundary conditions
   for(int i=1;i<=numBCs;i++) {
      std::stringstream ss;

      ss << condPrefix_ << i;
      std::string cond = pl->get<std::string>(ss.str());

      std::pair<std::string, unsigned int> matcherPair = getMatcherTypeAndDim(cond);
      std::string matcherType = matcherPair.first;
      unsigned int matcherDim = matcherPair.second;
      if(matcherType == "coord"){
        matchers_.push_back(buildMatcher(cond));
      }else if(matcherType == "edge")
        edgeMatchers_.push_back(buildMatcher(cond));
      else if(matcherType == "face")
        faceMatchers_.push_back(buildMatcher(cond));
      else if(matcherType == "all"){
        matchers_.push_back(buildMatcher(replaceMatcherType(cond,"coord"))); 
        edgeMatchers_.push_back(buildMatcher(replaceMatcherType(cond,"edge")));
        if(matcherDim > 2)
          faceMatchers_.push_back(buildMatcher(replaceMatcherType(cond,"face"))); 
      } 
   }

   // Order BCs with all coords first, followed by edges, then faces
   matchers_.insert(matchers_.end(),edgeMatchers_.begin(),edgeMatchers_.end());
   matchers_.insert(matchers_.end(),faceMatchers_.begin(),faceMatchers_.end());

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
                   "Boundary condition pairs formatted: <MatchCondition> <bndry1>;<bndry2>");
      // TODO swap this to change default search (see also STK_interface.hpp)
      pl->set<bool>(searchStr_,false,"Use bounding box search for GID match (requires STKSearch) or not");
   }

   return pl.getConst();
}

Teuchos::RCP<Teuchos::ParameterList> PeriodicBC_Parser::getValidParameters(int count) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(count<0,std::logic_error,
                      "PeriodicBC requires a positive number (or none) of periodic boundary conditions.");

   Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList);
   pl->set(countStr_,count);
   // TODO swap this to change default search (see also STK_interface.hpp)
   pl->set<bool>(searchStr_,false,"Use bounding box search for GID match (requires STKSearch) or not");

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

std::pair<std::string, unsigned int> PeriodicBC_Parser::getMatcherTypeAndDim(const std::string & buildStr) const
{
   std::string::size_type endMatch = buildStr.find_first_of(' ');

   std::string matcher = trim(buildStr.substr(0,endMatch));

   std::string::size_type hyphenMatch = matcher.find_last_of('-');

   TEUCHOS_TEST_FOR_EXCEPTION(hyphenMatch==std::string::npos,std::logic_error,
       "Failed parsing parameter list: could not find periodic boundary "
       "condition matcher \"" << matcher << "\" "
       "in string \"" << buildStr << "\n"
       "Matcher " << matcher << " requires a hyphen, e.g. x-coord, yz-edge\"");

   std::string matcherType = trim(matcher.substr(hyphenMatch+1,matcher.length()));

   TEUCHOS_TEST_FOR_EXCEPTION((matcherType != "coord") && (matcherType != "edge") && (matcherType != "face") && (matcherType != "all"),std::logic_error,
       "Failed parsing parameter list: could not find periodic boundary "
       "condition matcher \"" << matcher << "\" "
       "in string \"" << buildStr << "\n"
       "Type " << matcherType << " is not a valid boundary condition type. Must be coord, edge, face, or all\"");

   std::string matcherCoord = trim(matcher.substr(0,hyphenMatch));
   unsigned int matcherDim = 3;
   if((matcherCoord == "x") || (matcherCoord == "y") || (matcherCoord == "z"))
     matcherDim = 2;

   return std::make_pair(matcherType,matcherDim);
}

std::string PeriodicBC_Parser::replaceMatcherType(const std::string & buildStr, const std::string & matcherType) const
{
   std::string::size_type allPosition = buildStr.find("all");

   std::string beforeType = trim(buildStr.substr(0,allPosition));
   std::string afterType = trim(buildStr.substr(allPosition+3,buildStr.length()));

   return beforeType + matcherType + " " + afterType;
}

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
   std::string s_matcher, bndry1, bndry2;

   buildMatcher_Tokenize_withParams(buildStr,s_matcher,params,bndry1,bndry2);

   if(s_matcher=="x-coord") {
     panzer_stk::CoordMatcher matcher(0,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher);
   }

   if(s_matcher=="y-coord") {
     panzer_stk::CoordMatcher matcher(1,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher);
   }

   if(s_matcher=="z-coord") {
     panzer_stk::CoordMatcher matcher(2,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher);
   }

   if(s_matcher=="x-edge") {
     panzer_stk::CoordMatcher matcher(0,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher,"edge");
   }

   if(s_matcher=="y-edge") {
     panzer_stk::CoordMatcher matcher(1,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher,"edge");
   }

   if(s_matcher=="z-edge") {
     panzer_stk::CoordMatcher matcher(2,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher,"edge");
   }

   if(s_matcher=="xy-coord" || s_matcher=="yx-coord") {
     panzer_stk::PlaneMatcher matcher(0,1,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher);
   }

   if(s_matcher=="xz-coord" || s_matcher=="zx-coord") {
     panzer_stk::PlaneMatcher matcher(0,2,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher);
   }

   if(s_matcher=="yz-coord" || s_matcher=="zy-coord") {
     panzer_stk::PlaneMatcher matcher(1,2,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher);
   }

   if(s_matcher=="xy-edge" || s_matcher=="yx-edge") {
     panzer_stk::PlaneMatcher matcher(0,1,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher,"edge");
   }

   if(s_matcher=="xz-edge" || s_matcher=="zx-edge") {
     panzer_stk::PlaneMatcher matcher(0,2,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher,"edge");
   }

   if(s_matcher=="yz-edge" || s_matcher=="zy-edge") {
     panzer_stk::PlaneMatcher matcher(1,2,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher,"edge");
   }

   if(s_matcher=="xy-face" || s_matcher=="yx-face") {
     panzer_stk::PlaneMatcher matcher(0,1,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher,"face");
   }

   if(s_matcher=="xz-face" || s_matcher=="zx-face") {
     panzer_stk::PlaneMatcher matcher(0,2,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher,"face");
   }

   if(s_matcher=="yz-face" || s_matcher=="zy-face") {
     panzer_stk::PlaneMatcher matcher(1,2,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher,"face");
   }

   if(s_matcher=="wx-coord") {
     panzer_stk::WedgeMatcher matcher(WedgeMatcher::MirrorPlane::XZ_PLANE,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher);
   }

   if(s_matcher=="wy-coord") {
     panzer_stk::WedgeMatcher matcher(WedgeMatcher::MirrorPlane::YZ_PLANE,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher);
   }

   if(s_matcher=="wx-edge") {
     panzer_stk::WedgeMatcher matcher(WedgeMatcher::MirrorPlane::XZ_PLANE,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher,"edge");
   }

   if(s_matcher=="wy-edge") {
     panzer_stk::WedgeMatcher matcher(WedgeMatcher::MirrorPlane::YZ_PLANE,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher,"edge");
   }

   if(s_matcher=="wx-face") {
     panzer_stk::WedgeMatcher matcher(WedgeMatcher::MirrorPlane::XZ_PLANE,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher,"face");
   }

   if(s_matcher=="wy-face") {
     panzer_stk::WedgeMatcher matcher(WedgeMatcher::MirrorPlane::YZ_PLANE,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher,"face");
   }

   if(s_matcher=="(xy)z-quarter-coord") {
     panzer_stk::QuarterPlaneMatcher matcher(0,1,2,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher);
   }

   if(s_matcher=="(yx)z-quarter-coord") {
     panzer_stk::QuarterPlaneMatcher matcher(1,0,2,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher);
   }

   if(s_matcher=="(xz)y-quarter-coord") {
     panzer_stk::QuarterPlaneMatcher matcher(0,2,1,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher);
   }

   if(s_matcher=="(zx)y-quarter-coord") {
     panzer_stk::QuarterPlaneMatcher matcher(2,0,1,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher);
   }

   if(s_matcher=="(yz)x-quarter-coord") {
     panzer_stk::QuarterPlaneMatcher matcher(1,2,0,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher);
   }

   if(s_matcher=="(zy)x-quarter-coord") {
     panzer_stk::QuarterPlaneMatcher matcher(2,1,0,params);
     return panzer_stk::buildPeriodicBC_Matcher(bndry1,bndry2,matcher);
   }

   TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
       "Failed parsing parameter list: could not find periodic boundary "
       "condition matcher \"" << s_matcher << "\" "
       "in string \"" << buildStr << "\"");

   return Teuchos::null;
}

}
