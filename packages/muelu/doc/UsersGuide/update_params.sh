#!/bin/sh

xsltproc tex.xsl masterList.xml > paramlist.tex
xsltproc tex_hidden.xsl masterList.xml > paramlist_hidden.tex


code_file="../../src/MueCentral/MueLu_MasterList.cpp"

echo '// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include "MueLu_Exceptions.hpp"
#include "MueLu_MasterList.hpp"

namespace MueLu {

  Teuchos::RCP<const Teuchos::ParameterList> MasterList::List() {
    if (masterList_.is_null()) {
      masterList_ = Teuchos::getParametersFromXmlString(stringList_);
    }

    return masterList_;
  }

  Teuchos::RCP<Teuchos::ParameterList> MasterList::GetProblemSpecificList(std::string const & problemType) {

    if ( (problemType != problemType_) || problemSpecificList_.is_null() ) {
      if (DefaultProblemTypeLists_.find(problemType) != DefaultProblemTypeLists_.end()) {
        problemType_ = problemType;
        problemSpecificList_ = Teuchos::getParametersFromXmlString(DefaultProblemTypeLists_[problemType]);
      } else {
        //TODO provide valid problem types
        TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "Invalid problem type " << problemType << ".");
      }
    }
    return problemSpecificList_;
  }
 
   std::string MasterList::interpretParameterName(const std::string& name, const std::string& value) {

    // used to concatenate the return string
    std::stringstream ss;

    // put in short cuts here!

    if (name == "verbosity") { 
      std::string verb = "none";
      if (value == "\"0\"") verb = "none";
      if (value == "\"1\"" || value == "\"2\"" || value == "\"3\"") verb = "low";
      if (value == "\"4\"" || value == "\"5\"" || value == "\"6\"") verb = "medium";
      if (value == "\"7\"" || value == "\"8\"") verb = "high";
      if (value == "\"9\"") verb = "extreme";
      if (value == "\"10\"") verb = "test";  
      verb = "\"" + verb + "\"";
      ss << "<Parameter name=\"verbosity\" type=\"string\" value=" << verb << "/>"; 
      return ss.str(); 
    }
    
    if (name == "cycle type") {
      std::stringstream temp1; temp1 << "\"" << "MGV" << "\"";
      std::stringstream temp2; temp2 << "\"" << "MGV" << "\"";
      if (value == temp1.str() ) { ss << "<Parameter name=\"cycle type\" type=\"string\" value=\"V\"/>"; return ss.str(); }
      else if (value == temp2.str()) { ss << "<Parameter name=\"cycle type\" type=\"string\" value=\"W\"/>"; return ss.str(); }
      else TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "MasterList::interpretParameterName, Line " << __LINE__ << ". "
                                           << "The parameter " << value << " is not supported by MueLu.");
      return ss.str();
    }    

    // energy minimization is enabled
    if (name == "multigrid algorithm") {
      std::stringstream temp; temp << "\"" << "1" << "\"";
      if (value == temp.str() ) { ss << "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"pg\"/>"; return ss.str(); }
    }

    if (name == "repartition: enable") {
      std::stringstream temp1; temp1 << "\"" << "1" << "\"";
      if (value == temp1.str()) {
        RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
        *out << "WARNING: repartitioning in MueLu is different to ML's. Please refer to the MueLu user's Manual for more information." << std::endl;
      }
    }
    
    // put in auto-generated code here
' > $code_file
xsltproc gen_interpreter.xsl masterList.xml >> $code_file

echo '
    return "";
  }
 
  Teuchos::RCP<Teuchos::ParameterList> MasterList::masterList_ = Teuchos::null;
  Teuchos::RCP<Teuchos::ParameterList> MasterList::problemSpecificList_ = Teuchos::null;
  std::string                          MasterList::problemType_ = "unknown";
  const std::string                    MasterList::stringList_ =' >> $code_file

xsltproc paramlist.xsl masterList.xml >> $code_file

echo ';' >> $code_file

echo '  std::map<std::string,std::string> MasterList::DefaultProblemTypeLists_ = DefaultProblemStrings<std::string,std::string>' >> $code_file

PROBLEM_TYPES=( "Poisson-2D" "Poisson-3D" "Elasticity-2D" "Elasticity-3D" "MHD" "ConvectionDiffusion" )

for i in "${PROBLEM_TYPES[@]}"; do
  echo "(\"$i\"," >> $code_file
  xsltproc --stringparam prob_type "$i" probtypelist.xsl masterList.xml >> $code_file
  echo ')' >> $code_file
done

echo ";" >> $code_file

echo '  std::map<std::string,std::string> MasterList::ML2MueLuLists_ = DefaultProblemStrings<std::string,std::string>' >> $code_file
xsltproc ml2muelu.xsl masterList.xml >> $code_file
echo ';

}
' >> $code_file

# fix quotation
sed -i '/<Parameter/ s/\\""/\\"/g' $code_file
sed -i '/<Parameter/ s/"\\"/\\"/g' $code_file

# generate LaTeX files (MueLu options and ML compatibility)
SECTIONS=( "general" "smoothing_and_coarse" "aggregation" "misc" "multigrid" "rebalancing" "reuse" )
for i in "${SECTIONS[@]}"; do
  xsltproc --stringparam section "$i" options.xsl   masterList.xml > options_$i.tex
  xsltproc --stringparam section "$i" mloptions.xsl masterList.xml > mloptions_$i.tex
done
