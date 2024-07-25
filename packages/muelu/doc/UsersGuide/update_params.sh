#!/bin/bash

xsltproc tex.xsl masterList.xml > paramlist.tex
returncode=$?
if [ $returncode -ne 0 ]; then
    # There seems to be an issue with certain versions of xsltproc
    # suddenly segfaulting.
    echo "xsltproc exited unexpectedly. Please investigate!"
    exit returncode;
fi

xsltproc tex_hidden.xsl masterList.xml > paramlist_hidden.tex

if [ "$1" != "" ]; then
    code_file=$1
else
    code_file="../../src/MueCentral/MueLu_MasterList.cpp"
fi


echo '// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// clang-format off
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
      if (value == temp1.str() ) { ss << "<Parameter name=\"cycle type\" type=\"string\" value=\"V\"/>"; }
      else if (value == temp2.str()) { ss << "<Parameter name=\"cycle type\" type=\"string\" value=\"W\"/>"; }
      else TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "MasterList::interpretParameterName, Line " << __LINE__ << ". "
                                           << "The parameter " << value << " is not supported by MueLu.");
      return ss.str();
    }

    // energy minimization is enabled
    if (name == "multigrid algorithm") {
      std::stringstream temp; temp << "\"" << "1" << "\"";
      if (value == temp.str() ) { ss << "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"pg\"/>"; return ss.str(); }
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

PROBLEM_TYPES=( "Poisson-2D" "Poisson-2D-complex" "Poisson-3D" "Poisson-3D-complex" "Elasticity-2D" "Elasticity-2D-complex" "Elasticity-3D" "Elasticity-3D-complex" "MHD" "ConvectionDiffusion" )

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

# fix quotation using sed
# GH: similar to other instances in MueLu, we need to work around the GNU/BSD sed issue
#     the moral of the story is -i is not portable for testing!
sed '/<Parameter/ s/\\""/\\"/g' "$code_file" > "$code_file.tmp"
mv "$code_file.tmp" "$code_file"
sed '/<Parameter/ s/"\\"/\\"/g' "$code_file" > "$code_file.tmp"
mv "$code_file.tmp" "$code_file"

# generate LaTeX files (MueLu options and ML compatibility)
SECTIONS=( "general" "smoothing_and_coarse" "aggregation" "misc" "multigrid" "rebalancing" "reuse" "refmaxwell" )
for i in "${SECTIONS[@]}"; do
  xsltproc --stringparam section "$i" options.xsl   masterList.xml > options_$i.tex
  xsltproc --stringparam section "$i" mloptions.xsl masterList.xml > mloptions_$i.tex
done
