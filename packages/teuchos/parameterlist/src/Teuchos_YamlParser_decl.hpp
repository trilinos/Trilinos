// @HEADER
//
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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

#ifndef TEUCHOS_YAMLPARSER_DECL_H_
#define TEUCHOS_YAMLPARSER_DECL_H_

#include "yaml-cpp/yaml.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_PtrDecl.hpp"
#include "Teuchos_FileInputSource.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstring>

namespace Teuchos
{

#define MAKE_EXCEPTION_TYPE(Name) \
class Name : public Teuchos::ExceptionBase \
{ \
  public: \
    Name(const std::string& arg) : ExceptionBase(arg) {} \
};

MAKE_EXCEPTION_TYPE(YamlKeyError)
MAKE_EXCEPTION_TYPE(YamlScalarError)
MAKE_EXCEPTION_TYPE(YamlSequenceError)
MAKE_EXCEPTION_TYPE(YamlStructureError)
MAKE_EXCEPTION_TYPE(YamlUndefinedNodeError)

#undef MAKE_EXCEPTION_TYPE

std::string convertXmlToYaml(const std::string& xmlFileName); //returns filename of produced YAML file
void convertXmlToYaml(const std::string& xmlFileName, const std::string& yamlFileName); //writes to given filename
void convertXmlToYaml(std::istream& xmlStream, std::ostream& yamlStream);
bool haveSameValuesUnordered(const Teuchos::ParameterList& lhs, const Teuchos::ParameterList& rhs, bool verbose = false);

//Class modeled after Teuchos::XMLParameterListReader
namespace YAMLParameterList
{
  Teuchos::RCP<Teuchos::ParameterList> parseYamlText(const std::string& text);
  Teuchos::RCP<Teuchos::ParameterList> parseYamlText(const char* text);
  Teuchos::RCP<Teuchos::ParameterList> parseYamlFile(const std::string& yamlFile);
  Teuchos::RCP<Teuchos::ParameterList> parseYamlStream(std::istream& yaml);
  void writeYamlStream(std::ostream& yamlFile, const Teuchos::ParameterList& pl);
  void writeYamlFile(const std::string& yamlFile, const Teuchos::ParameterList& pl);
  Teuchos::RCP<Teuchos::ParameterList> readParams(std::vector<YAML::Node>& lists);
  //load all k-v pairs within node into param list (checks if node is map, and handles nesting)
  //topLevel means to put sub-pairs directly into parent and not create named sublists
  void processMapNode(const YAML::Node& node, Teuchos::ParameterList& parent, bool topLevel = false);
  void processKeyValueNode(const std::string& key, const YAML::Node& node, Teuchos::ParameterList& parent, bool topLevel = false);
  //  template<typename T> Teuchos::Array<T> getYamlArray(const YAML::Node& node);
  void writeParameterList(const Teuchos::ParameterList& pl, std::ostream& yaml, int indentLevel);
  void writeParameter(const std::string& paramName, const Teuchos::ParameterEntry& entry, std::ostream& yaml, int indentLevel);    //throws if the entry's type is not supported
  void generalWriteString(const std::string& str, std::ostream& yaml);
  void generalWriteDouble(double d, std::ostream& yaml);
  bool stringNeedsQuotes(const std::string& str);
  typedef Teuchos::ParameterList::ConstIterator PLIter;
}

} //namespace

#endif
