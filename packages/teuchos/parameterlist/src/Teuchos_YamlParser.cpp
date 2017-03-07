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

#ifndef TEUCHOS_YAMLPARSER_DEF_H_
#define TEUCHOS_YAMLPARSER_DEF_H_

#include "Teuchos_YamlParser_decl.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_TwoDArray.hpp"

namespace Teuchos
{

template<typename T> Teuchos::Array<T> getYamlArray(const YAML::Node& node)
{
  Teuchos::Array<T> arr;
  for(YAML::const_iterator it = node.begin(); it != node.end(); it++)
  {
    arr.push_back(it->as<T>());
  }
  return arr;
}

void checkYamlTwoDArray(const YAML::Node& node, const std::string& key)
{
  for (YAML::const_iterator it = node.begin(); it != node.end(); ++it)
  {
    if (it->size() != node.begin()->size())
    {
      throw YamlSequenceError(std::string("TwoDArray \"") + key + "\" has irregular sizes");
    }
  }
}

template<typename T> Teuchos::TwoDArray<T> getYamlTwoDArray(const YAML::Node& node)
{
  Teuchos::TwoDArray<T> arr;
  typename Teuchos::TwoDArray<T>::size_type i, j;
  arr.resizeRows(node.size());
  arr.resizeCols(node.begin()->size());
  i = 0;
  for (YAML::const_iterator rit = node.begin(); rit != node.end(); ++rit)
  {
    j = 0;
    for (YAML::const_iterator cit = rit->begin(); cit != rit->end(); ++cit)
    {
      arr(i, j) = cit->as<T>();
      ++j;
    }
    ++i;
  }
  return arr;
}

/* Helper functions */

void updateParametersFromYamlFile(const std::string& yamlFileName,
                                  const Teuchos::Ptr<Teuchos::ParameterList>& paramList)
{
  //load the YAML file in as a new param list
  Teuchos::RCP<Teuchos::ParameterList> updated = YAMLParameterList::parseYamlFile(yamlFileName);
  //now update the original list (overwriting values with same key)
  paramList->setParameters(*updated);
}

void updateParametersFromYamlCString(const char* const data,
                                     const Teuchos::Ptr<Teuchos::ParameterList>& paramList,
                                     bool overwrite)
{
  Teuchos::RCP<Teuchos::ParameterList> updated = YAMLParameterList::parseYamlText(data);
  if(overwrite)
  {
    paramList->setParameters(*updated);
  }
  else
  {
    paramList->setParametersNotAlreadySet(*updated);
  }
}

void updateParametersFromYamlString(const std::string& yamlData,
                                  const Teuchos::Ptr<Teuchos::ParameterList>& paramList,
                                  bool overwrite)
{
  Teuchos::RCP<Teuchos::ParameterList> updated = YAMLParameterList::parseYamlText(yamlData);
  if(overwrite)
  {
    paramList->setParameters(*updated);
  }
  else
  {
    paramList->setParametersNotAlreadySet(*updated);
  }
}

Teuchos::RCP<Teuchos::ParameterList> getParametersFromYamlFile(const std::string& yamlFileName)
{
  return YAMLParameterList::parseYamlFile(yamlFileName);
}


std::string convertXmlToYaml(const std::string& xmlFileName)
{
  //load the parameter list from xml
  Teuchos::RCP<Teuchos::ParameterList> toConvert = Teuchos::getParametersFromXmlFile(xmlFileName);
  //replace the file extension ".xml" with ".yaml", or append it if there was no extension
  std::string yamlFileName;
  if(xmlFileName.find(".xml") == std::string::npos)
  {
    yamlFileName = xmlFileName + ".yaml";
  }
  else
  {
    yamlFileName = xmlFileName.substr(0, xmlFileName.length() - 4) + ".yaml";
  }
  YAMLParameterList::writeYamlFile(yamlFileName, toConvert);
  return yamlFileName;
}

void convertXmlToYaml(const std::string& xmlFileName, const std::string& yamlFileName)
{
  //load the parameter list from xml
  Teuchos::RCP<Teuchos::ParameterList> toConvert = Teuchos::getParametersFromXmlFile(xmlFileName);
  //replace the file extension ".xml" with ".yaml", or append it if there was no extension
  YAMLParameterList::writeYamlFile(yamlFileName, toConvert);
}

bool haveSameValuesUnordered(const Teuchos::ParameterList& lhs, const Teuchos::ParameterList& rhs, bool verbose)
{
  typedef Teuchos::ParameterList::ConstIterator Iter;
  Iter i = lhs.begin();
  Iter j = rhs.begin();
  if(lhs.name() != rhs.name())
  {
    if(verbose)
    {
      std::cout << "Parameter list names: \"" << lhs.name() << "\" and \"" << rhs.name() << "\".\n";
    }
    return false;
  }
  for(; i != lhs.end(); i++)
  {
    const std::string& key = lhs.name(i);
    const Teuchos::ParameterEntry& val1 = lhs.entry(i);
    //check that rhs also contains this key
    if(!rhs.isParameter(key))
    {
      if(verbose)
      {
        std::cout << "One list is missing parameter: \"" << key << "\"\n";
      }
      return false;
    }
    const Teuchos::ParameterEntry& val2 = rhs.getEntry(key);
    const Teuchos::any& any1 = val1.getAny(false);
    const Teuchos::any& any2 = val2.getAny(false);
    //check that types match
    if(any1.type() != any2.type())
    {
      if(verbose)
      {
        std::cout << "Values for key \"" << key << "\" have different types.\n";
      }
      return false;
    }
    //check for parameter list special case (don't use operator==)
    if(any1.type() == typeid(Teuchos::ParameterList))
    {
      if(!haveSameValuesUnordered(Teuchos::any_cast<Teuchos::ParameterList>(any1), Teuchos::any_cast<Teuchos::ParameterList>(any2), verbose))
      {
        //Don't need to print message here, the deepest list not matching will do that
        return false;
      }
    }
    else
    {
      //otherwise, use == to compare the values
      if(!(val1 == val2))
      {
        if(verbose)
        {
          std::cout << "Values for key \"" << key << "\" are different.\n";
        }
        return false;
      }
    }
    j++;
  }
  //lists must have same # of entries
  if(j != rhs.end())
  {
    if(verbose)
    {
      std::cout << "Lists \"" << lhs.name() << "\" and \"" << rhs.name() << "\" have different number of parameters.\n";
    }
    return false;
  }
  return true;
}

namespace YAMLParameterList
{

Teuchos::RCP<Teuchos::ParameterList> parseYamlText(const std::string& text)
{
  Teuchos::ParameterList pl;
  std::vector<YAML::Node> baseMap = YAML::LoadAll(text);
  return readParams(baseMap);
}

Teuchos::RCP<Teuchos::ParameterList> parseYamlText(const char* text)
{
  Teuchos::ParameterList pl;
  std::vector<YAML::Node> baseMap = YAML::LoadAll(text);
  return readParams(baseMap);
}

Teuchos::RCP<Teuchos::ParameterList> parseYamlFile(const std::string& yamlFile)
{
  std::vector<YAML::Node> baseMap = YAML::LoadAllFromFile(yamlFile);
  return readParams(baseMap);
}

Teuchos::RCP<Teuchos::ParameterList> readParams(std::vector<YAML::Node>& lists)
{
  Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList); //pl is the root ParameterList to be returned
  //If there is exactly one element in "lists", assume it is the anonymous top-level parameter list
  //If there are more than one, place them all in the anonymous top-level list
  for(size_t i = 0; i < lists.size(); i++)
  {
    processMapNode(lists[i], *pl, true);
  }
  return pl;
}

void processMapNode(const YAML::Node& node, Teuchos::ParameterList& parent, bool topLevel)
{
  if(node.Type() != YAML::NodeType::Map)
  {
    throw YamlStructureError("All top-level elements of the YAML file must be maps.");
  }
  if(topLevel)
  {
    parent.setName("ANONYMOUS");
    processMapNode(node.begin()->second, parent);
  }
  else
  {
    for(YAML::const_iterator i = node.begin(); i != node.end(); i++)
    {
      //make sure the key type is a string
      if(i->first.Type() != YAML::NodeType::Scalar)
      {
        throw YamlKeyError("Keys must be plain strings");
      }
      //if this conversion fails and throws for any reason (shouldn't), let the caller handle it
      const std::string key = i->first.as<std::string>();
      processKeyValueNode(key, i->second, parent, topLevel);
    }
  }
}

void processKeyValueNode(const std::string& key, const YAML::Node& node, Teuchos::ParameterList& parent, bool topLevel)
{
  //node (value) type can be a map (for nested param lists),
  //a scalar (int, double, string), or a sequence of doubles (vector<double>)
  if(node.Type() == YAML::NodeType::Scalar)
  {
    try
    {
      parent.set(key, node.as<int>());
    }
    catch(...)
    {
      try
      {
        parent.set(key, node.as<double>());
      }
      catch(...)
      {
        try
        {
          std::string rawString = node.as<std::string>();
          if(rawString == "true")
          {
            parent.set<bool>(key, true);
          }
          else if(rawString == "false")
          {
            parent.set<bool>(key, false);
          }
          else
          {
            parent.set(key, rawString);
          }
        }
        catch(...)
        {
          throw YamlScalarError("YAML scalars must be int, double, bool or string.");
        }
      }
    }
  }
  else if(node.Type() == YAML::NodeType::Map)
  {
    if(topLevel)
    {
      processMapNode(node, parent);
    }
    else
    {
      Teuchos::ParameterList& sublist = parent.sublist(key);
      processMapNode(node, sublist);
    }
  }
  else if(node.Type() == YAML::NodeType::Sequence)
  {
    if (node.begin()->Type() == YAML::NodeType::Sequence) {
      checkYamlTwoDArray(node, key);
      try
      {
        node.begin()->begin()->as<int>();
        parent.set(key, getYamlTwoDArray<int>(node));
      }
      catch(...)
      {
        try
        {
          node.begin()->begin()->as<double>();
          parent.set(key, getYamlTwoDArray<double>(node));
        }
        catch(...)
        {
          try
          {
            node.begin()->begin()->as<std::string>();
            parent.set(key, getYamlTwoDArray<std::string>(node));
          }
          catch(...)
          {
            throw YamlSequenceError(std::string("TwoDArray \"") + key + "\" must contain int, double, bool or string");
          }
        }
      }
    } else {
      //typeString is used to provide a useful error message if types inconsistent
      try
      {
        node.begin()->as<int>();
        parent.set(key, getYamlArray<int>(node));
      }
      catch(...)
      {
        try
        {
          node.begin()->as<double>();
          parent.set(key, getYamlArray<double>(node));
        }
        catch(...)
        {
          try
          {
            node.begin()->as<std::string>();
            parent.set(key, getYamlArray<std::string>(node));
          }
          catch(...)
          {
            throw YamlSequenceError(std::string("Array \"") + key + "\" must contain int, double, bool or string");
          }
        }
      }
    }
  }
  else if(node.Type() == YAML::NodeType::Null)
  {
    //treat NULL as empty string (not an error)
    parent.set(key, std::string());
  }
  else
  {
    //Undefined
    throw YamlUndefinedNodeError("Value type in a key-value pair must be one of: int, double, string, array, sublist.");
  }
}



void writeYamlFile(const std::string& yamlFile, Teuchos::RCP<Teuchos::ParameterList>& pl)
{
  std::ofstream yaml(yamlFile);
  yaml << "%YAML 1.1\n---\n";
  yaml << "ANONYMOUS:";         //original top-level list name is not stored by ParameterList
  if(pl->numParams() == 0)
  {
    yaml << " { }\n";
  }
  else
  {
    writeParameterList(*pl, yaml, 2);
  }
  yaml << "...";
}

void writeParameterList(Teuchos::ParameterList& pl, std::ofstream& yaml, int indentLevel)
{
  if(pl.begin() == pl.end())
  {
    yaml << "{ }\n";
  }
  else
  {
    yaml << '\n';
    for(PLIter it = pl.begin(); it != pl.end(); it++)
    {
      writeParameter(pl.name(it), pl.entry(it), yaml, indentLevel);
    }
  }
}

void writeParameter(const std::string& paramName, const Teuchos::ParameterEntry& entry, std::ofstream& yaml, int indentLevel)
{
  for(int i = 0; i < indentLevel; i++)
  {
    yaml << ' ';
  }
  generalWriteString(paramName, yaml);
  yaml << ": ";
  if(entry.isList())
  {
    writeParameterList(Teuchos::getValue<Teuchos::ParameterList>(entry), yaml, indentLevel + 2);
    return;
  }
  else if(entry.isArray())
  {
    yaml << '[';
    if(entry.isType<Teuchos::Array<int> >())
    {
      Teuchos::Array<int>& arr = Teuchos::getValue<Teuchos::Array<int> >(entry);
      for(int i = 0; i < arr.size(); i++)
      {
        yaml << arr[i];
        if(i != arr.size() - 1)
          yaml << ", ";
      }
    }
    else if(entry.isType<Teuchos::Array<double> >())
    {
      Teuchos::Array<double>& arr = Teuchos::getValue<Teuchos::Array<double> >(entry);
      for(int i = 0; i < arr.size(); i++)
      {
        generalWriteDouble(arr[i], yaml);
        if(i != arr.size() - 1)
          yaml << ", ";
      }
    }
    else if(entry.isType<Teuchos::Array<std::string> >())
    {
      Teuchos::Array<std::string>& arr = Teuchos::getValue<Teuchos::Array<std::string> >(entry);
      for(int i = 0; i < arr.size(); i++)
      {
        generalWriteString(arr[i], yaml);
        if(i != arr.size() - 1)
          yaml << ", ";
      }
    }
    yaml << ']';
  }
  else if(entry.isType<int>())
  {
    yaml << Teuchos::getValue<int>(entry);
  }
  else if(entry.isType<double>())
  {
    generalWriteDouble(Teuchos::getValue<double>(entry), yaml);
  }
  else if(entry.isType<std::string>())
  {
    std::string& str = Teuchos::getValue<std::string>(entry);
    if(strchr(str.c_str(), '\n'))
    {
      //need explicit indentation so that indentation in the string is preserved
      yaml << "|2-\n";    
      //for each line, apply indent then print the line verbatim
      size_t index = 0;
      while(true)
      {
        size_t next = str.find('\n', index);
        for(int i = 0; i < indentLevel + 2; i++)
        {
          yaml << ' ';
        }
        if(next == std::string::npos)
        {
          yaml << str.substr(index, std::string::npos);
          break;
        }
        else
        {
          yaml << str.substr(index, next - index) << '\n';
        }
        index = next + 1;
      }
    }
    else
    {
      generalWriteString(str, yaml);
    }
  }
  else if(entry.isType<bool>())
  {
    yaml << (Teuchos::getValue<bool>(entry) ? "true" : "false");
  }
  yaml << '\n';
}

void generalWriteString(const std::string& str, std::ofstream& yaml)
{
  if(stringNeedsQuotes(str))
  {
    yaml << '\'' << str << '\'';
  }
  else
  {
    yaml << str;
  }
}

void generalWriteDouble(double d, std::ofstream& yaml)
{
  yaml << std::showpoint << std::setprecision(8);
  if(d < 1e6 && d > 1e-5)
  {
    //use regular notation
    yaml << d;
  }
  else
  {
    yaml << std::scientific << d;
  }
  yaml << std::fixed;
}

bool stringNeedsQuotes(const std::string& str)
{
  return strpbrk(str.c_str(), ":{}[],&*#?|-<>=!%@\\");
}

} //namespace YAMLParameterList

} //namespace Teuchos

#endif
