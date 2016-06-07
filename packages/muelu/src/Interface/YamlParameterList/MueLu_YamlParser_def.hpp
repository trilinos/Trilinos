// @HEADER
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

#ifndef MUELU_YAMLPARSER_DEF_H_
#define MUELU_YAMLPARSER_DEF_H_

#include "MueLu_YamlParser_decl.hpp"

namespace MueLu
{

/* Helper functions */

void updateParametersFromYamlFile(const std::string& yamlFileName,
                                  const Teuchos::Ptr<Teuchos::ParameterList>& paramList)
{
  //load the YAML file in as a new param list
  MueLu::YAMLParameterListReader reader;
  Teuchos::RCP<Teuchos::ParameterList> updated = reader.parseYamlFile(yamlFileName);
  //now update the original list (overwriting values with same key)
  paramList->setParameters(*updated);
}

void updateParametersFromYamlCString(const char* const data,
                                     const Teuchos::Ptr<Teuchos::ParameterList>& paramList,
                                     bool overwrite)
{
  MueLu::YAMLParameterListReader reader;
  Teuchos::RCP<Teuchos::ParameterList> updated = reader.parseYamlText(data);
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
  MueLu::YAMLParameterListReader reader;
  Teuchos::RCP<Teuchos::ParameterList> updated = reader.parseYamlText(yamlData);
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
  MueLu::YAMLParameterListReader reader;
  return reader.parseYamlFile(yamlFileName);
}

void updateParametersFromYamlFileAndBroadcast(const std::string &yamlFileName, const Teuchos::Ptr<Teuchos::ParameterList> &paramList, const Teuchos::Comm<int> &comm, bool overwrite)
{
  struct SafeFile
  {
    SafeFile(const char* fname, const char* options)
    {
      handle = fopen(fname, options);
    }
    ~SafeFile()
    {
      if(handle)
        fclose(handle);
    }
    FILE* handle;
  };
  //BMK note: see teuchos/comm/src/Teuchos_XMLParameterListHelpers.cpp
  if(comm.getSize() == 1)
  {
    updateParametersFromYamlFile(yamlFileName, paramList);
  }
  else
  {
    if(comm.getRank() == 0)
    {
      //BMK: TODO! //reader.setAllowsDuplicateSublists(false);
      //create a string and load file contents into it
      //C way for readability and speed, same thing with C++ streams is slow & ugly
      SafeFile yamlFile(yamlFileName.c_str(), "rb");
      if(!yamlFile.handle)
      {
        throw std::runtime_error(std::string("Failed to open YAML file \"") + yamlFileName + "\"for reading.");
      }
      fseek(yamlFile.handle, 0, SEEK_END);
      int strsize = ftell(yamlFile.handle) + 1;
      rewind(yamlFile.handle);
      //Make the array raii
      Teuchos::ArrayRCP<char> contents(new char[strsize], 0, strsize, true);
      fread((void*) contents.get(), strsize - 1, 1, yamlFile.handle);
      contents.get()[strsize - 1] = 0;
      Teuchos::broadcast<int, int>(comm, 0, &strsize);
      Teuchos::broadcast<int, char>(comm, 0, strsize, contents.get());
      updateParametersFromYamlCString(contents.get(), paramList, overwrite);
    }
    else
    {
      int strsize;
      Teuchos::broadcast<int, int>(comm, 0, &strsize);
      Teuchos::ArrayRCP<char> contents(new char[strsize], 0, strsize, true);
      Teuchos::broadcast<int, char>(comm, 0, strsize, contents.get());
      updateParametersFromYamlCString(contents.get(), paramList, overwrite);
    }
  }
}

/* YAMLParameterListReader functions */

Teuchos::RCP<Teuchos::ParameterList> YAMLParameterListReader::parseYamlText(const std::string& text)
{
  Teuchos::ParameterList pl;
  std::vector<YAML::Node> baseMap = YAML::LoadAll(text);
  return readParams(baseMap);
}

Teuchos::RCP<Teuchos::ParameterList> YAMLParameterListReader::parseYamlText(const char* text)
{
  Teuchos::ParameterList pl;
  std::vector<YAML::Node> baseMap = YAML::LoadAll(text);
  return readParams(baseMap);
}

Teuchos::RCP<Teuchos::ParameterList> YAMLParameterListReader::parseYamlFile(const std::string& yamlFile)
{
  std::vector<YAML::Node> baseMap = YAML::LoadAllFromFile(yamlFile);
  return readParams(baseMap);
}

Teuchos::RCP<Teuchos::ParameterList> YAMLParameterListReader::readParams(std::vector<YAML::Node>& lists)
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

void YAMLParameterListReader::processMapNode(const YAML::Node& node, Teuchos::ParameterList& parent, bool topLevel)
{
  if(node.Type() != YAML::NodeType::Map)
  {
    throw YamlStructureError("All top-level elements of the YAML file must be maps.");
  }
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

void YAMLParameterListReader::processKeyValueNode(const std::string& key, const YAML::Node& node, Teuchos::ParameterList& parent, bool topLevel)
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
          parent.set(key, node.as<bool>());
        }
        catch(...)
        {
          try
          {
            parent.set(key, node.as<std::string>());
          }
          catch(...)
          {
            throw YamlScalarError("YAML scalars must be int, double, bool or string.");
          }
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

template<typename T> Teuchos::Array<T> YAMLParameterListReader::getYamlArray(const YAML::Node& node)
{
  Teuchos::Array<T> arr;
  for(YAML::const_iterator it = node.begin(); it != node.end(); it++)
  {
    arr.push_back(it->as<T>());
  }
  return arr;
}

} //namespace MueLu

#endif
