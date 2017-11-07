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

/*! \file Teuchos_YamlParser_decl.hpp
    \brief Functions to convert between ParameterList and YAML

YAML is a human-readable data serialization format. Teuchos provides a
YAML parameter list interpreter. It produces Teuchos::ParameterList
objects equivalent to those produced by the Teuchos XML helper functions.

Here is a simple example XML parameter list:
\code{.xml}
<ParameterList>
  <ParameterList Input>
    <Parameter name="values" type="Array(double)" value="{54.3 -4.5 2.0}"/>
    <Parameter name="myfunc" type="string" value="
def func(a, b):
  return a * 2 - b"/>
  </ParameterList>
  <ParameterList Solver>
    <Parameter name="iterations" type="int" value="5"/>
    <Parameter name="tolerance" type="double" value="1e-7"/>
    <Parameter name="do output" type="bool" value="true"/>
    <Parameter name="output file" type="string" value="output.txt"/>
  </ParameterList>
</ParameterList>
\endcode

Here is an equivalent YAML parameter list:
\code{.yaml}
%YAML 1.1
---
ANONYMOUS:
  Input:
    values: [54.3, -4.5, 2.0]
    myfunc: |-

      def func(a, b):
        return a * 2 - b
  Solver:
    iterations: 5
    tolerance: 1e-7
    do output: yes
    output file: output.txt
...
\endcode

The nested structure and key-value pairs of these two lists are identical.
To a program querying them for settings, they are indistinguishable.

These are the general rules for creating a YAML parameter list:
<ol>
<li> First line must be <code>&#37YAML 1.1</code>, second must be <code>---</code>, and last must be <code>...</code> </li>
<li> Nested map structure is determined by indentation. SPACES ONLY, NO TABS! </li>
<li> As with the above example, for a top-level anonymous parameter list, <code>ANONYMOUS:</code> must be explicit </li>
<li> Type is inferred. <code>5</code> is an int, <code>5.0</code> is a double, and <code>'5.0'</code> is a string </li>
<li> Quotation marks (single or double) are optional for strings, but required for
     strings with special characters: <code>:{}[],&*#?|<>=!%@\</code>. </li>
<li> Quotation marks also turn non-string types into strings: <code>'3'</code> is a string </li>
<li> As with XML parameter lists, keys are regular strings </li>
<li> Even though YAML supports several names for bool true/false, only <code>true</code> and
     <code>false</code> are supported by the parameter list reader. </li>
<li> Arrays of int, double and string supported. <code>exampleArray: [hello, world, goodbye]</code> </li>
<li> <code>[3, 4, 5]</code> is an int array, <code>[3, 4, 5.0]</code> is a double array,
     and <code>[3, '4', 5.0]</code> is a string array </li>
<li> For multi-line strings, place <code>|</code> after the <code>key:</code> and
     then indent the string one level deeper than the key </li>
<li> To preserve indentation in a multiline string, place <code>|2</code> after the <code>key:</code> and
     then indent your string's content by 2 spaces relative to the key. </li>
</ol>

Note that the <code>|-</code> and empty line after <code>myfunc:</code> in this
example remove the trailing newline from the <code>myfunc</code> content and add
a newline to the start of it, respectively.
This is only done to mimic the fact that the XML value has a newline at
the beginning and not the end due to the way it was placed for readability in the
XML file. However, a more idiomatic way to include a multi-line string in YAML is:

\code{.yaml}
    myfunc: |
      def func(a, b):
        return a * 2 - b
\endcode

This will create a string with a newline at the end and not the beginning, which
is more natural, will print better to other file formats, looks more like the
contents of its own file from the perspective of the parser handling it (e.g. the RTCompiler)
and looks better in YAML itself.

*/

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_PtrDecl.hpp"
#include "Teuchos_FileInputSource.hpp"

#include <iostream>
#include <string>

namespace Teuchos
{

std::string convertXmlToYaml(const std::string& xmlFileName); //returns filename of produced YAML file
void convertXmlToYaml(const std::string& xmlFileName, const std::string& yamlFileName); //writes to given filename
void convertXmlToYaml(std::istream& xmlStream, std::ostream& yamlStream);

//Class modeled after Teuchos::XMLParameterListReader
namespace YAMLParameterList
{
  Teuchos::RCP<Teuchos::ParameterList> parseYamlText(const std::string& text,
      const std::string& name);
  Teuchos::RCP<Teuchos::ParameterList> parseYamlFile(const std::string& yamlFile);
  Teuchos::RCP<Teuchos::ParameterList> parseYamlStream(std::istream& yaml);
  void writeYamlStream(std::ostream& yamlFile, const Teuchos::ParameterList& pl);
  void writeYamlFile(const std::string& yamlFile, const Teuchos::ParameterList& pl);
  void writeParameterList(const Teuchos::ParameterList& pl, std::ostream& yaml, int indentLevel);
  void writeParameter(const std::string& paramName, const Teuchos::ParameterEntry& entry, std::ostream& yaml, int indentLevel);    //throws if the entry's type is not supported
  void generalWriteString(const std::string& str, std::ostream& yaml);
  void generalWriteDouble(double d, std::ostream& yaml);
  bool stringNeedsQuotes(const std::string& str);
  typedef Teuchos::ParameterList::ConstIterator PLIter;
}

} //namespace

#endif
