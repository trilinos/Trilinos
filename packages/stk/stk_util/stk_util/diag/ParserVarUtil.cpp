// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include "stk_util/diag/ParserVarUtil.hpp"
#include "stk_util/diag/StringUtil.hpp"
#include "stk_util/util/ReportHandler.hpp"

namespace stk {
namespace diag {


  bool remove_entry(const std::vector<sierra::String> &valid,  const sierra::String& input) 
  {
    for(const auto& vf : valid){
      if(vf == input){ 
        return true;
      }
    }
    return false;
  }

  std::string strip_entry(const std::vector<sierra::String> &valid,  const sierra::String& input) 
  {
    std::string val = input.s_str(); 
    for(auto& toRemove : valid) {
      if (val.rfind(toRemove, 0) == 0) {
        // Erase the substring from the original string
        val.erase(0, toRemove.s_str().length());
        return val;
      }
    }
    return val;
  }
  
  std::vector<sierra::String> remove_var_eq(const std::vector<sierra::String>& input) 
  {
    static const std::vector<sierra::String> validFirst{"=", "is", "are", "variable", "variable=","variables", "variables="};
    static const std::vector<sierra::String> validSecond{"=", "is", "are"};
    static const std::vector<sierra::String> stripEntries{"=", "variable=", "variables="};
    std::vector<sierra::String> output; 
    bool firstEntryRemoved = false;
    for(int i = 0; i < (int)input.size(); i++ ){
      if(i == 0) {
        firstEntryRemoved = remove_entry(validFirst,input[i]);
        if(!firstEntryRemoved){
          std::string variable = strip_entry(stripEntries, input[i]);
          output.emplace_back(variable);
        }
      } else if(i == 1 && firstEntryRemoved) {
        if(!remove_entry(validSecond,input[i])){
          std::string variable = strip_entry(stripEntries, input[i]);
          output.emplace_back(variable);
        }
      } else{
        output.emplace_back(input[i]);
      }
    }
    return output;
  }
  
  std::vector<sierra::String> join_paren_vars(const std::vector<sierra::String>& input, sierra::String tag) 
  {
    std::vector<sierra::String> varList;
    for(int i=0; i<(int)input.size(); i++) {
      sierra::String val = input[i];
      if(i+1<(int)input.size() && input[i+1].s_str().front() == '('){
        val = val + input[++i];
      }
      if(val.s_str().find('(') != std::string::npos) {
        while(val.s_str().back() != ')'){
          if(i+1 >= (int)input.size()){
            varList.push_back(val);
            return varList;
          }
          if(input[i+1].s_str().front() == ')' || val.s_str().back() == '('){
            val = val + input[++i];
          } else {
            val = val +  ","+ input[++i];
          }
        }
      }
      ParsedVariable var(val, tag);
      varList.push_back(var.name());
    }
    return varList;
  }


  
  //
  //  Convert range strings into enum values.  Note, cannot pick an exact index as the
  //  index might depend on the type.  (zy in a symmetric tensor is not the same as zy in a full tensor.)
  //
  //  Return:  The variable component.  Positive values are regular integer components, negative values are enums with
  //  special meanings
  //
  VariableComponent::VariableComponent(const std::string &rangeString, const sierra::String &varName)
                                      : m_savedRangeString(rangeString) {
    sierra::String range_lower = sierra::make_lower(rangeString);
    if (range_lower == "x") {
      m_componentEnum = VECTOR_X_COMPONENT;
    } else if (range_lower == "y") {
      m_componentEnum = VECTOR_Y_COMPONENT;
    } else if (range_lower == "z") {
      m_componentEnum = VECTOR_Z_COMPONENT;
    } else if (range_lower == "*" || range_lower == ":") {
      m_componentEnum = ALL_COMPONENTS;
    } else if (range_lower == "") {
      m_componentEnum = UNSPECIFIED_COMPONENT;
    } else if (range_lower == "xx") {
      m_componentEnum = TENSOR_XX_COMPONENT;
    } else if (range_lower == "yy") {
      m_componentEnum = TENSOR_YY_COMPONENT;
    } else if (range_lower == "zz") {
      m_componentEnum = TENSOR_ZZ_COMPONENT;
    } else if (range_lower == "xy") {
      m_componentEnum = TENSOR_XY_COMPONENT;
    } else if (range_lower == "xz") {
      m_componentEnum = TENSOR_XZ_COMPONENT;
    } else if (range_lower == "yz") {
      m_componentEnum = TENSOR_YZ_COMPONENT;
    } else if (range_lower == "yx") {
       m_componentEnum = TENSOR_YX_COMPONENT;
    } else if (range_lower == "zx") {
      m_componentEnum = TENSOR_ZX_COMPONENT;
    } else if (range_lower == "zy") {
      m_componentEnum = TENSOR_ZY_COMPONENT;
    } else {
      if (std::sscanf(range_lower.c_str(), "%d", &m_integerComponent) > 0 && std::stoi(range_lower.c_str()) > 0) {
        m_componentEnum = INTEGER_COMPONENT;
      }
    }
  }

  ParsedVariable::ParsedVariable(const sierra::String &var_name, sierra::String tag) {
    std::smatch matches;
    static const std::string beginWord = R"(^([A-Za-z0-9_\.\-\%\:\>\+]+)?\s*)";
    static const std::string openParen = R"(\()";
    static const std::string comma = R"(,)";
    static const std::string closeParen = R"(\)\s*$)";
    static const std::string componentOrIndex = R"(\s*(:|\*|x|y|z|xx|yy|zz|xy|xz|yz|yx|zx|zy|[0-9]+)\s*)";
    static const std::string indexOrAll = R"(\s*(:|\*|[0-9]+)\s*)";
    static const std::regex baseWithOneComponent{beginWord + openParen + componentOrIndex + closeParen};
    static const std::regex baseWithTwoComponentComma{beginWord + openParen + componentOrIndex + comma + indexOrAll + closeParen};
    static const std::regex baseWithTwoComponent{beginWord + openParen + componentOrIndex + indexOrAll + closeParen};
    static const std::regex baseNameOnly{beginWord + std::string(R"(\s*$)")};
    std::string lowerStr = sierra::make_lower(var_name).s_str();
    if (std::regex_search(lowerStr, matches, baseWithOneComponent)) {
      std::string tempName = var_name.s_str();
      tempName.resize(matches[1].str().size());
      tempName.erase(std::remove(tempName.begin(), tempName.end(), ' '), tempName.end());
      baseName = tempName;
      comp1 = VariableComponent(matches[2].str(), matches[1].str());
      comp2 = VariableComponent("", matches[1].str());
      explictNumComponents = 1;
    } else if (std::regex_search(lowerStr, matches, baseNameOnly)) {
      std::string tempName = var_name.s_str();
      tempName.erase(std::remove(tempName.begin(), tempName.end(), ' '), tempName.end());
      baseName = tempName;
      comp1 = VariableComponent(":", matches[1].str());
      comp2 = VariableComponent(":", matches[1].str());
      explictNumComponents = 0;
    } else if (std::regex_search(lowerStr, matches, baseWithTwoComponentComma) || 
               std::regex_search(lowerStr, matches, baseWithTwoComponent)) {
      std::string tempName = var_name.s_str();
      tempName.resize(matches[1].str().size());
      tempName.erase(std::remove(tempName.begin(), tempName.end(), ' '), tempName.end());
      baseName = tempName;
      comp1 = VariableComponent(matches[2].str(), matches[1].str());
      comp2 = VariableComponent(matches[3].str(), matches[1].str());
      explictNumComponents = 2;
    }

    if (comp1.UnknownComponent() || comp2.UnknownComponent() || baseName.empty()) {
      std::ostringstream error;
      error << "The variable '" << var_name << "' could not be parsed out of the variable list: "<< tag <<  ".\n"
            << "  Valid syntax takes one of 3 forms:\n"
            << "  variable_name, variable_name(<component>), or variable_name(<component>,<intg-pt>).\n"
            << "  <component> can be one of [:,*,x,y,z,xx,yy,zz,xy,xz,yx,yz,zx,xy] or a positive integer. \n"
            << "  <intg-pt>   can be one of [:,*] or a positive integer. \n"
            << ErrorTrace;
        throw std::runtime_error(error.str());
    }

  }

  sierra::String ParsedVariable::name() const {
    if (explictNumComponents == 0) {
      return baseName;
    } if (explictNumComponents == 1) {
      return baseName + "(" + comp1.rangeString() + ")";
    } if (explictNumComponents == 2) {
      return baseName + "(" + comp1.rangeString() + "," + comp2.rangeString() + ")";
    }
    return "";
  }

}
}