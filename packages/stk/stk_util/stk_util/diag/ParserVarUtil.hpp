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

#ifndef STK_UTIL_DIAG_ParserVarUtil_hpp
#define STK_UTIL_DIAG_ParserVarUtil_hpp

#include <string>  // for string
#include <vector>
#include "stk_util/diag/String.hpp"
#include <regex>
namespace stk {
namespace diag {
  /*!
   *  Parser reads a string like: 'stress(1,2)' as two strings,
   *  'stress(1' and '2)'.  This routines checks if two strings
   *  form a single variable name and if so concats them back
   *  together.  Returns true if the strings were concated,
   *  false otherwise.
   */

   bool remove_entry(const std::vector<sierra::String> &valid,  const sierra::String& input);
   std::string strip_entry(const std::vector<sierra::String> &valid,  const sierra::String& input);
   std::vector<sierra::String> remove_var_eq(const std::vector<sierra::String>& input);
   std::vector<sierra::String> join_paren_vars(const std::vector<sierra::String>& input, sierra::String tag  = "");

   class VariableComponent {
    public:
     VariableComponent() = default;
     VariableComponent(const std::string &rangeString, const sierra::String &varName);

     enum ComponentEnum {
       UNKNOWN_COMPONENT = 0,
       UNSPECIFIED_COMPONENT,
       ALL_COMPONENTS,
       VECTOR_X_COMPONENT,
       VECTOR_Y_COMPONENT,
       VECTOR_Z_COMPONENT,
       TENSOR_XX_COMPONENT,
       TENSOR_YY_COMPONENT,
       TENSOR_ZZ_COMPONENT,
       TENSOR_XY_COMPONENT,
       TENSOR_XZ_COMPONENT,
       TENSOR_YZ_COMPONENT,
       TENSOR_YX_COMPONENT,
       TENSOR_ZX_COMPONENT,
       TENSOR_ZY_COMPONENT,
       INTEGER_COMPONENT
     };

     inline bool UnknownComponent() const { return m_componentEnum == UNKNOWN_COMPONENT; }
     inline bool AllComponents() const { return m_componentEnum == ALL_COMPONENTS; }
     inline bool UnspecifiedComponent() const { return m_componentEnum == UNSPECIFIED_COMPONENT; }

     inline bool SpecifiedComponent() const { return m_componentEnum > ALL_COMPONENTS; }
     inline ComponentEnum ComponentEnumVal() const { return m_componentEnum; }

     sierra::String rangeString() const { return m_savedRangeString; }
     int getIntegerComponent() const { return m_integerComponent;}
    private:
     ComponentEnum m_componentEnum{UNKNOWN_COMPONENT};
     int m_integerComponent{-1};
     sierra::String m_savedRangeString{""};
   };

   // Extract a user specified variable and sub components
   struct ParsedVariable {
     ParsedVariable() = default;
     ParsedVariable( const ParsedVariable& ) = default;
     ParsedVariable(const sierra::String &var_name, sierra::String tag  = "");
     sierra::String name() const;

     sierra::String baseName{""};
     stk::diag::VariableComponent comp1;
     stk::diag::VariableComponent comp2;
     int explictNumComponents{-1};
   };

}
}

#endif
