// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_String_Utilities_hpp
#define Tempus_String_Utilities_hpp

#include <vector>
#include <string>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Tempus {

  //! Removes whitespace at beginning and end of string
  void trim(std::string& str);

  //! Tokenize a string, put tokens in a vector
  void StringTokenizer(std::vector<std::string>& tokens,
                       const std::string& str,
                       const std::string delimiter = ",",bool trim=false);

  //! Turn a vector of tokens into a vector of doubles
  void TokensToDoubles(std::vector<double> & values,
                   const std::vector<std::string> & tokens);

  //! Turn a vector of tokens into a vector of ints
  void TokensToInts(std::vector<int> & values,
                   const std::vector<std::string> & tokens);

  /** Read in a parameter field and return the correct scalar field. This parses
    * scalar type data
    */
  template <typename ScalarT>
  ScalarT getScalarParameter(const std::string & field,
                   const Teuchos::ParameterList & plist)
  {
      return plist.get<double>(field);
  }
}

#endif // Tempus_String_Utilities_hpp
