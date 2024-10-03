// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_STRING_UTILITIES_HPP
#define PANZER_STRING_UTILITIES_HPP

#include <vector>
#include <string>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace panzer {

  //! Removes whitespace at beginning and end of string
  void trim(std::string& str);
  
  //! Tokenize a string, put tokens in a vector
  void StringTokenizer(std::vector<std::string>& tokens,
		       const std::string& str,
		       const std::string delimiter = ",",bool trim=false);

  //! Turn a vector of tokens into a vector of doubles
  void TokensToDoubles(std::vector<double> & values,const std::vector<std::string> & tokens);

  //! Turn a vector of tokens into a vector of ints
  void TokensToInts(std::vector<int> & values,const std::vector<std::string> & tokens);

  /** Read in a parameter field and return the correct scalar field. This parses
    * scalar type data
    */
  template <typename ScalarT>
  ScalarT getScalarParameter(const std::string & field,const Teuchos::ParameterList & plist)
  {
      return plist.get<double>(field);
  }
}

#endif
