// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_CommaSeparatedEntryValidator_hpp__
#define __Panzer_CommaSeparatedEntryValidator_hpp__

#include <map>
#include <set>
#include <vector>
#include <string>
#include <iostream>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_ParameterEntryValidator.hpp"
#include "Teuchos_RCP.hpp"

namespace panzer {

/** This class validates a response type. Essentially
  * it is used to make sure the parameter value is correctly
  * formatted. 
  */
class CommaSeparatedEntryValidator : public Teuchos::ParameterEntryValidator {
public:
  /** A basic constructor. If <code>allowEmpty</code> is true then the 
    * empty string is a valid entry.
    */
  CommaSeparatedEntryValidator(bool allowEmpty=false) : allowEmpty_(allowEmpty) {}

  ValidStringsList validStringValues() const
  { return Teuchos::null; }

  void validate(const Teuchos::ParameterEntry & entry, 
                const std::string & paramName,
                const std::string & sublistName) const;

  const std::string getXMLTypeName() const
  { return "string-list"; }

  void printDoc(const std::string & docString, std::ostream &out) const;

  //! Utility function for tokenizing
  static void split(const std::string & str,
                    const std::string & delim,
                    std::vector<std::string> & tokens);
private:

  bool allowEmpty_; //! Is an empty string valid?
};

}

#endif
