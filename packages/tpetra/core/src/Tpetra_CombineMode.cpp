// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_CombineMode.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

namespace Tpetra {

  void
  setCombineModeParameter (Teuchos::ParameterList& plist,
                           const std::string& paramName)
  {
    typedef Tpetra::CombineMode enum_type;
    typedef Teuchos::StringToIntegralParameterEntryValidator<enum_type>
      validator_type;

    const std::string docString = "Tpetra::CombineMode: rule for combining "
      "entries that overlap across processes, when redistributing data via a "
      "Tpetra::Import or Tpetra::Export";
    const std::string defaultVal = "ADD";
    const bool caseSensitive = false;

    const Teuchos::Array<std::string>::size_type numParams = 6;
    Teuchos::Array<std::string> strs (numParams);
    Teuchos::Array<std::string> docs (numParams);
    Teuchos::Array<enum_type> vals (numParams);

    strs[0] = "ADD";
    strs[1] = "INSERT";
    strs[2] = "REPLACE";
    strs[3] = "ABSMAX";
    strs[4] = "ZERO";
    strs[5] = "ADD_ASSIGN";

    docs[0] = "Sum new values";
    docs[1] = "Insert new values that don't currently exist";
    docs[2] = "Replace existing values with new values";
    docs[3] = "Replace old value with maximum of magnitudes of old and new values";
    docs[4] = "Replace old values with zero";
    docs[5] = "Do addition assignment (+=) of new values into existing value; "
              "may not be supported by all classes";

    vals[0] = ADD;
    vals[1] = INSERT;
    vals[2] = REPLACE;
    vals[3] = ABSMAX;
    vals[4] = ZERO;
    vals[5] = ADD_ASSIGN;

    plist.set (paramName, defaultVal, docString,
               Teuchos::rcp (new validator_type (strs (), docs (), vals (),
                                                 defaultVal, caseSensitive)));
  }

  std::string combineModeToString (const CombineMode combineMode)
  {
    std::string combineModeStr;
    switch (combineMode) {
    case ADD:
      combineModeStr = "ADD";
      break;
    case REPLACE:
      combineModeStr = "REPLACE";
      break;
    case ABSMAX:
      combineModeStr = "ABSMAX";
      break;
    case INSERT:
      combineModeStr = "INSERT";
      break;
    case ZERO:
      combineModeStr = "ZERO";
      break;
    case ADD_ASSIGN:
      combineModeStr = "ADD_ASSIGN";
      break;
    default:
      combineModeStr = "INVALID";
    }
    return combineModeStr;
  }

} // namespace Tpetra
