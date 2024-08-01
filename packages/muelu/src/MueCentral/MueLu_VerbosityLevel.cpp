// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MueLu_VerbosityLevel.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Utilities.hpp"
#include <string>
#include <locale>

namespace MueLu {

VerbLevel toMueLuVerbLevel(const Teuchos::EVerbosityLevel verbLevel) {
  switch (verbLevel) {
    case Teuchos::VERB_NONE:
      return None;
    case Teuchos::VERB_DEFAULT:
      return Default;
    case Teuchos::VERB_LOW:
      return Low;
    case Teuchos::VERB_MEDIUM:
      return Medium;
    case Teuchos::VERB_HIGH:
      return High;
    case Teuchos::VERB_EXTREME:
      return Extreme;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Unknown enum value found.");
  }
}

std::string
lowerCase(const std::string& s) {
  typedef std::string::value_type char_t;
  typedef std::ctype<char_t> facet_type;
  const facet_type& facet = std::use_facet<facet_type>(std::locale());

  const std::string::size_type len = s.size();
  std::string s_lc(s);
  for (std::string::size_type k = 0; k < len; ++k) {
    s_lc[k] = facet.tolower(s[k]);
  }

  return s_lc;
}

MsgType toVerbLevel(const std::string& verbLevelStr) {
  std::map<std::string, MsgType> verbMap;
  // for developers
  verbMap["errors"]         = Errors;
  verbMap["warnings0"]      = Warnings0;
  verbMap["warnings00"]     = Warnings00;
  verbMap["warnings1"]      = Warnings1;
  verbMap["perfWarnings"]   = PerfWarnings;
  verbMap["runtime0"]       = Runtime0;
  verbMap["runtime1"]       = Runtime1;
  verbMap["runtimeTimings"] = RuntimeTimings;
  verbMap["noTimeReport"]   = NoTimeReport;
  verbMap["parameters0"]    = Parameters0;
  verbMap["parameters1"]    = Parameters1;
  verbMap["statistics0"]    = Statistics0;
  verbMap["statistics1"]    = Statistics1;
  verbMap["timings0"]       = Timings0;
  verbMap["timings1"]       = Timings1;
  verbMap["timingsByLevel"] = TimingsByLevel;
  verbMap["external"]       = External;
  verbMap["developer"]      = Developer;
  verbMap["debug"]          = Debug;
  verbMap["test"]           = Test;

  verbMap["warnings"]      = Warnings;
  verbMap["runtime"]       = Runtime;
  verbMap["parameters"]    = Parameters;
  verbMap["statistics"]    = Statistics;
  verbMap["timings"]       = Timings;
  verbMap["test"]          = Test;
  verbMap["interfacetest"] = InterfaceTest;
  // for users and developers
  verbMap["none"]    = None;
  verbMap["low"]     = Low;
  verbMap["medium"]  = Medium;
  verbMap["high"]    = High;
  verbMap["extreme"] = Extreme;

  std::string lcVerb = lowerCase(verbLevelStr);
  if (verbMap.find(lcVerb) != verbMap.end())
    return verbMap[lcVerb];
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::ParameterListInterpreter():: invalid verbosity level: " << verbLevelStr);
}

}  // namespace MueLu
