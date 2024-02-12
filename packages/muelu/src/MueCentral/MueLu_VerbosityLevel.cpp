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
