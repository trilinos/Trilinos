/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
// @HEADER
*/
#include "Tpetra_Details_Behavior.hpp"
#include "TpetraCore_config.h"
#include <algorithm> // std::transform
#include <cstdlib> // std::getenv
#include <cctype> // std::toupper
#include <string>
#include <map>
#include <vector>
#include <functional>
#include "Teuchos_TestForException.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include <stdexcept>

namespace Tpetra {
namespace Details {

namespace BehaviorDetails {
std::map<std::string, std::map<std::string, bool> > namedVariableMap_;
bool verboseDisabled_ = false;
bool timingDisabled_ = false;
}

namespace { // (anonymous)

  enum EnvironmentVariableState
  {
    EnvironmentVariableIsSet_ON,
    EnvironmentVariableIsSet_OFF,
    EnvironmentVariableIsSet,
    EnvironmentVariableIsNotSet
  };

  // See example here:
  //
  // http://en.cppreference.com/w/cpp/string/byte/toupper
  std::string stringToUpper (std::string s)
  {
    std::transform (s.begin (), s.end (), s.begin (),
                    [] (unsigned char c) { return std::toupper (c); });
    return s;
  }

  void
  split(const std::string& s,
        std::function<void(const std::string&)> f,
        const char sep=',')
  {
    typedef std::string::size_type size_type;
    size_type cur_pos, last_pos=0, length=s.length();
    while(last_pos < length + 1)
    {
      cur_pos = s.find_first_of(sep, last_pos);
      if(cur_pos == std::string::npos)
      {
        cur_pos = length;
      }
      if(cur_pos!=last_pos) {
        auto token = std::string(s.data()+last_pos, (size_type)cur_pos-last_pos);
        f(token);
      }
      last_pos = cur_pos + 1;
     }
    return;
  }

  EnvironmentVariableState
  environmentVariableState(const std::string& environmentVariableValue)
  {
    std::string v = stringToUpper(environmentVariableValue);
    if      (v == "1" || v == "YES" || v == "TRUE"  || v == "ON")
      // Environment variable is "ON"
      return EnvironmentVariableIsSet_ON;
    else if (v == "0" || v == "NO"  || v == "FALSE" || v == "OFF")
      // Environment variable is "OFF"
      return EnvironmentVariableIsSet_OFF;
    // Environment has some other non-boolean value
    return EnvironmentVariableIsSet;
  }

  void
  setEnvironmentVariableMap (const char environmentVariableName[],
                             std::map<std::string,std::map<std::string, bool> >& valsMap,
                             const bool defaultValue)
  {
    using std::map;
    using std::getenv;
    using std::string;
    using std::vector;

    // Set the default value for this variable
    valsMap[environmentVariableName] = map<string,bool>{{"DEFAULT", defaultValue}};

    const char* varVal = getenv (environmentVariableName);
    if (varVal == nullptr) {
      // Environment variable is not set, use the default value for any named
      // variants
      return;
    }

    // Variable is not empty.
    const string varStr(varVal);
    vector<string> names;
    split(varStr, [&](const string& x){names.push_back(x);});
    for (auto const& name: names) {
      auto state = environmentVariableState(name);
      if (state == EnvironmentVariableIsSet_ON) {
        // Environment variable was set as ENVAR_NAME=[1,YES,TRUE,ON]
        // Global value takes precedence
        valsMap[environmentVariableName]["DEFAULT"] = true;
      }
      else if (state == EnvironmentVariableIsSet_OFF) {
        // Environment variable was set as ENVAR_NAME=[0,NO,FALSE,OFF]
        // Global value takes precedence
        valsMap[environmentVariableName]["DEFAULT"] = false;
      }
      else {
        // Environment variable was set as ENVAR_NAME=...:name:...
        // So we set the mapping true for this named variant
        valsMap[environmentVariableName][name] = true;
      }
    }
    return;
  }

  bool
  getEnvironmentVariableAsBool (const char environmentVariableName[],
                                const bool defaultValue)
  {
    const char* varVal = std::getenv (environmentVariableName);

    bool retVal = defaultValue;
    if (varVal != nullptr) {
      auto state = environmentVariableState(std::string(varVal));
      if (state == EnvironmentVariableIsSet_ON) retVal = true;
      else if (state == EnvironmentVariableIsSet_OFF) retVal = false;
    }
    return retVal;
  }

  size_t
  getEnvironmentVariableAsSize(const char environmentVariableName[],
                               const size_t defaultValue)
  {
    const char prefix[] = "Tpetra::Details::Behavior: ";

    const char* varVal = std::getenv(environmentVariableName);
    if (varVal == nullptr) {
      return defaultValue;
    }
    else {
      // This could throw invalid_argument or out_of_range.
      // Go ahead and let it do so.
      long long val = std::stoll(stringToUpper(varVal));
      if (val < static_cast<long long>(0)) {
        // If negative - user has requested threshold be lifted
        return std::numeric_limits<size_t>::max();
      }
//      TEUCHOS_TEST_FOR_EXCEPTION
//        (val < static_cast<long long>(0), std::out_of_range,
//         prefix << "Environment variable \""
//         << environmentVariableName << "\" is supposed to be a size, "
//         "but it has a negative integer value " << val << ".");
      if (sizeof(long long) > sizeof(size_t)) {
        // It's hard to test this code, but I want to try writing it
        // at least, in case we ever have to run on 32-bit machines or
        // machines with sizeof(long long)=16 and sizeof(size_t)=8.
        constexpr long long maxSizeT =
          static_cast<long long>(std::numeric_limits<size_t>::max());
        TEUCHOS_TEST_FOR_EXCEPTION
          (val > maxSizeT, std::out_of_range, prefix << "Environment "
           "variable \"" << environmentVariableName << "\" has a "
           "value " << val << " larger than the largest size_t value "
           << maxSizeT << ".");
      }
      return static_cast<size_t>(val);
    }
  }

  bool
  idempotentlyGetEnvironmentVariableAsBool (bool& value,
                                            bool& initialized,
                                            const char environmentVariableName[],
                                            const bool defaultValue)
  {
    if (! initialized) {
      value = getEnvironmentVariableAsBool (environmentVariableName,
                                            defaultValue);
      initialized = true;
    }
    return value;
  }

  bool
  idempotentlyGetNamedEnvironmentVariableAsBool (const char name[],
                                                 bool& initialized,
                                                 const char environmentVariableName[],
                                                 const bool defaultValue)
  {
    using BehaviorDetails::namedVariableMap_;
    if (! initialized) {
      setEnvironmentVariableMap (environmentVariableName,
                                 namedVariableMap_,
                                 defaultValue);
      initialized = true;
    }
    auto thisEnvironmentVariableMap = namedVariableMap_[environmentVariableName];
    auto thisEnvironmentVariable = thisEnvironmentVariableMap.find(name);
    if (thisEnvironmentVariable != thisEnvironmentVariableMap.end())
      return thisEnvironmentVariable->second;
    return thisEnvironmentVariableMap["DEFAULT"];
  }

  size_t
  idempotentlyGetEnvironmentVariableAsSize
    (size_t& value,
     bool& initialized,
     const char environmentVariableName[],
     const size_t defaultValue)
  {
    if (! initialized) {
      value = getEnvironmentVariableAsSize(environmentVariableName,
                                           defaultValue);
      initialized = true;
    }
    return value;
  }

  constexpr bool debugDefault () {
#ifdef HAVE_TPETRA_DEBUG
    return true;
#else
    return false;
#endif // HAVE_TPETRA_DEBUG
  }

  constexpr bool verboseDefault () {
    return false;
  }

  constexpr bool timingDefault () {
    return false;
  }

  constexpr bool assumeMpiIsCudaAwareDefault () {
#ifdef TPETRA_ASSUME_CUDA_AWARE_MPI
    return true;
#else
    return false;
#endif // TPETRA_ASSUME_CUDA_AWARE_MPI
  }

  constexpr bool hierarchicalUnpackDefault () {
    return true;
  }

} // namespace (anonymous)

bool Behavior::debug ()
{
  constexpr char envVarName[] = "TPETRA_DEBUG";
  constexpr bool defaultValue = debugDefault ();

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariableAsBool (value_,
                                                   initialized_,
                                                   envVarName,
                                                   defaultValue);
}

bool Behavior::verbose ()
{
  if (BehaviorDetails::verboseDisabled_) return false;

  constexpr char envVarName[] = "TPETRA_VERBOSE";
  constexpr bool defaultValue = verboseDefault ();

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariableAsBool (value_,
                                                   initialized_,
                                                   envVarName,
                                                   defaultValue);
}

bool Behavior::timing ()
{
  if (BehaviorDetails::timingDisabled_) return false;

  constexpr char envVarName[] = "TPETRA_TIMING";
  constexpr bool defaultValue = timingDefault ();

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariableAsBool (value_,
                                                   initialized_,
                                                   envVarName,
                                                   defaultValue);
}

bool Behavior::assumeMpiIsCudaAware ()
{
  constexpr char envVarName[] = "TPETRA_ASSUME_CUDA_AWARE_MPI";
  constexpr bool defaultValue = assumeMpiIsCudaAwareDefault ();

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariableAsBool (value_,
                                                   initialized_,
                                                   envVarName,
                                                   defaultValue);
}

int Behavior::TAFC_OptimizationCoreCount ()
{
    // only call getenv once, save the value.
    static int savedval=-1;
    if(savedval!=-1) return savedval;
    const char* varVal = std::getenv ("MM_TAFC_OptimizationCoreCount");
    if (varVal == nullptr) {
        savedval = 3000;
        return savedval;
    }
    savedval = std::stoi(std::string(varVal));
    return savedval;
}

size_t Behavior::verbosePrintCountThreshold ()
{
  constexpr char envVarName[] = "TPETRA_VERBOSE_PRINT_COUNT_THRESHOLD";
  constexpr size_t defaultValue (200);

  static size_t value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariableAsSize
    (value_, initialized_, envVarName, defaultValue);
}

size_t Behavior::rowImbalanceThreshold ()
{
  constexpr char envVarName[] = "TPETRA_ROW_IMBALANCE_THRESHOLD";
  constexpr size_t defaultValue (256);

  static size_t value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariableAsSize
    (value_, initialized_, envVarName, defaultValue);
}

bool Behavior::useMergePathMultiVector()
{
  constexpr char envVarName[] = "TPETRA_MULTIVECTOR_USE_MERGE_PATH";
  constexpr bool defaultValue = false;

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariableAsBool
    (value_, initialized_, envVarName, defaultValue);
}

  size_t Behavior::multivectorKernelLocationThreshold ()
{
  constexpr char envVarName[] = "TPETRA_VECTOR_DEVICE_THRESHOLD";
  constexpr size_t defaultValue (10000);

  static size_t value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariableAsSize
    (value_, initialized_, envVarName, defaultValue);
}

size_t Behavior::hierarchicalUnpackBatchSize ()
{
  constexpr char envVarName[] = "TPETRA_HIERARCHICAL_UNPACK_BATCH_SIZE";

#ifdef HAVE_TPETRA_INST_CUDA
  constexpr size_t defaultValue (16);
#else
  constexpr size_t defaultValue (256);
#endif

  static size_t value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariableAsSize
    (value_, initialized_, envVarName, defaultValue);
}

size_t Behavior::hierarchicalUnpackTeamSize ()
{
  constexpr char envVarName[] = "TPETRA_HIERARCHICAL_UNPACK_TEAM_SIZE";
#ifdef HAVE_TPETRA_INST_CUDA
  const size_t defaultValue (16);
#else
  const size_t defaultValue (Teuchos::OrdinalTraits<size_t>::invalid ());
#endif

  static size_t value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariableAsSize
    (value_, initialized_, envVarName, defaultValue);
}

bool Behavior::profilingRegionUseTeuchosTimers ()
{
  constexpr char envVarName[] = "TPETRA_USE_TEUCHOS_TIMERS";
  constexpr bool defaultValue(false);

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariableAsBool
    (value_, initialized_, envVarName, defaultValue);
}

bool Behavior::profilingRegionUseKokkosProfiling ()
{
  constexpr char envVarName[] = "TPETRA_USE_KOKKOS_PROFILING";
  constexpr bool defaultValue(false);

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariableAsBool
    (value_, initialized_, envVarName, defaultValue);
}


bool Behavior::debug (const char name[])
{
  constexpr char envVarName[] = "TPETRA_DEBUG";
  constexpr bool defaultValue = false;

  static bool initialized_ = false;
  return idempotentlyGetNamedEnvironmentVariableAsBool (name,
                                                        initialized_,
                                                        envVarName,
                                                        defaultValue);
}

bool Behavior::verbose (const char name[])
{
  if (BehaviorDetails::verboseDisabled_) return false;

  constexpr char envVarName[] = "TPETRA_VERBOSE";
  constexpr bool defaultValue = false;

  static bool initialized_ = false;
  return idempotentlyGetNamedEnvironmentVariableAsBool (name,
                                                        initialized_,
                                                        envVarName,
                                                        defaultValue);
}

void Behavior::enable_verbose_behavior () {
  BehaviorDetails::verboseDisabled_ = false;
}

void Behavior::disable_verbose_behavior () {
  BehaviorDetails::verboseDisabled_ = true;
}

bool Behavior::timing (const char name[])
{
  if (BehaviorDetails::timingDisabled_) return false;

  constexpr char envVarName[] = "TPETRA_TIMING";
  constexpr bool defaultValue = false;

  static bool initialized_ = false;
  return idempotentlyGetNamedEnvironmentVariableAsBool (name,
                                                        initialized_,
                                                        envVarName,
                                                        defaultValue);
}

void Behavior::enable_timing() {
  BehaviorDetails::timingDisabled_ = false;
}

void Behavior::disable_timing() {
  BehaviorDetails::timingDisabled_ = true;
}

bool Behavior::hierarchicalUnpack ()
{
  constexpr char envVarName[] = "TPETRA_HIERARCHICAL_UNPACK";
  constexpr bool defaultValue = hierarchicalUnpackDefault();

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariableAsBool (value_,
                                                   initialized_,
                                                   envVarName,
                                                   defaultValue);
}

} // namespace Details
} // namespace Tpetra

