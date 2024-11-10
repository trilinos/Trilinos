// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// clang-format on
#include <algorithm> // std::transform
#include <array>
#include <cctype>  // std::toupper
#include <cstdlib> // std::getenv
#include <functional>
#include <map>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_TestForException.hpp"
#include "TpetraCore_config.h"
#include "Tpetra_Details_Behavior.hpp"
#include "KokkosKernels_config.h"  // for TPL enable macros

/*! \file Tpetra_Details_Behavior.cpp

This file is primarily concerned with providing programmatic access to certain
environment variables. The variable is evaluated on first access, and the data
is cached to provide the same result on future accesses even if the underlying
environment has changed.

There are boolean variables:

  TPETRA_ENV_VAR=1|ON|YES|TRUE     <lower or upper case>
  TPETRA_ENV_VAR=0|OFF|NO|FALSE    <lower or upper case>

There are size_t variables (negative values are max size_t)

  TPETRA_ENV_VAR=XXXXX  <XXXX is parsed as a long long, then cast to size_t>

There are "named" boolean variables, where multiple names are associated with
the variable.

  TPETRA_ENV_VAR=name1:name2:....

If a name is queried for a variable, and it was provided, the result is true.
Else false.

There are also two special cases:

  `CUDA_LAUNCH_BLOCKING` is not prefixed with "TPETRA_" but is otherwise a
boolean variable

  `MM_TAFC_OptimizationCoreCount` is not prefixed with "TPETRA_",
and is parsed as an int instead of a long long

*/

// environ should be available on posix platforms
#if not(defined(WIN) && (_MSC_VER >= 1900))
// needs to be in the global namespace
extern char **environ;
#endif

namespace Tpetra {
namespace Details {

namespace BehaviorDetails {

constexpr const std::string_view RESERVED_PREFIX = "TPETRA_";
constexpr const std::string_view ASSUME_GPU_AWARE_MPI =
    "TPETRA_ASSUME_GPU_AWARE_MPI";
constexpr const std::string_view CUDA_LAUNCH_BLOCKING = "CUDA_LAUNCH_BLOCKING";
constexpr const std::string_view MM_TAFC_OptimizationCoreCount =
    "MM_TAFC_OptimizationCoreCount";
constexpr const std::string_view VERBOSE_PRINT_COUNT_THRESHOLD =
    "TPETRA_VERBOSE_PRINT_COUNT_THRESHOLD";
constexpr const std::string_view ROW_IMBALANCE_THRESHOLD =
    "TPETRA_ROW_IMBALANCE_THRESHOLD";
constexpr const std::string_view MULTIVECTOR_USE_MERGE_PATH =
    "TPETRA_MULTIVECTOR_USE_MERGE_PATH";
constexpr const std::string_view VECTOR_DEVICE_THRESHOLD =
    "TPETRA_VECTOR_DEVICE_THRESHOLD";
constexpr const std::string_view HIERARCHICAL_UNPACK_BATCH_SIZE =
    "TPETRA_HIERARCHICAL_UNPACK_BATCH_SIZE";
constexpr const std::string_view HIERARCHICAL_UNPACK_TEAM_SIZE =
    "TPETRA_HIERARCHICAL_UNPACK_TEAM_SIZE";
constexpr const std::string_view USE_TEUCHOS_TIMERS =
    "TPETRA_USE_TEUCHOS_TIMERS";
constexpr const std::string_view USE_KOKKOS_PROFILING =
    "TPETRA_USE_KOKKOS_PROFILING";
constexpr const std::string_view DEBUG = "TPETRA_DEBUG";
constexpr const std::string_view VERBOSE = "TPETRA_VERBOSE";
constexpr const std::string_view TIMING = "TPETRA_TIMING";
constexpr const std::string_view HIERARCHICAL_UNPACK =
    "TPETRA_HIERARCHICAL_UNPACK";
constexpr const std::string_view SKIP_COPY_AND_PERMUTE =
    "TPETRA_SKIP_COPY_AND_PERMUTE";
constexpr const std::string_view FUSED_RESIDUAL = "TPETRA_FUSED_RESIDUAL";
constexpr const std::string_view OVERLAP = "TPETRA_OVERLAP";
constexpr const std::string_view SPACES_ID_WARN_LIMIT =
    "TPETRA_SPACES_ID_WARN_LIMIT";
constexpr const std::string_view TIME_KOKKOS_DEEP_COPY =
    "TPETRA_TIME_KOKKOS_DEEP_COPY";
constexpr const std::string_view TIME_KOKKOS_DEEP_COPY_VERBOSE1 =
    "TPETRA_TIME_KOKKOS_DEEP_COPY_VERBOSE1";
constexpr const std::string_view TIME_KOKKOS_DEEP_COPY_VERBOSE2 =
    "TPETRA_TIME_KOKKOS_DEEP_COPY_VERBOSE2";
constexpr const std::string_view TIME_KOKKOS_FENCE = "TPETRA_TIME_KOKKOS_FENCE";
constexpr const std::string_view TIME_KOKKOS_FUNCTIONS =
    "TPETRA_TIME_KOKKOS_FUNCTIONS";

// construct an std::array of string_view with any number of provided
// string_views
template <typename... Elems>
constexpr std::array<std::string_view, sizeof...(Elems)>
make_array(Elems &&... elems) {
  return {std::forward<Elems>(elems)...};
}

constexpr const auto RECOGNIZED_VARS = make_array(
    ASSUME_GPU_AWARE_MPI, CUDA_LAUNCH_BLOCKING, MM_TAFC_OptimizationCoreCount,
    VERBOSE_PRINT_COUNT_THRESHOLD, ROW_IMBALANCE_THRESHOLD,
    MULTIVECTOR_USE_MERGE_PATH, VECTOR_DEVICE_THRESHOLD,
    HIERARCHICAL_UNPACK_BATCH_SIZE, HIERARCHICAL_UNPACK_TEAM_SIZE,
    USE_TEUCHOS_TIMERS, USE_KOKKOS_PROFILING, DEBUG, VERBOSE, TIMING,
    HIERARCHICAL_UNPACK, SKIP_COPY_AND_PERMUTE, FUSED_RESIDUAL, OVERLAP,
    SPACES_ID_WARN_LIMIT, TIME_KOKKOS_DEEP_COPY, TIME_KOKKOS_DEEP_COPY_VERBOSE1,
    TIME_KOKKOS_DEEP_COPY_VERBOSE2, TIME_KOKKOS_FENCE, TIME_KOKKOS_FUNCTIONS);

// clang-format off
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
  split(const std::string_view s,
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
// clang-format on

template <typename T>
T getEnvironmentVariable(const std::string_view environmentVariableName,
                         const T defaultValue) {
  const char prefix[] = "Tpetra::Details::Behavior: ";

  const char *varVal = std::getenv(environmentVariableName.data());
  if (varVal == nullptr) {
    return defaultValue;
  } else {
    std::stringstream ss(varVal);
    T parsed;
    ss >> parsed;

    TEUCHOS_TEST_FOR_EXCEPTION(!ss, std::out_of_range,
                               prefix << "Environment "
                                         "variable \""
                                      << environmentVariableName
                                      << "\" has a "
                                         "value "
                                      << varVal
                                      << " that cannot be parsed as a "
                                      << typeid(T).name() << ".");

    return parsed;
  }
}

// full specialization of bool to preserve historical Tpetra parsing behavior
template <>
bool getEnvironmentVariable<bool>(
    const std::string_view environmentVariableName, const bool defaultValue) {
  const char *varVal = std::getenv(environmentVariableName.data());
  bool retVal = defaultValue;
  if (varVal != nullptr) {
    auto state = environmentVariableState(std::string(varVal));
    if (state == EnvironmentVariableIsSet_ON)
      retVal = true;
    else if (state == EnvironmentVariableIsSet_OFF)
      retVal = false;
  }
  return retVal;
}

/*! full specialization of size_t to preserve historical Tpetra parsing behavior

Parse it as a long long. If negative, return max size_t.
Else, return cast to size_t
*/
template <>
size_t
getEnvironmentVariable<size_t>(const std::string_view environmentVariableName,
                               const size_t defaultValue) {
  const char prefix[] = "Tpetra::Details::Behavior: ";

  const char *varVal = std::getenv(environmentVariableName.data());
  if (varVal == nullptr) {
    return defaultValue;
  } else {
    long long val = std::stoll(stringToUpper(varVal));
    if (val < static_cast<long long>(0)) {
      // If negative - user has requested threshold be lifted
      return std::numeric_limits<size_t>::max();
    }
    if (sizeof(long long) > sizeof(size_t)) {
      // It's hard to test this code, but I want to try writing it
      // at least, in case we ever have to run on 32-bit machines or
      // machines with sizeof(long long)=16 and sizeof(size_t)=8.
      constexpr long long maxSizeT =
          static_cast<long long>(std::numeric_limits<size_t>::max());
      TEUCHOS_TEST_FOR_EXCEPTION(
          val > maxSizeT, std::out_of_range,
          prefix << "Environment "
                    "variable \""
                 << environmentVariableName
                 << "\" has a "
                    "value "
                 << val << " larger than the largest size_t value " << maxSizeT
                 << ".");
    }
    return static_cast<size_t>(val);
  }
}

template <typename T>
T idempotentlyGetEnvironmentVariable(
    T &value, bool &initialized, const std::string_view environmentVariableName,
    const T defaultValue) {
  if (!initialized) {
    value = getEnvironmentVariable<T>(environmentVariableName, defaultValue);
    initialized = true;
  }
  return value;
}

// clang-format off
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

  constexpr bool assumeMpiIsGPUAwareDefault () {
#ifdef TPETRA_ASSUME_GPU_AWARE_MPI
    return true;
#else
    return false;
#endif // TPETRA_ASSUME_GPU_AWARE_MPI
  }

  constexpr bool cudaLaunchBlockingDefault () {
    return false;
  }

  constexpr bool hierarchicalUnpackDefault () {
    return true;
  }

} // namespace (anonymous)
// clang-format on

void Behavior::reject_unrecognized_env_vars() {

  static bool once = false;

  if (!once) {
    const char prefix[] = "Tpetra::Details::Behavior: ";
    char **env;
#if defined(WIN) && (_MSC_VER >= 1900)
    env = *__p__environ();
#else
    env = environ; // defined at the top of this file as extern char **environ;
#endif
    for (; *env; ++env) {

      std::string name;
      std::string value;
      const std::string_view ev(*env);

      // split name=value on the first =, everything before = is name
      split(
          ev,
          [&](const std::string &s) {
            if (name.empty()) {
              name = s;
            } else {
              value = s;
            }
          },
          '=');

      if (name.size() >= BehaviorDetails::RESERVED_PREFIX.size() &&
          name.substr(0, BehaviorDetails::RESERVED_PREFIX.size()) ==
              BehaviorDetails::RESERVED_PREFIX) {
        const auto it = std::find(BehaviorDetails::RECOGNIZED_VARS.begin(),
                                  BehaviorDetails::RECOGNIZED_VARS.end(), name);
        TEUCHOS_TEST_FOR_EXCEPTION(
            it == BehaviorDetails::RECOGNIZED_VARS.end(), std::out_of_range,
            prefix << "Environment "
                      "variable \""
                   << name << "\" (prefixed with \""
                   << BehaviorDetails::RESERVED_PREFIX
                   << "\") is not a recognized Tpetra variable.");
      }
    }

    once = true;
  }
}

bool Behavior::debug() {
  constexpr bool defaultValue = debugDefault();

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::DEBUG, defaultValue);
}

bool Behavior::verbose() {
  if (BehaviorDetails::verboseDisabled_)
    return false;

  constexpr bool defaultValue = verboseDefault();

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::VERBOSE, defaultValue);
}

bool Behavior::timing() {
  if (BehaviorDetails::timingDisabled_)
    return false;

  constexpr bool defaultValue = timingDefault();

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::TIMING, defaultValue);
}

bool Behavior::assumeMpiIsGPUAware() {
  constexpr bool defaultValue = assumeMpiIsGPUAwareDefault();

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::ASSUME_GPU_AWARE_MPI,
      defaultValue);
}

bool Behavior::cudaLaunchBlocking() {
  constexpr bool defaultValue = cudaLaunchBlockingDefault();

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::CUDA_LAUNCH_BLOCKING,
      defaultValue);
}

int Behavior::TAFC_OptimizationCoreCount() {
  constexpr int _default = 3000;
  static int value_ = _default;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::MM_TAFC_OptimizationCoreCount,
      _default);
}

size_t Behavior::verbosePrintCountThreshold() {
  constexpr size_t defaultValue(200);

  static size_t value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::VERBOSE_PRINT_COUNT_THRESHOLD,
      defaultValue);
}

size_t Behavior::rowImbalanceThreshold() {
  constexpr size_t defaultValue(256);

  static size_t value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::ROW_IMBALANCE_THRESHOLD,
      defaultValue);
}

bool Behavior::useMergePathMultiVector() {
  constexpr bool defaultValue = false;

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::MULTIVECTOR_USE_MERGE_PATH,
      defaultValue);
}

size_t Behavior::multivectorKernelLocationThreshold() {
  constexpr size_t defaultValue(22000);

  static size_t value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::VECTOR_DEVICE_THRESHOLD,
      defaultValue);
}

size_t Behavior::hierarchicalUnpackBatchSize() {

#ifdef HAVE_TPETRA_INST_CUDA
  constexpr size_t defaultValue(16);
#else
  constexpr size_t defaultValue(256);
#endif

  static size_t value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::HIERARCHICAL_UNPACK_BATCH_SIZE,
      defaultValue);
}

size_t Behavior::hierarchicalUnpackTeamSize() {
#ifdef HAVE_TPETRA_INST_CUDA
  const size_t defaultValue(16);
#else
  const size_t defaultValue(Teuchos::OrdinalTraits<size_t>::invalid());
#endif

  static size_t value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::HIERARCHICAL_UNPACK_TEAM_SIZE,
      defaultValue);
}

bool Behavior::profilingRegionUseTeuchosTimers() {
  constexpr bool defaultValue(false);

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::USE_TEUCHOS_TIMERS, defaultValue);
}

bool Behavior::profilingRegionUseKokkosProfiling() {
  constexpr bool defaultValue(false);

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::USE_KOKKOS_PROFILING,
      defaultValue);
}

bool Behavior::debug(const char name[]) {
  constexpr bool defaultValue = false;

  static bool initialized_ = false;
  return idempotentlyGetNamedEnvironmentVariableAsBool(
      name, initialized_, BehaviorDetails::DEBUG.data(), defaultValue);
}

bool Behavior::verbose(const char name[]) {
  if (BehaviorDetails::verboseDisabled_)
    return false;

  constexpr bool defaultValue = false;

  static bool initialized_ = false;
  return idempotentlyGetNamedEnvironmentVariableAsBool(
      name, initialized_, BehaviorDetails::VERBOSE.data(), defaultValue);
}

void Behavior::enable_verbose_behavior() {
  BehaviorDetails::verboseDisabled_ = false;
}

void Behavior::disable_verbose_behavior() {
  BehaviorDetails::verboseDisabled_ = true;
}

bool Behavior::timing(const char name[]) {
  if (BehaviorDetails::timingDisabled_)
    return false;

  constexpr bool defaultValue = false;

  static bool initialized_ = false;
  return idempotentlyGetNamedEnvironmentVariableAsBool(
      name, initialized_, BehaviorDetails::TIMING.data(), defaultValue);
}

void Behavior::enable_timing() { BehaviorDetails::timingDisabled_ = false; }

void Behavior::disable_timing() { BehaviorDetails::timingDisabled_ = true; }

bool Behavior::hierarchicalUnpack() {
  constexpr bool defaultValue = hierarchicalUnpackDefault();

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::HIERARCHICAL_UNPACK, defaultValue);
}

bool Behavior::skipCopyAndPermuteIfPossible() {
  constexpr bool defaultValue(false);

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::SKIP_COPY_AND_PERMUTE,
      defaultValue);
}

bool Behavior::fusedResidual() {
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) || \
    defined(KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE) || \
    defined(KOKKOSKERNELS_ENABLE_TPL_MKL)
  constexpr bool defaultValue(false);
#else
  constexpr bool defaultValue(true);
#endif

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::FUSED_RESIDUAL, defaultValue);
}

bool Behavior::overlapCommunicationAndComputation() {
  constexpr bool defaultValue(false);

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::OVERLAP, defaultValue);
}

size_t Behavior::spacesIdWarnLimit() {
  constexpr size_t defaultValue(16);

  static size_t value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::SPACES_ID_WARN_LIMIT,
      defaultValue);
}

bool Behavior::timeKokkosDeepCopy() {
  constexpr bool defaultValue(false);

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::TIME_KOKKOS_DEEP_COPY,
      defaultValue);
}

bool Behavior::timeKokkosDeepCopyVerbose1() {
  constexpr bool defaultValue(false);

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::TIME_KOKKOS_DEEP_COPY_VERBOSE1,
      defaultValue);
}

bool Behavior::timeKokkosDeepCopyVerbose2() {
  constexpr bool defaultValue(false);

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::TIME_KOKKOS_DEEP_COPY_VERBOSE2,
      defaultValue);
}

bool Behavior::timeKokkosFence() {
  constexpr bool defaultValue(false);

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::TIME_KOKKOS_FENCE, defaultValue);
}

bool Behavior::timeKokkosFunctions() {
  constexpr bool defaultValue(false);

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::TIME_KOKKOS_FUNCTIONS,
      defaultValue);
}

} // namespace Details
} // namespace Tpetra
