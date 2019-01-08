#include "Tpetra_Details_Behavior.hpp"
#include "TpetraCore_config.h"
#include <algorithm> // std::transform
#include <atomic> // std::atomic_thread_fence, std::memory_order_release
#include <cstdlib> // std::getenv
#include <cctype> // std::toupper
#include <mutex> // std::call_once, std::once_flag
#include <string>
#include <map>
#include <vector>
#include <functional>

namespace Tpetra {
namespace Details {

namespace BehaviorDetails {
std::map<std::string, std::map<std::string, bool> > namedVariableMap_;
bool verboseDisabled_ = false;
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
    if (varVal == NULL) {
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
    if (varVal != NULL) {
      auto state = environmentVariableState(std::string(varVal));
      if (state == EnvironmentVariableIsSet_ON) retVal = true;
      else if (state == EnvironmentVariableIsSet_OFF) retVal = false;
    }
    return retVal;
  }

  bool
  idempotentlyGetEnvironmentVariableAsBool (std::once_flag& once_flag,
                                            bool& value,
                                            bool& initialized,
                                            const char environmentVariableName[],
                                            const bool defaultValue)
  {
    // The extra "initialized" check avoids the cost of synchronizing
    // on the std::call_once for every call to this function.  We want
    // it to be cheap to get the Boolean value, so that users aren't
    // tempted to try to cache it themselves.
    if (! initialized) {
      std::call_once (once_flag, [&] () {
          value = getEnvironmentVariableAsBool (environmentVariableName,
                                                defaultValue);
          // http://preshing.com/20130922/acquire-and-release-fences/
          //
          // "A release fence prevents the memory reordering of any
          // read or write which precedes it in program order with any
          // write which follows it in program order."
          //
          // The point is to prevent the assignment to 'value' from
          // getting reordered after the assignment to 'initialized'
          // (the so-called "StoreStore" reordering).  That would be
          // bad in this case, because then other threads might read
          // 'initialized' as true, yet would fail to pick up the
          // change to 'value'.
          //
          // It's harmless if other threads don't see the write to
          // 'initialized', but did see the write to 'value'.  In that
          // case, they would just attempt and fail to enter the
          // std::call_once, and return (the correct value of)
          // 'value'.
          std::atomic_thread_fence (std::memory_order_release);

          initialized = true;
        });
    }
    return value;
  }

  bool
  idempotentlyGetNamedEnvironmentVariableAsBool (const char name[],
                                                 std::once_flag& once_flag,
                                                 bool& initialized,
                                                 const char environmentVariableName[],
                                                 const bool defaultValue)
  {
    using BehaviorDetails::namedVariableMap_;
    // The extra "initialized" check avoids the cost of synchronizing
    // on the std::call_once for every call to this function.  We want
    // it to be cheap to fill the namedVariableMap_, so that users aren't
    // tempted to try to cache it themselves.
    if (! initialized) {
      std::call_once (once_flag, [&] () {
          setEnvironmentVariableMap (environmentVariableName,
                                     namedVariableMap_,
                                     defaultValue);
          // http://preshing.com/20130922/acquire-and-release-fences/
          //
          // "A release fence prevents the memory reordering of any
          // read or write which precedes it in program order with any
          // write which follows it in program order."
          //
          // The point is to prevent the assignment to 'namedValueMap' from
          // getting reordered after the assignment to 'initialized' (the
          // so-called "StoreStore" reordering).  That would be bad in this
          // case, because then other threads might read 'initialized' as true,
          // yet would fail to pick up the change to 'namedValueMap'.
          //
          // It's harmless if other threads don't see the write to
          // 'initialized', but did see the write to 'namedValueMap'.  In that
          // case, they would just attempt and fail to enter the std::call_once,
          // and return (the correct value of) 'namedValueMap'.
          std::atomic_thread_fence (std::memory_order_release);

          initialized = true;
        });
    }
    auto thisEnvironmentVariableMap = namedVariableMap_[environmentVariableName];
    auto thisEnvironmentVariable = thisEnvironmentVariableMap.find(name);
    if (thisEnvironmentVariable != thisEnvironmentVariableMap.end())
      return thisEnvironmentVariable->second;
    return thisEnvironmentVariableMap["DEFAULT"];
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

  constexpr bool assumeMpiIsCudaAwareDefault () {
#ifdef TPETRA_ASSUME_CUDA_AWARE_MPI
    return true;
#else
    return false;
#endif // TPETRA_ASSUME_CUDA_AWARE_MPI
  }   

} // namespace (anonymous)

bool Behavior::debug ()
{
  constexpr char envVarName[] = "TPETRA_DEBUG";
  constexpr bool defaultValue = debugDefault ();

  static std::once_flag flag_;
  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariableAsBool (flag_,
                                                   value_,
                                                   initialized_,
                                                   envVarName,
                                                   defaultValue);
}

bool Behavior::verbose ()
{
  if (BehaviorDetails::verboseDisabled_) return false;

  constexpr char envVarName[] = "TPETRA_VERBOSE";
  constexpr bool defaultValue = verboseDefault ();

  static std::once_flag flag_;
  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariableAsBool (flag_,
                                                   value_,
                                                   initialized_,
                                                   envVarName,
                                                   defaultValue);
}

bool Behavior::assumeMpiIsCudaAware ()
{
  constexpr char envVarName[] = "TPETRA_ASSUME_CUDA_AWARE_MPI";
  constexpr bool defaultValue = assumeMpiIsCudaAwareDefault ();

  static std::once_flag flag_;
  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return idempotentlyGetEnvironmentVariableAsBool (flag_,
                                                   value_,
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
    if (varVal == NULL) {
        savedval = 3000; 
        return savedval; 
    }
    savedval = std::stoi(std::string(varVal));
    return savedval;
}


bool Behavior::debug (const char name[])
{
  constexpr char envVarName[] = "TPETRA_DEBUG";
  constexpr bool defaultValue = false;

  static std::once_flag flag_;
  static bool initialized_ = false;
  return idempotentlyGetNamedEnvironmentVariableAsBool (name,
                                                        flag_,
                                                        initialized_,
                                                        envVarName,
                                                        defaultValue);
}

bool Behavior::verbose (const char name[])
{
  if (BehaviorDetails::verboseDisabled_) return false;

  constexpr char envVarName[] = "TPETRA_VERBOSE";
  constexpr bool defaultValue = false;

  static std::once_flag flag_;
  static bool initialized_ = false;
  return idempotentlyGetNamedEnvironmentVariableAsBool (name,
                                                        flag_,
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

} // namespace Details
} // namespace Tpetra

