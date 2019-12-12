#include "Tpetra_Details_Behavior.hpp"
#include "TpetraCore_config.h"
#include <algorithm> // std::transform
#include <cstdlib> // std::getenv
#include <cctype> // std::toupper
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
    // only call getenv once, save the value.
    static int savedval=-1;
    if(savedval!=-1) return static_cast<size_t>(savedval);
    const char* varVal = std::getenv ("TPETRA_VERBOSE_PRINT_COUNT_THRESHOLD");
    if (varVal == nullptr) {
        savedval = 200; 
        return static_cast<size_t>(savedval); 
    }
    savedval = std::stoi(std::string(varVal));
    return static_cast<size_t>(savedval);
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

} // namespace Details
} // namespace Tpetra

