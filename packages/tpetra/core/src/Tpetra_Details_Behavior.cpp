#include "Tpetra_Details_Behavior.hpp"
#include "TpetraCore_config.h"
#include <algorithm> // std::transform
#include <atomic> // std::atomic_thread_fence, std::memory_order_release
#include <cstdlib> // std::getenv
#include <cctype> // std::toupper
#include <mutex> // std::call_once, std::once_flag
#include <string>

namespace Tpetra {
namespace Details {

namespace { // (anonymous)

  // See example here:
  //
  // http://en.cppreference.com/w/cpp/string/byte/toupper
  std::string stringToUpper (std::string s)
  {
    std::transform (s.begin (), s.end (), s.begin (),
                    [] (unsigned char c) { return std::toupper (c); });
    return s;
  }

  bool
  getEnvironmentVariableAsBool (const char environmentVariableName[],
                                const bool defaultValue)
  {
    const char* varVal = std::getenv (environmentVariableName);

    bool retVal = defaultValue;
    if (varVal != NULL) {
      const std::string varStr (stringToUpper (std::string (varVal)));

      if (varStr == "1" || varStr == "YES" || varStr == "TRUE") {
        retVal = true;
      }
      else if (varStr == "0" || varStr == "NO" || varStr == "FALSE") {
        retVal = false;
      }
      // Otherwise, use the default value.
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

} // namespace Details
} // namespace Tpetra

