#include "KokkosCompat_Details_KokkosInit.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Kokkos_Core.hpp"

#include <cstdlib> // std::atexit, setenv
#include <mutex> // std::call_once, std::once_flag
#include <stdexcept>
#include <vector>
#include <functional> // std::ref

namespace KokkosCompat {
namespace Details {

bool
isKokkosInitialized ()
{
  return Kokkos::is_initialized();
}

namespace { // (anonymous)

void
initializeKokkosWithCmdLineArgs (int& numArgs, char* theArgs[])
{
  if (! isKokkosInitialized ()) {
    Kokkos::initialize (numArgs, theArgs);
    // Since Tpetra initialized Kokkos, it can finalize all the
    // execution spaces.
    std::atexit (Kokkos::finalize_all);
  }
}

void
initializeKokkosWithStruct (const Kokkos::InitArguments& args)
{
  if (! isKokkosInitialized ()) {
    Kokkos::initialize (args);
    // Since Tpetra initialized Kokkos, it can finalize all the
    // execution spaces.
    std::atexit (Kokkos::finalize_all);
  }
}

void
initializeKokkosWithNoArgs ()
{
  if (! isKokkosInitialized ()) {
    std::vector<std::string> args =
      Teuchos::GlobalMPISession::getArgv ();
    int narg = static_cast<int> (args.size ());
    if (narg == 0) {
      initializeKokkosWithCmdLineArgs (narg, NULL);
    }
    else {
      std::vector<char*> args_c (narg);
      for (int k = 0; k < narg; ++k) {
        // mfh 25 Oct 2017: I feel a bit uncomfortable about this
        // const_cast, but there is no other way to pass
        // command-line arguments to Kokkos::initialize.
        args_c[k] = const_cast<char*> (args[k].c_str ());
      }
      initializeKokkosWithCmdLineArgs (narg, args_c.data ());
    }
  }
}

} // namespace (anonymous)

std::once_flag initKokkosWithCmdLineArgs_flag;

void
initializeKokkos (int& numArgs, char* theArgs[])
{
  std::call_once (initKokkosWithCmdLineArgs_flag,
                  initializeKokkosWithCmdLineArgs,
                  std::ref (numArgs), theArgs);
}

std::once_flag initKokkosWithStructArgs_flag;

void
initializeKokkos (const Kokkos::InitArguments& theArgs)
{
  std::call_once (initKokkosWithCmdLineArgs_flag,
                  initializeKokkosWithStruct, theArgs);
}

std::once_flag initKokkosWithNoArgs_flag;

void
initializeKokkos ()
{
  std::call_once (initKokkosWithNoArgs_flag,
                  initializeKokkosWithNoArgs);
}

void
setUpEnvironmentForCuda ()
{
#ifdef KOKKOS_ENABLE_CUDA
  // mfh 24 Oct 2017: This code used to call setenv to set the
  // CUDA_LAUNCH_BLOCKING environment variable to 1.  This doesn't
  // necessarily have an effect, since the CUDA run-time system may
  // not read environment variables after the process has started.
  // Thus, we instead read the environment variable and throw if it
  // is not set.
  const char varName[] = "CUDA_LAUNCH_BLOCKING";
  const char* varVal = std::getenv (varName);
  if (varVal == NULL) {
    throw std::runtime_error ("When running Tpetra with CUDA, Tpetra "
                              "requires that the environment variable "
                              "CUDA_LAUNCH_BLOCKING be set (e.g., to 1).");
  }
  //putenv("CUDA_LAUNCH_BLOCKING=1"); // see note below

  // mfh 22 Jul 2016: POSIX prefers use of setenv over putenv:
  //
  // http://pubs.opengroup.org/onlinepubs/009695399/functions/putenv.html
  //
  // There are two problems with putenv:
  //
  // 1. It takes its input as a nonconst pointer.  This causes a
  //    build warning (e.g., "conversion from a string literal to
  //    'char *' is deprecated") when using a string literal, in
  //    this case "CUDA_LAUNCH_BLOCKING=1", as the input argument.
  //
  // 2. putenv does not let us fix the above build warning by
  //    copying into a temporary nonconst char array and passing
  //    that into putenv.  This is because putenv is free to keep
  //    the original input pointer, rather than copying the input
  //    string.  Thus, if anyone tries to read the environment
  //    variable (via getenv) after the temporary array has been
  //    freed, they would be reading invalid memory.  The above
  //    POSIX standard entry alludes to this: "A potential error
  //    is to call putenv() with an automatic variable as the
  //    argument, then return from the calling function while
  //    string is still part of the environment."

#endif // KOKKOS_ENABLE_CUDA
}

/// \brief Get the value of the "Verbose" parameter as a \c bool.
///
/// This method lets the "Verbose" parameter have type either \c int
/// or \c bool, and returns its value as \c bool.  If the "Verbose"
/// parameter does not exist in the given list, return the default
/// value, which is \c false.
bool
getVerboseParameter (const Teuchos::ParameterList& params)
{
  const bool defaultValue = false; // default value of the parameter

  if (params.isParameter ("Verbose")) {
    if (params.isType<bool> ("Verbose")) { // is it a bool?
      return params.get<bool> ("Verbose");
    }
    else if (params.isType<int> ("Verbose")) { // is it an int?
      return params.get<int> ("Verbose");
    }
    // It might be polite to throw at this point with a helpful
    // message explaining that the parameter has the wrong type,
    // but that would change current behavior, so I'll just
    // leave it.
  }
  return defaultValue;
}

Teuchos::ParameterList getDefaultNodeParameters ()
{
  Teuchos::ParameterList params;
  params.set ("Verbose", 0);
  // -1 says "Let Kokkos pick"
  params.set ("Num Threads", -1);
  params.set ("Num NUMA", -1);
  params.set ("Num CoresPerNUMA", -1);
  params.set ("Device", 0);
  return params;
}

void
getNodeParameters (int& curNumThreads,
                   int& curNumNUMA,
                   int& curNumCoresPerNUMA,
                   int& curDevice,
                   bool& verbose,
                   const Teuchos::ParameterList& params)
{
  curNumThreads = -1; // -1 says "let Kokkos pick"
  if (params.isType<int> ("Num Threads")) {
    curNumThreads = params.get<int> ("Num Threads");
  }
  curNumNUMA = -1; // -1 says "let Kokkos pick"
  if (params.isType<int> ("Num NUMA")) {
    curNumNUMA = params.get<int> ("Num NUMA");
  }
  curNumCoresPerNUMA = -1; // -1 says "let Kokkos pick"
  if (params.isType<int> ("Num CoresPerNUMA")) {
    curNumCoresPerNUMA = params.get<int> ("Num CoresPerNUMA");
  }
  curDevice = 0; // -1 does NOT say "let Kokkos pick" for Cuda Devices
  if (params.isType<int> ("Device")) {
    curDevice = params.get<int> ("Device");
  }
  verbose = getVerboseParameter (params);
}

} // namespace Details
} // namespace KokkosCompat
