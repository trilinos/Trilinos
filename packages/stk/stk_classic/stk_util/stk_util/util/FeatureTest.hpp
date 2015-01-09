#ifndef STK_UTIL_UTIL_FeatureTest_h
#  define STK_UTIL_UTIL_FeatureTest_h

// #include <boost/config.hpp>

///
/// @addtogroup FeatureTestDetail
/// @{
///
/// The following feature test and bug macros been defined to conditionally compile in
/// specific features and compile workarounds for known bugs.
///
/// Bugs and workarounds:
///
/// @def SIERRA_TYPE_INFO_BEFORE_EQUALITY_BUG
/// SIERRA_TYPE_INFO_BEFORE_EQUALITY_BUG -- The before() function of the std::type_info
///   class has been implemented improperly, resulting in type_info which are equal to
///   always return true.  Equal std::type_info objects should alway return false.  Define
///   this is type_info::before() is not implemented properly.
///
/// @def SIERRA_AUTO_PTR_ASSIGNMENT_BUG
/// SIERRA_AUTO_PTR_ASSIGNMENT_BUG -- The assignment operator of an std::auto_ptr object
///   does not properly pass the ownership token to the lhs object.  This results in the
///   object being deleted, resulting on heap corruption.  The reset() function should be
///   used rather that assignement.  Define this if the auto_ptr assign function does not
///   work properly.
///
/// @def SIERRA_TEMPLATE_CALL_BUG
/// SIERRA_TEMPLATE_CALL_BUG -- The calling of a template member function erroneously
///   requires a "template" specified on GNU 3.2 and 3.3 compilers which is invalid for
///   other compilers.  Define this use the erroneous implementation.
///
/// @def SIERRA_TEMPLATE_FUNCTION_SELECT_BUG
/// SIERRA_TEMPLATE_FUNCTION_SELECT_BUG -- The function selection algorithm is broken and
///   special workaround code must included to allow the compile to complete.
///
/// @def SIERRA_SRAND_PARALLEL_IO_BUG
/// SIERRA_SRAND_PARALLEL_IO_BUG -- There is some kind of bug in the Redstorm I/O subsystem
///   that reveals itself when the random sequence is seeded.  Define this to disable the
///   seeding of the random number generator.
///
/// @def SIERRA_MPI_ALLREDUCE_USER_FUNCTION_BUG
/// SIERRA_MPI_ALLREDUCE_USER_FUNCTION_BUG -- The platforms MPI_Allreduce function does
///   not behave properly when using a user defined function.
///
/// @def SIERRA_SETVBUF_OUTPUT
/// SIERRA_SETVBUF_OUTPUT -- The output buffering on Redstorm significantly degrades
///   performance due to buffering issues.  This macro causes the opening of log files to
///   set the size of the C++ output buffer to be large.  Define this if you want the large
///   log output buffers defined.
///
/// @def SIERRA_DIAG_ENDL_NOFLUSH
/// SIERRA_DIAG_ENDL_NOFLUSH -- The output buffering on Redstorm significantly degrades
///   performance due to buffering issues.  This macro causes the diagnostic writer to not
///   flush the output buffer when writing newlines.  Define this if you want output buffer
///   flushing to not occur on every diagnostic output line.
///
/// @def SIERRA_USE_PLATFORM_DEMANGLER
/// SIERRA_USE_PLATFORM_DEMANGLER -- The platform type_info::name() function returns a
///   mangled name which needs to be demangled.  Implement the platform specific function
///   in Slib_EnvPlatform.C
///
///
/// Features:
///
/// @def SIERRA_DLOPEN_ENABLED
/// SIERRA_DLOPEN_ENABLED -- The dlopen functionality is built into the sierra
///   applications.  However it may not function on all platforms or be desired on all
///   distributions. Define this macro if you wish to implement this functionality.
///
/// @def SIERRA_MEMORY_INFO
/// SIERRA_MEMORY_INFO -- The platform supports memory usage information.
///
/// @def SIERRA_HEAP_INFO
/// SIERRA_HEAP_INFO -- The platform supports heap usage information.
///
/// @def SIERRA_MPI_ABORT_SIGNAL
/// SIERRA_MPI_ABORT_SIGNAL -- The MPI sends this signal to abort the processes.
///
/// @def SIERRA_USER_SHUTDOWN_SIGNAL
/// SIERRA_USER_SHUTDOWN_SIGNAL -- The user sends this signal to tell the application to
///   exit gracefully at it earliest convenience.
///
/// @def SIERRA_SHUTDOWN_SIGNAL
/// SIERRA_SHUTDOWN_SIGNAL -- The platform sends this signal to tell the application to
///   exit gracefully at it earliest convenience.
///
/// @def SIERRA_INCLUDE_LIBPAPI
/// SIERRA_INCLUDE_LIBPAPI -- The Performance API can be utilized for timing and operation
///   counting.  Define this macro if you wish to enable the high accuracy timing and
///   operation counting.  (NOT IMPLEMENTED)
///

// Platform/operating system based features and bugs
#if defined(REDS)				// Redstorm
#  define SIERRA_SETVBUF_OUTPUT 1
#  define SIERRA_DIAG_ENDL_NOFLUSH
#  define SIERRA_SRAND_PARALLEL_IO_BUG
#  define SIERRA_HEAP_INFO
#  define SIERRA_MEMORY_INFO
#  define SIERRA_USER_SHUTDOWN_SIGNAL SIGUSR1
#  define SIERRA_SHUTDOWN_SIGNAL SIGTERM

#elif defined(_CRAYXE)			        // Cray
#  define SIERRA_SETVBUF_OUTPUT 1
#  define SIERRA_DIAG_ENDL_NOFLUSH
#  define SIERRA_SRAND_PARALLEL_IO_BUG
#  define SIERRA_HEAP_INFO
#  define SIERRA_MEMORY_INFO
#  define SIERRA_SHUTDOWN_SIGNAL SIGTERM
#  define SIERRA_USER_SHUTDOWN_SIGNAL SIGURG

#elif defined(__sun)				// sass8000/sass9000
#  define SIERRA_HEAP_INFO
#  define SIERRA_MEMORY_INFO
#  define SIERRA_MPI_ABORT_SIGNAL SIGTERM
#  define SIERRA_USER_SHUTDOWN_SIGNAL SIGHUP

#elif defined(__linux__)	// Generic linux
#  define SIERRA_USE_PLATFORM_DEMANGLER
#  define SIERRA_HEAP_INFO
#  define SIERRA_MEMORY_INFO
#  define SIERRA_MPI_ABORT_SIGNAL SIGTERM
#  define SIERRA_USER_SHUTDOWN_SIGNAL SIGUSR1

#elif defined(__APPLE__)	// MacOS
#  define SIERRA_USE_PLATFORM_DEMANGLER
#  define SIERRA_HEAP_INFO
#  define SIERRA_MEMORY_INFO
#  define SIERRA_MPI_ABORT_SIGNAL SIGTERM
#  define SIERRA_USER_SHUTDOWN_SIGNAL SIGUSR1

#else // Unknown platform
#  warning Could not determine platform/operating system
#endif


// Compiler/runtime specific features and bugs

#if defined(__xlC__)				// IBM compiler
#  if __xlC__ < 0x0800
#    define SIERRA_TEMPLATE_FUNCTION_SELECT_BUG
#  endif
#  define SIERRA_TYPE_INFO_BEFORE_EQUALITY_BUG
#  define SIERRA_AUTO_PTR_ASSIGNMENT_BUG

#elif defined(__sun) || defined(__SUNPRO_CC)	// Sun spro compiler
#  define SIERRA_MPI_ALLREDUCE_USER_FUNCTION_BUG

#elif defined(__INTEL_COMPILER)			// Intel compiler
#  if __INTEL_COMPILER/100 == 10 && defined(__ia64) // Version 10 Intel compiler on ia64
#    define SIERRA_IA64_OPTIMIZER_FIX
#  elif __INTEL_COMPILER/100 == 11 && defined(__ia64) // Version 11 Intel compiler on ia64
#    define SIERRA_IA64_OPTIMIZER_FIX
#  elif __INTEL_COMPILER/100 >= 12 && defined(__ia64) // Version 12+ Intel compiler on ia64
#    define SIERRA_IA64_OPTIMIZER_WARN
#  endif

#elif defined(__PGI)				// PGI compiler

#elif defined(_CRAYC)				// CRAY compiler

#elif defined(__APPLE_CC__)			// Apple compiler

#elif defined(__PATHSCALE__)			// Pathscale compiler
#  if (__GNUC__ > 3) || defined(PATHSCALE_GNU4) // Only with gcc3 front-end
#    define NO_SIERRA_TEMPLATE_CALL_BUG
#  else
#    define SIERRA_TEMPLATE_CALL_BUG
#  endif

#elif defined(__GNUC__)				// GNU compiler, do this last since *every* compiler thinks its gcc
#  if __GNUC__ == 3 && __GNUC_MINOR__ < 4
#    define SIERRA_TEMPLATE_CALL_BUG
#  endif

#else // Unknown compiler
#  warning Could not determine compiler/runtime
#endif

///
/// @}
///

#endif // STK_UTIL_UTIL_FeatureTest_h
