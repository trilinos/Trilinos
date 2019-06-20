// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef STK_UTIL_UTIL_FeatureTest_h
#define STK_UTIL_UTIL_FeatureTest_h

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
/// @def SIERRA_TEMPLATE_FUNCTION_SELECT_BUG
/// SIERRA_TEMPLATE_FUNCTION_SELECT_BUG -- The function selection algorithm is broken and
///   special workaround code must included to allow the compile to complete.
///
/// @def SIERRA_USE_PLATFORM_DEMANGLER
/// SIERRA_USE_PLATFORM_DEMANGLER -- The platform type_info::name() function returns a
///   mangled name which needs to be demangled.
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

// Platform/operating system based features and bugs
#if defined(_CRAYXE)			        // Cray
#  define SIERRA_HEAP_INFO
#  define SIERRA_MEMORY_INFO
#  define SIERRA_SHUTDOWN_SIGNAL SIGTERM
#  define SIERRA_USER_SHUTDOWN_SIGNAL SIGURG

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

#elif defined(__GNUC__)				// GNU compiler, do this last since *every* compiler thinks its gcc

#else // Unknown compiler
#  warning Could not determine compiler/runtime
#endif

///
/// @}
///

#endif // STK_UTIL_UTIL_FeatureTest_h
