// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_TestHelpers.hpp
    \brief common code used by tests
*/
#ifndef ZOLTAN2_TESTHELPERS_HPP
#define ZOLTAN2_TESTHELPERS_HPP

#include <Zoltan2_Util.hpp>
#include <iostream>

#include <Tpetra_Map.hpp>
typedef Tpetra::Map<>::node_type znode_t;

// The path to the directory of test data

#define STR_VALUE(path) #path
#define PATH_NAME(path) STR_VALUE(path)

#ifdef Z2_DATA_DIR
  std::string testDataFilePath(PATH_NAME(Z2_DATA_DIR));
#else
  std::string testDataFilePath(".");
#endif

// The path to the Zoltan1 test directory.  We use
// some of their data for testing.

#ifdef Z1_TEST_DIR
  std::string zoltanTestDirectory(PATH_NAME(Z1_TEST_DIR));
#else
  std::string zoltanTestDirectory(".");
#endif

//////////////////////////////////////////////////////////////////////////
//
// If Tpetra is compiled with explicit instantiation,
// we have to choose data types that are compiled into Tpetra.
//
// Epetra uses (scalar/lno/gno) == (double/int/int) data types.  If we
// are using these data types, we can test Epetra user input.

// TODO:  KDD 8/13/14
// Global definitions of types gno_t, lno_t, zgid_t and
// scalar_t can cause bugs in the code.  If a class fails to define these
// types, but this file is included before the class file, the types
// from Zoltan2_TestHelpers.hpp will be used in the class.  Compilation in
// user programs (without Zoltan2_TestHelpers.hpp) would then fail.  An
// example of this bug was in the GeometricGenerator class, which used
// scalar_t without defining it.
// In this "fix," I changed gno_t, lno_t, zgid_t, scalar_t, and node_t to
// zgno_t, zlno_t, zzgid_t, zscalar_t and znode_t in Zoltan2_TestHelpers.hpp.
// This change is not the best fix; a better fix would remove the global
// definitions, but that would require more work.  (An even better change
// would use the Teuchos test framework to cycle through various options,
// but that would require even more work and should be delayed until we
// revamp the testing.)

#include <TpetraCore_config.h>

typedef int zpart_t; // added this for consistency but needs further discussion

typedef Tpetra::Map<>::local_ordinal_type zlno_t;
typedef Tpetra::Map<>::global_ordinal_type zgno_t;

#ifdef HAVE_TPETRA_DOUBLE
  typedef double zscalar_t;
# define HAVE_EPETRA_SCALAR_TYPE
#else
  typedef float zscalar_t;
#endif

#if defined HAVE_TPETRA_INT_INT
#  if defined HAVE_EPETRA_SCALAR_TYPE
#    define HAVE_EPETRA_DATA_TYPES
#  endif
#endif

#ifndef HAVE_ZOLTAN2_EPETRA
#  undef HAVE_EPETRA_SCALAR_TYPE
#  undef HAVE_EPETRA_DATA_TYPES
#endif

//////////////////////////////////////////////////////////////////////////

#define MEMORY_CHECK(iPrint, msg) \
  if (iPrint){ \
    long kb = Zoltan2::getProcessKilobytes(); \
    std::cout.width(10); \
    std::cout.fill('*'); \
    std::cout << kb << " KB, " << msg << std::endl; \
    std::cout.width(0); \
    std::cout.fill(' '); \
  }

#include <ErrorHandlingForTests.hpp>
#include <UserInputForTests.hpp>
#include <PrintData.hpp>

#endif
