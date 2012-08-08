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

typedef Kokkos::DefaultNode::DefaultNodeType node_t;

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

//
// If Tpetra is compiled with explicit instantiation,
// then we will use these types in our tests.
//
// Epetra uses int, int, double data types.  If we
// are using these data types, then we can test 
// cases of Epetra user input.
//


#if defined HAVE_ZOLTAN2_INST_FLOAT_INT_LONG

typedef int lno_t;
typedef long gno_t;
typedef float scalar_t;

#elif defined HAVE_ZOLTAN2_INST_DOUBLE_INT_LONG

typedef int lno_t;
typedef long gno_t;
typedef double scalar_t;

#elif defined HAVE_ZOLTAN2_INST_FLOAT_INT_INT

typedef int lno_t;
typedef int gno_t;
typedef float scalar_t;

#elif defined HAVE_ZOLTAN2_INST_DOUBLE_INT_INT

typedef int lno_t;
typedef int gno_t;
typedef double scalar_t;
#define HAVE_EPETRA_DATA_TYPES

#else

typedef int lno_t;
typedef int gno_t;
typedef double scalar_t;
#define HAVE_EPETRA_DATA_TYPES

#endif

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
