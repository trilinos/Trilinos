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


// Unit test for Zoltan2_TPLTraits.hpp
// Passes various zgno_t types to ASSIGN_TPL_T.
// Some combinations should work without error; 
// for these, this test FAILS if ASSIGN_TPL_T throws an error.
// Some combinations should throw an error; 
// for these, this test says it is GOOD if ASSIGN_TPL_T throws an error.

#include <Teuchos_GlobalMPISession.hpp>   
#include <Teuchos_RCP.hpp>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_TPLTraits.hpp>

#ifdef HAVE_ZOLTAN2_SCOTCH
// stdint.h for int64_t in scotch header
#include <stdint.h>
#ifndef HAVE_ZOLTAN2_MPI
#include "scotch.h"
#else
#include "ptscotch.h"
#endif
#endif

#ifdef HAVE_ZOLTAN2_PARMETIS

#ifndef HAVE_ZOLTAN2_MPI
// ParMETIS requires compilation with MPI.  
// If MPI is not available, make compilation fail.
#error "TPL ParMETIS requires compilation with MPI.  Configure with -DTPL_ENABLE_MPI:BOOL=ON or -DZoltan2_ENABLE_ParMETIS:BOOL=OFF"
  
#else

#include "parmetis.h"
#if (PARMETIS_MAJOR_VERSION < 4)
// Zoltan2 requires ParMETIS v4.x.  
// Make compilation fail for earlier versions of ParMETIS.
#error "Specified version of ParMETIS is not compatible with Zoltan2; upgrade to ParMETIS v4 or later, or build Zoltan2 without ParMETIS."
  
#else

// MPI and ParMETIS version requirements are met.  Proceed.
#define PARMETIS_IS_OK 1

#endif  // ParMETIS version check
#endif  // HAVE_ZOLTAN2_MPI
#endif  // HAVE_ZOLTAN2_PARMETIS



#define PRINTMSG(s) \
  std::cout << (s) << " " << __FILE__ << ":" << __LINE__ << std::endl

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment);

  int ierr = 0;

  ///////////////////////////////////////////////////////////
  // Test conversions into integers

  // Assignments that should always work 
  // (since the zgno value fits in an integer)
  int intIdx;
  try {
    int zgno = 123;
    Zoltan2::TPL_Traits<int,int>::ASSIGN_TPL_T(intIdx, zgno, env);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: int to int");
    ierr++;
  }

  try {
    unsigned int zgno = 123;
    Zoltan2::TPL_Traits<int,unsigned int>::ASSIGN_TPL_T(intIdx, zgno, env);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: unsigned int to int");
    ierr++;
  }

  try {
    long zgno = 123;
    Zoltan2::TPL_Traits<int,long>::ASSIGN_TPL_T(intIdx, zgno, env);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: long to int");
    ierr++;
  }
 
  try {
    size_t zgno = 123;
    Zoltan2::TPL_Traits<int,size_t>::ASSIGN_TPL_T(intIdx, zgno, env);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: size_t to int");
    ierr++;
  }

  // Assignments that should not work
  try {
    long long zgno = (long long)1 << 40;
    Zoltan2::TPL_Traits<int,long long>::ASSIGN_TPL_T(intIdx, zgno, env);
  }
  catch (std::exception &e) {
    PRINTMSG("GOOD: big long long to int throws exception");
  }

  try {
    size_t zgno = (size_t)1 << 40;
    Zoltan2::TPL_Traits<int,size_t>::ASSIGN_TPL_T(intIdx, zgno, env);
  }
  catch (std::exception &e) {
    PRINTMSG("GOOD: big size_t to int throws exception");
  }

  try {
    unsigned zgno = (1 << 31) + 1;
    Zoltan2::TPL_Traits<int,unsigned>::ASSIGN_TPL_T(intIdx, zgno, env);
  }
  catch (std::exception &e) {
    PRINTMSG("GOOD: huge unsigned to int throws exception");
  }

  ///////////////////////////////////////////////////////////
  // Test conversions into size_t

  // Assignments that should always work 

  size_t sizetIdx;
  try {
    long long zgno = (long long)1 << 40;
    Zoltan2::TPL_Traits<size_t,long long>::ASSIGN_TPL_T(sizetIdx, zgno, env);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: big long long to size_t");
    ierr++;
  }
 
  try {
    size_t zgno = (size_t)1 << 40;
    Zoltan2::TPL_Traits<size_t,size_t>::ASSIGN_TPL_T(sizetIdx, zgno, env);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: big size_t to size_t");
    ierr++;
  }

  ///////////////////////////////////////////////////////////
  // Test conversions into int64_t

  // Assignments that should always work 

  int64_t int64Idx;
  try {
    long long zgno = (long long)1 << 40;
    Zoltan2::TPL_Traits<int64_t,long long>::ASSIGN_TPL_T(int64Idx, zgno, env);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: big long long to int64_t");
    ierr++;
  }
 
  try {
    size_t zgno = (size_t)1 << 40;
    Zoltan2::TPL_Traits<int64_t,size_t>::ASSIGN_TPL_T(int64Idx, zgno, env);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: big size_t to int64_t");
    ierr++;
  }

  // Assignments that should not work
  try {
    size_t zgno = ((size_t)1 << 63) + 1 ;
    Zoltan2::TPL_Traits<int64_t,size_t>::ASSIGN_TPL_T(int64Idx, zgno, env);
  }
  catch (std::exception &e) {
    PRINTMSG("GOOD: huge size_t to int64_t threw exception");
  }

#ifdef HAVE_ZOLTAN2_SCOTCH
  ///////////////////////////////////////////////////////////
  // Test conversions into SCOTCH_Num

  SCOTCH_Num scotchIdx;

  // Assignments that should always work 
  // (since the zgno value fits in an integer)
  try {
    int zgno = 123;
    Zoltan2::TPL_Traits<SCOTCH_Num,int>::ASSIGN_TPL_T(scotchIdx, zgno, env);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: int to SCOTCH_Num");
    ierr++;
  }

  try {
    unsigned int zgno = 123;
    Zoltan2::TPL_Traits<SCOTCH_Num,unsigned int>::ASSIGN_TPL_T(scotchIdx, zgno, env);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: unsigned int to SCOTCH_Num");
    ierr++;
  }

  try {
    long zgno = 123;
    Zoltan2::TPL_Traits<SCOTCH_Num,long>::ASSIGN_TPL_T(scotchIdx, zgno, env);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: long to SCOTCH_Num");
    ierr++;
  }
 
  try {
    size_t zgno = 123;
    Zoltan2::TPL_Traits<SCOTCH_Num,size_t>::ASSIGN_TPL_T(scotchIdx, zgno, env);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: size_t to SCOTCH_Num");
    ierr++;
  }

  if (sizeof(SCOTCH_Num) == 8) {

    try {
      long long zgno = (long long)1 << 40;
      Zoltan2::TPL_Traits<SCOTCH_Num,long long>::ASSIGN_TPL_T(scotchIdx, zgno, env);
    }
    catch (std::exception &e) {
      PRINTMSG("FAIL: big unsigned int to SCOTCH_Num");
      ierr++;
    }

    try {
      size_t zgno = (size_t)1 << 40;
      Zoltan2::TPL_Traits<SCOTCH_Num,size_t>::ASSIGN_TPL_T(scotchIdx, zgno, env);
    }
    catch (std::exception &e) {
      PRINTMSG("FAIL: big size_t to SCOTCH_Num");
      ierr++;
    }
  }

  // Assignments that should not work
  if (sizeof(SCOTCH_Num) == 4) {
    try {
      long long zgno = (long long)1 << 40;
      Zoltan2::TPL_Traits<SCOTCH_Num,long long>::ASSIGN_TPL_T(scotchIdx, zgno, env);
    }
    catch (std::exception &e) {
      PRINTMSG("GOOD: big long long to 4-byte SCOTCH_Num throws exception");
    }

    try {
      size_t zgno = (size_t)1 << 40;
      Zoltan2::TPL_Traits<SCOTCH_Num,size_t>::ASSIGN_TPL_T(scotchIdx, zgno, env);
    }
    catch (std::exception &e) {
      PRINTMSG("GOOD: big size_t to 4-byte SCOTCH_Num throws exception");
    }
  }

#endif  // HAVE_ZOLTAN2_SCOTCH

#ifdef PARMETIS_IS_OK
  ///////////////////////////////////////////////////////////
  // Test conversions into ParMETIS' idx_t

  idx_t parmetisIdx;

  // Assignments that should always work 
  // (since the zgno value fits in an integer)
  try {
    int zgno = 123;
    Zoltan2::TPL_Traits<idx_t,int>::ASSIGN_TPL_T(parmetisIdx, zgno, env);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: int to ParMETIS' idx_t");
    ierr++;
  }

  try {
    unsigned int zgno = 123;
    Zoltan2::TPL_Traits<idx_t,unsigned int>::ASSIGN_TPL_T(parmetisIdx, zgno, env);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: unsigned int to ParMETIS' idx_t");
    ierr++;
  }

  try {
    long zgno = 123;
    Zoltan2::TPL_Traits<idx_t,long>::ASSIGN_TPL_T(parmetisIdx, zgno, env);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: long to ParMETIS' idx_t");
    ierr++;
  }
 
  try {
    size_t zgno = 123;
    Zoltan2::TPL_Traits<idx_t,size_t>::ASSIGN_TPL_T(parmetisIdx, zgno, env);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: size_t to ParMETIS' idx_t");
    ierr++;
  }

  if (sizeof(idx_t) == 8) {

    try {
      long long zgno = (long long)1 << 40;
      Zoltan2::TPL_Traits<idx_t,long long>::ASSIGN_TPL_T(parmetisIdx, zgno, env);
    }
    catch (std::exception &e) {
      PRINTMSG("FAIL: big unsigned int to ParMETIS' idx_t");
      ierr++;
    }

    try {
      size_t zgno = (size_t)1 << 40;
      Zoltan2::TPL_Traits<idx_t,size_t>::ASSIGN_TPL_T(parmetisIdx, zgno, env);
    }
    catch (std::exception &e) {
      PRINTMSG("FAIL: big size_t to ParMETIS' idx_t");
      ierr++;
    }
  }

  // Assignments that should not work
  if (sizeof(idx_t) == 4) {
    try {
      long long zgno = (long long)1 << 40;
      Zoltan2::TPL_Traits<idx_t,long long>::ASSIGN_TPL_T(parmetisIdx, zgno, env);
    }
    catch (std::exception &e) {
      PRINTMSG("GOOD: big long long to 4-byte ParMETIS' idx_t throws exception");
    }

    try {
      size_t zgno = (size_t)1 << 40;
      Zoltan2::TPL_Traits<idx_t,size_t>::ASSIGN_TPL_T(parmetisIdx, zgno, env);
    }
    catch (std::exception &e) {
      PRINTMSG("GOOD: big size_t to 4-byte ParMETIS' idx_t throws exception");
    }
  }
#endif
  ///////////////////////////////////////////////////////////

  if (ierr == 0)
    std::cout << "PASS" << std::endl;
  else
    std::cout << "FAIL" << std::endl;

  return 0;
}

