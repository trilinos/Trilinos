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
// Passes various zgno_t types to ASSIGN.
// Some combinations should work without error; 
// for these, this test FAILS if ASSIGN throws an error.
// Some combinations should throw an error; 
// for these, this test says it is GOOD if ASSIGN throws an error.

#include <Teuchos_RCP.hpp>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_TPLTraits.hpp>

#define MAX(a,b) ((a) > (b) ? (a) : (b))

#ifdef HAVE_ZOLTAN2_SCOTCH
// stdint.h for int64_t in scotch header
#include <stdint.h>
extern "C"{
#ifndef HAVE_ZOLTAN2_MPI
#include "scotch.h"
#else
#include "ptscotch.h"
#endif
}
#endif   // HAVE_ZOLTAN2_SCOTCH

#ifdef HAVE_ZOLTAN2_PARMETIS

extern "C"{
#include "parmetis.h"
}

#define PARMETIS_IS_OK 1

#endif  // HAVE_ZOLTAN2_MPI


#define PRINTMSG(s) \
  std::cout << (s) << " " << __FILE__ << ":" << __LINE__ << std::endl

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);

  int ierr = 0;

  ///////////////////////////////////////////////////////////
  // Test conversions into integers

  // Assignments that should always work 
  // (since the zgno value fits in an integer)
  int intIdx;
  try {
    int zgno = 123;
    Zoltan2::TPL_Traits<int,int>::ASSIGN(intIdx, zgno);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: int to int");
    ierr++;
  }

  try {
    unsigned int zgno = 123;
    Zoltan2::TPL_Traits<int,unsigned int>::ASSIGN(intIdx, zgno);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: unsigned int to int");
    ierr++;
  }

  try {
    long zgno = 123;
    Zoltan2::TPL_Traits<int,long>::ASSIGN(intIdx, zgno);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: long to int");
    ierr++;
  }
 
  try {
    size_t zgno = 123;
    Zoltan2::TPL_Traits<int,size_t>::ASSIGN(intIdx, zgno);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: size_t to int");
    ierr++;
  }

  // Assignments that should not work
  try {
    long long zgno = (long long)1 << 40;
    Zoltan2::TPL_Traits<int,long long>::ASSIGN(intIdx, zgno);
  }
  catch (std::exception &e) {
    PRINTMSG("GOOD: big long long to int throws exception");
  }

  try {
    size_t zgno = (size_t)1 << 40;
    Zoltan2::TPL_Traits<int,size_t>::ASSIGN(intIdx, zgno);
  }
  catch (std::exception &e) {
    PRINTMSG("GOOD: big size_t to int throws exception");
  }

  try {
    unsigned zgno = (1 << 31) + 1;
    Zoltan2::TPL_Traits<int,unsigned>::ASSIGN(intIdx, zgno);
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
    Zoltan2::TPL_Traits<size_t,long long>::ASSIGN(sizetIdx, zgno);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: big long long to size_t");
    ierr++;
  }
 
  try {
    size_t zgno = (size_t)1 << 40;
    Zoltan2::TPL_Traits<size_t,size_t>::ASSIGN(sizetIdx, zgno);
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
    Zoltan2::TPL_Traits<int64_t,long long>::ASSIGN(int64Idx, zgno);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: big long long to int64_t");
    ierr++;
  }
 
  try {
    size_t zgno = (size_t)1 << 40;
    Zoltan2::TPL_Traits<int64_t,size_t>::ASSIGN(int64Idx, zgno);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: big size_t to int64_t");
    ierr++;
  }

  // Assignments that should not work
  try {
    size_t zgno = ((size_t)1 << 63) + 1 ;
    Zoltan2::TPL_Traits<int64_t,size_t>::ASSIGN(int64Idx, zgno);
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
    Zoltan2::TPL_Traits<SCOTCH_Num,int>::ASSIGN(scotchIdx, zgno);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: int to SCOTCH_Num");
    ierr++;
  }

  try {
    unsigned int zgno = 123;
    Zoltan2::TPL_Traits<SCOTCH_Num,unsigned int>::ASSIGN(scotchIdx, zgno);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: unsigned int to SCOTCH_Num");
    ierr++;
  }

  try {
    long zgno = 123;
    Zoltan2::TPL_Traits<SCOTCH_Num,long>::ASSIGN(scotchIdx, zgno);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: long to SCOTCH_Num");
    ierr++;
  }
 
  try {
    size_t zgno = 123;
    Zoltan2::TPL_Traits<SCOTCH_Num,size_t>::ASSIGN(scotchIdx, zgno);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: size_t to SCOTCH_Num");
    ierr++;
  }

  if (sizeof(SCOTCH_Num) == 8) {

    try {
      long long zgno = (long long)1 << 40;
      Zoltan2::TPL_Traits<SCOTCH_Num,long long>::ASSIGN(scotchIdx, zgno);
    }
    catch (std::exception &e) {
      PRINTMSG("FAIL: big unsigned int to SCOTCH_Num");
      ierr++;
    }

    try {
      size_t zgno = (size_t)1 << 40;
      Zoltan2::TPL_Traits<SCOTCH_Num,size_t>::ASSIGN(scotchIdx, zgno);
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
      Zoltan2::TPL_Traits<SCOTCH_Num,long long>::ASSIGN(scotchIdx, zgno);
    }
    catch (std::exception &e) {
      PRINTMSG("GOOD: big long long to 4-byte SCOTCH_Num throws exception");
    }

    try {
      size_t zgno = (size_t)1 << 40;
      Zoltan2::TPL_Traits<SCOTCH_Num,size_t>::ASSIGN(scotchIdx, zgno);
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
    Zoltan2::TPL_Traits<idx_t,int>::ASSIGN(parmetisIdx, zgno);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: int to ParMETIS' idx_t");
    ierr++;
  }

  try {
    unsigned int zgno = 123;
    Zoltan2::TPL_Traits<idx_t,unsigned int>::ASSIGN(parmetisIdx, zgno);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: unsigned int to ParMETIS' idx_t");
    ierr++;
  }

  try {
    long zgno = 123;
    Zoltan2::TPL_Traits<idx_t,long>::ASSIGN(parmetisIdx, zgno);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: long to ParMETIS' idx_t");
    ierr++;
  }
 
  try {
    size_t zgno = 123;
    Zoltan2::TPL_Traits<idx_t,size_t>::ASSIGN(parmetisIdx, zgno);
  }
  catch (std::exception &e) {
    PRINTMSG("FAIL: size_t to ParMETIS' idx_t");
    ierr++;
  }

  if (sizeof(idx_t) == 8) {

    try {
      long long zgno = (long long)1 << 40;
      Zoltan2::TPL_Traits<idx_t,long long>::ASSIGN(parmetisIdx, zgno);
    }
    catch (std::exception &e) {
      PRINTMSG("FAIL: big unsigned int to ParMETIS' idx_t");
      ierr++;
    }

    try {
      size_t zgno = (size_t)1 << 40;
      Zoltan2::TPL_Traits<idx_t,size_t>::ASSIGN(parmetisIdx, zgno);
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
      Zoltan2::TPL_Traits<idx_t,long long>::ASSIGN(parmetisIdx, zgno);
    }
    catch (std::exception &e) {
      PRINTMSG("GOOD: big long long to 4-byte ParMETIS' idx_t throws exception");
    }

    try {
      size_t zgno = (size_t)1 << 40;
      Zoltan2::TPL_Traits<idx_t,size_t>::ASSIGN(parmetisIdx, zgno);
    }
    catch (std::exception &e) {
      PRINTMSG("GOOD: big size_t to 4-byte ParMETIS' idx_t throws exception");
    }
  }
#endif

  ///////////////////////////////////////////////////////////
  // Test conversions into and from ZOLTAN_ID_PTR

  ZOLTAN_ID_PTR zoltanGID = new ZOLTAN_ID_TYPE[4];

  {
    typedef char test_t;
    test_t zgno = 'a';
    zoltanGID[0] = 0; zoltanGID[1] = 0; zoltanGID[2] = 0; zoltanGID[3] = 0;

    int num_gid = MAX(1, sizeof(test_t) / sizeof(ZOLTAN_ID_TYPE));
    if (Zoltan2::TPL_Traits<ZOLTAN_ID_PTR,test_t>::NUM_ID != num_gid) {
      PRINTMSG("FAIL: NUM_ID wrong for char");
      ierr++;
    }

    Zoltan2::TPL_Traits<ZOLTAN_ID_PTR,test_t>::ASSIGN(zoltanGID, zgno);
    if (zoltanGID[0] != ZOLTAN_ID_TYPE(zgno) || zoltanGID[1] != 0 || 
        zoltanGID[2] != 0 || zoltanGID[3] != 0) {
      PRINTMSG("FAIL: char to ZOLTAN_ID_PTR");
      ierr++;
    }

    test_t back;
    Zoltan2::TPL_Traits<test_t,ZOLTAN_ID_PTR>::ASSIGN(back, zoltanGID);
    if (back != zgno) {
      PRINTMSG("FAIL: ZOLTAN_ID_PTR to char");
      ierr++;
    }
  }

  {
    typedef short test_t;
    test_t zgno = 63;
    zoltanGID[0] = 0; zoltanGID[1] = 0; zoltanGID[2] = 0; zoltanGID[3] = 0;

    int num_gid = MAX(1, sizeof(test_t) / sizeof(ZOLTAN_ID_TYPE));
    if (Zoltan2::TPL_Traits<ZOLTAN_ID_PTR,test_t>::NUM_ID != num_gid) {
      PRINTMSG("FAIL: NUM_ID wrong for short");
      ierr++;
    }

    Zoltan2::TPL_Traits<ZOLTAN_ID_PTR,test_t>::ASSIGN(zoltanGID, zgno);
    if (zoltanGID[0] != ZOLTAN_ID_TYPE(zgno) || zoltanGID[1] != 0 || 
        zoltanGID[2] != 0 || zoltanGID[3] != 0) {
      PRINTMSG("FAIL: short to ZOLTAN_ID_PTR");
      ierr++;
    }

    test_t back;
    Zoltan2::TPL_Traits<test_t,ZOLTAN_ID_PTR>::ASSIGN(back, zoltanGID);
    if (back != zgno) {
      PRINTMSG("FAIL: ZOLTAN_ID_PTR to short");
      ierr++;
    }
  }

  {
    typedef int test_t;
    test_t zgno = 123;
    zoltanGID[0] = 0; zoltanGID[1] = 0; zoltanGID[2] = 0; zoltanGID[3] = 0;

    int num_gid = MAX(1, sizeof(test_t) / sizeof(ZOLTAN_ID_TYPE));
    if (Zoltan2::TPL_Traits<ZOLTAN_ID_PTR,test_t>::NUM_ID != num_gid) {
      PRINTMSG("FAIL: NUM_ID wrong for int");
      ierr++;
    }

    Zoltan2::TPL_Traits<ZOLTAN_ID_PTR,test_t>::ASSIGN(zoltanGID, zgno);
    if (zoltanGID[0] != ZOLTAN_ID_TYPE(zgno) || zoltanGID[1] != 0 || 
        zoltanGID[2] != 0 || zoltanGID[3] != 0) {
      PRINTMSG("FAIL: int to ZOLTAN_ID_PTR");
      ierr++;
    }

    test_t back;
    Zoltan2::TPL_Traits<test_t,ZOLTAN_ID_PTR>::ASSIGN(back, zoltanGID);
    if (back != zgno) {
      PRINTMSG("FAIL: ZOLTAN_ID_PTR to int");
      ierr++;
    }
  }

  {
    typedef unsigned int test_t;
    test_t zgno = 456;
    zoltanGID[0] = 0; zoltanGID[1] = 0; zoltanGID[2] = 0; zoltanGID[3] = 0;

    int num_gid = MAX(1, sizeof(test_t) / sizeof(ZOLTAN_ID_TYPE));
    if (Zoltan2::TPL_Traits<ZOLTAN_ID_PTR,test_t>::NUM_ID != num_gid) {
      PRINTMSG("FAIL: NUM_ID wrong for unsigned int");
      ierr++;
    }

    Zoltan2::TPL_Traits<ZOLTAN_ID_PTR,test_t>::ASSIGN(zoltanGID, zgno);
    if (zoltanGID[0] != zgno || zoltanGID[1] != 0 || 
        zoltanGID[2] != 0 || zoltanGID[3] != 0) {
      PRINTMSG("FAIL: unsigned int to ZOLTAN_ID_PTR");
      ierr++;
    }

    test_t back;
    Zoltan2::TPL_Traits<test_t,ZOLTAN_ID_PTR>::ASSIGN(back, zoltanGID);
    if (back != zgno) {
      PRINTMSG("FAIL: ZOLTAN_ID_PTR to unsigned int");
      ierr++;
    }
  }

  {
    typedef long long test_t;
    test_t zgno = ((test_t)1 << 34) + (test_t)17;
    zoltanGID[0] = 0; zoltanGID[1] = 0; zoltanGID[2] = 0; zoltanGID[3] = 0;

    int num_gid = MAX(1, sizeof(test_t) / sizeof(ZOLTAN_ID_TYPE));
    if (Zoltan2::TPL_Traits<ZOLTAN_ID_PTR,test_t>::NUM_ID != num_gid) {
      PRINTMSG("FAIL: NUM_ID wrong for long long");
      ierr++;
    }

    Zoltan2::TPL_Traits<ZOLTAN_ID_PTR,test_t>::ASSIGN(zoltanGID, zgno);
    if (sizeof(ZOLTAN_ID_TYPE) == sizeof(unsigned int)) {
      if (zoltanGID[0] != 17 || zoltanGID[1] != 4 || 
          zoltanGID[2] != 0 || zoltanGID[3] != 0) {
        PRINTMSG("FAIL: long long to ZOLTAN_ID_PTR");
        ierr++;
      }
    }
    else if (sizeof(ZOLTAN_ID_TYPE) == sizeof(unsigned long long)) {
      if (test_t(zoltanGID[0]) != zgno || zoltanGID[1] != 0 ||
          zoltanGID[2] != 0 || zoltanGID[3] != 0) {
        PRINTMSG("FAIL: long long to ZOLTAN_ID_PTR");
        ierr++;
      }
    }
    else {
      // should never get here 
      PRINTMSG("FAIL: unknown sizeof(ZOLTAN_ID_TYPE)");
      ierr++;
    }


    test_t back;
    Zoltan2::TPL_Traits<test_t,ZOLTAN_ID_PTR>::ASSIGN(back, zoltanGID);
    if (back != zgno) {
      std::cout << "back " << back << " != zgno " << zgno << std::endl;
      PRINTMSG("FAIL: ZOLTAN_ID_PTR to long long");
      ierr++;
    }
  }

  {
    typedef unsigned long long test_t;
    test_t zgno = ((test_t)1 << 36) + (test_t)25;
    zoltanGID[0] = 0; zoltanGID[1] = 0; zoltanGID[2] = 0; zoltanGID[3] = 0;

    int num_gid = MAX(1, sizeof(test_t) / sizeof(ZOLTAN_ID_TYPE));
    if (Zoltan2::TPL_Traits<ZOLTAN_ID_PTR,test_t>::NUM_ID != num_gid) {
      PRINTMSG("FAIL: NUM_ID wrong for unsigned long long");
      ierr++;
    }

    Zoltan2::TPL_Traits<ZOLTAN_ID_PTR,test_t>::ASSIGN(zoltanGID, zgno);
    if (sizeof(ZOLTAN_ID_TYPE) == sizeof(unsigned int)) {
      if (zoltanGID[0] != 25 || zoltanGID[1] != 16 || 
          zoltanGID[2] != 0 || zoltanGID[3] != 0) {
        PRINTMSG("FAIL: unsigned long long to ZOLTAN_ID_PTR");
        ierr++;
      }
    }
    else if (sizeof(ZOLTAN_ID_TYPE) == sizeof(unsigned long long)) {
      if (zoltanGID[0] != zgno || zoltanGID[1] != 0 ||
          zoltanGID[2] != 0 || zoltanGID[3] != 0) {
        PRINTMSG("FAIL: long long to ZOLTAN_ID_PTR");
        ierr++;
      }
    }
    else {
      // should never get here 
      PRINTMSG("FAIL: unknown sizeof(ZOLTAN_ID_TYPE)");
      ierr++;
    }



    test_t back;
    Zoltan2::TPL_Traits<test_t,ZOLTAN_ID_PTR>::ASSIGN(back, zoltanGID);
    if (back != zgno) {
      std::cout << "back " << back << " != zgno " << zgno << std::endl;
      PRINTMSG("FAIL: ZOLTAN_ID_PTR to unsigned long long");
      ierr++;
    }
  }

  {
    typedef size_t test_t;
    test_t zgno = 0;
    for (size_t i = 0; i < 8*sizeof(test_t); i++) zgno += (test_t)1<<i;
    zoltanGID[0] = 0; zoltanGID[1] = 0; zoltanGID[2] = 0; zoltanGID[3] = 0;

    int num_gid = MAX(1, sizeof(test_t) / sizeof(ZOLTAN_ID_TYPE));
    if (Zoltan2::TPL_Traits<ZOLTAN_ID_PTR,test_t>::NUM_ID != num_gid) {
      PRINTMSG("FAIL: NUM_ID wrong for size_t");
      ierr++;
    }

    Zoltan2::TPL_Traits<ZOLTAN_ID_PTR,test_t>::ASSIGN(zoltanGID, zgno);
    for (int i = 0; i < num_gid; i++)
      if (zoltanGID[i] != std::numeric_limits<ZOLTAN_ID_TYPE>::max()) {
        PRINTMSG("FAIL: size_t to ZOLTAN_ID_PTR");
        ierr++;
      }
    for (int i = num_gid; i < 4; i++)
      if (zoltanGID[i] != 0) {
        PRINTMSG("FAIL: size_t to ZOLTAN_ID_PTR");
        ierr++;
      }

    test_t back;
    Zoltan2::TPL_Traits<test_t,ZOLTAN_ID_PTR>::ASSIGN(back, zoltanGID);
    if (back != zgno) {
      std::cout << "back " << back << " != zgno " << zgno << std::endl;
      PRINTMSG("FAIL: ZOLTAN_ID_PTR to size_t");
      ierr++;
    }
  }
  delete [] zoltanGID;

  ///////////////////////////////////////////////////////////

  if (ierr == 0)
    std::cout << "PASS" << std::endl;
  else
    std::cout << "FAIL" << std::endl;

  return 0;
}

