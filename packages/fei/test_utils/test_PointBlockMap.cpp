/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/

#include <fei_iostream.hpp>

#include <test_utils/test_PointBlockMap.hpp>
#include <snl_fei_PointBlockMap.hpp>

#undef fei_file
#define fei_file "test_PointBlockMap.cpp"
#include <fei_ErrMacros.hpp>

test_PointBlockMap::test_PointBlockMap(MPI_Comm comm)
 : tester(comm)
{
}

test_PointBlockMap::~test_PointBlockMap()
{
}

int test_PointBlockMap::runtests()
{
  CHK_ERR( test1() );
  CHK_ERR( test2() );
  CHK_ERR( test3() );
  CHK_ERR( test4() );
  return(0);
}

int test_PointBlockMap::test1()
{
  snl_fei::PointBlockMap pbMap;

  CHK_ERR( pbMap.setEqn(0, 0) );
  CHK_ERR( pbMap.setEqn(1, 0) );
  CHK_ERR( pbMap.setEqn(2, 0) );
  CHK_ERR( pbMap.setEqn(3, 1) );
  CHK_ERR( pbMap.setEqn(4, 1) );
  CHK_ERR( pbMap.setEqn(5, 2) );

  CHK_ERR( pbMap.setBlkEqnSize(0, 3) );
  CHK_ERR( pbMap.setBlkEqnSize(1, 2) );
  CHK_ERR( pbMap.setBlkEqnSize(2, 1) );

  int blkEqnOffset = pbMap.getBlkEqnOffset(1,3);
  if (blkEqnOffset != 0) {
    ERReturn(-1);
  }

  int size = pbMap.getBlkEqnSize(1);
  if (size != 2) {
    ERReturn(-1);
  }

  int blkSize, ptEqn;
  CHK_ERR( pbMap.getBlkEqnInfo(1, ptEqn, blkSize) );
  if (ptEqn != 3) {
    ERReturn(-1);
  }
  if (blkSize != 2) {
    ERReturn(-1);
  }

  int blkEqn = pbMap.eqnToBlkEqn(4);
  if (blkEqn != 1) {
    ERReturn(-1);
  }

  if (pbMap.isExactlyBlkEqn(4)) {
    ERReturn(-1);
  }

  int blkOffset;
  CHK_ERR( pbMap.getPtEqnInfo(4, blkEqn, blkOffset) );
  if (blkEqn != 1) {
    ERReturn(-1);
  }

  if (blkOffset != 1) {
    ERReturn(-1);
  }

  if (!pbMap.isExactlyBlkEqn(3)) {
    ERReturn(-1);
  }

  int maxBlkSize = pbMap.getMaxBlkEqnSize();
  if (maxBlkSize == 0) {
    ERReturn(-1);
  }

  return(0);
}

int test_PointBlockMap::test2()
{
  snl_fei::PointBlockMap pbMap;

  CHK_ERR( pbMap.setEqn(0, 0) );
  CHK_ERR( pbMap.setEqn(1, 0) );
  CHK_ERR( pbMap.setEqn(2, 0) );
  CHK_ERR( pbMap.setEqn(3, 1) );
  CHK_ERR( pbMap.setEqn(4, 1) );
  CHK_ERR( pbMap.setEqn(5, 2) );

  CHK_ERR( pbMap.setBlkEqnSize(0, 3) );
  CHK_ERR( pbMap.setBlkEqnSize(1, 2) );
  CHK_ERR( pbMap.setBlkEqnSize(2, 1) );

  std::map<int,std::pair<int,int> >* blkEqns = pbMap.getBlkEqns();

  std::map<int,std::pair<int,int> >::const_iterator
    iter = blkEqns->begin(),
    iter_end = blkEqns->end();

  for(; iter != iter_end; ++iter) {
    if ((*iter).first < 0) {
      CHK_ERR( -1 );
    }
  }

  return(0);
}

int test_PointBlockMap::test3()
{
  return(0);
}

int test_PointBlockMap::test4()
{
  return(0);
}
