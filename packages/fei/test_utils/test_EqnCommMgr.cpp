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

#include <fei_macros.hpp>

#include <test_utils/test_EqnCommMgr.hpp>
#include <fei_CommUtils.hpp>
#include <fei_defs.h>

#include <fei_ProcEqns.hpp>
#include <fei_EqnBuffer.hpp>
#include <fei_EqnCommMgr.hpp>

#undef fei_file
#define fei_file "test_EqnCommMgr.cpp"
#include <fei_ErrMacros.hpp>

test_EqnCommMgr::test_EqnCommMgr(MPI_Comm comm)
 : tester(comm)
{
}

test_EqnCommMgr::~test_EqnCommMgr()
{
}

int test_EqnCommMgr::runtests()
{
  CHK_ERR( test1() );
  CHK_ERR( test2() );
  CHK_ERR( test3() );
  CHK_ERR( test4() );
  return(0);
}

int test_EqnCommMgr::test1()
{
  EqnCommMgr* eqnCommMgr = new EqnCommMgr(comm_);

  int numProcs = fei::numProcs(comm_);
  int localProc = fei::localProc(comm_);

//  int numGlobalEqns = numProcs*5;
  int numLocalEqns = 5;
  int firstLocalEqn = localProc*numLocalEqns;
//  int lastLocalEqn = (localProc+1)*numLocalEqns - 1;

  if (numProcs > 1) {
    for(int p=0; p<numProcs; p++) {
      if (p == localProc) continue;

      for(int i=0; i<numLocalEqns; i++) {
	if (p != 2) eqnCommMgr->addLocalEqn(firstLocalEqn+i, p);
      }
    }
  }

  eqnCommMgr->setNumRHSs(1);

  int p;
  for(p=0; p<numProcs; p++) {
    if (p == localProc) continue;

    for(int i=0; i<numLocalEqns; i++) {
      int eqn = p*numLocalEqns + i;

      eqnCommMgr->addRemoteIndices(eqn, p, &eqn, 1);
    }
  }

  CHK_ERR( eqnCommMgr->exchangeIndices() );

  double zero = 0.0;
  for(p=0; p<numProcs; p++) {
    if (p == localProc) continue;

    for(int i=0; i<numLocalEqns; i++) {
      int eqn = p*numLocalEqns + i;

      eqnCommMgr->addSolnValues(&eqn, &zero, 1);
    }
  }

  EqnCommMgr* eCopy = eqnCommMgr->deepCopy();

  std::vector<int>& localEqns = eqnCommMgr->localEqnNumbers();
  std::vector<int>& localEqnsCopy = eCopy->localEqnNumbers();

  if (localEqns != localEqnsCopy) {
    ERReturn(-1);
  }

  eqnCommMgr->exchangeSoln();

  eqnCommMgr->resetCoefs();

  delete eqnCommMgr;
  delete eCopy;

  return(0);
}

int test_EqnCommMgr::test2()
{
  FEI_COUT << "testing ProcEqns...";

  ProcEqns procEqns;

  procEqns.addEqn(0, localProc_);
  procEqns.addEqn(1, localProc_);
  procEqns.addEqn(2, localProc_);

  procEqns.addEqn(3, 2, localProc_+1);
  procEqns.addEqn(4, 2, localProc_+1);
  procEqns.addEqn(5, 2, localProc_+1);

  ProcEqns* pCopy = procEqns.deepCopy();

  std::vector<int>& eqnsPerProc = procEqns.eqnsPerProcPtr();
  std::vector<int>& eqnsPerProcCopy = pCopy->eqnsPerProcPtr();

  if (eqnsPerProc != eqnsPerProcCopy) {
    ERReturn(-1);
  }

  delete pCopy;

  FEI_COUT << FEI_ENDL;
  return(0);
}

int test_EqnCommMgr::test3()
{
  return(0);
}

int test_EqnCommMgr::test4()
{
  return(0);
}
