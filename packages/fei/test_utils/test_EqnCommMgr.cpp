/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/


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
