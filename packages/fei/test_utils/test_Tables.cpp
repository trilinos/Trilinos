/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/


#include <fei_macros.hpp>

#include <test_utils/test_Tables.hpp>

#include <snl_fei_RaggedTable.hpp>
#include <fei_ProcEqns.hpp>
#include <fei_ctg_set.hpp>
#include <map>
#undef fei_file
#define fei_file "test_Tables.cpp"
#include <fei_ErrMacros.hpp>

test_Tables::test_Tables(MPI_Comm comm)
 : tester(comm)
{
}

test_Tables::~test_Tables()
{
}

int test_Tables::runtests()
{
  //This class doesn't have any parallel tests. The Table objects are
  //purely serial objects, so only run tests if numProcs_==1.
  if (numProcs_ > 1) {
    return(0);
  }

  CHK_ERR( test1() );
  CHK_ERR( test2() );
  CHK_ERR( test3() );
  CHK_ERR( test4() );
  return(0);
}

int test_Tables::test1()
{
  int len = 100;

  snl_fei::RaggedTable<std::map<int,fei::ctg_set<int>*>,fei::ctg_set<int> > ordTable(0,1);
  snl_fei::RaggedTable<std::map<int,fei::ctg_set<int>*>,fei::ctg_set<int> > ordTable2(0,1);

  std::vector<int> keys(len), values(len);
  int i;
  for(i=0; i<len; ++i) {
    keys[i] = i;
    values[i] = i;

    ordTable.addIndices(i, 1, &i);
  }

  keys.push_back(len);
  values.push_back(len);
  ++len;

  ordTable.addIndices(len, &keys[0], len, &values[0] );

  bool same = ordTable.equal(ordTable2, true);
  if (same) {
    return(-1);
  }

  for(i=0; i<len; ++i) {
    keys[i] = i;
    values[i] = i;

    ordTable2.addIndices(i, 1, &i);
  }

  same = ordTable.equal(ordTable2, true);

  ordTable2.addIndices(len, &keys[0], len, &values[0] );

  same = ordTable.equal(ordTable2, false);
  if (!same) {
    return(-2);
  }

  return(0);
}

int test_Tables::test2()
{
  return(0);
}

int test_Tables::test3()
{
  int len = 100;

  ProcEqns peqns;

  std::vector<int> keys(len), values(len);
  for(int i=0; i<len; ++i) {
    keys[i] = i;
    values[i] = i;

    peqns.addEqn(i, i);
    peqns.addEqn(i, i, len-i);
  }

  return(0);
}

int test_Tables::test4()
{
  return(0);
}
