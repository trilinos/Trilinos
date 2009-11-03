/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_mpi.h>

#include <test_utils/test_Set.hpp>

#include <snl_fei_Utils.hpp>
#include <fei_ctg_set.hpp>

#undef fei_file
#define fei_file "test_Set.cpp"
#include <fei_ErrMacros.hpp>

test_Set::test_Set(MPI_Comm comm)
 : tester(comm)
{
}

test_Set::~test_Set()
{
}

template<typename SET_TYPE>
void set_test1(SET_TYPE& set_obj)
{
  if (set_obj.size() < 1) {
    typename SET_TYPE::const_iterator
      s_beg = set_obj.begin(),
      s_end = set_obj.end();

    if (s_beg != s_end) {
      throw std::runtime_error("failed test 1");
    }
  }
  else set_obj.clear();

  std::pair<typename SET_TYPE::const_iterator,bool> result = set_obj.insert(5);

  if (!result.second) {
    throw std::runtime_error("failed test 2");
  }

  result = set_obj.insert(4);

  if (!result.second) {
    throw std::runtime_error("failed test 3");
  }

  result = set_obj.insert(7);

  ++(result.first);
  if (result.first != set_obj.end()) {
    throw std::runtime_error("failed test 4");
  }

  result = set_obj.insert(6);

  if (!result.second) {
    throw std::runtime_error("failed test 5");
  }

  ++(result.first);

  if (*(result.first) != 7) {
    throw std::runtime_error("failed test 6");
  }

  ++(result.first);
  if (result.first != set_obj.end()) {
    throw std::runtime_error("failed test 7");
  }

  result = set_obj.insert(2);
  result = set_obj.insert(3);

  ++(result.first);
  if (*(result.first) != 4) {
    throw std::runtime_error("failed test 8");
  }

  SET_TYPE set_copy(set_obj);

  if (set_copy.size() != set_obj.size()) {
    throw std::runtime_error("failed test 9");
  }

  typename SET_TYPE::const_iterator
    s_iter = set_obj.begin(),
    s_end = set_obj.end();

  typename SET_TYPE::const_iterator
    c_iter = set_copy.begin(),
    c_end = set_copy.end();

  for(; s_iter != s_end; ++s_iter) {
    if (*s_iter != *c_iter) {
      throw std::runtime_error("failed test 10");
    }
    ++c_iter;
  }

  if (c_iter != c_end) {
    throw std::runtime_error("failed test 11");
  }
}

template<typename SET_TYPE>
void set_test2(SET_TYPE& set_obj)
{
  if (set_obj.size() < 1) {
    typename SET_TYPE::const_iterator
      s_beg = set_obj.begin(),
      s_end = set_obj.end();

    if (s_beg != s_end) {
      throw std::runtime_error("failed test2 1");
    }
  }
  else set_obj.clear();

  set_obj.insert2(5);
  set_obj.insert2(4);
  set_obj.insert2(7);
  set_obj.insert2(6);
  set_obj.insert2(2);
  set_obj.insert2(3);

  SET_TYPE set_copy(set_obj);

  if (set_copy.size() != set_obj.size()) {
    throw std::runtime_error("failed test2 2");
  }

  typename SET_TYPE::const_iterator
    s_iter = set_obj.begin(),
    s_end = set_obj.end();

  typename SET_TYPE::const_iterator
    c_iter = set_copy.begin(),
    c_end = set_copy.end();

  for(; s_iter != s_end; ++s_iter) {
    if (*s_iter != *c_iter) {
      throw std::runtime_error("failed test2 3");
    }
    ++c_iter;
  }

  if (c_iter != c_end) {
    throw std::runtime_error("failed test2 4");
  }
}

int test_Set::runtests()
{
  if (numProcs_ > 1) return(0);

  CHK_ERR( test1() );
  CHK_ERR( test2() );
  CHK_ERR( test3() );
  CHK_ERR( test4() );
  CHK_ERR( test5() );
  CHK_ERR( test6() );
  CHK_ERR( test7() );
  CHK_ERR( test8() );
  CHK_ERR( test9() );

  return(0);
}

int test_Set::test1()
{
  return(0);
}

int test_Set::test2()
{
  FEI_COUT << "testing fei::ctg_set<int> insert,insert2,find,iterate...";
  fei::ctg_set<int> sset2;

  sset2.insert(5);
  sset2.insert(8);
  sset2.insert(3);
  sset2.insert(0);
  sset2.insert(4);
  sset2.insert(1);
  sset2.insert(2);

  fei::ctg_set<int>::const_iterator
    ss2_iter = sset2.begin(),
    ss2_end = sset2.end();

  int i=0;
  for(; ss2_iter != ss2_end; ++ss2_iter) {
    if (*ss2_iter != i && *ss2_iter != 8) {
      return(-1);
    }
    ++i;
  }

  int size2 = sset2.size();
  if (size2 != 7) {
    return(-1);
  }

  fei::ctg_set<int>::const_iterator iter4 = sset2.find(4);
  if (*iter4 != 4) {
    return(-2);
  }

  ++iter4;
  if (*iter4 != 5) {
    return(-3);
  }

  fei::ctg_set<int>::const_iterator iter8 = sset2.find(8);
  if (*iter8 != 8) {
    return(-4);
  }

  set_test2(sset2);

  fei::ctg_set<int> sset3;

  sset3.insert2(1);
  sset3.insert2(3);
  sset3.insert2(6);
  sset3.insert2(8);
  sset3.insert2(0);
  sset3.insert2(2);
  sset3.insert2(9);
  sset3.insert2(11);
  sset3.insert2(4);
  sset3.insert2(10);

  int size3 = sset3.size();
  if (size3 != 10) {
    return(-1);
  }

  fei::ctg_set<int>::const_iterator ss3_iter4 = sset3.find(4);
  if (*ss3_iter4 != 4) {
    return(-2);
  }

  ++ss3_iter4;
  if (*ss3_iter4 != 6) {
    return(-3);
  }

  fei::ctg_set<int>::const_iterator ss3_iter8 = sset3.find(8);
  if (*ss3_iter8 != 8) {
    return(-4);
  }

  FEI_COUT << "ok"<<FEI_ENDL;
  return(0);
}

int test_Set::test3()
{
  return(0);
}

int test_Set::test4()
{

  return(0);
}

int test_Set::test5()
{
  FEI_COUT << "testing fei::binarySearch(...,start,end,...)...";

  std::vector<int> array;
  for(int i=0; i<10; ++i) array.push_back(i);

  int start = 2;
  int end = 6;
  int insertPoint = -1;
  int offset = fei::binarySearch(9, &array[0], array.size(),
				     start, end, insertPoint);
  if (offset >= 0) {
    return(-1);
  }

  offset = fei::binarySearch(4, &array[0], array.size(),
				 start, end, insertPoint);

  if (offset < 0) {
    return(-1);
  }

  fei::ctg_set<int> sset;
  sset.insert2(1);
  sset.insert2(5);
  sset.insert2(9);
  sset.insert2(0);
  sset.insert2(4);
  sset.insert2(8);

  if (sset.size() != 6) {
    ERReturn(-1);
  }

  if (sset.find(0) == sset.end()) {
    ERReturn(-1);
  }

  FEI_COUT << "ok"<<FEI_ENDL;

  return(0);
}

int test_Set::test6()
{
  return(0);
}

int test_Set::test7()
{
  return(0);
}


int test_Set::test8()
{
  return(0);
}

int test_Set::test9()
{

  return(0);
}
