/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_mpi.h>
#include <test_utils/test_Database.hpp>
#include <snl_fei_Utils.hpp>
#include <snl_fei_MapContig.hpp>

#undef fei_file
#define fei_file "test_Database.cpp"

#include <fei_ErrMacros.hpp>

test_Database::test_Database(MPI_Comm comm)
 : tester(comm)
{
}

test_Database::~test_Database()
{
}

void test_MapContig_1()
{
  FEI_COUT << "testing snl_fei::MapContig...";

  snl_fei::MapContig<int> mc(0, 3);

  std::pair<snl_fei::MapContig<int>::iterator,bool> mpair = mc.insert(std::pair<int,int>(1, 2));

  snl_fei::MapContig<int>::iterator miter = mpair.first;

  if ((*miter).second != 2) {
    throw std::runtime_error("MapContig insert iter test 1 failed.");
  }

  mc.insert(std::pair<int,int>(0,1));

  snl_fei::MapContig<int>::iterator
    m_iter = mc.begin(),
    m_end = mc.end();

  if ((*m_iter).first != 0) {
    throw std::runtime_error("MapContig iter test 1 failed.");
  }

  if ((*m_iter).second != 1) {
    throw std::runtime_error("MapContig iter test 2 failed.");
  }

  ++m_iter;

  if ((*m_iter).first != 1) {
    throw std::runtime_error("MapContig iter test 3 failed.");
  }

  if ((*m_iter).second != 2) {
    throw std::runtime_error("MapContig iter test 4 failed.");
  }

  FEI_COUT << "ok"<<FEI_ENDL;
}

int test_Database::runtests()
{
  if (numProcs_ > 1) return(0);

  test_MapContig_1();

  CHK_ERR( test1() );
  CHK_ERR( test2() );
  CHK_ERR( test3() );
  CHK_ERR( test4() );
  CHK_ERR( test5() );
  CHK_ERR( test6() );
  CHK_ERR( test7() );
  CHK_ERR( test8() );

  return(0);
}

int test_Database::test1()
{
  return(0);
}

int test_Database::test2()
{
  return(0);
}

int test_Database::test3()
{
  return(0);
}

int test_Database::test4()
{

  return(0);
}

int test_Database::test5()
{

  return(0);
}

int test_Database::test6()
{

  return(0);
}

int test_Database::test7()
{

  return(0);
}

int test_Database::test8()
{

  return(0);
}
