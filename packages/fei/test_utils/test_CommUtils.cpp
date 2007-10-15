/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <test_utils/test_CommUtils.hpp>

#include <snl_fei_CommUtils.hpp>
#include <feiArray.hpp>
#include <cmath>
#undef fei_file
#define fei_file "test_CommUtils.cpp"
#include <fei_ErrMacros.hpp>

test_CommUtils::test_CommUtils(MPI_Comm comm)
 : tester(comm)
{
}

test_CommUtils::~test_CommUtils()
{
}

int test_CommUtils::runtests()
{
  CHK_ERR( test1() );
  CHK_ERR( test2() );
  CHK_ERR( test3() );
  CHK_ERR( test4() );
  return(0);
}

int test_CommUtils::test1()
{
  std::vector<int> send(8, 1), recv;
  std::vector<int> recvLengths;

  snl_fei::CommUtils<int> commUtils(comm_);

  CHK_ERR( commUtils.Allgatherv(send, recvLengths, recv) );

  if ((int)recvLengths.size() != numProcs_) ERReturn(-1);
  if ((int)recv.size() != 8*numProcs_) ERReturn(-1);

  for(unsigned i=0; i<recv.size(); i++) {
    if (recv[i] != 1) ERReturn(-1);
  }

  //use a zero-length send-buffer on odd-numbered processors
  if (localProc_%2 != 0) send.resize(0);

  CHK_ERR( commUtils.Allgatherv(send, recvLengths, recv) );

  int expectedLength = 0;
  for(int p=0; p<numProcs_; p++) {
    if (p%2 == 0) expectedLength += 8;
  }

  if ((int)recvLengths.size() != numProcs_) ERReturn(-1);
  if ((int)recv.size() != expectedLength) ERReturn(-1);

  for(unsigned j=0; j<recv.size(); j++) {
    if (recv[j] != 1) ERReturn(-1);
  }

  std::vector<int> local(5, localProc_), global;

  CHK_ERR( commUtils.GlobalMax(local, global) );

  if (global.size() != local.size()) ERReturn(-1);

  for(unsigned i=0; i<global.size(); i++) {
    if (global[i] != numProcs_-1) ERReturn(-1);
  }

  return(0);
}

int test_CommUtils::test2()
{
  std::vector<double> send(8, 1.0), recv;
  std::vector<int> recvLengths;

  snl_fei::CommUtils<double> commUtils(comm_);

  CHK_ERR( commUtils.Allgatherv(send, recvLengths, recv) );

  if ((int)recvLengths.size() != numProcs_) ERReturn(-1);
  if ((int)recv.size() != 8*numProcs_) ERReturn(-1);

  for(unsigned i=0; i<recv.size(); i++) {
    if (std::abs(recv[i] - 1.0) > 1.e-49) ERReturn(-1);
  }

  //use a zero-length send-buffer on odd-numbered processors
  if (localProc_%2 != 0) send.resize(0);

  CHK_ERR( commUtils.Allgatherv(send, recvLengths, recv) );

  int expectedLength = 0;
  for(int p=0; p<numProcs_; p++) {
    if (p%2 == 0) expectedLength += 8;
  }

  if ((int)recvLengths.size() != numProcs_) ERReturn(-1);
  if ((int)recv.size() != expectedLength) ERReturn(-1);

  for(unsigned j=0; j<recv.size(); j++) {
    if (std::abs(recv[j] - 1.0) > 1.e-49) ERReturn(-1);
  }

  std::vector<double> local(5, localProc_), global;

  CHK_ERR( commUtils.GlobalMax(local, global) );

  if (global.size() != local.size()) ERReturn(-1);

  for(unsigned i=0; i<global.size(); i++) {
    if (std::abs(global[i] - numProcs_+1) > 1.e-49) ERReturn(-1);
  }

  return(0);
}

int test_CommUtils::test3()
{
  return(0);
}

int test_CommUtils::test4()
{
  return(0);
}
