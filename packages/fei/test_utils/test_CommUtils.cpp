/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_macros.hpp"

#include "test_utils/test_CommUtils.hpp"
#include "fei_TemplateUtils.hpp"

#include "snl_fei_CommUtils.hpp"
#include "feiArray.hpp"
#include <cmath>
#undef fei_file
#define fei_file "test_CommUtils.cpp"
#include "fei_ErrMacros.hpp"

test_CommUtils::test_CommUtils(MPI_Comm comm)
 : tester(comm)
{
}

test_CommUtils::~test_CommUtils()
{
}

void test_fei_Allgatherv(MPI_Comm comm_);

int test_CommUtils::runtests()
{
  test_fei_Allgatherv(comm_);

  CHK_ERR( test1() );
  CHK_ERR( test2() );
  CHK_ERR( test3() );
  CHK_ERR( test4() );

  return(0);
}

int test_CommUtils::test1()
{
  FEI_COUT << "testing CommUtils.Allgatherv...";

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

  FEI_COUT << "ok" << FEI_ENDL;

  return(0);
}

int test_CommUtils::test2()
{
  return(0);
}

void test_fei_Allgatherv(MPI_Comm comm_)
{
  FEI_COUT << "testing fei::Allgatherv...";

  std::vector<double> send(8, 1.0), recv;
  std::vector<double> send2(8, 1.0), recv2;
  std::vector<int> recvLengths;
  std::vector<int> recvLengths2;

  snl_fei::CommUtils<double> commUtils(comm_);

  if (commUtils.Allgatherv(send, recvLengths, recv) != 0) {
    throw fei::Exception("commUtils.Allgatherv test failed 1.");
  }

  if ((int)recvLengths.size() != commUtils.numProcs()) {
    throw fei::Exception("commUtils.Allgatherv test failed 2.");
  }
  if ((int)recv.size() != 8*commUtils.numProcs()) {
    throw fei::Exception("commUtils.Allgatherv test failed 3.");
  }

  for(unsigned i=0; i<recv.size(); i++) {
    if (std::abs(recv[i] - 1.0) > 1.e-49) {
      throw fei::Exception("commUtils.Allgatherv test failed 4.");
    }
  }

  if (fei::Allgatherv<double>(comm_, send2, recvLengths2, recv2) != 0) {
    throw fei::Exception("fei::Allgatherv test failed 1.");
  }

  if (recvLengths2 != recvLengths || recv2 != recv) {
    throw fei::Exception("fei::Allgatherv test failed 2.");
  }

  //use a zero-length send-buffer on odd-numbered processors
  if (commUtils.localProc()%2 != 0) send.resize(0);

  if (commUtils.Allgatherv(send, recvLengths, recv) != 0) {
    throw fei::Exception("commUtils.Allgatherv test failed 5.");
  }

  int expectedLength = 0;
  for(int p=0; p<commUtils.numProcs(); p++) {
    if (p%2 == 0) expectedLength += 8;
  }

  if ((int)recvLengths.size() != commUtils.numProcs()) {
    throw fei::Exception("commUtils.Allgatherv test failed 6.");
  }
  if ((int)recv.size() != expectedLength) {
    throw fei::Exception("commUtils.Allgatherv test failed 7.");
  }

  for(unsigned j=0; j<recv.size(); j++) {
    if (std::abs(recv[j] - 1.0) > 1.e-49) {
      throw fei::Exception("commUtils.Allgatherv test failed 8.");
    }
  }

  std::vector<double> local(5, commUtils.localProc()), global;

  if (commUtils.GlobalMax(local, global) != 0) {
    throw fei::Exception("commUtils.Allgatherv test failed 9.");
  }

  if (global.size() != local.size()) {
    throw fei::Exception("commUtils.Allgatherv test failed 10.");
  }

  for(unsigned i=0; i<global.size(); i++) {
    if (std::abs(global[i] - commUtils.numProcs()+1) > 1.e-49) {
      throw fei::Exception("commUtils.Allgatherv test failed 11.");
    }
  }

  FEI_COUT << "ok" << FEI_ENDL;
}

int test_CommUtils::test3()
{
  return(0);
}

int test_CommUtils::test4()
{
  return(0);
}
