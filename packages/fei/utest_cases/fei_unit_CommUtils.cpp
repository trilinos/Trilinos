
#include <fei_iostream.hpp>
#include <fei_Exception.hpp>
#include "fei_TemplateUtils.hpp"

#include "snl_fei_CommUtils.hpp"
#include "feiArray.hpp"

#include <fei_unit_CommUtils.hpp>

#include <vector>
#include <cmath>

#undef fei_file
#define fei_file "fei_unit_CommUtils.cpp"
#include "fei_ErrMacros.hpp"

void test_fei_Allgatherv(MPI_Comm comm)
{
  FEI_COUT << "testing fei::Allgatherv...";

  std::vector<double> send(8, 1.0), recv;
  std::vector<double> send2(8, 1.0), recv2;
  std::vector<int> recvLengths;
  std::vector<int> recvLengths2;

  snl_fei::CommUtils<double> commUtils(comm);

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

  if (fei::Allgatherv<double>(comm, send2, recvLengths2, recv2) != 0) {
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

int test_CommUtils_test1(MPI_Comm comm)
{
  FEI_COUT << "testing CommUtils.Allgatherv...";

  std::vector<int> send(8, 1), recv;
  std::vector<int> recvLengths;

  snl_fei::CommUtils<int> commUtils(comm);

  CHK_ERR( commUtils.Allgatherv(send, recvLengths, recv) );

  if ((int)recvLengths.size() != commUtils.numProcs()) ERReturn(-1);
  if ((int)recv.size() != 8*commUtils.numProcs()) ERReturn(-1);

  for(unsigned i=0; i<recv.size(); i++) {
    if (recv[i] != 1) ERReturn(-1);
  }

  //use a zero-length send-buffer on odd-numbered processors
  if (commUtils.localProc()%2 != 0) send.resize(0);

  CHK_ERR( commUtils.Allgatherv(send, recvLengths, recv) );

  int expectedLength = 0;
  for(int p=0; p<commUtils.numProcs(); p++) {
    if (p%2 == 0) expectedLength += 8;
  }

  if ((int)recvLengths.size() != commUtils.numProcs()) ERReturn(-1);
  if ((int)recv.size() != expectedLength) ERReturn(-1);

  for(unsigned j=0; j<recv.size(); j++) {
    if (recv[j] != 1) ERReturn(-1);
  }

  std::vector<int> local(5, commUtils.localProc()), global;

  CHK_ERR( commUtils.GlobalMax(local, global) );

  if (global.size() != local.size()) ERReturn(-1);

  for(unsigned i=0; i<global.size(); i++) {
    if (global[i] != commUtils.numProcs()-1) ERReturn(-1);
  }

  FEI_COUT << "ok" << FEI_ENDL;

  return(0);
}

bool test_CommUtils::run(MPI_Comm comm)
{
  test_fei_Allgatherv(comm);

  if (test_CommUtils_test1(comm) != 0) {
    throw fei::Exception("test_CommUtils_test1 failed.");
  }

  return true;
}

