
#include <fei_iostream.hpp>
#include "fei_CommUtils.hpp"

#include <fei_unit_CommUtils.hpp>

#include <vector>
#include <cmath>

#undef fei_file
#define fei_file "fei_unit_CommUtils.cpp"
#include "fei_ErrMacros.hpp"

void test_fei_Allgatherv(MPI_Comm comm)
{
  FEI_COUT << "testing fei::Allgatherv...";

  int num_procs = fei::numProcs(comm);
  int local_proc = fei::localProc(comm);

  std::vector<double> send(8, 1.0), recv;
  std::vector<int> recvLengths;
  std::vector<int> recvLengths2;

  if (fei::Allgatherv(comm, send, recvLengths, recv) != 0) {
    throw std::runtime_error("fei::Allgatherv test failed 1.");
  }

  if ((int)recvLengths.size() != num_procs) {
    throw std::runtime_error("fei::Allgatherv test failed 2.");
  }
  if ((int)recv.size() != 8*num_procs) {
    throw std::runtime_error("fei::Allgatherv test failed 3.");
  }

  for(unsigned i=0; i<recv.size(); i++) {
    if (std::abs(recv[i] - 1.0) > 1.e-49) {
      throw std::runtime_error("fei::Allgatherv test failed 4.");
    }
  }

  //use a zero-length send-buffer on odd-numbered processors
  if (local_proc%2 != 0) send.resize(0);

  if (fei::Allgatherv(comm, send, recvLengths, recv) != 0) {
    throw std::runtime_error("fei::Allgatherv test failed 5.");
  }

  int expectedLength = 0;
  for(int p=0; p<num_procs; p++) {
    if (p%2 == 0) expectedLength += 8;
  }

  if ((int)recvLengths.size() != num_procs) {
    throw std::runtime_error("fei::Allgatherv test failed 6.");
  }
  if ((int)recv.size() != expectedLength) {
    throw std::runtime_error("fei::Allgatherv test failed 7.");
  }

  for(unsigned j=0; j<recv.size(); j++) {
    if (std::abs(recv[j] - 1.0) > 1.e-49) {
      throw std::runtime_error("fei::Allgatherv test failed 8.");
    }
  }

  std::vector<double> local(5, fei::localProc(comm)), global;

  if (fei::GlobalMax(comm, local, global) != 0) {
    throw std::runtime_error("fei::Allgatherv test failed 9.");
  }

  if (global.size() != local.size()) {
    throw std::runtime_error("fei::Allgatherv test failed 10.");
  }

  for(unsigned i=0; i<global.size(); i++) {
    if (std::abs(global[i] - fei::numProcs(comm)+1) > 1.e-49) {
      throw std::runtime_error("fei::Allgatherv test failed 11.");
    }
  }

  FEI_COUT << "ok" << FEI_ENDL;
}

int test_CommUtils_test1(MPI_Comm comm)
{
  FEI_COUT << "testing CommUtils.Allgatherv...";

  std::vector<int> send(8, 1), recv;
  std::vector<int> recvLengths;

  CHK_ERR( fei::Allgatherv(comm, send, recvLengths, recv) );

  if ((int)recvLengths.size() != fei::numProcs(comm)) ERReturn(-1);
  if ((int)recv.size() != 8*fei::numProcs(comm)) ERReturn(-1);

  for(unsigned i=0; i<recv.size(); i++) {
    if (recv[i] != 1) ERReturn(-1);
  }

  //use a zero-length send-buffer on odd-numbered processors
  if (fei::localProc(comm)%2 != 0) send.resize(0);

  CHK_ERR( fei::Allgatherv(comm, send, recvLengths, recv) );

  int expectedLength = 0;
  for(int p=0; p<fei::numProcs(comm); p++) {
    if (p%2 == 0) expectedLength += 8;
  }

  if ((int)recvLengths.size() != fei::numProcs(comm)) ERReturn(-1);
  if ((int)recv.size() != expectedLength) ERReturn(-1);

  for(unsigned j=0; j<recv.size(); j++) {
    if (recv[j] != 1) ERReturn(-1);
  }

  std::vector<int> local(5, fei::localProc(comm)), global;

  CHK_ERR( fei::GlobalMax(comm, local, global) );

  if (global.size() != local.size()) ERReturn(-1);

  for(unsigned i=0; i<global.size(); i++) {
    if (global[i] != fei::numProcs(comm)-1) ERReturn(-1);
  }

  FEI_COUT << "ok" << FEI_ENDL;

  return(0);
}

bool test_CommUtils::run(MPI_Comm comm)
{
  test_fei_Allgatherv(comm);

  if (test_CommUtils_test1(comm) != 0) {
    throw std::runtime_error("test_CommUtils_test1 failed.");
  }

  return true;
}

