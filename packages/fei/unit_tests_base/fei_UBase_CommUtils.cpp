
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <fei_iostream.hpp>
#include "fei_CommUtils.hpp"

#include <vector>
#include <cmath>

#undef fei_file
#define fei_file "fei_unit_CommUtils.cpp"
#include "fei_ErrMacros.hpp"

namespace {

TEUCHOS_UNIT_TEST(CommUtils, fei_Allgatherv)
{
  MPI_Comm comm = MPI_COMM_WORLD;

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
}

TEUCHOS_UNIT_TEST(CommUtils, test_CommUtils_test1)
{
  MPI_Comm comm = MPI_COMM_WORLD;

  std::vector<int> send(8, 1), recv;
  std::vector<int> recvLengths;

  fei::Allgatherv(comm, send, recvLengths, recv);

  TEUCHOS_TEST_EQUALITY((int)recvLengths.size(), fei::numProcs(comm), out, success);
  TEUCHOS_TEST_EQUALITY((int)recv.size(), 8*fei::numProcs(comm), out, success);

  for(unsigned i=0; i<recv.size(); i++) {
    TEUCHOS_TEST_EQUALITY(recv[i], 1, out, success);
  }

  //use a zero-length send-buffer on odd-numbered processors
  if (fei::localProc(comm)%2 != 0) send.resize(0);

  fei::Allgatherv(comm, send, recvLengths, recv);

  int expectedLength = 0;
  for(int p=0; p<fei::numProcs(comm); p++) {
    if (p%2 == 0) expectedLength += 8;
  }

  TEUCHOS_TEST_EQUALITY((int)recvLengths.size(), fei::numProcs(comm), out, success);
  TEUCHOS_TEST_EQUALITY((int)recv.size(), expectedLength, out, success);

  for(unsigned j=0; j<recv.size(); j++) {
    TEUCHOS_TEST_EQUALITY(recv[j], 1, out, success);
  }

  std::vector<int> local(5, fei::localProc(comm)), global;

  fei::GlobalMax(comm, local, global);

  TEUCHOS_TEST_EQUALITY(global.size(), local.size(), out, success);

  for(unsigned i=0; i<global.size(); i++) {
    TEUCHOS_TEST_EQUALITY(global[i], fei::numProcs(comm)-1, out, success);
  }
}

}//namespace <anonymous>
