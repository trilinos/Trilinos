
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <fei_CommMap.hpp>
#include <fei_CommUtils.hpp>
#include <iostream>

namespace {

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(exchangeCommMapData, test0, T)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int numProcs = fei::numProcs(comm);
  if (numProcs == 1) return;

  int localProc = fei::localProc(comm);

  typename fei::CommMap<T>::Type send_comm_map;
  typename fei::CommMap<T>::Type recv_comm_map;

  std::vector<T> send_items(localProc+1);
  for(int i=0; i<localProc+1; ++i) send_items[i] = static_cast<T>(i);

  for(int i=0; i<numProcs; ++i) {
    if (i == localProc) continue;

    fei::addItemsToCommMap<T>(i, send_items.size(), &send_items[0], send_comm_map);
  }

  fei::exchangeCommMapData<T>(comm, send_comm_map, recv_comm_map);

  int rsize = recv_comm_map.size();
  TEUCHOS_TEST_EQUALITY(rsize, numProcs-1, out, success);

  typename fei::CommMap<T>::Type::iterator
    r_iter = recv_comm_map.begin(), r_end = recv_comm_map.end();

  for(; r_iter != r_end; ++r_iter) {
    int rproc = r_iter->first;
    TEUCHOS_TEST_EQUALITY(rproc!=localProc, true, out, success);

    std::vector<T>& rvec = r_iter->second;
    int rvsize = rvec.size();
    TEUCHOS_TEST_EQUALITY(rvsize, rproc+1, out, success);
  }
}

#define UNIT_TEST_GROUP(TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(exchangeCommMapData,test0,TYPE)

typedef long int longint;
UNIT_TEST_GROUP(int)
UNIT_TEST_GROUP(longint)
UNIT_TEST_GROUP(float)
UNIT_TEST_GROUP(double)

}//namespace <anonymous>

