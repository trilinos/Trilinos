
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <snl_fei_RecordCollection.hpp>
#include <fei_SharedIDs.hpp>
#include <fei_CommUtils.hpp>
#include <fei_set_shared_ids.hpp>
#include <fei_FieldMask.hpp>
#include <iostream>

namespace {

//The fei::set_shared_ids function (which this test exercises) is not
//currently a template. But it will be in the future, which is why this
//test uses the template mechanism.

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(set_shared_ids, test0, T)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int numProcs = fei::numProcs(comm);
  if (numProcs == 1) return;

  int localProc = fei::localProc(comm);

  snl_fei::RecordCollection records(localProc);
  int local_id = localProc;

  std::vector<fei::FieldMask*> fieldMasks;

  records.initRecords(1, &local_id, fieldMasks);

  if (localProc > 0) {
    local_id = localProc-1;
    records.initRecords(1, &local_id, fieldMasks);
  }

  if (localProc < numProcs-1) {
    local_id = localProc+1;
    records.initRecords(1, &local_id, fieldMasks);
  }

  //this template parameter will someday be T instead of int:
  fei::SharedIDs<int> shared;

  fei::set_shared_ids(comm, records, shared, records.getMinID(), records.getMaxID());

  fei::SharedIDs<int>::map_type& shID_map = shared.getSharedIDs();

  fei::SharedIDs<int>::map_type::iterator
    s_iter = shID_map.begin(), s_end = shID_map.end();

  if (localProc == 0) {
    size_t expected_num_shIDs = 2;
    TEUCHOS_TEST_EQUALITY(shID_map.size(), expected_num_shIDs, out, success);

    int ID0 = s_iter->first;
    TEUCHOS_TEST_EQUALITY(ID0, 0, out, success);
    std::set<int>& procs0 = s_iter->second;
    size_t expected_numprocs0 = 2;
    TEUCHOS_TEST_EQUALITY(procs0.size(), expected_numprocs0, out, success);

    ++s_iter;

    int ID1 = s_iter->first;
    TEUCHOS_TEST_EQUALITY(ID1, 1, out, success);
    std::set<int>& procs1 = s_iter->second;
    size_t expected_numprocs1 = numProcs>2 ? 3 : 2;
    TEUCHOS_TEST_EQUALITY(procs1.size(), expected_numprocs1, out, success);
  }
  else if (localProc == numProcs-1) {
    size_t expected_num_shIDs = 2;
    TEUCHOS_TEST_EQUALITY(shID_map.size(), expected_num_shIDs, out, success);

    int ID0 = s_iter->first;
    TEUCHOS_TEST_EQUALITY(ID0, localProc-1, out, success);
    std::set<int>& procs0 = s_iter->second;
    size_t expected_numprocs0 = numProcs>2 ? 3 : 2;
    TEUCHOS_TEST_EQUALITY(procs0.size(), expected_numprocs0, out, success);

    ++s_iter;

    int ID1 = s_iter->first;
    TEUCHOS_TEST_EQUALITY(ID1, localProc, out, success);
    std::set<int>& procs1 = s_iter->second;
    size_t expected_numprocs1 = 2;
    TEUCHOS_TEST_EQUALITY(procs1.size(), expected_numprocs1, out, success);
  }
  else {
    size_t expected_num_shIDs = 3;
    TEUCHOS_TEST_EQUALITY(shID_map.size(), expected_num_shIDs, out, success);

    int ID0 = s_iter->first;
    TEUCHOS_TEST_EQUALITY(ID0, localProc-1, out, success);
    std::set<int>& procs0 = s_iter->second;
    TEUCHOS_TEST_EQUALITY(procs0.size(), 3, out, success);

    ++s_iter;

    int ID1 = s_iter->first;
    TEUCHOS_TEST_EQUALITY(ID1, localProc, out, success);
    std::set<int>& procs1 = s_iter->second;
    TEUCHOS_TEST_EQUALITY(procs1.size(), 3, out, success);

    ++s_iter;

    int ID2 = s_iter->first;
    TEUCHOS_TEST_EQUALITY(ID2, localProc+1, out, success);
    std::set<int>& procs2 = s_iter->second;
    TEUCHOS_TEST_EQUALITY(procs2.size(), 3, out, success);
  }
}

#define UNIT_TEST_GROUP(TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(set_shared_ids,test0,TYPE)

UNIT_TEST_GROUP(int)
// typedef long int longint;
// UNIT_TEST_GROUP(longint)

}//namespace <anonymous>

