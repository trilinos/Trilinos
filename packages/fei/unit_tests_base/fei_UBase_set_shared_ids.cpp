/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


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

  fei::set_shared_ids(comm, records, shared);

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

