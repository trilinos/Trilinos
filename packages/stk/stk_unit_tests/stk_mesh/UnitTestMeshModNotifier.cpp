// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 



#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <ostream>                      // for operator<<, ostream
#include <stk_mesh/base/ModificationNotifier.hpp>
#include <stk_mesh/base/ModificationObserver.hpp>
#include <stk_mesh/base/Types.hpp>      // for EntityRank
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_util/util/ReportHandler.hpp>  // for ThrowRequireMsg
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for all_reduce_max
#include <vector>                       // for vector
#include "mpi.h"                        // for ompi_communicator_t, etc

namespace {

void reduce_max(std::vector<size_t> &localVector, stk::ParallelMachine comm)
{
  std::vector<size_t> globalVector(localVector.size());
  stk::all_reduce_max(comm, localVector.data(), globalVector.data(), globalVector.size());
  localVector.swap(globalVector);
}

class TestListener : public stk::mesh::ModificationObserver {
public:
  TestListener(stk::ParallelMachine comm)
    : stk::mesh::ModificationObserver(stk::mesh::ModificationObserverPriority::APPLICATION),
      m_comm(comm),
      m_local_entities_created_or_deleted(5, 0),
      m_entity_comm_info_changed(5, 0),
      m_buckets_changed(5, 0)
  {}

  virtual void local_entities_created_or_deleted_notification(stk::mesh::EntityRank rank)
  {
    STK_ThrowRequireMsg(m_local_entities_created_or_deleted.size() > static_cast<unsigned>(rank), "TestListener::local_entities_created_or_deleted ERROR, rank ("<<rank<<") out of range.");
    m_local_entities_created_or_deleted[rank]++;
  }

  virtual void local_entity_comm_info_changed_notification(stk::mesh::EntityRank rank)
  {
    STK_ThrowRequireMsg(m_entity_comm_info_changed.size() > static_cast<unsigned>(rank), "TestListener::global_entity_comm_info_changed ERROR, rank ("<<rank<<") out of range.");
    m_entity_comm_info_changed[rank]++;
  }

  virtual void local_buckets_changed_notification(stk::mesh::EntityRank rank)
  {
    STK_ThrowRequireMsg(m_buckets_changed.size() > static_cast<unsigned>(rank), "TestListener::global_buckets_changed ERROR, rank ("<<rank<<") out of range.");
    m_buckets_changed[rank]++;
  }

  virtual void finished_modification_end_notification()
  {
    reduce_max(m_local_entities_created_or_deleted, m_comm);
    reduce_max(m_entity_comm_info_changed, m_comm);
    reduce_max(m_buckets_changed, m_comm);
  }


  size_t get_local_entities_created_or_deleted(stk::mesh::EntityRank rank) const
  {
    return m_local_entities_created_or_deleted[rank];
  }

  size_t get_global_entity_comm_info_changed(stk::mesh::EntityRank rank) const
  {
    return m_entity_comm_info_changed[rank];
  }

  size_t get_global_buckets_changed(stk::mesh::EntityRank rank) const
  {
    return m_buckets_changed[rank];
  }

private:
  stk::ParallelMachine m_comm;
  std::vector<size_t> m_local_entities_created_or_deleted;
  std::vector<size_t> m_entity_comm_info_changed;
  std::vector<size_t> m_buckets_changed;
};

}

TEST(MeshModNotifier, testLocalEvents)
{
  std::shared_ptr<TestListener> listener = std::make_shared<TestListener>(MPI_COMM_WORLD);
  stk::mesh::ModificationNotifier notifier;
  notifier.register_observer(listener);

  EXPECT_EQ(0u, listener->get_local_entities_created_or_deleted(stk::topology::NODE_RANK));

  notifier.notify_local_entities_created_or_deleted(stk::topology::NODE_RANK);

  EXPECT_EQ(1u, listener->get_local_entities_created_or_deleted(stk::topology::NODE_RANK));
}

TEST(MeshModNotifier, testGlobalEvents)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(comm) == 2)
  {
    int procId = stk::parallel_machine_rank(comm);

    std::shared_ptr<TestListener> listener = std::make_shared<TestListener>(comm);
    stk::mesh::ModificationNotifier notifier;
    notifier.register_observer(listener);

    EXPECT_EQ(0u, listener->get_global_entity_comm_info_changed(stk::topology::NODE_RANK));
    EXPECT_EQ(0u, listener->get_global_buckets_changed(stk::topology::NODE_RANK));

    if (procId == 0) {
      notifier.notify_local_entity_comm_info_changed(stk::topology::NODE_RANK);
      notifier.notify_local_buckets_changed(stk::topology::NODE_RANK);
    }
    notifier.notify_finished_modification_end(comm);

    EXPECT_EQ(1u, listener->get_global_entity_comm_info_changed(stk::topology::NODE_RANK));
    EXPECT_EQ(1u, listener->get_global_buckets_changed(stk::topology::NODE_RANK));
  }
}

