// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_PARALLELCOMMHELPERS_H_
#define KRINO_INCLUDE_AKRI_PARALLELCOMMHELPERS_H_

#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/EntityLess.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <algorithm>
#include <vector>

namespace krino {

// PACK_OP shoudl have an operator()(stk::CommBuffer & buffer, const stk::mesh::EntityKey & entity_key)
// UNPACK_OP should have an operator()(CommBuffer & buffer)
template <class PACK_OP, class UNPACK_OP>
void communicate_with_shared_procs(const stk::mesh::BulkData& mesh,
    const std::vector<stk::mesh::Entity> & local_entities,
    PACK_OP pack_op,
    UNPACK_OP unpack_op)
{
  // The input vector, local_entities, must only contain entities that are locally_owned or globally_shared (or both).

  stk::CommSparse comm_spec(mesh.parallel());
  std::vector<int> sharing_procs;

  for (int phase=0;phase<2;++phase)
  {
    for (auto&& local_entity : local_entities)
    {
      ThrowAssert(mesh.bucket(local_entity).owned() || mesh.bucket(local_entity).shared());
      mesh.comm_shared_procs(local_entity, sharing_procs);
      for (auto&& sharing_proc : sharing_procs)
      {
        if (sharing_proc != mesh.parallel_rank())
        {
          pack_op(comm_spec.send_buffer(sharing_proc), mesh.entity_key(local_entity));
        }
      }
    }

    if ( phase == 0 )
    {
      comm_spec.allocate_buffers();
    }
    else
    {
      comm_spec.communicate();
    }
  }

  for(int i = 0; i < mesh.parallel_size(); ++i)
  {
    if(i != mesh.parallel_rank())
    {
      auto & recv_buffer = comm_spec.recv_buffer(i);
      while(recv_buffer.remaining())
      {
        unpack_op(recv_buffer);
      }
    }
  }
}

}

#endif /* KRINO_INCLUDE_AKRI_PARALLELCOMMHELPERS_H_ */
