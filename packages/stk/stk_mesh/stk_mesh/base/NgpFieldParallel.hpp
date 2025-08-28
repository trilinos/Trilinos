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

#ifndef STK_MESH_NGPFIELDPARALLEL_HPP
#define STK_MESH_NGPFIELDPARALLEL_HPP

#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_mesh/base/Ngp.hpp"
#include "stk_mesh/base/NgpField.hpp"
#include "stk_mesh/base/NgpParallelComm.hpp"
#include "stk_mesh/base/NgpParallelDataExchange.hpp"
#include <vector>

namespace stk {
namespace mesh {

template <typename T>
void parallel_sum(const stk::mesh::BulkData & bulk, const std::vector<stk::mesh::NgpField<T> *> & ngpFields,
                  bool doFinalSyncBackToDevice = true)
{
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  const std::vector<stk::mesh::FieldBase *> & allStkFields = meta.get_fields();

  std::vector<const stk::mesh::FieldBase *> stkFields;
  for (stk::mesh::NgpField<T> * ngpField : ngpFields) {
    stkFields.push_back(allStkFields[ngpField->get_ordinal()]);
    ngpField->sync_to_host();
  }

  stk::mesh::parallel_sum(bulk, stkFields);

  for (stk::mesh::NgpField<T> * ngpField : ngpFields) {
    ngpField->modify_on_host();

    if (doFinalSyncBackToDevice) {
      ngpField->sync_to_device();
    }
  }
}

template <typename T>
void parallel_sum_including_ghosts(const stk::mesh::BulkData & bulk,
                                   const std::vector<stk::mesh::NgpField<T> *> & ngpFields,
                                   bool doFinalSyncBackToDevice = true,
                                   bool deterministic = true)
{
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  const std::vector<stk::mesh::FieldBase *> & allStkFields = meta.get_fields();

  std::vector<const stk::mesh::FieldBase *> stkFields;
  for (stk::mesh::NgpField<T> * ngpField : ngpFields) {
    stkFields.push_back(allStkFields[ngpField->get_ordinal()]);
    ngpField->sync_to_host();
  }

  stk::mesh::parallel_sum_including_ghosts(bulk, stkFields, deterministic);

  for (stk::mesh::NgpField<T> * ngpField : ngpFields) {
    ngpField->modify_on_host();

    if (doFinalSyncBackToDevice) {
      ngpField->sync_to_device();
    }
  }
}

template <typename T>
void do_final_sync_to_device(const std::vector<NgpField<T>*>& ngpFields)
{
  for(auto ngpField : ngpFields) {
    ngpField->sync_to_device();
  }
}

template <typename T>
void copy_owned_to_shared(const stk::mesh::BulkData & bulk,
                          const std::vector<stk::mesh::NgpField<T> *> & ngpFields,
                          bool doFinalSyncBackToDevice = true)
{
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  const std::vector<stk::mesh::FieldBase *> & allStkFields = meta.get_fields();

  std::vector<const stk::mesh::FieldBase *> stkFields;
  for (stk::mesh::NgpField<T> * ngpField : ngpFields) {
    stkFields.push_back(allStkFields[ngpField->get_ordinal()]);
  }

  stk::mesh::copy_owned_to_shared(bulk, stkFields);

  if(doFinalSyncBackToDevice) {
    do_final_sync_to_device(ngpFields);
  }
}

template <typename T>
void communicate_field_data(const stk::mesh::Ghosting & ghosting,
                            const std::vector<stk::mesh::NgpField<T> *> & ngpFields,
                            bool doFinalSyncBackToDevice = true)
{
  const stk::mesh::MetaData & meta = ghosting.mesh().mesh_meta_data();
  const std::vector<stk::mesh::FieldBase *> & allStkFields = meta.get_fields();

  std::vector<const stk::mesh::FieldBase *> stkFields;
  for (stk::mesh::NgpField<T> * ngpField : ngpFields) {
    stkFields.push_back(allStkFields[ngpField->get_ordinal()]);
  }

  stk::mesh::communicate_field_data(ghosting, stkFields);

  if(doFinalSyncBackToDevice) {
    do_final_sync_to_device(ngpFields);
  }
}

template <typename T>
void communicate_field_data(const stk::mesh::BulkData & bulk,
                            const std::vector<stk::mesh::NgpField<T> *> & ngpFields,
                            bool doFinalSyncBackToDevice = true)
{
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  const std::vector<stk::mesh::FieldBase *> & allStkFields = meta.get_fields();

  std::vector<const stk::mesh::FieldBase *> stkFields;
  for (stk::mesh::NgpField<T> * ngpField : ngpFields) {
    stkFields.push_back(allStkFields[ngpField->get_ordinal()]);
  }

  stk::mesh::communicate_field_data(bulk, stkFields);

  if(doFinalSyncBackToDevice) {
    do_final_sync_to_device(ngpFields);
  }
}

template <typename T>
void parallel_sum_device_mpi(const stk::mesh::NgpMesh& ngpMesh, const std::vector<stk::mesh::NgpField<T> *> & ngpFields)
{
  const stk::mesh::BulkData& bulk = ngpMesh.get_bulk_on_host();

  for (stk::mesh::NgpField<T> * ngpField : ngpFields) {
    ngpField->sync_to_device();
  }

  ParallelSumDataExchangeSymPackUnpackHandler<stk::mesh::NgpMesh,stk::mesh::NgpField<T>> exchangeHandler(ngpMesh, ngpFields);

  const bool includeGhosts = false;
  const bool deterministic = false;
  stk::mesh::ngp_parallel_data_exchange_sym_pack_unpack<double,Operation::SUM>(bulk.parallel(),
                                                                exchangeHandler,
                                                                includeGhosts,
                                                                deterministic);

  for (stk::mesh::NgpField<T> * ngpField : ngpFields) {
    ngpField->modify_on_device();
  }
}

template <Operation OP, typename NGPMESH, typename NGPFIELD>
void parallel_op_including_ghosts_device_mpi(const NGPMESH& ngpMesh,
                                  const std::vector<NGPFIELD*> & ngpFields,
                                             bool deterministic)
{
  const stk::mesh::BulkData& bulk = ngpMesh.get_bulk_on_host();

  STK_ThrowRequireMsg(bulk.has_symmetric_ghost_info() && bulk.in_synchronized_state(),
     "parallel_sym_including_ghosts_device_mpi requires symmetric ghost info, and the mesh also can not be in a modifiable state.");

  using ThisExecSpace = typename NGPMESH::MeshExecSpace;
  constexpr bool onDevice = !Kokkos::SpaceAccessibility<ThisExecSpace, stk::ngp::HostMemSpace>::accessible;
  for (NGPFIELD * ngpField : ngpFields) {
    if constexpr (onDevice) {
      ngpField->sync_to_device();
    }
    else {
      ngpField->sync_to_host();
    }
  }

  using ExchangeHandler = ParallelSumDataExchangeSymPackUnpackHandler<NGPMESH,NGPFIELD>;
  ExchangeHandler exchangeHandler(ngpMesh, ngpFields);

  const bool includeGhosts = true;
  using T = typename NGPFIELD::value_type;
  stk::mesh::ngp_parallel_data_exchange_sym_pack_unpack<T,OP,ExchangeHandler>(bulk.parallel(),
                                                                           exchangeHandler,
                                                                           includeGhosts,
                                                                           deterministic);
  for (NGPFIELD * ngpField : ngpFields) {
    if constexpr (onDevice) {
      ngpField->modify_on_device();
    }
    else {
      ngpField->modify_on_host();
    }
  }
}

template <typename NGPMESH, typename NGPFIELD>
void parallel_sum_including_ghosts_device_mpi(const NGPMESH& ngpMesh,
                                  const std::vector<NGPFIELD*> & ngpFields,
                                              bool deterministic = true)
{
  parallel_op_including_ghosts_device_mpi<Operation::SUM>(ngpMesh, ngpFields, deterministic);
}

template <typename NgpMesh, typename NgpField, typename MemSpace = stk::ngp::MemSpace>
void parallel_sum(NgpMesh const& ngpMesh, std::vector<NgpField*> const& ngpFields, bool doFinalSyncBackToDevice = true)
{
  if constexpr (!std::is_same_v<Kokkos::DefaultHostExecutionSpace, Kokkos::DefaultExecutionSpace> &&
                Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace, MemSpace>::accessible) {
    parallel_sum_device_mpi(ngpMesh, ngpFields);
  } else {
    parallel_sum(ngpMesh.get_bulk_on_host(), ngpFields, doFinalSyncBackToDevice);
  }
}

template <typename NGPMESH, typename NGPFIELD, typename MemSpace = stk::ngp::MemSpace>
void parallel_sum_including_ghosts(NGPMESH const& ngpMesh, std::vector<NGPFIELD*> const& ngpFields, bool deterministic = true)
{
  STK_ThrowRequireMsg((Kokkos::SpaceAccessibility<typename NGPMESH::MeshExecSpace, MemSpace>::accessible), "parallel_sum_including_ghosts MemSpace not accessible from NGPMESH::MeshExecSpace");

  parallel_op_including_ghosts_device_mpi<Operation::SUM>(ngpMesh, ngpFields, deterministic);
}

template <typename NGPMESH, typename NGPFIELD, typename MemSpace = stk::ngp::MemSpace>
void parallel_max_including_ghosts(NGPMESH const& ngpMesh, std::vector<NGPFIELD*> const& ngpFields, bool deterministic = true)
{
  STK_ThrowRequireMsg((Kokkos::SpaceAccessibility<typename NGPMESH::MeshExecSpace, MemSpace>::accessible), "parallel_max_including_ghosts MemSpace not accessible from NGPMESH::MeshExecSpace");

  parallel_op_including_ghosts_device_mpi<Operation::MAX>(ngpMesh, ngpFields, deterministic);
}

template <typename NGPMESH, typename NGPFIELD, typename MemSpace = stk::ngp::MemSpace>
void parallel_min_including_ghosts(NGPMESH const& ngpMesh, std::vector<NGPFIELD*> const& ngpFields, bool deterministic = true)
{
  STK_ThrowRequireMsg((Kokkos::SpaceAccessibility<typename NGPMESH::MeshExecSpace, MemSpace>::accessible), "parallel_min_including_ghosts MemSpace not accessible from NGPMESH::MeshExecSpace");

  parallel_op_including_ghosts_device_mpi<Operation::MIN>(ngpMesh, ngpFields, deterministic);
}

}
}

#endif

