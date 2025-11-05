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

#ifndef stk_mesh_FieldParallel_hpp
#define stk_mesh_FieldParallel_hpp

#include <stk_util/stk_config.h>
#include <stk_mesh/base/Types.hpp>      // for EntityProc
#include <stk_mesh/base/FieldBase.hpp>  // for FieldBase
#include <stk_mesh/base/BulkData.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_mesh/base/NgpParallelComm.hpp>
#include <stk_util/util/ReportHandler.hpp>  // for ThrowRequireMsg

#include <stddef.h>                     // for size_t
#include <vector>                       // for vector

namespace stk { namespace mesh { class Ghosting; } }

namespace stk::mesh {

/**
 * This file contains some helper functions that are part of the Field API.
 */

/** Send field-data from entities to their ghosts, for a specified 'ghosting'.
 * For entities that are ghosted, this function updates field-data from the
 * original entity to the ghosts.
 */
void communicate_field_data(const Ghosting& ghosts, const std::vector<const FieldBase*>& fields);

void communicate_field_data(const BulkData& mesh, const std::vector<const FieldBase*>& fields);

/** Copy data for the given fields, from owned entities to shared-but-not-owned entities.
 * I.e., shared-but-not-owned entities get an update of the field-data from the owned entity.
*/
inline
void copy_owned_to_shared( const BulkData& mesh,
                           const std::vector< const FieldBase *> & fields)
{
  communicate_field_data(*mesh.ghostings()[BulkData::SHARED], fields);
}

//----------------------------------------------------------------------

/** Sum/Max/Min (assemble) field-data for the specified fields on shared entities such that each shared entity
 * will have the same field values on each sharing proc.
 */
void parallel_sum(const BulkData& mesh, const std::vector<const FieldBase*>& fields, bool deterministic = true);
void parallel_max(const BulkData& mesh, const std::vector<const FieldBase*>& fields);
void parallel_min(const BulkData& mesh, const std::vector<const FieldBase*>& fields);

void parallel_sum_including_ghosts(const BulkData & mesh, const std::vector<const FieldBase *> & fields, bool deterministic = true);
void parallel_max_including_ghosts(const BulkData & mesh, const std::vector<const FieldBase *> & fields, bool deterministic = true);
void parallel_min_including_ghosts(const BulkData & mesh, const std::vector<const FieldBase *> & fields, bool deterministic = true);

template <typename NgpSpace, Operation OP, bool includeGhosts>
void ngp_parallel_op(const BulkData & bulk,
                     const std::vector<const FieldBase*> & fields,
                     bool deterministic)
{
  if (bulk.parallel_size() == 1 || fields.empty()) return;

  if (fields[0]->type_is<double>()) {
    ngp_parallel_data_excahnge_sym_pack_unpack<double, OP, NgpSpace>(bulk.parallel(), bulk, fields, includeGhosts, deterministic);
  }
  else if (fields[0]->type_is<int>()) {
    ngp_parallel_data_excahnge_sym_pack_unpack<int, OP, NgpSpace>(bulk.parallel(), bulk, fields, includeGhosts, deterministic);
  }
  else if (fields[0]->type_is<float>()) {
    ngp_parallel_data_excahnge_sym_pack_unpack<float, OP, NgpSpace>(bulk.parallel(), bulk, fields, includeGhosts, deterministic);
  }
  else if (fields[0]->type_is<long double>()) {
    ngp_parallel_data_excahnge_sym_pack_unpack<long double, OP, NgpSpace>(bulk.parallel(), bulk, fields, includeGhosts, deterministic);
  }
  else if (fields[0]->type_is<unsigned long>()) {
    ngp_parallel_data_excahnge_sym_pack_unpack<unsigned long, OP, NgpSpace>(bulk.parallel(), bulk, fields, includeGhosts, deterministic);
  }
}

template <typename NgpSpace>
void parallel_sum(const BulkData & bulk,
                  const std::vector<const FieldBase*> & fields,
                  bool deterministic = true)
{
  if constexpr (std::is_same_v<NgpSpace, stk::ngp::HostSpace>) {
    parallel_sum(bulk, fields, deterministic);
  }
  else {
    constexpr bool includeGhosts = false;
    ngp_parallel_op<NgpSpace, Operation::SUM, includeGhosts>(bulk, fields, deterministic);
  }
}

template <typename NgpSpace>
void parallel_sum_including_ghosts(const BulkData & bulk,
                                   const std::vector<const FieldBase*> & fields,
                                   bool deterministic = true)
{
  if constexpr (std::is_same_v<NgpSpace, stk::ngp::HostSpace>) {
    parallel_sum_including_ghosts(bulk, fields, deterministic);
  }
  else {
    constexpr bool includeGhosts = true;
    ngp_parallel_op<NgpSpace, Operation::SUM, includeGhosts>(bulk, fields, deterministic);
  }
}

template <typename NgpSpace>
void parallel_min(const BulkData & bulk, const std::vector<const FieldBase*> & fields, bool deterministic = true)
{
  if constexpr (std::is_same_v<NgpSpace, stk::ngp::HostSpace>) {
    parallel_min(bulk, fields);
  }
  else {
    constexpr bool includeGhosts = false;
    ngp_parallel_op<NgpSpace, Operation::MIN, includeGhosts>(bulk, fields, deterministic);
  }
}

template <typename NgpSpace>
void parallel_min_including_ghosts(const BulkData & bulk,
                                   const std::vector<const FieldBase*> & fields,
                                   bool deterministic = true)
{
  if constexpr (std::is_same_v<NgpSpace, stk::ngp::HostSpace>) {
    parallel_min_including_ghosts(bulk, fields, deterministic);
  }
  else {
    constexpr bool includeGhosts = true;
    ngp_parallel_op<NgpSpace, Operation::MIN, includeGhosts>(bulk, fields, deterministic);
  }
}

template <typename NgpSpace>
void parallel_max(const BulkData & bulk,
                  const std::vector<const FieldBase*> & fields,
                  bool deterministic = true)
{
  if constexpr (std::is_same_v<NgpSpace, stk::ngp::HostSpace>) {
    parallel_max(bulk, fields);
  }
  else {
    constexpr bool includeGhosts = false;
    ngp_parallel_op<NgpSpace, Operation::MAX, includeGhosts>(bulk, fields, deterministic);
  }
}

template <typename NgpSpace>
void parallel_max_including_ghosts(const BulkData & bulk,
                                   const std::vector<const FieldBase*> & fields,
                                   bool deterministic = true)
{
  if constexpr (std::is_same_v<NgpSpace, stk::ngp::HostSpace>) {
    parallel_max_including_ghosts(bulk, fields, deterministic);
  }
  else {
    constexpr bool includeGhosts = true;
    ngp_parallel_op<NgpSpace, Operation::MAX, includeGhosts>(bulk, fields, deterministic);
  }
}
} // namespace stk::mesh

#endif

