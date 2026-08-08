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

#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_mesh/base/Ngp.hpp"
#include "stk_mesh/base/NgpMesh.hpp"
#include "stk_mesh/base/NgpField.hpp"
#include "stk_mesh/base/NgpParallelComm.hpp"
#include "stk_mesh/base/NgpParallelDataExchange.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"
#include <vector>

namespace stk::mesh {

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

} // namespace stk::mesh

#endif

