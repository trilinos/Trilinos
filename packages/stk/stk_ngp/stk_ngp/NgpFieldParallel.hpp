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

#ifndef NGPFIELDPARALLEL_HPP
#define NGPFIELDPARALLEL_HPP

#include <vector>
#include "stk_ngp/Ngp.hpp"
#include "stk_ngp/NgpField.hpp"
#include "stk_mesh/base/FieldParallel.hpp"

namespace ngp {

template <typename T>
void parallel_sum(stk::mesh::BulkData & bulk, const std::vector<ngp::Field<T> *> & ngpFields)
{
  stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  const std::vector<stk::mesh::FieldBase *> & allStkFields = meta.get_fields();

  std::vector<const stk::mesh::FieldBase *> stkFields;
  for (ngp::Field<T> * ngpField : ngpFields) {
    stkFields.push_back(allStkFields[ngpField->get_ordinal()]);
    ngpField->sync_to_host();
  }

  stk::mesh::parallel_sum(bulk, stkFields);

  for (ngp::Field<T> * ngpField : ngpFields) {
    ngpField->modify_on_host();
    ngpField->sync_to_device();
  }
}

}

#endif

