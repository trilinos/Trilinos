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

#ifndef STK_MESH_NGPSTATEDFIELD_HPP
#define STK_MESH_NGPSTATEDFIELD_HPP

#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_util/stk_config.h>

namespace stk {
namespace mesh {

const unsigned MAX_NUM_FIELD_STATES = 6;

template<typename T>
class NgpMultistateField
{
public:
  NgpMultistateField() {}

  NgpMultistateField(stk::mesh::BulkData& bulk, const stk::mesh::FieldBase &stkField) :
    numOldStates(stkField.number_of_states() - 1)
  {
    stk::mesh::FieldBase * fieldOfState = stkField.field_state(stk::mesh::StateNew);
    fieldNewState = stk::mesh::NgpField<T>(bulk, *fieldOfState);
    for(unsigned i = 0; i < numOldStates; i++)
    {
      fieldOfState = stkField.field_state(get_old_state_for_index(i));
      fieldOldStates[i] = stk::mesh::NgpConstField<T>(bulk, *fieldOfState);
    }
  }

  STK_FUNCTION
  virtual ~NgpMultistateField() {}

  STK_FUNCTION
  unsigned get_num_states() const
  {
    return numOldStates + 1;
  }

  STK_FUNCTION
  stk::mesh::NgpField<T> get_field_new_state() const
  {
    return fieldNewState;
  }

  STK_FUNCTION
  stk::mesh::NgpConstField<T> get_field_old_state(stk::mesh::FieldState state) const
  {
    NGP_ThrowRequire(state != stk::mesh::StateNew);
    return fieldOldStates[state-1];
  }

  STK_FUNCTION
  void increment_state()
  {
    for(unsigned i=numOldStates-1; i>0; i--) {
      fieldOldStates[i].swap(fieldOldStates[i-1]);
    }
    fieldOldStates[0].swap(fieldNewState);
  }

  void sync_to_host()
  {
    copy_device_to_host();
  }

  void sync_to_device()
  {
    copy_host_to_device();
  }

  void modify_on_host()
  {
    fieldNewState.modify_on_host();
  }

  void modify_on_device()
  {
    fieldNewState.modify_on_device();
  }

private:
  void clear_sync_state()
  {
    fieldNewState.clear_sync_state();
    for(unsigned i=0; i<numOldStates; i++)
    {
      fieldOldStates[i].clear_sync_state();
    }
  }

  void copy_device_to_host()
  {
    fieldNewState.sync_to_host();
    for(unsigned i=0; i<numOldStates; i++)
    {
      fieldOldStates[i].sync_to_host();
    }
  }

  void copy_host_to_device()
  {
    fieldNewState.sync_to_device();
    for(unsigned i=0; i<numOldStates; i++)
    {
      fieldOldStates[i].sync_to_device();
    }
  }

  stk::mesh::FieldState get_old_state_for_index(unsigned i)
  {
    return static_cast<stk::mesh::FieldState>(i+1);
  }

  unsigned numOldStates;
  stk::mesh::NgpField<T> fieldNewState;
  stk::mesh::NgpConstField<T> fieldOldStates[MAX_NUM_FIELD_STATES];
};

template<typename T>
class ConvenientNgpMultistateField : public NgpMultistateField<T>
{
public:
  ConvenientNgpMultistateField(stk::mesh::BulkData& bulk, const stk::mesh::FieldBase &stkField) :
    NgpMultistateField<T>(bulk, stkField)
  {}

  STK_FUNCTION
  T& get_new(const stk::mesh::NgpMesh& ngpMesh, stk::mesh::Entity entity, int component) const
  {
    return NgpMultistateField<T>::get_field_new_state().get(ngpMesh, entity, component);
  }

  STK_FUNCTION
  T get_old(stk::mesh::FieldState state, const stk::mesh::NgpMesh& ngpMesh, stk::mesh::Entity entity, int component) const
  {
    return NgpMultistateField<T>::get_field_old_state(state).get(ngpMesh, entity, component);
  }

  STK_FUNCTION
  T& get_new(stk::mesh::NgpMesh::MeshIndex entity, int component) const
  {
    return NgpMultistateField<T>::get_field_new_state().get(entity, component);
  }

  STK_FUNCTION
  T get_old(stk::mesh::FieldState state, stk::mesh::NgpMesh::MeshIndex entity, int component) const
  {
    return NgpMultistateField<T>::get_field_old_state(state).get(entity, component);
  }
};

}
}

#endif
