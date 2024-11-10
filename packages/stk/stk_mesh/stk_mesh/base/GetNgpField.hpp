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

#ifndef GETNGPFIELD_HPP
#define GETNGPFIELD_HPP

#include "stk_mesh/base/NgpField.hpp"
#include "stk_mesh/base/NgpFieldSyncDebugger.hpp"
#include "stk_mesh/base/FieldBase.hpp"

namespace stk {
namespace mesh {

template <typename T, template <typename> class NgpDebugger = DefaultNgpFieldSyncDebugger>
NgpField<T, NgpDebugger> & get_updated_ngp_field_async(const FieldBase & stkField, const stk::ngp::ExecSpace& execSpace)
{
  NgpFieldBase * ngpField = impl::get_ngp_field(stkField);

  if (ngpField == nullptr) {
    ngpField = new NgpField<T, NgpDebugger>(stkField.get_mesh(), stkField, true);
    ngpField->update_field(execSpace);
    ngpField->debug_initialize_debug_views();
    impl::set_ngp_field(stkField, ngpField);
    ngpField->clear_host_sync_state();
  }
  else {
    if (stkField.get_mesh().synchronized_count() != ngpField->synchronized_count()) {
      ngpField->update_field(execSpace);
    }
  }

  return dynamic_cast< NgpField<T, NgpDebugger>& >(*ngpField);
}

template <typename T, template <typename> class NgpDebugger = DefaultNgpFieldSyncDebugger>
NgpField<T, NgpDebugger> & get_updated_ngp_field_async(const FieldBase & stkField, stk::ngp::ExecSpace&& execSpace)
{
  NgpFieldBase * ngpField = impl::get_ngp_field(stkField);

  if (ngpField == nullptr) {
    ngpField = new NgpField<T, NgpDebugger>(stkField.get_mesh(), stkField, true);
    ngpField->update_field(std::forward<stk::ngp::ExecSpace>(execSpace));
    ngpField->debug_initialize_debug_views();
    impl::set_ngp_field(stkField, ngpField);
    ngpField->clear_host_sync_state();
  }
  else {
    if (stkField.get_mesh().synchronized_count() != ngpField->synchronized_count()) {
      ngpField->update_field(std::forward<stk::ngp::ExecSpace>(execSpace));
    }
  }

  return dynamic_cast< NgpField<T, NgpDebugger>& >(*ngpField);
}

template <typename T, template <typename> class NgpDebugger = DefaultNgpFieldSyncDebugger>
NgpField<T, NgpDebugger> & get_updated_ngp_field(const FieldBase & stkField)
{
  auto& ngpFieldRef = get_updated_ngp_field_async<T, NgpDebugger>(stkField, Kokkos::DefaultExecutionSpace());
  ngpFieldRef.fence();
  return ngpFieldRef;
}

}}

#endif // GETNGPFIELD_HPP

