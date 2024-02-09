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

#ifndef STK_MESH_BASE_NGPFIELDBLAS_HPP
#define STK_MESH_BASE_NGPFIELDBLAS_HPP

#include <stk_util/stk_config.h>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>

#include <complex>
#include <string>
#include <iostream>
#include <algorithm>

namespace stk {
namespace mesh {

namespace impl {
//************ implementation detail, not for public use ********************
//************ public functions start below after impl namespace ************

template<class Scalar, class EXEC_SPACE>
inline
void field_fill_no_sync_or_mark(const Scalar alpha, DeviceField<Scalar>& ngpField, const EXEC_SPACE& execSpace)
{
  auto ngpView = impl::get_device_data(ngpField);
  Kokkos::deep_copy(execSpace, ngpView, alpha);
}

template<class Scalar, class EXEC_SPACE>
inline
void field_fill_no_sync_or_mark(const Scalar alpha, HostField<Scalar>& ngpField, const EXEC_SPACE& execSpace)
{
  stk::mesh::field_fill(alpha, *ngpField.get_field_base());
}

//************ end of implementation detail *********************************
} // namespace impl

template<class Scalar, typename EXEC_SPACE>
inline
void field_fill(const Scalar alpha, const FieldBase& field, const EXEC_SPACE& execSpace,
                bool IsDeviceExecSpaceUserOverride = (!std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>))
{
  bool isActuallyDeviceExecSpace = !std::is_same_v<stk::ngp::HostExecSpace,EXEC_SPACE>;
  bool operateOnDevice = isActuallyDeviceExecSpace;
#ifndef STK_USE_DEVICE_MESH
  operateOnDevice = false;
#endif

  field.clear_sync_state();

  if (operateOnDevice) {
    NgpField<Scalar>& ngpField = get_updated_ngp_field<Scalar>(field);
    impl::field_fill_no_sync_or_mark(alpha, ngpField, execSpace);
  }
  else {
    stk::mesh::field_fill(alpha, field);
  }

  if (IsDeviceExecSpaceUserOverride) {
    field.modify_on_device();
  }
  else {
    field.modify_on_host();
  }
}

} // mesh
} // stk

#endif // STK_MESH_BASE_NGPFIELDBLAS_HPP

