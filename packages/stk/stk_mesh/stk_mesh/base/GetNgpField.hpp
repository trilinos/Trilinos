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
#include "stk_mesh/base/FieldBase.hpp"

namespace stk {
namespace mesh {

template <typename T>
NgpField<T> & get_updated_ngp_field(const FieldBase & stkField)
{
  NgpFieldBase * ngpField = stkField.get_ngp_field();

  if (ngpField == nullptr) {
    NgpField<T>* ngpFieldT = new NgpField<T>(stkField.get_mesh(), stkField);
    ngpField = ngpFieldT;
    stkField.set_ngp_field(ngpField);
    const unsigned numStates = stkField.number_of_states();
    NgpField<T>* stateNgpFields[MaximumFieldStates];
    
    if (numStates > 1) {
      FieldState requestedState = stkField.state();

      for (unsigned state = 0; state < numStates; ++state) {
        FieldState fieldState = static_cast<FieldState>(state);
        if (fieldState != requestedState) {
          FieldBase* stateStkField = stkField.field_state(fieldState);
          NgpField<T>* stateNgpField = new NgpField<T>(stateStkField->get_mesh(), *stateStkField);
          stateStkField->set_ngp_field(stateNgpField);
          stateNgpFields[state] = stateNgpField;
        }
        else {
          stateNgpFields[state] = ngpFieldT;
        }
      }   
    }   
    else {
      stateNgpFields[0] = ngpFieldT;
    }   

    for (unsigned state = 0; state < numStates; ++state) {
      stateNgpFields[state]->set_field_states(stateNgpFields);
    }
  }
  else {
    if (stkField.get_mesh().synchronized_count() != ngpField->synchronized_count()) {
      ngpField->update_field();
    }
  }

  return dynamic_cast< NgpField<T>& >(*ngpField);
}

}}

#endif // GETNGPFIELD_HPP
