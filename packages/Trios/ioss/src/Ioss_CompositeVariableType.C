// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <Ioss_CompositeVariableType.h>
#include <assert.h>
#include <cstdio>
#include <Ioss_VariableType.h>
#include <Ioss_Utils.h>
#include <string>

namespace Ioss {
std::string CompositeVariableType::composite_name(const std::string &base, int copies)
{
  static std::string SEPARATOR("*");
  std::string name = base;
  name += SEPARATOR;
  name += Utils::to_string(copies);
  return name;
}

VariableType*
CompositeVariableType::composite_variable_type(const VariableType *inst, int copies)
{
  VariableType* comp_inst = NULL;

  // See if we already constructed this composite type...
  std::string composite_type = CompositeVariableType::composite_name(inst->name(),copies);

  VariableTypeMap::iterator iter = registry().find(composite_type);
  if (iter == registry().end()) {
    // Not found, construct new type...
    comp_inst = new CompositeVariableType(inst, copies, true);
  } else {
    comp_inst = (*iter).second;
  }
  return comp_inst;
}

CompositeVariableType::CompositeVariableType(const VariableType *base_type, int copies, bool delete_me)
  : VariableType(composite_name(base_type->name(),copies),
		 base_type->component_count()*copies, delete_me),
    baseType(base_type), copies_(copies)
{}

CompositeVariableType::CompositeVariableType(const std::string& my_name, int number_components, bool delete_me)
  : VariableType(my_name, number_components, delete_me), baseType(NULL), copies_(0)
{}

std::string CompositeVariableType::label(int which, const char suffix_sep) const
{
  static char tmp_sep[2];

  // NOTE: 'which' is 1-based
  assert(which > 0 && which <= component_count());

  int base_comp = baseType->component_count();
  int which_instance = (which-1) / base_comp;
  int which_base = (which-1) % base_comp;
  
  std::string my_label = baseType->label(which_base+1, suffix_sep);
  if (suffix_sep != 0 && base_comp > 1) {
    tmp_sep[0] = suffix_sep;
    my_label += tmp_sep;
  }
  my_label += VariableType::numeric_label(which_instance+1, copies_, name());
  return my_label;
}
}
