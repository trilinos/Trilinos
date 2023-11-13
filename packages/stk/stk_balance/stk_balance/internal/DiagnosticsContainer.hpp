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
#ifndef DIAGNOSTICSCONTAINER_HPP
#define DIAGNOSTICSCONTAINER_HPP

#include "stk_util/diag/StringUtil.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_balance/internal/Diagnostics.hpp"
#include <vector>
#include <typeinfo>
#include <utility>
#include <cstddef>
#include <type_traits>
#include <algorithm>
#include <string>

namespace stk::balance {

class DiagnosticsContainer
{
public:
  DiagnosticsContainer() = default;
  ~DiagnosticsContainer();

  using TypeVector = std::vector<const std::type_info *>;
  using ValueVector = std::vector<Diagnostic *>;

  template <typename T, typename... Args>
  void register_diagnostic(Args&&... args);

  template <typename T>
  T * get() const;

  size_t size() const { return m_types.size(); }
  void clear();

  ValueVector::const_iterator begin() const { return m_values.begin(); }
  ValueVector::const_iterator end() const { return m_values.end(); }

private:
  TypeVector m_types;
  ValueVector m_values;
};


template <typename T, typename... Args>
void DiagnosticsContainer::register_diagnostic(Args&&... args)
{
  static_assert(std::is_base_of_v<Diagnostic, T>, "Registered type must have a base class of stk::balance::Diagnostic");
  if constexpr (std::is_base_of_v<Diagnostic, T>)
  {
    const std::type_info * newType = &typeid(T);
    STK_ThrowRequireMsg(std::none_of(m_types.begin(), m_types.end(),
                                 [newType](const std::type_info * existingType) { return *newType == *existingType; }),
                    std::string("Can only register each Diagnostic once: ") + sierra::demangle(newType->name()));
    m_types.push_back(newType);
    m_values.push_back(new T{std::forward<Args>(args)...});
  }
}

template <typename T>
T *
DiagnosticsContainer::get() const
{
  const std::type_info * requestedType = &typeid(T);
  auto it = std::find_if(m_types.begin(), m_types.end(),
                         [requestedType](const std::type_info * existingType) { return *requestedType == *existingType; });
  if (it != m_types.end()) {
    return static_cast<T*>(m_values[it - m_types.begin()]);
  }
  else {
    return nullptr;
  }
}


namespace impl {
extern DiagnosticsContainer g_diagnosticsContainer;
}


template <typename T, typename... Args>
void register_diagnostic(Args&&... args) {
  impl::g_diagnosticsContainer.register_diagnostic<T>(std::forward<Args>(args)...);
}

template <typename T>
T * get_diagnostic() {
  return impl::g_diagnosticsContainer.get<T>();
}

} // namespace stk::balance

#endif // DIAGNOSTICSCONTAINER_HPP
