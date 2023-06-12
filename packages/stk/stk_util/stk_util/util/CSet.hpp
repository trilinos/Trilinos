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

#ifndef stk_util_util_CSet_hpp
#define stk_util_util_CSet_hpp

#include "stk_util/util/ReportHandler.hpp"
#include <typeinfo>  // for type_info
#include <vector>    // for vector
#include <iterator>   // for advance
#include <memory>

namespace stk {

class CSet {
public:
  CSet() = default;
  ~CSet() = default;

  CSet(const CSet & rhs) = default;
  CSet(CSet && rhs) = default;
  CSet & operator=(const CSet & rhs) = default;
  CSet & operator=(CSet && rhs) = default;

  template<class T> const T * get() const;
  template<class T> const T * insert_with_delete(const T * value);
  template<class T> const T * insert_no_delete(const T * value);
  template<class T> bool remove(const T * value);

  using TypeVector = std::vector<const std::type_info *>;
  using ValueVector = std::vector<std::shared_ptr<const void>>;

private:

  template <typename T>
  class Deleter
  {
  public:
    Deleter()
      : m_doDelete(true)
    { }
    ~Deleter() = default;
    void operator()(const void * v) {
      if (m_doDelete) {
        delete reinterpret_cast<const T*>(v);
      }
    }
    void cancel_delete() { m_doDelete = false; }

  private:
    bool m_doDelete;
  };

  class InactiveDeleter
  {
  public:
    InactiveDeleter() = default;
    ~InactiveDeleter() = default;
    void operator()(const void * ) { }
  };

  const void * p_get(const std::type_info & type) const;

  template <typename DELETER, typename T>
  const T * p_insert(const DELETER & deleter, const T * value);

  template <typename T>
  bool p_remove(const std::type_info & manager, const void * value);

  TypeVector m_type;
  ValueVector m_value;
};


template<class T>
inline
const T * CSet::get() const
{
  return static_cast<const T*>(p_get(typeid(T)));
}

template<class T>
inline
const T * CSet::insert_with_delete(const T * arg_value)
{
  return p_insert(Deleter<T>(), arg_value);
}

template<class T>
inline
const T * CSet::insert_no_delete(const T * arg_value)
{
  return p_insert(InactiveDeleter(), arg_value);
}

template<class T>
inline
bool CSet::remove(const T * value)
{
  return p_remove<T>(typeid(T), value);
}


namespace cset {

struct less_cset {
  bool operator()(const std::type_info * lhs,
                  const std::type_info * rhs) const
  {
    return lhs->before(*rhs);
  }
};

CSet::TypeVector::iterator
lower_bound(CSet::TypeVector & v, const std::type_info * t);

}

template <typename DELETER, typename T>
inline
const T * CSet::p_insert(const DELETER & deleter, const T * value)
{
  const std::type_info * type = &typeid(T);

  auto im = cset::lower_bound(m_type, type);
  const size_t offset = im - m_type.begin();

  STK_ThrowAssert(m_value.size() == m_type.size());
  auto iv = m_value.begin();
  std::advance(iv, offset);

  if (im == m_type.end() || *type != **im) {
    im = m_type.insert(im, type);
    iv = m_value.insert(iv , std::shared_ptr<const void>(value, deleter));
  }

  STK_ThrowAssert(iv != m_value.end());
  return static_cast<const T*>((*iv).get());
}

template <typename T>
bool CSet::p_remove(const std::type_info & type, const void * value)
{
  bool result = false;
  const auto im = cset::lower_bound(m_type, &type);

  if (im != m_type.end()) {
    const size_t offset = im - m_type.begin();

    if (offset <= m_value.size()) {
      auto iv = m_value.begin();
      std::advance(iv, offset);

      result = (type == **im) && (value == (*iv).get());

      if (result) {
        m_type.erase(im);

        if (Deleter<T> * deleter = std::get_deleter<Deleter<T>>(*iv)) {
          deleter->cancel_delete();
        }
        m_value.erase(iv);
      }
    }
  }

  return result;
}

}

#endif // stk_util_util_CSet_hpp
