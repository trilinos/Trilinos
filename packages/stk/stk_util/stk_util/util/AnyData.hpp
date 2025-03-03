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

#ifndef SIERRA_UTILITY_ANY_DATA_H
#define SIERRA_UTILITY_ANY_DATA_H

#include "stk_util/util/Writer.hpp"     // for operator"", Writer
#include "stk_util/diag/StringUtil.hpp" // for demangle
#include "stk_util/util/Writer_fwd.hpp" // for Writer
#include <typeinfo>
#include <stdexcept>
#include <string>
#include <ostream>
#include <vector>

namespace sierra {

std::ostream& operator<<(std::ostream& os,
                         const std::vector<double>& t); /// WORKAROUND: Intel 10.0 cannot find symbol for some reason
std::ostream&
operator<<(std::ostream& os, const std::vector<int>& t); /// WORKAROUND: Intel 10.0 cannot find symbol for some reason

/**
 * @brief Class <b>bad_any_data_cast</b> is a bad cast exception thrown when attempting
 * top cast a resource to a different type from the resource.
 */
class bad_any_data_cast : public std::bad_cast {
 public:
  bad_any_data_cast(const std::string& message) noexcept : m_message(message) {}

  virtual const char* what() const noexcept override { return m_message.c_str(); }

 private:
  std::string m_message;
};

/**
 * @brief Function <b>cast_error</b> creates a <b>bad_any_data_cast</b>, message
 * describing the invalid conversion.
 */
inline bad_any_data_cast cast_error(const std::type_info& from_type, const std::type_info& to_type)
{
  std::ostringstream strout;
  strout << "Cannot cast resource data of type " << demangle(from_type.name()) << " to type "
         << demangle(to_type.name());
  return bad_any_data_cast(strout.str());
}

/**
 * @brief Interface class <b>Data&lt;void&gt;</b> is used as the base class of all resources.  Only
 * object derived from <b>Data%lt;void&gt;</b> may be inserted into the <b>Resources.</b>
 */
class AnyData {
 public:
  AnyData() {}

  virtual ~AnyData() {}

  template <typename T>
  const T& value() const;

  /**
   * @brief Member function <b>value</b> attempts to cast the data to the
   * specified type.  If the data is not of the specified type, an exception is thrown.
   *
   * @return a <b>T</b> reference to the data.
   */
  template <typename T>
  T& value();

  /**
   * @brief Pure virtual member function <b>type</b> returns the type of the value
   * stored in the <b>data</b> object.
   *
   * @return a <b>std::type_info</b> const reference to the type
   */
  virtual const std::type_info& type() const = 0;

  virtual std::ostream& dump(std::ostream& os) const = 0;
};

/**
 * @brief template<b>Data&lt;T&gt;</b> defines the interface for a resource.
 * Every resource derived class must implement a <b>getValue()</b> and a const
 * <b>getValue()</b> virtual member functions.
 *
 * These functions are designed for the usage of properties to be implementation
 * independent from the storage of the actual data for the resource.
 */
template <typename T>
class Data : public AnyData {
 public:
  Data() {}

  virtual ~Data() {}

  AnyData& operator=(const T& t)
  {
    value() = t;
    return *this;
  }

  /**
   * @brief Member function <b>getValue</b> returns a const reference to the
   * resource storage.
   *
   * @return a <b>T</b> const reference to the resource storage.
   */
  virtual const T& value() const = 0;

  /**
   * @brief Member function <b>getValue</b> returns a reference to the resource
   * storage.
   *
   * @return a <b>T</b> reference to the resource storge
   */
  virtual T& value() = 0;

  /**
   * @brief Pure virtual member function <b>type</b> returns the type of the value
   * stored.
   *
   * @return a <b>std::type_info</b> const reference to the type.
   */
  virtual const std::type_info& type() const override { return typeid(T); }

  virtual std::ostream& dump(std::ostream& os) const override
  {
    T v = value();
    return os << v;
  }
};

template <typename T>
inline const T& AnyData::value() const
{
  if(typeid(T) != type()) throw cast_error(type(), typeid(T));
  const Data<T>* t = static_cast<const Data<T>*>(this);
  return t->value();
}

/**
 * @brief Member function <b>value</b> attempts to cast the data to the
 * specified type.  If the data is not of the specified type, an exception is thrown.
 *
 * @return a <b>T</b> reference to the data.
 *
 */
template <typename T>
inline T& AnyData::value()
{
  if(typeid(T) != type()) throw cast_error(type(), typeid(T));
  Data<T>* t = static_cast<Data<T>*>(this);
  return t->value();
}

/**
 * @brief template<b>Value</b> defines a resource which stores it's data by * value.
 */
template <typename T>
class Value : public Data<T> {
 public:
  Value(const T& t)
    : Data<T>()
    , m_t(t)
  {
  }

  virtual ~Value() {}

  virtual const T& value() const override { return m_t; }

  virtual T& value() override { return m_t; }

  Value& operator=(const T& t)
  {
    m_t = t;
    return *this;
  }

 private:
  T m_t; ///< Data value
};

/**
 * @brief template<b>Reference</b> defines a resource which stores it's data by reference.
 */
template <typename T>
class Reference : public Data<T> {
 public:
  Reference(T& t)
    : Data<T>()
    , m_t(t)
  {
  }

  virtual ~Reference() {}

  virtual const T& value() const override { return m_t; }

  virtual T& value() override { return m_t; }

  Reference& operator=(const T& t)
  {
    m_t = t;
    return *this;
  }

 private:
  T& m_t; ///< Data reference
};

Diag::Writer& operator<<(Diag::Writer& dout, const AnyData& data);

inline std::ostream& operator<<(std::ostream& os, const std::vector<double>& t)
{
  for(unsigned i = 0; i < t.size(); i++)
    os << t[i] << " ";
  return os;
}
inline std::ostream& operator<<(std::ostream& os, const std::vector<int>& t)
{
  for (unsigned i = 0; i < t.size(); i++) os << t[i] << " ";
  return os;
}

inline Diag::Writer& operator<<(Diag::Writer& dout, const AnyData& data)
{
  std::ostringstream ss;
  data.dump(ss);
  return dout << ss.str();
}

} // namespace sierra

#endif // SIERRA_UTILITY_ANY_DATA_H
