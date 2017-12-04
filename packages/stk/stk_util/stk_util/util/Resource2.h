//   ------------------------------------------------------------
//   Copyright 2006 - 2010 Sandia Corporation.
//   Under the terms of Contract DE-AC04-94AL85000, there is a
//   non-exclusive license for use of this work by or on behalf
//   of the U.S. Government.  Export of this program may require
//   a license from the United States Government.
//   ------------------------------------------------------------

#ifndef SIERRA_SLIB_RESOURCE2_H
#define SIERRA_SLIB_RESOURCE2_H

#include "stk_util/diag/String.hpp"       // for String, operator<, etc
#include "stk_util/util/AnyData.hpp"      // for Value, AnyData
#include "stk_util/util/Writer_fwd.hpp"   // for Writer
#include <algorithm>                      // for find_if, lower_bound
#include <stddef.h>                       // for size_t
#include <typeinfo>                       // for type_info
#include <vector>                         // for vector, etc

namespace sierra {
namespace Rsrc2 {
class Resource;
class ResourceList;
struct Resource_;
}
}

namespace sierra {
namespace Rsrc2 {

static const char SEP = '.'; ///< Resource path name separator

/**
 * @brief Function <b>cast_error</b> creates a <b>bad_any_data_cast</b>, assembling a
 * message describing the name and invalid conversion.
 */
bad_any_data_cast cast_error(const Resource& resource_name, const AnyData* data, const std::type_info& to_type);

Resource_* find(Resource_* resource_, const String& resource_name);

Diag::Writer& operator<<(Diag::Writer& dout, const Resource_& resource_);

class Resource {
  friend struct Resource_;
  friend Resource_* find(Resource_* resource_, const String& resource_name);

 public:
  typedef std::vector<Resource>::iterator iterator;
  typedef std::vector<Resource>::const_iterator const_iterator;

  Resource()
    : m_resource_(nullptr)
  {
  }

  explicit Resource(Resource_* resource_)
    : m_resource_(resource_)
  {
  }

  ~Resource() {}

  const String& name() const;

  String path() const;

  const AnyData* data() const;

  AnyData* data();

  template <typename T>
  bool isType() const
  {
    return data() && typeid(T) == data()->type();
  }

  const std::type_info& type() const { return data()->type(); }

  template <typename T>
  const T& value() const
  {
    if (data() && data()->type() == typeid(T)) {
      const Data<T>* data_t = static_cast<const Data<T>*>(data());
      return data_t->value();
    }
    throw cast_error(*this, data(), typeid(T));
  }

  template <typename T>
  T& value()
  {
    if (data() && data()->type() == typeid(T)) {
      Data<T>* data_t = static_cast<Data<T>*>(data());
      return data_t->value();
    }
    throw cast_error(*this, data(), typeid(T));
  }

  Resource get(const String& resource_name) const;

  void match(const String& resource_name, ResourceList& result_list) const;

  Resource create(const String& resource_name, AnyData* data = nullptr);

  template <typename T>
  Resource create_value(const String& resource_name, T t)
  {
    return create(resource_name, new Value<T>(t));
  }

  void destroy();

  bool exists(const String& resource_name) const;

  template <typename T>
  T& value(const String& resource_name)
  {
    return get(resource_name).value<T>();
  }

  template <typename T>
  const T& value(const String& resource_name) const
  {
    return get(resource_name).value<T>();
  }

  iterator begin();
  const_iterator begin() const;

  iterator end();
  const_iterator end() const;

  Diag::Writer& verbose_print(Diag::Writer& dout) const;

 protected:
  Resource_* m_resource_;
};

inline bool operator<(const Resource& left, const Resource& right) { return left.name() < right.name(); }

class ResourceRoot : public Resource {
 public:
  explicit ResourceRoot(const String& resource_name);

  ~ResourceRoot();
};

class ResourceList {
 public:
  typedef std::vector<Resource>::iterator iterator;
  typedef std::vector<Resource>::const_iterator const_iterator;

 private:
  struct find_pred {
    explicit find_pred(const String& resource_name)
      : m_name(resource_name)
    {
    }

    bool operator()(const Resource& resource) { return resource.name() == m_name; }

    const String& m_name;
  };

 public:
  size_t size() const { return m_resourceVector.size(); }

  bool empty() const { return m_resourceVector.empty(); }

  Resource front() const { return m_resourceVector.front(); }

  void insert(Resource resource)
  {
    m_resourceVector.insert(std::lower_bound(m_resourceVector.begin(), m_resourceVector.end(), resource), resource);
  }

  iterator find(const String& resource_name) { return std::find_if(begin(), end(), find_pred(resource_name)); }

  void erase(iterator it) { m_resourceVector.erase(it); }

  void match(const String& resource_name, ResourceList& result_list) const;

  iterator begin() { return m_resourceVector.begin(); }

  const_iterator begin() const { return m_resourceVector.begin(); }

  iterator end() { return m_resourceVector.end(); }

  const_iterator end() const { return m_resourceVector.end(); }

 private:
  std::vector<Resource> m_resourceVector;
};

inline Diag::Writer& operator<<(Diag::Writer& dout, const Resource& resource) { return resource.verbose_print(dout); }

Diag::Writer& operator<<(Diag::Writer& dout, const ResourceList& resource_list);

} // namespace Rsrc
} // namespace sierra

#endif // SIERRA_SLIB_RESOURCE2_H
