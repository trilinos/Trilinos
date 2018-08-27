//   ------------------------------------------------------------
//   Copyright 2006 - 2010 Sandia Corporation.
//   Under the terms of Contract DE-AC04-94AL85000, there is a
//   non-exclusive license for use of this work by or on behalf
//   of the U.S. Government.  Export of this program may require
//   a license from the United States Government.
//   ------------------------------------------------------------

#include "stk_util/diag/StringUtil.hpp"   // for demangle, case_strcmp
#include "stk_util/diag/WriterExt.hpp"    // for operator<<
#include "stk_util/util/AnyData.hpp"      // for Value, AnyData
#include "stk_util/util/Writer.hpp"       // for Writer, operator<<, dendl
#include <stk_util/util/Resource2.h>
#include <algorithm> // for find
#include <sstream>   // for operator<<, basic_ostream, etc
#include <stdexcept> // for runtime_error
#include <stk_util/util/ReportHandler.hpp>
#include <string> // for char_traits, operator<<, etc

namespace sierra {
namespace Rsrc2 {

bad_any_data_cast cast_error(const Resource& resource, const AnyData* data, const std::type_info& to_type)
{
  std::ostringstream strout;
  strout << "Cannot cast resource '" << resource.path();
  if (data) strout << "' of type " << demangle(data->type().name());
  strout << "' to type " << demangle(to_type.name());
  if (data) strout << " because it has no value assigned";
  return bad_any_data_cast(strout.str());
}

std::runtime_error not_found(const String& resource_name, const String& parent_path)
{
  std::ostringstream strout;

  strout << "Resource '" << resource_name << "' not found in '" << parent_path << "'";

  return std::runtime_error(strout.str());
}

struct Resource_ {
  Resource_(Resource_* parent, const String& resource_name, AnyData* value)
    : m_parent(parent)
    , m_name(resource_name)
    , m_value(value)
    , m_resourceList()
  {
  }

  ~Resource_()
  {
    delete m_value;

    for (auto it = m_resourceList.begin(); it != m_resourceList.end(); ++it)
      if ((*it).m_resource_->m_parent && (*it).m_resource_->m_parent == this) delete (*it).m_resource_;
  }

  Resource_* m_parent;
  String m_name;
  AnyData* m_value;
  ResourceList m_resourceList;

  String path() const
  {
    std::string s(m_name.c_str());
    auto parent = m_parent;
    while (parent) {
      s = std::string(parent->m_name.c_str()).append(".").append(s);
      parent = parent->m_parent;
    }
    return s.c_str();
  }
};

ResourceList& getResourceMap(Resource_* resource) { return resource->m_resourceList; }

ResourceRoot::ResourceRoot(const String& resource_name)
  : Resource(0)
{
  m_resource_ = new Resource_(0, resource_name, 0);
}

ResourceRoot::~ResourceRoot() { destroy(); }

const String& Resource::name() const { return m_resource_->m_name; }

String Resource::path() const { return m_resource_->path(); }

AnyData* Resource::data() { return m_resource_->m_value; }

const AnyData* Resource::data() const { return m_resource_->m_value; }

Diag::Writer& Resource::verbose_print(Diag::Writer& dout) const
{
  if (dout.shouldPrint()) {
    if (m_resource_) {
      dout << *m_resource_;
    }
    else {
      dout << "<empty>";
    }
  }

  return dout;
}

Resource::iterator Resource::begin()
{
  auto& resource_list = getResourceMap(m_resource_);
  return resource_list.begin();
}

Resource::const_iterator Resource::begin() const
{
  const auto& resource_list = getResourceMap(m_resource_);
  return resource_list.begin();
}

Resource::iterator Resource::end()
{
  auto& resource_list = getResourceMap(m_resource_);
  return resource_list.end();
}

Resource::const_iterator Resource::end() const
{
  const auto& resource_list = getResourceMap(m_resource_);
  return resource_list.end();
}

Resource Resource::create(const String& resource_name, AnyData* any_data)
{
  auto& resource_list = getResourceMap(m_resource_);

  ThrowRequireMsg(resource_list.find(resource_name) == resource_list.end(),
                  "Resource " << resource_name << " already exists in resource group ");

  auto new_resource_ = new Resource_(m_resource_, resource_name, any_data);

  resource_list.insert(Resource(new_resource_));

  return Resource(new_resource_);
}

void Resource::destroy()
{
  delete m_resource_;
  m_resource_ = nullptr;
}

Resource_* find(Resource_* resource_, const String& resource_name)
{
  ResourceList* resource_list = &getResourceMap(resource_);

  auto c = resource_name.begin();

  while (1) {
    auto d = std::find(c, resource_name.end(), SEP);
    String t(c, d);

    auto it = resource_list->find(t);
    if (it == resource_list->end())
      return 0;
    else if (d == resource_name.end())
      return (*it).m_resource_;
    else {
      c = ++d;
      resource_list = &getResourceMap((*it).m_resource_);
    }
  }
}

Resource Resource::get(const String& resource_name) const
{
  auto resource_ = find(m_resource_, resource_name);
  if (!resource_) throw not_found(resource_name, path());
  return Resource(resource_);
}

bool Resource::exists(const String& resource_name) const { return find(m_resource_, resource_name); }

namespace {

bool match2(const String& resource_name, const String& s)
{
  if (resource_name.size() < s.size()) return false;

  const char* n_start = resource_name.c_str();
  const char* n = n_start + resource_name.size() - s.size();
  while (n != n_start && *n != SEP)
    --n;
  if (*n == SEP) ++n;

  return case_strcmp(n, s.c_str()) == 0;
}

} // namespace <unnamed>

void Resource::match(const String& resource_name, ResourceList& result_list) const
{
  if (match2(path(), resource_name)) result_list.insert(*this);

  ResourceList& resource_list = getResourceMap(m_resource_);
  resource_list.match(resource_name, result_list);
}

void ResourceList::match(const String& resource_name, ResourceList& result_list) const
{
  for (auto it = begin(); it != end(); ++it) {
    (*it).match(resource_name, result_list);
  }
}

Diag::Writer& operator<<(Diag::Writer& dout, const ResourceList& resource_list)
{
  for (auto it = resource_list.begin(); it != resource_list.end(); ++it) {
    dout << (*it);
  }
  return dout;
}

Diag::Writer& operator<<(Diag::Writer& dout, const Resource_& resource_)
{
  if (dout.shouldPrint()) {
    dout << "[" << resource_.path() << "]";

    if (resource_.m_value) dout << " " << *resource_.m_value;

    dout << Diag::dendl << resource_.m_resourceList;
  }

  return dout;
}

} // namespace Rsrc2
} // namespace sierra
