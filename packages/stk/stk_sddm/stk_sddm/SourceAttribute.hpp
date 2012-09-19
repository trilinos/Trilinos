#ifndef SDDM_SOURCE_ATTRIBUTE_HPP
#define SDDM_SOURCE_ATTRIBUTE_HPP

#include <ostream>
#include <string>

#include <stk_sddm/AttributeInterface.hpp>

namespace stk {
namespace sddm {

class SourceAttribute : public AttributeInterface<SourceAttribute>
{
public:
  SourceAttribute()
  {}

  const char *name() const;
  
  virtual std::ostream &dump(std::ostream &os) const = 0;
  
  virtual std::ostream &xml(std::ostream &os, size_t indent) const = 0;

  virtual std::ostream &describe(std::ostream &oss) const = 0;
};
 

class FileSourceAttribute : public SourceAttribute
{
public:
  FileSourceAttribute(const std::string &path, size_t line)
    : SourceAttribute(),
      m_path(path),
      m_line(line)
  {}

  virtual std::ostream &dump(std::ostream &os) const;
  
  virtual std::ostream &xml(std::ostream &os, size_t indent) const;

  virtual std::ostream &describe(std::ostream &oss) const;
  
protected:
  const std::string     m_path;
  const size_t          m_line;
};


class MethodSourceAttribute : public SourceAttribute
{
public:
  MethodSourceAttribute(const char *method)
    : SourceAttribute(),
      m_method(method)
  {}

  virtual std::ostream &dump(std::ostream &os) const;
  
  virtual std::ostream &xml(std::ostream &os, size_t indent) const;

  virtual std::ostream &describe(std::ostream &oss) const;
  
protected:
  const std::string     m_method;
};
 
} // namespace sddm
} // namespace stk

#endif // SDDM_SOURCE_ATTRIBUTE_HPP
