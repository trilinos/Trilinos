#include <sstream>
#include <iomanip>
#include <stk_sddm/SourceAttribute.hpp>

namespace stk {
namespace sddm {

const char *
SourceAttribute::name() const {
  return "Source";
}


std::ostream &
FileSourceAttribute::dump(std::ostream &os) const {
  os << "@" << name() << "(Path=\"" << m_path << ", Line=" << m_line << ")";
  return os;
}
  

std::ostream &
FileSourceAttribute::xml(std::ostream &os, size_t indent) const {
  os << std::setw(indent*2) << "" << "<Attribute type=\"" << name() << "\">" << std::endl;
  os << std::setw(indent*2) << "" << "  <Path>" << m_path << "</Path>" << std::endl;
  os << std::setw(indent*2) << "" << "  <Line>" << m_line << "</Line>" << std::endl;
  os << std::setw(indent*2) << "" << "</Attribute>" << std::endl;
  
  return os;
}


std::ostream &
FileSourceAttribute::describe(
  std::ostream &        oss) const
{
  oss << m_path << ":" << m_line;

  return oss;
}

std::ostream &
MethodSourceAttribute::dump(std::ostream &os) const {
  os << "@" << name() << "(Method=\"" << m_method << ")";
  return os;
}
  

std::ostream &
MethodSourceAttribute::xml(std::ostream &os, size_t indent) const {
  os << std::setw(indent*2) << "" << "<Attribute type=\"" << name() << "\">" << std::endl;
  os << std::setw(indent*2) << "" << "  <Method>" << m_method << "</Method>" << std::endl;
  os << std::setw(indent*2) << "" << "</Attribute>" << std::endl;
  
  return os;
}


std::ostream &
MethodSourceAttribute::describe(
  std::ostream &        oss) const
{
  oss << m_method;

  return oss;
}

} // namespace sddm
} // namespace stk

