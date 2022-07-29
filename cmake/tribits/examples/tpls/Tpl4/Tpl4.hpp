#ifndef TPL4_HPP
#define TPL4_HPP

#include "Tpl2a.hpp"
#include "Tpl2b.hpp"
#include "Tpl3.hpp"

#include <string>

namespace Tpl4 {

/** \brief . */
inline std::string itsme()
{
  return "Tpl4";
}

/** \brief . */
inline std::string deps()
{
  std::ostringstream oss;
  oss
    << Tpl3::itsme() << "{" << Tpl3::deps() << "}"
    << ", "
    << Tpl2::a_itsme() << "{" << Tpl2::a_deps() << "}"
    << ", "
    << Tpl2::b_itsme() << "{" << Tpl2::b_deps() << "}";
  return oss.str();
}

} // namespace Tpl4

#endif // TPL4_HPP
