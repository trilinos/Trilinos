#ifndef HEADER_ONLY_TPL_STUFF_HPP
#define HEADER_ONLY_TPL_STUFF_HPP

#include <string>

namespace HeaderOnlyTpl {

/** \brief . */
template <typename T>
T  sqr(const T &v)
{
  return v*v;
}

/** \brief . */
inline std::string itsme()
{
  return "headeronlytpl";
}

} // namespace HeaderOnlyTpl


#endif // HEADER_ONLY_TPL_STUFF_HPP
