// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "BelosSolverFactory.hpp"
#include <locale>

namespace Belos {
namespace Impl {

void
printStringArray (std::ostream& out,
                  const Teuchos::ArrayView<const std::string>& array)
{
  typedef Teuchos::ArrayView<std::string>::const_iterator iter_type;

  out << "[";
  for (iter_type iter = array.begin(); iter != array.end(); ++iter) {
    out << "\"" << *iter << "\"";
    if (iter + 1 != array.end()) {
      out << ", ";
    }
  }
  out << "]";
}

void
printStringArray (std::ostream& out,
                  const std::vector<std::string>& array)
{
  Teuchos::ArrayView<const std::string> av;
  if (array.size () == 0) {
    av = Teuchos::ArrayView<const std::string> (NULL, 0);
  }
  else {
    av = Teuchos::ArrayView<const std::string> (&array[0], array.size ());
  }
  printStringArray (out, av);
}

std::string
upperCase (const std::string& s)
{
  typedef std::string::value_type char_t;
  typedef std::ctype<char_t> facet_type;
  const facet_type& facet = std::use_facet<facet_type> (std::locale ());

  const std::string::size_type len = s.size ();
  std::string s_uc (s);
  for (std::string::size_type k = 0; k < len; ++k) {
    s_uc[k] = facet.toupper (s[k]);
  }

  return s_uc;
}

} // namespace Impl
} // namespace Belos


