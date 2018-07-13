//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

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


