// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_GetParameter.cpp
 *  \brief Convenience methods for working with the parameter list.
 */

#include <Zoltan2_GetParameter.hpp>

namespace Zoltan2{

/*! \brief Return a sublist of the given parameter list.
 *   \param pl  a Teuchos::ParameterList.
 *   \param listName  the name of a parameter list that may be a sublist
 *                          of pl.
 *   \return the requested parameter list if it exists, an empty list
 *                    otherwise.
 *
 *  If the sublist does not exist, an empty parameter list named "emptyList"
 *  is returned.  If the input list \c pl is such an empty list, then
 *  the empty list is returned.  In this way getList() can be nested when
 *  it is not known if the intermediate lists exist.  For example:
 *
   \code
         getList(getList(getParameters(), "partitioning"), "geometric")
   \endcode
 *
 * will work (by returning an empty list) even if there is
 * no "partitioning" list.
 */

const Teuchos::ParameterList & getParameterList(
  const Teuchos::ParameterList &superList, const char *listName)
{
  static Teuchos::ParameterList emptyList("emptyList");

  if (superList.name() == std::string("emptyList"))
    return superList;

  const Teuchos::ParameterEntry *sublist = superList.getEntryPtr(listName);

  if (!sublist || !sublist->isList()){
    return emptyList;
  }

  return superList.sublist(listName);
}

} // namespace Zoltan2

