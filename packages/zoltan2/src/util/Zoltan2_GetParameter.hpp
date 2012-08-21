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

/*! \file Zoltan2_GetParameter.hpp
 *  \brief Convenience methods for working with the parameter list.
 */

#ifndef ZOLTAN2_PARAMETERS_HPP
#define ZOLTAN2_PARAMETERS_HPP

#include <Zoltan2_Standards.hpp>
#include <Teuchos_ParameterList.hpp>

namespace Zoltan2{

const Teuchos::ParameterList & getParameterList(
  const Teuchos::ParameterList &superList, const char *listName);


/*! \brief Find the value of the named parameter in the list.
 *  \param pl A parameter list that may contain the parameter.
 *  \param name  The name of the parameter entry.
 *  \param set  On return, true if the parameter is set and false if
 *               it is not set (does not appear in the parameter list).
 *  \param value On return, if the entry was found, this will be set
 *                    to the value of the entry.  Otherwise it is
 *                    untouched.
 */

template <typename T>
  void getParameterValue(const Teuchos::ParameterList &pl,
    const char *name, bool &isSet, T &value)
{
  isSet = false;

  if (pl.name() == std::string("emptyList"))
    return;

  const Teuchos::ParameterEntry *pe = pl.getEntryPtr(name);

  if (!pe)
    return;

  isSet = true;
  value = pe->getValue<T>(&value);
}

/*! \brief Find the value of a second level parameter.
 *  \param name1  The sublist containing the parameter.
 *  \param name2  The name of the parameter entry.
 *  \param set  On return, true if the parameter is set and false if
 *               it is not set (does not appear in the parameter list).
 *  \param value On return, if the entry was found, this will be set
 *                    to the value of the entry.  Otherwise it is
 *                    untouched.
 */

template <typename T>
  void getParameterValue(const Teuchos::ParameterList &pl,
    const char *name1, const char *name2, 
    bool &isSet, T &value)
{
  getParameterValue(
    getParameterList(pl, name1), name2, 
    isSet, value);
}

/*! \brief Find the value of a third level parameter.
 *  \param name1  The top level sublist.
 *  \param name2  The sublist in \c name1 that contains the parameter.
 *  \param name3  The name of the parameter entry.
 *  \param set  On return, true if the parameter is set and false if
 *               it is not set (does not appear in the parameter list).
 *  \param value On return, if the entry was found, this will be set
 *                    to the value of the entry.  Otherwise it is
 *                    untouched.
 */

template <typename T>
  void getParameterValue(const Teuchos::ParameterList &pl,
    const char *name1, const char *name2, const char *name3,
    bool &isSet, T &value)
{
  getParameterValue(
    getParameterList(getParameterList(pl, name1), name2), 
    name3, isSet, value);
}

} // namespace Zoltan2

#endif
