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

#ifndef BELOS_DETAILS_EBELOSSOLVERTYPE_HPP
#define BELOS_DETAILS_EBELOSSOLVERTYPE_HPP

/// \file Belos_Details_EBelosSolverType.hpp
/// \brief Declaration of alias functions for solver names.

#include <string>
#include <utility>
#include <vector>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Teuchos {
  // forward declaration (avoid unnecessary include)
  class ParameterList;
} // namespace Teuchos
#endif // NOT DOXYGEN_SHOULD_SKIP_THIS

namespace Belos {
namespace Details {

/// \brief Number of Belos solvers supported for any linear algebra
///   implementation ("generically").
///
/// This may differ from the number of supported solver
/// <i>names</i>, since we may accept multiple names ("aliases") for
/// some solvers.
int numSupportedSolvers ();

/// \brief Get the candidate canonical name for a given candidate
///   alias.
///
/// \param candidateAlias [in] Must be upper case.
///
/// \return The candidate canonical name, and whether the candidate
///   canonical name was recognized as an alias.  If it was NOT
///   recognized as an alias, this function just returns its input.
///
/// The keys of this map do not necessarily include canonical solver
/// names.  If a candidate name isn't a key in this map, then it must
/// be a canonical name in order to be valid.  There doesn't need to
/// be an alias for each solver.
///
/// \note To Belos developers: If you want to add a new alias, first
///   add the mapping from alias to canonical solver name in the
///   SolverFactory constructor.  Then, edit
///   Details::reviseParameterListForAlias to do any modifications of
///   the input ParameterList associated with that alias.
std::pair<std::string, bool>
getCanonicalNameFromAlias (const std::string& candidateAlias);

/// \brief List of supported aliases (to canonical solver names).
///
/// This list does _not_ include canonical names.
std::vector<std::string>
solverNameAliases ();

/// \brief List of canonical solver names.
///
/// This list does _not_ include aliases.
std::vector<std::string>
canonicalSolverNames ();

/// \brief Modify the input ParameterList appropriately for the given
///   solver alias.
///
/// Some aliases include modifications or special checking of the
/// input ParameterList.  All alias-related ParameterList revision
/// happens in this method.
void
reviseParameterListForAlias (const std::string& aliasName,
                             Teuchos::ParameterList& solverParams);

} // namespace Details
} // namespace Belos

#endif // BELOS_DETAILS_EBELOSSOLVERTYPE_HPP
