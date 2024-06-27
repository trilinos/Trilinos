// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
