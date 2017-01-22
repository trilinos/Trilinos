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

#include "Belos_Details_EBelosSolverType.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Belos {
namespace Details {

EBelosSolverType
getEnumFromCanonicalName (const std::string& canonicalName)
{
  // Mapping from canonical solver name (a string) to its
  // corresponding enum value.  This mapping is one-to-one.
  //
  // NOTE (mfh 12 Aug 2015) The keys need to be all uppercase.
  if (canonicalName == "BLOCK GMRES") {
    return SOLVER_TYPE_BLOCK_GMRES;
  }
  else if (canonicalName == "PSEUDOBLOCK GMRES") {
    return SOLVER_TYPE_PSEUDO_BLOCK_GMRES;
  }
  else if (canonicalName == "BLOCK CG") {
    return SOLVER_TYPE_BLOCK_CG;
  }
  else if (canonicalName == "PSEUDOBLOCK CG") {
    return SOLVER_TYPE_PSEUDO_BLOCK_CG;
  }
  else if (canonicalName == "PSEUDOBLOCK STOCHASTIC CG") {
    return SOLVER_TYPE_STOCHASTIC_CG;
  }
  else if (canonicalName == "GCRODR") {
    return SOLVER_TYPE_GCRODR;
  }
  else if (canonicalName == "RCG") {
    return SOLVER_TYPE_RCG;
  }
  else if (canonicalName == "MINRES") {
    return SOLVER_TYPE_MINRES;
  }
  else if (canonicalName == "LSQR") {
    return SOLVER_TYPE_LSQR;
  }
  else if (canonicalName == "TFQMR") {
    return SOLVER_TYPE_TFQMR;
  }
  else if (canonicalName == "PSEUDOBLOCK TFQMR") {
    return SOLVER_TYPE_PSEUDO_BLOCK_TFQMR;
  }
  else if (canonicalName == "HYBRID BLOCK GMRES") {
    return SOLVER_TYPE_GMRES_POLY;
  }
  else if (canonicalName == "PCPG") {
    return SOLVER_TYPE_PCPG;
  }
  else if (canonicalName == "FIXED POINT") {
    return SOLVER_TYPE_FIXED_POINT;
  }
  else if (canonicalName == "BICGSTAB") {
    return SOLVER_TYPE_BICGSTAB;
  }
  else {
    // An error code; not actually a valid value.
    return SOLVER_TYPE_UPPER_BOUND;
  }
}

std::pair<std::string, bool>
getCanonicalNameFromAlias (const std::string& candidateAlias)
{
  // NOTE (mfh 29 Nov 2011) Accessing the flexible capability requires
  // setting a parameter in the solver's parameter list.  This affects
  // the SolverFactory's interface, since using the "Flexible GMRES"
  // alias requires modifying the user's parameter list if necessary.
  // This is a good idea because users may not know about the
  // parameter, or may have forgotten.
  //
  // NOTE (mfh 12 Aug 2015) The keys and values need to be all uppercase.

  if (candidateAlias == "GMRES") {
    return std::make_pair (std::string ("PSEUDOBLOCK GMRES"), true);
  }
  else if (candidateAlias == "BLOCK GMRES") {
    return std::make_pair (std::string ("BLOCK GMRES"), true);
  }
  else if (candidateAlias == "FLEXIBLE GMRES") {
    return std::make_pair (std::string ("BLOCK GMRES"), true);
  }
  else if (candidateAlias == "CG") {
    return std::make_pair (std::string ("PSEUDOBLOCK CG"), true);
  }
  else if (candidateAlias == "PSEUDOBLOCKCG") {
    return std::make_pair (std::string ("PSEUDOBLOCK CG"), true);
  }
  else if (candidateAlias == "STOCHASTIC CG") {
    return std::make_pair (std::string ("PSEUDOBLOCK STOCHASTIC CG"), true);
  }
  else if (candidateAlias == "RECYCLING CG") {
    return std::make_pair (std::string ("RCG"), true);
  }
  else if (candidateAlias == "RECYCLING GMRES") {
    return std::make_pair (std::string ("GCRODR"), true);
  }
  // For compatibility with Stratimikos' Belos adapter.
  else if (candidateAlias == "PSEUDO BLOCK GMRES") {
    return std::make_pair (std::string ("PSEUDOBLOCK GMRES"), true);
  }
  else if (candidateAlias == "PSEUDOBLOCKGMRES") {
    return std::make_pair (std::string ("PSEUDOBLOCK GMRES"), true);
  }
  else if (candidateAlias == "PSEUDO BLOCK CG") {
    return std::make_pair (std::string ("PSEUDOBLOCK CG"), true);
  }
  else if (candidateAlias == "PSEUDOBLOCKCG") {
    return std::make_pair (std::string ("PSEUDOBLOCK CG"), true);
  }
  else if (candidateAlias == "TRANSPOSE-FREE QMR") {
    return std::make_pair (std::string ("TFQMR"), true);
  }
  else if (candidateAlias == "PSEUDO BLOCK TFQMR") {
    return std::make_pair (std::string ("PSEUDOBLOCK TFQMR"), true);
  }
  else if (candidateAlias == "PSEUDO BLOCK TRANSPOSE-FREE QMR") {
    return std::make_pair (std::string ("PSEUDOBLOCK TFQMR"), true);
  }
  else if (candidateAlias == "GMRESPOLY") {
    return std::make_pair (std::string ("HYBRID BLOCK GMRES"), true);
  }
  else if (candidateAlias == "SEED GMRES") {
    return std::make_pair (std::string ("HYBRID BLOCK GMRES"), true);
  }
  else if (candidateAlias == "CGPOLY") {
    return std::make_pair (std::string ("PCPG"), true);
  }
  else if (candidateAlias == "SEED CG") {
    return std::make_pair (std::string ("PCPG"), true);
  }
  else if (candidateAlias == "FIXED POINT") {
    return std::make_pair (std::string ("FIXED POINT"), true);
  }
  else if (candidateAlias == "BICGSTAB") {
    return std::make_pair (std::string ("BICGSTAB"), true);
  }
  else { // not a known alias
    return std::make_pair (candidateAlias, false);
  }
}

std::vector<std::string>
solverNameAliases ()
{
#ifdef HAVE_TEUCHOSCORE_CXX11
  return {
    {"GMRES"},
    {"BLOCK GMRES"},
    {"FLEXIBLE GMRES"},
    {"CG"},
    {"PSEUDOBLOCKCG"},
    {"STOCHASTIC CG"},
    {"RECYCLING CG"},
    {"RECYCLING GMRES"},
    {"PSEUDO BLOCK GMRES"},
    {"PSEUDOBLOCKGMRES"},
    {"PSEUDO BLOCK CG"},
    {"PSEUDOBLOCKCG"},
    {"TRANSPOSE-FREE QMR"},
    {"PSEUDO BLOCK TFQMR"},
    {"PSEUDO BLOCK TRANSPOSE-FREE QMR"},
    {"GMRESPOLY"},
    {"SEED GMRES"},
    {"CGPOLY"},
    {"SEED CG"},
    {"FIXED POINT"},
    {"BICGSTAB"}
  };
#else // NOT HAVE_TEUCHOSCORE_CXX11
  std::vector<std::string> names;

  names.push_back ("GMRES");
  names.push_back ("BLOCK GMRES");
  names.push_back ("FLEXIBLE GMRES");
  names.push_back ("CG");
  names.push_back ("PSEUDOBLOCKCG");
  names.push_back ("STOCHASTIC CG");
  names.push_back ("RECYCLING CG");
  names.push_back ("RECYCLING GMRES");
  names.push_back ("PSEUDO BLOCK GMRES");
  names.push_back ("PSEUDOBLOCKGMRES");
  names.push_back ("PSEUDO BLOCK CG");
  names.push_back ("PSEUDOBLOCKCG");
  names.push_back ("TRANSPOSE-FREE QMR");
  names.push_back ("PSEUDO BLOCK TFQMR");
  names.push_back ("PSEUDO BLOCK TRANSPOSE-FREE QMR");
  names.push_back ("GMRESPOLY");
  names.push_back ("SEED GMRES");
  names.push_back ("CGPOLY");
  names.push_back ("SEED CG");
  names.push_back ("FIXED POINT");
  names.push_back ("BICGSTAB");

  return names;
#endif // HAVE_TEUCHOSCORE_CXX11
}

int numSupportedSolvers () {
  return static_cast<int> (SOLVER_TYPE_UPPER_BOUND);
}

std::vector<std::string>
canonicalSolverNames ()
{
#ifdef HAVE_TEUCHOSCORE_CXX11
  return {
    {"BLOCK GMRES"},
    {"PSEUDOBLOCK GMRES"},
    {"BLOCK CG"},
    {"PSEUDOBLOCK CG"},
    {"PSEUDOBLOCK STOCHASTIC CG"},
    {"GCRODR"},
    {"RCG"},
    {"MINRES"},
    {"LSQR"},
    {"TFQMR"},
    {"PSEUDOBLOCK TFQMR"},
    {"HYBRID BLOCK GMRES"},
    {"PCPG"},
    {"FIXED POINT"},
    {"BICGSTAB"}
  };
#else
  std::vector<std::string> canonicalNames (numSupportedSolvers ());

  canonicalNames.push_back ("BLOCK GMRES");
  canonicalNames.push_back ("PSEUDOBLOCK GMRES");
  canonicalNames.push_back ("BLOCK CG");
  canonicalNames.push_back ("PSEUDOBLOCK CG");
  canonicalNames.push_back ("PSEUDOBLOCK STOCHASTIC CG");
  canonicalNames.push_back ("GCRODR");
  canonicalNames.push_back ("RCG");
  canonicalNames.push_back ("MINRES");
  canonicalNames.push_back ("LSQR");
  canonicalNames.push_back ("TFQMR");
  canonicalNames.push_back ("PSEUDOBLOCK TFQMR");
  canonicalNames.push_back ("HYBRID BLOCK GMRES");
  canonicalNames.push_back ("PCPG");
  canonicalNames.push_back ("FIXED POINT");
  canonicalNames.push_back ("BICGSTAB");

  return canonicalNames;
#endif // HAVE_TEUCHOSCORE_CXX11
}

void
reviseParameterListForAlias (const std::string& aliasName,
                             Teuchos::ParameterList& solverParams)
{
  if (aliasName == "FLEXIBLE GMRES") {
    // "Gmres" uses title case in this solver's parameter list.  For
    // our alias, we prefer the all-capitals "GMRES" that the
    // algorithm's authors (Saad and Schultz) used.
    solverParams.set ("Flexible Gmres", true);
  }
}

} // namespace Details
} // namespace Belos

