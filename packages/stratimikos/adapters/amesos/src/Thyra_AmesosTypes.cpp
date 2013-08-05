/*
// @HEADER
// ***********************************************************************
// 
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#include "Thyra_AmesosTypes.hpp"

namespace Thyra {

const Amesos::ESolverType Amesos::solverTypeValues[Amesos::numSolverTypes] =
{
  Amesos::LAPACK
#ifdef HAVE_AMESOS_KLU
  ,Amesos::KLU
#endif
#ifdef HAVE_AMESOS_UMFPACK
  ,Amesos::UMFPACK
#endif
#ifdef HAVE_AMESOS_SUPERLU
  ,Amesos::SUPERLU
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
  ,Amesos::SUPERLUDIST
#endif
#ifdef HAVE_AMESOS_TAUCS
  ,Amesos::TAUCS
#endif
#ifdef HAVE_AMESOS_PARDISO
  ,Amesos::PARDISO
#endif
#ifdef HAVE_AMESOS_PASTIX
  ,Amesos::PASTIX
#endif
#ifdef HAVE_AMESOS_PARAKLETE
  ,Amesos::PARAKLETE
#endif
#ifdef HAVE_AMESOS_MUMPS
  ,Amesos::MUMPS
#endif
#ifdef HAVE_AMESOS_SCALAPACK
  ,Amesos::SCALAPACK
#endif
#ifdef HAVE_AMESOS_DSCPACK
  ,Amesos::DSCPACK
#endif
};

const char* Amesos::solverTypeNames[Amesos::numSolverTypes] =
{
  "Lapack"
#ifdef HAVE_AMESOS_KLU
  ,"Klu"
#endif
#ifdef HAVE_AMESOS_UMFPACK
  ,"Umfpack"
#endif
#ifdef HAVE_AMESOS_SUPERLU
  ,"Superlu"
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
  ,"Superludist"
#endif
#ifdef HAVE_AMESOS_TAUCS
  ,"Taucs"
#endif
#ifdef HAVE_AMESOS_PARDISO
  ,"Pardiso"
#endif
#ifdef HAVE_AMESOS_PASTIX
  ,"Pastix"
#endif
#ifdef HAVE_AMESOS_PARAKLETE
  ,"Paraklete"
#endif
#ifdef HAVE_AMESOS_MUMPS
  ,"Mumps"
#endif
#ifdef HAVE_AMESOS_SCALAPACK
  ,"Scalapack"
#endif
#ifdef HAVE_AMESOS_DSCPACK
  ,"Dscpack"
#endif
};

const bool Amesos::supportsUnsymmetric[Amesos::numSolverTypes] =
{
  true
#ifdef HAVE_AMESOS_KLU
  ,true
#endif
#ifdef HAVE_AMESOS_UMFPACK
  ,true
#endif
#ifdef HAVE_AMESOS_SUPERLU
  ,true
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
  ,true
#endif
#ifdef HAVE_AMESOS_TAUCS
  ,false
#endif
#ifdef HAVE_AMESOS_PARDISO
  ,true
#endif
#ifdef HAVE_AMESOS_PASTIX
  ,true
#endif
#ifdef HAVE_AMESOS_PARAKLETE
  ,true
#endif
#ifdef HAVE_AMESOS_MUMPS
  ,true
#endif
#ifdef HAVE_AMESOS_SCALAPACK
  ,true
#endif
#ifdef HAVE_AMESOS_DSCPACK
  ,false
#endif
};

Teuchos::StringToIntMap
Amesos::solverTypeNameToEnumMap(
  "Amesos::SolverType"
  ,Amesos::numSolverTypes
  ,Amesos::solverTypeNames
  );

const Amesos::ERefactorizationPolicy Amesos::refactorizationPolicyValues[Amesos::numRefactorizationPolices] =
{
  Amesos::REPIVOT_ON_REFACTORIZATION
  ,Amesos::NO_PIVOT_ON_REFACTORIZATION
};

const char* Amesos::refactorizationPolicyNames[Amesos::numRefactorizationPolices] =
{
  "RepivotOnRefactorization"
  ,"NoPivotOnRefactorization"
};

Teuchos::StringToIntMap
Amesos::refactorizationPolicyNameToEnumMap(
  "Amesos::RefactorizationPolices"
  ,Amesos::numRefactorizationPolices
  ,Amesos::refactorizationPolicyNames
  );

} // namespace Thyra
