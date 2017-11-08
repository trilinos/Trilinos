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

#include "Thyra_Amesos2Types.hpp"

namespace Thyra {

const Amesos2::ESolverType Amesos2::solverTypeValues[Amesos2::numSolverTypes] =
{
  Amesos2::KLU2
#ifdef HAVE_AMESOS2_LAPACK
  ,Amesos2::LAPACK
#endif
#ifdef HAVE_AMESOS2_SUPERLU
  ,Amesos2::SUPERLU
#endif
#ifdef HAVE_AMESOS2_SUPERLUMT
  ,Amesos2::SUPERLUMT
#endif
#ifdef HAVE_AMESOS2_SUPERLUDIST
  ,Amesos2::SUPERLUDIST
#endif
#ifdef HAVE_AMESOS2_PARDISO_MKL
  ,Amesos2::PARDISO_MKL
#endif
#ifdef HAVE_AMESOS2_CHOLMOD
  ,Amesos2::CHOLMOD
#endif
#ifdef HAVE_AMESOS2_BASKER
  ,Amesos2::BASKER
#endif
#ifdef HAVE_AMESOS2_MUMPS
  ,Amesos2::MUMPS
#endif
};

const char* Amesos2::solverTypeNames[Amesos2::numSolverTypes] =
{
  "KLU2"
#ifdef HAVE_AMESOS2_LAPACK
  ,"LAPACK"
#endif
#ifdef HAVE_AMESOS2_SUPERLU
  ,"SuperLU"
#endif
#ifdef HAVE_AMESOS2_SUPERLUMT
  ,"SuperLU_MT"
#endif
#ifdef HAVE_AMESOS2_SUPERLUDIST
  ,"SuperLU_DIST"
#endif
#ifdef HAVE_AMESOS2_PARDISO_MKL
  ,"PARDISOMKL"
#endif
#ifdef HAVE_AMESOS2_CHOLMOD
  ,"Cholmod"
#endif
#ifdef HAVE_AMESOS2_BASKER
  ,"Basker"
#endif
#ifdef HAVE_AMESOS2_MUMPS
  ,"MUMPS"
#endif
};

const bool Amesos2::supportsUnsymmetric[Amesos2::numSolverTypes] =
{
  true
#ifdef HAVE_AMESOS2_LAPACK
  ,true
#endif
#ifdef HAVE_AMESOS2_SUPERLU
  ,true
#endif
#ifdef HAVE_AMESOS2_SUPERLUMT
  ,true
#endif
#ifdef HAVE_AMESOS2_SUPERLUDIST
  ,true
#endif
#ifdef HAVE_AMESOS2_PARDISO_MKL
  ,true
#endif
#ifdef HAVE_AMESOS2_CHOLMOD
  ,false //don't know, being conservative
#endif
#ifdef HAVE_AMESOS2_BASKER
  ,false //don't know, being conservative
#endif
#ifdef HAVE_AMESOS2_MUMPS
  ,true
#endif
};

Teuchos::StringToIntMap
Amesos2::solverTypeNameToEnumMap(
  "Amesos2::SolverType"
  ,Amesos2::numSolverTypes
  ,Amesos2::solverTypeNames
  );

const Amesos2::ERefactorizationPolicy Amesos2::refactorizationPolicyValues[Amesos2::numRefactorizationPolices] =
{
  Amesos2::REPIVOT_ON_REFACTORIZATION
  ,Amesos2::NO_PIVOT_ON_REFACTORIZATION
};

const char* Amesos2::refactorizationPolicyNames[Amesos2::numRefactorizationPolices] =
{
  "RepivotOnRefactorization"
  ,"NoPivotOnRefactorization"
};

Teuchos::StringToIntMap
Amesos2::refactorizationPolicyNameToEnumMap(
  "Amesos2::RefactorizationPolices"
  ,Amesos2::numRefactorizationPolices
  ,Amesos2::refactorizationPolicyNames
  );

} // namespace Thyra
