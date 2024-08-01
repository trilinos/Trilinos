// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
