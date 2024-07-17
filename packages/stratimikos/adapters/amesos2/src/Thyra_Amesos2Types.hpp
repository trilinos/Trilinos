// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_AMESOS2_TYPES_HPP
#define THYRA_AMESOS2_TYPES_HPP

#include "Amesos2_config.h"
#include "Teuchos_StringToIntMap.hpp"

namespace Thyra {

namespace Amesos2 {

/** \brief .
\ingroup Amesos2_Thyra_adapters_grp
*/
enum ESolverType {
  KLU2
#ifdef HAVE_AMESOS2_LAPACK
  ,LAPACK
#endif
#ifdef HAVE_AMESOS2_SUPERLU
  ,SUPERLU
#endif
#ifdef HAVE_AMESOS2_SUPERLUMT
  ,SUPERLUMT
#endif
#ifdef HAVE_AMESOS2_SUPERLUDIST
  ,SUPERLUDIST
#endif
#ifdef HAVE_AMESOS2_PARDISO_MKL
  ,PARDISO_MKL
#endif
#ifdef HAVE_AMESOS2_CHOLMOD
  ,CHOLMOD
#endif
#ifdef HAVE_AMESOS2_BASKER
  ,BASKER
#endif
#ifdef HAVE_AMESOS2_MUMPS
  ,MUMPS
#endif
}; 

/** \brief .
\ingroup Amesos2_Thyra_adapters_grp
*/
const int numSolverTypes = 1
#ifdef HAVE_AMESOS2_LAPACK
+1
#endif
#ifdef HAVE_AMESOS2_SUPERLU
+1
#endif
#ifdef HAVE_AMESOS2_SUPERLUMT
+1
#endif
#ifdef HAVE_AMESOS2_SUPERLUDIST
+1
#endif
#ifdef HAVE_AMESOS2_PARDISO_MKL
+1
#endif
#ifdef HAVE_AMESOS2_CHOLMOD
+1
#endif
#ifdef HAVE_AMESOS2_BASKER
+1
#endif
#ifdef HAVE_AMESOS2_MUMPS
+1
#endif
; 

/** \brief .
\ingroup Amesos2_Thyra_adapters_grp
*/
extern const ESolverType solverTypeValues[numSolverTypes];

/** \brief .
\ingroup Amesos2_Thyra_adapters_grp
*/
extern const char* solverTypeNames[numSolverTypes];

/** \brief .
\ingroup Amesos2_Thyra_adapters_grp
*/
extern const bool supportsUnsymmetric[numSolverTypes];

/** \brief .
\ingroup Amesos2_Thyra_adapters_grp
*/
inline const char* toString(const ESolverType solverType)
{ return solverTypeNames[solverType]; }

/** \brief .
\ingroup Amesos2_Thyra_adapters_grp
*/
extern Teuchos::StringToIntMap solverTypeNameToEnumMap;

/** \brief The policy used on refactoring a matrix.
\ingroup Amesos2_Thyra_adapters_grp
*/
enum ERefactorizationPolicy {
  REPIVOT_ON_REFACTORIZATION     ///< Completely new pivoting will be used on refactorizations!
  ,NO_PIVOT_ON_REFACTORIZATION   ///< No piviting, or only minor repivoting, will be used on refactorizations!
};

/** \brief .
\ingroup Amesos2_Thyra_adapters_grp
*/
const int numRefactorizationPolices = 2;

/** \brief .
\ingroup Amesos2_Thyra_adapters_grp
*/
extern const ERefactorizationPolicy refactorizationPolicyValues[numRefactorizationPolices];

/** \brief . */
extern const char* refactorizationPolicyNames[numRefactorizationPolices];

/** \brief .
\ingroup Amesos2_Thyra_adapters_grp
*/
inline const char* toString(const ERefactorizationPolicy refactorizationPolicy)
{ return refactorizationPolicyNames[refactorizationPolicy]; }

/** \brief .
\ingroup Amesos2_Thyra_adapters_grp
*/
extern Teuchos::StringToIntMap refactorizationPolicyNameToEnumMap;

}// namespace Amesos2

} // namespace Thyra

#endif // THYRA_AMESOS2_TYPES_HPP
