/*
// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#ifndef THYRA_AMESOS_TYPES_HPP
#define THYRA_AMESOS_TYPES_HPP

#include "Amesos_ConfigDefs.h"
#include "Teuchos_StringToIntMap.hpp"

namespace Thyra {

namespace Amesos {

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
enum ESolverType {
  LAPACK
#ifdef HAVE_AMESOS_KLU
  ,KLU
#endif
#ifdef HAVE_AMESOS_UMFPACK
  ,UMFPACK
#endif
#ifdef HAVE_AMESOS_SUPERLU
  ,SUPERLU
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
  ,SUPERLUDIST
#endif
#ifdef HAVE_AMESOS_TAUCS
  ,TAUCS
#endif
#ifdef HAVE_AMESOS_PARDISO
  ,PARDISO
#endif
#ifdef HAVE_AMESOS_PASTIX
  ,PASTIX
#endif
#ifdef HAVE_AMESOS_PARAKLETE
  ,PARAKLETE
#endif
#ifdef HAVE_AMESOS_MUMPS
  ,MUMPS
#endif
#ifdef HAVE_AMESOS_SCALAPACK
  ,SCALAPACK
#endif
#ifdef HAVE_AMESOS_DSCPACK
  ,DSCPACK
#endif
}; 

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
const int numSolverTypes = 1 // LAPACK 
#ifdef HAVE_AMESOS_KLU
+1
#endif
#ifdef HAVE_AMESOS_UMFPACK
+1
#endif
#ifdef HAVE_AMESOS_SUPERLU
+1
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
+1
#endif
#ifdef HAVE_AMESOS_TAUCS
+1
#endif
#ifdef HAVE_AMESOS_PARDISO
+1
#endif
#ifdef HAVE_AMESOS_PASTIX
+1
#endif
#ifdef HAVE_AMESOS_PARAKLETE
+1
#endif
#ifdef HAVE_AMESOS_MUMPS
+1
#endif
#ifdef HAVE_AMESOS_SCALAPACK
+1
#endif
#ifdef HAVE_AMESOS_DSCPACK
+1
#endif
;

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
extern const ESolverType solverTypeValues[numSolverTypes];

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
extern const char* solverTypeNames[numSolverTypes];

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
extern const bool supportsUnsymmetric[numSolverTypes];

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
inline const char* toString(const ESolverType solverType)
{ return solverTypeNames[solverType]; }

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
extern Teuchos::StringToIntMap solverTypeNameToEnumMap;

/** \brief The policy used on refactoring a matrix.
\ingroup Amesos_Thyra_adapters_grp
*/
enum ERefactorizationPolicy {
  REPIVOT_ON_REFACTORIZATION     ///< Completely new pivoting will be used on refactorizations!
  ,NO_PIVOT_ON_REFACTORIZATION   ///< No piviting, or only minor repivoting, will be used on refactorizations!
};

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
const int numRefactorizationPolices = 2;

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
extern const ERefactorizationPolicy refactorizationPolicyValues[numRefactorizationPolices];

/** \brief . */
extern const char* refactorizationPolicyNames[numRefactorizationPolices];

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
inline const char* toString(const ERefactorizationPolicy refactorizationPolicy)
{ return refactorizationPolicyNames[refactorizationPolicy]; }

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
extern Teuchos::StringToIntMap refactorizationPolicyNameToEnumMap;

} // namespace Amesos

} // namespace Thyra

#endif // THYRA_AMESOS_TYPES_HPP
