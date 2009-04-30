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
