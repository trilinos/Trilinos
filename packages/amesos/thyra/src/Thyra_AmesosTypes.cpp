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

const Amesos::ESolverType Amesos::SolverTypeValues[Amesos::numSolverTypes] =
{
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

const char* Amesos::SolverTypeNames[Amesos::numSolverTypes] =
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

const Amesos::ERefactorizationPolicy Amesos::RefactorizationPolicyValues[Amesos::numRefactorizationPolices] =
{
  REPIVOT_ON_REFACTORIZATION
  ,NO_PIVOT_ON_REFACTORIZATION
};

const char* Amesos::RefactorizationPolicyNames[Amesos::numRefactorizationPolices] =
{
  "RepivotOnRefactorization"
  ,"NoPivotOnRefactorization"
};

} // namespace Thyra
