// @HEADER
// ***********************************************************************
//
//                Amesos2: Direct Sparse Solver Package
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

/**
   \file   Amesos2_Factory_def.hpp
   \author Eric Bavier <etbavier@sandia.gov>
   \date   Thu Jul  8 22:57:26 2010

   \brief  Templated Amesos2 solver Factory. Definitions.
*/


#ifndef AMESOS2_FACTORY_DEF_HPP
#define AMESOS2_FACTORY_DEF_HPP


#include <Teuchos_RCP.hpp>
using Teuchos::rcp;
using Teuchos::RCP;

namespace Amesos {


template <class Matrix,
          class Vector >
RCP<SolverBase>
Factory<Matrix,Vector>::create(Matrix* A, Vector* X, Vector* B)
{
  std::string solver = "Klu2";
  // Pass non-owning RCP objects to other factory method
  return( create(solver, rcp(A,false), rcp(X,false), rcp(B,false)) );
}


template <class Matrix,
          class Vector >
RCP<SolverBase>
Factory<Matrix,Vector>::create(RCP<Matrix> A, RCP<Vector> X, RCP<Vector> B)
{
  std::string solver = "Klu2";
  return( create(solver, A, X, B) );
}


template <class Matrix,
          class Vector >
RCP<SolverBase>
Factory<Matrix,Vector>::create(const char* solverName, Matrix* A, Vector* X, Vector* B)
{
  std::string solver = solverName;
  // Pass non-owning RCP objects to other factory method
  return( create(solver, rcp(A,false), rcp(X,false), rcp(B,false)) );
}


template <class Matrix,
          class Vector >
RCP<SolverBase>
Factory<Matrix,Vector>::create(
  const char* solverName,
  const RCP<Matrix> A,
  const RCP<Vector> X,
  const RCP<Vector> B)
{
  std::string solver = solverName;
  return( create(solver, A, X, B) );
}


template <class Matrix,
          class Vector >
RCP<SolverBase>
Factory<Matrix,Vector>::create(const std::string solverName, Matrix* A, Vector* X, Vector* B)
{
  // Pass non-owning RCP objects to other factory method
  return( create(solverName, rcp(A,false), rcp(X,false), rcp(B,false)) );
}


template <class Matrix,
          class Vector >
RCP<SolverBase>
Factory<Matrix,Vector>::create(
  const std::string solverName,
  const RCP<Matrix> A,
  const RCP<Vector> X,
  const RCP<Vector> B)
{
  // Check for our native solver first.
  // 
  // Remove compiler guards once interface is finalized, since we will always include it.
#ifdef HAVE_AMESOS2_KLU2
  if((solverName == "Amesos2_Klu2") || (solverName == "Klu2") || (solverName == "KLU2")){
    return( rcp(new Klu2<Matrix,Vector>(A, X, B)) );
  }
#endif

#ifdef HAVE_AMESOS2_KLU
  if((solverName == "Amesos2_Klu") || (solverName == "Klu")){
    return( rcp(new Klu<Matrix,Vector>(A, X, B)) );
  }
#endif

#ifdef HAVE_AMESOS2_LAPACK
  if((solverName == "Amesos2_Lapack") || (solverName == "Lapack")){
    return( rcp(new Lapack<Matrix,Vector>(A, X, B)) );
  }
#endif

#ifdef HAVE_AMESOS2_MUMPS
  if((solverName == "Amesos2_Mumps") || (solverName == "Mumps")){
    return( rcp(new Mumps<Matrix,Vector>(A, X, B)) );
  }
#endif

#ifdef HAVE_AMESOS2_SCALAPACK
  if((solverName == "Amesos2_Scalapack") || (solverName == "Scalapack")){
    return( rcp(new Scalapack<Matrix,Vector>(A, X, B)) );
  }
#endif

#ifdef HAVE_AMESOS2_UMFPACK
  if((solverName == "Amesos2_Umfpack") || (solverName == "Umfpack")){
    return( rcp(new Umfpack<Matrix,Vector>(A, X, B)) );
  }
#endif

#ifdef HAVE_AMESOS2_SUPERLUDIST
  if((solverName == "Amesos2_Superludist") || (solverName == "Superludist")){
    return( rcp(new Superludist<Matrix,Vector>(A, X, B)) );
  }
#endif

#ifdef HAVE_AMESOS2_SUPERLU
  if((solverName == "Amesos2_Superlu") || (solverName == "Superlu")){
    return( rcp(new Superlu<Matrix,Vector>(A, X, B)) );
  }
#endif

#ifdef HAVE_AMESOS2_DSCPACK
  if((solverName == "Amesos2_Dscpack") || (solverName == "Dscpack")){
    return( rcp(new Dscpack<Matrix,Vector>(A, X, B)) );
  }
#endif

#ifdef HAVE_AMESOS2_PARDISO
  if((solverName == "Amesos2_Pardiso") || (solverName == "Pardiso")){
    return( rcp(new Pardiso<Matrix,Vector>(A, X, B)) );
  }
#endif

#ifdef HAVE_AMESOS2_TAUCS
  if((solverName == "Amesos2_Taucs") || (solverName == "Taucs")){
    return( rcp(new Taucs<Matrix,Vector>(A, X, B)) );
  }
#endif

#ifdef HAVE_AMESOS2_PARAKLETE
  if((solverName == "Amesos2_Paraklete") || (solverName == "Paraklete")){
    return( rcp(new Paraklete<Matrix,Vector>(A, X, B)) );
  }
#endif

  /* If none of the above conditionals are satisfied, then the solver
   * requested is not yet supported.  We throw a runtime exception stating
   * this, and return null.
   */
  std::string err_msg = solverName + " is not implemented";
  TEST_FOR_EXCEPTION(true, std::invalid_argument, err_msg);
  return( Teuchos::null );
}


template <class Matrix,
          class Vector >
bool Factory<Matrix,Vector>::query(const char* solverName){
  std::string solver = solverName;
  return( query(solver) );
}


template <class Matrix,
          class Vector >
bool Factory<Matrix,Vector>::query(const std::string solverName){
#ifdef HAVE_AMESOS2_KLU
  if((solverName == "Amesos2_Klu") || (solverName == "Klu")){
    return( true );
  }
#endif

#ifdef HAVE_AMESOS2_LAPACK
  if((solverName == "Amesos2_Lapack") || (solverName == "Lapack")){
    return( true );
  }
#endif

#ifdef HAVE_AMESOS2_MUMPS
  if((solverName == "Amesos2_Mumps") || (solverName == "Mumps")){
    return( true );
  }
#endif

#ifdef HAVE_AMESOS2_SCALAPACK
  if((solverName == "Amesos2_Scalapack") || (solverName == "Scalapack")){
    return( true );
  }
#endif

#ifdef HAVE_AMESOS2_UMFPACK
  if((solverName == "Amesos2_Umfpack") || (solverName == "Umfpack")){
    return( true );
  }
#endif

#ifdef HAVE_AMESOS2_SUPERLUDIST
  if((solverName == "Amesos2_Superludist") || (solverName == "Superludist")){
    return( true );
  }
#endif

#ifdef HAVE_AMESOS2_SUPERLU
  if((solverName == "Amesos2_Superlu") || (solverName == "Superlu")){
    return( true );
  }
#endif

#ifdef HAVE_AMESOS2_DSCPACK
  if((solverName == "Amesos2_Dscpack") || (solverName == "Dscpack")){
    return( true );
  }
#endif

#ifdef HAVE_AMESOS2_PARDISO
  if((solverName == "Amesos2_Pardiso") || (solverName == "Pardiso")){
    return( true );
  }
#endif

#ifdef HAVE_AMESOS2_SUPERLU
  if((solverName == "Amesos2_Taucs") || (solverName == "Taucs")){
    return( true );
  }
#endif

#ifdef HAVE_AMESOS2_PARAKLETE
  if((solverName == "Amesos2_Paraklete") || (solverName == "Paraklete")){
    return( true );
  }
#endif

  // Otherwise, the solver is not available
  return( false );
}

// TODO: Here in Amesos.cpp there is a function defined
// getValidParameters.  I wonder if it would be more appropriate
// to define this function in the base Solver and concrete Solver
// classes instead.


} // end namespace Amesos

#endif	// AMESOS2_FACTORY_DEF_HPP
