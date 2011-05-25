// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package 
//                  Copyright 2010 Sandia Corporation
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
// ***********************************************************************
//
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

  /*
   * Utility function to transform a string into all lowercase
   */
  std::string tolower(const std::string& s);


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
  Factory<Matrix,Vector>::create(const char* solverName,
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
  Factory<Matrix,Vector>::create(const std::string solver_name,
				 const RCP<Matrix> A,
				 const RCP<Vector> X,
				 const RCP<Vector> B)
  {
    std::string solverName = tolower(solver_name); // for easy string checking
    // Check for our native solver first.
    // 
    // Remove compiler guards once interface is finalized, since we will always include it?
#ifdef HAVE_AMESOS2_KLU2
    if((solverName == "amesos2_klu2") || (solverName == "klu2")){
      return( rcp(new Klu2<Matrix,Vector>(A, X, B)) );
    }
#endif

#ifdef HAVE_AMESOS2_KLU
    if((solverName == "amesos2_klu") || (solverName == "klu")){
      return( rcp(new Klu<Matrix,Vector>(A, X, B)) );
    }
#endif

#ifdef HAVE_AMESOS2_LAPACK
    if((solverName == "amesos2_lapack") || (solverName == "lapack")){
      return( rcp(new Lapack<Matrix,Vector>(A, X, B)) );
    }
#endif

#ifdef HAVE_AMESOS2_MUMPS
    if((solverName == "amesos2_mumps") || (solverName == "mumps")){
      return( rcp(new Mumps<Matrix,Vector>(A, X, B)) );
    }
#endif

#ifdef HAVE_AMESOS2_SCALAPACK
    if((solverName == "amesos2_scalapack") || (solverName == "scalapack")){
      return( rcp(new Scalapack<Matrix,Vector>(A, X, B)) );
    }
#endif

#ifdef HAVE_AMESOS2_UMFPACK
    if((solverName == "amesos2_umfpack") || (solverName == "umfpack")){
      return( rcp(new Umfpack<Matrix,Vector>(A, X, B)) );
    }
#endif

#ifdef HAVE_AMESOS2_SUPERLUDIST
    if((solverName == "amesos2_superludist") ||
       (solverName == "superludist") ||
       (solverName == "amesos2_superlu_dist") ||
       (solverName == "superlu_dist")){
      return( rcp(new Superludist<Matrix,Vector>(A, X, B)) );
    }
#endif

#ifdef HAVE_AMESOS2_SUPERLUMT
    if((solverName == "amesos2_superlumt") ||
       (solverName == "superlumt") ||
       (solverName == "amesos2_superlu_mt") ||
       (solverName == "superlu_mt")){
      return( rcp(new Superlumt<Matrix,Vector>(A, X, B)) );
    }
#endif

#ifdef HAVE_AMESOS2_SUPERLU
    if((solverName == "amesos2_superlu") ||
       (solverName == "superlu")){
      return( rcp(new Superlu<Matrix,Vector>(A, X, B)) );
    }
#endif

#ifdef HAVE_AMESOS2_DSCPACK
    if((solverName == "amesos2_dscpack") || (solverName == "dscpack")){
      return( rcp(new Dscpack<Matrix,Vector>(A, X, B)) );
    }
#endif

#ifdef HAVE_AMESOS2_PARDISO
    if((solverName == "amesos2_pardiso") || (solverName == "pardiso")){
      return( rcp(new Pardiso<Matrix,Vector>(A, X, B)) );
    }
#endif

#ifdef HAVE_AMESOS2_TAUCS
    if((solverName == "amesos2_taucs") || (solverName == "taucs")){
      return( rcp(new Taucs<Matrix,Vector>(A, X, B)) );
    }
#endif

#ifdef HAVE_AMESOS2_PARAKLETE
    if((solverName == "amesos2_paraklete") || (solverName == "paraklete")){
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
#ifdef HAVE_AMESOS2_KLU2
    if((solverName == "amesos2_klu2") || (solverName == "klu2")){
      return( true );
    }
#endif

#ifdef HAVE_AMESOS2_KLU
    if((solverName == "amesos2_klu") || (solverName == "klu")){
      return( true );
    }
#endif

#ifdef HAVE_AMESOS2_LAPACK
    if((solverName == "amesos2_lapack") || (solverName == "lapack")){
      return( true );
    }
#endif

#ifdef HAVE_AMESOS2_MUMPS
    if((solverName == "amesos2_mumps") || (solverName == "mumps")){
      return( true );
    }
#endif

#ifdef HAVE_AMESOS2_SCALAPACK
    if((solverName == "amesos2_scalapack") || (solverName == "scalapack")){
      return( true );
    }
#endif

#ifdef HAVE_AMESOS2_UMFPACK
    if((solverName == "amesos2_umfpack") || (solverName == "umfpack")){
      return( true );
    }
#endif

#ifdef HAVE_AMESOS2_SUPERLUDIST
    if((solverName == "amesos2_superludist") ||
       (solverName == "superludist") ||
       (solverName == "amesos2_superlu_dist") ||
       (solverName == "superlu_dist")){
      return( true );
    }
#endif

#ifdef HAVE_AMESOS2_SUPERLUMT
    if((solverName == "amesos2_superlumt") ||
       (solverName == "superlumt") ||
       (solverName == "amesos2_superlu_mt") ||
       (solverName == "superlu_mt")){
      return( true );
    }
#endif

#ifdef HAVE_AMESOS2_SUPERLU
    if((solverName == "amesos2_superlu") ||
       (solverName == "superlu")){
      return( true );
    }
#endif

#ifdef HAVE_AMESOS2_DSCPACK
    if((solverName == "amesos2_dscpack") || (solverName == "dscpack")){
      return( true );
    }
#endif

#ifdef HAVE_AMESOS2_PARDISO
    if((solverName == "amesos2_pardiso") || (solverName == "pardiso")){
      return( true );
    }
#endif

#ifdef HAVE_AMESOS2_TAUCS
    if((solverName == "amesos2_taucs") || (solverName == "taucs")){
      return( true );
    }
#endif

#ifdef HAVE_AMESOS2_PARAKLETE
    if((solverName == "amesos2_paraklete") || (solverName == "paraklete")){
      return( true );
    }
#endif

    // Otherwise, the solver is not available
    return( false );
  }


  std::string tolower(const std::string& s)
  {
    std::locale loc;
    std::string rtn = s;
    for (size_t i=0; i<rtn.length(); ++i)
      {
	rtn[i] = tolower(rtn[i],loc);
      }
    return rtn;
  }


  // Here in Amesos.cpp there is a function defined getValidParameters.
  // We decided it would be best to have such functionality be part of
  // the Solver class, with some functionality being delegated to the
  // solver implementations


} // end namespace Amesos

#endif	// AMESOS2_FACTORY_DEF_HPP
