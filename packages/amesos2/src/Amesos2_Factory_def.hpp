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

   \brief  Templated Amesos2 solver factory methods. Definitions.
*/


#ifndef AMESOS2_FACTORY_DEF_HPP
#define AMESOS2_FACTORY_DEF_HPP

#include "Amesos2_config.h"

#include <Teuchos_RCP.hpp>

namespace Amesos2 {

  using Teuchos::rcp;
  using Teuchos::RCP;

  /*
   * Utility function to transform a string into all lowercase
   */
  std::string tolower(const std::string& s);


  template <class Matrix,
	    class Vector >
  Solver<Matrix,Vector>*
  create(Matrix* A, Vector* X, Vector* B)
  {
    std::string solver = "Klu2";
    // Pass non-owning RCP objects to other factory method
    return( create(solver, rcp(A,false), rcp(X,false), rcp(B,false)).getRawPtr() );
  }


  template <class Matrix,
	    class Vector >
  RCP<Solver<Matrix,Vector> >
  create(RCP<const Matrix> A, RCP<Vector> X, RCP<const Vector> B)
  {
    std::string solver = "Klu2";
    return( create(solver, A, X, B) );
  }


  template <class Matrix,
	    class Vector >
  Solver<Matrix,Vector>*
  create(const char* solverName, const Matrix* A, Vector* X, const Vector* B)
  {
    std::string solver = solverName;
    // Pass non-owning RCP objects to other factory method
    return( create(solver, rcp(A,false), rcp(X,false), rcp(B,false)).getRawPtr() );
  }


  template <class Matrix,
	    class Vector >
  RCP<Solver<Matrix,Vector> >
  create(const char* solverName,
	 const RCP<const Matrix> A, const RCP<Vector> X, const RCP<const Vector> B)
  {
    std::string solver = solverName;
    return( create(solver, A, X, B) );
  }


  template <class Matrix,
	    class Vector >
  Solver<Matrix,Vector>*
  create(const std::string solverName, const Matrix* A){
    return( create(solverName, rcp(A,false), RCP<Vector>(), RCP<Vector>()).getRawPtr() );
  }


  template <class Matrix,
	    class Vector >
  RCP<Solver<Matrix,Vector> >
  create(const std::string solverName, const RCP<const Matrix> A){
    return( create(solverName, A, RCP<Vector>(), RCP<Vector>()) );
  }


  template <class Matrix,
	    class Vector >
  RCP<Solver<Matrix,Vector> >
  create(const std::string solverName, const Matrix* A, Vector* X, const Vector* B)
  {
    // Pass non-owning RCP objects to other factory method
    return( create(solverName, rcp(A,false), rcp(X,false), rcp(B,false)) );
  }

  
  template <class Matrix,
	    class Vector >
  RCP<Solver<Matrix,Vector> >
  create(const std::string solver_name,
	 const RCP<const Matrix> A, const RCP<Vector> X, const RCP<const Vector> B)
  {
    std::string solverName = tolower(solver_name); // for easy string checking

    // Check for our native solver first.  Treat KLU and KLU2 as equals.
    // 
    // We use compiler guards in case a user does want to disable KLU2
#ifdef HAVE_AMESOS2_KLU2
    if((solverName == "amesos2_klu2") || (solverName == "klu2") ||
       (solverName == "amesos2_klu")  || (solverName == "klu")){
      return handle_solver_type_support<Klu2,Matrix,Vector>::apply(A, X, B);
      // return( rcp(new Klu2<Matrix,Vector>(A, X, B)) );
    }
#endif

#ifdef HAVE_AMESOS2_MUMPS
    if((solverName == "amesos2_mumps") || (solverName == "mumps")){
      return handle_solver_type_support<Mumps,Matrix,Vector>::apply(A, X, B);
      // return( rcp(new Mumps<Matrix,Vector>(A, X, B)) );
    }
#endif

#ifdef HAVE_AMESOS2_UMFPACK
    if((solverName == "amesos2_umfpack") || (solverName == "umfpack")){
      return handle_solver_type_support<Umfpack,Matrix,Vector>::apply(A, X, B);
      // return( rcp(new Umfpack<Matrix,Vector>(A, X, B)) );
    }
#endif

#ifdef HAVE_AMESOS2_SUPERLUDIST
    if((solverName == "amesos2_superludist") ||
       (solverName == "superludist") ||
       (solverName == "amesos2_superlu_dist") ||
       (solverName == "superlu_dist")){
      return handle_solver_type_support<Superludist,Matrix,Vector>::apply(A, X, B);
      // return( rcp(new Superludist<Matrix,Vector>(A, X, B)) );
    }
#endif

#ifdef HAVE_AMESOS2_SUPERLUMT
    if((solverName == "amesos2_superlumt") ||
       (solverName == "superlumt") ||
       (solverName == "amesos2_superlu_mt") ||
       (solverName == "superlu_mt")){
      return handle_solver_type_support<Superlumt,Matrix,Vector>::apply(A, X, B);
      // return( rcp(new Superlumt<Matrix,Vector>(A, X, B)) );
    }
#endif

#ifdef HAVE_AMESOS2_SUPERLU
    if((solverName == "amesos2_superlu") ||
       (solverName == "superlu")){
      return handle_solver_type_support<Superlu,Matrix,Vector>::apply(A, X, B);
      // return( rcp(new Superlu<Matrix,Vector>(A, X, B)) );
    }
#endif

#ifdef HAVE_AMESOS2_DSCPACK
    if((solverName == "amesos2_dscpack") || (solverName == "dscpack")){
      return handle_solver_type_support<Dscpack,Matrix,Vector>::apply(A, X, B);
      // return( rcp(new Dscpack<Matrix,Vector>(A, X, B)) );
    }
#endif

#ifdef HAVE_AMESOS2_PARDISO
    if((solverName == "amesos2_pardiso") || (solverName == "pardiso")){
      return handle_solver_type_support<Pardiso,Matrix,Vector>::apply(A, X, B);
      // return( rcp(new Pardiso<Matrix,Vector>(A, X, B)) );
    }
#endif

#ifdef HAVE_AMESOS2_TAUCS
    if((solverName == "amesos2_taucs") || (solverName == "taucs")){
      return handle_solver_type_support<Taucs,Matrix,Vector>::apply(A, X, B);
      // return( rcp(new Taucs<Matrix,Vector>(A, X, B)) );
    }
#endif

    /* If none of the above conditionals are satisfied, then the solver
     * requested is not yet supported.  We throw a runtime exception stating
     * this, and return null.
     */
    std::string err_msg = solver_name + " is not enabled or is not supported";
    TEST_FOR_EXCEPTION(true, std::invalid_argument, err_msg);
    return( Teuchos::null );
  }


  /**********************
   *   QUERY function   *
   **********************/

  bool query(const char* solverName){
    std::string solver = solverName;
    return( query(solver) );
  }

  bool query(const std::string solver_name){
    std::string solverName = tolower(solver_name); // for easier string checking
#ifdef HAVE_AMESOS2_KLU2
    if((solverName == "amesos2_klu2") || (solverName == "klu2") ||
       (solverName == "amesos2_klu")  || (solverName == "klu")){
      return( true );
    }
#endif

#ifdef HAVE_AMESOS2_KLU
    if((solverName == "amesos2_klu") || (solverName == "klu")){
      return( true );
    }
#endif

#ifdef HAVE_AMESOS2_MUMPS
    if((solverName == "amesos2_mumps") || (solverName == "mumps")){
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
} // end namespace Amesos2

#endif  // AMESOS2_FACTORY_DEF_HPP
