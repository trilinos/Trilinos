// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
 * \file   Amesos2_Factory.hpp
 * \author Eric Bavier <etbavie@sandia.gov>
 * \date   Tue Jun  7 16:25:16 2011
 *
 * \brief  Contains declarations for Amesos2::create and Amesos2::query.
 *
 * Amesos2 defines the nonmember factory method \c Amesos2::create for
 * creating instances of Amesos2 solvers.  If a users asks
 * Amesos2::create to create a solver with a matrix whose scalar type
 * is not supported by that solver, then a runtime
 * std::invalid_argument exception will be thrown.
 *
 * The \c Amesos2::query function can be used to ask Amesos2 at
 * runtime whether a particular solver is supported.
 *
 * \attention
 * Users should favor these factory methods for creating Amesos2 solver
 * instances over explicitly instantiating their own.
 *
 * \note A solver's third-party library must be enabled in the
 * Trilinos build, and Amesos2 must be also told to enable it.  Put
 * <tt>Amesos2_ENABLE_<i>SOLVERNAME</i>:BOOL=ON</tt> in your Trilinos
 * configuration script to do this, where <i>SOLVERNAME</i> is the
 * name of the solver you would like to enable.
 *
 * \section usage Example Usage
 *
 * \code
 * typedef Tpetra::CrsMatrix<double,int> MAT;
 * typedef Tpetra::MultiVector<double,int> VEC;
 * // ... Create A of type RCP<MAT>, and X and B of type RCP<VEC> ...
 * RCP<Amesos2::Solver<MAT,VEC> > solver = Amesos2::create<MAT,VEC>("SuperLU", A, X, B);
 * \endcode
 *
 * \ingroup amesos2_solver_framework
 */

#ifndef AMESOS2_FACTORY_HPP
#define AMESOS2_FACTORY_HPP

#include "Amesos2_config.h"

#include "Amesos2_Solver.hpp"
#include "Amesos2_SolverTraits.hpp"

#include "Teuchos_ScalarTraits.hpp"
#include "Amesos2_MultiVecAdapter.hpp"
#include "Amesos2_MatrixTraits.hpp"
#include "Amesos2_ctassert.hpp"

#ifdef HAVE_AMESOS2_BASKER
#include "Amesos2_Basker.hpp"
#endif

#ifdef HAVE_AMESOS2_SHYLU_NODEBASKER
#include "Amesos2_ShyLUBasker.hpp"
#endif

#if defined(HAVE_AMESOS2_KLU2)
#include "Amesos2_KLU2.hpp"
#endif

#ifdef HAVE_AMESOS2_SUPERLUDIST // Distributed-memory SuperLU
#include "Amesos2_Superludist.hpp"
#endif

#ifdef HAVE_AMESOS2_SUPERLUMT   // Multi-threaded SuperLU
#include "Amesos2_Superlumt.hpp"
#endif

#ifdef HAVE_AMESOS2_UMFPACK     // Umfpack
#include "Amesos2_Umfpack.hpp"
#endif

#ifdef HAVE_AMESOS2_SHYLU_NODETACHO       // Tacho
#include "Amesos2_Tacho.hpp"
#endif

#ifdef HAVE_AMESOS2_SUPERLU     // Sequential SuperLU
#include "Amesos2_Superlu.hpp"
#endif

#ifdef HAVE_AMESOS2_PARDISO_MKL // MKL version of Pardiso
#include "Amesos2_PardisoMKL.hpp"
#endif

#ifdef HAVE_AMESOS2_CSS_MKL // Cluster-Sparse solver from MKL
#include "Amesos2_CssMKL.hpp"
#endif

#ifdef HAVE_AMESOS2_LAPACK
#include "Amesos2_Lapack.hpp"
#endif

#if defined (HAVE_AMESOS2_CHOLMOD) && defined (HAVE_AMESOS2_EXPERIMENTAL)
#include "Amesos2_Cholmod.hpp"
#endif

#if defined (HAVE_AMESOS2_CUSOLVER) && defined (HAVE_AMESOS2_CUSPARSE)
#include "Amesos2_cuSOLVER.hpp"
#endif

#ifdef HAVE_AMESOS2_MUMPS
#include "Amesos2_MUMPS.hpp"
#endif

#ifdef HAVE_AMESOS2_STRUMPACK
#include "Amesos2_STRUMPACK.hpp"
#endif


namespace Amesos2 {

  template <class,class> class Solver;

  /*
   * Utility function to transform a string into all lowercase
   */
  std::string tolower(const std::string& s);


  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
   * and RHS vector B.
   *
   * This is the default option which interfaces with the native KLU2 solver.
   *
   * \param [in] A pointer to a matrix of coefficients
   * \param [in] X pointer to LHS solution vector
   * \param [in] B pointer to RHS vector
   *
   * \return A \c C pointer to a KLU2 solver interface.
   *
   * \relatesalso Amesos2::Solver
   */
  template < class Matrix,
             class Vector >
  Solver<Matrix,Vector>*
  create(const Matrix* A, Vector* X, const Vector* B);


  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
   * and RHS vector B.
   *
   * This is the default option which interfaces with the native KLU2 solver.
   *
   * \param [in] A <tt>Teuchos::RCP</tt> to the coefficient matrix
   * \param [in] X <tt>Teuchos::RCP</tt> to LHS solution vector
   * \param [in] B <tt>Teuchos::RCP</tt> to RHS vector
   *
   * \return A <tt>Teuchos::RCP</tt> to a KLU2 solver interface.
   *
   * \relatesalso Amesos2::Solver
   */
  template < class Matrix,
             class Vector >
  Teuchos::RCP<Solver<Matrix,Vector> >
  create(Teuchos::RCP<const Matrix> A,
         Teuchos::RCP<Vector>       X,
         Teuchos::RCP<const Vector> B);


  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
   * and RHS vector B.
   *
   * \param [in] solverName A \c C character string with the name of the
   *                        underlying third-party solver desired.
   * \param [in] A pointer to a matrix of coefficients
   * \param [in] X pointer to LHS solution vector
   * \param [in] B pointer to RHS vector
   *
   * \return A \c C pointer to an Amesos2 solver interface.
   *
   * \throw std::invalid_argument The third-party solver named by \c
   * solverName is not supported.
   *
   * \relatesalso Amesos2::Solver
   */
  template < class Matrix,
             class Vector >
  Solver<Matrix,Vector>*
  create(const char* solverName, const Matrix* A, Vector* X, const Vector* B);


  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
   * and RHS vector B.
   *
   * \param [in] solverName The name of the desired third-party solver
   * \param [in] A <tt>Teuchos::RCP</tt> to the coefficient matrix
   * \param [in] X <tt>Teuchos::RCP</tt> to LHS solution vector
   * \param [in] B <tt>Teuchos::RCP</tt> to RHS vector
   *
   * \return A <tt>Teuchos::RCP</tt> to an Amesos2 solver interface.
   *
   * \throw std::invalid_argument The third-party solver named by \c
   * solverName is not supported.
   *
   * \relatesalso Amesos2::Solver
   */
  template < class Matrix,
             class Vector >
  Teuchos::RCP<Solver<Matrix,Vector> >
  create(const char* solverName,
         const Teuchos::RCP<const Matrix> A,
         const Teuchos::RCP<Vector>       X,
         const Teuchos::RCP<const Vector> B);


  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
   * and RHS vector B.
   *
   * \param [in] solverName The name of the desired third-party solver
   * \param [in] A \c C pointer to the coefficient matrix
   * \param [in] X \c C pointer to LHS solution vector
   * \param [in] B \c C pointer to RHS vector
   *
   * \return A \c C pointer to an Amesos2 solver interface.
   *
   * \throw std::invalid_argument The third-party solver named by \c
   * solverName is not supported.
   *
   * \relatesalso Amesos2::Solver
   */
  template < class Matrix,
             class Vector >
  Solver<Matrix,Vector>*
  create(const std::string& solverName, const Matrix* A, Vector* X, const Vector* B);


  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
   * and RHS vector B.
   *
   * \param [in] solverName The name of the desired third-party solver
   * \param [in] A <tt>Teuchos::RCP</tt> to the coefficient matrix
   * \param [in] X <tt>Teuchos::RCP</tt> to LHS solution vector
   * \param [in] B <tt>Teuchos::RCP</tt> to RHS vector
   *
   * \return A <tt>Teuchos::RCP</tt> to an Amesos2 solver interface.
   *
   * \throw std::invalid_argument The third-party solver named by \c
   * solverName is not supported.
   *
   * \relatesalso Amesos2::Solver
   */
  template < class Matrix,
             class Vector >
  Teuchos::RCP<Solver<Matrix,Vector> >
  create(const std::string& solverName,
         const Teuchos::RCP<const Matrix> A,
         const Teuchos::RCP<Vector>       X,
         const Teuchos::RCP<const Vector> B);


  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A.
   *
   * Suitable for cases where numeric factorization must be performed
   * before the X and B vectors are known.  Before a solve, the \c
   * setX() and \c setB() functions should be used to set X and B, or
   * the overloaded solve(X,B) method should be used.
   *
   * \param [in] solverName The name of the desired third-party solver
   * \param [in] A \c C pointer to the coefficient matrix
   *
   * \return A \c C pointer to an Amesos2 solver interface.
   *
   * \throw std::invalid_argument The third-party solver named by \c
   * solverName is not supported.
   *
   * \relatesalso Amesos2::Solver
   */
  template < class Matrix,
             class Vector >
  Solver<Matrix,Vector>*
  create(const std::string& solverName, const Matrix* A);


  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A.
   *
   * Suitable for cases where numeric factorization must be performed
   * before the X and B vectors are known.  Before a solve, the \c
   * setX() and \c setB() functions should be used to set X and B, or
   * the overloaded solve(X,B) method should be used.
   *
   * \param [in] solverName The name of the desired third-party solver
   * \param [in] A <tt>Teuchos::RCP</tt> to the coefficient matrix
   *
   * \return A <tt>Teuchos::RCP</tt> to an Amesos2 solver interface.
   *
   * \throw std::invalid_argument The third-party solver named by \c
   * solverName is not supported.
   *
   * \relatesalso Amesos2::Solver
   */
  template < class Matrix,
             class Vector >
  Teuchos::RCP<Solver<Matrix,Vector> >
  create(const std::string& solverName,
         const Teuchos::RCP<const Matrix> A);


  /////////////////////////////////////////////////////
  // Meta-functions to help with creation of solvers //
  /////////////////////////////////////////////////////

  template < template <class,class> class ConcreteSolver,
             class Matrix,
             class Vector >
  struct create_solver_with_supported_type {
    static Teuchos::RCP<Solver<Matrix,Vector> > apply(Teuchos::RCP<const Matrix> A,
                                                      Teuchos::RCP<Vector>       X,
                                                      Teuchos::RCP<const Vector> B )
    {
      ctassert<
        std::is_same_v<
          typename MatrixTraits<Matrix>::scalar_t,
          typename MultiVecAdapter<Vector>::scalar_t
        >
      > same_scalar_assertion;
      (void)same_scalar_assertion; // This stops the compiler from warning about unused declared variables

      // If our assertion did not fail, then create and return a new solver
      return rcp( new ConcreteSolver<Matrix,Vector>(A, X, B) );
    }
  };

/**
 * \internal
 *
 * If the apply method of this struct is ever called, then it means
 * that the user requested to create a concrete solver interface for
 * a matrix whose scalar type is not supported by the solver.  In
 * such a case we throw a runtime exception.
 */
template < template <class,class> class ConcreteSolver,
           class Matrix,
           class Vector >
struct throw_no_scalar_support_exception {
  static Teuchos::RCP<Solver<Matrix,Vector> > apply(Teuchos::RCP<const Matrix> A,
                                                    Teuchos::RCP<Vector>       X,
                                                    Teuchos::RCP<const Vector> B )
  {
    // typedef ConcreteSolver<Matrix,Vector> concretesolver_matrix_vector;
    typedef typename MatrixTraits<Matrix>::scalar_t scalar_t;
    TEUCHOS_TEST_FOR_EXCEPTION( true,
                        std::invalid_argument,
                        "The requested Amesos2 "
                        // << concretesolver_matrix_vector::name <<
                        " solver interface does not support the " <<
                        Teuchos::ScalarTraits<scalar_t>::name() <<
                        " scalar type." );
  }
};

template < template <class,class> class ConcreteSolver,
           class Matrix,
           class Vector >
struct throw_no_matrix_support_exception {
  static Teuchos::RCP<Solver<Matrix,Vector> > apply(Teuchos::RCP<const Matrix> A,
                                                    Teuchos::RCP<Vector>       X,
                                                    Teuchos::RCP<const Vector> B )
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true,
                        std::invalid_argument,
                        "This solver does not support the kokkos adapter." );
  }
};

  /**
   * \internal
   *
   * Utility meta-function which binds to an exception-throwing
   * runtime function if the solver does not support the scalar type
   * of the matrix.  Otherwise, if the scalar type is supported, then
   * this returns an RCP to a new concrete Amesos2 solver of the given
   * type.
   */
  template < template <class,class> class ConcreteSolver,
             class Matrix,
             class Vector >
  struct handle_solver_scalar_type_support {
    static Teuchos::RCP<Solver<Matrix,Vector> > apply(Teuchos::RCP<const Matrix> A,
                                                      Teuchos::RCP<Vector>       X,
                                                      Teuchos::RCP<const Vector> B )
    {
      return std::conditional_t<
      solver_supports_scalar<ConcreteSolver, typename MatrixTraits<Matrix>::scalar_t>::value,
        create_solver_with_supported_type<ConcreteSolver,Matrix,Vector>,
        throw_no_scalar_support_exception<ConcreteSolver,Matrix,Vector> >::apply(A, X, B);
    }
  };

  /**
   * \internal
   *
   * Utility meta-function which binds to an exception-throwing
   * runtime function if the solver does not support the scalar type
   * of the matrix or the matrix adapter.  Otherwise, if the scalar type and
   * adapter are supported, then this returns an RCP to a new concrete Amesos2
   * solver of the given type.
   */
  template < template <class,class> class ConcreteSolver,
             class Matrix,
             class Vector >
  struct handle_solver_matrix_and_type_support {
    static Teuchos::RCP<Solver<Matrix,Vector> > apply(Teuchos::RCP<const Matrix> A,
                                                      Teuchos::RCP<Vector>       X,
                                                      Teuchos::RCP<const Vector> B )
    {
      return std::conditional_t<
        solver_supports_matrix<ConcreteSolver, Matrix>::value,
        handle_solver_scalar_type_support<ConcreteSolver,Matrix,Vector>,
        throw_no_matrix_support_exception<ConcreteSolver,Matrix,Vector> >::apply(A, X, B);
    }
  };

  /////////////////////
  // Query Functions //
  /////////////////////

  /**
   * \brief Queries the Factory for support of the named third-party library.
   *
   * \return \c true if the solver is supported.
   *
   * \relatesalso Amesos2::Solver
   */
  bool query(const char* solverName);


  /**
   * \brief Queries the Factory for support of the named third-party library.
   *
   * \return \c true if the solver is supported.
   *
   * \relatesalso Amesos2::Solver
   */
  bool query(const std::string& solverName);


  /////////////////
  // Definitions //
  /////////////////

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
  Teuchos::RCP<Solver<Matrix,Vector> >
  create(Teuchos::RCP<const Matrix> A,
         Teuchos::RCP<Vector>       X,
         Teuchos::RCP<const Vector> B)
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
    // Pass non-owning Teuchos::RCP objects to other factory method
    return( create(solver, rcp(A,false), rcp(X,false), rcp(B,false)).getRawPtr() );
  }


  template <class Matrix,
            class Vector >
  Teuchos::RCP<Solver<Matrix,Vector> >
  create(const char* solverName,
         const Teuchos::RCP<const Matrix> A,
         const Teuchos::RCP<Vector>       X,
         const Teuchos::RCP<const Vector> B)
  {
    std::string solver = solverName;
    return( create(solver, A, X, B) );
  }


  template <class Matrix,
            class Vector >
  Solver<Matrix,Vector>*
  create(const std::string& solverName, const Matrix* A){
    return( create(solverName, rcp(A,false),
                   Teuchos::RCP<Vector>(),
                   Teuchos::RCP<const Vector>()).getRawPtr() );
  }


  template <class Matrix,
            class Vector >
  Teuchos::RCP<Solver<Matrix,Vector> >
  create(const std::string& solverName, const Teuchos::RCP<const Matrix> A){
    return( create(solverName, A, Teuchos::RCP<Vector>(), Teuchos::RCP<const Vector>()) );
  }


  template <class Matrix,
            class Vector >
  Teuchos::RCP<Solver<Matrix,Vector> >
  create(const std::string& solverName, const Matrix* A, Vector* X, const Vector* B)
  {
    // Pass non-owning Teuchos::RCP objects to other factory method
    return( create(solverName, rcp(A,false), rcp(X,false), rcp(B,false)) );
  }


  template <class Matrix,
            class Vector >
  Teuchos::RCP<Solver<Matrix,Vector> >
  create(const std::string& solver_name,
         const Teuchos::RCP<const Matrix> A,
         const Teuchos::RCP<Vector>       X,
         const Teuchos::RCP<const Vector> B)
  {
    std::string solverName = tolower(solver_name); // for easy string checking

    // Check for our native solver first.  Treat KLU and KLU2 as equals.
    //
    // We use compiler guards in case a user does want to disable KLU2
#ifdef HAVE_AMESOS2_SHYLU_NODEBASKER
    if((solverName == "ShyLUBasker") || (solverName == "shylubasker") || (solverName == "amesos2_shylubasker"))
    {
      return handle_solver_matrix_and_type_support<ShyLUBasker, Matrix,Vector>::apply(A,X,B);
    }
#endif

#ifdef HAVE_AMESOS2_BASKER
    if((solverName == "Basker") || (solverName == "basker") || (solverName == "amesos2_basker"))
    {
      return handle_solver_matrix_and_type_support<Basker, Matrix,Vector>::apply(A,X,B);
    }
#endif


#ifdef HAVE_AMESOS2_KLU2
    if((solverName == "amesos2_klu2") || (solverName == "klu2") ||
        (solverName == "amesos2_klu")  || (solverName == "klu")){
      return handle_solver_matrix_and_type_support<KLU2,Matrix,Vector>::apply(A, X, B);
    }
#endif

#ifdef HAVE_AMESOS2_SUPERLUDIST
    if((solverName == "amesos2_superludist") ||
       (solverName == "superludist") ||
       (solverName == "amesos2_superlu_dist") ||
       (solverName == "superlu_dist")){
      return handle_solver_matrix_and_type_support<Superludist,Matrix,Vector>::apply(A, X, B);
    }
#endif

#ifdef HAVE_AMESOS2_SUPERLUMT
    if((solverName == "amesos2_superlumt") ||
       (solverName == "superlumt") ||
       (solverName == "amesos2_superlu_mt") ||
       (solverName == "superlu_mt")){
      return handle_solver_matrix_and_type_support<Superlumt,Matrix,Vector>::apply(A, X, B);
    }
#endif

#ifdef HAVE_AMESOS2_UMFPACK
    if((solverName == "amesos2_umfpack") ||
       (solverName == "umfpack")){
      return handle_solver_matrix_and_type_support<Umfpack,Matrix,Vector>::apply(A, X, B);
    }
#endif

#ifdef HAVE_AMESOS2_SHYLU_NODETACHO
    if((solverName == "amesos2_tacho") ||
       (solverName == "tacho")){
      return handle_solver_matrix_and_type_support<TachoSolver,Matrix,Vector>::apply(A, X, B);
    }

#endif

#ifdef HAVE_AMESOS2_SUPERLU
    if((solverName == "amesos2_superlu") ||
       (solverName == "superlu")){
      return handle_solver_matrix_and_type_support<Superlu,Matrix,Vector>::apply(A, X, B);
    }
#endif

#ifdef HAVE_AMESOS2_PARDISO_MKL
    if((solverName == "amesos2_pardiso_mkl") ||
       (solverName == "pardiso_mkl") ||
       (solverName == "amesos2_pardisomkl")  ||
       (solverName == "pardisomkl")){
      return handle_solver_matrix_and_type_support<PardisoMKL,Matrix,Vector>::apply(A, X, B);
    }
#endif
#ifdef HAVE_AMESOS2_CSS_MKL
    if((solverName == "amesos2_css_mkl") ||
       (solverName == "css_mkl") ||
       (solverName == "amesos2_cssmkl")  ||
       (solverName == "cssmkl")){
      return handle_solver_matrix_and_type_support<CssMKL,Matrix,Vector>::apply(A, X, B);
    }
#endif

#ifdef HAVE_AMESOS2_LAPACK
    if((solverName == "amesos2_lapack") ||
       (solverName == "lapack")){
      return handle_solver_matrix_and_type_support<Lapack,Matrix,Vector>::apply(A, X, B);
    }
#endif


#ifdef HAVE_AMESOS2_MUMPS
    if((solverName == "MUMPS") || (solverName == "mumps") ||
       (solverName == "amesos2_MUMPS") || (solverName == "amesos2_mumps"))
      {
        return handle_solver_matrix_and_type_support<MUMPS,Matrix,Vector>::apply(A,X,B);
      }
#endif
              
#ifdef HAVE_AMESOS2_STRUMPACK
    if((solverName == "STRUMPACK") || (solverName == "strumpack") ||
       (solverName == "amesos2_STRUMPACK") || (solverName == "amesos2_strumpack"))
      {
        return handle_solver_matrix_and_type_support<STRUMPACK,Matrix,Vector>::apply(A,X,B);
      }
#endif

#if defined (HAVE_AMESOS2_CHOLMOD) && defined (HAVE_AMESOS2_EXPERIMENTAL)
    if(solverName == "amesos2_cholmod" || solverName == "cholmod")
      return handle_solver_matrix_and_type_support<Cholmod,Matrix,Vector>::apply(A, X, B);
#endif

#if defined (HAVE_AMESOS2_CUSOLVER) && defined (HAVE_AMESOS2_CUSPARSE)
    if(solverName == "amesos2_cusolver" || solverName == "cusolver")
      return handle_solver_matrix_and_type_support<cuSOLVER,Matrix,Vector>::apply(A, X, B);
#endif

    /* If none of the above conditionals are satisfied, then the solver
     * requested is not yet supported.  We throw a runtime exception stating
     * this, and return null.
     */
    std::string err_msg = solver_name + " is not enabled or is not supported";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, err_msg);
    //return( Teuchos::null ); // unreachable
  }

} // end namespace Amesos2

#endif  // AMESOS2_FACTORY_HPP
