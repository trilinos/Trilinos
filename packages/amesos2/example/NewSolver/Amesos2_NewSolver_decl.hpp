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
   \file   NewSolver_decl.hpp
   \author John Doe <jd@sandia.gov>
   \date   Thu Jul  8 16:08:10 2010
   
   \brief  A template class that does nothing useful besides show developers
           what, in general, needs to be done to add a new solver interface to
           the Amesos2 collection.
*/


#ifndef AMESOS2_NEWSOLVER_DECL_HPP
#define AMESOS2_NEWSOLVER_DECL_HPP

#include "Amesos2_Solver.hpp"
#include "Amesos2_NewSolver_FunctionMap.hpp"
#include "Amesos2_NewSolver_TypeMap.hpp"


namespace Amesos {


/** \brief Amesos2 interface to the NewSolver package.
 *
 * \section newsolver_options Supported Options
 *
 * Currently, the following parameters/options are recognized:
 *
 * <ul>
 *   <li> \c "Equilibrate" : { \c true | \c false } </li>
 *   <li> \c "PivotThreshold" : \c double value </li>
 *   <li> \c "Transpose" : { \c "Transpose" | \c "NoTranspose" |
 *     \c "Conjugate" } </li>
 *   <li>etc...</li>
 * </ul>
 */
template <class Matrix,
          class Vector>
class NewSolver : public Solver<Amesos::NewSolver, Matrix, Vector> 
{
  friend class Solver<Amesos::NewSolver,Matrix,Vector>; // Give our base access
                                					              // to our private
				                                                // implementation funcs
public:
  static const char* name;	// declaration. Initialization outside.

  typedef NewSolver<Matrix,Vector>                                 solver_type;
  typedef Matrix                                                   matrix_type;
  typedef Vector                                                   vector_type;
  typedef typename MatrixAdapter<matrix_type>::scalar_type         scalar_type;
  typedef typename MatrixAdapter<matrix_type>::local_ordinal_type  local_ordinal_type;
  typedef typename MatrixAdapter<matrix_type>::global_ordinal_type global_ordinal_type;
  typedef typename MatrixAdapter<matrix_type>::global_size_type    global_size_type;
  typedef typename MatrixAdapter<matrix_type>::node_type           node_type;

  /// \name Constructor/Destructor methods
  //@{ 

  /**
   * \brief Initialize from Teuchos::RCP.
   *
   * \warning Should not be called directly!  Use instead
   * Amesos::Factory::create() to initialize a NewSolver interface.
   */
  NewSolver(
    Teuchos::RCP<Matrix> A,
    Teuchos::RCP<Vector> X,
    Teuchos::RCP<Vector> B);
  

  /// Destructor
  ~NewSolver( );

  //@}

private:
  
  /**
   * \brief Performs pre-ordering on the matrix to increase efficiency.
   */
  int preOrdering_impl();


  /**
   * \brief Perform symbolic factorization of the matrix using NewSolver.
   *
   * Called first in the sequence before numericFactorization.
   *
   * \throw std::runtime_error NewSolver is not able to factor the matrix.
   */
  int symbolicFactorization_impl();


  /**
   * \brief NewSolver specific numeric factorization
   * 
   * \throw std::runtime_error NewSolver is not able to factor the matrix
   */
  int numericFactorization_impl();


  /**
   * \brief NewSolver specific solve.
   *
   * Uses the symbolic and numeric factorizations, along with the RHS vector
   * \c B to solve the sparse system of equations.
   *
   * \throw std::runtime_error NewSolver is not able to solve the system.
   */
  int solve_impl();


  /**
   * \brief Determines whether the shape of the matrix is OK for this solver.
   */
  bool matrixShapeOK_impl() const;


  /**
   * This method is hooked in by our Amesos2::Solver parent class, which
   * handles the status and control methods, and this method handles
   * solver-specific parameters.
   *
   * See also: \ref newsolver_options 
   */
  void setParameters_impl(
    const Teuchos::RCP<Teuchos::ParameterList> & parameterList );


  /* Declare private variable necessary for interaction with the NewSolver
   * TPL.
   *
   * For example, the following Arrays are persisting storage arrays for A, X,
   * and B that can be used with solvers expecting a compressed-column
   * representation of the matrix A.
   */
   
  /// Stores the values of the nonzero entries for NewSolver
  Teuchos::Array<typename TypeMap<Amesos::NewSolver,scalar_type>::type> nzvals_;
  /// Stores the location in \c Ai_ and Aval_ that starts row j
  Teuchos::Array<int> rowind_;
  /// Stores the row indices of the nonzero entries
  Teuchos::Array<int> colptr_;
  /// Persisting 1D store for X
  Teuchos::Array<typename TypeMap<Amesos::NewSolver,scalar_type>::type> xvals_;  int ldx_;
  /// Persisting 1D store for B
  Teuchos::Array<typename TypeMap<Amesos::NewSolver,scalar_type>::type> bvals_;  int ldb_;

  /* Note: In the above, must use "Amesos::NewSolver" rather than "NewSolver"
   * because otherwise the compiler references the specialized type of the
   * class, and not the templated type that is required for Amesos::TypeMap.
   * This is only when referencing Amesos::NewSolver from within the
   * Amesos::NewSolver class.
   */
  
};                              // End class NewSolver


} // end namespace Amesos

#endif	// AMESOS2_NEWSOLVER_DECL_HPP
