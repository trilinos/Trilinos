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
   \file   NewSolver_decl.hpp
   \author Eric Bavier <etbavier@sandia.gov>
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
