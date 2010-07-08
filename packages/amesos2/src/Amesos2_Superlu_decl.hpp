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

/*
 * symbolicFactorization_impl() performs no action other than making sure that
 * Factorization is performed.  numericFactorization_impl() performs the
 * factorization but no solve (right hand side is set to 0 vectors) reuses
 * factors only if ReuseFactorization_ is set.
 *  
 *  If FactorizationOK_() && ReuseSymbolic_ 
 *     reFactor()
 *  else
 *     factor() 
 *
 * solve() 
 *
 * factor() does everything from scratch:
 *  - Redistributes the data if necessary
 *  - Deletes any data structures left over from the previous call to
 *    Factor()
 *  - Copies the data into the format that SuperLU wants it
 *  - Calls dgssvx to factor the matrix with factor set to true
 * reFactor()
 *  - Redistributes the data if necessary
 *    + Attempting to check to make sure that the non-zero structure is
 *      unchanged
 *  - Copies the data into the format that SuperLU already has it
 *    FIRST PASS - assert( false ) 
 */

#ifndef AMESOS2_SUPERLU_DECL_HPP
#define AMESOS2_SUPERLU_DECL_HPP

#include "Amesos2_Solver.hpp"
#include "Amesos2_Superlu_FunctionMap.hpp"
#include "Amesos2_Superlu_TypeMap.hpp"


namespace Amesos {


template <class Matrix,
          class Vector>
class Superlu : public Solver<Amesos::Superlu, Matrix, Vector> 
{
  friend class Solver<Amesos::Superlu,Matrix,Vector>; // Give our base access
                              					              // to our private
				                                              // implementation funcs
public:
  static const char* name;	// declaration. Initialization outside.

  typedef Superlu<Matrix,Vector>                                   solver_type;
  typedef Matrix                                                   matrix_type;
  typedef Vector                                                   vector_type;
  typedef typename MatrixAdapter<matrix_type>::scalar_type         scalar_type;
  typedef typename MatrixAdapter<matrix_type>::local_ordinal_type  local_ordinal_type;
  typedef typename MatrixAdapter<matrix_type>::global_ordinal_type global_ordinal_type;
  typedef typename MatrixAdapter<matrix_type>::global_size_type    global_size_type;
  typedef typename MatrixAdapter<matrix_type>::node_type           node_type;

  //@{ \name Constructor/Destructor methods
    
  Superlu(
    Teuchos::RCP<Matrix> A,
    Teuchos::RCP<Vector> X,
    Teuchos::RCP<Vector> B);
  

  // Destructor
  ~Superlu( );

private:
  
  /** \brief Performs pre-ordering on the matrix to increase efficiency.
   *
   * Superlu does not support pre-ordering, so this method does nothing.
   */
  int preOrdering_impl();


  /** \brief Perform symbolic factorization of the matrix using Superlu.
   *
   * Called first in the sequence before numericFactorization.
   *
   * \throw std::runtime_error Superlu is not able to factor the matrix.
   */
  int symbolicFactorization_impl();


  /** \brief Superlu specific numeric factorization
   * 
   * \throw std::runtime_error Superlu is not able to factor the matrix
   */
  int numericFactorization_impl();


  /** \brief Superlu specific solve.
   *
   * Uses the symbolic and numeric factorizations, along with the RHS vector
   * \c B to solve the sparse system of equations.
   *
   * \throw std::runtime_error Superlu is not able to solve the system.
   */
  int solve_impl();


  /** \brief Determines whether the shape of the matrix is OK for this solver.
   */
  bool matrixShapeOK_impl() const;


  /* This method is hooked in by our Amesos2::Solver parent class, which
   * handles the status and control methods, and this method handles
   * solver-specific parameters.
   */
  void setParameters_impl(
    const Teuchos::RCP<Teuchos::ParameterList> & parameterList );


  /* Convert the matrix to the form that is expected by SuperLU (Ap,
   * Ai, Aval).  We assume that the matrix has already been formed
   * on processor 0 as serialMatrixA_
   */
  void factor();

    
  /* Convert matrix to the form that SuperLU expects (Ap, Ai, Aval).
   * We assume that convertToSerial() has already been called and
   * that serialMatrixA_ is stored as it should be.  Only processor
   * 0 does anything useful here.
   */
  void refactor();

  // Not used at the moment
  bool recreateSLUMatrixAndCheckSameStructure();
  

  typedef typename TypeMap<Amesos::Superlu,scalar_type>::magnitude_type magnitude_type;
    
  // struct holds all data necessary to make a superlu factorization or solve call
  struct SLUData {
    SLU::SuperMatrix A, B, X, L, U;
#ifdef USE_DGSTRF
    SLU::SuperMatrix AC;
#endif
    SLU::superlu_options_t options;
    SLU::mem_usage_t mem_usage;
    SLU::SuperLUStat_t stat;

    Teuchos::Array<magnitude_type> berr; ///<  forward error bounds
    Teuchos::Array<magnitude_type> ferr; ///< backward error bounds
    Teuchos::Array<int> perm_r;
    Teuchos::Array<int> perm_c;
    Teuchos::Array<int> etree;
    Teuchos::Array<magnitude_type> R;
    Teuchos::Array<magnitude_type> C;

    char equed;                  // Flags whether to equilibrate the matrix
                                 // before solve
    int relax;
    int panel_size;
  } data_;

  // The following Arrays are persisting storage arrays for A, X, and B
  /// Stores the values of the nonzero entries for SuperLU
  Teuchos::Array<typename TypeMap<Amesos::Superlu,scalar_type>::type> nzvals_;
  /// Stores the location in \c Ai_ and Aval_ that starts row j
  Teuchos::Array<int> rowind_;
  /// Stores the row indices of the nonzero entries
  Teuchos::Array<int> colptr_;
  /// Persisting 1D store for X
  Teuchos::Array<typename TypeMap<Amesos::Superlu,scalar_type>::type> xvals_;  int ldx_;
  /// Persisting 1D store for B
  Teuchos::Array<typename TypeMap<Amesos::Superlu,scalar_type>::type> bvals_;  int ldb_;

  /* Note: In the above, must use "Amesos::Superlu" rather than "Superlu"
   * because otherwise the compiler references the specialized type of the
   * class, and not the templated type that is required for Amesos::TypeMap
   */
  
  /* SuperLU can accept input in either compressed-row or compressed-column
   * storage.  We will store and pass matrices in *compressed-row* format
   * because that is the format Amesos used.
   */
  

  bool factorizationDone_;
  
};				// End class Superlu


} // end namespace Amesos

#endif	// AMESOS2_SUPERLU_DECL_HPP
