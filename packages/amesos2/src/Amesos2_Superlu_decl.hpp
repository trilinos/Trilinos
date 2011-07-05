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
  \file   Amesos2_Superlu_decl.hpp
  \author Eric Bavier <etbavier@sandia.gov>
  \date   Thu Jul  8 22:52:34 2010

  \brief  Amesos2 Superlu declarations.
*/


#ifndef AMESOS2_SUPERLU_DECL_HPP
#define AMESOS2_SUPERLU_DECL_HPP

#include "Amesos2_SolverCore.hpp"
#include "Amesos2_Superlu_MatrixHelper.hpp"
#include "Amesos2_Superlu_FunctionMap.hpp"
#include "Amesos2_Superlu_TypeMap.hpp"


namespace Amesos {


/** \brief Amesos2 interface to the SuperLU package.
 *
 * See the \ref superlu_parameters "summary of SuperLU parameters"
 * supported by this Amesos2 interface.
 *
 * \ingroup amesos2_solver_interfaces
 */
template <class Matrix,
          class Vector>
class Superlu : public SolverCore<Amesos::Superlu, Matrix, Vector>
{
  friend class SolverCore<Amesos::Superlu,Matrix,Vector>; // Give our base access
                                                          // to our private
                                                          // implementation funcs
public:

  /// Name of this solver interface.
  static const char* name;	// declaration. Initialization outside.

  typedef Superlu<Matrix,Vector>                                       type;
  typedef SolverCore<Amesos::Superlu,Matrix,Vector>              super_type;

  // Since typedef's are not inheritted, go grab them
  typedef typename super_type::scalar_type                      scalar_type;
  typedef typename super_type::local_ordinal_type        local_ordinal_type;
  typedef typename super_type::global_ordinal_type      global_ordinal_type;
  typedef typename super_type::global_size_type            global_size_type;

  /*
   * The SuperLU interface will need two other typedef's, which are:
   * - the superlu type that corresponds to scalar_type and
   * - the corresponding type to use for magnitude
   */
  typedef typename TypeMap<Amesos::Superlu,scalar_type>::type           slu_type;
  typedef typename TypeMap<Amesos::Superlu,scalar_type>::magnitude_type magnitude_type;

  typedef FunctionMap<Amesos::Superlu,scalar_type>             function_map;
  typedef MatrixHelper<Amesos::Superlu>                       matrix_helper;

  /// \name Constructor/Destructor methods
  //@{

  /**
   * \brief Initialize from Teuchos::RCP.
   *
   * \warning Should not be called directly!  Use instead
   * Amesos::create() to initialize a Superlu interface.
   */
  Superlu(
    Teuchos::RCP<Matrix> A,
    Teuchos::RCP<Vector> X,
    Teuchos::RCP<Vector> B);


  /// Destructor
  ~Superlu( );

  //@}

private:

  /**
   * \brief Performs pre-ordering on the matrix to increase efficiency.
   *
   * Superlu does not support pre-ordering, so this method does nothing.
   */
  int preOrdering_impl();


  /**
   * \brief Perform symbolic factorization of the matrix using Superlu.
   *
   * Called first in the sequence before numericFactorization.
   *
   * \throw std::runtime_error Superlu is not able to factor the matrix.
   */
  int symbolicFactorization_impl();


  /**
   * \brief Superlu specific numeric factorization
   *
   * \throw std::runtime_error Superlu is not able to factor the matrix
   */
  int numericFactorization_impl();


  /**
   * \brief Superlu specific solve.
   *
   * Uses the symbolic and numeric factorizations, along with the RHS
   * vector \c B to solve the sparse system of equations.  The
   * solution is placed in X.
   *
   * \throw std::runtime_error Superlu is not able to solve the system.
   *
   * \callgraph
   */
  int solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> > X,
		 const Teuchos::Ptr<MultiVecAdapter<Vector> > B) const;
  // TODO: B should technically be declared const across the board


  /**
   * \brief Determines whether the shape of the matrix is OK for this solver.
   */
  bool matrixShapeOK_impl() const;


  /**
   * Currently, the following SuperLU parameters/options are
   * recognized and acted upon:
   *
   * <ul>
   *   <li> \c "Trans" : { \c "NOTRANS" | \c "TRANS" |
   *     \c "CONJ" }.  Will also recognize the \c "Transpose" : { \c true
   *     | \c false } option which is equivalent to \c "TRANS" and
   *     \c "NOTRANS" , respectively.</li>
   *   <li> \c "Equil" : { \c "YES" | \c "NO" } or, equivalently, { \c true | \c false }.
   *     Specifies whether the solver to equilibrate the matrix before solving.</li>
   *   <li> \c "IterRefine" : { \c "NO" | \c "SINGLE" | \c "DOUBLE" | \c "EXTRA"
   *     }. Specifies whether to perform iterative refinement, and in
   *     what precision to compute the residual.</li>
   *   <li> \c "SymmetricMode" : { \c "Yes" | \c "NO" } or, equivalently,
   *     { \c true | \c false }.
   *   <li> \c "DiagPivotThresh" : \c double value. Specifies the threshold
   *     used for a diagonal to be considered an acceptable pivot.</li>
   *   <li> \c "ColPerm" which takes one of the following:
   *     <ul>
   *     <li> \c "NATURAL" : natural ordering.</li>
   *     <li> \c "MMD_AT_PLUS_A" : minimum degree ordering on the structure of
   *       \f$ A^T + A\f$ .</li>
   *     <li> \c "MMD_ATA" : minimum degree ordering on the structure of
   *       \f$ A T A \f$ .</li>
   *     <li> \c "COLAMD" : approximate minimum degree column ordering.
   *       (default)</li>
   *     </ul>
   * </ul>
   */
  void setParameters_impl(
    const Teuchos::RCP<Teuchos::ParameterList> & parameterList );


  /**
   * Hooked in by Amesos::SolverCore parent class.
   *
   * \return a const Teuchos::ParameterList of all valid parameters for this
   * solver.
   */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters_impl() const;


  // struct holds all data necessary to make a superlu factorization or solve call
  mutable struct SLUData {
    SLU::SuperMatrix A, B, X, L, U; // matrix A in NCformat
    SLU::SuperMatrix AC;	// permuted matrix A in NCPformat

    SLU::superlu_options_t options;
    SLU::mem_usage_t mem_usage;
    SLU::SuperLUStat_t stat;

    Teuchos::Array<magnitude_type> berr; ///< backward error bounds
    Teuchos::Array<magnitude_type> ferr; ///<  forward error bounds
    Teuchos::Array<int> perm_r;
    Teuchos::Array<int> perm_c;
    Teuchos::Array<int> etree;
    Teuchos::Array<magnitude_type> R;
    Teuchos::Array<magnitude_type> C;

    char equed;
    bool rowequ, colequ;	// flags what type of equilibration
				// has been performed

    int relax;
    int panel_size;
  } data_;

  // The following Arrays are persisting storage arrays for A, X, and B
  /// Stores the values of the nonzero entries for SuperLU
  Teuchos::Array<slu_type> nzvals_;
  /// Stores the location in \c Ai_ and Aval_ that starts row j
  Teuchos::Array<int> rowind_;
  /// Stores the row indices of the nonzero entries
  Teuchos::Array<int> colptr_;

  /// Persisting 1D store for X
  Teuchos::Array<typename TypeMap<Amesos::Superlu,scalar_type>::type> xvals_;  int ldx_;
  /// Persisting 1D store for B
  Teuchos::Array<typename TypeMap<Amesos::Superlu,scalar_type>::type> bvals_;  int ldb_;

  /* Note: In the above, must use "Amesos::Superlu" rather than
   * "Superlu" because otherwise the compiler references the
   * specialized type of the class, and not the templated type that is
   * required for Amesos::TypeMap
   */

  /* SuperLU can accept input in either compressed-row or
   * compressed-column storage.  We will store and pass matrices in
   * *compressed-column* format.
   */

  /*
   * Internal flag that is used for the numericFactorization_impl
   * routine.  If true, then the superlu gstrf routine should have
   * SamePattern_SameRowPerm in its options.  Otherwise, it should
   * factor from scratch.
   *
   * This is somewhat of a kludge to get around the fact that the
   * superlu routines all expect something different from the options
   * struct.  The big issue is that we don't want gstrf doing the
   * symbolic factorization if it doesn't need to.  On the other hand,
   * we can't leave options.Fact set to SamePattern_SameRowPerm
   * because the solver driver needs it to be set at FACTORED.  But
   * having it set at FACTORED upon re-entrance into
   * numericFactorization prompts gstrf to redo the symbolic
   * factorization.
   */
  bool same_symbolic_;

};				// End class Superlu


} // end namespace Amesos

#endif	// AMESOS2_SUPERLU_DECL_HPP
