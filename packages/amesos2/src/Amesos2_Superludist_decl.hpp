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
  \file   Amesos2_Superludist_decl.hpp
  \author Eric Bavier <etbavie@sandia.gov>
  \date   Tue Jun 21 13:32:31 MDT 2011

  \brief  Amesos2 SuperLU_Dist declarations.
*/


#ifndef AMESOS2_SUPERLUDIST_DECL_HPP
#define AMESOS2_SUPERLUDIST_DECL_HPP

#include "Amesos2_SolverCore.hpp"
#include "Amesos2_Superludist_MatrixHelper.hpp"
#include "Amesos2_Superludist_FunctionMap.hpp"


namespace Amesos {


/** \brief Amesos2 interface to the distributed memory version of SuperLU.
 *
 * The distributed memory version of SuperLU, SuperLU_DIST, is supported by
 * this Amesos2 interface.  Currently support is for the SuperLU_DIST
 * 2.5 version.
 *
 * \section slu_dist_options Supported Options
 *
 * Currently, the following parameters/options are recognized:
 *
 * <ul>
 *   <li> \c "nprocs"(int) : Specifies the number of threads to be spawned.
 *     Default is 1.</li>
 *   <li> \c "Equil" : { \c "YES" | \c "NO" } or, equivalently, { \c true | \c false }.
 *     Specifies whether the solver to equilibrate the matrix before solving.</li>
 *   <li> \c "IterRefine" : { \c "NO" | \c "SINGLE" | \c "DOUBLE" | \c "EXTRA"
 *     }. Specifies whether to perform iterative refinement, and in
 *     what precision to compute the residual. (Not currently supported)</li>
 *   <li> \c "ColPerm" which takes one of the following:
 *     <ul>
 *     <li> \c "NATURAL" : natural ordering.</li>
 *     <li> \c "MMD_AT_PLUS_A" : minimum degree ordering on the structure of
 *       \f$ A^T + A\f$ .</li>
 *     <li> \c "MMD_ATA" : minimum degree ordering on the structure of
 *       \f$ A T A \f$ .</li>
 *     <li> \c "COLAMD" : approximate minimum degree column ordering.
 *       (default)</li>
 *     <li> \c "MY_PERMC" : use the ordering given in the "perm_c" parameter
 *       given by the user.  The value of the "perm_c" parameter should be a
 *       Teuchos::Array<int> with length equal to the number of global columns.
 *       </li>
 *     </ul>
 *   <li> \c "perm_c" : a Teuchos::Array<int> with length equal to the number
 *     of global columns.  Will assume <tt> ColPerm = MY_PERMC</tt>.</li>
 * </ul>
 *
 * \warning After creation, the size of the matrix should not change
 * (i.e. when using setA())
 */
template <class Matrix,
          class Vector>
class Superludist : public SolverCore<Amesos::Superludist, Matrix, Vector>
{
  friend class SolverCore<Amesos::Superludist,Matrix,Vector>; // Give our base access
                                                              // to our private
                                                              // implementation funcs
public:

  /// Name of this solver interface.
  static const char* name;	// declaration. Initialization outside.

  typedef Superludist<Matrix,Vector>                                   type;
  typedef SolverCore<Amesos::Superludist,Matrix,Vector>          super_type;

  typedef Matrix                                                matrix_type;
  typedef Vector                                                vector_type;
  
  // Since typedef's are not inheritted, go grab them
  typedef typename super_type::scalar_type                      scalar_type;
  typedef typename super_type::local_ordinal_type        local_ordinal_type;
  typedef typename super_type::global_ordinal_type      global_ordinal_type;
  typedef typename super_type::global_size_type            global_size_type;
  typedef typename super_type::node_type                          node_type;

  typedef TypeMap<Amesos::Superludist,scalar_type>                 type_map;

  typedef typename type_map::type                                  slu_type;
  typedef typename type_map::magnitude_type                  magnitude_type;

  typedef FunctionMap<Amesos::Superludist,slu_type>            function_map;
  typedef MatrixHelper<Amesos::Superludist>                   matrix_helper;

  
  /// \name Constructor/Destructor methods
  //@{

  /**
   * \brief Initialize from Teuchos::RCP.
   *
   * \warning Should not be called directly!  Use instead
   * Amesos::create() to initialize a SuperLU_DIST interface.
   */
  Superludist(
    Teuchos::RCP<Matrix> A,
    Teuchos::RCP<Vector> X,
    Teuchos::RCP<Vector> B);


  /// Destructor
  ~Superludist( );

  //@}

private:

  /**
   * \brief Performs pre-ordering on the matrix to increase efficiency.
   *
   * SuperLU_DIST supports several forms of column permutations.  Refer
   * to \ref slu_mt_options for the available \c ColPerm options.
   */
  int preOrdering_impl();


  /**
   * \brief Perform symbolic factorization of the matrix using SuperLU_DIST.
   *
   * Called second in the sequence before numericFactorization.
   *
   * \throw std::runtime_error SuperLU_DIST is not able to factor the matrix.
   */
  int symbolicFactorization_impl();


  /**
   * \brief SuperLU_DIST specific numeric factorization
   *
   * SuperLU_DIST factors the matrix in a shared memory environment
   * using \c nprocs threads, where nprocs defaults to \c 1 if it is
   * not changed through \c setParameters().
   *
   * \throw std::runtime_error SuperLU_DIST is not able to factor the matrix
   */
  int numericFactorization_impl();


  /**
   * \brief SuperLU_DIST specific solve.
   *
   * Uses the symbolic and numeric factorizations, along with the RHS
   * vector \c B to solve the sparse system of equations.  The
   * solution is placed in X.
   *
   * \throw std::runtime_error SuperLU_DIST is not able to solve the system.
   *
   * \callgraph
   */
  int solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> > X,
		 const Teuchos::Ptr<MultiVecAdapter<Vector> > B) const;


  /**
   * \brief Determines whether the shape of the matrix is OK for this solver.
   *
   * SuperLU_DIST supports square matrices.
   */
  bool matrixShapeOK_impl() const;


  /**
   * This method is hooked in by our Amesos::Solver parent class, which
   * handles the status and control methods, and this method handles
   * solver-specific parameters.
   *
   * See also: \ref slu_mt_options
   */
  void setParameters_impl(
    const Teuchos::RCP<Teuchos::ParameterList> & parameterList );


  /**
   * Hooked in by Amesos::Solver parent class.
   *
   * \return a const Teuchos::ParameterList of all valid parameters for this
   * solver.
   */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters_impl() const;


  /**
   * Calculates a SuperLU_DIST grid size of \c nprow by \c npcol
   * processes which will try to utilize all \c nprocs available
   * processes, but in case of failure, will return a square grid that
   * may not use all \c nprocs processes.
   *
   * If you're ever not pleased with how the algorithm's heuristics
   * treat prime numbers, don't give a prime for \c nprocs.
   *
   * \note the \c nprocs , \c nprow and \c npcol parameters may be set
   * together directly with setParameters()
   */
  void get_default_grid_size(int nprocs, SLUD::int_t& nprow, SLUD::int_t& npcol) const;


  // struct holds all data necessary to make a superlu factorization or solve call
  mutable struct SLUData {
    SLUD::SuperMatrix A,;
    SLUD::SuperMatrix AC; ///< The column-permuted matrix which will be factored
    typename type_map::LUstruct_t LU; ///< storage for L and U factors
    SLUD::Glu_freeable_t glu_freeable; ///< freeable storage used during symbolic fact

    /// Communicator for parallel column-ordering and symbolic fact.
    /// The number of processors in this communicator shall be the
    /// next power of 2 less than grid->nprow * grid->npcol.
    int                            domains;
    MPI_Comm                       symb_comm;
    SLUD::int_t                   *sizes, *fstVtxSep; // memory allocated by get_perm_c_parmetis
    SLUD::Pslu_freeable_t          pslu_freeable;

    SLUD::superlu_options_t          options;
    SLUD::mem_usage_t                mem_usage;
    SLUD::gridinfo_t                 grid;
    MPI_Comm                         grid_comm;
    typename type_map::LUstruct_t    lu; ///< stores the L and U factors
    SLUD::SuperLUStat_t              stat;
    typename type_map::SOLVEstruct_t solve_struct;

    Teuchos::Array<magnitude_type> berr; ///< backward error bounds
    Teuchos::Array<magnitude_type> ferr; ///< forward error bounds

    SLUD::ScalePermstruct_t        scale_perm; // R, C, perm_r, and perm_c found in here
    Teuchos::Array<magnitude_type> R, C;       // equilibration scalings
    Teuchos::Array<magnitude_type> R1, C1;     // row-permutation scalings
    Teuchos::Array<SLUD::int_t>    perm_r, perm_c;

    SLUD::DiagScale_t equed;	///< Whether/what kind of equilibration to use/has been used
    bool rowequ, colequ;	///< whether row/col equilibration has been applied to AC
    magnitude_type rowcnd, colcnd, amax;
  } data_;

  // The following Arrays are persisting storage arrays for A, X, and B
  /// Stores the values of the nonzero entries for SuperLU
  Teuchos::Array<slu_type> nzvals_;
  /// Stores the row indices of the nonzero entries
  Teuchos::Array<int> colind_;
  /// Stores the location in \c Ai_ and Aval_ that starts row j
  Teuchos::Array<int> rowptr_;
  /// 1D store for B and X values
  mutable Teuchos::Array<slu_type> bxvals_;
  mutable size_t ldbx_;

  /// \c true if this processor is in SuperLU_DISTS's 2D process grid
  bool in_grid_;
  bool same_symbolic_;
  mutable bool same_solve_struct_; // may be modified in solve_impl, but still `logically const'

  /// Maps rows of the matrix to processors in the SuperLU_DIST processor grid
  Teuchos::RCP<const Tpetra::Map<local_ordinal_type,
				 global_ordinal_type,
				 node_type> > superlu_rowmap_;
  
};				// End class Superludist


} // end namespace Amesos

#endif	// AMESOS2_SUPERLUDIST_DECL_HPP
