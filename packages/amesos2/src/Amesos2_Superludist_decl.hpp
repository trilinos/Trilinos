// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
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

#include "Amesos2_SolverTraits.hpp"
#include "Amesos2_SolverCore.hpp"
#include "Amesos2_Superludist_FunctionMap.hpp"

namespace Amesos2 {


/** \brief Amesos2 interface to the distributed memory version of SuperLU.
 *
 * The distributed memory version of SuperLU, SuperLU_DIST, is
 * supported by this Amesos2 interface.  Currently support is for the
 * SuperLU_DIST 2.5 version.
 *
 * This interface to SuperLU_DIST currently does not support row
 * permutations due to a sequential bottleneck present in SuperLU_DIST
 * when calculating such a row permutation (i.e. the matrix must be
 * brought to the root processor, which then finds the row permutation
 * and broadcasts this permutation to the other processors).  In the
 * future we may support row permutations through Zoltan.  By not
 * supporting this option, we make use of the entirely distributed
 * interface to SuperLU_DIST.  On the other hand, if you absolutely
 * need row permutations and your matrix will fit on a single node,
 * then you may consider using Amesos2's SuperLU_MT interface instead.
 *
 * SuperLU_DIST does not provide a means to directly get the number of
 * non-zeros in the L and U factors.
 *
 * \warning After creation, the size of the matrix should not change
 * (i.e. when using setA())
 *
 * \ingroup amesos2_solver_interfaces
 */
template <class Matrix,
          class Vector>
class Superludist : public SolverCore<Amesos2::Superludist, Matrix, Vector>
{
  friend class SolverCore<Amesos2::Superludist,Matrix,Vector>; // Give our base access
                                                              // to our private
                                                              // implementation funcs
public:

  /// Name of this solver interface.
  static const char* name;      // declaration. Initialization outside.

  typedef Superludist<Matrix,Vector>                                   type;
  typedef SolverCore<Amesos2::Superludist,Matrix,Vector>         super_type;

  typedef Matrix                                                matrix_type;
  typedef Vector                                                vector_type;

  // Since typedef's are not inheritted, go grab them
  typedef typename super_type::scalar_type                      scalar_type;
  typedef typename super_type::local_ordinal_type        local_ordinal_type;
  typedef typename super_type::global_ordinal_type      global_ordinal_type;
  typedef typename super_type::global_size_type            global_size_type;
  typedef typename super_type::node_type                          node_type;

  typedef TypeMap<Amesos2::Superludist,scalar_type>                type_map;

  typedef typename type_map::type                                  slu_type;
  typedef typename type_map::magnitude_type                  magnitude_type;

  typedef FunctionMap<Amesos2::Superludist,slu_type>           function_map;


  /// \name Constructor/Destructor methods
  //@{

  /**
   * \brief Initialize from Teuchos::RCP.
   *
   * \warning Should not be called directly!  Use instead
   * Amesos2::create() to initialize a SuperLU_DIST interface.
   */
  Superludist(Teuchos::RCP<const Matrix> A,
              Teuchos::RCP<Vector>       X,
              Teuchos::RCP<const Vector> B);


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
                 const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const;


  /**
   * \brief Determines whether the shape of the matrix is OK for this solver.
   *
   * SuperLU_DIST supports square matrices.
   */
  bool matrixShapeOK_impl() const;


  /**
   * Currently, the following SuperLU_DIST parameters/options are recognized:
   *
   * <ul>
   *   <li> \c "npcol"(int) and "nprow"(int) : Specified together, these parameters
   *     set the size of the SuperLU_DIST processor grid to \c nprow rows by
   *     \c npcol columns.  If these parameters are not set, the SuperLU_DIST
   *     interface uses a heuristic to pick the grid dimensions based on the
   *     number of processors in the matrix' communicator.</li>
   *   <li> \c "ColPerm" which takes one of the following:
   *     <ul>
   *     <li> \c "NATURAL" : natural column ordering.</li>
   *     <li> \c "PARMETIS" : use the ParMETIS TPL to order the columns. (default)</li>
   *     </ul>
   *   </li>
   *   <li> \c "ReplaceTinyPivot" : { \c true | \c false }.  Specifies
   *     whether to replace tiny diagonals with
   *     \f$\sqrt{\epsilon}\cdot\| A \|\f$ during LU
   *     factorization. (default: \c true)</li>
   * </ul>
   */
  /*
   * The following options could be supported in the future:
   *
   *   <li> \c "Equil" : { \c "YES" | \c "NO" } or, equivalently, { \c true | \c false }.
   *     Specifies whether the solver to equilibrate the matrix before solving.</li>
   *   <li> \c "IterRefine" : { \c "NO" | \c "SINGLE" | \c "DOUBLE" | \c "EXTRA"
   *     }. Specifies whether to perform iterative refinement, and in
   *     what precision to compute the residual. (Not currently supported)</li>
   */
  void setParameters_impl(
    const Teuchos::RCP<Teuchos::ParameterList> & parameterList );


  /**
   * Hooked in by Amesos2::Solver parent class.
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


  /**
   * \brief Reads matrix data into internal solver structures.
   *
   * Loads data from the matrix A into the internal SuperLU_DIST
   * matrix structure.  This function requires communication accross
   * all processors as the matrix is redistributed as necessary to the
   * processors in the SuperLU_DIST process grid.
   *
   * \post
   * - nzvals_, colind_, and rowptr_ arrays are sized to match the portion
   *   of the matrix on this processor.
   *
   * \return \c true if the matrix was loaded, \c false if not
   */
  bool loadA_impl(EPhase current_phase);


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
    MPI_Comm                         mat_comm; ///< Raw communicator used by the matrix A
    typename type_map::LUstruct_t    lu; ///< stores the L and U factors
    SLUD::SuperLUStat_t              stat;
    typename type_map::SOLVEstruct_t solve_struct;

    Teuchos::Array<magnitude_type> berr; ///< backward error bounds
    Teuchos::Array<magnitude_type> ferr; ///< forward error bounds

    SLUD::ScalePermstruct_t        scale_perm; // R, C, perm_r, and perm_c found in here
    Teuchos::Array<magnitude_type> R, C;       // equilibration scalings
    Teuchos::Array<magnitude_type> R1, C1;     // row-permutation scalings
    Teuchos::Array<SLUD::int_t>    perm_r, perm_c;

    SLUD::DiagScale_t equed;    ///< Whether/what kind of equilibration to use/has been used
    bool rowequ, colequ;        ///< whether row/col equilibration has been applied to AC
    magnitude_type rowcnd, colcnd, amax;
  } data_;

  // The following Arrays are persisting storage arrays for A, X, and B
  /// Stores the values of the nonzero entries for SuperLU_DIST
  Teuchos::Array<slu_type> nzvals_;
  /// Stores the row indices of the nonzero entries
  Teuchos::Array<int> colind_;
  /// Stores the location in \c Ai_ and Aval_ that starts row j
  Teuchos::Array<int> rowptr_;
  /// 1D store for B values
  mutable Teuchos::Array<slu_type> bvals_;
  /// 1D store for X values
  mutable Teuchos::Array<slu_type> xvals_;

  /// \c true if this processor is in SuperLU_DISTS's 2D process grid
  bool in_grid_;
  bool same_symbolic_;
  mutable bool same_solve_struct_; // may be modified in solve_impl, but still `logically const'

  /// Maps rows of the matrix to processors in the SuperLU_DIST processor grid
  Teuchos::RCP<const Tpetra::Map<local_ordinal_type,
                                 global_ordinal_type,
                                 node_type> > superlu_rowmap_;

};                              // End class Superludist


// Specialize the solver_traits template for SuperLU_DIST
template <>
struct solver_traits<Superludist> {
#ifdef HAVE_TEUCHOS_COMPLEX
  typedef Meta::make_list3<double, std::complex<double>, SLUD::Z::doublecomplex> supported_scalars;
#else
  typedef Meta::make_list1<double> supported_scalars;
#endif
};

} // end namespace Amesos2

#endif  // AMESOS2_SUPERLUDIST_DECL_HPP
