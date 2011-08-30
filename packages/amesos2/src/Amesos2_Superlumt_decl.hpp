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
  \file   Amesos2_Superlumt_decl.hpp
  \author Eric Bavier <etbavie@sandia.gov>
  \date   Mon May 23 17:06:53 MDT 2011

  \brief  Amesos2 SuperLU_MT declarations.
*/


#ifndef AMESOS2_SUPERLUMT_DECL_HPP
#define AMESOS2_SUPERLUMT_DECL_HPP

#include "Amesos2_SolverTraits.hpp"
#include "Amesos2_SolverCore_decl.hpp"
#include "Amesos2_Superlumt_FunctionMap.hpp"


namespace Amesos2 {


/** \brief Amesos2 interface to the Multi-threaded version of SuperLU.
 *
 * The multi-threaded version of SuperLU, SuperLU_MT, is supported by
 * this Amesos2 interface.  Currently support is for the SuperLU_MT
 * 2.0 version.
 *
 * See the \ref superlu_mt_parameters "summary of SuperLU_MT
 * parameters" supported this Amesos2 interface
 *
 * \ingroup amesos2_solver_interfaces
 */
template <class Matrix,
          class Vector>
class Superlumt : public SolverCore<Amesos2::Superlumt, Matrix, Vector>
{
  friend class SolverCore<Amesos2::Superlumt,Matrix,Vector>; // Give our base access
                                                            // to our private
                                                            // implementation funcs
public:

  /// Name of this solver interface.
  static const char* name;      // declaration. Initialization outside.

  typedef Superlumt<Matrix,Vector>                                     type;
  typedef SolverCore<Amesos2::Superlumt,Matrix,Vector>           super_type;

  // Since typedef's are not inheritted, go grab them
  typedef typename super_type::scalar_type                      scalar_type;
  typedef typename super_type::local_ordinal_type        local_ordinal_type;
  typedef typename super_type::global_ordinal_type      global_ordinal_type;
  typedef typename super_type::global_size_type            global_size_type;

  typedef TypeMap<Amesos2::Superlumt,scalar_type>                  type_map;

  typedef typename type_map::type                                  slu_type;
  typedef typename type_map::magnitude_type                  magnitude_type;

  typedef FunctionMap<Amesos2::Superlumt,slu_type>             function_map;


  /// \name Constructor/Destructor methods
  //@{

  /**
   * \brief Initialize from Teuchos::RCP.
   *
   * \warning Should not be called directly!  Use instead
   * Amesos2::create() to initialize a SuperLU_MT interface.
   */
  Superlumt(Teuchos::RCP<const Matrix> A,
            Teuchos::RCP<Vector>       X,
            Teuchos::RCP<const Vector> B);


  /// Destructor
  ~Superlumt( );

  //@}

private:

  /**
   * \brief Performs pre-ordering on the matrix to increase efficiency.
   *
   * SuperLU_MT supports several forms of column permutations.  Refer
   * to \ref slu_mt_options for the available \c ColPerm options.
   */
  int preOrdering_impl();


  /**
   * \brief Perform symbolic factorization of the matrix using SuperLU_MT.
   *
   * Called second in the sequence before numericFactorization.
   *
   * \throw std::runtime_error SuperLU_MT is not able to factor the matrix.
   */
  int symbolicFactorization_impl();


  /**
   * \brief SuperLU_MT specific numeric factorization
   *
   * SuperLU_MT factors the matrix in a shared memory environment
   * using \c nprocs threads, where nprocs defaults to \c 1 if it is
   * not changed through \c setParameters().
   *
   * \throw std::runtime_error SuperLU_MT is not able to factor the matrix
   */
  int numericFactorization_impl();


  /**
   * \brief SuperLU_MT specific solve.
   *
   * Uses the symbolic and numeric factorizations, along with the RHS
   * vector \c B to solve the sparse system of equations.  The
   * solution is placed in X.
   *
   * \throw std::runtime_error SuperLU_MT is not able to solve the system.
   *
   * \callgraph
   */
  int solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> > X,
                 const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const;


  /**
   * \brief Determines whether the shape of the matrix is OK for this solver.
   *
   * SuperLU_MT supports square matrices.
   */
  bool matrixShapeOK_impl() const;


  /**
   * The following SuperLU_MT parameters/options are recognized:
   *
   * <ul>
   *   <li> \c "nprocs" (int) : Specifies the number of threads to be spawned.
   *     Default is 1.</li>
   *   <li> \c "trans" : { \c "NOTRANS" | \c "TRANS" |
   *     \c "CONJ" }.  Will also recognize the \c "Transpose" : { \c true
   *     | \c false } option which is equivalent to \c "TRANS" and
   *     \c "NOTRANS" , respectively.</li>
   *   <li> \c "panel_size" (int) : Specifies the number of consecutive
   *     columns to be treated as a unit of task.</li>
   *   <li> \c "relax" (int) : Specifies the number of columns to be grouped as a relaxed
   *     supernode.</li>
   *   <li> \c "Equil" : { \c true | \c false }.  Specifies whether
   *     the solver to equilibrate the matrix before solving.</li>
   *   <li> \c "SymmetricMode" : { \c true | \c false }.  Specifies
   *   whether to use the symmetric mode.</li>
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
   *
   * Note that the \c nprocs, \c panel_size, and \c relax options are
   * recognized by SuperLU_MT but not by SuperLU.  Note also that it
   * is no typo in "trans", it really is lower-case (as opposed to
   * upper-case in SuperLU)
   */
  void setParameters_impl(
    const Teuchos::RCP<Teuchos::ParameterList> & parameterList );
  /* Parameters to support in the future:
   *
   *  <li> \c "IterRefine" : { \c "NO" | \c "SINGLE" | \c "DOUBLE" | \c "EXTRA"
   *     }. Specifies whether to perform iterative refinement, and in
   *     what precision to compute the residual. (Not currently supported)</li>
   */

  /**
   * Hooked in by Amesos2::SolverCore parent class.
   *
   * \return a const Teuchos::ParameterList of all valid parameters for this
   * solver.
   */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters_impl() const;


  /**
   * \brief Reads matrix data into internal structures
   *
   * \param [in] current_phase an indication of which solution phase this
   *                           load is being performed for.
   *
   * \return \c true if the matrix was loaded, \c false if not
   */
  bool loadA_impl(EPhase current_phase);


  // struct holds all data necessary to make a superlu factorization or solve call
  mutable struct SLUData {
    SLUMT::SuperMatrix A, BX, L, U;
    SLUMT::SuperMatrix AC; ///< The column-permuted matrix which will be factored

    SLUMT::superlumt_options_t options;
    SLUMT::superlu_memusage_t mem_usage;
    SLUMT::Gstat_t stat;

    Teuchos::Array<magnitude_type> berr; ///< backward error bounds
    Teuchos::Array<magnitude_type> ferr; ///< forward error bounds
    Teuchos::Array<int> perm_r;
    Teuchos::Array<int> perm_c;
    Teuchos::Array<magnitude_type> R;
    Teuchos::Array<magnitude_type> C;

    // in contrast to SuperLU, memory for etree will be allocated by
    // pxgssvx and the pointer will be stored in `options'

    SLUMT::equed_t equed;       ///< Whether/what kind of equilibration to use
    bool rowequ, colequ;        ///< whether row/col equilibration has been applied to AC
  } data_;

  // The following Arrays are persisting storage arrays for A, X, and B
  /// Stores the values of the nonzero entries for SuperLU
  Teuchos::Array<typename TypeMap<Amesos2::Superlumt,scalar_type>::type> nzvals_;
  /// Stores the row indices of the nonzero entries
  Teuchos::Array<int> rowind_;
  /// Stores the location in \c Ai_ and Aval_ that starts row j
  Teuchos::Array<int> colptr_;

  /* Note: In the above, must use "Amesos2::Superlumt" rather than
   * "Superlumt" because otherwise the compiler references the
   * specialized type of the class, and not the templated type that is
   * required for Amesos2::TypeMap
   */

  /* SuperLU can accept input in either compressed-row or
   * compressed-column storage.  We will store and pass matrices in
   * *compressed-row* format because that is the format Amesos used.
   */

};                              // End class Superlumt


// Specialize the solver_traits struct for SuperLU_MT
template <>
struct solver_traits<Superlumt> {
#ifdef HAVE_TEUCHOS_COMPLEX
  typedef Meta::make_list6<float,
                           double,
                           std::complex<float>,
                           std::complex<double>,
                           SLUMT::C::complex,
                           SLUMT::Z::doublecomplex> supported_scalars;
#else
  typedef Meta::make_list2<float, double> supported_scalars;
#endif
};

} // end namespace Amesos2

#endif  // AMESOS2_SUPERLUMT_DECL_HPP
