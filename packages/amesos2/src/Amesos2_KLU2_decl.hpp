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
  \file   Amesos2_KLU2_decl.hpp
  \author Siva Rajamanickam <srajama@sandia.gov>

  \brief  Amesos2 KLU2 declarations.
*/


#ifndef AMESOS2_KLU2_DECL_HPP
#define AMESOS2_KLU2_DECL_HPP

#include "Amesos2_SolverTraits.hpp"
#include "Amesos2_SolverCore.hpp"
#include "Amesos2_KLU2_FunctionMap.hpp"


namespace Amesos2 {


/** \brief Amesos2 interface to the KLU2 package.
 *
 * See the \ref KLU2_parameters "summary of KLU2 parameters"
 * supported by this Amesos2 interface.
 *
 * \ingroup amesos2_solver_interfaces
 */
template <class Matrix,
          class Vector>
class KLU2 : public SolverCore<Amesos2::KLU2, Matrix, Vector>
{
  friend class SolverCore<Amesos2::KLU2,Matrix,Vector>; // Give our base access
                                                          // to our private
                                                          // implementation funcs
public:

  /// Name of this solver interface.
  static const char* name;      // declaration. Initialization outside.

  typedef KLU2<Matrix,Vector>                                       type;
  typedef SolverCore<Amesos2::KLU2,Matrix,Vector>             super_type;

  // Since typedef's are not inheritted, go grab them
  typedef typename super_type::scalar_type                      scalar_type;
  typedef typename super_type::local_ordinal_type        local_ordinal_type;
  typedef typename super_type::global_ordinal_type      global_ordinal_type;
  typedef typename super_type::global_size_type            global_size_type;

  typedef TypeMap<Amesos2::KLU2,scalar_type>                    type_map;

  /*
   * The KLU2 interface will need two other typedef's, which are:
   * - the KLU2 type that corresponds to scalar_type and
   * - the corresponding type to use for magnitude
   */
  typedef typename type_map::type                                 klu2_type;
  typedef typename type_map::dtype                               klu2_dtype;

  typedef FunctionMap<Amesos2::KLU2,klu2_type>                 function_map;

  typedef Matrix                                                matrix_type;
  typedef MatrixAdapter<matrix_type>                    matrix_adapter_type;

  /// \name Constructor/Destructor methods
  //@{

  /**
   * \brief Initialize from Teuchos::RCP.
   *
   * \warning Should not be called directly!  Use instead
   * Amesos2::create() to initialize a KLU2 interface.
   */
  KLU2(Teuchos::RCP<const Matrix> A,
          Teuchos::RCP<Vector>       X,
          Teuchos::RCP<const Vector> B);


  /// Destructor
  ~KLU2( );

  //@}

private:

 /**
  * \brief can we optimize size_type and ordinal_type for straight pass through,
  * also check that is_contiguous_ flag set to true
  */
  bool single_proc_optimization() const;

  /**
   * \brief Performs pre-ordering on the matrix to increase efficiency.
   *
   * KLU2 does not support pre-ordering, so this method does nothing.
   */
  int preOrdering_impl();


  /**
   * \brief Perform symbolic factorization of the matrix using KLU2.
   *
   * Called first in the sequence before numericFactorization.
   *
   * \throw std::runtime_error KLU2 is not able to factor the matrix.
   */
  int symbolicFactorization_impl();


  /**
   * \brief KLU2 specific numeric factorization
   *
   * \throw std::runtime_error KLU2 is not able to factor the matrix
   */
  int numericFactorization_impl();


  /**
   * \brief KLU2 specific solve.
   *
   * Uses the symbolic and numeric factorizations, along with the RHS
   * vector \c B to solve the sparse system of equations.  The
   * solution is placed in X.
   *
   * \throw std::runtime_error KLU2 is not able to solve the system.
   *
   * \callgraph
   */
  int solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
                 const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const;


  /**
   * \brief Determines whether the shape of the matrix is OK for this solver.
   */
  bool matrixShapeOK_impl() const;


  /**
   * Currently, the following KLU2 parameters/options are
   * recognized and acted upon:
   *
   * <ul>
   *   <li> \c "Trans" : { \c "NOTRANS" | \c "TRANS" |
   *     \c "CONJ" }.  Specifies whether to solve with the transpose system.</li>
   *   <li> \c "Equil" : { \c true | \c false }.  Specifies whether
   *     the solver to equilibrate the matrix before solving.</li>
   *   <li> \c "IterRefine" : { \c "NO" | \c "SINGLE" | \c "DOUBLE" | \c "EXTRA"
   *     }. Specifies whether to perform iterative refinement, and in
   *     what precision to compute the residual.</li>
   *   <li> \c "SymmetricMode" : { \c true | \c false }.</li>
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

  // struct holds all data necessary for KLU2 factorization or solve call
  mutable struct KLU2Data {
      ::KLU2::klu_symbolic<klu2_dtype, local_ordinal_type> *symbolic_;
      ::KLU2::klu_numeric<klu2_dtype, local_ordinal_type> *numeric_;
      ::KLU2::klu_common<klu2_dtype, local_ordinal_type> common_;
  } data_ ;

  typedef Kokkos::DefaultHostExecutionSpace HostSpaceType;
  typedef Kokkos::View<local_ordinal_type*, HostSpaceType> host_ordinal_type_array;

  typedef Kokkos::View<klu2_type*, HostSpaceType>     host_value_type_array;

  // The following Views are persisting storage arrays for A, X, and B
  /// Stores the values of the nonzero entries for CHOLMOD
  host_value_type_array host_nzvals_view_;

  /// Stores the location in \c Ai_ and Aval_ that starts row j
  host_ordinal_type_array host_rows_view_;
  /// Stores the row indices of the nonzero entries
  host_ordinal_type_array host_col_ptr_view_;

  typedef typename Kokkos::View<klu2_type**, Kokkos::LayoutLeft, HostSpaceType>
    host_solve_array_t;

  /// Persisting 1D store for X
  mutable host_solve_array_t xValues_;

  /// Persisting 1D store for B
  mutable host_solve_array_t bValues_;

  /// Transpose flag
  /// 0: Non-transpose, 1: Transpose, 2: Conjugate-transpose
  int transFlag_;

  bool is_contiguous_;
};                              // End class KLU2


// Specialize solver_traits struct for KLU2
template <>
struct solver_traits<KLU2> {
#ifdef HAVE_TEUCHOS_COMPLEX
  typedef Meta::make_list6<float,
                           double,
                           Kokkos::complex<float>,
                           Kokkos::complex<double>,
                           std::complex<float>,
                           std::complex<double> > supported_scalars;
#else
  typedef Meta::make_list2<float, double> supported_scalars;
#endif
};

template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
struct solver_supports_matrix<KLU2,
  KokkosSparse::CrsMatrix<Scalar, LocalOrdinal, ExecutionSpace>> {
  static const bool value = true;
};

} // end namespace Amesos2

#endif  // AMESOS2_KLU2_DECL_HPP
