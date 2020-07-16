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
  \file   Amesos2_Cholmod_decl.hpp
  \author Kevin Deweese <kdewees@software.sandia.gov> 
  \date   Tue Aug 27 17:06:53 2013

  \brief  Amesos2 CHOLMOD declarations.
*/


#ifndef AMESOS2_CHOLMOD_DECL_HPP
#define AMESOS2_CHOLMOD_DECL_HPP

#include "Amesos2_SolverTraits.hpp"
#include "Amesos2_SolverCore.hpp"
#include "Amesos2_Cholmod_FunctionMap.hpp"

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV) && defined(KOKKOSKERNELS_ENABLE_TPL_CHOLMOD)
#include "KokkosKernels_Handle.hpp"
#endif

namespace Amesos2 {


/** \brief Amesos2 interface to the CHOLMOD package.
 *
 * See the \ref cholmod_parameters "summary of CHOLMOD parameters"
 * supported by this Amesos2 interface.
 *
 * \ingroup amesos2_solver_interfaces
 */
template <class Matrix,
          class Vector>
class Cholmod : public SolverCore<Amesos2::Cholmod, Matrix, Vector>
{
  friend class SolverCore<Amesos2::Cholmod,Matrix,Vector>; // Give our base access
                                                          // to our private
                                                          // implementation funcs
public:

  /// Name of this solver interface.
  static const char* name;      // declaration. Initialization outside.

  typedef Cholmod<Matrix,Vector>                                       type;
  typedef SolverCore<Amesos2::Cholmod,Matrix,Vector>             super_type;

  // Since typedef's are not inheritted, go grab them
  typedef typename super_type::scalar_type                    scalar_type;
  typedef typename super_type::local_ordinal_type      local_ordinal_type;
  typedef typename super_type::global_ordinal_type    global_ordinal_type;
  typedef typename super_type::global_size_type          global_size_type;
  typedef typename super_type::node_type                        node_type;

  typedef TypeMap<Amesos2::Cholmod,scalar_type>                    type_map;

  /*
   * The CHOLMOD interface will need two other typedef's, which are:
   * - the CHOLMOD type that corresponds to scalar_type and
   * - the corresponding type to use for magnitude
   */
  typedef typename type_map::type                                 chol_type;
  typedef typename type_map::magnitude_type                  magnitude_type;

  typedef FunctionMap<Amesos2::Cholmod,chol_type>              function_map;

  /// \name Constructor/Destructor methods
  //@{

  /**
   * \brief Initialize from Teuchos::RCP.
   *
   * \warning Should not be called directly!  Use instead
   * Amesos2::create() to initialize a CHOLMOD interface.
   */
  Cholmod(Teuchos::RCP<const Matrix> A,
          Teuchos::RCP<Vector>       X,
          Teuchos::RCP<const Vector> B);


  /// Destructor
  ~Cholmod( );

  //@}

private:

  /**
   * \brief Performs pre-ordering on the matrix to increase efficiency.
   */
  int preOrdering_impl();


  /**
   * \brief Perform symbolic factorization of the matrix using CHOLMOD.
   *
   * Called first in the sequence before numericFactorization.
   *
   * \throw std::runtime_error CHOLMOD is not able to factor the matrix.
   */
  int symbolicFactorization_impl();


  /**
   * \brief CHOLMOD specific numeric factorization
   *
   * \throw std::runtime_error CHOLMOD is not able to factor the matrix
   */
  int numericFactorization_impl();


  /**
   * \brief CHOLMOD specific solve.
   *
   * Uses the symbolic and numeric factorizations, along with the RHS
   * vector \c B to solve the sparse system of equations.  The
   * solution is placed in X.
   *
   * \throw std::runtime_error CHOLMOD is not able to solve the system.
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
   * Currently, the following CHOLMOD parameters/options are
   * recognized and acted upon:
   *
   * <ul>
   *   <li> \c "Trans" : { \c "NOTRANS" | \c "TRANS" |
   *     \c "CONJ" }.  Specifies whether to solve with the transpose system.</li>
   *   <li> \c "Equil" : { \c true | \c false }.  Specifies whether
   *     the solver to equilibrate the matrix before solving.</li>
   *   <li> \c "IterRefine" : { \c "NO" | \c "SLU_SINGLE" | \c "SLU_DOUBLE" | \c "EXTRA"
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


  // struct holds all data necessary to make a cholmod factorization or solve call
  mutable struct CholData {
    cholmod_sparse A;
    cholmod_dense x, b;
    cholmod_dense *Y, *E;
    cholmod_factor *L;
    cholmod_common c;
  } data_;

  typedef Kokkos::DefaultHostExecutionSpace                   HostExecSpaceType;
  typedef typename HostExecSpaceType::memory_space             HostMemSpaceType;

  typedef long size_type;
  typedef long ordinal_type;
  typedef Kokkos::View<size_type*, HostExecSpaceType>       host_size_type_array;
  typedef Kokkos::View<ordinal_type*, HostExecSpaceType> host_ordinal_type_array;

  typedef Kokkos::View<chol_type*, HostExecSpaceType>      host_value_type_array;

  // The following Views are persisting storage arrays for A, X, and B
  /// Stores the values of the nonzero entries for CHOLMOD
  host_value_type_array host_nzvals_view_;
  /// Stores the location in \c Ai_ and Aval_ that starts row j
  host_size_type_array host_rows_view_;
  /// Stores the row indices of the nonzero entries
  host_ordinal_type_array host_col_ptr_view_;

  typedef typename Kokkos::View<chol_type**, Kokkos::LayoutLeft, HostExecSpaceType>
    host_solve_array_t;

  /// Persisting 1D store for X
  mutable host_solve_array_t host_xValues_;

  /// Persisting 1D store for B
  mutable host_solve_array_t host_bValues_;

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV) && defined(KOKKOSKERNELS_ENABLE_TPL_CHOLMOD)

  typedef Kokkos::DefaultExecutionSpace DeviceExecSpaceType;

  #ifdef KOKKOS_ENABLE_CUDA
    // solver will be UVM off even though Tpetra is CudaUVMSpace
    typedef typename Kokkos::CudaSpace DeviceMemSpaceType;
  #else
    typedef typename DeviceExecSpaceType::memory_space DeviceMemSpaceType;
  #endif

  typedef Kokkos::View<chol_type**, Kokkos::LayoutLeft, DeviceMemSpaceType>
    device_solve_array_t;
  // For triangular solves we have both host and device versions of xValues and
  // bValues because a parameter can turn it on or off.
  mutable device_solve_array_t device_xValues_;
  mutable device_solve_array_t device_bValues_;
  typedef Kokkos::View<int*, HostMemSpaceType>                 host_int_array;
  typedef Kokkos::View<int*, DeviceMemSpaceType>              device_int_array;
  host_int_array host_trsv_etree_;
  host_int_array host_trsv_perm_;
  device_int_array device_trsv_perm_;
  mutable device_solve_array_t device_trsv_rhs_;
  mutable device_solve_array_t device_trsv_sol_;
  typedef KokkosKernels::Experimental::KokkosKernelsHandle <size_type, ordinal_type, chol_type,
    DeviceExecSpaceType, DeviceMemSpaceType, DeviceMemSpaceType> kernel_handle_type;
  mutable kernel_handle_type device_khL_;
  mutable kernel_handle_type device_khU_;
#endif

  bool firstsolve;
  
  // Used as a hack around cholmod doing ordering and symfact together
  bool skip_symfact;

  Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> > map;

  bool is_contiguous_;
  bool use_triangular_solves_;

  void triangular_solve_symbolic();
  void triangular_solve_numeric();

public: // for GPU
  void triangular_solve() const; // Only for internal use - public to support kernels
};                              // End class Cholmod


/* Specialize solver_traits struct for Cholmod
 *
 * Based on the CHOLMOD documentation, the support for
 * single-precision complex numbers is unclear.  Much of the
 * discussion of complex types only makes explicit mention of 'double'
 * types.  So, be pessimistic for now and don't declare
 * single-precision complex support
 */
template <>
struct solver_traits<Cholmod> {

// Cholmod does not yet support float.
#ifdef HAVE_TEUCHOS_COMPLEX
  typedef Meta::make_list3<double, std::complex<double>,
                           Kokkos::complex<double>> supported_scalars;
#else
  typedef Meta::make_list1<double> supported_scalars;
#endif
};

template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
struct solver_supports_matrix<Cholmod,
  KokkosSparse::CrsMatrix<Scalar, LocalOrdinal, ExecutionSpace>> {
  static const bool value = true;
};

} // end namespace Amesos2

#endif  // AMESOS2_CHOLMOD_DECL_HPP
