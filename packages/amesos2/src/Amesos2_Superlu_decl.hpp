// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
  \file   Amesos2_Superlu_decl.hpp
  \author Eric Bavier <etbavier@sandia.gov>
  \date   Thu Jul  8 22:52:34 2010

  \brief  Amesos2 Superlu declarations.
*/


#ifndef AMESOS2_SUPERLU_DECL_HPP
#define AMESOS2_SUPERLU_DECL_HPP

#include "Amesos2_SolverTraits.hpp"
#include "Amesos2_SolverCore.hpp"
#include "Amesos2_Superlu_FunctionMap.hpp"

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV) && defined(KOKKOSKERNELS_ENABLE_TPL_SUPERLU)
#include "KokkosKernels_Handle.hpp"
#endif

namespace Amesos2 {


/** \brief Amesos2 interface to the SuperLU package.
 *
 * See the \ref superlu_parameters "summary of SuperLU parameters"
 * supported by this Amesos2 interface.
 *
 * \ingroup amesos2_solver_interfaces
 */
template <class Matrix,
          class Vector>
class Superlu : public SolverCore<Amesos2::Superlu, Matrix, Vector>
{
  friend class SolverCore<Amesos2::Superlu,Matrix,Vector>; // Give our base access
                                                          // to our private
                                                          // implementation funcs
public:

  /// Name of this solver interface.
  static const char* name;      // declaration. Initialization outside.

  typedef Superlu<Matrix,Vector>                                       type;
  typedef SolverCore<Amesos2::Superlu,Matrix,Vector>             super_type;

  // Since typedef's are not inheritted, go grab them
  typedef typename super_type::scalar_type                      scalar_type;
  typedef typename super_type::local_ordinal_type        local_ordinal_type;
  typedef typename super_type::global_ordinal_type      global_ordinal_type;
  typedef typename super_type::global_size_type            global_size_type;

  typedef TypeMap<Amesos2::Superlu,scalar_type>                    type_map;

  /*
   * The SuperLU interface will need two other typedef's, which are:
   * - the superlu type that corresponds to scalar_type and
   * - the corresponding type to use for magnitude
   */
  typedef typename type_map::type                                  slu_type;
  typedef typename type_map::convert_type                  slu_convert_type;
  typedef typename type_map::magnitude_type                  magnitude_type;

  typedef FunctionMap<Amesos2::Superlu,slu_type>               function_map;

  /// \name Constructor/Destructor methods
  //@{

  /**
   * \brief Initialize from Teuchos::RCP.
   *
   * \warning Should not be called directly!  Use instead
   * Amesos2::create() to initialize a Superlu interface.
   */
  Superlu(Teuchos::RCP<const Matrix> A,
          Teuchos::RCP<Vector>       X,
          Teuchos::RCP<const Vector> B);


  /// Destructor
  ~Superlu( );

  //@}

  /// Returns a short description of this Solver
  std::string description() const override;

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
  int solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
                 const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const;


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
   *   <li> \c "ILU_Flag" : { \c false | \c true }. If true, run ILU routines.</li>
   *   <li> \c "ILU_DropTol" : \c double value. ILUT drop tolerance.</li>
   *   <li> \c "ILU_FillFactor" : \c double value. ILUT fill factor.</li>
   *   <li> \c "ILU_Norm" : { \c "INF_NORM" | \c "ONE_NORM" | \c "TWO_NORM"}.</li>
   *   <li> \c "ILU_MILU" : { \c "SILU" | \c "SMILU_1" | \c "SMILU_2" | \c "SMILU_3"}.  Type of modified ILU algorithm to use.</li>
   *   <li> \c "ILU_FillTol" : \c double value. ILUT fill tolerance.</li>
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

  typedef Kokkos::DefaultHostExecutionSpace HostExecSpaceType;

  // struct holds all data necessary to make a superlu factorization or solve call
  mutable struct SLUData {
    SLU::SuperMatrix A, B, X, L, U; // matrix A in NCformat
    SLU::SuperMatrix AC;        // permuted matrix A in NCPformat

    SLU::superlu_options_t options;
    SLU::mem_usage_t mem_usage;
#ifdef HAVE_AMESOS2_SUPERLU5_API
    SLU::GlobalLU_t lu;      // Use for gssvx and gsisx in SuperLU 5.0
#endif
    SLU::SuperLUStat_t stat;



    typedef Kokkos::View<magnitude_type*, HostExecSpaceType>    host_mag_array;
    typedef Kokkos::View<int*, HostExecSpaceType>               host_int_array;
    host_mag_array berr; ///< backward error bounds
    host_mag_array ferr; ///<  forward error bounds
    host_int_array perm_r;
    host_int_array perm_c;
    host_int_array etree;
    host_mag_array R;
    host_mag_array C;

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV) && defined(KOKKOSKERNELS_ENABLE_TPL_SUPERLU)
    host_int_array parents;
#endif

    char equed;
    bool rowequ, colequ;        // flags what type of equilibration
                                // has been performed
    magnitude_type anorm, rcond; // condition number estimate

    int relax;
    int panel_size;
  } data_;

  typedef int size_type;
  typedef int ordinal_type;
  typedef Kokkos::View<size_type*, HostExecSpaceType>       host_size_type_array;
  typedef Kokkos::View<ordinal_type*, HostExecSpaceType> host_ordinal_type_array;
  typedef Kokkos::View<slu_type*, HostExecSpaceType>       host_value_type_array;

  // The following Arrays are persisting storage arrays for A, X, and B
  /// Stores the values of the nonzero entries for SuperLU
  host_value_type_array host_nzvals_view_;
  Teuchos::Array<slu_convert_type> convert_nzvals_; // copy to SuperLU native array before calling SuperLU

  /// Stores the location in \c Ai_ and Aval_ that starts row j
  host_size_type_array host_rows_view_;
  /// Stores the row indices of the nonzero entries
  host_ordinal_type_array host_col_ptr_view_;

  typedef typename Kokkos::View<slu_type**, Kokkos::LayoutLeft, HostExecSpaceType>
    host_solve_array_t;

  /// Persisting 1D store for X
  mutable host_solve_array_t host_xValues_;
  mutable Teuchos::Array<slu_convert_type> convert_xValues_; // copy to SuperLU native array before calling SuperLU

  /// Persisting 1D store for B
  mutable host_solve_array_t host_bValues_;
  mutable Teuchos::Array<slu_convert_type> convert_bValues_; // copy to SuperLU native array before calling SuperLU

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV) && defined(KOKKOSKERNELS_ENABLE_TPL_SUPERLU)
  typedef Kokkos::DefaultExecutionSpace DeviceExecSpaceType;

  #ifdef KOKKOS_ENABLE_CUDA
    // solver will be UVM off even though Tpetra is CudaUVMSpace
    typedef typename Kokkos::CudaSpace DeviceMemSpaceType;
  #else
    typedef typename DeviceExecSpaceType::memory_space DeviceMemSpaceType;
  #endif

  typedef Kokkos::View<slu_type**, Kokkos::LayoutLeft, DeviceMemSpaceType>
    device_solve_array_t;
  // For triangular solves we have both host and device versions of xValues and
  // bValues because a parameter can turn it on or off.
  mutable device_solve_array_t device_xValues_;
  mutable device_solve_array_t device_bValues_;
  typedef Kokkos::View<int*,            DeviceMemSpaceType>      device_int_array;
  typedef Kokkos::View<magnitude_type*, DeviceMemSpaceType>      device_mag_array;
  device_int_array device_trsv_perm_r_;
  device_int_array device_trsv_perm_c_;
  device_mag_array device_trsv_R_;
  device_mag_array device_trsv_C_;
  mutable device_solve_array_t device_trsv_rhs_;
  mutable device_solve_array_t device_trsv_sol_;
  typedef KokkosKernels::Experimental::KokkosKernelsHandle <size_type, ordinal_type, slu_type,
    DeviceExecSpaceType, DeviceMemSpaceType, DeviceMemSpaceType> kernel_handle_type;
  mutable kernel_handle_type device_khL_;
  mutable kernel_handle_type device_khU_;
  /* parameters for SpTRSV */
  bool sptrsv_invert_diag_;
  bool sptrsv_invert_offdiag_;
  bool sptrsv_u_in_csr_;
  bool sptrsv_merge_supernodes_;
  bool sptrsv_use_spmv_;
#endif

  /* Note: In the above, must use "Amesos2::Superlu" rather than
   * "Superlu" because otherwise the compiler references the
   * specialized type of the class, and not the templated type that is
   * required for Amesos2::TypeMap
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
  bool ILU_Flag_;

  bool is_contiguous_;
  bool use_triangular_solves_;

  void triangular_solve_factor();

  /* call metis before SuperLU */
  bool use_metis_;
  bool symmetrize_metis_;

  public: // for GPU
    void triangular_solve() const; // Only for internal use - public to support kernels
};                              // End class Superlu


// Specialize solver_traits struct for SuperLU
template <>
struct solver_traits<Superlu> {
#ifdef HAVE_TEUCHOS_COMPLEX
  typedef Meta::make_list6<float, double,
                           std::complex<float>, std::complex<double>,
                           Kokkos::complex<float>, Kokkos::complex<double>>
                           supported_scalars;
#else
  typedef Meta::make_list2<float, double> supported_scalars;
#endif
};

template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
struct solver_supports_matrix<Superlu,
  KokkosSparse::CrsMatrix<Scalar, LocalOrdinal, ExecutionSpace>> {
  static const bool value = true;
};

} // end namespace Amesos2

#endif  // AMESOS2_SUPERLU_DECL_HPP
