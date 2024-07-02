// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

  using type                = Cholmod<Matrix,Vector>;
  using super_type          = SolverCore<Amesos2::Cholmod,Matrix,Vector>;

  // Since using's are not inheritted, go grab them
  using scalar_type         = typename super_type::scalar_type;
  using local_ordinal_type  = typename super_type::local_ordinal_type;
  using global_ordinal_type = typename super_type::global_ordinal_type;
  using global_size_type    = typename super_type::global_size_type;
  using node_type           = typename super_type::node_type;

  using type_map            = TypeMap<Amesos2::Cholmod,scalar_type>;

  /*
   * The CHOLMOD interface will need two other using's, which are:
   * - the CHOLMOD type that corresponds to scalar_type and
   * - the corresponding type to use for magnitude
   */
  using chol_type           = typename type_map::type;
  using magnitude_type      = typename type_map::magnitude_type;

  using function_map        = FunctionMap<Amesos2::Cholmod,chol_type>;

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
   *   <li> \c "nmethods" : { \c int value }.  Specifies the number of different ordering methods to try.</li>
   *   <li> \c "print" : \c int value. Specifies the verbosity of the print statements.</li>
   *   <li> \c "dbound" : \c int value. Specifies the smallest absolute value on the diagonal D for the LDL factorization.</li>
   *   <li> \c "PreferUpper" : \c int value. Specifies whether the matrix will be stored in upper triangular form.</li>
   *   <li> \c "useGPU" : \c int value. 1: Use GPU is 1, 0: Do not use GPU, -1: ENV CHOLMOD_USE_GPU set GPU usage..</li>
   *   <li> \c "Enable_KokkosKernels_TriangularSolves" : \c bool value. Whether to use triangular solves.</li>
   *   <li> \c "CholmodInt" : \c bool value. Whether to use cholmod int form.</li>Whether GIDs contiguous the
   *              smallest absolute value on the diagonal D for the LDL factorization.</li>
   *   <li> \c "SuperNodal" : \c bool value. Whether to use super nodal.</li>
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

  // use_cholmod_int_type controls whether we use CHOLMOD_INT or CHOLMOD_LONG.
  // To preserve a simple interface for the user where this can be picked
  // simply by setting a parameter, we prepare both types of arrays and just
  // one will actually be used.
  typedef int size_int_type;
  typedef int ordinal_int_type;

  typedef long size_long_type;
  typedef long ordinal_long_type;

  typedef Kokkos::View<size_long_type*, HostExecSpaceType>       host_size_long_type_array;
  typedef Kokkos::View<ordinal_long_type*, HostExecSpaceType> host_ordinal_long_type_array;

  typedef Kokkos::View<size_int_type*, HostExecSpaceType>       host_size_int_type_array;
  typedef Kokkos::View<ordinal_int_type*, HostExecSpaceType> host_ordinal_int_type_array;

  typedef Kokkos::View<chol_type*, HostExecSpaceType>      host_value_type_array;

  // The following Views are persisting storage arrays for A, X, and B
  /// Stores the values of the nonzero entries for CHOLMOD
  host_value_type_array host_nzvals_view_;
  /// Stores the location in \c Ai_ and Aval_ that starts row j
  host_size_int_type_array host_rows_int_view_;
  host_size_long_type_array host_rows_long_view_;
  /// Stores the row indices of the nonzero entries
  host_ordinal_int_type_array host_col_ptr_int_view_;
  host_ordinal_long_type_array host_col_ptr_long_view_;

  typedef typename Kokkos::View<chol_type**, Kokkos::LayoutLeft, HostExecSpaceType>
    host_solve_array_t;

  /// Persisting 1D store for X
  mutable host_solve_array_t host_xValues_;

  /// Persisting 1D store for B
  mutable host_solve_array_t host_bValues_;

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV) && defined(KOKKOSKERNELS_ENABLE_TPL_CHOLMOD)

  using DeviceExecSpaceType= Kokkos::DefaultExecutionSpace;

  #ifdef KOKKOS_ENABLE_CUDA
    // solver will be UVM off even though Tpetra is CudaUVMSpace
    using DeviceMemSpaceType = typename Kokkos::CudaSpace;
  #elif KOKKOS_ENABLE_HIP
    // same as above, make the solver UVM off
    using DeviceMemSpaceType = typename Kokkos::HIPSpace;
  #else
    using DeviceMemSpaceType = typename DeviceExecSpaceType::memory_space;
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
  typedef KokkosKernels::Experimental::KokkosKernelsHandle <size_int_type, ordinal_int_type, chol_type,
    DeviceExecSpaceType, DeviceMemSpaceType, DeviceMemSpaceType> kernel_handle_int_type;
  typedef KokkosKernels::Experimental::KokkosKernelsHandle <size_long_type, ordinal_long_type, chol_type,
    DeviceExecSpaceType, DeviceMemSpaceType, DeviceMemSpaceType> kernel_handle_long_type;
  mutable kernel_handle_int_type device_int_khL_;
  mutable kernel_handle_int_type device_int_khU_;
  mutable kernel_handle_long_type device_long_khL_;
  mutable kernel_handle_long_type device_long_khU_;
#endif

  bool firstsolve;
  
  // Used as a hack around cholmod doing ordering and symfact together
  bool skip_symfact;

  Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> > map;

  bool is_contiguous_;
  bool use_triangular_solves_;
  bool use_cholmod_int_type_; // controls if Cholmod is int or long

  void triangular_solve_symbolic();
  void triangular_solve_numeric();

public: // for GPU
  void triangular_solve() const; // Only for internal use - public to support kernels
};                              // End class Cholmod

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
