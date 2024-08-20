// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef AMESOS2_CUSOLVER_DECL_HPP
#define AMESOS2_CUSOLVER_DECL_HPP

#include "Amesos2_SolverTraits.hpp"
#include "Amesos2_SolverCore.hpp"
#include "Amesos2_cuSOLVER_FunctionMap.hpp"

namespace Amesos2 {

/** \brief Amesos2 interface to cuSOLVER.
 *
 * \ingroup amesos2_solver_interfaces
 */
template <class Matrix,
          class Vector>
class cuSOLVER : public SolverCore<Amesos2::cuSOLVER, Matrix, Vector>
{
  friend class SolverCore<Amesos2::cuSOLVER,Matrix,Vector>;

public:

  /// Name of this solver interface.
  static const char* name;      // declaration. Initialization outside.

  typedef cuSOLVER<Matrix,Vector>                                        type;
  typedef SolverCore<Amesos2::cuSOLVER,Matrix,Vector>              super_type;

  typedef typename super_type::scalar_type                        scalar_type;
  typedef typename super_type::local_ordinal_type          local_ordinal_type;
  typedef typename super_type::global_ordinal_type        global_ordinal_type;
  typedef typename super_type::global_size_type              global_size_type;
  typedef typename super_type::node_type                            node_type;

  typedef TypeMap<Amesos2::cuSOLVER,scalar_type>                     type_map;

  typedef typename type_map::type                               cusolver_type;
  typedef typename type_map::magnitude_type                    magnitude_type;

  typedef FunctionMap<Amesos2::cuSOLVER,cusolver_type>           function_map;

  #ifdef KOKKOS_ENABLE_CUDA
    // solver will be UVM off
    typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>       device_type;
  #else
    typedef Kokkos::DefaultExecutionSpace::device_type            device_type;
  #endif

  typedef int size_type;
  typedef int ordinal_type;
  typedef Kokkos::View<size_type*, device_type>        device_size_type_array;
  typedef Kokkos::View<ordinal_type*, device_type>  device_ordinal_type_array;
  typedef Kokkos::View<cusolver_type*, device_type>   device_value_type_array;

  /// \name Constructor/Destructor methods
  //@{

  /**
   * \brief Initialize from Teuchos::RCP.
   *
   * \warning Should not be called directly!  Use instead
   * Amesos2::create() to initialize a cuSOLVER interface.
   */
  cuSOLVER(Teuchos::RCP<const Matrix> A,
          Teuchos::RCP<Vector>       X,
          Teuchos::RCP<const Vector> B);


  /// Destructor
  ~cuSOLVER( );

  //@}

private:

  /**
   * \brief Performs pre-ordering on the matrix to increase efficiency.
   */
  int preOrdering_impl();

  /**
   * \brief Perform symbolic factorization of the matrix using cuSOLVER.
   *
   * Called first in the sequence before numericFactorization.
   *
   * \throw std::runtime_error cuSOLVER is not able to factor the matrix.
   */
  int symbolicFactorization_impl();

  /**
   * \brief cuSOLVER specific numeric factorization
   *
   * \throw std::runtime_error cuSOLVER is not able to factor the matrix
   */
  int numericFactorization_impl();

  /**
   * \brief cuSOLVER specific solve.
   *
   * Uses the symbolic and numeric factorizations, along with the RHS
   * vector \c B to solve the sparse system of equations.  The
   * solution is placed in X.
   *
   * \throw std::runtime_error cuSOLVER is not able to solve the system.
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
   * Currently, the following cuSOLVER parameters/options are
   * recognized and acted upon:
   *
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
   * \param [in] cur rent_phase an indication of which solution phase this
   *                           load is being performed for.
   *
   * \return \c true if the matrix was loaded, \c false if not
   */
  bool loadA_impl(EPhase current_phase);

  /**
   * \brief can we optimize size_type and ordinal_type for straight pass through
   */
  bool do_optimization() const;

  // struct holds all data necessary to make a superlu factorization or solve call
  mutable struct cuSolverData {
    cusolverSpHandle_t handle;
    csrcholInfo_t chol_info;
    cusparseMatDescr_t desc;
    bool bReorder;
  } data_;

  typedef Kokkos::View<cusolver_type**, Kokkos::LayoutLeft, device_type> device_solve_array_t;

  mutable device_solve_array_t xValues_;
  mutable device_solve_array_t bValues_;
  mutable device_value_type_array buffer_;

  device_value_type_array device_nzvals_view_;
  device_size_type_array device_row_ptr_view_;
  device_ordinal_type_array device_cols_view_;
  size_t sorted_nnz;

  // data for reordering
  typedef Kokkos::View<ordinal_type*, device_type> permute_array_t;
  permute_array_t device_perm_;
  permute_array_t device_peri_;
  mutable device_solve_array_t permute_result_;
};                              // End class cuSOLVER

template <>
struct solver_traits<cuSOLVER> {
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
struct solver_supports_matrix<cuSOLVER,
  KokkosSparse::CrsMatrix<Scalar, LocalOrdinal, ExecutionSpace>> {
  static const bool value = true;
};

} // end namespace Amesos2

#endif  // AMESOS2_CUSOLVER_DECL_HPP
