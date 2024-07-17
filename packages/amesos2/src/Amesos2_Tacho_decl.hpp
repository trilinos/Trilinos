// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef AMESOS2_TACHO_DECL_HPP
#define AMESOS2_TACHO_DECL_HPP

#include <Kokkos_Core.hpp>

#include "Amesos2_SolverTraits.hpp"
#include "Amesos2_SolverCore.hpp"
#include "Amesos2_Tacho_FunctionMap.hpp"

#include "Tacho.hpp"
#include "Tacho_Solver.hpp"

namespace Amesos2 {

/** \brief Amesos2 interface to the Tacho package.
 *
 * See the \ref tacho_parameters "summary of Tacho parameters"
 * supported by this Amesos2 interface.
 *
 * \ingroup amesos2_solver_interfaces
 */
template <class Matrix,
          class Vector>
class TachoSolver : public SolverCore<Amesos2::TachoSolver, Matrix, Vector>
{
  friend class SolverCore<Amesos2::TachoSolver,Matrix,Vector>; // Give our base access
                                                          // to our private
                                                          // implementation funcs
public:

  /// Name of this solver interface.
  static const char* name;      // declaration. Initialization outside.

  typedef TachoSolver<Matrix,Vector>                                   type;
  typedef SolverCore<Amesos2::TachoSolver,Matrix,Vector>         super_type;

  // Since typedef's are not inheritted, go grab them
  typedef typename super_type::scalar_type                      scalar_type;
  typedef typename super_type::local_ordinal_type        local_ordinal_type;
  typedef typename super_type::global_size_type            global_size_type;

  typedef TypeMap<Amesos2::TachoSolver,scalar_type>                type_map;

  /*
   * The Tacho interface will need two other typedef's, which are:
   * - the tacho type that corresponds to scalar_type and
   * - the corresponding type to use for magnitude
   */
  typedef typename type_map::type                                tacho_type;
  typedef typename type_map::magnitude_type                  magnitude_type;

  typedef FunctionMap<Amesos2::TachoSolver,tacho_type>         function_map;

  // TODO - Not sure yet best place for organizing these typedefs
  typedef Tacho::ordinal_type                                  ordinal_type;
  typedef Tacho::size_type                                        size_type;

  typedef Kokkos::DefaultExecutionSpace                     exec_space_type;
  typedef Kokkos::DefaultHostExecutionSpace            host_exec_space_type;

  typedef typename
    Tacho::UseThisDevice<exec_space_type>::device_type          device_type;
  typedef typename
    Tacho::UseThisDevice<host_exec_space_type>::device_type host_device_type;

  typedef Tacho::DummyTaskScheduler<exec_space_type>          scheduler_type;

  typedef Kokkos::View<size_type*, device_type>       device_size_type_array;
  typedef Kokkos::View<ordinal_type*, device_type> device_ordinal_type_array;
  typedef Kokkos::View<tacho_type*, device_type>     device_value_type_array;

  // also work with host space - right now symbolic requires host space so we
  // do everything in device space if source was device, then deep_copy
  typedef Kokkos::View<size_type*, host_device_type>       host_size_type_array;
  typedef Kokkos::View<ordinal_type*, host_device_type> host_ordinal_type_array;

  /// \name Constructor/Destructor methods
  //@{

  /**
   * \brief Initialize from Teuchos::RCP.
   *
   * \warning Should not be called directly!  Use instead
   * Amesos2::create() to initialize a Tacho interface.
   */
  TachoSolver(Teuchos::RCP<const Matrix> A,
          Teuchos::RCP<Vector>       X,
          Teuchos::RCP<const Vector> B);


  /// Destructor
  ~TachoSolver( );

  //@}

  /// Returns a short description of this Solver
  std::string description() const override;

private:

  /**
   * \brief Performs pre-ordering on the matrix to increase efficiency.
   *
   * Tacho does not support pre-ordering, so this method does nothing.
   */
  int preOrdering_impl();


  /**
   * \brief Perform symbolic factorization of the matrix using Tacho.
   *
   * Called first in the sequence before numericFactorization.
   *
   * \throw std::runtime_error Tacho is not able to factor the matrix.
   */
public: // made this public for CUDA parallel_for usage
  int symbolicFactorization_impl();

private:

  /**
   * \brief Tacho specific numeric factorization
   *
   * \throw std::runtime_error Tacho is not able to factor the matrix
   */
  int numericFactorization_impl();

  /**
   * \brief Tacho specific solve.
   *
   * Uses the symbolic and numeric factorizations, along with the RHS
   * vector \c B to solve the sparse system of equations.  The
   * solution is placed in X.
   *
   * \throw std::runtime_error Tacho is not able to solve the system.
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


  /**
   * \brief can we optimize size_type and ordinal_type for straight pass through
   */
  bool do_optimization() const;

  // struct holds all data necessary to make a tacho factorization or solve call
  mutable struct TACHOData {
    typename Tacho::Solver<tacho_type, scheduler_type> solver;

    // TODO: Implement the paramter options - confirm which we want and which have been implemented
    int method;
    int variant;
    int small_problem_threshold_size;
    // int num_kokkos_threads;
    // int max_num_superblocks;
  } data_;

  typedef typename Tacho::Solver<tacho_type, scheduler_type>::value_type_matrix
    device_solve_array_t;

  // used as an internal workspace - possibly we can store this better in TACHOData
  mutable device_solve_array_t workspace_;

  // x and b for solve - only allocated first time
  mutable device_solve_array_t xValues_;
  mutable device_solve_array_t bValues_;

  // numeric is on device
  device_value_type_array device_nzvals_view_;

  // symbolic is done on host for Tacho so store these versions as well
  host_size_type_array host_row_ptr_view_;
  host_ordinal_type_array host_cols_view_;
};                              // End class Tacho


// Specialize solver_traits struct for Tacho
template <>
struct solver_traits<TachoSolver> {
#ifdef HAVE_TEUCHOS_COMPLEX
  typedef Meta::make_list6<float,
                           double,
                           std::complex<float>,
                           std::complex<double>,
                           Kokkos::complex<float>,
                           Kokkos::complex<double>
                           >supported_scalars;
#else
  typedef Meta::make_list2<float,
                           double
                           >supported_scalars;
#endif
};

template <typename Scalar, typename LocalOrdinal, typename ExecutionSpace>
struct solver_supports_matrix<TachoSolver,
  KokkosSparse::CrsMatrix<Scalar, LocalOrdinal, ExecutionSpace>> {
  static const bool value = true;
};

} // end namespace Amesos2

#endif  // AMESOS2_TACHO_DECL_HPP
