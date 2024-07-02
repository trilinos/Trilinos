// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef AMESOS2_UMFPACK_DECL_HPP
#define AMESOS2_UMFPACK_DECL_HPP

#include "Amesos2_SolverTraits.hpp"
#include "Amesos2_SolverCore.hpp"
#include "Amesos2_Umfpack_FunctionMap.hpp"

namespace Amesos2 {


/** \brief Amesos2 interface to the Umfpack package.
 *
 * See the \ref umfpack_parameters "summary of Umfpack parameters"
 * supported by this Amesos2 interface.
 *
 * \ingroup amesos2_solver_interfaces
 */
template <class Matrix,
          class Vector>
class Umfpack : public SolverCore<Amesos2::Umfpack, Matrix, Vector>
{
  friend class SolverCore<Amesos2::Umfpack,Matrix,Vector>; // Give our base access
                                                          // to our private
                                                          // implementation funcs
public:

  /// Name of this solver interface.
  static const char* name;      // declaration. Initialization outside.

  typedef Umfpack<Matrix,Vector>                                       type;
  typedef SolverCore<Amesos2::Umfpack,Matrix,Vector>             super_type;

  // Since typedef's are not inheritted, go grab them
  typedef typename super_type::scalar_type                      scalar_type;
  typedef typename super_type::local_ordinal_type        local_ordinal_type;
  typedef typename super_type::global_ordinal_type      global_ordinal_type;
  typedef typename super_type::global_size_type            global_size_type;

  typedef TypeMap<Amesos2::Umfpack,scalar_type>                    type_map;

  /*
   * The Umfpack interface will need two other typedef's, which are:
   * - the umfpack type that corresponds to scalar_type and
   * - the corresponding type to use for magnitude
   */
  typedef typename type_map::type                              umfpack_type;
  typedef typename type_map::magnitude_type                  magnitude_type;

  typedef FunctionMap<Amesos2::Umfpack,umfpack_type>           function_map;

  typedef Kokkos::DefaultHostExecutionSpace              HostExecSpaceType;
  typedef Kokkos::View<int*, HostExecSpaceType>          host_size_type_array;
  typedef Kokkos::View<int*, HostExecSpaceType>          host_ordinal_type_array;
  typedef Kokkos::View<umfpack_type*, HostExecSpaceType> host_value_type_array;

  /// \name Constructor/Destructor methods
  //@{

  /**
   * \brief Initialize from Teuchos::RCP.
   *
   * \warning Should not be called directly!  Use instead
   * Amesos2::create() to initialize a Umfpack interface.
   */
  Umfpack(Teuchos::RCP<const Matrix> A,
          Teuchos::RCP<Vector>       X,
          Teuchos::RCP<const Vector> B);


  /// Destructor
  ~Umfpack( );

  //@}

  /// Returns a short description of this Solver
  std::string description() const;

private:

  /**
   * \brief Performs pre-ordering on the matrix to increase efficiency.
   *
   * Umfpack does not support pre-ordering, so this method does nothing.
   */
  int preOrdering_impl();


  /**
   * \brief Perform symbolic factorization of the matrix using Umfpack.
   *
   * Called first in the sequence before numericFactorization.
   *
   * \throw std::runtime_error Umfpack is not able to factor the matrix.
   */
  int symbolicFactorization_impl();


  /**
   * \brief Umfpack specific numeric factorization
   *
   * \throw std::runtime_error Umfpack is not able to factor the matrix
   */
  int numericFactorization_impl();

  /**
   * \brief Umfpack specific solve.
   *
   * Uses the symbolic and numeric factorizations, along with the RHS
   * vector \c B to solve the sparse system of equations.  The
   * solution is placed in X.
   *
   * \throw std::runtime_error Umfpack is not able to solve the system.
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

  // struct holds all data necessary to make a umfpack factorization or solve call
  mutable struct UMFPACKData {
    // Umfpack internal opaque object
    void *Symbolic;
    void *Numeric;

    // Info and Control state
    double Info[UMFPACK_INFO];
    double Control[UMFPACK_CONTROL];
  } data_;

  // The following Kokkos::View's are persisting storage for A's CCS arrays
  /// Stores the values of the nonzero entries for Umfpack
  host_value_type_array nzvals_view_;
  /// Stores the location in \c Ai_ and Aval_ that starts row j
  host_ordinal_type_array rowind_view_;
  /// Stores the row indices of the nonzero entries
  host_size_type_array colptr_view_;

  // The following Arrays are persisting storage arrays for X, and B
  /// Persisting 1D store for X
  Teuchos::Array<umfpack_type> xvals_;  int ldx_;
  /// Persisting 1D store for B
  Teuchos::Array<umfpack_type> bvals_;  int ldb_;

  bool is_contiguous_;
};                              // End class Umfpack


// Specialize solver_traits struct for Umfpack
template <>
struct solver_traits<Umfpack> {
#ifdef HAVE_TEUCHOS_COMPLEX
  typedef Meta::make_list4<float,
                           double,
                           std::complex<float>,
                           std::complex<double>
                           >supported_scalars;
#else
  typedef Meta::make_list2<float,
                           double
                           >supported_scalars;
#endif
};

} // end namespace Amesos2

#endif  // AMESOS2_UMFPACK_DECL_HPP
