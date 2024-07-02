// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
   \file   Amesos2_ShyLUBasker_decl.hpp
   \author Joshua Dennis Booth <jdbooth@sandia.gov>
           Siva Rajamanickam <srajama@sandia.gov>
           Nathan Ellingwood <ndellin@sandia.gov>

   \brief  Amesos2 ShyLUBasker declarations.
*/


#ifndef AMESOS2_SHYLUBASKER_DECL_HPP
#define AMESOS2_SHYLUBASKER_DECL_HPP

#include "Amesos2_SolverTraits.hpp"
#include "Amesos2_SolverCore.hpp"
#include "Amesos2_ShyLUBasker_FunctionMap.hpp"

#include "shylubasker_decl.hpp"
#include "shylubasker_def.hpp"
#include "shylubasker_trilinos_decl.hpp"


namespace Amesos2 {

/** \brief Amesos2 interface to the Baker package.
 *
 * See the \ref ShyLUBasker_parameters "summary of ShyLUBasker parameters"
 * supported by this Amesos2 interface.
 *
 * \ingroup amesos2_solver_interfaces
 */
template <class Matrix,class Vector>
class ShyLUBasker : public SolverCore<Amesos2::ShyLUBasker, Matrix, Vector>
{
  friend class SolverCore<Amesos2::ShyLUBasker,Matrix,Vector>; // Give our base access
                                                          // to our private
                                                          // implementation funcs
public:

  /// Name of this solver interface.
  static const char* name;      // declaration. Initialization outside.


  typedef ShyLUBasker<Matrix,Vector>                           type;

  typedef SolverCore<Amesos2::ShyLUBasker,Matrix,Vector>       super_type;

  // Since typedef's are not inheritted, go grab them
  typedef typename super_type::scalar_type                     scalar_type;
  typedef typename super_type::local_ordinal_type              local_ordinal_type;
  typedef typename super_type::global_ordinal_type             global_ordinal_type;
  typedef typename super_type::global_size_type                global_size_type;
  typedef typename super_type::node_type                       node_type;

  typedef TypeMap<Amesos2::ShyLUBasker,scalar_type>            type_map;
  typedef typename type_map::type                              shylubasker_type;
  typedef typename type_map::dtype                             shylubasker_dtype;

  typedef typename type_map::type                              slu_type;

  typedef FunctionMap<Amesos2::ShyLUBasker,shylubasker_type>   function_map;

  typedef Matrix                                               matrix_type;
  typedef MatrixAdapter<matrix_type>                           matrix_adapter_type;

  typedef Kokkos::DefaultHostExecutionSpace                    HostExecSpaceType;
//  typedef Kokkos::View<local_ordinal_type*, HostExecSpaceType> host_size_type_array;
  typedef Kokkos::View<local_ordinal_type*, HostExecSpaceType> host_ordinal_type_array;
  typedef Kokkos::View<shylubasker_type*, HostExecSpaceType>   host_value_type_array;


  ShyLUBasker( Teuchos::RCP<const Matrix> A,
          Teuchos::RCP<Vector>       X,
          Teuchos::RCP<const Vector> B);
   ~ShyLUBasker( );


private:

 /**
  * \brief can we optimize size_type and ordinal_type for straight pass through,
  * also check that is_contiguous_ flag set to true
  */
  bool single_proc_optimization() const;


  /**
   * \brief Performs pre-ordering on the matrix to increase efficiency.
   *
   *   Come back to add support to Amesos for preordering
   */
  int preOrdering_impl();


  int symbolicFactorization_impl();


  /**
   * \brief ShyLUBasker specific numeric factorization
   *
   * \throw std::runtime_error ShyLUBasker is not able to factor the matrix
   */
  int numericFactorization_impl();


  /**
   * \brief ShyLUBasker specific solve.
   *
   * Uses the symbolic and numeric factorizations, along with the RHS
   * vector \c B to solve the sparse system of equations.  The
   * solution is placed in X.
   *
   * \throw std::runtime_error ShyLUBasker is not able to solve the system.
   *
   * \callgraph
   */
  int solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
                 const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const;


  /**
   * \brief Determines whether the shape of the matrix is OK for this solver.
   */
  bool matrixShapeOK_impl() const;


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


  // Members
  int num_threads;

  // The following Kokkos::View's are persisting storage for A's CCS arrays
  /// Stores the values of the nonzero entries for Umfpack
  host_value_type_array nzvals_view_;
  /// Stores the location in \c Ai_ and Aval_ that starts row j
  host_ordinal_type_array rowind_view_;
  /// Stores the row indices of the nonzero entries
  host_ordinal_type_array colptr_view_;


  bool is_contiguous_;

  typedef typename Kokkos::View<shylubasker_type**, Kokkos::LayoutLeft, 
                                typename HostExecSpaceType::memory_space> host_solve_array_t;

  /// Persisting 1D store for X
  mutable host_solve_array_t xValues_;
  int ldx_;

  /// Persisting 1D store for B
  mutable host_solve_array_t bValues_;
  int ldb_;

    /*Handle for ShyLUBasker object*/
 
#if defined( HAVE_AMESOS2_KOKKOS ) && defined( KOKKOS_ENABLE_OPENMP )
  /*
  typedef typename node_type::device_type  kokkos_device;
  typedef typename kokkos_device::execution_space kokkos_exe;
  static_assert(std::is_same<kokkos_exe,Kokkos::OpenMP>::value,
  "Kokkos node type not support by experimental ShyLUBasker Amesos2");
  */
  typedef Kokkos::OpenMP Exe_Space;
   ::BaskerNS::BaskerTrilinosInterface<local_ordinal_type, shylubasker_dtype, Exe_Space> *ShyLUbasker;
#else
  #pragma message("Amesos_ShyLUBasker_decl Error: ENABLED SHYLU_NODEBASKER BUT NOT KOKKOS or NOT OPENMP!")
#endif


}; // End class ShyLUBasker


// Specialize solver_traits struct for ShyLUBasker
// TODO
template <>
struct solver_traits<ShyLUBasker> {
#ifdef HAVE_TEUCHOS_COMPLEX
/*
  typedef Meta::make_list4<float,
                           double,
                           std::complex<float>,
                           std::complex<double> > supported_scalars;
*/
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
struct solver_supports_matrix<ShyLUBasker,
  KokkosSparse::CrsMatrix<Scalar, LocalOrdinal, ExecutionSpace>> {
  static const bool value = true;
};

} // end namespace Amesos2

#endif  // AMESOS2_SHYLUBASKER_DECL_HPP
