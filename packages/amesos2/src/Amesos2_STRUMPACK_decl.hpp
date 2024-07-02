// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef AMESOS2_STRUMPACK_DECL_HPP
#define AMESOS2_STRUMPACK_DECL_HPP

#include "Amesos2_SolverTraits.hpp"
#include "Amesos2_SolverCore.hpp"

#ifdef HAVE_MPI
#include "StrumpackSparseSolverMPIDist.hpp"
#else
#include "StrumpackSparseSolver.hpp"
#endif

namespace Amesos2 {


/** \brief Amesos2 interface to STRUMPACK direct solver and preconditioner.
 *
 * Currently support is for the STRUMPACK 2.1.0 version.
 *
 * \warning After creation, the size of the matrix should not change
 * (i.e. when using setA())
 *
 * \ingroup amesos2_solver_interfaces
 */
template <class Matrix,
          class Vector>
class STRUMPACK : public SolverCore<Amesos2::STRUMPACK, Matrix, Vector>
{
  friend class SolverCore<Amesos2::STRUMPACK,Matrix,Vector>; // Give our base access
                                                             // to our private
                                                             // implementation funcs
public:

  /// Name of this solver interface.
  static const char* name;      // declaration. Initialization outside.

  typedef STRUMPACK<Matrix,Vector>                                     type;
  typedef SolverCore<Amesos2::STRUMPACK,Matrix,Vector>           strum_type;

  typedef Matrix                                                matrix_type;
  typedef Vector                                                vector_type;

  // Since typedef's are not inheritted, go grab them
  typedef typename strum_type::scalar_type                      scalar_type;
  typedef typename strum_type::local_ordinal_type        local_ordinal_type;
  typedef typename strum_type::global_ordinal_type      global_ordinal_type;
  typedef typename strum_type::global_size_type            global_size_type;
  typedef typename strum_type::node_type                          node_type;

  typedef Kokkos::DefaultHostExecutionSpace                     HostExecSpaceType;
  typedef Kokkos::View<global_ordinal_type*, HostExecSpaceType> host_ordinal_type_array;
  typedef Kokkos::View<scalar_type*, HostExecSpaceType>         host_value_type_array;

  /// \name Constructor/Destructor methods
  //@{

  /**
   * \brief Initialize from Teuchos::RCP.
   *
   * \warning Should not be called directly!  Use instead
   * Amesos2::create() to initialize a STRUMPACK interface.
   */
  STRUMPACK(Teuchos::RCP<const Matrix> A,
            Teuchos::RCP<Vector>       X,
            Teuchos::RCP<const Vector> B);


  /// Destructor
  ~STRUMPACK( );

  //@}

private:

  /**
   * \brief Performs pre-ordering on the matrix to increase efficiency.
   *
   * STRUMPACK supports several forms of column permutations.
   * TODO document different ordering?
   */
  int preOrdering_impl();


  /**
   * \brief Perform symbolic factorization of the matrix using STRUMPACK.
   *
   * Called second in the sequence before numericFactorization.
   *
   * \throw std::runtime_error STRUMPACK is not able to factor the matrix.
   */
  int symbolicFactorization_impl();


  /**
   * \brief STRUMPACK specific numeric factorization
   *
   * \throw std::runtime_error STRUMPACK is not able to factor the matrix
   */
  int numericFactorization_impl();


  /**
   * \brief STRUMPACK specific solve.
   *
   * Uses the symbolic and numeric factorizations, along with the RHS
   * vector \c B to solve the sparse system of equations.  The
   * solution is placed in X.
   *
   * \throw std::runtime_error STRUMPACK is not able to solve the system.
   *
   * \callgraph
   */
  int solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> > X,
                 const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const;


  /**
   * \brief Determines whether the shape of the matrix is OK for this solver.
   *
   * STRUMPACK supports square matrices.
   */
  bool matrixShapeOK_impl() const;


  /**
   * Currently, the following STRUMPACK parameters/options are recognized:
   *
   * <ul>
   *   <li> \c "Ordering" which takes one of the following:
   *     <ul>
   *     <li> \c "NATURAL" : natural reordering.</li>
   *     <li> \c "METIS" : use METIS nested-dissection reordering.</li>
   *     <li> \c "PARMETIS" : use the ParMETIS nested-dissection reordering. (default)</li>
   *     <li> \c "SCOTCH" : use Scotch nested-dissection reordering.</li>
   *     <li> \c "PTSCOTCH" : use PT-Scotch nested-dissection reordering.</li>
   *     <li> \c "RCM" : use RCM reordering.</li>
   *     <li> \c "GEOMETRIC" : a simple geometric nested dissection code that only works for regular meshes.</li>
   *     </ul>
   *   <li> \c "Matching" which takes one of the following:
   *     <ul>
   *     <li> \c "NONE" : don't do anything. (default)</li>
   *     <li> \c "MAX_CARDINALITY" : Maximum cardinality.</li>
   *     <li> \c "MAX_SMALLEST_DIAGONAL" : Maximum smallest diagonal value.</li>
   *     <li> \c "MAX_SMALLEST_DIAGONAL_2" : Same as MAX_SMALLEST_DIAGONAL, but different algorithm.</li>
   *     <li> \c "MAX_DIAGONAL_SUM" : Maximum sum of diagonal values.</li>
   *     <li> \c "MAX_DIAGONAL_PRODUCT_SCALING" : Maximum product of diagonal values and row and column scaling.</li>
   *     <li> \c "COMBBLAS" : Use AWPM from CombBLAS.</li>
   *     </ul>
   *   <li> \c "ReplaceTinyPivot" which takes one of the following:
   *     <ul>
   *     <li> \c true : replace tiny diagonals during LU factorization.</li>
   *     <li> \c false : do not replace tiny diagonals during LU factorization.</li>
   *     </ul>
   *   </li>
   *
   * </ul>
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
   * \brief Reads matrix data into internal solver structures.
   *
   * Loads data from the matrix A into the internal STRUMPACK matrix
   * structure.  This function requires communication accross all
   * processors as the matrix is redistributed as necessary in
   * STRUMPACK internally.
   *
   * \post
   * - nzvals_, colind_, and rowptr_ arrays are sized to match the portion
   *   of the matrix on this processor.
   *
   * \return \c true if the matrix was loaded, \c false if not
   */
  bool loadA_impl(EPhase current_phase);

#ifdef HAVE_MPI
  Teuchos::RCP<strumpack::StrumpackSparseSolverMPIDist<scalar_type,global_ordinal_type>> sp_;
#else
  Teuchos::RCP<strumpack::StrumpackSparseSolver<scalar_type,global_ordinal_type>> sp_;
#endif


  // The following Arrays are persisting storage arrays for A, X, and B
  /// Stores the values of the nonzero entries for STRUMPACK
  host_value_type_array nzvals_view_;
  /// Stores the row indices of the nonzero entries
  host_ordinal_type_array colind_view_;
  /// Stores the location in \c Ai_ and Aval_ that starts row j
  host_ordinal_type_array rowptr_view_;
  // /// 1D store for B values
   mutable Teuchos::Array<scalar_type> bvals_;
  // /// 1D store for X values
   mutable Teuchos::Array<scalar_type> xvals_;

#ifdef HAVE_MPI
  /// Maps rows of the matrix to processors in the STRUMPACK processor grid
  Teuchos::RCP<const Tpetra::Map<local_ordinal_type,
                                 global_ordinal_type,
                                 node_type> > strumpack_rowmap_;
#endif
};                              // End class STRUMPACK


// Specialize the solver_traits template for STRUMPACK
template <>
struct solver_traits<STRUMPACK> {
#if defined(HAVE_TEUCHOS_COMPLEX) && !defined(__clang__)
  typedef Meta::make_list4<float, double, std::complex<float>, std::complex<double>> supported_scalars;
#else
  typedef Meta::make_list2<float, double> supported_scalars;
#endif
};

} // end namespace Amesos2

#endif  // AMESOS2_STRUMPACK_DECL_HPP
