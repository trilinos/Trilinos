// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/**
   \file   Amesos2_CssMKL_decl.hpp
   \author Eric Bavier <etbavie@sandia.gov>
   \date   Wed Jul 27 12:52:30 MDT 2011

   \brief  A template class that does nothing useful besides show developers
           what, in general, needs to be done to add a new solver interface to
           the Amesos2 collection.
*/


#ifndef AMESOS2_CSSMKL_DECL_HPP
#define AMESOS2_CSSMKL_DECL_HPP

#include <map>

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "Amesos2_SolverTraits.hpp"
#include "Amesos2_SolverCore.hpp"
#include "Amesos2_CssMKL_FunctionMap.hpp"


namespace Amesos2 {


  /** \brief Amesos2 interface to the CssMKL package.
   *
   * This class provides access to the MKL's distributed-memory
   * sparse direct solver, called Cluster Sparse Solver (CSS).
   * Access is provided for \c float and \c double scalar types, in
   * both real and complex.  Access to 32-bit or 64-bit integer routines
   * are provided based on the global_ordinal_type.
   *
   * \ingroup amesos2_solver_interfaces
   */
  template <class Matrix,
            class Vector>
  class CssMKL : public SolverCore<Amesos2::CssMKL, Matrix, Vector>
  {
    friend class SolverCore<Amesos2::CssMKL,Matrix,Vector>; // Give our base access
                                                            // to our private
                                                            // implementation funcs
  public:

    /// The name of this solver interface
    static const char* name;    // declaration. Initialization outside.

    typedef CssMKL<Matrix,Vector>                                       type;
    typedef SolverCore<Amesos2::CssMKL,Matrix,Vector>             super_type;

    // Since typedef's are not inheritted, go grab them
    typedef typename super_type::scalar_type                         scalar_type;
    typedef typename super_type::local_ordinal_type           local_ordinal_type;
    typedef typename super_type::global_ordinal_type         global_ordinal_type;
    typedef typename super_type::global_size_type               global_size_type;
    typedef typename super_type::node_type                             node_type;
    typedef Tpetra::Map<local_ordinal_type,
                                   global_ordinal_type,
                                   node_type>                           map_type;

    typedef TypeMap<Amesos2::PardisoMKL,scalar_type>                    type_map;

    typedef typename type_map::type                           solver_scalar_type;
    typedef typename type_map::magnitude_type              solver_magnitude_type;

    // This may be PMKL::_INTEGER_t or long long int depending on the
    // mapping and input ordinal
    typedef typename TypeMap<Amesos2::PardisoMKL,global_ordinal_type>::type int_t;

    /* For CssMKL we dispatch based on the integer type instead of
     * the scalar type:
     *   - _INTEGER_t    => use the cluster_sparse_solver(...)    method
     *   - long long int => use the cluster_sparse_solver_64(...) method
     */
    typedef FunctionMap<Amesos2::CssMKL,int_t>                  function_map;

    typedef Kokkos::DefaultHostExecutionSpace HostExecSpaceType;
    typedef Kokkos::View<int_t*,              HostExecSpaceType>    host_size_type_array;
    typedef Kokkos::View<int_t*,              HostExecSpaceType>    host_ordinal_type_array;
    typedef Kokkos::View<solver_scalar_type*, HostExecSpaceType>    host_value_type_array;

    /// \name Constructor/Destructor methods
    //@{

    /**
     * \brief Initialize from Teuchos::RCP.
     *
     * \warning Should not be called directly!  Use instead
     * Amesos2::create() to initialize a CssMKL interface.
     */
    CssMKL(Teuchos::RCP<const Matrix> A,
               Teuchos::RCP<Vector>       X,
               Teuchos::RCP<const Vector> B);


    /// Destructor
    ~CssMKL( );

    //@}

  private:

    /**
     * \brief Performs pre-ordering on the matrix to increase efficiency.
     *
     * CssMKL does reordering internally during symbolic
     * factorization.  Please refer to the \c "IPARM(2)" parameter for
     * some reordering options.
     */
    int preOrdering_impl();


    /**
     * \brief Perform symbolic factorization of the matrix using CssMKL.
     *
     * Called the sequence before numericFactorization.
     *
     * \throw std::runtime_error CssMKL is not able to factor the matrix.
     */
    int symbolicFactorization_impl();


    /**
     * \brief CssMKL specific numeric factorization
     *
     * \throw std::runtime_error CssMKL is not able to factor the matrix
     */
    int numericFactorization_impl();


    /**
     * \brief CssMKL specific solve.
     *
     * Uses the symbolic and numeric factorizations, along with the RHS vector
     * \c B to solve the sparse system of equations.
     *
     * The solution of the system is placed in X.
     *
     * \throw std::runtime_error CssMKL is not able to solve the system.
     */
    int solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
                   const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const;


    /**
     * \brief Determines whether the shape of the matrix is OK for this solver.
     *
     * CSS handles square matrices.
     */
    bool matrixShapeOK_impl() const;


    /**
     * The MKL CSS parameters that are currently recognized are:
     *
     * <ul>
     *  <li> \c "IPARM(2)"</li>: Fill-in reducing ordering for the input matrix.
     *                           2*: nested dissection algorithm from the METIS package,
     *                           3: thread-parallel version of the nested dissection algorithm, and
     *                           10: MPI-parallel version of the nested dissection.
     *  <li> \c "IPARM(8)"</li>: Iterative refinement step.
     *                           0*: performs two steps of iterative refinement when perturbed pivots are obtained during the numerical factorization.
     *                           >0: Maximum number of iterative refinement steps that the solver performs.
     *                           <0: Same as above, but the accumulation of the residuum uses extended precision.
     *  <li> \c "IPARM(10)"</li>: Pivoting perturbation.
     *                            8: default value for symmetric indefinite matrices
     *                            13: default value for nonsymmetric matrices
     *  <li> \c "IPARM(12)"</li>: Solve with transposed or conjugate transposed matrix A.
     *                            0*: Solve a linear system AX = B.
     *                            1: Solve a conjugate transposed system A^HX = B 
     *                            2: Solve a transposed system A^TX = B
     *  <li> \c "IPARM(13)"</li>: Improved accuracy using (non-) symmetric weighted matching.
     *                            0: Disable matching. Default for symmetric indefinite matrices.
     *                            1: nable matching. Default for nonsymmetric matrices.
     *  <li> \c "IPARM(18)"</li>: Report the number of non-zero elements in the factors.
     *                            <0*: Enable reporting.
     *                            >=0: Disable reporting.
     *  <li> \c "IPARM(27)"</li>: Matrix checker.
     *                            0*:Do not check the sparse matrix representation for error.
     *                            1: Check integer arrays ia and ja. In particular, check whether the column indices are sorted in increasing order within each row.
     * </ul>
     *
     */
    void setParameters_impl(const Teuchos::RCP<Teuchos::ParameterList> & parameterList );


    /**
     * \return a const Teuchos::ParameterList of all valid parameters
     * (set to their default values) for this solver.
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


    ////////// Internal routines (not called from outside) //////////

    /** \internal
     *
     * \brief Throws an appropriate runtime error in the event that
     * error < 0 .
     *
     * \param phase the phase for which this error is being checked.
     *              The meaning of a particular error value may depend
     *              on which phase was last performed
     *
     * \param error the error value returned by CssMKL for the
     *              given phase.
     *
     * We broadcast the input value from the rank=0 image to all
     * others before checking the value.  Before doing this we convert
     * the error into an \c int value which allow us to easily
     * broadcast its value to all process images without having to
     * enable Teuchos long long support in the case where the user is
     * making use of pardiso_64.  The valid values of error certainly
     * fit within an int.
     */
    void check_css_mkl_error(EPhase phase, int_t error) const;

    /** \internal
     *
     * Sets the internal mtype_ member.  Errors are thrown for
     * unacceptable scalar/mtype combinations.
     *
     * \param mtype the type of the matrix.  This may come as input
     * from the interface user, or may be set to the default value in
     * case mtype == 0 on entry to this function.
     */
    void set_css_mkl_matrix_type(int_t mtype = 0);
    void set_css_mkl_default_parameters(void* pt[], int_t iparm[]) const;


    /* Declare private variables necessary for interaction with the
     * CssMKL TPL.
     *
     * For example, the following Arrays are persisting storage arrays
     * for A, X, and B that can be used with solvers expecting a
     * compressed-row representation of the matrix A.
     */

    /// Stores the values of the nonzero entries for CssMKL
    host_value_type_array nzvals_view_;
    host_value_type_array nzvals_temp_;
    /// Stores the location in \c Ai_ and Aval_ that starts row j
    host_ordinal_type_array colind_view_;
    /// Stores the row indices of the nonzero entries
    host_size_type_array rowptr_view_;
    /// Persisting, contiguous, 1D store for X
    mutable Teuchos::Array<solver_scalar_type> xvals_;
    /// Persisting, contiguous, 1D store for B
    mutable Teuchos::Array<solver_scalar_type> bvals_;

    /// CssMKL internal data address pointer
    mutable void* pt_[64];
    /// The matrix type.  We deal only with unsymmetrix matrices
    int_t mtype_;
    /// Number of equations in the sparse linear system
    int_t n_;
    /// Permutation vector
    Teuchos::Array<int_t> perm_;
    /// number of righthand-side vectors
    mutable int_t nrhs_;

    bool css_initialized_;
    bool is_contiguous_;

    /// The messaging level.  Set to 1 if you wish for CSS to print statistical info
    int_t msglvl_;

    /// CssMKL parameter vector.  Note that the documentation uses
    /// 1-based indexing, but our interface must use 0-based indexing
    int_t iparm_[64];

    // We will deal with 1 factor at a time
    static const int_t maxfct_;
    static const int_t mnum_;


    static const bool complex_
    = Meta::or_<std::is_same_v<solver_scalar_type, PMKL::_MKL_Complex8>,
                std::is_same_v<solver_scalar_type, PMKL::_DOUBLE_COMPLEX_t>>::value;

    MPI_Fint CssComm_;
    Teuchos::RCP<const map_type> css_rowmap_;
    Teuchos::RCP<const map_type> css_contig_rowmap_;
    Teuchos::RCP<const map_type> css_contig_colmap_;

};                              // End class CssMKL


// Specialize the solver_traits struct for CssMKL.
template <>
struct solver_traits<CssMKL> {
#ifdef HAVE_TEUCHOS_COMPLEX
  typedef Meta::make_list6<float,
                           double,
                           std::complex<float>,
                           std::complex<double>,
                           PMKL::_MKL_Complex8,
                           PMKL::_DOUBLE_COMPLEX_t> supported_scalars;
#else
typedef Meta::make_list2<float,
                         double> supported_scalars;
#endif
};

} // end namespace Amesos

#endif  // AMESOS2_CSSMKL_DECL_HPP
