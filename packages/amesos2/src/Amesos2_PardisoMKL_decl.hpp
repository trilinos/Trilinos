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
   \file   Amesos2_PardisoMKL_decl.hpp
   \author Eric Bavier <etbavie@sandia.gov>
   \date   Wed Jul 27 12:52:30 MDT 2011

   \brief  A template class that does nothing useful besides show developers
           what, in general, needs to be done to add a new solver interface to
           the Amesos2 collection.
*/


#ifndef AMESOS2_PARDISOMKL_DECL_HPP
#define AMESOS2_PARDISOMKL_DECL_HPP

#include <map>

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "Amesos2_SolverTraits.hpp"
#include "Amesos2_SolverCore.hpp"
#include "Amesos2_PardisoMKL_FunctionMap.hpp"


namespace Amesos2 {


  /** \brief Amesos2 interface to the PardisoMKL package.
   *
   * This class provides access to the Pardiso (MKL version 10.3 and
   * compatible) sparse direct solver with out-of-core solve support.
   * Access is provided for \c float and \c double scalar types, in
   * both real and complex.  Access to to Pardiso's 64-bit integer
   * routines is also provided.
   *
   * \ingroup amesos2_solver_interfaces
   */
  template <class Matrix,
            class Vector>
  class PardisoMKL : public SolverCore<Amesos2::PardisoMKL, Matrix, Vector>
  {
    friend class SolverCore<Amesos2::PardisoMKL,Matrix,Vector>; // Give our base access
                                                                // to our private
                                                                // implementation funcs
  public:

    /// The name of this solver interface
    static const char* name;    // declaration. Initialization outside.

    typedef PardisoMKL<Matrix,Vector>                                       type;
    typedef SolverCore<Amesos2::PardisoMKL,Matrix,Vector>             super_type;

    // Since typedef's are not inheritted, go grab them
    typedef typename super_type::scalar_type                         scalar_type;
    typedef typename super_type::local_ordinal_type           local_ordinal_type;
    typedef typename super_type::global_ordinal_type         global_ordinal_type;
    typedef typename super_type::global_size_type               global_size_type;

    typedef TypeMap<Amesos2::PardisoMKL,scalar_type>                    type_map;

    typedef typename type_map::type                           solver_scalar_type;
    typedef typename type_map::magnitude_type              solver_magnitude_type;

    // This may be PMKL::_INTEGER_t or long long int depending on the
    // mapping and input ordinal
    typedef typename TypeMap<Amesos2::PardisoMKL,local_ordinal_type>::type int_t;

    /* For PardisoMKL we dispatch based on the integer type instead of
     * the scalar type:
     *   - _INTEGER_t    => use the pardiso(...)    method
     *   - long long int => use the pardiso_64(...) method
     */
    typedef FunctionMap<Amesos2::PardisoMKL,int_t>                  function_map;


    /// \name Constructor/Destructor methods
    //@{

    /**
     * \brief Initialize from Teuchos::RCP.
     *
     * \warning Should not be called directly!  Use instead
     * Amesos2::create() to initialize a PardisoMKL interface.
     */
    PardisoMKL(Teuchos::RCP<const Matrix> A,
               Teuchos::RCP<Vector>       X,
               Teuchos::RCP<const Vector> B);


    /// Destructor
    ~PardisoMKL( );

    //@}

  private:

    /**
     * \brief Performs pre-ordering on the matrix to increase efficiency.
     *
     * PardisoMKL does reordering internally during symbolic
     * factorization.  Please refer to the \c "IPARM(2)" parameter for
     * some reordering options.
     */
    int preOrdering_impl();


    /**
     * \brief Perform symbolic factorization of the matrix using PardisoMKL.
     *
     * Called the sequence before numericFactorization.
     *
     * \throw std::runtime_error PardisoMKL is not able to factor the matrix.
     */
    int symbolicFactorization_impl();


    /**
     * \brief PardisoMKL specific numeric factorization
     *
     * \throw std::runtime_error PardisoMKL is not able to factor the matrix
     */
    int numericFactorization_impl();


    /**
     * \brief PardisoMKL specific solve.
     *
     * Uses the symbolic and numeric factorizations, along with the RHS vector
     * \c B to solve the sparse system of equations.
     *
     * The solution of the system is placed in X.
     *
     * \throw std::runtime_error PardisoMKL is not able to solve the system.
     */
    int solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
                   const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const;


    /**
     * \brief Determines whether the shape of the matrix is OK for this solver.
     *
     * Pardiso MKL handles square matrices.
     */
    bool matrixShapeOK_impl() const;


    /**
     * The Pardiso MKL parameters that are currently recognized are:
     *
     * <ul>
     *  <li> \c "IPARM(2)"</li>
     *  <li> \c "IPARM(4)"</li>
     *  <li> \c "IPARM(8)"</li>
     *  <li> \c "IPARM(10)"</li>
     *  <li> \c "IPARM(18)"</li>
     *  <li> \c "IPARM(24)"</li>
     *  <li> \c "IPARM(25)"</li>
     *  <li> \c "IPARM(60)"</li>
     * </ul>
     *
     * Please see the Pardiso MKL documentation for a summary of the
     * meaning and valid values for each parameter.
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
     * \param error the error value returned by PardisoMKL for the
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
    void check_pardiso_mkl_error(EPhase phase, int_t error) const;

    /** \internal
     *
     * Sets the internal mtype_ member.  Errors are thrown for
     * unacceptable scalar/mtype combinations.
     *
     * \param mtype the type of the matrix.  This may come as input
     * from the interface user, or may be set to the default value in
     * case mtype == 0 on entry to this function.
     */
    void set_pardiso_mkl_matrix_type(int_t mtype = 0);


    /* Declare private variables necessary for interaction with the
     * PardisoMKL TPL.
     *
     * For example, the following Arrays are persisting storage arrays
     * for A, X, and B that can be used with solvers expecting a
     * compressed-row representation of the matrix A.
     */

    /// Stores the values of the nonzero entries for PardisoMKL
    Teuchos::Array<solver_scalar_type> nzvals_;
    /// Stores the location in \c Ai_ and Aval_ that starts row j
    Teuchos::Array<int_t> colind_;
    /// Stores the row indices of the nonzero entries
    Teuchos::Array<int_t> rowptr_;
    /// Persisting, contiguous, 1D store for X
    mutable Teuchos::Array<solver_scalar_type> xvals_;
    /// Persisting, contiguous, 1D store for B
    mutable Teuchos::Array<solver_scalar_type> bvals_;

    /// PardisoMKL internal data address pointer
    mutable void* pt_[64];
    /// The matrix type.  We deal only with unsymmetrix matrices
    int_t mtype_;
    /// Number of equations in the sparse linear system
    int_t n_;
    /// Permutation vector
    Teuchos::Array<int_t> perm_;
    /// number of righthand-side vectors
    mutable int_t nrhs_;

    /// PardisoMKL parameter vector.  Note that the documentation uses
    /// 1-based indexing, but our interface must use 0-based indexing
    int_t iparm_[64];

    /// The messaging level.  Set to 1 if you wish for Pardiso MKL to print statistical info
    static const int_t msglvl_;

    // We will deal with 1 factor at a time
    static const int_t maxfct_;
    static const int_t mnum_;


    static const bool complex_
    = Meta::or_<Meta::is_same<solver_scalar_type, PMKL::_MKL_Complex8>::value,
                Meta::is_same<solver_scalar_type, PMKL::_DOUBLE_COMPLEX_t>::value>::value;

    mutable std::map<int,Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> > > validators;

};                              // End class PardisoMKL


// Specialize the solver_traits struct for PardisoMKL.
template <>
struct solver_traits<PardisoMKL> {
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

#endif  // AMESOS2_PARDISOMKL_DECL_HPP
