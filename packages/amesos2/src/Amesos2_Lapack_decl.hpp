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
   \file   Amesos2_Lapack_decl.hpp
   \author Eric Bavier <etbavie@sandia.gov>
   \date   Fri Jul 29 15:43:28 MDT 2011
   
   \brief  Declarations for the Amesos2 interface to LAPACK
*/


#ifndef AMESOS2_LAPACK_DECL_HPP
#define AMESOS2_LAPACK_DECL_HPP

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_SerialDenseSolver.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>

#include "Amesos2_SolverTraits.hpp"
#include "Amesos2_SolverCore.hpp"


namespace Amesos2 {


  /** \brief Amesos2 interface to the LAPACK.
   *
   * This Amesos2 interface uses a sequential, dense solver based on
   * the LAPACK routines.  This interface exists mostly for
   * pedagogical purposes (to contrast sparse and dense solver
   * performance), but could be useful in cases where the input matrix
   * is sparse but "almost dense".
   *
   * \ingroup amesos2_solver_interfaces 
   */
  template <class Matrix,
	    class Vector>
  class Lapack : public SolverCore<Amesos2::Lapack, Matrix, Vector> 
  {
    friend class SolverCore<Amesos2::Lapack,Matrix,Vector>; // Give our base access
                                                            // to our private
                                                            // implementation funcs
  public:

    /// The name of this solver interface
    static const char* name;	// declaration. Initialization outside.

    typedef Lapack<Matrix,Vector>                                                type;
    typedef SolverCore<Amesos2::Lapack,Matrix,Vector>                      super_type;

    // Since typedef's are not inheritted, go grab them
    typedef typename super_type::scalar_type                              scalar_type;
    typedef typename super_type::local_ordinal_type                local_ordinal_type;
    typedef typename super_type::global_ordinal_type              global_ordinal_type;
    typedef typename super_type::global_size_type                    global_size_type;

    typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;



    /// \name Constructor/Destructor methods
    //@{ 

    /**
     * \brief Initialize from Teuchos::RCP.
     *
     * \warning Should not be called directly!  Use instead
     * Amesos2::create() to initialize a Lapack interface.
     */
    Lapack(Teuchos::RCP<const Matrix> A,
	   Teuchos::RCP<Vector>       X,
	   Teuchos::RCP<const Vector> B);
    

    /// Destructor
    ~Lapack( );

    //@}

  private:
  
    /**
     * \brief No-op.
     */
    int preOrdering_impl();


    /**
     * \brief No-op
     */
    int symbolicFactorization_impl();


    /**
     * \brief Perform numeric factorization using LAPACK
     * 
     * \throw std::runtime_error Lapack is not able to factor the matrix
     */
    int numericFactorization_impl();


    /**
     * \brief Lapack solve.
     *
     * Uses the numeric factorization, along with the RHS vector \c B
     * to solve the system of equations.
     *
     * The solution of the system is placed in X.
     *
     * \throw std::runtime_error Lapack is not able to solve the system.
     */
    int solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
		   const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const;


    /**
     * \brief Determines whether the shape of the matrix is OK for this solver.
     */
    bool matrixShapeOK_impl() const;


    /**
     * The LAPACK interface supports the following parameters:
     *
     * <ul>
     *   <li> \c "Equilibrate" (bool) : equilibrate the matrix before
     *   factorization.  The default is \c true .</li>
     * </ul>
     */
    void setParameters_impl(const Teuchos::RCP<Teuchos::ParameterList> & parameterList );

  
    /**
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


    // We keep some `temporary' storage for when we retrive values
    // from the input matrix
    //
    // NOTE:: It seems that the Teuchos::LAPACK wrappers only really
    // work for the OrdinalType=int case, so that's what we'll work
    // with, but ideally we should just be able to give
    // global_ordinal_type.

    /// Stores the values of the nonzero entries
    Teuchos::Array<scalar_type> nzvals_;
    /// Stores the column indices of the nonzero entries
    // Teuchos::Array<global_ordinal_type> rowind_;
    Teuchos::Array<int> rowind_;
    /// Stores the location in \c Ai_ and Aval_ that starts col j
    // Teuchos::Array<global_size_type> colptr_;
    Teuchos::Array<int> colptr_;
    /// Store for RHS and Solution values
    mutable Teuchos::Array<scalar_type> rhsvals_;

    /// L and U storage
    // Teuchos::SerialDenseMatrix<global_ordinal_type,scalar_type> lu_;
    Teuchos::SerialDenseMatrix<int,scalar_type> lu_;

    /// The serial solver
    // Teuchos::SerialDenseSolver<global_ordinal_type,scalar_type> solver_;
    mutable Teuchos::SerialDenseSolver<int,scalar_type> solver_;

  };                              // End class Lapack

  
  // Specialize the solver_traits struct for Lapack.
  //
  // Specializations of Teuchos::LAPACK only exist for real and
  // complex float and double types.
  //
  // TODO: Reinstate the complex support once the bug in
  // Teuchos::SerialDenseSolver has been fixed
  template <>
  struct solver_traits<Lapack> {
// #ifdef HAVE_TEUCHOS_COMPLEX
//     typedef Meta::make_list4<float,
// 			     double,
// 			     std::complex<float>,
// 			     std::complex<double> >supported_scalars;
// #else
    typedef Meta::make_list2<float, double> supported_scalars;
// #endif
};
  
} // end namespace Amesos2

#endif	// AMESOS2_NEWSOLVER_DECL_HPP
