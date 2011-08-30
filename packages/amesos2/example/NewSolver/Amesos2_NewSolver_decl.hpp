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
   \file   Amesos2_NewSolver_decl.hpp
   \author John Doe <jd@sandia.gov>
   \date   Thu Jul  8 16:08:10 2010
   
   \brief  A template class that does nothing useful besides show developers
           what, in general, needs to be done to add a new solver interface to
           the Amesos2 collection.
*/


#ifndef AMESOS2_NEWSOLVER_DECL_HPP
#define AMESOS2_NEWSOLVER_DECL_HPP

#include "Amesos2_SolverTraits.hpp"
#include "Amesos2_SolverCore.hpp"
#include "Amesos2_NewSolver_FunctionMap.hpp"
#include "Amesos2_NewSolver_TypeMap.hpp"


namespace Amesos2 {


  /** \brief Amesos2 interface to the NewSolver package.
   *
   * Describe here anything you would like about the NewSolver solver:
   * what sorts of problems it is good at solving, other documentation
   * about the solver itself, what scalar and ordinal types are
   * supported, etc.
   *
   * \ingroup amesos2_solver_interfaces 
   */
  template <class Matrix,
	    class Vector>
  class NewSolver : public SolverCore<Amesos2::NewSolver, Matrix, Vector> 
  {
    friend class SolverCore<Amesos2::NewSolver,Matrix,Vector>; // Give our base access
    // to our private
    // implementation funcs
  public:

    /// The name of this solver interface
    static const char* name;	// declaration. Initialization outside.

    typedef NewSolver<Matrix,Vector>                                     type;
    typedef SolverCore<Amesos2::NewSolver,Matrix,Vector>           super_type;

    // Since typedef's are not inheritted, go grab them
    typedef typename super_type::scalar_type                      scalar_type;
    typedef typename super_type::local_ordinal_type        local_ordinal_type;
    typedef typename super_type::global_ordinal_type      global_ordinal_type;
    typedef typename super_type::global_size_type            global_size_type;

    typedef TypeMap<Amesos2::NewSolver,scalar_type>                  type_map;

    typedef typename type_map::type                        solver_scalar_type;
    typedef typename type_map::magnitude_type           solver_magnitude_type;

    typedef FunctionMap<Amesos2::NewSolver,slu_type>             function_map;


    /// \name Constructor/Destructor methods
    //@{ 

    /**
     * \brief Initialize from Teuchos::RCP.
     *
     * \warning Should not be called directly!  Use instead
     * Amesos2::create() to initialize a NewSolver interface.
     */
    NewSolver(Teuchos::RCP<const Matrix> A,
	      Teuchos::RCP<Vector>       X,
	      Teuchos::RCP<const Vector> B);
  

    /// Destructor
    ~NewSolver( );

    //@}

  private:
  
    /**
     * \brief Performs pre-ordering on the matrix to increase efficiency.
     */
    int preOrdering_impl();


    /**
     * \brief Perform symbolic factorization of the matrix using NewSolver.
     *
     * Called first in the sequence before numericFactorization.
     *
     * \throw std::runtime_error NewSolver is not able to factor the matrix.
     */
    int symbolicFactorization_impl();


    /**
     * \brief NewSolver specific numeric factorization
     * 
     * \throw std::runtime_error NewSolver is not able to factor the matrix
     */
    int numericFactorization_impl();


    /**
     * \brief NewSolver specific solve.
     *
     * Uses the symbolic and numeric factorizations, along with the RHS vector
     * \c B to solve the sparse system of equations.
     *
     * The solution of the system is placed in X.
     *
     * \throw std::runtime_error NewSolver is not able to solve the system.
     */
    int solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
		   const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const;


    /**
     * \brief Determines whether the shape of the matrix is OK for this solver.
     */
    bool matrixShapeOK_impl() const;


    /**
     * Describe here the options that this interface supports.  The body
     * of this comment should be copied into the amesos2_parameters
     * section in Amesos2.hpp.
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



    /* Declare private variables necessary for interaction with the
     * NewSolver TPL.
     *
     * For example, the following Arrays are persisting storage arrays
     * for A, X, and B that can be used with solvers expecting a
     * compressed-row representation of the matrix A.
     */
   
    /// Stores the values of the nonzero entries for NewSolver
    Teuchos::Array<solver_scalar_type> nzvals_;
    /// Stores the location in \c Ai_ and Aval_ that starts row j
    Teuchos::Array<int> colind_;
    /// Stores the row indices of the nonzero entries
    Teuchos::Array<int> rowptr_;
    /// Persisting 1D store for X
    Teuchos::Array<solver_scalar_type> xvals_;  int ldx_;
    /// Persisting 1D store for B
    Teuchos::Array<solver_scalar_type> bvals_;  int ldb_;

  };                              // End class NewSolver

  
  // Specialize the solver_traits struct for NewSolver.
  //
  // The `supported_scalars' member type will be used by
  // Amesos2::create to control instantiation of the NewSolver
  // interface.  That is, create() will only instantiate the NewSolver
  // interface for types that it supports, and will throw a runtime
  // error for any other types.
  //
  // For example, if NewSolver supports float, double, complex float,
  // and complex double scalar types.
  template <>
  struct solver_traits<NewSolver> {
    typedef Meta::make_list4<float,
			     double,
			     std::complex<float>,
			     std::complex<double> > supported_scalars;
};
  
} // end namespace Amesos2

#endif	// AMESOS2_NEWSOLVER_DECL_HPP
