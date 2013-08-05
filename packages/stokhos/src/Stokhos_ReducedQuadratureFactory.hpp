// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_REDUCED_QUADRATURE_FACTORY_HPP
#define STOKHOS_REDUCED_QUADRATURE_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Stokhos_UserDefinedQuadrature.hpp"

namespace Stokhos {

  /*! 
   * \brief Generate a basis from a given set of PCE expansions that is 
   * orthogonal with respect to the product measure induced by these expansions.
   */
  /*!
   * Given the PCE expansions, first build a non-orthogonal monomial basis.  
   * Orthogonalize this basis using Gram-Schmidt, then build a quadrature rule
   * using the simplex method.
   */
  template <typename ordinal_type, typename value_type>
  class ReducedQuadratureFactory {
  public:

    //! Constructor
    /*!
     * \param params Parameter dictating choice of reduction method
     */
    ReducedQuadratureFactory(const Teuchos::ParameterList& params);

    //! Destructor
    virtual ~ReducedQuadratureFactory() {}

    //! Get reduced quadrature object
    virtual Teuchos::RCP<const Stokhos::UserDefinedQuadrature<ordinal_type, value_type> >
    createReducedQuadrature(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q2,
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
      const Teuchos::Array<value_type>& weights) const;

  protected:

    void reducedQuadrature_Q_Squared(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
      const Teuchos::Array<value_type>& weights,
      Teuchos::RCP< Teuchos::Array<value_type> >& red_weights,
      Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& red_points,
      Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& red_values
      ) const;

    void reducedQuadrature_Q_Squared_CPQR(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
      const Teuchos::Array<value_type>& weights,
      Teuchos::RCP< Teuchos::Array<value_type> >& red_weights,
      Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& red_points,
      Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& red_values
      ) const;

     void reducedQuadrature_Q_Squared_CPQR2(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
      const Teuchos::Array<value_type>& weights,
      Teuchos::RCP< Teuchos::Array<value_type> >& red_weights,
      Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& red_points,
      Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& red_values
      ) const;

    void reducedQuadrature_Q2(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q2,
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
      const Teuchos::Array<value_type>& weights,
      Teuchos::RCP< Teuchos::Array<value_type> >& red_weights,
      Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& red_points,
      Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& red_values
      ) const;

    void reducedQuadrature_Q2_CPQR(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q2,
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
      const Teuchos::Array<value_type>& weights,
      Teuchos::RCP< Teuchos::Array<value_type> >& red_weights,
      Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& red_points,
      Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& red_values
      ) const;

    void underdetermined_solver(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A,
      const Teuchos::SerialDenseVector<ordinal_type, value_type>& b,
      Teuchos::SerialDenseVector<ordinal_type, value_type>& x,
      Teuchos::ETransp transa, Teuchos::EUplo uplo) const;

    void solver_TRSM(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A,
      const Teuchos::SerialDenseVector<ordinal_type, value_type>& b,
      Teuchos::SerialDenseVector<ordinal_type, value_type>& x,
      Teuchos::ETransp transa, Teuchos::EUplo uplo) const;

    void solver_GLPK(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A,
      const Teuchos::SerialDenseVector<ordinal_type, value_type>& b,
      Teuchos::SerialDenseVector<ordinal_type, value_type>& x,
      Teuchos::ETransp transa, Teuchos::EUplo uplo) const;

    void solver_CLP(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A,
      const Teuchos::SerialDenseVector<ordinal_type, value_type>& b,
      Teuchos::SerialDenseVector<ordinal_type, value_type>& x,
      Teuchos::ETransp transa, Teuchos::EUplo uplo) const;

    void solver_CLP_IP(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A,
      const Teuchos::SerialDenseVector<ordinal_type, value_type>& b,
      Teuchos::SerialDenseVector<ordinal_type, value_type>& x,
      Teuchos::ETransp transa, Teuchos::EUplo uplo) const;

    void solver_qpOASES(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A,
      const Teuchos::SerialDenseVector<ordinal_type, value_type>& b,
      Teuchos::SerialDenseVector<ordinal_type, value_type>& x,
      Teuchos::ETransp transa, Teuchos::EUplo uplo) const;

    void solver_CompressedSensing(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A,
      const Teuchos::SerialDenseVector<ordinal_type, value_type>& b,
      Teuchos::SerialDenseVector<ordinal_type, value_type>& x,
      Teuchos::ETransp transa, Teuchos::EUplo uplo) const;

    ordinal_type computeRank(
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& R,
      const value_type tol) const;

    /*!
     * \brief Compute bionomial coefficient (n ; k) = n!/( k! (n-k)! )
     */
    ordinal_type n_choose_k(const ordinal_type& n, const ordinal_type& k) const;

  private:

    // Prohibit copying
    ReducedQuadratureFactory(const ReducedQuadratureFactory&);

    // Prohibit Assignment
    ReducedQuadratureFactory& operator=(const ReducedQuadratureFactory&);
    
  protected:

    //! Algorithm parameters
    mutable Teuchos::ParameterList params;

    //! Reduction method
    std::string reduction_method;

    //! Underdetermined solver method
    std::string solver_method;

    //! Whether to eliminate dependent rows in constraints
    bool eliminate_dependent_rows;

    //! Whether to print a bunch of stuff out
    bool verbose;

    //! Dimension reduction tolerance
    value_type reduction_tol;

    //! Value used in LP-based objective function
    value_type objective_value;

    Teuchos::LAPACK<ordinal_type,value_type> lapack;
    Teuchos::BLAS<ordinal_type,value_type> blas;

  }; // class ReducedQuadratureFactory

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_ReducedQuadratureFactoryImp.hpp"

#endif
