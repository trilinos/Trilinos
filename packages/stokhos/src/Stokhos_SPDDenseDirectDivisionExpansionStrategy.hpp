// $Id$ 
// $Source$ 
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

#ifndef STOKHOS_SPD_DENSE_DIRECT_DIVISION_EXPANSION_STRATEGY_HPP
#define STOKHOS_SPD_DENSE_DIRECT_DIVISION_EXPANSION_STRATEGY_HPP

#include "Stokhos_DivisionExpansionStrategy.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"

#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"
#include "Teuchos_SerialSpdDenseSolver.hpp"

namespace Stokhos {

  //! Strategy interface for computing PCE of a/b using only b[0]
  /*!
   * Such a strategy is only useful when the division occurs in a preconditioner
   */
  template <typename ordinal_type, typename value_type, typename node_type> 
  class SPDDenseDirectDivisionExpansionStrategy :
    public DivisionExpansionStrategy<ordinal_type,value_type,node_type> {
  public:

    //! Constructor
    SPDDenseDirectDivisionExpansionStrategy(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis_,
      const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk_);

    //! Destructor
    virtual ~SPDDenseDirectDivisionExpansionStrategy() {}
 
    // Division operation:  c = \alpha*(a/b) + beta*c
    virtual void divide(
      Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c,
      const value_type& alpha,
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b,
      const value_type& beta);

  private:

    // Prohibit copying
    SPDDenseDirectDivisionExpansionStrategy(
      const SPDDenseDirectDivisionExpansionStrategy&);

    // Prohibit Assignment
    SPDDenseDirectDivisionExpansionStrategy& operator=(
      const SPDDenseDirectDivisionExpansionStrategy& b);

  protected:

    //! Basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> > basis;

    //! Short-hand for Cijk
    typedef Stokhos::Sparse3Tensor<ordinal_type, value_type> Cijk_type;

    //! Triple product
    Teuchos::RCP<const Cijk_type> Cijk;

    //! Dense matrices for linear system
    Teuchos::RCP< Teuchos::SerialSymDenseMatrix<ordinal_type,value_type> > A;

    Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type,value_type> > X, B;

    //! Serial dense solver
    Teuchos::SerialSpdDenseSolver<ordinal_type,value_type> solver;
    
  }; // class SPDDenseDirectDivisionExpansionStrategy

} // namespace Stokhos

template <typename ordinal_type, typename value_type, typename node_type> 
Stokhos::SPDDenseDirectDivisionExpansionStrategy<ordinal_type,value_type,node_type>::
SPDDenseDirectDivisionExpansionStrategy(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis_,
  const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk_) :
  basis(basis_),
  Cijk(Cijk_),
  solver()
{
  ordinal_type sz = basis->size();
  A = Teuchos::rcp(new Teuchos::SerialSymDenseMatrix<ordinal_type,value_type>(
		     sz, sz));
  B = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type,value_type>(
		     sz, 1));
  X = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type,value_type>(
		     sz, 1));
}

template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::SPDDenseDirectDivisionExpansionStrategy<ordinal_type,value_type,node_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c,
       const value_type& alpha,
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b,
       const value_type& beta)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::SPDDenseDirectDivisionStrategy::divide()");
#endif

  ordinal_type sz = basis->size();
  ordinal_type pa = a.size();
  ordinal_type pb = b.size();
  ordinal_type pc;
  if (pb > 1)
    pc = sz;
  else
    pc = pa;
  if (c.size() != pc)
    c.resize(pc);

  const value_type* ca = a.coeff();
  const value_type* cb = b.coeff();
  value_type* cc = c.coeff();

  if (pb > 1) {
    // Compute A
    A->putScalar(0.0);
    typename Cijk_type::k_iterator k_begin = Cijk->k_begin();
    typename Cijk_type::k_iterator k_end = Cijk->k_end();
    if (pb < Cijk->num_k())
      k_end = Cijk->find_k(pb);
    value_type cijk;
    ordinal_type i,j,k;
    for (typename Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
      k = index(k_it);
      for (typename Cijk_type::kj_iterator j_it = Cijk->j_begin(k_it); 
	   j_it != Cijk->j_end(k_it); ++j_it) {
	j = index(j_it);
	for (typename Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it);
	     i_it != Cijk->i_end(j_it); ++i_it) {
	  i = index(i_it);
	  cijk = value(i_it);
	  (*A)(i,j) += cijk*cb[k];
	}
      }
    }
    A->setUpper();

    // Compute B
    B->putScalar(0.0);
    for (ordinal_type i=0; i<pa; i++)
      (*B)(i,0) = ca[i]*basis->norm_squared(i);

    // Setup solver
    solver.setMatrix(A);
    solver.setVectors(X, B);
    if (solver.shouldEquilibrate()) {
      solver.factorWithEquilibration(true);
      solver.equilibrateMatrix();
    }

    // Compute X = A^{-1}*B
    solver.solve();

    // value_type kappa;
    // solver.reciprocalConditionEstimate(kappa);
    // std::cout << "kappa = " << 1.0/kappa << std::endl;

    // Compute c
    for (ordinal_type i=0; i<pc; i++)
      cc[i] = alpha*(*X)(i,0) + beta*cc[i];
  }
  else {
    for (ordinal_type i=0; i<pc; i++)
      cc[i] = alpha*ca[i]/cb[0] + beta*cc[i];
  }
}

#endif // STOKHOS_DIVISION_EXPANSION_STRATEGY_HPP
