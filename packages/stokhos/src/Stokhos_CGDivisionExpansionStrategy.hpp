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
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_CG_DIVISION_EXPANSION_STRATEGY_HPP
#define STOKHOS_CG_DIVISION_EXPANSION_STRATEGY_HPP

#include "Stokhos_DivisionExpansionStrategy.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "Stokhos_DiagPreconditioner.hpp"
#include "Stokhos_JacobiPreconditioner.hpp"
#include "Stokhos_GSPreconditioner.hpp"
#include "Stokhos_SchurPreconditioner.hpp"
#include "Stokhos_InversePreconditioner.hpp"
#include "Stokhos_BlockPreconditioner.hpp"

#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"

#include <iostream>
#include <fstream>

namespace Stokhos {

  //! Strategy interface for computing PCE of a/b using only b[0]
  /*!
   * Such a strategy is only useful when the division occurs in a preconditioner
   */
  template <typename ordinal_type, typename value_type, typename node_type> 
  class CGDivisionExpansionStrategy :
    public DivisionExpansionStrategy<ordinal_type,value_type,node_type> {
  public:

    //! Constructor
    CGDivisionExpansionStrategy(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis_,
      const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk_, const int prec_iter_, const double tol_, const int PrecNum_, const int max_it_, const int linear_, const int diag_, const int equil_);

    //! Destructor
    virtual ~CGDivisionExpansionStrategy() {}
 
    // Division operation:  c = \alpha*(a/b) + beta*c
    virtual void divide(
      Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c,
      const value_type& alpha,
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b,
      const value_type& beta);

  private:

    // Prohibit copying
    CGDivisionExpansionStrategy(
      const CGDivisionExpansionStrategy&);

    // Prohibit Assignment
    CGDivisionExpansionStrategy& operator=(
      const CGDivisionExpansionStrategy& b);

      int CG(const Teuchos::SerialDenseMatrix<int, double> &  A, Teuchos::SerialDenseMatrix<int,double> &  X, const Teuchos::SerialDenseMatrix<int,double> &   B, int max_iter, double tolerance, int prec_iter, int order, int dim, int PrecNum, const Teuchos::SerialDenseMatrix<int, double> &  M, int diag);

  protected:

    //! Basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> > basis;

    //! Short-hand for Cijk
    typedef Stokhos::Sparse3Tensor<ordinal_type, value_type> Cijk_type;

    //! Triple product
    Teuchos::RCP<const Cijk_type> Cijk;

    //! Dense matrices for linear system
    Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type,value_type> > A, X, B, M;
    
    //! Tolerance for CG
    int prec_iter;

    double tol;
 
    int PrecNum;

    int max_it;

    int linear;

    int diag;

    int equil;
       
  }; // class CGDivisionExpansionStrategy

} // namespace Stokhos

template <typename ordinal_type, typename value_type, typename node_type> 
Stokhos::CGDivisionExpansionStrategy<ordinal_type,value_type,node_type>::
CGDivisionExpansionStrategy(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis_,
  const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk_,
  const int prec_iter_, const double tol_, const int PrecNum_, const int max_it_, const int linear_, const int diag_, const int equil_): 
  basis(basis_),
  Cijk(Cijk_),
  prec_iter(prec_iter_),
  tol(tol_),
  PrecNum(PrecNum_),
  max_it(max_it_),
  linear(linear_),
  diag(diag_),
  equil(equil_)
  
{
  ordinal_type sz = basis->size();
  A = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type,value_type>(
		     sz, sz));
  B = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type,value_type>(
		     sz, 1));
  X = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type,value_type>(
		     sz, 1));
   M = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type,value_type>(
                     sz, sz));

}


template <typename ordinal_type, typename value_type, typename node_type> 
void
Stokhos::CGDivisionExpansionStrategy<ordinal_type,value_type,node_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c,
       const value_type& alpha,
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b,
       const value_type& beta)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::CGDivisionStrategy::divide()");
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
	     i_it  != Cijk->i_end(j_it); ++i_it) {       
	  i = index(i_it);
	  cijk = value(i_it);
	  (*A)(i,j) += cijk*cb[k];
	}
      }
    }
// Compute B
    B->putScalar(0.0);
    for (ordinal_type i=0; i<pa; i++)
        (*B)(i,0) = ca[i]*basis->norm_squared(i);
Teuchos::SerialDenseMatrix<int,double> D(sz, 1);
//Equilibrate the linear system
if (equil == 1){
	//Create diag mtx of max row entries
	for (int i=0; i<sz; i++){
		Teuchos::SerialDenseMatrix<int, double> r(Teuchos::View, *A, 1, sz, i, 0);
		D(i,0)=sqrt(r.normOne());
	}


	//Compute inv(D)*A*inv(D)
	for (int i=0; i<sz; i++){
		for (int j=0; j<sz; j++){
			(*A)(i,j)=(*A)(i,j)/(D(i,0)*D(j,0));
	}}

	//Scale b by inv(D)
	for (int i=0; i<sz; i++){
		(*B)(i,0)=(*B)(i,0)/D(i,0);
      	}

}

if (linear == 1){
//Compute M, the linear matrix to be used in the preconditioner

 pb = basis->dimension()+1;

 M->putScalar(0.0);
    if (pb < Cijk->num_k())
      k_end = Cijk->find_k(pb);
    for (typename Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
      k = index(k_it);
      for ( typename Cijk_type::kj_iterator j_it = Cijk->j_begin(k_it);
           j_it != Cijk->j_end(k_it); ++j_it) {
        j = index(j_it);
        for ( typename Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it);
             i_it  != Cijk->i_end(j_it); ++i_it) {
          i = index(i_it);
          cijk = value(i_it);
          (*M)(i,j) += cijk*cb[k];
        }
      }
    }

 //Scale M
 if (equil == 1){
     //Compute inv(D)*M*inv(D)
	for (int i=0; i<sz; i++){
        	for (int j=0; j<sz; j++){
            		(*M)(i,j)=(*M)(i,j)/(D(i,0)*D(j,0));
            }
 	}
  }
         CG(*A,*X,*B, max_it, tol, prec_iter, basis->order(), basis->dimension(), PrecNum, *M, diag);
}

else{
  
	CG(*A,*X,*B, max_it, tol, prec_iter, basis->order(), basis->dimension(), PrecNum, *A, diag);
}
  
if (equil == 1 ) {
	//Rescale X 
	for (int i=0; i<sz; i++){
        	(*X)(i,0)=(*X)(i,0)/D(i,0);
        }
}
   
// Compute c
    for (ordinal_type i=0; i<pc; i++)
      cc[i] = alpha*(*X)(i,0) + beta*cc[i];
  }
  else {
    for (ordinal_type i=0; i<pc; i++)
      cc[i] = alpha*ca[i]/cb[0] + beta*cc[i];
  }
}
 

template <typename ordinal_type, typename value_type, typename node_type>
int
Stokhos::CGDivisionExpansionStrategy<ordinal_type,value_type,node_type>::
CG(const Teuchos::SerialDenseMatrix<int, double> &  A, Teuchos::SerialDenseMatrix<int,double> &  X, const Teuchos::SerialDenseMatrix<int,double> &   B, int max_iter, double tolerance, int prec_iter, int order , int m, int PrecNum, const Teuchos::SerialDenseMatrix<int, double> &  M, int diag)

{
  int n = A.numRows();
  int k=0;
  double resid;
  Teuchos::SerialDenseMatrix<int, double> Ax(n,1);
  Ax.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0, A, X, 0.0);
  Teuchos::SerialDenseMatrix<int, double> r(B);
  r-=Ax;
  resid=r.normFrobenius();
  Teuchos::SerialDenseMatrix<int, double> p(r);
  Teuchos::SerialDenseMatrix<int, double> rho(1,1);
  Teuchos::SerialDenseMatrix<int, double> oldrho(1,1);
  Teuchos::SerialDenseMatrix<int, double> pAp(1,1);
  Teuchos::SerialDenseMatrix<int, double> Ap(n,1);
  double b;
  double a;
  while (resid > tolerance && k < max_iter){
    Teuchos::SerialDenseMatrix<int, double> z(r);
    //Solve Mz=r
    if (PrecNum != 0){
        if (PrecNum == 1){
               Stokhos::DiagPreconditioner precond(M);
               precond.ApplyInverse(r,z,prec_iter);
          }
  	else if (PrecNum == 2){
        	Stokhos::JacobiPreconditioner precond(M);
        	precond.ApplyInverse(r,z,2);
          }
  	else if (PrecNum == 3){
        	Stokhos::GSPreconditioner precond(M,0);
        	precond.ApplyInverse(r,z,1);
          }
        else if (PrecNum == 4){
      		Stokhos::SchurPreconditioner precond(M, order, m, diag);
                precond.ApplyInverse(r,z,prec_iter);            
          }
    }
    rho.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0, r, z, 0.0);
  

    if (k==0){
       p.assign(z);
       rho.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, r, z, 0.0);  
       }
    else {
       b=rho(0,0)/oldrho(0,0);
       p.scale(b);
       p+=z; 
    }
    Ap.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0, A, p, 0.0);
    pAp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0, p, Ap, 0.0);
    a=rho(0,0)/pAp(0,0);
    Teuchos::SerialDenseMatrix<int, double> scalep(p);
    scalep.scale(a);
    X+=scalep;
    Ap.scale(a);
    r-=Ap;
    oldrho.assign(rho);
    resid=r.normFrobenius();
    k++;
    }                      
 
  std::cout << "iteration count  " << k << std::endl;
  return 0; 
}

 #endif // STOKHOS_DIVISION_EXPANSION_STRATEGY_HPP
