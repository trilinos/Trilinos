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

#ifndef STOKHOS_GMRES_DIVISION_EXPANSION_STRATEGY_HPP
#define STOKHOS_GMRES_DIVISION_EXPANSION_STRATEGY_HPP

#include "Stokhos_DivisionExpansionStrategy.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "Stokhos_DiagPreconditioner.hpp"
#include "Stokhos_JacobiPreconditioner.hpp"
#include "Stokhos_GSPreconditioner.hpp"
#include "Stokhos_SchurPreconditioner.hpp"
#include "Stokhos_BlockPreconditioner.hpp"

#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"

namespace Stokhos {

  //! Strategy interface for computing PCE of a/b using only b[0]
  /*!
   * Such a strategy is only useful when the division occurs in a preconditioner
   */
  template <typename ordinal_type, typename value_type, typename node_type> 
  class GMRESDivisionExpansionStrategy :
    public DivisionExpansionStrategy<ordinal_type,value_type,node_type> {
  public:

    //! Constructor
    GMRESDivisionExpansionStrategy(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis_,
      const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk_, const int prec_iter_, const double tol_, const int PrecNum_, const int max_it_, const int linear_, const int diag_, const int equil_);

    //! Destructor
    virtual ~GMRESDivisionExpansionStrategy() {}
 
    // Division operation:  c = \alpha*(a/b) + beta*c
    virtual void divide(
      Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c,
      const value_type& alpha,
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b,
      const value_type& beta);

  private:

    // Prohibit copying
    GMRESDivisionExpansionStrategy(
      const GMRESDivisionExpansionStrategy&);

    // Prohibit Assignment
    GMRESDivisionExpansionStrategy& operator=(
      const GMRESDivisionExpansionStrategy& b);

      int GMRES(const Teuchos::SerialDenseMatrix<int, double> &  A, Teuchos::SerialDenseMatrix<int,double> &  X, const Teuchos::SerialDenseMatrix<int,double> &   B, int max_iter, double tolerance, int prec_iter, int order, int dim,int PrecNum, const Teuchos::SerialDenseMatrix<int, double>& M, int diag);

  protected:

    //! Basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> > basis;

    //! Short-hand for Cijk
    typedef Stokhos::Sparse3Tensor<ordinal_type, value_type> Cijk_type;

    //! Triple product
    Teuchos::RCP<const Cijk_type> Cijk;

    //! Dense matrices for linear system
    Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type,value_type> > A, X, B, M;
    
    //! Tolerance for GMRES

    int prec_iter;

    double tol;

    int PrecNum;

    int max_it;

    int linear;

    int diag;

    int equil;
   
    
  }; // class GMRESDivisionExpansionStrategy

} // namespace Stokhos

template <typename ordinal_type, typename value_type, typename node_type> 
Stokhos::GMRESDivisionExpansionStrategy<ordinal_type,value_type,node_type>::
GMRESDivisionExpansionStrategy(
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
Stokhos::GMRESDivisionExpansionStrategy<ordinal_type,value_type,node_type>::
divide(Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c,
       const value_type& alpha,
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b,
       const value_type& beta)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::GMRESDivisionStrategy::divide()");
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

   Teuchos::SerialDenseMatrix<ordinal_type,value_type> D(sz, 1);
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


   	// Compute X = A^{-1}*B  
    	GMRES(*A,*X,*B, max_it, tol, prec_iter, basis->order(), basis->dimension(), PrecNum, *M, diag);
  }

  else{
	GMRES(*A,*X,*B, max_it, tol, prec_iter, basis->order(), basis->dimension(), PrecNum, *A, diag);
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
Stokhos::GMRESDivisionExpansionStrategy<ordinal_type,value_type,node_type>::
GMRES(const Teuchos::SerialDenseMatrix<int, double> &  A, Teuchos::SerialDenseMatrix<int,double> &  X, const Teuchos::SerialDenseMatrix<int,double> &   B, int max_iter, double tolerance, int prec_iter, int order, int dim, int PrecNum, const Teuchos::SerialDenseMatrix<int, double> & M, int diag)
{
 int n = A.numRows();
 int k = 1;
 double resid;
 Teuchos::SerialDenseMatrix<int, double> P(n,n);
 Teuchos::SerialDenseMatrix<int, double> Ax(n,1);
 Ax.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0, A, X, 0.0);
 Teuchos::SerialDenseMatrix<int, double> r0(B);
 r0-=Ax;
 resid=r0.normFrobenius();
 //define vector v=r/norm(r) where r=b-Ax
 Teuchos::SerialDenseMatrix<int, double> v(n,1);
 r0.scale(1/resid);
 Teuchos::SerialDenseMatrix<int, double> h(1,1);
 //Matrix of orthog basis vectors V
 Teuchos::SerialDenseMatrix<int, double> V(n,1);
 //Set v=r0/norm(r0) to be 1st col of V
 for (int i=0; i<n; i++){
    V(i,0)=r0(i,0);
  }
  //right hand side
  Teuchos::SerialDenseMatrix<int, double> bb(1,1);
  bb(0,0)=resid;
  Teuchos::SerialDenseMatrix<int, double> w(n,1);
  Teuchos::SerialDenseMatrix<int, double> c;
  Teuchos::SerialDenseMatrix<int, double> s;
  while (resid > tolerance && k < max_iter){
    h.reshape(k+1,k);
    //Arnoldi iteration(Gram-Schmidt )
    V.reshape(n,k+1);
    //set vk to be kth col of V
    Teuchos::SerialDenseMatrix<int, double> vk(Teuchos::Copy, V, n,1,0,k-1);
    //Preconditioning step: solve Mz=vk
    Teuchos::SerialDenseMatrix<int, double> z(vk);
        if (PrecNum == 1){
               Stokhos::DiagPreconditioner precond(M);
               precond.ApplyInverse(vk,z,prec_iter);
        }
        else if (PrecNum == 2){
	       	Stokhos::JacobiPreconditioner precond(M);
        	precond.ApplyInverse(vk,z,2);
        }
 	else if (PrecNum == 3){
        	Stokhos::GSPreconditioner precond(M,1);
        	precond.ApplyInverse(vk,z,1);
        }
  	else if (PrecNum == 4){
        	Stokhos::SchurPreconditioner precond(M, order, dim, diag);
	        precond.ApplyInverse(vk,z,prec_iter);
	}

    w.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1, A, z, 0.0);
    Teuchos::SerialDenseMatrix<int, double> vi(n,1);
    Teuchos::SerialDenseMatrix<int, double> ip(1,1);
    for (int i=0; i<k; i++){
       //set vi to be ith col of V
       Teuchos::SerialDenseMatrix<int, double> vi(Teuchos::Copy, V, n,1,0,i);
       //Calculate inner product
       ip.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, vi, w, 0.0);
       h(i,k-1)= ip(0,0);
       //scale vi by h(i,k-1)
       vi.scale(ip(0,0));
       w-=vi;
     }
     h(k,k-1)=w.normFrobenius(); 
     w.scale(1.0/h(k,k-1));
     //add column vk+1=w to V
     for (int i=0; i<n; i++){
        V(i,k)=w(i,0);
        }
     //Solve upper hessenberg least squares problem via Givens rotations
     //Compute previous Givens rotations
     for (int i=0; i<k-1; i++){
        double q=c(i,0)*h(i,k-1)+s(i,0)*h(i+1,k-1);
        h(i+1,k-1)=-1*s(i,0)*h(i,k-1)+c(i,0)*h(i+1,k-1);
        h(i,k-1)=q;

     }
     //Compute next Givens rotations
     c.reshape(k,1);
     s.reshape(k,1);
     bb.reshape(k+1,1);
     double l = sqrt(h(k-1,k-1)*h(k-1,k-1)+h(k,k-1)*h(k,k-1));
     c(k-1,0)=h(k-1,k-1)/l;
     s(k-1,0)=h(k,k-1)/l;
    
     // Givens rotation on h and bb
     h(k-1,k-1)=l;
     h(k,k-1)=0;
    
     bb(k,0)=-s(k-1,0)*bb(k-1,0);
     bb(k-1,0)=c(k-1,0)*bb(k-1,0);
    
     //Determine residual    
     resid = fabs(bb(k,0));
     k++;
   }
   //Extract upper triangular square matrix
   bb.reshape(h.numRows()-1 ,1);
   //Solve linear system
   int info;
   Teuchos::LAPACK<int, double> lapack;
   lapack.TRTRS('U', 'N', 'N', h.numRows()-1, 1, h.values(), h.stride(), bb.values(), bb.stride(),&info);
   Teuchos::SerialDenseMatrix<int, double> ans(X);
   V.reshape(n,k-1);
   ans.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, V, bb, 0.0);
          if (PrecNum == 1){
               Stokhos::DiagPreconditioner precond(M);
               precond.ApplyInverse(ans,ans,prec_iter);
          }
          else if (PrecNum == 2){
                Stokhos::JacobiPreconditioner precond(M);
                precond.ApplyInverse(ans,ans,2);
           }
          else if (PrecNum == 3){
                Stokhos::GSPreconditioner precond(M,1);
                precond.ApplyInverse(ans,ans,1);
          }
          else if (PrecNum == 4){
                Stokhos::SchurPreconditioner precond(M, order, dim, diag);
                precond.ApplyInverse(ans,ans,prec_iter);}
   X+=ans;

   std::cout << "iteration count=  " << k-1 << std::endl;        


return 0;
}



 #endif // STOKHOS_DIVISION_EXPANSION_STRATEGY_HPP
