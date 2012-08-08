// Id$ 
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
// @HEADER  Implements GMRES for scalar divison

#include "Stokhos_DivisionExpansionStrategy.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"

#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Stokhos.hpp"
#include "Stokhos_Sacado.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"


int main(int argc, char *argv[]) 
{
  Teuchos::SerialDenseMatrix<int, double> A(4,4);
  Teuchos::SerialDenseMatrix<int, double> B(4,1);
  Teuchos::SerialDenseMatrix<int, double> X(4,1);
  My_Matrix.random();
  My_Vector.random();
  X.putScalar(0.0);
  
  int maxit=15;
  double tolerance=1e-6;
   
  gmres(A,B,X,maxit,tolerance); 
  
  return 0;
}
//Gmres 
void gmres(const Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type,value_type>&  A, const Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type,value_type>&  X,const Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type,value_type>&  B, int max_iter, double tolerance)
{
  
  int n; 
  int k;
  double resid;
  k=1;

  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Ax(n,1);
  Ax.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0, A, X, 0.0);
  
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> r0(B);
  r0-=Ax;
  
  resid=B.normFrobenius();
  
  //define vector v=r/norm(r) where r=b-Ax
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> v(n,1);
  v=r0.scale(1.0/resid);
  

  Teuchos::SerialDenseMatrix<ordinal_type, value_type> h(1,1);

  //Matrix of orthog basis vectors V
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> V(n,n);
  
  //vector vj that is jth col of V
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> vj(n,1);
  //Set v to be 1st col of V
  for (int i=0; i<n; i++){
          V(i,0)=v(i,0);
         }
  
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> bb(1,1);
  bb(1,1)=resid;


  //Matrices needed in least squares problem 
  
  
  
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Gk();  //current step rotation mtx
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Htemp(); //G*H 
  // Teuchos::SerialDenseMatrix<ordinal_type, value_type> R(n,1);   //upper triangular mtx R of QR factorization
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> c();  

  Teuchos::SerialDenseMatrix<ordinal_type, value_type> w(n,1);
  while (resid > tolerance && k < max_iter){
    h.reshape(k+1,k);
    //Arnoldi iteration(via G-S)
    w.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0, A, v, 0.0);
    
    for (int j=0; j<k; j++){
      for (int i=0; i<n; i++){
         //set vector vj to be jth col of V
         vj(i,0)=V(i,j);
         }
       //Calculate inner product
       Teuchos::SerialDenseMatrix<ordinal_type, value_type> ip(1,1);
       ip.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, vj, w, 0.0);
       double* ipvalue = ip.values();      
       h(j,k-1)= ipvalue;
       //scale vj by h(j,k-1)
       vj.scale(ipvalue);
       w-=vj;
    }
    h(k;k-1)=w.normFrobenius();     
    v=w.scale(w.normFrobenius());   
    //add column v to V
    for (int i=0; i<n; i++){
          V(i,k)=v(i,0);
         }
    //Solve upper hessenberg least squares problem via Givens rotations
    if (k==1){
      Teuchos::SerialDenseMatrix<ordinal_type, value_type> Htemp(h);
      Teuchos::SerialDenseMatrix<ordinal_type, value_type> G(2,2);   //product of previous Gks (rotation mtxs)
      G(0,0)=1.0;  // 2x2 identity mtx 
      G(1,1)=1.0;
      G(0,1)=0.0;
      G(1,0)=0.0;
    }
    else {
      //set Htemp=G*H where trans(Htemp)=Q approx
      G.reshape(k+1,k+1);
      G(k,k)=1;  //need to fill non-diag entries of G with 0???
      Htemp.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, G, h, 0.0);
    }
    //form Givens rotation matrix Gk to eliminate the h(k,k-1) entry ???

    Gk.shape(k+1,k+1);
    for (int l=0; l<k; k++){
       Gk(l,l)=1;
       }
    Gk(k-1,k-1)=(Htemp[k-1,k-1])/sqrt(Htemp[k-1,k-1]*Htemp[k-1,k-1]+Htemp[k,k-1]*Htemp[k,k-1]); 
    Gk(k-1,k)=(Htemp[k,k-1])/sqrt(Htemp[k-1,k-1]*Htemp[k-1,k-1]+Htemp[k,k-1]*Htemp[k,k-1]);
    Gk(k,k-1)=(-Htemp[k,k-1])/sqrt(Htemp[k-1,k-1]*Htemp[k-1,k-1]+Htemp[k,k-1]*Htemp[k,k-1]);
    Gk(k,k)=(Htemp[k-1,k-1])/sqrt(Htemp[k-1,k-1]*Htemp[k-1,k-1]+Htemp[k,k-1]*Htemp[k,k-1]);
    //Combine with previous Givens rotations
    G.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Gk, G, 0.0);
    R.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, G, H, 0.0);
    bb=resize(k+1,1);
    c.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, G, bb, 0.0);
    //Determine residual    
    resid = fabs(c[0,k]);
    
          
  } 
   //Extract upper triangular square matrix
   R.reshape(R.numRows()-1, R.numCols());
   c.reshape(R.numRows(),1);
   //Solve linear system
   int info;
   //TRTRS (const char UPLO, const char TRANS, const char DIAG, const OrdinalType n, const OrdinalType nrhs, const ScalarType *A, const OrdinalType lda, ScalarType *B, const OrdinalType ldb, OrdinalType *info)
   Teuchos::LAPACK<ordinal_type, scalar_type>::TRTRS('U', Teuchos::NO_TRANS, 'N', R.numRows()-1, 1, R.values(), R.stride(), c.values(), c.stride(),&info); 

return 0;
}






/*
  template <typename ordinal_type, typename value_type, typename node_type> 
  class GMRESDivisionExpansionStrategy :
    public DivisionExpansionStrategy<ordinal_type,value_type,node_type> {
  public:

    //! Constructor
    GMRESDivisionExpansionStrategy(
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis_,
      const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk_);

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

  protected:

    //! Basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> > basis;

    //! Short-hand for Cijk
    typedef Stokhos::Sparse3Tensor<ordinal_type, value_type> Cijk_type;

    //! Triple product
    Teuchos::RCP<const Cijk_type> Cijk;

    //! Dense matrices for linear system
    Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type,value_type> > A, X, B;

    //! Serial dense solver
    Teuchos::SerialDenseSolver<ordinal_type,value_type> solver;
    //Create GMRES alg=solver=gmres
  }; // class DenseDirectDivisionExpansionStrategy

} // namespace Stokhos

template <typename ordinal_type, typename value_type, typename node_type> 
Stokhos::GMRESDivisionExpansionStrategy<ordinal_type,value_type,node_type>::
GMRESDivisionExpansionStrategy(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis_,
  const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk_) :
  basis(basis_),
  Cijk(Cijk_),
  gmres()
{
  ordinal_type sz = basis->size();
  A = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type,value_type>(
		     sz, sz));
  B = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type,value_type>(
		     sz, 1));
  X = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type,value_type>(
		     sz, 1));
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
    A->putScalar(0.0);//sets all values of A to 0.0
    typename Cijk_type::k_iterator k_begin = Cijk->k_begin();
    typename Cijk_type::k_iterator k_end = Cijk->k_end();
    if (pb < Cijk->num_k())
      k_end = Cijk->find_k(pb);
    value_type tmp, cijk;
    ordinal_type i,j,k;
    for (typename Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
      k = index(k_it);
      tmp = value_type(0.0);
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
   if (pb>1)  {
     //Compute A
     A->putScalar(0.0);
     typename Cijk_type::k_iterator k_begin = Cijk->k_begin();
     typename Cijk_type::k_iterator k_end = Cijk->k_end();
     if (pb < Cijk->num_k())
     k_end = Cijk->find_k(pb);
     value_type tmp, cijk;
     ordinal_type i,j,k;
     for (typename Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
       k = index(k_it);
       tmp = value_type(0.0);
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
    solver.solve();*/



#endif // STOKHOS_DIVISION_EXPANSION_STRATEGY_HPP  
