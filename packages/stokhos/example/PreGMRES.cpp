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
#include "Teuchos_TestForException.hpp"


//Mean-Based Preconditioned GMRES  
int pregmres(const Teuchos::SerialDenseMatrix<int, double> &  A, const Teuchos::SerialDenseMatrix<int,double> &  X,const Teuchos::SerialDenseMatrix<int,double> &   B, int max_iter, double tolerance)
{
  int n; 
  int k;
  double resid;
  k=1;
  n=A.numRows();
  std::cout << A << std::endl;
  Teuchos::SerialDenseMatrix<int, double> D(n,1);

  //Get diagonal entries of A 
  for (int i=0; i<n; i++){
    D(i,0)=A(i,i);
  }
  

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
    
    std::cout << "k = " << k << std::endl;
    h.reshape(k+1,k);
    //Arnoldi iteration(Gram-Schmidt )
    V.reshape(n,k+1);    
    //set vk to be kth col of V
    Teuchos::SerialDenseMatrix<int, double> vk(Teuchos::Copy, V, n,1,0,k-1);
    //Preconditioning step w=AMj(-1)vj
    w.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1/D(k-1,0), A, vk, 0.0);  

   
    
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
    
    w.scale(1.0/w.normFrobenius());   
    //add column vk+1=w to V
    for (int i=0; i<n; i++){
          V(i,k)=w(i,0);
         } 
   //Solve upper hessenberg least squares problem via Givens rotations
   //Compute previous Givens rotations
    for (int i=0; i<k-1; i++){
       h(i,k-1)=c(i,0)*h(i,k-1)+s(i,0)*h(i+1,k-1);
       h(i+1,k-1)=-s(i,0)*h(i,k-1)+c(i,0)*h(i+1,k-1);
     }  
     //Compute next Givens rotations
     c.reshape(k,1);
     s.reshape(k,1); 
     bb.reshape(k+1,1);
     double l = sqrt(h(k-1,k-1)*h(k-1,k-1)+h(k,k-1)*h(k,k-1));
     c(k-1,0)=h(k-1,k-1)/l;
     s(k-1,0)=h(k,k-1)/l;
     std::cout <<" h(k,k-1)= " << h(k,k-1) << std::endl;
     // Givens rotation on h and bb
     h(k-1,k-1)=l;
     h(k,k-1)=0;
     bb(k-1,0)=c(k-1,0)*bb(k-1,0);
     bb(k,0)=-s(k-1,0)*bb(k-1,0);

    //Determine residual    
    resid = fabs(bb(k,0));

    std::cout << "resid = " << resid << std::endl;
    k=k+1;
  } 
   //Extract upper triangular square matrix
   bb.reshape(h.numRows()-1 ,1);

   //Solve linear system
   int info;

   Teuchos::LAPACK<int, double> lapack;
   lapack.TRTRS('U', 'N', 'N', h.numRows()-1, 1, h.values(), h.stride(), bb.values(), bb.stride(),&info); 
   //Found y=Mx
   for (int i=0; i<k-1; i++){
      bb(i,0)=bb(i,0)/D(i,0);
   }

   V.reshape(n,k-1);
   Teuchos::SerialDenseMatrix<int, double> ans(X);
   ans.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, V, bb, 1.0);
   std::cout << "ans= " << ans << std::endl;

   std::cout << "h= " << h << std::endl;



return 0;
}

int main()
{
  Teuchos::SerialDenseMatrix<int, double> A(4,4);
  Teuchos::SerialDenseMatrix<int, double> B(4,1);
  Teuchos::SerialDenseMatrix<int, double> X(4,1);

r(5.0);
  X.putScalar(0.0);

  int maxit=20;
  double tolerance=1e-6;

  pregmres(A,X,B,maxit,tolerance);

  return 0;
}





