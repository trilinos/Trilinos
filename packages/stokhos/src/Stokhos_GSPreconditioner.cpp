#include "Stokhos_GSPreconditioner.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_LAPACK.hpp"


Stokhos::GSPreconditioner::
GSPreconditioner(
  const Teuchos::SerialDenseMatrix<int, double>& A_, const int sym_) :
       A(A_),
       sym(sym_)
{
}


Stokhos::GSPreconditioner::
~GSPreconditioner()
{
}
//m is the number of GS iterations to solve Mz=r
int
Stokhos::GSPreconditioner::
ApplyInverse(const Teuchos::SerialDenseMatrix<int, double>& Input,
                             Teuchos::SerialDenseMatrix<int, double>& Result, int m) const   //If sym=0 then do symmetric Gauss Seidel
{
  Result.assign(Input);
  int n=A.numRows();
  int info;
  Teuchos::LAPACK<int, double> lapack;
  
  //Get lower triangular part of A  
  Teuchos::SerialDenseMatrix<int, double> L(A);
  for (int i=0; i<n; i++){
     for (int j=0; j<n; j++){
         if (j>i) 
             L(i,j)=0;
         }
  }
  
  if (sym==0){
  //Get inv(diag(A))=D
  Teuchos::SerialDenseMatrix<int, double> D(n,n);
  for (int i=0; i<n; i++){
     for (int j=0; j<n; j++){
         D(i,i)=1/A(i,i);
         }
  }

  //Get upper triangular part of A  
  Teuchos::SerialDenseMatrix<int, double> U(A);
  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
       if (i>j)
            U(i,j)=0;
        }
       }
  

  Result.assign(Input);
  for (int j=0; j<m; j++){
    //compute M=(L+D)inv(diagA)
    Teuchos::SerialDenseMatrix<int, double> M(n,n);
    M.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0, L, D, 0.0);
           
    //Forward solve Lz=r
    lapack.TRTRS('L', 'N', 'N', M.numRows(), 1, M.values(), M.stride(), Result.values(), Result.stride(),&info);
    
    //Backward solve Uw=z
    lapack.TRTRS('U', 'N', 'N', U.numRows(), 1, U.values(), U.stride(), Result.values(), Result.stride(),&info);
   }
   
    }
   else{
    lapack.TRTRS('L', 'N', 'N', L.numRows(), 1, L.values(), L.stride(), Result.values(), Result.stride(),&info);
    }
  return 0;
}
