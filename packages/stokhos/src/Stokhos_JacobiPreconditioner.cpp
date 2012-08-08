#include "Stokhos_JacobiPreconditioner.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

Stokhos::JacobiPreconditioner::
JacobiPreconditioner(
  const Teuchos::SerialDenseMatrix<int, double>& A_) :
       A(A_)
{
}


Stokhos::JacobiPreconditioner::
~JacobiPreconditioner()
{
}

int
Stokhos::JacobiPreconditioner::
ApplyInverse(const Teuchos::SerialDenseMatrix<int, double>& Input,
                             Teuchos::SerialDenseMatrix<int, double>& Result, int m) const
{
  int n=Input.numRows();
  Teuchos::SerialDenseMatrix<int, double> G(A);
  Teuchos::SerialDenseMatrix<int, double> z(n,1);
  Teuchos::SerialDenseMatrix<int, double> invDz(n,1);


  for (int j=0; j<m; j++){

    if (j==0){  // Compute z=D-1r
      for (int i=0; i<n; i++)
         z(i,0)=Input(i,0)/A(i,i);

    }
    else {
      //Compute G=invD(-L-U)=I-inv(D)A 
      
      for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
           if (j==i)
             G(i,j)=0;
           else 
             G(i,j)=-A(i,j)/A(i,i);
         }
      }
      //z=Gz+inv(D)r
      //do inv(D)r 
      for (int i=0; i<n; i++){
        invDz(i,0)=Input(i,0)/A(i,i);
      }
      
      //Teuchos::SerialDenseMatrix<int, double> Gz(invDz);
      invDz.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0, G, z, 1.0);
      z.assign(invDz);
     
 }

      
  }
 Result.assign(z);
 return 0;
}
