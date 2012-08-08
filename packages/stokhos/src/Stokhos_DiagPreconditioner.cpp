#include "Stokhos_DiagPreconditioner.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

Stokhos::DiagPreconditioner::
DiagPreconditioner(
  const Teuchos::SerialDenseMatrix<int, double>& A_) :
       A(A_)
{
}


Stokhos::DiagPreconditioner::
~DiagPreconditioner()
{
}

int
Stokhos::DiagPreconditioner::
ApplyInverse(const Teuchos::SerialDenseMatrix<int, double>& Input,
                             Teuchos::SerialDenseMatrix<int, double>& Result, int m) const
{
  int n=Input.numRows();
  for (int i=0; i<n; i++){
      Result(i,0)=Input(i,0)/A(i,i);
      }
  return 0;
}
