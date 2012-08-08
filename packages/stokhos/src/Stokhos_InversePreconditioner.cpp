#include "Stokhos_InversePreconditioner.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"


Stokhos::InversePreconditioner::
InversePreconditioner(
  const Teuchos::SerialDenseMatrix<int, double>& A_) :
       A(A_)
{
}


Stokhos::InversePreconditioner::
~InversePreconditioner()
{
}

int
Stokhos::InversePreconditioner::
ApplyInverse(const Teuchos::SerialDenseMatrix<int, double>& Input,
                             Teuchos::SerialDenseMatrix<int, double>& Result, int m) const
{
      Teuchos::RCP< Teuchos::SerialDenseMatrix<int, double> > AA, UU, RR;
      AA = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int, double> (A));
      UU = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int, double> (Result));
      RR = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int, double> (Input));

      
      // Setup solver
      Teuchos::SerialDenseSolver<int, double> solver;
      solver.setMatrix(AA);
      solver.setVectors(UU, RR);
      //Solve A*Result=Input
      if (solver.shouldEquilibrate()) {
         solver.factorWithEquilibration(true);
         solver.equilibrateMatrix();
         }
         solver.solve();
      

      for (int i=0; i<A.numRows(); i++)
          Result(i,0)=(*UU)(i,0);
      
      std::cout << "Result = " << Result << std::endl;  
 return 0;
}
