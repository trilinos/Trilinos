#include "Stokhos_BlockPreconditioner.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"

//Computes the exact Schur complement block LU decomposition


Stokhos::BlockPreconditioner::
BlockPreconditioner(
  const Teuchos::SerialDenseMatrix<int, double>& K_,const int p_, const int m_) : 
       K(K_),
       p(p_),
       m(m_)
{
}


Stokhos::BlockPreconditioner::
~BlockPreconditioner()
{
}

int facto(int n)

{
  if (n > 1)
   return (n * facto(n-1));
  else
   return (1);
}

int siz (int n, int m) //n is the polynomial order and m is the number of random variables
{
    return (facto(n+m)/(facto(n)*facto(m)));
 }


int
Stokhos::BlockPreconditioner::
ApplyInverse(const Teuchos::SerialDenseMatrix<int, double>& Input,
                             Teuchos::SerialDenseMatrix<int, double>& Result, int n) const
{ //Solve M*Result=Input
     int c=siz(p,m);
     int s = siz(p-1,m);
  
     //Split residual
     Teuchos::SerialDenseMatrix<int, double> r1(Teuchos::Copy, Input, s, 1);
     Teuchos::SerialDenseMatrix<int, double> r2(Teuchos::Copy, Input, c-s, 1, s, 0);
     
     //Split Result     
     Teuchos::SerialDenseMatrix<int, double> u1(Teuchos::Copy, Result, s, 1);
     Teuchos::SerialDenseMatrix<int, double> u2(Teuchos::Copy, Result, c-s, 1, s, 0);

     Teuchos::SerialDenseMatrix<int, double> B(Teuchos::View, K, s, c-s, 0, s);
     Teuchos::SerialDenseMatrix<int, double> D(Teuchos::View, K, c-s, c-s, s,s);

     //rD=inv(D)r2

     Teuchos::SerialDenseMatrix<int, double> Dr(c-s,1);
      
     for (int i=0; i<c-s; i++)
         Dr(i,0)=r2(i,0)/D(i,i);
     
     int ret = r1.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS, -1.0, B, Dr, 1.0);
     TEUCHOS_ASSERT(ret == 0);

      //Compute S=A-B*inv(D)*Bt
      Teuchos::SerialDenseMatrix<int, double> S(Teuchos::Copy, K, s, s);
      //Compute B*inv(D)
      Teuchos::SerialDenseMatrix<int, double> BinvD(s,c-s);
      for (int i=0; i<c-s; i++) //col
          for (int j=0; j<s; j++) //row
              BinvD(j,i)=B(j,i)/D(i,i);
     
      S.multiply(Teuchos::NO_TRANS,Teuchos::TRANS, -1.0, BinvD, B, 1.0);
    
      Teuchos::RCP< Teuchos::SerialDenseMatrix<int, double> > SS, w, rr;
      SS = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int, double> (S));
      w = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int, double> (s,1));
      rr = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int, double> (r1));
     
      
      // Setup solver
      Teuchos::SerialDenseSolver<int, double> solver;
      solver.setMatrix(SS);
      solver.setVectors(w, rr);
      //Solve S*w=r1
      if (solver.shouldEquilibrate()) {
         solver.factorWithEquilibration(true);
         solver.equilibrateMatrix();
         }
      solver.solve();
     
      for (int i=0; i<s; i++)
          Result(i,0)=(*w)(i,0);
     
       

      ret = r2.multiply(Teuchos::TRANS,Teuchos::NO_TRANS, -1.0, B, *w, 1.0);
      TEUCHOS_ASSERT(ret == 0);
          
      for (int i=s; i<c; i++)
          Result(i,0)=r2(-s+i,0)/D(-s+i, -s+i);

   
  
 return 0;
}
