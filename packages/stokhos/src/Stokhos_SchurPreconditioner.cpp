#include "Stokhos_SchurPreconditioner.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_RCP.hpp"


Stokhos::SchurPreconditioner::
SchurPreconditioner(
  const Teuchos::SerialDenseMatrix<int, double>& K_, const int p_, const int m_, const int diag_) :
       K(K_),
       p(p_),
       m(m_),
       diag(diag_)
{
}


Stokhos::SchurPreconditioner::
~SchurPreconditioner()
{
}

int fact (int n)
{
  if (n > 1)
   n= n * fact (n-1);
  else
   n = 1;
  return n;
 
}

int size (int n, int m) //n is the polynomial order and m is the number of random variables
{
 // return (fact(n+m)/(fact(n)*fact(m)));
    int min, max;
    if (n == 0 ) 
        return 1;
    else {
	if (n<=m){
	   max = m;
           min = n;
        }
        else {
           max = n;
           min = m;
        }
    
    int num = n+m;
    for (int i=1; i<=min-1; i++)
        num = num*(n+m-i);
    return num/fact(min); 
 }
}
int
Stokhos::SchurPreconditioner::
ApplyInverse(const Teuchos::SerialDenseMatrix<int, double>& Input,
                             Teuchos::SerialDenseMatrix<int, double>& U, const int n) const   //p: polynomial order; m: number of random variables; n: level to stop and form actual Schur Complement
{
  Teuchos::SerialDenseMatrix<int, double> Resid(Input);
  int ret;
  int lmin;
  

  if (n<=1)
      lmin=1;  
  else 
      lmin=n;

   //Create array of solvers to solve Dr=r
/*   int sz=p-lmin+1;
   Array< RCP<const Teuchos::SerialDenseSolver<int,double> > > Solvers(sz);
    for (int i=0; i<sz; i++) {
      Solvers[i] = rcp(new Stokhos::LegendreBasis<int,double>(p));
    }
*/          
   Teuchos::SerialDenseSolver<int, double> solver;
   for (int l=p; l>=lmin; l--){
//   Teuchos::SerialDenseSolver<int, double> solver;    
   int c=size(l,m);
   int s = size(l-1,m);   
//   std::cout << " c = " << c << " s = " << s << std::endl;
  
   //split residual
   Teuchos::SerialDenseMatrix<int, double> rminus(Teuchos::View, Resid, s, 1);
   Teuchos::SerialDenseMatrix<int, double> r(Teuchos::Copy, Resid, c-s, 1, s, 0);
   
   //Compute pre-correction

   Teuchos::SerialDenseMatrix<int, double> B(Teuchos::View, K, s, c-s, 0, s);   
   Teuchos::SerialDenseMatrix<int, double> D(Teuchos::View, K, c-s, c-s, s,s); 

 
   //Computing inv(D)r
   if (diag == 0){
      //For D diagonal
      for (int i=0; i<c-s; i++)
	r(i,0)=r(i,0)/D(i,i);

   }
   else{
   //For full D block
      Teuchos::RCP< Teuchos::SerialDenseMatrix<int, double> > DD, RR;
      DD = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int, double> (D));
      
      RR = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int, double> (r));
     

     // Setup solver
     solver.setMatrix(DD);
     solver.setVectors(RR, RR);
     //Solve D*r=r
     if (solver.shouldEquilibrate()) {
        solver.factorWithEquilibration(true);
        solver.equilibrateMatrix();
        }
        solver.solve();
    
  
      for (int i=0; i<c-s; i++){
             r(i,0)=(*RR)(i,0);
       }
 
   }


   Teuchos::SerialDenseMatrix<int, double> rm(rminus);  
   ret = rm.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS, -1.0, B, r, 1.0);
   TEUCHOS_ASSERT(ret == 0);
   
 
   //set r(l-1)=g(l-1)
   if (l>lmin){
     for (int i=0; i<s; i++)
        Resid(i,0)=rm(i,0); 
   

   }

   else {
     //For l=1, solve A0u0=g0
     if (n<1){

        U(0,0)=rm(0,0)/K(0,0);
     }
     //Compute the Schur completement
     else {
      for (int i=0; i<s; i++)
         Resid(i,0)=rm(i,0);
   //   Teuchos::SerialDenseSolver<int, double> solver;
      Teuchos::SerialDenseMatrix<int, double> S(Teuchos::Copy, K, s, s);    
      Teuchos::SerialDenseMatrix<int, double> BinvD(s,c-s);   
     
	    

        for (int i=0; i<c-s; i++)
          for (int j=0; j<s; j++)
              BinvD(j,i)=B(j,i)/D(i,i);





 
      S.multiply(Teuchos::NO_TRANS,Teuchos::TRANS, -1.0, BinvD, B, 1.0);
 
      Teuchos::SerialDenseMatrix<int, double> Ul(Teuchos::View, U, s, 1);
      Teuchos::RCP< Teuchos::SerialDenseMatrix<int, double> > SS, UU, RR; 
      SS = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int, double> (S));
      UU = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int, double> (Ul));
      RR = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int, double> (rm));
      Teuchos::SerialDenseSolver<int, double> solver2;

//      std::cout << "S  = " << S << std::endl;
      // Setup solver
      solver2.setMatrix(SS);
      solver2.setVectors(UU, RR);
      //Solve S*u=rm
      if (solver2.shouldEquilibrate()) {
        solver2.factorWithEquilibration(true);
        solver2.equilibrateMatrix();
         }
         solver2.solve();
   
      for (int i=0; i<s; i++)
         U(i,0)=(*UU)(i,0);
      }

}
  

}


   
  for (int l=lmin; l<=p; l++){
   //compute post-correction
   int c = size(l,m);
   int s = size(l-1,m);
 
   Teuchos::SerialDenseMatrix<int, double> B(Teuchos::View, K, s, c-s, 0, s);
   Teuchos::SerialDenseMatrix<int, double> D(Teuchos::Copy, K, c-s, c-s, s,s);
   //split residual
   Teuchos::SerialDenseMatrix<int, double> r(Teuchos::Copy, Resid, c-s, 1, s, 0);
   
   //Get first s values of U
   Teuchos::SerialDenseMatrix<int, double> Ul(Teuchos::View, U, s, 1);
   
   ret = r.multiply(Teuchos::TRANS,Teuchos::NO_TRANS, -1.0, B, Ul, 1.0);

   if (diag == 0) {
   //For diag D
   //Concatenate U
     for (int i=s; i<c; i++){
        U(i,0)=r(-s+i,0)/D(-s+i,-s+i);
     }
   }
   else {
   //For block D
//      Teuchos::SerialDenseSolver<int, double> solver2;
      Teuchos::RCP< Teuchos::SerialDenseMatrix<int, double> > DD, RR;
      DD = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int, double> (D));
      RR = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int, double> (r));
     // Setup solver
     solver.setMatrix(DD);
     solver.setVectors(RR, RR);
     //Solve D*r=r
     if (solver.shouldEquilibrate()) {
          solver.factorWithEquilibration(true);
          solver.equilibrateMatrix();
          }
          solver.solve();

     for (int i=s; i<c; i++){
       U(i,0)=(*RR)(-s+i,0); 
       } 
   }

 }
  
  return 0;
}

 
