#include<iostream>
#include<iomanip>
#include"numerics/NodalBasis.hpp"
#include"numerics/InnerProductMatrix.hpp" 


typedef double RealT;

int main(int argc, char* argv[]) {
  
    const Teuchos::LAPACK<int,RealT> * const lapack = new Teuchos::LAPACK<int,RealT>();
    NodalBasis<RealT> nb(lapack,3,6);
    InnerProductMatrix<RealT> mass(nb.L_,nb.L_,nb.wq_);
    InnerProductMatrix<RealT> stiffness(nb.Lp_,nb.Lp_,nb.wq_);

    return 0;
}
