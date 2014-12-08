#include<iostream>
#include<vector>
#include"Teuchos_LAPACK.hpp"

#ifndef __INNER_PRODUCT_MATRIX__
#define __INNER_PRODUCT_MATRIX__

template<class Real>
class InnerProductMatrix{

    public:
        InnerProductMatrix(const std::vector<Real> &U,
                           const std::vector<Real> &V,
                           const std::vector<Real> &w);

        ~InnerProductMatrix(); 
        void apply(const std::vector<Real> &x, std::vector<Real> &b);          

    private:
        const int nq_;
        const int ni_;
        std::vector<Real> M_;     


};

template<class Real>
class InnerProductMatrixSolver : public InnerProductMatrix<Real> {
    

}; 






template<class Real>
InnerProductMatrix<Real>::InnerProductMatrix( const std::vector<Real> &U, 
                                              const std::vector<Real> &V,  
                                              const std::vector<Real> &w) : 
                                              nq_(w.size()), ni_(U.size()/nq_), M_(ni_*ni_,0) {
    for(int i=0;i<ni_;++i) {
        for(int j=0;j<ni_;++j) {
            for(int k=0;k<nq_;++k) { 
                M_[i+j*ni_] += w[k]*U[k+i*nq_]*V[k+j*nq_];
            }
        }
    }
}

template<class Real>
InnerProductMatrix<Real>::~InnerProductMatrix(){
}

template<class Real>
void InnerProductMatrix<Real>::apply(const std::vector<Real> &x, std::vector<Real> &b) {

    for(int i=0;i<ni_;++i) {
        b[i] = 0;
        for(int j=0;i<ni_;++j ) {
            b[i] += M_[i*ni_+j]*x[j]; 
        }
    } 
}         






#endif

