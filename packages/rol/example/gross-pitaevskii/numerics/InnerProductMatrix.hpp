#include<iostream>
#include<vector>
#include"ROL_StdVector.hpp"
#include"Teuchos_LAPACK.hpp"

#ifndef __INNER_PRODUCT_MATRIX__
#define __INNER_PRODUCT_MATRIX__


template<class Real>
class InnerProductMatrix{
    public:
        InnerProductMatrix(const std::vector<Real> &U,
                           const std::vector<Real> &V,
                           const std::vector<Real> &w,
                           int a);

        InnerProductMatrix(const std::vector<Real> &U,
                           const std::vector<Real> &V,
                           const std::vector<Real> &w,
                           const std::vector<Real> &a);

        InnerProductMatrix(InnerProductMatrix<Real> *ipm);

        void update(const std::vector<Real> &a);

        virtual ~InnerProductMatrix(); 

        void apply(Teuchos::RCP<const std::vector<Real> > xp,
                   Teuchos::RCP<std::vector<Real> > bp);          
   
        void applyadd(Teuchos::RCP<const std::vector<Real> > xp,
                      Teuchos::RCP<std::vector<Real> > bp);          

        void applyaddtimes(Teuchos::RCP<const std::vector<Real> > xp,
                           Teuchos::RCP<std::vector<Real> > bp, Real factor);          


        Real inner(Teuchos::RCP<const std::vector<Real> > up,
                   Teuchos::RCP<const std::vector<Real> > vp); 

        // This method does nothing in the base class
        virtual void solve(Teuchos::RCP<const std::vector<Real> > bp, 
                           Teuchos::RCP<std::vector<Real> > xp){};    

        // This method does nothing in the base class
        virtual Real inv_inner(Teuchos::RCP<const std::vector<Real> > up,
                               Teuchos::RCP<const std::vector<Real> > vp){
                               return 0;} 

    protected:
        const int nq_;
        const int ni_;
        const std::vector<Real> U_;
        const std::vector<Real> V_;
        const std::vector<Real> w_;  

        std::vector<Real> M_;     


};



template<class Real>
InnerProductMatrix<Real>::InnerProductMatrix( const std::vector<Real> &U, 
                                              const std::vector<Real> &V,  
                                              const std::vector<Real> &w,
                                              const int a=1 ) : 
                                              nq_(w.size()), ni_(U.size()/nq_), 
                                              U_(U),V_(V),w_(w),M_(ni_*ni_,0) {
    for(int i=0;i<ni_;++i) {
        for(int j=0;j<ni_;++j) {
            for(int k=0;k<nq_;++k) { 
                M_[i+j*ni_] += a*w_[k]*U_[k+i*nq_]*V_[k+j*nq_];
            }
        }
    }
}

template<class Real>
InnerProductMatrix<Real>::InnerProductMatrix( const std::vector<Real> &U, 
                                              const std::vector<Real> &V,  
                                              const std::vector<Real> &w,
                                              const std::vector<Real> &a ) : 
                                              nq_(w.size()), ni_(U.size()/nq_),
                                              U_(U),V_(V),w_(w),M_(ni_*ni_,0) {
    for(int i=0;i<ni_;++i) {
        for(int j=0;j<ni_;++j) {
            for(int k=0;k<nq_;++k) { 
                M_[i+j*ni_] += a[k]*w_[k]*U_[k+i*nq_]*V_[k+j*nq_];
            }
        }
    }
}

template<class Real>
InnerProductMatrix<Real>::~InnerProductMatrix(){
}

template<class Real>
void InnerProductMatrix<Real>::apply(Teuchos::RCP<const std::vector<Real> > xp,
                                     Teuchos::RCP<std::vector<Real> > bp ) {
    for(int i=0;i<ni_;++i) {
        (*bp)[i] = 0;
        for(int j=0;j<ni_;++j ) {
            (*bp)[i] += M_[i+ni_*j]*(*xp)[j]; 
        }
    } 
}         

template<class Real>
void InnerProductMatrix<Real>::applyadd(Teuchos::RCP<const std::vector<Real> > xp,
                                        Teuchos::RCP<std::vector<Real> > bp ) {
    for(int i=0;i<ni_;++i) {
        for(int j=0;j<ni_;++j ) {
            (*bp)[i] += M_[i+ni_*j]*(*xp)[j]; 
        }
    } 
} 

template<class Real>
void InnerProductMatrix<Real>::applyaddtimes(Teuchos::RCP<const std::vector<Real> > xp,
                                             Teuchos::RCP<std::vector<Real> > bp, Real factor ) {
    for(int i=0;i<ni_;++i) {
        for(int j=0;j<ni_;++j ) {
            (*bp)[i] += factor*M_[i+ni_*j]*(*xp)[j]; 
        }
    } 
} 

template<class Real>
void InnerProductMatrix<Real>::update( const std::vector<Real> &a ){

    std::fill(M_.begin(),M_.end(),0);
    for(int i=0;i<ni_;++i) {
        for(int j=0;j<ni_;++j) {
            for(int k=0;k<nq_;++k) { 
                M_[i+j*ni_] += a[k]*w_[k]*U_[k+i*nq_]*V_[k+j*nq_];
            }
        }
    }
}





//! Compute the inner product \f$u^\top M v\f$
template<class Real>
Real InnerProductMatrix<Real>::inner( Teuchos::RCP<const std::vector<Real> > up,
                                      Teuchos::RCP<const std::vector<Real> > vp ) {
    Real J = 0;
    Teuchos::RCP<std::vector<Real> > Mvp = Teuchos::rcp( new std::vector<Real> (ni_,0) );
    this->apply(vp,Mvp);
    for(int i=0;i<ni_;++i) {
        J += (*up)[i]*(*Mvp)[i];
    } 
    return J;
}





//! \brief This class adds a solve method
template<class Real>
class InnerProductMatrixSolver : public InnerProductMatrix<Real> {

    private:
        Teuchos::RCP<Teuchos::LAPACK<int,Real> > lapack_; 
        const int         ni_;
        const int         nq_;
        std::vector<Real> M_; 
        const char        TRANS_;
        std::vector<int>  ipiv_;
        std::vector<Real> PLU_;     
        const int         nrhs_; 
        int               info_;

    // Solve the system Ax=b for x
    public:
       InnerProductMatrixSolver(Teuchos::RCP<Teuchos::LAPACK<int,Real> > lapack,
                                const std::vector<Real> &U, 
                                const std::vector<Real> &V,  
                                const std::vector<Real> &w,
                                const int a );

       InnerProductMatrixSolver(Teuchos::RCP<Teuchos::LAPACK<int,Real> > lapack,
                                const std::vector<Real> &U, 
                                const std::vector<Real> &V,  
                                const std::vector<Real> &w,
                                const std::vector<Real> &a);

       void solve(Teuchos::RCP<const std::vector<Real> > bp, 
                  Teuchos::RCP<std::vector<Real> > xp);    

       Real inv_inner(Teuchos::RCP<const std::vector<Real> > up, 
                      Teuchos::RCP<const std::vector<Real> > vp); 
}; 


template<class Real>
InnerProductMatrixSolver<Real>::InnerProductMatrixSolver(Teuchos::RCP<Teuchos::LAPACK<int,Real> > lapack,  
                                                         const std::vector<Real> &U=std::vector<Real>(), 
                                                         const std::vector<Real> &V=std::vector<Real>(),  
                                                         const std::vector<Real> &w=std::vector<Real>(),
                                                         const int a=1 ):
                                                         InnerProductMatrix<Real>(U,V,w,a),
                                                         lapack_(lapack),   
                                                         ni_(InnerProductMatrix<Real>::ni_),
                                                         nq_(InnerProductMatrix<Real>::nq_),
                                                         M_(InnerProductMatrix<Real>::M_),  
                                                         TRANS_('N'), ipiv_(ni_,0), PLU_(ni_*ni_,0), 
                                                         nrhs_(1),info_(0){
    PLU_ = M_;

    // Do matrix factorization
    lapack->GETRF(ni_,ni_,&PLU_[0],ni_,&ipiv_[0],&info_);

} 

template<class Real>
InnerProductMatrixSolver<Real>::InnerProductMatrixSolver(Teuchos::RCP<Teuchos::LAPACK<int,Real> > lapack,  
                                                         const std::vector<Real> &U=std::vector<Real>(), 
                                                         const std::vector<Real> &V=std::vector<Real>(),  
                                                         const std::vector<Real> &w=std::vector<Real>(),
                                                         const std::vector<Real> &a=std::vector<Real>()):
                                                         InnerProductMatrix<Real>(U,V,w,a),
                                                         lapack_(lapack),   
                                                         ni_(InnerProductMatrix<Real>::ni_),
                                                         nq_(InnerProductMatrix<Real>::nq_),
                                                         M_(InnerProductMatrix<Real>::M_),  
                                                         TRANS_('N'), ipiv_(ni_,0), PLU_(ni_*ni_,0), 
                                                         nrhs_(1),info_(0){
    PLU_ = M_;

    // Do matrix factorization
    lapack->GETRF(ni_,ni_,&PLU_[0],ni_,&ipiv_[0],&info_);

} 

//! \brief solve \f$Mx=b\f$ for \f$x\f$ 
template<class Real>
void InnerProductMatrixSolver<Real>::solve(Teuchos::RCP<const std::vector<Real> > bp, 
                                           Teuchos::RCP<std::vector<Real> > xp){

    int nrhs = bp->size()/ni_;
 
    *xp = *bp; 
 
    // Solve LU-factored system 
    lapack_->GETRS(TRANS_,ni_,nrhs,&PLU_[0],ni_,&ipiv_[0],&(*xp)[0],ni_,&info_);

}

//! Compute the inner product \f$u^\top M^{-1} v\f$
template<class Real>
Real InnerProductMatrixSolver<Real>::inv_inner( Teuchos::RCP<const std::vector<Real> > up,
                                                Teuchos::RCP<const std::vector<Real> > vp ) {
    Real J = 0;
    Teuchos::RCP<std::vector<Real> > Mivp = Teuchos::rcp( new std::vector<Real> (ni_,0) );
    this->solve(vp,Mivp);
    for(int i=0;i<ni_;++i) {
        J += (*up)[i]*(*Mivp)[i];
    } 
    return J;
}

#endif

