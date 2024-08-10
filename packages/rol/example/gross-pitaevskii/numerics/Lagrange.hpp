// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include<vector>
#include<iostream>
#ifndef __LAGRANGE__
#define __LAGRANGE__

template<class Real>
class Lagrange{
    public:
        Lagrange(const std::vector<Real> &xin, const std::vector<Real> &xev );
        ~Lagrange();

        // Interpolate from the interpolation to the evaluation points
        void interp(const std::vector<Real> &f, std::vector<Real> &p);

        void dinterp(const std::vector<Real> &f, std::vector<Real> &p);

        // Evaluate the kth interpolant on the evaluation points
        void interpolant(const int k, std::vector<Real> &l);

        // Derivative of the interpolating polynomial
        void derivative(const int k, std::vector<Real> &d);

        /* Implement sum formulas as found in equation 4.2 of Trefethen 
           and Berrut SIAM Review, Vol. 46, No. 3, pp.501-517 */
        void bi_sum(const std::vector<Real> &f, std::vector<Real> &y);


    private:
        //! \param xin_ Vector of interpolation points  
        const std::vector<Real> xin_; 

        // \param xwv_ Vector of evaluation points
        const std::vector<Real> xev_;

        // \param nin_ Number of interpolation points
        const int nin_;

        // \param nev_ Number of evaluation points
        const int nev_;

        // \param w_ Vector of interpolation weights
        std::vector<Real> w_;

        // \param ell_ Vector conatining barycentric multiplicative polynomial on evaluation points
        std::vector<Real> ell_;

        // Pseudospectral differentiation matrix on interpolation points (column stacked)
        std::vector<Real> D_; 
};



/** \brief Interpolation object which interpolates from to the grid xin to xev  
    @param[in]  xin  vector of interpolation points
    @param[in]  xev  vector of evaluation points */ 
template<class Real>
Lagrange<Real>::Lagrange(const std::vector<Real> &xin, const std::vector<Real> &xev):
    xin_(xin), xev_(xev), nin_(xin.size()), nev_(xev.size()), w_(nin_,0), ell_(nev_,0), D_(nin_*nin_,0) {

    // For storing displacements between interpolation points
    double d;

   /* Compute the weights using as slightly modified version of the 
       algorithm on page 504 in the Trefethen & Berrut paper */
    
    w_[0] = 1;
    
    for(int j=1;j<nin_;++j)
    {
        w_[j] = 1.0;

        for(int k=0;k<j;++k)
        {
            d = xin_[k]-xin_[j];
            w_[k] *= d;
            w_[j] *= -d;
        }
    }

    for(int j=0;j<nin_;++j)
    {
        w_[j] = 1.0/w_[j];
    }
  
    std::vector<Real> ones(nin_,1.0);    // Create vector of ones
  
    this->bi_sum(ones,ell_);    // Compute the \f$\ell(x)\f$ polynomial

    for(int j=0;j<nev_;++j)
    {
        ell_[j] = 1.0/ell_[j];
    }

    for(int k=0;k<nin_;++k) {
        for(int i=0;i<nin_;++i) {
            if(i!=k) {  
                D_[i+k*nin_] = (w_[k]/w_[i])/(xin_[i]-xin_[k]); 
            }
        }       

   }
   for(int k=0;k<nin_;++k){ 
        for(int i=0;i<nin_;++i) {
            if(i!=k){
                D_[k+k*nin_] -= D_[k+i*nin_];  
            }
        }
    }
} 




template<class Real>
Lagrange<Real>::~Lagrange(){}

/** \brief This routine evaluates sums of the form shown in equation (4.2) in 
            the paper by J-P Berrut and L.N. Trefethen.  
    @param[in]   f  vector of values appearing in the sum
    @param[out]  y  the result */
template<class Real>
void Lagrange<Real>::bi_sum(const std::vector<Real> &f, std::vector<Real> &y){

    for(int j=0;j<nev_;++j)
    {
        y[j] = 0;
        for(int k=0;k<nin_;++k)
        {
            if(xev_[j] == xin_[k])
            {
                y[j] = f[j];
                break;
            }
            else
            {
                y[j] += w_[k]*f[k]/(xev_[j]-xin_[k]);
            }
        }    
    }    
}

/** \brief Given the values of a function on the interpolation points xin, stored 
           in f, evaluate the interpolating polynomial on the evaluation points xev
    @param[in]  f  vector of function values sampled at xin 
    @param[out] y  vector of interpolating polynomial values evaluated ay xev */
template<class Real>
void Lagrange<Real>::interp(const std::vector<Real> &f, std::vector<Real> &p){ 
    this->bi_sum(f,p);

    for(int j=0;j<nev_;++j)
    {
        p[j] *= ell_[j];
    }    
}

template<class Real>
void Lagrange<Real>::dinterp(const std::vector<Real> &f, std::vector<Real> &p){ 
    std::vector<Real> fb(nin_,0);

    std::copy(f.begin(),f.end(),fb.begin()+1);

    this->bi_sum(fb,p);

    for(int j=0;j<nev_;++j)
    {
        p[j] *= ell_[j];
    }    
}



/** \brief Evaluate the kth interpolant on the evaluation points 
    @param[in]  k  vector of function values sampled at xin 
    @param[out] y  vector of interpolating polynomial values evaluated ay xev */
template<class Real>
void Lagrange<Real>::interpolant(const int k, std::vector<Real> &l){ 
    std::vector<Real> f(nin_,0);
    std::fill(l.begin(),l.end(),0);
    f[k] = 1.0; 
    this->bi_sum(f,l);

    for(int j=0;j<nev_;++j)
    {
        l[j] = ell_[j]*w_[k]/(xev_[j]-xin_[k]);
    }    
}


/** \brief Derivative of the \f$k\f$th interpolant on the interpolation points */
template<class Real>
void Lagrange<Real>::derivative(const int k, std::vector<Real> &d ) {
    std::vector<Real> lp(nin_,0);
    std::fill(d.begin(),d.end(),0);
    std::copy(D_.begin()+k*nin_,D_.begin()+(k+1)*nin_,lp.begin());

    // Interpolate the derivative
    this->interp(lp,d);
}

#endif
