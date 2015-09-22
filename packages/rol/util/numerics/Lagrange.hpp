// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#include "ROL_StdVector.hpp"

#ifndef __LAGRANGE__
#define __LAGRANGE__

template<class Real>
class Lagrange{
    public:
        Lagrange(const ROL::Vector<Real> &xin, const ROL::Vector<Real> &xev );

        // Interpolate from the interpolation to the evaluation points
        void interp(const ROL::Vector<Real> &f, ROL::Vector<Real> &p);

        // Interpolate a function specified only on interior points and
        // assumed zero on the boundary
        void dinterp(const ROL::Vector<Real> &f, ROL::Vector<Real> &p);

        // Evaluate the kth interpolating polynomial on the evaluation points
        void interpolant(const int k, ROL::Vector<Real> &l);

        // Derivative of the kth interpolating polynomial
        void derivative(const int k, ROL::Vector<Real> &d);

        // Differentiate a function on the interpolation points
        // void diff(const ROL::Vector<Real> &f, ROL::Vector<Real> &df);


    private:
        //! \param xin_ Vector of interpolation points  
        Teuchos::RCP<const std::vector<Real> > xip_; 

        // \param xwv_ Vector of evaluation points
        Teuchos::RCP<const std::vector<Real> > xep_;

        // \param nin_ Number of interpolation points
        int nin_;

        // \param nev_ Number of evaluation points
        int nev_;

        // \param w_ Vector of interpolation weights
        Teuchos::RCP<std::vector<Real> > wp_;

        // \param ell_ Vector conatining barycentric multiplicative polynomial on evaluation points
        Teuchos::RCP<std::vector<Real> > ellp_;

        // Pseudospectral differentiation matrix on interpolation points (column stacked)
        Teuchos::RCP<std::vector<Real> > Dp_; 

        /* Implement sum formulas as found in equation 4.2 of Trefethen 
           and Berrut SIAM Review, Vol. 46, No. 3, pp.501--517 */
        void bi_sum( Teuchos::RCP<const std::vector<Real> > &fp, 
                     Teuchos::RCP<std::vector<Real> > &yp );

};



/** \brief Interpolation object which interpolates from to the grid xin to xev  
    @param[in]  xin  vector of interpolation points
    @param[in]  xev  vector of evaluation points */ 
template<class Real>
Lagrange<Real>::Lagrange(const ROL::Vector<Real> &xin, const ROL::Vector<Real> &xev) :
    xip_((Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(xin))).getVector()),
    xep_((Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(xev))).getVector()),
    nin_(xip_->size()),
    nev_(xep_->size()),
    wp_(Teuchos::rcp( new std::vector<Real>(nin_,0) )),     
    ellp_(Teuchos::rcp( new std::vector<Real>(nev_,0) )),
    Dp_(Teuchos::rcp( new std::vector<Real>(nin_*nin_,0) )) {


    // For storing displacements between interpolation points
    Real d;

    // Begin computing weights using as slightly modified version of the 
    // algorithm on page 504 in the Trefethen & Berrut paper 

    (*wp_)[0] = 1.0;
    
    for(int j=1;j<nin_;++j)
    {
        (*wp_)[j] = 1.0;

        for(int k=0;k<j;++k)
        {
            d = (*xip_)[k]-(*xip_)[j];
            (*wp_)[k] *=  d;
            (*wp_)[j] *= -d;
        }
    }

    for(int j=0;j<nin_;++j)
    {
        (*wp_)[j] = 1.0/(*wp_)[j];
    }

    // ---- End computing weights ---

    // Create vector of ones
    Teuchos::RCP<const std::vector<Real> > onesp = Teuchos::rcp( new std::vector<Real>(nin_,1.0) ); 

    this->bi_sum(onesp,ellp_);    // Compute the \f$\ell(x)\f$ polynomial

    for(int j=0;j<nev_;++j)
    {
        (*ellp_)[j] = 1.0/(*ellp_)[j];
    }

 
    for(int k=0;k<nin_;++k) {
        for(int i=0;i<nin_;++i) {
            if(i!=k) {  
                (*Dp_)[i+k*nin_] = ((*wp_)[k]/(*wp_)[i])/((*xip_)[i]-(*xip_)[k]); 
            }
        }       

   }
   for(int k=0;k<nin_;++k){ 
        for(int i=0;i<nin_;++i) {
            if(i!=k){
                (*Dp_)[k+k*nin_] -= (*Dp_)[k+i*nin_];  
            }
        }
    }
} 





/** \brief This routine evaluates sums of the form shown in equation (4.2) in 
            the paper by J-P Berrut and L.N. Trefethen.  
    @param[in]   f  vector of values appearing in the sum
    @param[out]  y  the result */
template<class Real>
void Lagrange<Real>::bi_sum( Teuchos::RCP<const std::vector<Real> >&fp, 
                             Teuchos::RCP<std::vector<Real> > &yp ){

    for(int j=0;j<nev_;++j)
    {
        (*yp)[j] = 0;
        for(int k=0;k<nin_;++k)
        {
            if((*xep_)[j] == (*xip_)[k])
            {
                (*yp)[j] = (*fp)[j];
                break;
            }
            else
            {
                (*yp)[j] += (*wp_)[k]*(*fp)[k]/((*xep_)[j]-(*xip_)[k]);
            }
        }    
    }    
}

/** \brief Given the values of a function on the interpolation points xin, stored 
           in f, evaluate the interpolating polynomial on the evaluation points xev
    @param[in]  f  vector of function values sampled at xin 
    @param[out] y  vector of interpolating polynomial values evaluated ay xev */
template<class Real>
void Lagrange<Real>::interp(const ROL::Vector<Real> &f, ROL::Vector<Real> &p){ 
    Teuchos::RCP<const std::vector<Real> > fp = 
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(f))).getVector();

    Teuchos::RCP<std::vector<Real> > pp = 
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(p)).getVector()); 

    this->bi_sum(fp,pp);

    for(int j=0;j<nev_;++j)
    {
        (*pp)[j] *= (*ellp_)[j];
    }    
}

template<class Real>
void Lagrange<Real>::dinterp(const ROL::Vector<Real> &f, ROL::Vector<Real> &p){ 

    Teuchos::RCP<const std::vector<Real> > fp = 
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(f))).getVector();

    Teuchos::RCP<std::vector<Real> > pp = 
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(p)).getVector()); 

    this->bi_sum(fp,pp);

    for(int j=0;j<nev_;++j)
    {
        (*pp)[j] *= (*ellp_)[j];
    }    
}



/** \brief Evaluate the kth interpolant on the evaluation points 
    @param[in]  k  vector of function values sampled at xin 
    @param[out] y  vector of interpolating polynomial values evaluated ay xev */
template<class Real>
void Lagrange<Real>::interpolant(const int k, ROL::Vector<Real> &l){ 

    Teuchos::RCP<std::vector<Real> > lp = 
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(l)).getVector()); 

    std::fill(lp->begin(),lp->end(),0);

    Teuchos::RCP<std::vector<Real> > ekp = Teuchos::rcp(new std::vector<Real>(nin_,0));
    (*ekp)[k] = 1.0; 

    ROL::StdVector<Real> ek(ekp);
    this->interp(ek,l);

}


/** \brief Derivative of the \f$k\f$th interpolant on the evaluation points */
template<class Real>
void Lagrange<Real>::derivative(const int k, ROL::Vector<Real> &d ) {

    Teuchos::RCP<std::vector<Real> > dp = 
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(d)).getVector()); 
   
    Teuchos::RCP<std::vector<Real> > lp = Teuchos::rcp( new std::vector<Real>(nin_,0) );
    std::fill(dp->begin(),dp->end(),0);
    std::copy(Dp_->begin()+k*nin_,Dp_->begin()+(k+1)*nin_,lp->begin());

    ROL::StdVector<Real> l(lp);

    // Interpolate the derivative
    this->interp(l,d);
}



#endif
