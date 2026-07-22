// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include<iostream>
#include"ROL_StdVector.hpp"
#include"OrthogonalPolynomials.hpp"
#include"Lagrange.hpp"
#include"LinearAlgebra.hpp"

#ifndef __NODAL_BASIS__
#define __NODAL_BASIS__


template<class Real>
struct NodalBasis{

        //! \param lapack_ pointer to LAPACK interface
        ROL::Ptr<Teuchos::LAPACK<int,Real> > lapack_;
 
        //! \param ni_ Number of interpolation points
        const int ni_;

        //! \param nq_ Number of quadrature points
        const int nq_;

        NodalBasis(ROL::Ptr<Teuchos::LAPACK<int,Real> > lapack, const int ni, const int nq);    
        ~NodalBasis();

         //! \param xi_ Vector of interpolation points
        std::vector<Real> xi_;

        //! \param xq_ Vector of quadrature points
        std::vector<Real> xq_;

        //! \param wq_ Vector of quadrature weights
        std::vector<Real> wq_;

        //! \param L_ Column stacked vector of vectors of interpolating functions  
        std::vector<Real> L_;
        std::vector<Real> Lp_;
        
        //! Object for working with Lagrange polynomials and their derivatives 
        ROL::Ptr<Lagrange<Real> > lagrange_; 

         
};

/** \brief Set up quantities we will need repeatedly
*/
template<class Real>
NodalBasis<Real>::NodalBasis(ROL::Ptr<Teuchos::LAPACK<int,Real> > const lapack, const int ni, const int nq): 
    lapack_(lapack), ni_(ni), nq_(nq), xi_(ni_,0), xq_(nq_,0), wq_(nq_,0), L_(ni_*nq_,0), Lp_(ni_*nq_,0)
    {

    // Generate Legendre-Gauss-Lobatto nodes and weights (unused)
    std::vector<Real> ai(ni_,0);
    std::vector<Real> bi(ni_,0);
    std::vector<Real> wi(ni_,0);
    rec_jacobi(this->lapack_,0,0,ai,bi); 
    rec_lobatto(this->lapack_,-1.0,1.0,ai,bi);
    gauss(this->lapack_,ai,bi,xi_,wi);

    // Generate Legendre-Gauss nodes and weights
    std::vector<Real> aq(nq_,0);
    std::vector<Real> bq(nq_,0);
    rec_jacobi(this->lapack_,0,0,aq,bq); 
    gauss(this->lapack_,aq,bq,xq_,wq_);

    std::vector<Real> e(ni_,0);
    std::vector<Real> ell(nq,0);

    lagrange_ = ROL::makePtr<Lagrange<Real>>(xi_,xq_);

    // Loop over canonical vectors
    for(int i=0;i<ni_;++i) {

        lagrange_->interpolant(i,ell);

        std::copy(ell.begin(),ell.end(),L_.begin()+i*nq_);
      

        lagrange_->derivative(i,ell);
        std::copy(ell.begin(),ell.end(),Lp_.begin()+i*nq_);

        // Rezero the vector 
        std::fill(e.begin(),e.end(),0);
    }     

}

template<class Real>
NodalBasis<Real>::~NodalBasis(){
}

#endif
