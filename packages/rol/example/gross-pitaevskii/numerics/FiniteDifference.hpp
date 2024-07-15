// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ROL_StdVector.hpp"
#include "Teuchos_LAPACK.hpp"

using namespace ROL;

template<class Real>
class FiniteDifference {
    
    private:
        int n_;
        double dx2_;
        Teuchos::LAPACK<int,Real> lapack_;

        // subdiagonal is -1/dx^2 
        std::vector<Real> dl_;

        // diagonal is 2/dx^2
        std::vector<Real> d_;

        // superdiagonal is -1/dx^2 
        std::vector<Real> du_;

        // second superdiagonal (du2 = 0)
        std::vector<Real> du2_;

            // Pivot indices
        std::vector<int> ipiv_; 
              
        int info_;
   

    public: 
 
        FiniteDifference(int n, double dx) : n_(n),dx2_(dx*dx),
            dl_(n_-1,-1.0/dx2_),
            d_(n_,2.0/dx2_),
            du_(n_-1,-1.0/dx2_),
            du2_(n_-2,0.0),
            ipiv_(n_,0) {

            lapack_.GTTRF(n_,&dl_[0],&d_[0],&du_[0],&du2_[0],&ipiv_[0],&info_);            
        }    

        //! Given f, compute -u''=f 
        void solve(ROL::Ptr<const std::vector<Real> > fp, ROL::Ptr<std::vector<Real> > up) {            
            for(int i=0;i<n_;++i) {
                (*up)[i] = (*fp)[i];
            }
            lapack_.GTTRS('n',n_,1,&dl_[0],&d_[0],&du_[0],&du2_[0],&ipiv_[0],&(*up)[0],n_,&info_);
        }

        //! Same as above but with overwrite in place
        void solve(ROL::Ptr<std::vector<Real> > up) {            
            lapack_.GTTRS('n',n_,1,&dl_[0],&d_[0],&du_[0],&du2_[0],&ipiv_[0],&(*up)[0],n_,&info_);
        }
 
        //! Given u, compute f = -u'' 
        void apply(ROL::Ptr<const std::vector<Real> > up,  ROL::Ptr<std::vector<Real> > fp) {
            (*fp)[0] = (2.0*(*up)[0]-(*up)[1])/dx2_;
  
            for(int i=1;i<n_-1;++i) {
                (*fp)[i] = (2.0*(*up)[i]-(*up)[i-1]-(*up)[i+1])/dx2_;
            } 
            (*fp)[n_-1] = (2.0*(*up)[n_-1]-(*up)[n_-2])/dx2_;
        }    

        //! Same as above but with overwrite in place
        void apply(ROL::Ptr<std::vector<Real> > fp) {

            ROL::Ptr<std::vector<Real> > up = ROL::makePtr<std::vector<Real>>(n_, 0.0);
             for(int i=0;i<n_;++i) {
                (*up)[i] = (*fp)[i];
            }

            (*fp)[0] = (2.0*(*up)[0]-(*up)[1])/dx2_;
            for(int i=1;i<n_-1;++i) {
                (*fp)[i] = (2.0*(*up)[i]-(*up)[i-1]-(*up)[i+1])/dx2_;
            } 
            (*fp)[n_-1] = (2.0*(*up)[n_-1]-(*up)[n_-2])/dx2_;

        }
};



