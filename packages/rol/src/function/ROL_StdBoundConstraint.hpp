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

/** \file
    \brief  Contains definitions for std::vector bound constraints.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef ROL_STDBOUNDCONSTRAINT_HPP
#define ROL_STDBOUNDCONSTRAINT_HPP

#include "ROL_StdVector.hpp"
#include "ROL_BoundConstraint.hpp"

namespace ROL {

  template<class Real>
  class StdBoundConstraint : public BoundConstraint<Real> {
  private: 
    int dim_;
    std::vector<Real> x_lo_;
    std::vector<Real> x_up_;
    Real min_diff_;
  public:
    StdBoundConstraint(std::vector<Real> &l, std::vector<Real> &u) : x_lo_(l), x_up_(u) {
      dim_ = x_lo_.size();
      for ( int i = 0; i < dim_; i++ ) { 
        if ( i == 0 ) {
          min_diff_ = x_up_[i] - x_lo_[i];
        }
        else {
          min_diff_ = ( (min_diff_ < (x_up_[i] - x_lo_[i])) ? min_diff_ : (x_up_[i] - x_lo_[i]) );
        }
      }
      min_diff_ *= 0.5;
    }

    bool isFeasible( const Vector<Real> &x ) {
      Teuchos::RCP<const std::vector<Real> > ex =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      bool val = true;
      int  cnt = 1;                                                     
      for ( int i = 0; i < this->dim_; i++ ) {
        if ( (*ex)[i] >= this->x_lo_[i] && (*ex)[i] <= this->x_up_[i] ) { cnt *= 1; }
        else                                                            { cnt *= 0; }
      } 
      if ( cnt == 0 ) { val = false; }
      return val; 
    }  

    void project( Vector<Real> &x ) {
      Teuchos::RCP<std::vector<Real> > ex =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(x)).getVector());
      for ( int i = 0; i < this->dim_; i++ ) {
        (*ex)[i] = std::max(this->x_lo_[i],std::min(this->x_up_[i],(*ex)[i]));
      }
    }

    void pruneLowerActive(Vector<Real> &v, const Vector<Real> &x, Real eps) {
      Teuchos::RCP<const std::vector<Real> > ex = 
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<std::vector<Real> > ev =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(v)).getVector());
      Real epsn = std::min(eps,this->min_diff_);
      for ( int i = 0; i < this->dim_; i++ ) {
        if ( ((*ex)[i] <= this->x_lo_[i]+epsn) ) {
          (*ev)[i] = 0.0;
        }
      }
    }           

    void pruneUpperActive(Vector<Real> &v, const Vector<Real> &x, Real eps) {
      Teuchos::RCP<const std::vector<Real> > ex = 
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<std::vector<Real> > ev =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(v)).getVector());
      Real epsn = std::min(eps,this->min_diff_);
      for ( int i = 0; i < this->dim_; i++ ) {
        if ( ((*ex)[i] >= this->x_up_[i]-epsn) ) {
          (*ev)[i] = 0.0;
        }
      }
    }           

    void pruneActive(Vector<Real> &v, const Vector<Real> &x, Real eps) {
      Teuchos::RCP<const std::vector<Real> > ex = 
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<std::vector<Real> > ev =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(v)).getVector());
      Real epsn = std::min(eps,this->min_diff_);
      for ( int i = 0; i < this->dim_; i++ ) {
        if ( ((*ex)[i] <= this->x_lo_[i]+epsn) || 
             ((*ex)[i] >= this->x_up_[i]-epsn) ) {
          (*ev)[i] = 0.0;
        }
      }
    }           

    void pruneLowerActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps) {
      Teuchos::RCP<const std::vector<Real> > ex = 
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > eg =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(g))).getVector();
      Teuchos::RCP<std::vector<Real> > ev =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(v)).getVector());
      Real epsn = std::min(eps,this->min_diff_);
      for ( int i = 0; i < this->dim_; i++ ) {
        if ( ((*ex)[i] <= this->x_lo_[i]+epsn && (*eg)[i] > 0.0) ) {
          (*ev)[i] = 0.0;
        }
      }
    }           

    void pruneUpperActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps) {
      Teuchos::RCP<const std::vector<Real> > ex = 
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > eg =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(g))).getVector();
      Teuchos::RCP<std::vector<Real> > ev =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(v)).getVector());
      Real epsn = std::min(eps,this->min_diff_);
      for ( int i = 0; i < this->dim_; i++ ) {
        if ( ((*ex)[i] >= this->x_up_[i]-epsn && (*eg)[i] < 0.0) ) {
          (*ev)[i] = 0.0;
        }
      }
    }           

    void pruneActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps) {
      Teuchos::RCP<const std::vector<Real> > ex = 
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > eg =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(g))).getVector();
      Teuchos::RCP<std::vector<Real> > ev =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(v)).getVector());
      Real epsn = std::min(eps,this->min_diff_);
      for ( int i = 0; i < this->dim_; i++ ) {
        if ( ((*ex)[i] <= this->x_lo_[i]+epsn && (*eg)[i] > 0.0) || 
             ((*ex)[i] >= this->x_up_[i]-epsn && (*eg)[i] < 0.0) ) {
          (*ev)[i] = 0.0;
        }
      }
    }        

    void setVectorToUpperBound( ROL::Vector<Real> &u ) {
      Teuchos::RCP<std::vector<Real> > us = Teuchos::rcp( new std::vector<Real>(this->dim_,0.0) );
      us->assign(this->x_up_.begin(),this->x_up_.end());
      Teuchos::RCP<ROL::Vector<Real> > up = Teuchos::rcp( new ROL::StdVector<Real>(us) );
      u.set(*up);
    }

    void setVectorToLowerBound( ROL::Vector<Real> &l ) {
      Teuchos::RCP<std::vector<Real> > ls = Teuchos::rcp( new std::vector<Real>(this->dim_,0.0) );
      ls->assign(this->x_lo_.begin(),this->x_lo_.end());
      Teuchos::RCP<ROL::Vector<Real> > lp = Teuchos::rcp( new ROL::StdVector<Real>(ls) );
      l.set(*lp);
    }   
  };  

}// End ROL Namespace

#endif
