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
    \brief  Contains definitions for Poisson material inversion.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_POISSONINVERSION_HPP
#define ROL_POISSONINVERSION_HPP

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_HelperFunctions.hpp"

#include "Teuchos_LAPACK.hpp"

namespace ROL {
namespace ZOO {

  /** \brief Poisson material inversion.
   */
  template<class Real>
  class Objective_PoissonInversion : public Objective<Real> {
  private:
    int nu_;
    int nz_;

    Real hu_;
    Real hz_;

    Real alpha_;

    Real eps_;
    int  reg_type_;

  public:

    /* CONSTRUCTOR */
    Objective_PoissonInversion(int nz = 32, Real alpha = 1.e-4) : nz_(nz), alpha_(alpha) {
      nu_       = nz_-1;
      hu_       = 1.0/((Real)nu_+1.0);
      hz_       = hu_; 
      eps_      = 1.e-4;
      reg_type_ = 2;
    }

    /* REGULARIZATION DEFINITIONS */
    Real reg_value(const Vector<Real> &z) {
      Teuchos::RCP<const std::vector<Real> > zp = ROL::StdVector_Helper::constDownCast(z);
//      Teuchos::RCP<const std::vector<Real> > zp =
//        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(z))).getVector();
      Real val = 0.0;
      for (int i = 0; i < this->nz_; i++) {
        if ( this->reg_type_ == 2 ) {
          val += this->alpha_/2.0 * this->hz_ * (*zp)[i]*(*zp)[i];
        }
        else if ( this->reg_type_ == 1 ) {
          val += this->alpha_ * this->hz_ * std::sqrt((*zp)[i]*(*zp)[i] + this->eps_);
        }
        else if ( this->reg_type_ == 0 ) {
          if ( i < this->nz_-1 ) {
            val += this->alpha_ * std::sqrt(std::pow((*zp)[i]-(*zp)[i+1],2.0)+this->eps_);
          }
        }
      }
      return val;
    }

    void reg_gradient(Vector<Real> &g, const Vector<Real> &z) {
      if ( this->reg_type_ == 2 ) {
        g.set(z);
        g.scale(this->alpha_*this->hz_);    
      } 
      else if ( this->reg_type_ == 1 ) {
        Teuchos::RCP<const std::vector<Real> > zp = ROL::StdVector_Helper::constDownCast(z);
        Teuchos::RCP<std::vector<Real> >       gp = ROL::StdVector_Helper::downCast(g);
//        Teuchos::RCP<const std::vector<Real> > zp =
//          (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(z))).getVector();
//        Teuchos::RCP<std::vector<Real> > gp =
//          Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(g)).getVector());
        for (int i = 0; i < this->nz_; i++) {
          (*gp)[i] = this->alpha_ * this->hz_ * (*zp)[i]/std::sqrt(std::pow((*zp)[i],2.0)+this->eps_);
        }
      }
      else if ( this->reg_type_ == 0 ) {
        Teuchos::RCP<const std::vector<Real> > zp = ROL::StdVector_Helper::constDownCast(z);
        Teuchos::RCP<std::vector<Real> >       gp = ROL::StdVector_Helper::downCast(g);
//        Teuchos::RCP<const std::vector<Real> > zp =
//          (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(z))).getVector();
//        Teuchos::RCP<std::vector<Real> > gp =
//          Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(g)).getVector());
        Real diff = 0.0;
        for (int i = 0; i < this->nz_; i++) {
          if ( i == 0 ) {
            diff     = (*zp)[i]-(*zp)[i+1];
            (*gp)[i] = this->alpha_ * diff/std::sqrt(std::pow(diff,2.0)+this->eps_);
          }
          else if ( i == this->nz_-1 ) {
            diff     = (*zp)[i-1]-(*zp)[i];
            (*gp)[i] = -this->alpha_ * diff/std::sqrt(std::pow(diff,2.0)+this->eps_);
          }
          else {
            diff      = (*zp)[i]-(*zp)[i+1];
            (*gp)[i]  = this->alpha_ * diff/std::sqrt(std::pow(diff,2.0)+this->eps_);
            diff      = (*zp)[i-1]-(*zp)[i];
            (*gp)[i] -= this->alpha_ * diff/std::sqrt(std::pow(diff,2.0)+this->eps_);
          }
        }
      }
    }

    void reg_hessVec(Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z) {
      if ( this->reg_type_ == 2 ) {
        hv.set(v);
        hv.scale(this->alpha_*this->hz_);
      }
      else if ( this->reg_type_ == 1 ) {
        Teuchos::RCP<const std::vector<Real> > zp  = ROL::StdVector_Helper::constDownCast(z);
        Teuchos::RCP<const std::vector<Real> > vp  = ROL::StdVector_Helper::constDownCast(v);
        Teuchos::RCP<std::vector<Real> >       hvp = ROL::StdVector_Helper::downCast(hv);
//        Teuchos::RCP<const std::vector<Real> > zp =
//          (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(z))).getVector();
//        Teuchos::RCP<const std::vector<Real> > vp =
//          (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
//        Teuchos::RCP<std::vector<Real> > hvp =
//          Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());
        for (int i = 0; i < this->nz_; i++) {
          (*hvp)[i] = this->alpha_*this->hz_*(*vp)[i]*this->eps_/std::pow(std::pow((*zp)[i],2.0)+this->eps_,3.0/2.0);
        }
      }
      else if ( this->reg_type_ == 0 ) {
        Teuchos::RCP<const std::vector<Real> > zp  = ROL::StdVector_Helper::constDownCast(z);
        Teuchos::RCP<const std::vector<Real> > vp  = ROL::StdVector_Helper::constDownCast(v);
        Teuchos::RCP<std::vector<Real> >       hvp = ROL::StdVector_Helper::downCast(hv);
//        Teuchos::RCP<const std::vector<Real> > zp =
//          (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(z))).getVector();
//        Teuchos::RCP<const std::vector<Real> > vp =
//          (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
//        Teuchos::RCP<std::vector<Real> > hvp =
//          Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());
        Real diff1 = 0.0;
        Real diff2 = 0.0;
        for (int i = 0; i < this->nz_; i++) {
          if ( i == 0 ) {
            diff1 = (*zp)[i]-(*zp)[i+1];
            diff1 = this->eps_/std::pow(std::pow(diff1,2.0)+this->eps_,3.0/2.0);
            (*hvp)[i] = this->alpha_* ((*vp)[i]*diff1 - (*vp)[i+1]*diff1);
          }
          else if ( i == this->nz_-1 ) {
            diff2 = (*zp)[i-1]-(*zp)[i];
            diff2 = this->eps_/std::pow(std::pow(diff2,2.0)+this->eps_,3.0/2.0);
            (*hvp)[i] = this->alpha_* (-(*vp)[i-1]*diff2 + (*vp)[i]*diff2);
          }
          else {
            diff1 = (*zp)[i]-(*zp)[i+1];
            diff1 = this->eps_/std::pow(std::pow(diff1,2.0)+this->eps_,3.0/2.0);
            diff2 = (*zp)[i-1]-(*zp)[i];
            diff2 = this->eps_/std::pow(std::pow(diff2,2.0)+this->eps_,3.0/2.0);
            (*hvp)[i] = this->alpha_* (-(*vp)[i-1]*diff2 + (*vp)[i]*(diff1 + diff2) - (*vp)[i+1]*diff1);
          }
        }
      }
    }

    /* FINITE ELEMENT DEFINTIONS */
    void apply_mass(Vector<Real> &Mf, const Vector<Real> &f ) {
      Teuchos::RCP<const std::vector<Real> > fp  = ROL::StdVector_Helper::constDownCast(f);
      Teuchos::RCP<std::vector<Real> >       Mfp = ROL::StdVector_Helper::downCast(Mf);
//      Teuchos::RCP<const std::vector<Real> > fp =
//        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(f))).getVector();
//      Teuchos::RCP<std::vector<Real> > Mfp =
//        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(Mf)).getVector());

      for (int i = 0; i < this->nu_; i++) {
        if ( i == 0 ) {
          (*Mfp)[i] = this->hu_/6.0*(4.0*(*fp)[i] + (*fp)[i+1]);
        }
        else if ( i == this->nu_-1 ) {
          (*Mfp)[i] = this->hu_/6.0*((*fp)[i-1] + 4.0*(*fp)[i]);
        }
        else {
          (*Mfp)[i] = this->hu_/6.0*((*fp)[i-1] + 4.0*(*fp)[i] + (*fp)[i+1]);
        }
      }
    }

    void solve_poisson(Vector<Real> &u, const Vector<Real> &z, Vector<Real> &b) {
      Teuchos::RCP<const std::vector<Real> > zp = ROL::StdVector_Helper::constDownCast(z);
      Teuchos::RCP<std::vector<Real> >       up = ROL::StdVector_Helper::downCast(u);
      Teuchos::RCP<std::vector<Real> >       bp = ROL::StdVector_Helper::downCast(b);
//      Teuchos::RCP<const std::vector<Real> > zp =
//        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(z))).getVector();
//      Teuchos::RCP<std::vector<Real> > up =
//        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(u)).getVector());
//      Teuchos::RCP<std::vector<Real> > bp =
//        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(b)).getVector());

      // Get Diagonal and Off-Diagonal Entries of PDE Jacobian
      std::vector<Real> d(this->nu_,1.0);
      std::vector<Real> o(this->nu_-1,1.0);
      for ( int i = 0; i < this->nu_; i++ ) {
        d[i] = (std::exp((*zp)[i]) + std::exp((*zp)[i+1]))/this->hu_;
        if ( i < this->nu_-1 ) {
          o[i] *= -std::exp((*zp)[i+1])/this->hu_;
        }
      }

      // Solve Tridiagonal System Using LAPACK's SPD Tridiagonal Solver
      Teuchos::LAPACK<int,Real> lp;
      int info;
      int ldb  = this->nu_;
      int nhrs = 1;
      lp.PTTRF(this->nu_,&d[0],&o[0],&info);
      lp.PTTRS(this->nu_,nhrs,&d[0],&o[0],&(*bp)[0],ldb,&info);
      u.set(b);
    }

    Real evaluate_target(Real x) {
      return x*(1.0-x);
    }

    void apply_linearized_control_operator( Vector<Real> &Bd, const Vector<Real> &z, 
                                      const Vector<Real> &d,  const Vector<Real> &u ) {
      Teuchos::RCP<const std::vector<Real> > zp  = ROL::StdVector_Helper::constDownCast(z);
      Teuchos::RCP<const std::vector<Real> > up  = ROL::StdVector_Helper::constDownCast(u);
      Teuchos::RCP<const std::vector<Real> > dp  = ROL::StdVector_Helper::constDownCast(d);
      Teuchos::RCP<std::vector<Real> >       Bdp = ROL::StdVector_Helper::downCast(Bd);
//      Teuchos::RCP<const std::vector<Real> > zp =
//        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(z))).getVector();
//      Teuchos::RCP<const std::vector<Real> > up =
//        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(u))).getVector();
//      Teuchos::RCP<const std::vector<Real> > dp =
//        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(d))).getVector();
//      Teuchos::RCP<std::vector<Real> > Bdp =
//        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(Bd)).getVector());

      for (int i = 0; i < this->nu_; i++) {
        if ( i == 0 ) {
          (*Bdp)[i] = 1.0/this->hu_*( std::exp((*zp)[i])*(*up)[i]*(*dp)[i] 
                                    + std::exp((*zp)[i+1])*((*up)[i]-(*up)[i+1])*(*dp)[i+1] );
        }
        else if ( i == this->nu_-1 ) {
          (*Bdp)[i] = 1.0/this->hu_*( std::exp((*zp)[i])*((*up)[i]-(*up)[i-1])*(*dp)[i] 
                                    + std::exp((*zp)[i+1])*(*up)[i]*(*dp)[i+1] );
        }
        else {
          (*Bdp)[i] = 1.0/this->hu_*( std::exp((*zp)[i])*((*up)[i]-(*up)[i-1])*(*dp)[i] 
                                    + std::exp((*zp)[i+1])*((*up)[i]-(*up)[i+1])*(*dp)[i+1] );
        }
      }
    }

    void apply_transposed_linearized_control_operator( Vector<Real> &Bd, const Vector<Real> &z,
                                                 const Vector<Real> &d,  const Vector<Real> &u ) {
      Teuchos::RCP<const std::vector<Real> > zp  = ROL::StdVector_Helper::constDownCast(z);
      Teuchos::RCP<const std::vector<Real> > up  = ROL::StdVector_Helper::constDownCast(u);
      Teuchos::RCP<const std::vector<Real> > dp  = ROL::StdVector_Helper::constDownCast(d);
      Teuchos::RCP<std::vector<Real> >       Bdp = ROL::StdVector_Helper::downCast(Bd);
//      Teuchos::RCP<const std::vector<Real> > zp =
//        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(z))).getVector();
//      Teuchos::RCP<const std::vector<Real> > up =
//        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(u))).getVector();
//      Teuchos::RCP<const std::vector<Real> > dp =
//        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(d))).getVector();
//      Teuchos::RCP<std::vector<Real> > Bdp =
//        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(Bd)).getVector());

      for (int i = 0; i < this->nz_; i++) {
        if ( i == 0 ) {
          (*Bdp)[i] = std::exp((*zp)[i])/this->hu_*(*up)[i]*(*dp)[i];
        }
        else if ( i == this->nz_-1 ) {
          (*Bdp)[i] = std::exp((*zp)[i])/this->hu_*(*up)[i-1]*(*dp)[i-1];
        }
        else {
          (*Bdp)[i] = std::exp((*zp)[i])/this->hu_*( ((*up)[i]-(*up)[i-1])*((*dp)[i]-(*dp)[i-1]) );
        }
      }
    }
    
    void apply_transposed_linearized_control_operator_2( Vector<Real> &Bd, const Vector<Real> &z, const Vector<Real> &v,
                                                   const Vector<Real> &d,  const Vector<Real> &u ) {
      Teuchos::RCP<const std::vector<Real> > zp  = ROL::StdVector_Helper::constDownCast(z);
      Teuchos::RCP<const std::vector<Real> > vp  = ROL::StdVector_Helper::constDownCast(v);
      Teuchos::RCP<const std::vector<Real> > up  = ROL::StdVector_Helper::constDownCast(u);
      Teuchos::RCP<const std::vector<Real> > dp  = ROL::StdVector_Helper::constDownCast(d);
      Teuchos::RCP<std::vector<Real> >       Bdp = ROL::StdVector_Helper::downCast(Bd);
//      Teuchos::RCP<const std::vector<Real> > zp =
//        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(z))).getVector();
//      Teuchos::RCP<const std::vector<Real> > vp =
//        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
//      Teuchos::RCP<const std::vector<Real> > up =
//        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(u))).getVector();
//      Teuchos::RCP<const std::vector<Real> > dp =
//        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(d))).getVector();
//      Teuchos::RCP<std::vector<Real> > Bdp =
//        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(Bd)).getVector());

      for (int i = 0; i < this->nz_; i++) {
        if ( i == 0 ) {
          (*Bdp)[i] = (*vp)[i]*std::exp((*zp)[i])/this->hu_*(*up)[i]*(*dp)[i];
        }
        else if ( i == this->nz_-1 ) {
          (*Bdp)[i] = (*vp)[i]*std::exp((*zp)[i])/this->hu_*(*up)[i-1]*(*dp)[i-1];
        }
        else {
          (*Bdp)[i] = (*vp)[i]*std::exp((*zp)[i])/this->hu_*( ((*up)[i]-(*up)[i-1])*((*dp)[i]-(*dp)[i-1]) );
        }
      }
    }

    /* STATE AND ADJOINT EQUATION DEFINTIONS */
    void solve_state_equation(Vector<Real> &u, const Vector<Real> &z) {
      Real k1 = 1.0;
      Real k2 = 2.0;
      // Right Hand Side
      Teuchos::RCP<std::vector<Real> > bp = Teuchos::rcp( new std::vector<Real> (this->nu_, 0.0) );
      for ( int i = 0; i < this->nu_; i++ ) {
        if ( (Real)(i+1)*this->hu_ < 0.5 ) {
         (*bp)[i] = 2.0*k1*this->hu_;
        }
        else if ( std::abs((Real)(i+1)*this->hu_ - 0.5) < ROL_EPSILON ) {
         (*bp)[i] = (k1+k2)*this->hu_;
        }
        else if ( (Real)(i+1)*this->hu_ > 0.5 ) {
         (*bp)[i] = 2.0*k2*this->hu_;
        }
      }
     
      StdVector<Real> b(bp);
      // Solve Equation
      this->solve_poisson(u,z,b);
    }

    void solve_adjoint_equation(Vector<Real> &p, const Vector<Real> &u, const Vector<Real> &z) {
      Teuchos::RCP<const std::vector<Real> > up = ROL::StdVector_Helper::constDownCast(u);
      StdVector<Real> res( Teuchos::rcp( new std::vector<Real>(this->nu_,0.0) ) );
      Teuchos::RCP<std::vector<Real> >       rp = ROL::StdVector_Helper::downCast(res);
//      Teuchos::RCP<const std::vector<Real> > up =
//        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(u))).getVector();
//      StdVector<Real> res( Teuchos::rcp( new std::vector<Real>(this->nu_,0.0) ) );
//      Teuchos::RCP<std::vector<Real> > rp =
//        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(res)).getVector());
      for (int i = 0; i < this->nu_; i++) {
        (*rp)[i] = -((*up)[i]-this->evaluate_target((Real)(i+1)*this->hu_));
      }
      StdVector<Real> Mres( Teuchos::rcp( new std::vector<Real>(this->nu_,0.0) ) );
      this->apply_mass(Mres,res);
      this->solve_poisson(p,z,Mres);
    }

    void solve_state_sensitivity_equation(Vector<Real> &w, const Vector<Real> &v, 
                                          const Vector<Real> &u, const Vector<Real> &z) {
      StdVector<Real> b( Teuchos::rcp( new std::vector<Real>(this->nu_,0.0) ) );
      this->apply_linearized_control_operator(b,z,v,u);
      this->solve_poisson(w,z,b);
    }

    void solve_adjoint_sensitivity_equation(Vector<Real> &q, const Vector<Real> &w, const Vector<Real> &v,
                                            const Vector<Real> &p, const Vector<Real> &u, const Vector<Real> &z) {
      StdVector<Real> res( Teuchos::rcp( new std::vector<Real>(this->nu_,0.0) ) );
      this->apply_mass(res,w);
      StdVector<Real> res1( Teuchos::rcp( new std::vector<Real>(this->nu_,0.0) ) );
      this->apply_linearized_control_operator(res1,z,v,p);
      res.axpy(-1.0,res1);
      this->solve_poisson(q,z,res);
    }

    /* OBJECTIVE FUNCTION DEFINITIONS */
    Real value( const Vector<Real> &z, Real &tol ) {
      // SOLVE STATE EQUATION
      StdVector<Real> u( Teuchos::rcp( new std::vector<Real>(this->nu_,0.0) ) );
      Teuchos::RCP<std::vector<Real> > up = ROL::StdVector_Helper::downCast(u);
//      Teuchos::RCP<std::vector<Real> > up =
//        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(u)).getVector());
      this->solve_state_equation(u,z);

      // COMPUTE MISFIT
      StdVector<Real> res( Teuchos::rcp( new std::vector<Real>(this->nu_,0.0) ) );
      Teuchos::RCP<std::vector<Real> > rp = ROL::StdVector_Helper::downCast(res);
//      Teuchos::RCP<std::vector<Real> > rp =
//        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(res)).getVector());
      for (int i = 0; i < this->nu_; i++) {
        (*rp)[i] = ((*up)[i]-this->evaluate_target((Real)(i+1)*this->hu_));
      }
      StdVector<Real> Mres( Teuchos::rcp( new std::vector<Real>(this->nu_,0.0) ) );
      this->apply_mass(Mres,res);
      return 0.5*Mres.dot(res) + this->reg_value(z);
    } 

    void gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
      // SOLVE STATE EQUATION
      StdVector<Real> u( Teuchos::rcp( new std::vector<Real>(this->nu_,0.0) ) );
      this->solve_state_equation(u,z);

      // SOLVE ADJOINT EQUATION
      StdVector<Real> p( Teuchos::rcp( new std::vector<Real>(this->nu_,0.0) ) );
      this->solve_adjoint_equation(p,u,z);

      // Apply Transpose of Linearized Control Operator
      this->apply_transposed_linearized_control_operator(g,z,p,u);
     
      // Regularization gradient
      StdVector<Real> g_reg( Teuchos::rcp( new std::vector<Real>(this->nz_,0.0) ) );
      this->reg_gradient(g_reg,z); 

      // Build Gradient
      g.plus(g_reg);
    }
#if USE_HESSVEC
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
      // SOLVE STATE EQUATION
      StdVector<Real> u( Teuchos::rcp( new std::vector<Real>(this->nu_,0.0) ) );
      this->solve_state_equation(u,z);

      // SOLVE ADJOINT EQUATION
      StdVector<Real> p( Teuchos::rcp( new std::vector<Real>(this->nu_,0.0) ) );
      this->solve_adjoint_equation(p,u,z);

      // SOLVE STATE SENSITIVITY EQUATION
      StdVector<Real> w( Teuchos::rcp( new std::vector<Real>(this->nu_,0.0) ) );
      this->solve_state_sensitivity_equation(w,v,u,z);

      // SOLVE ADJOINT SENSITIVITY EQUATION
      StdVector<Real> q( Teuchos::rcp( new std::vector<Real>(this->nu_,0.0) ) );
      this->solve_adjoint_sensitivity_equation(q,w,v,p,u,z);

      // Apply Transpose of Linearized Control Operator
      this->apply_transposed_linearized_control_operator(hv,z,q,u);
    
      // Apply Transpose of Linearized Control Operator
      StdVector<Real> tmp( Teuchos::rcp( new std::vector<Real>(this->nz_,0.0) ) );
      this->apply_transposed_linearized_control_operator(tmp,z,w,p);
      hv.axpy(-1.0,tmp); 

      // Apply Transpose of 2nd Derivative of Control Operator
      tmp.zero();
      this->apply_transposed_linearized_control_operator_2(tmp,z,v,p,u);
      hv.plus(tmp);

      // Regularization hessVec
      StdVector<Real> hv_reg( Teuchos::rcp( new std::vector<Real>(this->nz_,0.0) ) );
      this->reg_hessVec(hv_reg,v,z);

      // Build hessVec
      hv.plus(hv_reg);
    }
#endif

    void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

      // Cast hv and v vectors to std::vector.
      Teuchos::RCP<std::vector<Real> > hvp = ROL::StdVector_Helper::downCast(hv);
      Teuchos::RCP<std::vector<Real> >  vp = ROL::StdVector_Helper::downCast(v);
//      Teuchos::RCP<std::vector<Real> > hvp =
//        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());
//      Teuchos::RCP<std::vector<Real> > vp =
//        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector());

      int dim = vp->size();

      // Compute dense Hessian.
      Teuchos::SerialDenseMatrix<int, Real> H(dim, dim);
      Objective_PoissonInversion<Real> & obj = *this; 
      H = computeDenseHessian<Real>(obj, x);

      // Compute eigenvalues, sort real part.
      std::vector<std::vector<Real> > eigenvals = computeEigenvalues<Real>(H);
      std::sort((eigenvals[0]).begin(), (eigenvals[0]).end());

      // Perform 'inertia' correction.
      Real inertia = (eigenvals[0])[0];
      Real correction = 0.0;
      if ( inertia <= 0.0 ) {
        correction = 2.0*std::abs(inertia);
        if ( inertia == 0.0 ) {
          int cnt = 0;
          while ( eigenvals[0][cnt] == 0.0 ) {
            cnt++;
          }
          correction = 0.5*eigenvals[0][cnt];
          if ( cnt == dim-1 ) {
            correction = 1.0;
          }
        }
        for (int i=0; i<dim; i++) {
          H(i,i) += correction;
        }
      }

      // Compute dense inverse Hessian.
      Teuchos::SerialDenseMatrix<int, Real> invH = computeInverse<Real>(H);

      // Apply dense inverse Hessian.
      Teuchos::SerialDenseVector<int, Real> hv_teuchos(Teuchos::View, &((*hvp)[0]), dim);
      Teuchos::SerialDenseVector<int, Real> v_teuchos(Teuchos::View, &((*vp)[0]), dim);
      hv_teuchos.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, invH, v_teuchos, 0.0);
    }

  };

  template<class Real>
  void getPoissonInversion( Teuchos::RCP<Objective<Real> > &obj, Vector<Real> &x0, Vector<Real> &x ) {
    // Cast Initial Guess and Solution Vectors
    Teuchos::RCP<std::vector<Real> > x0p = ROL::StdVector_Helper::downCast(x0);
    Teuchos::RCP<std::vector<Real> >  xp = ROL::StdVector_Helper::downCast(x);
//    Teuchos::RCP<std::vector<Real> > x0p =
//      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(x0)).getVector());
//    Teuchos::RCP<std::vector<Real> > xp =
//      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(x)).getVector());
    int n = xp->size();
    // Resize Vectors
    n = 128;
    x0p->resize(n);
    xp->resize(n);
    // Instantiate Objective Function
    obj = Teuchos::rcp( new Objective_PoissonInversion<Real>(n,1.e-6) );
    // Get Initial Guess
    for (int i=0; i<n; i++) {
      (*x0p)[i] = 1.5;
    }
    // Get Solution
    Real h  = 1.0/((Real)n+1);
    Real pt = 0.0;
    Real k1 = 1.0;
    Real k2 = 2.0;
    for( int i=0; i<n; i++ ) {
      pt = (Real)(i+1)*h;
      if ( pt >= 0.0 && pt < 0.5 ) {
        (*xp)[i] = std::log(k1);
      }
      else if ( pt >= 0.5 && pt < 1.0 ) {
        (*xp)[i] = std::log(k2); 
      }
    }
  }

} // End ZOO Namespace
} // End ROL Namespace

#endif
