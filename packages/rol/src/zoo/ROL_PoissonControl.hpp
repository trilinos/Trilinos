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
    \brief  Contains definitions for Poisson optimal control.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_POISSONCONTROL_HPP
#define ROL_POISSONCONTROL_HPP

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"

namespace ROL {

  /** \brief Poisson distributed control.
   */
  template<class Real>
  class Objective_PoissonControl : public Objective<Real> {
  private:
    Real alpha_;

  public:

    Objective_PoissonControl(Real alpha = 1.e-4) : alpha_(alpha) {}

    void apply_mass(Vector<Real> &Mz, const Vector<Real> &z ) {
      Teuchos::RCP<const std::vector<Real> > zp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(z))).getVector();
      Teuchos::RCP<std::vector<Real> > Mzp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(Mz)).getVector());

      int  n = zp->size();
      Real h = 1.0/((Real)n+1.0);
      for (int i=0; i<n; i++) {
        if ( i == 0 ) {
          (*Mzp)[i] = h/6.0*(4.0*(*zp)[i] + (*zp)[i+1]);
        }
        else if ( i == n-1 ) {
          (*Mzp)[i] = h/6.0*((*zp)[i-1] + 4.0*(*zp)[i]);
        }
        else {
          (*Mzp)[i] = h/6.0*((*zp)[i-1] + 4.0*(*zp)[i] + (*zp)[i+1]);
        }
      }
    }

    void solve_poisson(Vector<Real> & u, const Vector<Real> & z) {
      Teuchos::RCP<std::vector<Real> > up =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(u)).getVector());

      int  n = up->size();
      Real h = 1.0/((Real)n+1.0);
      StdVector<Real> b( Teuchos::rcp( new std::vector<Real>(n,0.0) ) );
      Teuchos::RCP<std::vector<Real> > bp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(b)).getVector());
      this->apply_mass(b,z);

      Real d   =  2.0/h;
      Real o   = -1.0/h;
      Real m   = 0.0;
      std::vector<Real> c(n,o);
      c[0]     = c[0]/d;
      (*up)[0] = (*bp)[0]/d;
      for ( int i = 1; i < n; i++ ) {
        m        = 1.0/(d - o*c[i-1]);
        c[i]     = c[i]*m;
        (*up)[i] = ( (*bp)[i] - o*(*up)[i-1] )*m;
      }
      for ( int i = n-1; i > 0; i-- ) {
        (*up)[i-1] = (*up)[i-1] - c[i-1]*(*up)[i];
      }
    }

    Real evaluate_target(Real x) {
      Real val = 1.0/3.0*std::pow(x,4.0) - 2.0/3.0*std::pow(x,3.0) + 1.0/3.0*x + 8.0*this->alpha_;
      return val;
    }

    Real value( const Vector<Real> &z, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > zp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(z))).getVector();
      int  n    = zp->size();
      Real h    = 1.0/((Real)n+1.0);
      // SOLVE STATE EQUATION
      StdVector<Real> u( Teuchos::rcp( new std::vector<Real>(n,0.0) ) );
      this->solve_poisson(u,z);
      Teuchos::RCP<std::vector<Real> > up =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(u)).getVector());

      Real val  = 0.0;
      Real res  = 0.0;
      Real res1 = 0.0;
      Real res2 = 0.0;
      Real res3 = 0.0;
      for (int i=0; i<n; i++) {
        res = this->alpha_*(*zp)[i];
        if ( i == 0 ) {
          res *= h/6.0*(4.0*(*zp)[i] + (*zp)[i+1]);
          res1 = (*up)[i]-evaluate_target((Real)(i+1)*h);
          res2 = (*up)[i+1]-evaluate_target((Real)(i+2)*h);
          res += h/6.0*(4.0*res1 + res2)*res1;
        }
        else if ( i == n-1 ) {
          res *= h/6.0*((*zp)[i-1] + 4.0*(*zp)[i]);
          res1 = (*up)[i-1]-evaluate_target((Real)(i)*h);
          res2 = (*up)[i]-evaluate_target((Real)(i+1)*h);
          res += h/6.0*(res1 + 4.0*res2)*res2;
        }
        else {
          res *= h/6.0*((*zp)[i-1] + 4.0*(*zp)[i] + (*zp)[i+1]);
          res1 = (*up)[i-1]-evaluate_target((Real)(i)*h);
          res2 = (*up)[i]-evaluate_target((Real)(i+1)*h);
          res3 = (*up)[i+1]-evaluate_target((Real)(i+2)*h);
          res += h/6.0*(res1 + 4.0*res2 + res3)*res2;
        }
        val += 0.5*res;
      }
      return val;
   }

    void gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > zp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(z))).getVector();
      Teuchos::RCP<std::vector<Real> > gp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(g)).getVector());
      int  n = zp->size();
      Real h = 1.0/((Real)n+1.0);

      // SOLVE STATE EQUATION
      StdVector<Real> u( Teuchos::rcp( new std::vector<Real>(n,0.0) ) );
      this->solve_poisson(u,z);
      Teuchos::RCP<std::vector<Real> > up =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(u)).getVector());

      // SOLVE ADJOINT EQUATION
      StdVector<Real> res( Teuchos::rcp( new std::vector<Real>(n,0.0) ) );
      Teuchos::RCP<std::vector<Real> > rp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(res)).getVector());
      for (int i=0; i<n; i++) {
        (*rp)[i] = -((*up)[i]-evaluate_target((Real)(i+1)*h));
      }
      StdVector<Real> p( Teuchos::rcp( new std::vector<Real>(n,0.0) ) );
      this->solve_poisson(p,res);
      Teuchos::RCP<std::vector<Real> > pp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(p)).getVector());

      Real res1 = 0.0;
      Real res2 = 0.0;
      Real res3 = 0.0;
      for (int i=0; i<n; i++) {
        if ( i == 0 ) {
          res1 = this->alpha_*(*zp)[i] - (*pp)[i];
          res2 = this->alpha_*(*zp)[i+1] - (*pp)[i+1];
          (*gp)[i] = h/6.0*(4.0*res1 + res2);
        }
        else if ( i == n-1 ) {
          res1 = this->alpha_*(*zp)[i-1] - (*pp)[i-1];
          res2 = this->alpha_*(*zp)[i] - (*pp)[i];
          (*gp)[i] = h/6.0*(res1 + 4.0*res2);
        }
        else {
          res1 = this->alpha_*(*zp)[i-1] - (*pp)[i-1];
          res2 = this->alpha_*(*zp)[i] - (*pp)[i];
          res3 = this->alpha_*(*zp)[i+1] - (*pp)[i+1];
          (*gp)[i] = h/6.0*(res1 + 4.0*res2 + res3);
        }
      }
    }
#if USE_HESSVEC
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > zp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(z))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());

      int  n = zp->size();
      Real h = 1.0/((Real)n+1.0);

      // SOLVE STATE EQUATION
      StdVector<Real> u( Teuchos::rcp( new std::vector<Real>(n,0.0) ) );
      this->solve_poisson(u,v);
      Teuchos::RCP<std::vector<Real> > up =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(u)).getVector());

      // SOLVE ADJOINT EQUATION
      StdVector<Real> p( Teuchos::rcp( new std::vector<Real>(n,0.0) ) );
      this->solve_poisson(p,u);
      Teuchos::RCP<std::vector<Real> > pp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(p)).getVector());

      Real res1 = 0.0;
      Real res2 = 0.0;
      Real res3 = 0.0;
      for (int i=0; i<n; i++) {
        if ( i == 0 ) {
          res1 = this->alpha_*(*vp)[i] + (*pp)[i];
          res2 = this->alpha_*(*vp)[i+1] + (*pp)[i+1];
          (*hvp)[i] = h/6.0*(4.0*res1 + res2);
        }
        else if ( i == n-1 ) {
          res1 = this->alpha_*(*vp)[i-1] + (*pp)[i-1];
          res2 = this->alpha_*(*vp)[i] + (*pp)[i];
          (*hvp)[i] = h/6.0*(res1 + 4.0*res2);
        }
        else {
          res1 = this->alpha_*(*vp)[i-1] + (*pp)[i-1];
          res2 = this->alpha_*(*vp)[i] + (*pp)[i];
          res3 = this->alpha_*(*vp)[i+1] + (*pp)[i+1];
          (*hvp)[i] = h/6.0*(res1 + 4.0*res2 + res3);
        }
      }
    }
#endif
  };

  template<class Real>
  void getPoissonControl( Teuchos::RCP<Objective<Real> > &obj, Vector<Real> &x0, Vector<Real> &x ) {
    // Cast Initial Guess and Solution Vectors
    Teuchos::RCP<std::vector<Real> > x0p =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(x0)).getVector());
    Teuchos::RCP<std::vector<Real> > xp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(x)).getVector());
    int n = xp->size();
    // Resize Vectors
    n = 512;
    x0p->resize(n);
    xp->resize(n);
    // Instantiate Objective Function
    obj = Teuchos::rcp( new Objective_PoissonControl<Real> );
    // Get Initial Guess
    for (int i=0; i<n; i++) {
      (*x0p)[i] = 0.0;
    }
    // Get Solution
    Real h  = 1.0/((Real)n+1.0);
    Real pt = 0.0;
    for( int i=0; i<n; i++ ) {
      pt = (Real)(i+1)*h;
      (*xp)[i] = 4.0*pt*(1.0-pt);
    }
  }

}// End ROL Namespace

#endif
