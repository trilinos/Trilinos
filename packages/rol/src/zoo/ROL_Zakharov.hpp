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
    \brief  Contains definitions for the Zakharov function.

    Objective function: 
    \f[f(\mathbf{x}) = \mathbf{x}^\top\mathbf{x} + \frac{1}{4}(\mathbf{k}^\top \mathbf{x})^2 +
                                                   \frac{1}{16}(\mathbf{k}^\top \mathbf{x})^4 \f]
    Where \f$\mathbf{k}=(1,\cdots,n)\f$
 
    Gradient:
    \f[
    g=\nabla f(\mathbf{x}) = 2\mathbf{x} + 
                           \frac{1}{4}\left(2(\mathbf{k}^\top\mathbf{x})+(\mathbf{k}^\top\mathbf{x})^3\right)\mathbf{k} 
    \f]

    Hessian: 
    \f[
    H=\nabla^2 f(\mathbf{x}) = 2 I + \frac{1}{4}[2+3(\mathbf{k}^\top\mathbf{x})^2]\mathbf{kk}^\top
    \f]
 
    The Hessian is a multiple of the identity plus a rank one symmetric 
    matrix, therefore the action of the inverse Hessian can be 
    performed using the Sherman-Morrison formula.

    \f[
    H^{-1}\mathbf{v} = \frac{1}{2}\mathbf{v}-\frac{(\mathbf{k}^\top\mathbf{v})}
                                             {\frac{16}{2+3(\mathbf{k}^\top\mathbf{x})^2}+2\mathbf{k^\top}\mathbf{k}}\mathbf{k}
    \f]

    \author Created by G. von Winckel
**/

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_ZAKHAROV_HPP
#define ROL_ZAKHAROV_HPP

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"

namespace ROL {
namespace ZOO {

  /** \brief Zakharov function.
   */
  template<class Real>
  class Objective_Zakharov : public Objective<Real> {
  private:
    

  public:
    Objective_Zakharov(void) {}

    Real value( const Vector<Real> &x, Real &tol ) {
      StdVector<Real> & ex =
        Teuchos::dyn_cast<StdVector<Real> >(const_cast <Vector<Real> &>(x));
      Teuchos::RCP<const std::vector<Real> > xp = ex.getVector();
      int n = xp->size();

      Real val = 0;
      Real xdotx = 0; 
      Real kdotx = 0; 

      for(int i=0; i<n; i++ ) {
          xdotx += pow((*xp)[i],2);
          kdotx += double(i+1)*((*xp)[i]);
      }
      val = xdotx + pow(kdotx,2)/4.0 + pow(kdotx,4)/16.0;
      return val;
    }

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<std::vector<Real> > gp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(g)).getVector());

      int n = xp->size();

      Real kdotx = 0;
      for( int i=0; i<n; i++ ) {
        kdotx += double(i+1)*((*xp)[i]);
      }

      Real coeff = (2.0*kdotx+pow(kdotx,3))/4.0;  

      for( int i=0; i<n; i++ ) {
        (*gp)[i]   =  2.0*((*xp)[i]) + double(i+1)*coeff;
      }
    }

#if USE_HESSVEC
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());

      int n = xp->size();

      Real kdotk = 0;
      Real kdotx = 0;
      Real kdotv = 0;

      for( int i=0; i<n; i++ ) {
        kdotx += double(i+1)*((*xp)[i]);
        kdotk += pow(double(i+1),2);
        kdotv += double(i+1)*((*vp)[i]); 
      }
      
      Real coeff = (2.0+3.0*pow(kdotx,2))*kdotv/4.0;

      for( int i=0; i<n; i++ ) {
        (*hvp)[i] = 2.0*(*vp)[i] + coeff*double(i+1);
      }

    }
#endif
    void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());

      int n = xp->size();
      Real kdotv = 0;
      Real kdotx = 0; 
      Real kdotk = 0; 

      for( int i=0; i<n; i++) {
          kdotv += double(i+1)*((*vp)[i]); 
          kdotx += double(i+1)*((*xp)[i]); 
          kdotk += pow(double(i+1),2);
      }

      Real coeff = -kdotv/(16.0/(2.0+3.0*pow(kdotx,2))+2.0*kdotk);

      for( int i=0; i<n; i++) {
          (*hvp)[i] = 0.5*((*vp)[i]) + coeff*(i+1); 
      }       

    }
  };








  template<class Real>
  void getZakharov( Teuchos::RCP<Objective<Real> > &obj, Vector<Real> &x0, Vector<Real> &x ) {
    // Cast Initial Guess and Solution Vectors
    Teuchos::RCP<std::vector<Real> > x0p =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(x0)).getVector());
    Teuchos::RCP<std::vector<Real> > xp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(x)).getVector());
    int n = xp->size();
    // Resize Vectors
    n = 10;
    x0p->resize(n);
    xp->resize(n);
    // Instantiate Objective Function
    obj = Teuchos::rcp( new Objective_Zakharov<Real> );
    // Get Initial Guess
    (*x0p)[0] =  3.0;
    (*x0p)[1] =  3.0;
    // Get Solution
    (*xp)[0] = 0.0;
    (*xp)[1] = 0.0;
  }


}// End ZOO Namespace
}// End ROL Namespace

#endif
