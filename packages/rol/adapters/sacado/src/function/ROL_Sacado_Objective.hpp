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

#ifndef ROL_SACADO_OBJECTIVE
#define ROL_SACADO_OBJECTIVE

#include "Sacado.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"

using namespace ROL;

/** \brief Generic objective wrapper class for class that uses Sacado */
template<class Real, template<class> class Obj>
class Sacado_Objective : public Objective<Real> {

     protected:

        Obj<Real> obj_;

    /* Evaluate the gradient at x */
    template<class ScalarT>
    void gradientAD( Vector<ScalarT> &g, const Vector<ScalarT> &x, Real &tol );

    /* Compute the action of the Hessian evaluated at x on a vector v */
    template<class ScalarT>
    void hessVecAD( Vector<ScalarT> &hv, const Vector<ScalarT> &v, const Vector<ScalarT> &x, Real &tol );


    public:

    Sacado_Objective() : obj_(Obj<Real>()) {}
    Sacado_Objective(const Obj<Real> &obj) : obj_(obj) {}

    /* Evaluate the objective function at x */
    Real value( const Vector<Real> &x, Real &tol ) {
      return obj_.value(x,tol);
    }

    /* Evaluate the gradient at x */
    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
        this->gradientAD(g,x,tol);
    }

    /* Compute the action of the Hessian evaluated at x on a vector v */
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
        this->hessVecAD(hv,v,x,tol);
    }
};



template<class Real, template<class> class Obj>
template<class ScalarT>
void Sacado_Objective<Real,Obj>::gradientAD(Vector<ScalarT> &g, const Vector<ScalarT> &x, Real &tol) {

    // Data type which supports automatic differentiation
    typedef Sacado::Fad::DFad<ScalarT> FadType;
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

    using Teuchos::RCP;       using Teuchos::rcp;
    using Teuchos::dyn_cast;

    // Get a pointer to the optimization vector
    RCP<const vector> xp = dyn_cast<const SV>(x).getVector();

    // Get a pointer to the gradient vector
    RCP<vector> gp = dyn_cast<SV>(g).getVector();

    int n = xp->size();

    // Create a vector of independent variables
    RCP<Fadvector> x_fad_rcp = rcp( new Fadvector );
    x_fad_rcp->reserve(n);

    // Initialize constructor for each element
    for(int i=0; i<n; ++i) {
        x_fad_rcp->push_back(FadType(n,i,(*xp)[i]));
    }

    StdVector<FadType> x_fad(x_fad_rcp);

    // AD access to objective function
    FadType J_fad = obj_.value(x_fad,tol);

    // Evaluate gradient
    for(int i=0; i<n; ++i) {
        (*gp)[i] = J_fad.dx(i);
    }
}



template <class Real, template<class> class Obj>
template <class ScalarT>
void Sacado_Objective<Real,Obj>::hessVecAD( Vector<ScalarT> &hv, const Vector<ScalarT> &v,
                                            const Vector<ScalarT> &x, Real &tol ) {

    // Data type which supports automatic differentiation
    typedef Sacado::Fad::SFad<ScalarT,1> FadType;
    typedef std::vector<FadType>         Fadvector;
    typedef std::vector<ScalarT>         vector;
    typedef StdVector<ScalarT>           SV;

    using Teuchos::RCP;       using Teuchos::rcp;
    using Teuchos::dyn_cast;

    // Get a pointer to the optimization vector
    RCP<const vector> xp = dyn_cast<const SV>(x).getVector();

    // Get a pointer to the direction vector
    RCP<const vector> vp = dyn_cast<const SV>(v).getVector();

    RCP<vector> hvp = dyn_cast<SV>(hv).getVector();


    int n = xp->size();

    // Create a vector of independent variables
    RCP<Fadvector> x_fad_rcp = rcp( new Fadvector );
    x_fad_rcp->reserve(n);

    // Allocate for gradient
    RCP<Fadvector> g_fad_rcp = rcp( new Fadvector );
    g_fad_rcp->resize(n);

    for(int i=0; i<n; ++i) {
        x_fad_rcp->push_back(FadType(1,(*xp)[i]));
    }

    // Set directional derivative
    for(int i=0; i<n; ++i) {
        (*x_fad_rcp)[i].fastAccessDx(0) = (*vp)[i];
    }

    StdVector<FadType> x_fad(x_fad_rcp);
    StdVector<FadType> g_fad(g_fad_rcp);

    this->gradientAD(g_fad,x_fad,tol);

    for(int i=0; i<n; ++i) {
        (*hvp)[i] = (*g_fad_rcp)[i].dx(0);
    }
}


#endif
