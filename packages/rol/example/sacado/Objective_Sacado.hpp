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

#include "Sacado.hpp"

/** \brief Generic objective wrapper class for class that uses Sacado */
template<class Real, class Obj>
class ROL_Objective_Sacado : public Objective<Real> {
    private:
        Obj obj_;
    public:
    /* Evaluate the objective function at x */
    Real value( const Vector<Real> &x, Real &tol ) {
      return obj_.value(x,tol);
    }

    /* Evaluate the gradient at x */
    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
        obj_.gradient(g,x,tol); 
    }

    /* Compute the action of the Hessian evaluated at x on a vector v */
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
        obj_.hessVec(hv,v,x,tol);
    }
};



template<class Real, class Obj>
class Objective_Sacado {
    private:
        Obj obj_;
    public:

    /* Evaluate the objective function at x */
    template<class ScalarT> 
    ScalarT value( const Vector<ScalarT> &x, Real &tol ) {
      return obj_.value(x,tol);
    }

    /* Evaluate the gradient at x */
    template<class ScalarT> 
    void gradient( Vector<ScalarT> &g, const Vector<ScalarT> &x, Real &tol );

    /* Compute the action of the Hessian evaluated at x on a vector v */
    template<class ScalarT> 
    void hessVec( Vector<ScalarT> &hv, const Vector<ScalarT> &v, const Vector<ScalarT> &x, Real &tol ); 
};



/** \brief A Sacado-accesible function which computes the gradient of the Zakharov function 
    using automatic differentiation.
    @param[in]  x is the optimization vector
    @param[out] J is the value of the objective function
    @param[out] g is the gradient of the objective function    
*/
template<class Real, class Obj>
template<class ScalarT>
void Objective_Sacado<Real,Obj>::gradient(Vector<ScalarT> &g, const Vector<ScalarT> &x, Real &tol) { 

    // Data type which supports automatic differentiation 
    typedef Sacado::Fad::DFad<ScalarT> FadType;

    // Get a pointer to the optimization vector
    Teuchos::RCP<const std::vector<ScalarT> > xp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(x))).getVector();

    // Get a pointer to the gradient vector
    Teuchos::RCP<std::vector<ScalarT> > gp =
        Teuchos::rcp_const_cast<std::vector<ScalarT> > ((Teuchos::dyn_cast<StdVector<ScalarT> > (g)).getVector());
    
    int n = xp->size();
 
    // Create a vector of independent variables
    Teuchos::RCP<std::vector<FadType> > x_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    x_fad_rcp->reserve(n);
   
    // Initialize constructor for each element
    for(int i=0; i<n; ++i) {
        x_fad_rcp->push_back(FadType(n,i,(*xp)[i])); 
    }

    StdVector<FadType> x_fad(x_fad_rcp);

    // AD access to objective function
    FadType J_fad = this->value(x_fad,tol);

    // Evaluate gradient
    for(int i=0; i<n; ++i) {
        (*gp)[i] = J_fad.dx(i);
    }
}



/** \brief A Sacado-accesible function which computes the action of the Hessian on a given direction vector
    @param[in]   x is the optimization vector
    @param[in]   v is a given direction vector 
    @param[out]  J is the value of the objective function
    @param[out]  g is the gradient of the objective function
    @param[out]  hv is the action of the Hessian the direction v 
*/
template <class Real, class Obj>
template <class ScalarT>
void Objective_Sacado<Real,Obj>::hessVec( Vector<ScalarT> &hv, const Vector<ScalarT> &v, const Vector<ScalarT> &x, Real &tol ) {
  
    // Get a pointer to the optimization vector 
    Teuchos::RCP<const std::vector<ScalarT> > xp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(x))).getVector();

    // Get a pointer to the direction vector
    Teuchos::RCP<const std::vector<ScalarT> > vp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(v))).getVector();

    Teuchos::RCP<std::vector<ScalarT> > hvp =
        Teuchos::rcp_const_cast<std::vector<ScalarT> >((Teuchos::dyn_cast<StdVector<ScalarT> >(hv)).getVector());

    typedef Sacado::Fad::SFad<ScalarT,1> FadType;

    int n = xp->size();
   
    // Create a vector of independent variables
    Teuchos::RCP<std::vector<FadType> > x_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    x_fad_rcp->reserve(n);
 
    // Allocate for gradient   
    Teuchos::RCP<std::vector<FadType> > g_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    g_fad_rcp->reserve(n);

    for(int i=0; i<n; ++i) {
        x_fad_rcp->push_back(FadType(1,(*xp)[i]));
    }

    // Set directional derivative    
    for(int i=0; i<n; ++i) {
        (*x_fad_rcp)[i].fastAccessDx(0) = (*vp)[i];
    }
    
    StdVector<FadType> x_fad(x_fad_rcp);
    StdVector<FadType> g_fad(g_fad_rcp);

    this->gradient(g_fad,x_fad,tol);

    for(int i=0; i<n; ++i) {
        (*hvp)[i] = (*g_fad_rcp)[i].dx(0);            
    }
}


