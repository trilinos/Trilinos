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


#ifndef ROL_SACADO_EQUALITYCONSTRAINT
#define ROL_SACADO_EQUALITYCONSTRAINT

#include "Sacado.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Constraint.hpp"

namespace ROL {

//! \brief ROL interface wrapper for Sacado Constraint
template<class Real, template<class> class Constr>
class Sacado_Constraint : public Constraint<Real> {

    protected:

        Constr<Real> constr_;

        int dim_;

        template<class ScalarT>
        void applyJacobianAD( Vector<ScalarT> &jv, const Vector<ScalarT> &v,
                              const Vector<ScalarT> &x, Real &tol );   
 
        template<class ScalarT>
        void applyAdjointJacobianAD( Vector<ScalarT> &aju, const Vector<ScalarT> &u,
                                     const Vector<ScalarT> &x, Real &tol);

        template<class ScalarT>
        void applyAdjointHessianAD( Vector<ScalarT> &ahuv, const Vector<ScalarT> &u,
                                    const Vector<ScalarT> &v, const Vector<ScalarT> &x, Real &tol);


    public: 
 
    Sacado_Constraint(int dim) : dim_(dim) {}
    Sacado_Constraint(const Constr<Real> & constr) : constr_(constr) {}

    virtual void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
        constr_.value(c,x,tol);
    }
    
    virtual void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, 
                       const Vector<Real> &x, Real &tol) {
        this->applyJacobianAD(jv,v,x,tol);
    }

    virtual void applyAdjointJacobian(Vector<Real> &aju, const Vector<Real> &u,
                              const Vector<Real> &x, Real &tol) {
        this->applyAdjointJacobianAD(aju,u,x,tol);
    } 

    void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u,
                             const Vector<Real> &v, const Vector<Real> &x, Real &tol){
        this->applyAdjointHessianAD(ahuv,u,v,x,tol);
    }

};


template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_Constraint<Real,Constr>::applyJacobianAD(Vector<ScalarT> &jv, const Vector<ScalarT> &v, 
                                                             const Vector<ScalarT> &x, Real &tol) {

    // Data type which supports automatic differentiation 
    typedef Sacado::Fad::DFad<ScalarT> FadType;
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;
  
    using Teuchos::RCP;       using Teuchos::rcp;
    using Teuchos::dyn_cast; 

    RCP<const vector> xp = dyn_cast<const SV>(x).getVector();

    int n = xp->size();

    // Get a pointer to the direction vector
    RCP<const vector> vp = dyn_cast<const SV>(v).getVector();

    RCP<vector> jvp = dyn_cast<SV>(jv).getVector();

    // Create a vector of independent variables
    RCP<Fadvector> x_fad_rcp = rcp( new Fadvector );
    x_fad_rcp->reserve(n);

    // Initialize constructor for each element
    for(int i=0; i<n; ++i) {
        x_fad_rcp->push_back(FadType(n,i,(*xp)[i])); 
    }

    // Create a vector of independent variables
    RCP<Fadvector> c_fad_rcp = rcp( new Fadvector );
    c_fad_rcp->reserve(dim_);

    for(int j=0; j<dim_; ++j) {
        c_fad_rcp->push_back(0);  
    }

    StdVector<FadType> x_fad(x_fad_rcp);
    StdVector<FadType> c_fad(c_fad_rcp);

    // Evaluate constraint     
    constr_.value(c_fad,x_fad,tol);

    for(int i=0; i<dim_; ++i) {
        (*jvp)[i] = 0;
        for(int j=0; j<n; ++j) {
            (*jvp)[i] += (*vp)[j]*(*c_fad_rcp)[i].dx(j); 
        }   
    }       
}



template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_Constraint<Real,Constr>::applyAdjointJacobianAD(Vector<ScalarT> &aju, const Vector<ScalarT> &u, 
                                                                    const Vector<ScalarT> &x, Real &tol) {

   // Data type which supports automatic differentiation 
    typedef Sacado::Fad::DFad<ScalarT> FadType;
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;
  
    using Teuchos::RCP;       using Teuchos::rcp;
    using Teuchos::dyn_cast; 

    // Get a pointer to the optimization vector
    RCP<const vector> xp = dyn_cast<const SV>(x).getVector();

    // Get a pointer to the direction vector
    RCP<const vector> up = dyn_cast<const SV>(u).getVector();

    RCP<vector> ajup = dyn_cast<SV>(aju).getVector();

    int n = xp->size();

    // Create a vector of independent variables
    RCP<Fadvector> x_fad_rcp = rcp( new Fadvector );
    x_fad_rcp->reserve(n);

    // Initialize constructor for each element
    for(int i=0; i<n; ++i) {
        x_fad_rcp->push_back(FadType(n,i,(*xp)[i])); 
    }

    RCP<Fadvector> c_fad_rcp = rcp( new Fadvector );
    c_fad_rcp->reserve(dim_);
    for(int j=0; j<dim_; ++j) {
        c_fad_rcp->push_back(0);  
    }

    StdVector<FadType> x_fad(x_fad_rcp);
    StdVector<FadType> c_fad(c_fad_rcp);
 
    // Evaluate constraint
    constr_.value(c_fad,x_fad,tol);
    
    FadType udotc = 0;
    
    for(int j=0;j<dim_;++j){ 
        udotc += (*c_fad_rcp)[j]*(*up)[j];
    } 

    for(int i=0;i<n;++i){
        (*ajup)[i] = udotc.dx(i);
    } 
}


template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_Constraint<Real,Constr>::applyAdjointHessianAD(Vector<ScalarT> &ahuv, const Vector<ScalarT> &u,
                                                                   const Vector<ScalarT> &v, const Vector<ScalarT> &x, 
                                                                   Real &tol){

    // Data type which supports automatic differentiation 
    typedef Sacado::Fad::SFad<ScalarT,1> FadType;
    typedef std::vector<FadType>         Fadvector;
    typedef std::vector<ScalarT>         vector;
    typedef StdVector<ScalarT>           SV;
  
    using Teuchos::RCP;       using Teuchos::rcp;
    using Teuchos::dyn_cast;

    // Get a pointer to the optimization vector
    RCP<const vector> xp = dyn_cast<const SV>(x).getVector();

    // Get a pointer to the dual constraint vector
    RCP<const vector> up = dyn_cast<const SV>(u).getVector();

    // Get a pointer to the direction vector
    RCP<const vector> vp = dyn_cast<const SV>(v).getVector();

    // Get a pointer to the directional adjoint Hessian 
    RCP<vector> ahuvp = dyn_cast<SV>(ahuv).getVector();

    // Number of optimization variables
    int n = xp->size();

    // Create a vector of independent variables
    RCP<Fadvector> x_fad_rcp = rcp( new Fadvector );
    x_fad_rcp->reserve(n);

    // Allocate for directional adjoint Jacobian
    RCP<Fadvector> aju_fad_rcp = rcp( new Fadvector );
    aju_fad_rcp->reserve(n);

    for(int i=0; i<n; ++i) {
        x_fad_rcp->push_back(FadType(1,(*xp)[i]));

        // Set derivative direction
        (*x_fad_rcp)[i].fastAccessDx(0) = (*vp)[i];     

        aju_fad_rcp->push_back(0);
    }

    // Allocate for constraint vector
    RCP<Fadvector> c_fad_rcp = rcp( new Fadvector );
    c_fad_rcp->reserve(dim_);

    // Allocate for dual constraint vector
    RCP<Fadvector> u_fad_rcp = rcp( new Fadvector );
    u_fad_rcp->reserve(dim_);

     for(int j=0; j<dim_; ++j) {
        u_fad_rcp->push_back((*up)[j]);
    }

    StdVector<FadType> x_fad(x_fad_rcp);
    StdVector<FadType> u_fad(u_fad_rcp);
    StdVector<FadType> aju_fad(aju_fad_rcp);

    // Evaluate constraint adjoint Jacobian direction
    this->applyAdjointJacobianAD( aju_fad, u_fad, x_fad, tol);

    for(int i=0; i<n; ++i) {
        (*ahuvp)[i] = (*aju_fad_rcp)[i].dx(0);             
    }
}


}
#endif
