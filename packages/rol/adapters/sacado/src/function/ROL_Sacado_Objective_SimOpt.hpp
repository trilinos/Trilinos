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
#include "ROL_Objective_SimOpt.hpp"

using namespace ROL;

template <class Real, template<class> class Obj>
class Sacado_Objective_SimOpt : public Objective_SimOpt<Real> {
    private:
        Obj<Real> obj_;

    template<class ScalarT>
    void gradient_1AD(Vector<ScalarT> &g, const Vector<ScalarT> &u, const Vector<ScalarT> &z, Real &tol);

    template<class ScalarT>
    void gradient_2AD(Vector<ScalarT> &g, const Vector<ScalarT> &u, const Vector<ScalarT> &z, Real &tol);

    template<class ScalarT>
    void hessVec_11AD(Vector<ScalarT> &hv, const Vector<ScalarT> &v, const Vector<ScalarT> &u, const Vector<ScalarT> &z, Real &tol);

    template<class ScalarT>
    void hessVec_12AD(Vector<ScalarT> &hv, const Vector<ScalarT> &v, const Vector<ScalarT> &u, const Vector<ScalarT> &z, Real &tol);

    template<class ScalarT>
    void hessVec_21AD(Vector<ScalarT> &hv, const Vector<ScalarT> &v, const Vector<ScalarT> &u, const Vector<ScalarT> &z, Real &tol);

    template<class ScalarT>
    void hessVec_22AD(Vector<ScalarT> &hv, const Vector<ScalarT> &v, const Vector<ScalarT> &u, const Vector<ScalarT> &z, Real &tol);




    public:
    
    Real value(const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
        return obj_.value(u,z,tol);
    }
    
    void gradient_1(Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
        this->gradient_1AD(g,u,z,tol);
    }

    void gradient_2(Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
        this->gradient_2AD(g,u,z,tol);
    }

    void hessVec_11(Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
        this->hessVec_11AD(hv,v,u,z,tol);  
    }

    void hessVec_12(Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
        this->hessVec_12AD(hv,v,u,z,tol);  
    }

    void hessVec_21(Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
        this->hessVec_21AD(hv,v,u,z,tol);  
    }

    void hessVec_22(Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
        this->hessVec_22AD(hv,v,u,z,tol);  
    }

};


template<class Real, template<class> class Obj>
template<class ScalarT>
void Sacado_Objective_SimOpt<Real,Obj>::gradient_1AD(Vector<ScalarT> &g, const Vector<ScalarT> &u, 
                                                    const Vector<ScalarT> &z, Real &tol {
    typedef Sacado::Fad::DFad<ScalarT> FadType;

    Teuchos::RCP<const std::vector<ScalarT> > up = 
        (Teuchos::dyn_cas<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(u))).getVector();

    Teuchos::RCP<const std::vector<ScalarT> > zp = 
        (Teuchos::dyn_cas<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(z))).getVector();

    // Get a pointer to the gradient vector
    Teuchos::RCP<std::vector<ScalarT> > gp =
        Teuchos::rcp_const_cast<std::vector<ScalarT> > ((Teuchos::dyn_cast<StdVector<ScalarT> > (g)).getVector());

    int m = zp->size();
    int n = up->size();

    Teuchos::RCP<std::vector<FadType> > u_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    Teuchos::RCP<std::vector<FadType> > z_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    
    u_fad_rcp->reserve(n);
    z_fad_rcp->reserve(m);

    // Initialize constructor for each element of u
    for(int i=0; i<n; ++i) {
        u_fad_rcp->push_back(FadType(n,i,(*up)[i]));
    }

    // Initialize constructor for each element of z
    for(int j=0; j<m; ++j) {
        z_fad_rcp->push_back((*zp)[j]);
    }

    StdVector<FadType> u_fad(u_fad_rcp);
    StdVector<FadType> z_fad(z_fad_rcp);
  
    FadType J_fad = obj_.value(u_fad,z_fad,tol);

    // Evaluate gradient
    for(int i=0; i<n; ++i) {
        (*gp)[i] = J_fad.dx(i);
    }
}


template<class Real, template<class> class Obj>
template<class ScalarT>
void Sacado_Objective_SimOpt<Real,Obj>::gradient_2AD(Vector<ScalarT> &g, const Vector<ScalarT> &u, 
                                                    const Vector<ScalarT> &z, Real &tol {
    typedef Sacado::Fad::DFad<ScalarT> FadType;

    Teuchos::RCP<const std::vector<ScalarT> > up = 
        (Teuchos::dyn_cas<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(u))).getVector();

    Teuchos::RCP<const std::vector<ScalarT> > zp = 
        (Teuchos::dyn_cas<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(z))).getVector();

    // Get a pointer to the gradient vector
    Teuchos::RCP<std::vector<ScalarT> > gp =
        Teuchos::rcp_const_cast<std::vector<ScalarT> > ((Teuchos::dyn_cast<StdVector<ScalarT> > (g)).getVector());

    int m = zp->size();
    int n = up->size();

    Teuchos::RCP<std::vector<FadType> > u_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    Teuchos::RCP<std::vector<FadType> > z_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    
    u_fad_rcp->reserve(n);
    z_fad_rcp->reserve(m);

    // Initialize constructor for each element of u
    for(int i=0; i<n; ++i) {
        u_fad_rcp->push_back((*up)[i]);
    }

    // Initialize constructor for each element of z
    for(int j=0; j<m; ++j) {
        z_fad_rcp->push_back(FadType(m,j,(*zp)[j]));
    }

    StdVector<FadType> u_fad(u_fad_rcp);
    StdVector<FadType> z_fad(z_fad_rcp);
  
    FadType J_fad = obj_.value(u_fad,z_fad,tol);

    // Evaluate gradient
    for(int j=0; j<m; ++j) {
        (*gp)[j] = J_fad.dx(j);
    }
}



template <class Real, template<class> class Obj>
template <class ScalarT>
void Sacado_Objective_SimOpt<Real,Obj>::hessVec_11AD(Vector<ScalarT> &hv, const Vector<ScalarT> &v, 
                                                     const Vector<ScalarT> &u, const Vector<ScalarT> &z, Real &tol) {

    typedef Sacado::Fad::SFad<ScalarT,1> FadType;

    Teuchos::RCP<const std::vector<ScalarT> > up = 
        (Teuchos::dyn_cas<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(u))).getVector();

    Teuchos::RCP<const std::vector<ScalarT> > zp = 
        (Teuchos::dyn_cas<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(z))).getVector();

    Teuchos::RCP<const std::vector<ScalarT> > vp = 
        (Teuchos::dyn_cas<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(v))).getVector();
 
    Teuchos::RCP<std::vector<ScalarT> > hvp =
        Teuchos::rcp_const_cast<std::vector<ScalarT> > ((Teuchos::dyn_cast<StdVector<ScalarT> > (hv)).getVector());

    int n = up->size(); // vp and hvp have this size also
    int m = zp->size();

    Teuchos::RCP<std::vector<FadType> > u_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    Teuchos::RCP<std::vector<FadType> > z_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    Teuchos::RCP<std::vector<FadType> > g_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    
    u_fad_rcp->reserve(n);
    z_fad_rcp->reserve(m);
    g_fad_rcp->reserve(n);

    // Initialize constructor for each element of u
    for(int i=0; i<n; ++i) {
        u_fad_rcp->push_back(FadType(n,i,(*up)[i]));
        g_fad_rcp->push_back(0);
    }

    // Initialize constructor for each element of z
    for(int j=0; j<m; ++j) {
        z_fad_rcp->push_back((*zp)[j]);
    }
    
    StdVector<FadType> u_fad(u_fad_rcp);
    StdVector<FadType> z_fad(z_fad_rcp);
    StdVector<FadType> g_fad(g_fad_rcp);
 
    this->gradient_1AD(g_fad,u_fad,z_fad,tol);

    for(int i=0; i<n; ++i) {
        (*hvp)[i] = (*g_fad_rcp)[i].dx(0); 
    }

}


template <class Real, template<class> class Obj>
template <class ScalarT>
void Sacado_Objective_SimOpt<Real,Obj>::hessVec_12AD(Vector<ScalarT> &hv, const Vector<ScalarT> &v, 
                                                     const Vector<ScalarT> &u, const Vector<ScalarT> &z, Real &tol) {
    typedef Sacado::Fad::DFad<ScalarT> FadType;

    Teuchos::RCP<const std::vector<ScalarT> > up = 
        (Teuchos::dyn_cas<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(u))).getVector();

    Teuchos::RCP<const std::vector<ScalarT> > zp = 
        (Teuchos::dyn_cas<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(z))).getVector();

    Teuchos::RCP<const std::vector<ScalarT> > vp = 
        (Teuchos::dyn_cas<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(v))).getVector();
 
    Teuchos::RCP<std::vector<ScalarT> > hvp =
        Teuchos::rcp_const_cast<std::vector<ScalarT> > ((Teuchos::dyn_cast<StdVector<ScalarT> > (hv)).getVector());

    int n = up->size(); 
    int m = zp->size();

    Teuchos::RCP<std::vector<FadType> > u_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    Teuchos::RCP<std::vector<FadType> > z_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    Teuchos::RCP<std::vector<FadType> > g_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    
    u_fad_rcp->reserve(n);
    z_fad_rcp->reserve(m);
    g_fad_rcp->reserve(n);

     // Initialize constructor for each element of u
    for(int i=0; i<n; ++i) {
        u_fad_rcp->push_back((*up)[i]);
        g_fad_rcp->push_back(0);
    }

    // Initialize constructor for each element of z
    for(int j=0; j<m; ++j) {
        z_fad_rcp->push_back(FadType(m,j,(*zp)[j]));
    }
    
    StdVector<FadType> u_fad(u_fad_rcp);
    StdVector<FadType> z_fad(z_fad_rcp);
    StdVector<FadType> g_fad(g_fad_rcp);
 
    this->gradient_1AD(g_fad,u_fad,z_fad,tol);
   
    FadType vdotg = 0;

    for(int i=0; i<n; ++i) {
        vdotg += (*vp)[i]*(*g_fad_rcp)[i];
    }

    for(int j=0; j<m; ++j) {
        (*hvp)[j] = vdotg.dx(j); 
    }
}



template <class Real, template<class> class Obj>
template <class ScalarT>
void Sacado_Objective_SimOpt<Real,Obj>::hessVec_21AD(Vector<ScalarT> &hv, const Vector<ScalarT> &v, 
                                                     const Vector<ScalarT> &u, const Vector<ScalarT> &z, Real &tol) {
    typedef Sacado::Fad::DFad<ScalarT> FadType;

    Teuchos::RCP<const std::vector<ScalarT> > up = 
        (Teuchos::dyn_cas<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(u))).getVector();

    Teuchos::RCP<const std::vector<ScalarT> > zp = 
        (Teuchos::dyn_cas<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(z))).getVector();

    Teuchos::RCP<const std::vector<ScalarT> > vp = 
        (Teuchos::dyn_cas<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(v))).getVector();
 
    Teuchos::RCP<std::vector<ScalarT> > hvp =
        Teuchos::rcp_const_cast<std::vector<ScalarT> > ((Teuchos::dyn_cast<StdVector<ScalarT> > (hv)).getVector());

    int n = up->size(); 
    int m = zp->size();

    Teuchos::RCP<std::vector<FadType> > u_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    Teuchos::RCP<std::vector<FadType> > z_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    Teuchos::RCP<std::vector<FadType> > g_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    
    u_fad_rcp->reserve(n);
    z_fad_rcp->reserve(m);
    g_fad_rcp->reserve(m);

     // Initialize constructor for each element of u
    for(int i=0; i<n; ++i) {
        u_fad_rcp->push_back(FadType(n,i,(*up)[i]));
    }

    // Initialize constructor for each element of z
    for(int j=0; j<m; ++j) {
        z_fad_rcp->push_back((*zp)[j]);
        g_fad_rcp->push_back(0);
    }
    
    StdVector<FadType> u_fad(u_fad_rcp);
    StdVector<FadType> z_fad(z_fad_rcp);
    StdVector<FadType> g_fad(g_fad_rcp);
 
    this->gradient_2AD(g_fad,u_fad,z_fad,tol);
   
    FadType vdotg = 0;

    for(int j=0; i<m; ++j) {
        vdotg += (*vp)[j]*(*g_fad_rcp)[j];
    }

    for(int i=0; i<n; ++i) {
        (*hvp)[i] = vdotg.dx(i); 
    }
}





template <class Real, template<class> class Obj>
template <class ScalarT>
void Sacado_Objective_SimOpt<Real,Obj>::hessVec_22AD(Vector<ScalarT> &hv, const Vector<ScalarT> &v, 
                                                     const Vector<ScalarT> &u, const Vector<ScalarT> &z, Real &tol) {

    typedef Sacado::Fad::SFad<ScalarT,1> FadType;

    Teuchos::RCP<const std::vector<ScalarT> > up = 
        (Teuchos::dyn_cas<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(u))).getVector();

    Teuchos::RCP<const std::vector<ScalarT> > zp = 
        (Teuchos::dyn_cas<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(z))).getVector();

    Teuchos::RCP<const std::vector<ScalarT> > vp = 
        (Teuchos::dyn_cas<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(v))).getVector();
 
    Teuchos::RCP<std::vector<ScalarT> > hvp =
        Teuchos::rcp_const_cast<std::vector<ScalarT> > ((Teuchos::dyn_cast<StdVector<ScalarT> > (hv)).getVector());

    int n = up->size();
    int m = zp->size(); // vp and hvp have this size also

    Teuchos::RCP<std::vector<FadType> > u_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    Teuchos::RCP<std::vector<FadType> > z_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    Teuchos::RCP<std::vector<FadType> > g_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    
    u_fad_rcp->reserve(n);
    z_fad_rcp->reserve(m);
    g_fad_rcp->reserve(m);

    // Initialize constructor for each element of u
    for(int i=0; i<n; ++i) {
        u_fad_rcp->push_back((*up)[i]);
    }

    // Initialize constructor for each element of z
    for(int j=0; j<m; ++j) {
        z_fad_rcp->push_back(FadType(m,j,(*zp)[j]));
        g_fad_rcp->push_back(0);
    }
    
    StdVector<FadType> u_fad(u_fad_rcp);
    StdVector<FadType> z_fad(z_fad_rcp);
    StdVector<FadType> g_fad(g_fad_rcp);
 
    this->gradient_2AD(g_fad,u_fad,z_fad,tol);

    for(int j=0; j<m; ++j) {
        (*hvp)[j] = (*g_fad_rcp)[j].dx)(0); 
    }

}



#endif
