// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SACADO_OBJECTIVE
#define ROL_SACADO_OBJECTIVE

#include "Sacado.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Objective_SimOpt.hpp"

using namespace ROL;

template <class Real, template<class> class Obj>
class Sacado_Objective_SimOpt : public Objective_SimOpt<Real> {

    protected:
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
     
    
    Sacado_Objective_SimOpt() : obj_(Obj<Real>()) {}
    Sacado_Objective_SimOpt(const Obj<Real> &obj) : obj_(obj) {}    

    
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
                                                    const Vector<ScalarT> &z, Real &tol) {
    typedef Sacado::Fad::DFad<ScalarT> FadType;
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

           
    

    ROL::Ptr<const vector> up = dynamic_cast<const SV&>(u).getVector();
    ROL::Ptr<const vector> zp = dynamic_cast<const SV&>(z).getVector();

    // Get a pointer to the gradient vector
    ROL::Ptr<vector> gp = dynamic_cast<SV&>(g).getVector();

    int m = zp->size();
    int n = up->size();

    ROL::Ptr<Fadvector> u_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> z_fad_ptr = ROL::makePtr<Fadvector>();
    
    u_fad_ptr->reserve(n);
    z_fad_ptr->reserve(m);

    // Initialize constructor for each element of u
    for(int i=0; i<n; ++i) {
        u_fad_ptr->push_back(FadType(n,i,(*up)[i]));
    }

    // Initialize constructor for each element of z
    for(int j=0; j<m; ++j) {
        z_fad_ptr->push_back((*zp)[j]);
    }

    StdVector<FadType> u_fad(u_fad_ptr);
    StdVector<FadType> z_fad(z_fad_ptr);
  
    FadType J_fad = obj_.value(u_fad,z_fad,tol);

    // Evaluate gradient
    for(int i=0; i<n; ++i) {
        (*gp)[i] = J_fad.dx(i);
    }
}


template<class Real, template<class> class Obj>
template<class ScalarT>
void Sacado_Objective_SimOpt<Real,Obj>::gradient_2AD(Vector<ScalarT> &g, const Vector<ScalarT> &u, 
                                                    const Vector<ScalarT> &z, Real &tol) {
    typedef Sacado::Fad::DFad<ScalarT> FadType;
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

           
     

    ROL::Ptr<const vector> up = dynamic_cast<const SV&>(u).getVector();
    ROL::Ptr<const vector> zp = dynamic_cast<const SV&>(z).getVector();

    // Get a pointer to the gradient vector
    ROL::Ptr<vector> gp = dynamic_cast<SV&>(g).getVector();

    int m = zp->size();
    int n = up->size();

    ROL::Ptr<Fadvector> u_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> z_fad_ptr = ROL::makePtr<Fadvector>();
    
    u_fad_ptr->reserve(n);
    z_fad_ptr->reserve(m);

    // Initialize constructor for each element of u
    for(int i=0; i<n; ++i) {
        u_fad_ptr->push_back((*up)[i]);
    }

    // Initialize constructor for each element of z
    for(int j=0; j<m; ++j) {
        z_fad_ptr->push_back(FadType(m,j,(*zp)[j]));
    }

    StdVector<FadType> u_fad(u_fad_ptr);
    StdVector<FadType> z_fad(z_fad_ptr);
  
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
    typedef std::vector<FadType>         Fadvector;
    typedef std::vector<ScalarT>         vector;
    typedef StdVector<ScalarT>           SV;

           
     

    ROL::Ptr<const vector> up = dynamic_cast<const SV&>(u).getVector();
    ROL::Ptr<const vector> zp = dynamic_cast<const SV&>(z).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    ROL::Ptr<vector> hvp = dynamic_cast<SV&>(hv).getVector();

    int n = up->size(); // vp and hvp have this size also
    int m = zp->size();

    ROL::Ptr<Fadvector> u_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> z_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> g_fad_ptr = ROL::makePtr<Fadvector>();
    
    u_fad_ptr->reserve(n);
    z_fad_ptr->reserve(m);
    g_fad_ptr->reserve(n);

    // Initialize constructor for each element of u
    for(int i=0; i<n; ++i) {
        u_fad_ptr->push_back(FadType(1,(*up)[i]));
        (*u_fad_ptr)[i].fastAccessDx(0) = (*vp)[i];
        g_fad_ptr->push_back(0);
    }

    // Initialize constructor for each element of z
    for(int j=0; j<m; ++j) {
        z_fad_ptr->push_back((*zp)[j]);
    }
    
    StdVector<FadType> u_fad(u_fad_ptr);
    StdVector<FadType> z_fad(z_fad_ptr);
    StdVector<FadType> g_fad(g_fad_ptr);
 
    this->gradient_1AD(g_fad,u_fad,z_fad,tol);

    for(int i=0; i<n; ++i) {
        (*hvp)[i] = (*g_fad_ptr)[i].dx(0); 
    }

}


template <class Real, template<class> class Obj>
template <class ScalarT>
void Sacado_Objective_SimOpt<Real,Obj>::hessVec_12AD(Vector<ScalarT> &hv, const Vector<ScalarT> &v, 
                                                     const Vector<ScalarT> &u, const Vector<ScalarT> &z, Real &tol) {
    typedef Sacado::Fad::DFad<ScalarT> FadType;
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

           
      

    ROL::Ptr<const vector> up = dynamic_cast<const SV&>(u).getVector();
    ROL::Ptr<const vector> zp = dynamic_cast<const SV&>(z).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    ROL::Ptr<vector> hvp = dynamic_cast<SV&>(hv).getVector();

    int n = up->size(); 
    int m = zp->size();

    ROL::Ptr<Fadvector> u_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> z_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> g_fad_ptr = ROL::makePtr<Fadvector>();
    
    u_fad_ptr->reserve(n);
    z_fad_ptr->reserve(m);
    g_fad_ptr->reserve(n);

     // Initialize constructor for each element of u
    for(int i=0; i<n; ++i) {
        u_fad_ptr->push_back((*up)[i]);
        g_fad_ptr->push_back(0);
    }

    // Initialize constructor for each element of z
    for(int j=0; j<m; ++j) {
        z_fad_ptr->push_back(FadType(m,j,(*zp)[j]));
    }
    
    StdVector<FadType> u_fad(u_fad_ptr);
    StdVector<FadType> z_fad(z_fad_ptr);
    StdVector<FadType> g_fad(g_fad_ptr);
 
    this->gradient_1AD(g_fad,u_fad,z_fad,tol);
   
    FadType vdotg = 0;

    for(int i=0; i<n; ++i) {
        vdotg += (*vp)[i]*(*g_fad_ptr)[i];
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
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

           
     

    ROL::Ptr<const vector> up = dynamic_cast<const SV&>(u).getVector();
    ROL::Ptr<const vector> zp = dynamic_cast<const SV&>(z).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    ROL::Ptr<vector> hvp = dynamic_cast<SV&>(hv).getVector();

    int n = up->size(); 
    int m = zp->size();

    ROL::Ptr<Fadvector> u_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> z_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> g_fad_ptr = ROL::makePtr<Fadvector>();
    
    u_fad_ptr->reserve(n);
    z_fad_ptr->reserve(m);
    g_fad_ptr->reserve(m);

     // Initialize constructor for each element of u
    for(int i=0; i<n; ++i) {
        u_fad_ptr->push_back(FadType(n,i,(*up)[i]));
    }

    // Initialize constructor for each element of z
    for(int j=0; j<m; ++j) {
        z_fad_ptr->push_back((*zp)[j]);
        g_fad_ptr->push_back(0);
    }
    
    StdVector<FadType> u_fad(u_fad_ptr);
    StdVector<FadType> z_fad(z_fad_ptr);
    StdVector<FadType> g_fad(g_fad_ptr);
 
    this->gradient_2AD(g_fad,u_fad,z_fad,tol);
   
    FadType vdotg = 0;

    for(int j=0; j<m; ++j) {
        vdotg += (*vp)[j]*(*g_fad_ptr)[j];
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
    typedef std::vector<FadType>         Fadvector;
    typedef std::vector<ScalarT>         vector;
    typedef StdVector<ScalarT>           SV;

           
    

    ROL::Ptr<const vector> up = dynamic_cast<const SV&>(u).getVector();
    ROL::Ptr<const vector> zp = dynamic_cast<const SV&>(z).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    ROL::Ptr<vector> hvp = dynamic_cast<SV&>(hv).getVector();
  
    int n = up->size();
    int m = zp->size(); // vp and hvp have this size also

    
    ROL::Ptr<Fadvector> u_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> z_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> g_fad_ptr = ROL::makePtr<Fadvector>();
    
    u_fad_ptr->reserve(n);
    z_fad_ptr->reserve(m);
    g_fad_ptr->reserve(m);

    // Initialize constructor for each element of u
    for(int i=0; i<n; ++i) {
        u_fad_ptr->push_back((*up)[i]);
    }

    // Initialize constructor for each element of z
    for(int j=0; j<m; ++j) {
        z_fad_ptr->push_back(FadType(1,(*zp)[j]));
        (*z_fad_ptr)[j].fastAccessDx(0) = (*vp)[j];
        g_fad_ptr->push_back(0);
    }
    
    StdVector<FadType> u_fad(u_fad_ptr);
    StdVector<FadType> z_fad(z_fad_ptr);
    StdVector<FadType> g_fad(g_fad_ptr);
 
    this->gradient_2AD(g_fad,u_fad,z_fad,tol);

    for(int j=0; j<m; ++j) {
        (*hvp)[j] = (*g_fad_ptr)[j].dx(0); 
    }

}



#endif
