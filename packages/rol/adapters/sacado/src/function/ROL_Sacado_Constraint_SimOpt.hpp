// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SACADO_EQUALITYCONSTRAINT_SIMOPT
#define ROL_SACADO_EQUALITYCONSTRAINT_SIMOPT

#include "Sacado.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Constraint_SimOpt.hpp"

using namespace ROL;

//! \brief ROL interface wrapper for Sacado SimOpt Constraint
template<class Real, template<class> class Constr>
class Sacado_Constraint_SimOpt : public Constraint_SimOpt<Real> {
 

   
    protected:     
        Constr<Real> constr_;

        template<class ScalarT>
        void applyJacobian_1AD(Vector<ScalarT> &jv, const Vector<ScalarT> &v, const Vector<ScalarT> &u, 
                               const Vector<ScalarT> &z, Real &tol);

        template<class ScalarT>
        void applyJacobian_2AD(Vector<ScalarT> &jv, const Vector<ScalarT> &v, const Vector<ScalarT> &u, 
                               const Vector<ScalarT> &z, Real &tol);

        template<class ScalarT>
        void applyAdjointJacobian_1AD(Vector<ScalarT> &ajv, const Vector<ScalarT> &v, const Vector<ScalarT> &u, 
                                      const Vector<ScalarT> &z, Real &tol);

        template<class ScalarT>
        void applyAdjointJacobian_2AD(Vector<ScalarT> &ajv, const Vector<ScalarT> &v, const Vector<ScalarT> &u, 
                                      const Vector<ScalarT> &z, Real &tol);
         
        template<class ScalarT>
        void applyAdjointHessian_11AD(Vector<ScalarT> &ahwv, const Vector<ScalarT> &w, 
                                      const Vector<ScalarT> &v, const Vector<ScalarT> &u,
                                      const Vector<ScalarT> &z, Real &tol);         
        template<class ScalarT>
        void applyAdjointHessian_12AD(Vector<ScalarT> &ahwv, const Vector<ScalarT> &w, 
                                      const Vector<ScalarT> &v, const Vector<ScalarT> &u,
                                      const Vector<ScalarT> &z, Real &tol);         
        template<class ScalarT>
        void applyAdjointHessian_21AD(Vector<ScalarT> &ahwv, const Vector<ScalarT> &w, 
                                      const Vector<ScalarT> &v, const Vector<ScalarT> &u,
                                      const Vector<ScalarT> &z, Real &tol);         
        template<class ScalarT>
        void applyAdjointHessian_22AD(Vector<ScalarT> &ahwv, const Vector<ScalarT> &w, 
                                      const Vector<ScalarT> &v, const Vector<ScalarT> &u,
                                      const Vector<ScalarT> &z, Real &tol);         
 
    public:
        Sacado_Constraint_SimOpt() : constr_(Constr<Real>()) {}
        Sacado_Constraint_SimOpt(Constr<Real> constr) : constr_(constr) { }

        void value(Vector<Real> &c, const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
            constr_.value(c,u,z,tol);
        }

        void applyJacobian_1(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &u, 
                             const Vector<Real> &z, Real &tol) {
            this->applyJacobian_1AD(jv,v,u,z,tol);
        } 

        void applyJacobian_2(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &u, 
                             const Vector<Real> &z, Real &tol) {
            this->applyJacobian_2AD(jv,v,u,z,tol);
        } 

        void applyAdjointJacobian_1(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &u, 
                                    const Vector<Real> &z, Real &tol) {
            this->applyAdjointJacobian_1AD(ajv,v,u,z,tol);
        } 

        void applyAdjointJacobian_2(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &u, 
                                    const Vector<Real> &z, Real &tol) {
            this->applyAdjointJacobian_2AD(ajv,v,u,z,tol);
        } 

        void applyAdjointHessian_11(Vector<Real> &ahwv, const Vector<Real> &w, const Vector<Real> &v,
                                    const Vector<Real> &u, const Vector<Real> &z, Real &tol){
            this->applyAdjointHessian_11AD(ahwv,w,v,u,z,tol);
        }
 
        void applyAdjointHessian_12(Vector<Real> &ahwv, const Vector<Real> &w, const Vector<Real> &v,
                                    const Vector<Real> &u, const Vector<Real> &z, Real &tol){
            this->applyAdjointHessian_12AD(ahwv,w,v,u,z,tol);
        }
 
        void applyAdjointHessian_21(Vector<Real> &ahwv, const Vector<Real> &w, const Vector<Real> &v,
                                    const Vector<Real> &u, const Vector<Real> &z, Real &tol){
            this->applyAdjointHessian_21AD(ahwv,w,v,u,z,tol);
        }

        void applyAdjointHessian_22(Vector<Real> &ahwv, const Vector<Real> &w, const Vector<Real> &v,
                                    const Vector<Real> &u, const Vector<Real> &z, Real &tol){
            this->applyAdjointHessian_22AD(ahwv,w,v,u,z,tol);
        }
};



template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_Constraint_SimOpt<Real,Constr>::applyJacobian_1AD(Vector<ScalarT> &jv, const Vector<ScalarT> &v, const Vector<ScalarT> &u, 
                                                                      const Vector<ScalarT> &z, Real &tol) {
    //    v in U (n), jv in C (n)    

    typedef Sacado::Fad::DFad<ScalarT> FadType;
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

           
     

    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    ROL::Ptr<const vector> up = dynamic_cast<const SV&>(u).getVector();
    ROL::Ptr<const vector> zp = dynamic_cast<const SV&>(z).getVector();
    ROL::Ptr<vector> jvp = dynamic_cast<SV&>(jv).getVector(); 

    int n = up->size();
    int m = zp->size();
    
    ROL::Ptr<Fadvector> c_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> u_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> z_fad_ptr = ROL::makePtr<Fadvector>();

    c_fad_ptr->reserve(n);
    u_fad_ptr->reserve(n);
    z_fad_ptr->reserve(m);

    for(int i=0; i<n; ++i) {
        c_fad_ptr->push_back(0);
        u_fad_ptr->push_back(FadType(n,i,(*up)[i]));
    }

    for(int j=0; j<m; ++j) {
        z_fad_ptr->push_back((*zp)[j]); 
    }

    StdVector<FadType> c_fad(c_fad_ptr);
    StdVector<FadType> u_fad(u_fad_ptr);
    StdVector<FadType> z_fad(z_fad_ptr);
  
    // Evaluate constraint 
    constr_.value(c_fad,u_fad,z_fad,tol);

    for(int i=0; i<n; ++i) {
        (*jvp)[i] = 0; 
        for(int j=0; j<n; ++j) {
            (*jvp)[i] += (*vp)[j]*(*c_fad_ptr)[i].dx(j);
        }
    } 
}


template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_Constraint_SimOpt<Real,Constr>::applyJacobian_2AD(Vector<ScalarT> &jv, const Vector<ScalarT> &v, const Vector<ScalarT> &u, 
                                                                      const Vector<ScalarT> &z, Real &tol) {
    // v in Z (m), jv in C (n)

    typedef Sacado::Fad::DFad<ScalarT> FadType;
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

           
    

    ROL::Ptr<vector> jvp = dynamic_cast<SV&>(jv).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    ROL::Ptr<const vector> up = dynamic_cast<const SV&>(u).getVector();
    ROL::Ptr<const vector> zp = dynamic_cast<const SV&>(z).getVector();

    int n = up->size();
    int m = zp->size(); 
    
    ROL::Ptr<Fadvector> c_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> u_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> z_fad_ptr = ROL::makePtr<Fadvector>();

    c_fad_ptr->reserve(n);
    u_fad_ptr->reserve(n);
    z_fad_ptr->reserve(m);

    for(int i=0; i<n; ++i) {
        c_fad_ptr->push_back(0);
        u_fad_ptr->push_back((*up)[i]); 
    }

    for(int j=0; j<m; ++j) {
        z_fad_ptr->push_back(FadType(n,j,(*zp)[j]));
    }

    StdVector<FadType> c_fad(c_fad_ptr);
    StdVector<FadType> u_fad(u_fad_ptr);
    StdVector<FadType> z_fad(z_fad_ptr);
  
    // Evaluate constraint 
    constr_.value(c_fad,u_fad,z_fad,tol);

    for(int i=0; i<n; ++i) {
        (*jvp)[i] = 0; 
        for(int j=0; j<n; ++j) {
            (*jvp)[i] += (*vp)[j]*(*c_fad_ptr)[i].dx(j);
        }
    } 
 
}


template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_Constraint_SimOpt<Real,Constr>::applyAdjointJacobian_1AD(Vector<ScalarT> &ajv, const Vector<ScalarT> &v,
                                                                             const Vector<ScalarT> &u, const Vector<ScalarT> &z, Real &tol) {

    // v in C* (n), ajv in U* (n)

    typedef Sacado::Fad::DFad<ScalarT> FadType;
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

           
      

    ROL::Ptr<vector> ajvp = dynamic_cast<SV&>(ajv).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    ROL::Ptr<const vector> up = dynamic_cast<const SV&>(u).getVector();
    ROL::Ptr<const vector> zp = dynamic_cast<const SV&>(z).getVector();
    
    int n = up->size();
    int m = zp->size();
    
    ROL::Ptr<Fadvector> c_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> u_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> z_fad_ptr = ROL::makePtr<Fadvector>();

    c_fad_ptr->reserve(n);
    u_fad_ptr->reserve(n);
    z_fad_ptr->reserve(m);

    for(int i=0; i<n; ++i) {
        c_fad_ptr->push_back(0);
        u_fad_ptr->push_back(FadType(n,i,(*up)[i]));
    }

    for(int j=0; j<m; ++j) {
        z_fad_ptr->push_back((*zp)[j]); 
    }

    StdVector<FadType> c_fad(c_fad_ptr);
    StdVector<FadType> u_fad(u_fad_ptr);
    StdVector<FadType> z_fad(z_fad_ptr);
  
    // Evaluate constraint 
    constr_.value(c_fad,u_fad,z_fad,tol);

    FadType vdotc = 0;

    for(int i=0;i<n;++i) {
        vdotc += (*c_fad_ptr)[i]*(*vp)[i]; 
    } 

    for(int i=0;i<n;++i) {
        (*ajvp)[i] = vdotc.dx(i);
    }

}



template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_Constraint_SimOpt<Real,Constr>::applyAdjointJacobian_2AD(Vector<ScalarT> &ajv, const Vector<ScalarT> &v, 
                                                                             const Vector<ScalarT> &u, const Vector<ScalarT> &z, 
                                                                             Real &tol) {
    // v in C* (n), ajv in Z* (m)

    typedef Sacado::Fad::DFad<ScalarT> FadType;
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

           
      

    ROL::Ptr<vector> ajvp = dynamic_cast<SV&>(ajv).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();   
    ROL::Ptr<const vector> up = dynamic_cast<const SV&>(u).getVector();   
    ROL::Ptr<const vector> zp = dynamic_cast<const SV&>(z).getVector();   

    int n = up->size();
    int m = zp->size();
    
    ROL::Ptr<Fadvector> c_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> u_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> z_fad_ptr = ROL::makePtr<Fadvector>();

    c_fad_ptr->reserve(n);
    u_fad_ptr->reserve(n);
    z_fad_ptr->reserve(m);

    for(int i=0; i<n; ++i) {
        c_fad_ptr->push_back(0);
        u_fad_ptr->push_back((*up)[i]); 
    }

    for(int j=0; j<m; ++j) {
        z_fad_ptr->push_back(FadType(n,j,(*zp)[j]));
    }

    StdVector<FadType> c_fad(c_fad_ptr);
    StdVector<FadType> u_fad(u_fad_ptr);
    StdVector<FadType> z_fad(z_fad_ptr);
  
    // Evaluate constraint 
    constr_.value(c_fad,u_fad,z_fad,tol);

    FadType vdotc = 0;

    for(int i=0;i<n;++i) {
        vdotc += (*c_fad_ptr)[i]*(*vp)[i]; 
    } 

    for(int j=0;j<m;++j) {
        (*ajvp)[j] = vdotc.dx(j);
    }
}



template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_Constraint_SimOpt<Real,Constr>::applyAdjointHessian_11AD(Vector<ScalarT> &ahwv, const Vector<ScalarT> &w, 
                                                                             const Vector<ScalarT> &v, const Vector<ScalarT> &u,
                                                                             const Vector<ScalarT> &z, Real &tol) {
    // w in C* (n), v in U* (n), ahwv in U* (n)

    typedef Sacado::Fad::DFad<ScalarT> FadType;    
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

           
     

    ROL::Ptr<vector> ahwvp = dynamic_cast<SV&>(ahwv).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    ROL::Ptr<const vector> wp = dynamic_cast<const SV&>(w).getVector();
    ROL::Ptr<const vector> up = dynamic_cast<const SV&>(u).getVector();
    ROL::Ptr<const vector> zp = dynamic_cast<const SV&>(z).getVector();

    int n = up->size();
    int m = zp->size();

    ROL::Ptr<Fadvector> v_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> u_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> z_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> jv_fad_ptr = ROL::makePtr<Fadvector>();

    v_fad_ptr->reserve(n);
    u_fad_ptr->reserve(n);
    z_fad_ptr->reserve(m);
    jv_fad_ptr->reserve(n);

    for(int i=0; i<n; ++i) {
        v_fad_ptr->push_back((*vp)[i]);
        u_fad_ptr->push_back(FadType(n,i,(*up)[i]));
        jv_fad_ptr->push_back(0);
    }

    for(int j=0; j<m; ++j) {
        z_fad_ptr->push_back((*zp)[j]);    
    }

    StdVector<FadType> v_fad(v_fad_ptr);     
    StdVector<FadType> u_fad(u_fad_ptr);     
    StdVector<FadType> z_fad(z_fad_ptr);     
    StdVector<FadType> jv_fad(jv_fad_ptr);     

    this->applyJacobian_1AD(jv_fad,v_fad,u_fad,z_fad,tol);
 
    FadType wjv_fad = 0;
   
    for(int i=0; i<n; ++i) {
        wjv_fad += (*wp)[i]*(*jv_fad_ptr)[i];
    }

    for(int i=0; i<n; ++i) {
        (*ahwvp)[i] = wjv_fad.dx(i);
    }
}




template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_Constraint_SimOpt<Real,Constr>::applyAdjointHessian_12AD(Vector<ScalarT> &ahwv, const Vector<ScalarT> &w, 
                                                                             const Vector<ScalarT> &v, const Vector<ScalarT> &u,
                                                                             const Vector<ScalarT> &z, Real &tol) {
    // w in C* (n), v in U* (n), ahwv in Z* (m)

    typedef Sacado::Fad::DFad<ScalarT> FadType;    
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

           
      

    ROL::Ptr<vector> ahwvp = dynamic_cast<SV&>(ahwv).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    ROL::Ptr<const vector> wp = dynamic_cast<const SV&>(w).getVector();
    ROL::Ptr<const vector> up = dynamic_cast<const SV&>(u).getVector();
    ROL::Ptr<const vector> zp = dynamic_cast<const SV&>(z).getVector();

    int n = up->size();
    int m = zp->size();

    ROL::Ptr<Fadvector> v_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> u_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> z_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> jv_fad_ptr = ROL::makePtr<Fadvector>();

    v_fad_ptr->reserve(n);
    u_fad_ptr->reserve(n);
    z_fad_ptr->reserve(m);
    jv_fad_ptr->reserve(n);

    for(int i=0; i<n; ++i) {
        v_fad_ptr->push_back((*vp)[i]);
        u_fad_ptr->push_back((*up)[i]);
        jv_fad_ptr->push_back(0);
    }

    for(int j=0; j<m; ++j) {
        z_fad_ptr->push_back(FadType(m,j,(*zp)[j]));    
    }

    StdVector<FadType> v_fad(v_fad_ptr);     
    StdVector<FadType> u_fad(u_fad_ptr);     
    StdVector<FadType> z_fad(z_fad_ptr);     
    StdVector<FadType> jv_fad(jv_fad_ptr);     

    this->applyJacobian_1AD(jv_fad,v_fad,u_fad,z_fad,tol);
    
    FadType wjv_fad = 0;

    for(int i=0; i<n; ++i) {
        wjv_fad += (*wp)[i]*(*jv_fad_ptr)[i]; 
    }
 
    for(int j=0; j<m; ++j) {
        (*ahwvp)[j] = wjv_fad.dx(j);
    }
}


template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_Constraint_SimOpt<Real,Constr>::applyAdjointHessian_21AD(Vector<ScalarT> &ahwv, const Vector<ScalarT> &w, 
                                                                             const Vector<ScalarT> &v, const Vector<ScalarT> &u,
                                                                             const Vector<ScalarT> &z, Real &tol) {
    // w in C* (n), v in Z* (m), ahwv in U* (n)

    typedef Sacado::Fad::DFad<ScalarT> FadType;    
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

           
     

    ROL::Ptr<vector> ahwvp = dynamic_cast<SV&>(ahwv).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    ROL::Ptr<const vector> wp = dynamic_cast<const SV&>(w).getVector();
    ROL::Ptr<const vector> up = dynamic_cast<const SV&>(u).getVector();
    ROL::Ptr<const vector> zp = dynamic_cast<const SV&>(z).getVector();

    int n = up->size();
    int m = zp->size();

    ROL::Ptr<Fadvector> v_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> u_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> z_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> jv_fad_ptr = ROL::makePtr<Fadvector>();

    v_fad_ptr->reserve(m);
    u_fad_ptr->reserve(n);
    z_fad_ptr->reserve(m);
    jv_fad_ptr->reserve(n);

    for(int i=0; i<n; ++i) {
        u_fad_ptr->push_back(FadType(1,(*up)[i]));
        jv_fad_ptr->push_back(0);
    }

    for(int j=0; j<m; ++j) {
        v_fad_ptr->push_back((*vp)[j]);
        z_fad_ptr->push_back((*zp)[j]);    
    }

    StdVector<FadType> v_fad(v_fad_ptr);     
    StdVector<FadType> u_fad(u_fad_ptr);     
    StdVector<FadType> z_fad(z_fad_ptr);     
    StdVector<FadType> jv_fad(jv_fad_ptr);     

    this->applyJacobian_2AD(jv_fad,v_fad,u_fad,z_fad,tol);

    FadType wjv_fad = 0;

    for(int i=0; i<n; ++i) {
        wjv_fad += (*wp)[i]*(*jv_fad_ptr)[i];
    }

    for(int i=0; i<n; ++i) {
        (*ahwvp)[i] = wjv_fad.dx(i);
    }
}


template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_Constraint_SimOpt<Real,Constr>::applyAdjointHessian_22AD(Vector<ScalarT> &ahwv, const Vector<ScalarT> &w, 
                                                                             const Vector<ScalarT> &v, const Vector<ScalarT> &u,
                                                                             const Vector<ScalarT> &z, Real &tol) {
    // w in C* (n), v in Z* (m), ahwv in Z* (m)

    typedef Sacado::Fad::DFad<ScalarT> FadType;    
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

           
     

    ROL::Ptr<vector> ahwvp = dynamic_cast<SV&>(ahwv).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    ROL::Ptr<const vector> wp = dynamic_cast<const SV&>(w).getVector();
    ROL::Ptr<const vector> up = dynamic_cast<const SV&>(u).getVector();
    ROL::Ptr<const vector> zp = dynamic_cast<const SV&>(z).getVector();
 
    int n = up->size();
    int m = zp->size();

    ROL::Ptr<Fadvector> v_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> u_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> z_fad_ptr = ROL::makePtr<Fadvector>();
    ROL::Ptr<Fadvector> jv_fad_ptr = ROL::makePtr<Fadvector>();

    v_fad_ptr->reserve(m);
    u_fad_ptr->reserve(n);
    z_fad_ptr->reserve(m);
    jv_fad_ptr->reserve(n);

    for(int i=0; i<n; ++i) {
        u_fad_ptr->push_back((*up)[i]);
        jv_fad_ptr->push_back(0);
    }

    for(int j=0; j<m; ++j) {
        v_fad_ptr->push_back((*vp)[j]);
        z_fad_ptr->push_back(FadType(m,j,(*zp)[j]));    
    }

    StdVector<FadType> v_fad(v_fad_ptr);     
    StdVector<FadType> u_fad(u_fad_ptr);     
    StdVector<FadType> z_fad(z_fad_ptr);     
    StdVector<FadType> jv_fad(jv_fad_ptr);     

    this->applyJacobian_2AD(jv_fad,v_fad,u_fad,z_fad,tol);
 
    FadType wjv_fad = 0;

    for(int i=0; i<n; ++i) {
        wjv_fad += (*wp)[i]*(*jv_fad_ptr)[i];
    }

    for(int j=0; j<m; ++j) {
        (*ahwvp)[j] = wjv_fad.dx(j);

    }
}
#endif
