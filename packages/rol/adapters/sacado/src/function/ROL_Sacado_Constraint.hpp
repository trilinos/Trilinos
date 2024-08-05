// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  
           
     

    ROL::Ptr<const vector> xp = dynamic_cast<const SV&>(x).getVector();

    int n = xp->size();

    // Get a pointer to the direction vector
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();

    ROL::Ptr<vector> jvp = dynamic_cast<SV&>(jv).getVector();

    // Create a vector of independent variables
    ROL::Ptr<Fadvector> x_fad_ptr = ROL::makePtr<Fadvector>();
    x_fad_ptr->reserve(n);

    // Initialize constructor for each element
    for(int i=0; i<n; ++i) {
        x_fad_ptr->push_back(FadType(n,i,(*xp)[i])); 
    }

    // Create a vector of independent variables
    ROL::Ptr<Fadvector> c_fad_ptr = ROL::makePtr<Fadvector>();
    c_fad_ptr->reserve(dim_);

    for(int j=0; j<dim_; ++j) {
        c_fad_ptr->push_back(0);  
    }

    StdVector<FadType> x_fad(x_fad_ptr);
    StdVector<FadType> c_fad(c_fad_ptr);

    // Evaluate constraint     
    constr_.value(c_fad,x_fad,tol);

    for(int i=0; i<dim_; ++i) {
        (*jvp)[i] = 0;
        for(int j=0; j<n; ++j) {
            (*jvp)[i] += (*vp)[j]*(*c_fad_ptr)[i].dx(j); 
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
  
           
     

    // Get a pointer to the optimization vector
    ROL::Ptr<const vector> xp = dynamic_cast<const SV&>(x).getVector();

    // Get a pointer to the direction vector
    ROL::Ptr<const vector> up = dynamic_cast<const SV&>(u).getVector();

    ROL::Ptr<vector> ajup = dynamic_cast<SV&>(aju).getVector();

    int n = xp->size();

    // Create a vector of independent variables
    ROL::Ptr<Fadvector> x_fad_ptr = ROL::makePtr<Fadvector>();
    x_fad_ptr->reserve(n);

    // Initialize constructor for each element
    for(int i=0; i<n; ++i) {
        x_fad_ptr->push_back(FadType(n,i,(*xp)[i])); 
    }

    ROL::Ptr<Fadvector> c_fad_ptr = ROL::makePtr<Fadvector>();
    c_fad_ptr->reserve(dim_);
    for(int j=0; j<dim_; ++j) {
        c_fad_ptr->push_back(0);  
    }

    StdVector<FadType> x_fad(x_fad_ptr);
    StdVector<FadType> c_fad(c_fad_ptr);
 
    // Evaluate constraint
    constr_.value(c_fad,x_fad,tol);
    
    FadType udotc = 0;
    
    for(int j=0;j<dim_;++j){ 
        udotc += (*c_fad_ptr)[j]*(*up)[j];
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
  
           
    

    // Get a pointer to the optimization vector
    ROL::Ptr<const vector> xp = dynamic_cast<const SV&>(x).getVector();

    // Get a pointer to the dual constraint vector
    ROL::Ptr<const vector> up = dynamic_cast<const SV&>(u).getVector();

    // Get a pointer to the direction vector
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();

    // Get a pointer to the directional adjoint Hessian 
    ROL::Ptr<vector> ahuvp = dynamic_cast<SV&>(ahuv).getVector();

    // Number of optimization variables
    int n = xp->size();

    // Create a vector of independent variables
    ROL::Ptr<Fadvector> x_fad_ptr = ROL::makePtr<Fadvector>();
    x_fad_ptr->reserve(n);

    // Allocate for directional adjoint Jacobian
    ROL::Ptr<Fadvector> aju_fad_ptr = ROL::makePtr<Fadvector>();
    aju_fad_ptr->reserve(n);

    for(int i=0; i<n; ++i) {
        x_fad_ptr->push_back(FadType(1,(*xp)[i]));

        // Set derivative direction
        (*x_fad_ptr)[i].fastAccessDx(0) = (*vp)[i];     

        aju_fad_ptr->push_back(0);
    }

    // Allocate for constraint vector
    ROL::Ptr<Fadvector> c_fad_ptr = ROL::makePtr<Fadvector>();
    c_fad_ptr->reserve(dim_);

    // Allocate for dual constraint vector
    ROL::Ptr<Fadvector> u_fad_ptr = ROL::makePtr<Fadvector>();
    u_fad_ptr->reserve(dim_);

     for(int j=0; j<dim_; ++j) {
        u_fad_ptr->push_back((*up)[j]);
    }

    StdVector<FadType> x_fad(x_fad_ptr);
    StdVector<FadType> u_fad(u_fad_ptr);
    StdVector<FadType> aju_fad(aju_fad_ptr);

    // Evaluate constraint adjoint Jacobian direction
    this->applyAdjointJacobianAD( aju_fad, u_fad, x_fad, tol);

    for(int i=0; i<n; ++i) {
        (*ahuvp)[i] = (*aju_fad_ptr)[i].dx(0);             
    }
}


}
#endif
