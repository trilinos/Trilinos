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

           
    

    // Get a pointer to the optimization vector
    ROL::Ptr<const vector> xp = dynamic_cast<const SV&>(x).getVector();

    // Get a pointer to the gradient vector
    ROL::Ptr<vector> gp = dynamic_cast<SV&>(g).getVector();

    int n = xp->size();

    // Create a vector of independent variables
    ROL::Ptr<Fadvector> x_fad_ptr = ROL::makePtr<Fadvector>();
    x_fad_ptr->reserve(n);

    // Initialize constructor for each element
    for(int i=0; i<n; ++i) {
        x_fad_ptr->push_back(FadType(n,i,(*xp)[i]));
    }

    StdVector<FadType> x_fad(x_fad_ptr);

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

           
    

    // Get a pointer to the optimization vector
    ROL::Ptr<const vector> xp = dynamic_cast<const SV&>(x).getVector();

    // Get a pointer to the direction vector
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();

    ROL::Ptr<vector> hvp = dynamic_cast<SV&>(hv).getVector();


    int n = xp->size();

    // Create a vector of independent variables
    ROL::Ptr<Fadvector> x_fad_ptr = ROL::makePtr<Fadvector>();
    x_fad_ptr->reserve(n);

    // Allocate for gradient
    ROL::Ptr<Fadvector> g_fad_ptr = ROL::makePtr<Fadvector>();
    g_fad_ptr->resize(n);

    for(int i=0; i<n; ++i) {
        x_fad_ptr->push_back(FadType(1,(*xp)[i]));
    }

    // Set directional derivative
    for(int i=0; i<n; ++i) {
        (*x_fad_ptr)[i].fastAccessDx(0) = (*vp)[i];
    }

    StdVector<FadType> x_fad(x_fad_ptr);
    StdVector<FadType> g_fad(g_fad_ptr);

    this->gradientAD(g_fad,x_fad,tol);

    for(int i=0; i<n; ++i) {
        (*hvp)[i] = (*g_fad_ptr)[i].dx(0);
    }
}


#endif
