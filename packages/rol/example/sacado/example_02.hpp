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

#include <iostream>

#include "Sacado.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_EqualityConstraint.hpp"
#include "ROL_CompositeStepSQP.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

using namespace ROL;

// ---- Begin Sacado functions ----

/** \brief Sacado-accsible objective function:
           f(x) = exp(x1*x2*x3*x4*x5) + 0.5*(x1^3+x2^3+1)^2
*/
template <typename ScalarT>
ScalarT objfun(const Vector<ScalarT>& x) {
    Teuchos::RCP<const std::vector<ScalarT> > xp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(x))).getVector();
    int n = xp->size();

    // std::cout << "\033[1;37m xp->size() = " << n << "\033[1;38m" << std::endl; 
    TEUCHOS_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (objfun): "
                                                                 "Primal vector x must be of length 5.");

    ScalarT x1 = (*xp)[0];
    ScalarT x2 = (*xp)[1];
    ScalarT x3 = (*xp)[2];
    ScalarT x4 = (*xp)[3];
    ScalarT x5 = (*xp)[4];

    ScalarT J = exp(x1*x2*x3*x4*x5) - 0.5 * pow( (pow(x1,3)+pow(x2,3)+1.0), 2);
    return J;  
}


/** \brief A Sacado-accesible function which computes the gradient of the objective function 
    using automatic differentiation.
    @param[in]  x is the optimization vector
    @param[out] J is the value of the objective function
    @param[out] g is the gradient of the objective function    
*/
template <typename ScalarT>
void objgrad(const Vector<ScalarT>& x, ScalarT& J, Vector<ScalarT>& g) {

    // Data type which supports automatic differentiation 
    typedef Sacado::Fad::DFad<ScalarT> FadType;

    // Get a pointer to the optimization vector
    Teuchos::RCP<const std::vector<ScalarT> > xp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(x))).getVector();

    // Get a pointer to the gradient vector
    Teuchos::RCP<std::vector<ScalarT> > gp =
        Teuchos::rcp_const_cast<std::vector<ScalarT> > ((Teuchos::dyn_cast<StdVector<ScalarT> > (g)).getVector());
    
    int n = xp->size();

    TEUCHOS_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (objgrad): "
                                                                 "Primal vector x must be of length 5.");


    // Create a vector of independent variables
    Teuchos::RCP<std::vector<FadType> > x_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    x_fad_rcp->reserve(n);
   
    // Initialize constructor for each element
    for(int i=0; i<n; ++i) {
        x_fad_rcp->push_back(FadType(n,i,(*xp)[i])); 
    }

    StdVector<FadType> x_fad(x_fad_rcp);

    // AD access to objective function
    FadType J_fad = objfun(x_fad);

    // Evaluate objective function
    J = J_fad.val();

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
template <typename ScalarT>
void applyHessian(const Vector<ScalarT>& x, const Vector<ScalarT>& v,
             ScalarT& J, Vector<ScalarT>& g, Vector<ScalarT>& hv) {
  
    // Get a pointer to the optimization vector 
    Teuchos::RCP<const std::vector<ScalarT> > xp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(x))).getVector();

    // Get a pointer to the gradient vector
    Teuchos::RCP<std::vector<ScalarT> > gp =
        Teuchos::rcp_const_cast<std::vector<ScalarT> > ((Teuchos::dyn_cast<StdVector<ScalarT> > (g)).getVector());

    // Get a pointer to the direction vector
    Teuchos::RCP<const std::vector<ScalarT> > vp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(v))).getVector();

    Teuchos::RCP<std::vector<ScalarT> > hvp =
        Teuchos::rcp_const_cast<std::vector<ScalarT> >((Teuchos::dyn_cast<StdVector<ScalarT> >(hv)).getVector());

    typedef Sacado::Fad::SFad<ScalarT,1> FadType;

    int n = xp->size();

    TEUCHOS_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (applyHessian): "
                                                                   "Primal vector x must be of length 5."); 

    TEUCHOS_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (applyHessian): "
                                                                   "Input vector v must be of length 5.");    


    FadType J_fad; 

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

    objgrad(x_fad,J_fad,g_fad);

    J = J_fad.val();

    for(int i=0; i<n; ++i) {
        (*gp)[i] = (*g_fad_rcp)[i].val();
        (*hvp)[i] = (*g_fad_rcp)[i].dx(0);            
    }
}


/** \brief Sacado-accsible equality constraint function with c_i(x) = 0, where :
           c_1(x) = x1^2+x2^2+x3^2+x4^2+x5^2 - 10
           c_2(x) = x2*x3-5*x4*x5
           c_3(x) = x1^3 + x2^3 + 1
*/
template <typename ScalarT>
void confun(const Vector<ScalarT>& x, Vector<ScalarT>& c) {
    Teuchos::RCP<std::vector<ScalarT> > cp = 
        Teuchos::rcp_const_cast<std::vector<ScalarT> >((Teuchos::dyn_cast<StdVector<ScalarT> >(c)).getVector());
    Teuchos::RCP<const std::vector<ScalarT> > xp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(x))).getVector();

    int n = xp->size();
    int m = cp->size();

    TEUCHOS_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (confun): "
                                                                 "Primal vector x must be of length 5.");
    TEUCHOS_TEST_FOR_EXCEPTION( (m != 3), std::invalid_argument, ">>> ERROR (confun): "
                                                                 "Constraint vector c must be of length 3.");

    ScalarT x1 = (*xp)[0];
    ScalarT x2 = (*xp)[1];
    ScalarT x3 = (*xp)[2];
    ScalarT x4 = (*xp)[3];
    ScalarT x5 = (*xp)[4];

    (*cp)[0] = x1*x1+x2*x2+x3*x3+x4*x4+x5*x5 - 10.0;
    (*cp)[1] = x2*x3 - 5.0*x4*x5;
    (*cp)[2] = x1*x1*x1 + x2*x2*x2 + 1.0;

}

/** 
    \brief Compute the action of the constraint Jacobian on a direction vector 
    \f[ \nabla\mathbf{c}(\mathbf{x})\mathbf{v} = \sum\limits_{i=1}^n 
         \frac{\partial \mathbf{c}(\mathbf{x})}{\partial x_i} v_i\f]
*/
template <typename ScalarT>
void conjac( const Vector<ScalarT> &x,  const Vector<ScalarT> &v, 
                          Vector<ScalarT> &c,        Vector<ScalarT> &jv) {

    // Data type which supports automatic differentiation 
    typedef Sacado::Fad::DFad<ScalarT> FadType;

    Teuchos::RCP<std::vector<ScalarT> > cp = 
        Teuchos::rcp_const_cast<std::vector<ScalarT> >((Teuchos::dyn_cast<StdVector<ScalarT> >(c)).getVector());

    Teuchos::RCP<const std::vector<ScalarT> > xp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(x))).getVector();

    // Get a pointer to the direction vector
    Teuchos::RCP<const std::vector<ScalarT> > vp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(v))).getVector();

    Teuchos::RCP<std::vector<ScalarT> > jvp = 
        Teuchos::rcp_const_cast<std::vector<ScalarT> >((Teuchos::dyn_cast<StdVector<ScalarT> >(jv)).getVector());

    int n = xp->size();
    TEUCHOS_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (conjac): "
                                                                 "Primal vector x must be of length 5.");
    n = vp->size();
    TEUCHOS_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (conjac): "
                                                                 "Input vector v must be of length 5.");
    int m = cp->size();
    TEUCHOS_TEST_FOR_EXCEPTION( (m != 3), std::invalid_argument, ">>> ERROR (conjac): "
                                                                 "Constraint vector c must be of length 3.");
    m = jvp->size();
    TEUCHOS_TEST_FOR_EXCEPTION( (m != 3), std::invalid_argument, ">>> ERROR (conjac): "
                                                                 "Output vector jv must be of length 3.");

    // Create a vector of independent variables
    Teuchos::RCP<std::vector<FadType> > x_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    x_fad_rcp->reserve(n);

    // Initialize constructor for each element
    for(int i=0; i<n; ++i) {
        x_fad_rcp->push_back(FadType(n,i,(*xp)[i])); 
    }

    // Create a vector of independent variables
    Teuchos::RCP<std::vector<FadType> > c_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    c_fad_rcp->reserve(m);

    for(int j=0; j<m; ++j) {
        c_fad_rcp->push_back(0);  
    }

    StdVector<FadType> x_fad(x_fad_rcp);
    StdVector<FadType> c_fad(c_fad_rcp);

    // Evaluate constraint     
    confun(x_fad,c_fad);

    for(int i=0; i<m; ++i) {
        (*jvp)[i] = 0;
        for(int j=0; j<n; ++j) {
            (*jvp)[i] += (*vp)[j]*(*c_fad_rcp)[i].dx(j); 
        }   
    }       
}  

/** 
    \brief Compute the action of the transpose of the constraint Jacobian on a direction vector 
    \f[ [\nabla \mathbf{c}(\mathbf{x})]^\top\mathbf{u} =
              \left(\frac{\partial [ \mathbf{c}(\mathbf{x})]^\top\mathbf{u}}{\partial x_1} \cdots 
              \frac{\partial [ \mathbf{c}(\mathbf{x})]^\top\mathbf{u}}{\partial x_n} \right) \f]
*/
template <typename ScalarT>
void conadjjac( const Vector<ScalarT> &x, const Vector<ScalarT> &u, 
                      Vector<ScalarT> &c,       Vector<ScalarT> &aju) {

    // Data type which supports automatic differentiation 
    typedef Sacado::Fad::DFad<ScalarT> FadType;

    // Get a pointer to the constraint vector
    Teuchos::RCP<std::vector<ScalarT> > cp = 
        Teuchos::rcp_const_cast<std::vector<ScalarT> >((Teuchos::dyn_cast<StdVector<ScalarT> >(c)).getVector());

    // Get a pointer to the optimization vector
    Teuchos::RCP<const std::vector<ScalarT> > xp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(x))).getVector();

    // Get a pointer to the direction vector
    Teuchos::RCP<const std::vector<ScalarT> > up =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(u))).getVector();

    Teuchos::RCP<std::vector<ScalarT> > ajup = 
        Teuchos::rcp_const_cast<std::vector<ScalarT> >((Teuchos::dyn_cast<StdVector<ScalarT> >(aju)).getVector());

    int n = xp->size();
    TEUCHOS_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (conadjjac): "
                                                                 "Primal vector x must be of length 5.");
    int m = up->size();

    TEUCHOS_TEST_FOR_EXCEPTION( (m != 3), std::invalid_argument, ">>> ERROR (conadjjac): "
                                                                 "Constraint dual vector u must be of length 3.");
    m = cp->size();
    TEUCHOS_TEST_FOR_EXCEPTION( (m != 3), std::invalid_argument, ">>> ERROR (conadjjac): "
                                                                 "Constraint vector c must be of length 3.");
    n = ajup->size();
    TEUCHOS_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (conadjjac): "
                                                                 "Output vector aju must be of length 5.");

    // Create a vector of independent variables
    Teuchos::RCP<std::vector<FadType> > x_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    x_fad_rcp->reserve(n);

    // Initialize constructor for each element
    for(int i=0; i<n; ++i) {
        x_fad_rcp->push_back(FadType(n,i,(*xp)[i])); 
    }

    Teuchos::RCP<std::vector<FadType> > c_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    c_fad_rcp->reserve(m);
    for(int j=0; j<m; ++j) {
        c_fad_rcp->push_back(0);  
    }

    StdVector<FadType> x_fad(x_fad_rcp);
    StdVector<FadType> c_fad(c_fad_rcp);
 
    // Evaluate constraint
    confun(x_fad,c_fad);

    
    FadType udotc = 0;
    
    for(int j=0;j<m;++j){ 
        (*cp)[j] = (*c_fad_rcp)[j].val();
        udotc += (*c_fad_rcp)[j]*(*up)[j];
    } 

    for(int i=0;i<n;++i){
        (*ajup)[i] = udotc.dx(i);
    } 
}


template <typename ScalarT>
void conadjhess( const Vector<ScalarT> &x,   const Vector<ScalarT> &u,
                 const Vector<ScalarT> &v,         Vector<ScalarT> &c,
                       Vector<ScalarT> &aju,       Vector<ScalarT> &ahuv ) {
      
     // Data type which supports automatic differentiation 
    typedef Sacado::Fad::SFad<ScalarT,1> FadType;

    // Get a pointer to the optimization vector
    Teuchos::RCP<const std::vector<ScalarT> > xp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(x))).getVector();

    // Get a pointer to the dual constraint vector
    Teuchos::RCP<const std::vector<ScalarT> > up =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(u))).getVector();

    // Get a pointer to the direction vector
    Teuchos::RCP<const std::vector<ScalarT> > vp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(v))).getVector();

    // Get a pointer to the constraint vector
    Teuchos::RCP<std::vector<ScalarT> > cp = 
        Teuchos::rcp_const_cast<std::vector<ScalarT> >((Teuchos::dyn_cast<StdVector<ScalarT> >(c)).getVector());

    // Get a pointer to the directional adjoint Jacobian 
    Teuchos::RCP<std::vector<ScalarT> > ajup = 
        Teuchos::rcp_const_cast<std::vector<ScalarT> >((Teuchos::dyn_cast<StdVector<ScalarT> >(aju)).getVector());

    // Get a pointer to the directional adjoint Hessian 
    Teuchos::RCP<std::vector<ScalarT> > ahuvp = 
        Teuchos::rcp_const_cast<std::vector<ScalarT> >((Teuchos::dyn_cast<StdVector<ScalarT> >(ahuv)).getVector());

    // Number of optimization variables
    int n = xp->size();

    TEUCHOS_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (conadjhess): "
                                                                 "Primal vector x must be of length 5.");

    n = vp->size();

    TEUCHOS_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (conadjhess): "
                                                                 "Direction vector v must be of length 5.");
    int m = cp->size();

    TEUCHOS_TEST_FOR_EXCEPTION( (m != 3), std::invalid_argument, ">>> ERROR (conadjhess): "
                                                                 "Constraint vector c must be of length 3.");
    m = up->size();


    TEUCHOS_TEST_FOR_EXCEPTION( (m != 3), std::invalid_argument, ">>> ERROR (conadjhess): "
                                                                 "Constraint dual vector u must be of length 3.");
    n = ajup->size();
    TEUCHOS_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (conadjhess): "
                                                                 "Constraint adjoint Jacobian direction aju must be of length 5.");
    n = ahuvp->size();

    TEUCHOS_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (conadjhess): "
                                                                 "Output vector ahuv must be of length 5.");
    
    // Create a vector of independent variables
    Teuchos::RCP<std::vector<FadType> > x_fad_rcp =  Teuchos::rcp( new std::vector<FadType> );
    x_fad_rcp->reserve(n);

    for(int i=0; i<n; ++i) {
        x_fad_rcp->push_back(FadType(1,(*xp)[i]));
    }

    // Set directional derivative
    for(int i=0; i<n; ++i) {
        (*x_fad_rcp)[i].fastAccessDx(0) = (*vp)[i];     
    }

    // Allocate for constraint vector
    Teuchos::RCP<std::vector<FadType> > c_fad_rcp =  Teuchos::rcp( new std::vector<FadType> );
    c_fad_rcp->reserve(m);
    for(int j=0; j<m; ++j) {
        c_fad_rcp->push_back(0);  
    }

    // Allocate for dual constraint vector
    Teuchos::RCP<std::vector<FadType> > u_fad_rcp =  Teuchos::rcp( new std::vector<FadType> );
    u_fad_rcp->reserve(m);
   

    // Allocate for directional adjoint Jacobian
    Teuchos::RCP<std::vector<FadType> > aju_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    aju_fad_rcp->reserve(n);

    for(int i=0;i<n;++i) {
        aju_fad_rcp->push_back(0);
    } 

    for(int j=0;j<m;++j) {
        (*u_fad_rcp).push_back((*up)[j]);
    } 

    StdVector<FadType> x_fad(x_fad_rcp);
    StdVector<FadType> u_fad(u_fad_rcp);
    StdVector<FadType> c_fad(c_fad_rcp);
    StdVector<FadType> aju_fad(aju_fad_rcp);

    // Evaluate constraint adjoint Jacobian direction
    conadjjac( x_fad, u_fad, c_fad, aju_fad);

    for(int i=0; i<n; ++i) {
        (*ajup)[i]  = (*aju_fad_rcp)[i].val();
        (*ahuvp)[i] = (*aju_fad_rcp)[i].dx(0);             
    }
}

// ---- End Sacado functions ----




/** \brief Objective Class */
template<class Real>
class Objective_Sacado : public Objective<Real> {

    public:
    Objective_Sacado(void) { }

    /** \brief Evaluate the objective function at x */
    Real value( const Vector<Real> &x, Real &tol ) {

        Teuchos::RCP<const std::vector<Real> > xp =
            (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();

        Real J = objfun(x);

        return J;
    }

    /** \brief Evaluate the gradient at x */
    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {

        Teuchos::RCP<const std::vector<Real> > xp =
            (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();

        Teuchos::RCP<std::vector<Real> > gp =
            Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(g)).getVector());

        Real J;

        objgrad(x,J,g);

    }

    /** \brief Compute the action of the Hessian evaluated at x on a vector v */
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
        Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());

        int n = xp->size();

        Real J;

        Teuchos::RCP<std::vector<Real> > gp = Teuchos::rcp( new std::vector<Real>(n,0.0) );
        StdVector<Real> g(gp);

        applyHessian(x,v,J,g,hv);
    }
};


template<class Real>
class Constraint_Sacado : public EqualityConstraint<Real> {
    private:
        //! \var m_ Number of constraint equations
        int m_;

        //! \var n_ Number of optimization variables
        int n_;

    public:
        Constraint_Sacado(int m,int n) : m_(m), n_(n) {        
   }

    void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {

        Teuchos::RCP<std::vector<Real> > cp = 
            Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(c)).getVector());

        confun(x,c);

    }


    void applyJacobian( Vector<Real> &jv, const Vector<Real> &v, 
                        const Vector<Real> &x, Real &tol) {

        Teuchos::RCP<std::vector<Real> > cp = Teuchos::rcp( new std::vector<Real>(m_,0.0) );
        StdVector<Real> c(cp);

        conjac(x,v,c,jv);

    }


    void applyAdjointJacobian( Vector<Real> &aju, const Vector<Real> &u, 
                               const Vector<Real> &x, Real &tol) {

        Teuchos::RCP<std::vector<Real> > cp = Teuchos::rcp( new std::vector<Real>(m_,0.0) );
        StdVector<Real> c(cp);

        conadjjac(x,u,c,aju);
 
    } 


    void applyAdjointHessian( Vector<Real> &ahuv, const Vector<Real> &u, 
                              const Vector<Real> &v, const Vector<Real> &x, Real &tol) {

        Teuchos::RCP<std::vector<Real> > cp = Teuchos::rcp( new std::vector<Real>(m_,0.0) );
        StdVector<Real> c(cp);

        Teuchos::RCP<std::vector<Real> > ajup = Teuchos::rcp( new std::vector<Real>(n_,0.0) );
        StdVector<Real> aju(ajup);

        conadjhess(x,u,v,c,aju,ahuv);
   
    } 
};

