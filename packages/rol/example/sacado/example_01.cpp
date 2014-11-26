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
    \brief An example combining ROL and Sacado to mimize the Zakharov 
           function. The gradient and the action of the Hessian on a given
           vector are computed by Sacado using automatic differentiation.    
                    
           This implementation is far from optimal as vectors of AD type
           are being created repeatedly. A more efficient implementation 
           would encapsulate the functions Zakharov, objgrad, and hessVec
           in an object so that the AD vectors can be instantiated once. 

    \author Created by G. von Winckel
**/

#include <iostream>
#include <iomanip>

#include "Sacado.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_Algorithm.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"


using namespace ROL;


// The function to differentiate
template <typename ScalarT>
ScalarT Zakharov(const Vector<ScalarT>& x) {
    /** 
    Zakharov objective function
    \f[f(\mathbf{x}) = \mathbf{x}^\top\mathbf{x} 
                     + \frac{1}{4}(\mathbf{k}^\top \mathbf{x})^2 +
                       \frac{1}{16}(\mathbf{k}^\top \mathbf{x})^4 \f]
    Where \f$\mathbf{k}=(1,\cdots,n)\f$
    **/

    Teuchos::RCP<const std::vector<ScalarT> > xp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(x))).getVector();

    int n = xp->size();
    ScalarT xdotx = 0;
    ScalarT kdotx = 0;
    ScalarT J = 0;
   
    // Compute dot products 
    for(int i=0; i<n; ++i) {
        xdotx += pow((*xp)[i],2);       // (k,x)
        kdotx += double(i+1)*(*xp)[i];  // (x,x)
    }

    // Objective function
    J = xdotx + pow(kdotx,2)/4.0 + pow(kdotx,4)/16.0;
    
    return J;
}



/* Evaluate the objective function and its gradient at x */
template <typename ScalarT>
void objgrad(const Vector<ScalarT>& x, ScalarT& J, Vector<ScalarT>& g) {

    typedef Sacado::Fad::DFad<ScalarT> FadType;

    Teuchos::RCP<const std::vector<ScalarT> > xp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(x))).getVector();
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
    FadType J_fad = Zakharov(x_fad);

    // Evaluate objective function
    J = J_fad.val();

    // Evaluate gradient
    for(int i=0; i<n; ++i) {
        (*gp)[i] = J_fad.dx(i);
    }
}

/*
Evaluate the objection function and its gradient at x and evaluate
the action of the Hessian at x on a vector v to produce a direction Hv.
Essentially, we are computing the v-directional derivative of the gradient.
*/

template <typename ScalarT>
void applyHessian(const Vector<ScalarT>& x, const Vector<ScalarT>& v,
             ScalarT& J, Vector<ScalarT>& g, Vector<ScalarT>& hv) {
   
    Teuchos::RCP<const std::vector<ScalarT> > xp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(x))).getVector();
    Teuchos::RCP<std::vector<ScalarT> > gp =
        Teuchos::rcp_const_cast<std::vector<ScalarT> > ((Teuchos::dyn_cast<StdVector<ScalarT> > (g)).getVector());
    Teuchos::RCP<const std::vector<ScalarT> > vp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(v))).getVector();
    Teuchos::RCP<std::vector<ScalarT> > hvp =
        Teuchos::rcp_const_cast<std::vector<ScalarT> >((Teuchos::dyn_cast<StdVector<ScalarT> >(hv)).getVector());

    typedef Sacado::Fad::SFad<ScalarT,1> FadType;

    int n = xp->size();
   
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


/** Objective Class Definition begin **/


template<class Real>
class Objective_Sacado : public Objective<Real> {
    private:

    public:
    Objective_Sacado(void) {}
 
    /* Evaluate the objective function at x */
    Real value( const Vector<Real> &x, Real &tol ) {
      StdVector<Real> & ex =
        Teuchos::dyn_cast<StdVector<Real> >(const_cast <Vector<Real> &>(x));
      Teuchos::RCP<const std::vector<Real> > xp = ex.getVector();

      Real J = Zakharov(x);

      return J;
    }

    /* Evaluate the gradient at x */
    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
        Teuchos::RCP<const std::vector<Real> > xp =
            (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
        Teuchos::RCP<std::vector<Real> > gp =
            Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(g)).getVector());

        Real J;
        objgrad(x,J,g);

    }

    /* Compute the action of the Hessian evaluated at x on a vector v */
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

/** Objective Class Definition end **/



typedef double RealT;

int main(int argc, char **argv)
{
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);

    // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
    int iprint     = argc - 1;
    Teuchos::RCP<std::ostream> outStream;
    Teuchos::oblackholestream bhs; // outputs nothing
    if (iprint > 0)
        outStream = Teuchos::rcp(&std::cout, false);
    else
        outStream = Teuchos::rcp(&bhs, false);

    int errorFlag  = 0;

    // *** Example body.

    try {

        Objective_Sacado<RealT> obj;
    
        int dim = 10; // Set problem dimension. 

        // Load optimizer parameters form XML file
        Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp(new Teuchos::ParameterList());
        std::string paramfile = "parameters.xml";
        Teuchos::updateParametersFromXmlFile(paramfile,Teuchos::Ptr<Teuchos::ParameterList>(&*parlist));

        // Define Step
        LineSearchStep<RealT> step(*parlist);

        // Define Status Test
        RealT gtol  = 1e-12;  // norm of gradient tolerance
        RealT stol  = 1e-14;  // norm of step tolerance
        int   maxit = 100;    // maximum number of iterations
        StatusTest<RealT> status(gtol, stol, maxit);    

        // Define Algorithm
        DefaultAlgorithm<RealT> algo(step,status,false);

        // Iteration Vector
        Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
        // Set Initial Guess
        for (int i=0; i<dim; i++) {
            (*x_rcp)[i]   = 2;
        }

        StdVector<RealT> x(x_rcp);

        // Run Algorithm
        std::vector<std::string> output = algo.run(x, obj, false);
        for ( unsigned i = 0; i < output.size(); i++ ) {
            std::cout << output[i];
        }

        // Get True Solution
        Teuchos::RCP<std::vector<RealT> > xtrue_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
        StdVector<RealT> xtrue(xtrue_rcp);

        
        // Compute Error
        x.axpy(-1.0, xtrue);
        RealT abserr = x.norm();
        *outStream << std::scientific << "\n   Absolute Error: " << abserr;
        if ( abserr > sqrt(ROL_EPSILON) ) {
            errorFlag += 1;
        }
    }
    catch (std::logic_error err) {
        *outStream << err.what() << "\n";
        errorFlag = -1000;
    }; // end try

    if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
    else
        std::cout << "End Result: TEST PASSED\n";

    return 0;

}
