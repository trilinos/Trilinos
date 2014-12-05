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

/** \file   example_02.cpp
    \brief  Minimize the Gross-Pitaevskii functional and demonstrate 
            the effect of choice of function space of the Gradient on
            convergence. In this version we implement the correct Sobolev
            inner products and Reisz mapping using a finite difference 
            approximation of the derivative terms.          
                

    \details Minimize the one-dimensional Gross-Pitaevskii (GP) energy 
             functional
             \f[ J[\psi] = \int \frac{1}{2} |\nabla\psi|^2 + V(x)|\psi|^2 
                           +g|\psi|^4 \,\mathrm{d}x \f]
             Subject to the equality constraint that the particle density be
             normalized. 
             \f[ e(\psi) = \int |\psi|^2\,\mathrm{d}x - 1 = 0 \f]
             For simplicity, we will assume the wavefunction \f$\psi\f$ to 
             be real-valued, the potential function \f$ V(x)\geq 0\f$,
             the computational domain is the interval \f$[0,1]\f$, and that
             \f$\psi(0)=\psi(1)=0\f$. We also discretize the problem using
             second-order centered finite differences on a uniform grid. 

             \f[
             \psi''(x_i) \approx = \frac{\psi(x_{i-1})-2\psi(x_i)+\psi(x_{i+1})}{\Delta x^2}
             \f]

             The gradient with respect to the \f$L^2\f$ inner product is actually 
             an element of \f$H^{-1}\f$, so if we search in this direction, we
             are actually looking for a solution \f$\psi\f$ in a larger space than
             we should. If we compute the gradient with respect to the \f$H^1\f$ 
             inner product, by solving a Poisson equation, we search in the right space
             and the optimizer converges much faster than example_01.cpp which  
             does not do this.  
             

    \author Greg von Winckel
    \date   Wed Dec  3 16:40:45 MST 2014
*/

#include<algorithm>
#include<string>
#include"example_02.hpp"

typedef double RealT;

int main(int argc, char **argv) {

    
    // Set up MPI
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);

    // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
    int iprint     = argc - 1;
    Teuchos::RCP<std::ostream> outStream;
    Teuchos::oblackholestream bhs; // outputs nothing
    if (iprint > 0)
        outStream = Teuchos::rcp(&std::cout, false);
    else
        outStream = Teuchos::rcp(&bhs, false);

    int errorFlag = 0;
 

    Teuchos::ParameterList parlist;
    std::string paramfile = "parameters.xml";
    Teuchos::updateParametersFromXmlFile(paramfile,Teuchos::Ptr<Teuchos::ParameterList>(&parlist));
       
    int    nx         = parlist.get("Interior Grid Points",100);
    double gnl        = parlist.get("Nonlinearity Coefficient g",50.0);
    bool   exactsolve = parlist.get("Solve Exact Augmented System",false);

    // Command line option to override parameters.xml for solving the exact augmented system
    if(argc > 1) {
        std::string input = argv[1];
        std::transform(input.begin(), input.end(), input.begin(), ::tolower);  
        if(input=="exactsolve") {
            exactsolve = true;
        }
    }

 
    // Grid spacing
    RealT dx = 1.0/(nx+1);

    // Finite difference class
    FiniteDifference<RealT> *fd = new FiniteDifference<RealT>(nx,dx);

    // Pointer to linspace type vector \f$x_i = \frac{i+1}{n_x+1}\f$ where \f$i=0,\hdots,n_x\f$
    Teuchos::RCP<std::vector<RealT> > xi_rcp = Teuchos::rcp( new std::vector<RealT> (nx, 0.0) );
    
    for(int i=0; i<nx; ++i) {
        (*xi_rcp)[i] = RealT(i+1)/(nx+1);
    }
    
    // Pointer to potential vector (quadratic centered at x=0.5)
    Teuchos::RCP<std::vector<RealT> > V_rcp = Teuchos::rcp( new std::vector<RealT> (nx, 0.0) );
    for(int i=0; i<nx; ++i) {
       (*V_rcp)[i] = 100.0*pow((*xi_rcp)[i]-0.5,2);
    }

    StdVector<RealT> V(V_rcp);
        
    // Iteration Vector (pointer to optimzation vector)
    Teuchos::RCP<std::vector<RealT> > psi_rcp = Teuchos::rcp( new std::vector<RealT> (nx, 0.0) );
    OptStdVector<RealT> psi(psi_rcp,fd);
       
    // Set Initial Guess (normalized)
    double sqrt30 = sqrt(30);

    for (int i=0; i<nx; i++) {
        (*psi_rcp)[i]   = sqrt30*(*xi_rcp)[i]*(1.0-(*xi_rcp)[i]);
    }


    // Equality constraint value (scalar)  
    Teuchos::RCP<std::vector<RealT> > c_rcp = Teuchos::rcp( new std::vector<RealT> (1, 0.0) );
    ConStdVector<RealT> c(c_rcp);

    // Lagrange multiplier value (scalar)   
    Teuchos::RCP<std::vector<RealT> > lam_rcp = Teuchos::rcp( new std::vector<RealT> (1, 0.0) );
    ConDualStdVector<RealT> lam(lam_rcp);

    // Gradient   
    Teuchos::RCP<std::vector<RealT> > g_rcp = Teuchos::rcp( new std::vector<RealT> (nx, 0.0) );
    OptDualStdVector<RealT> g(g_rcp,fd);

    // Instantiate objective function  
    Objective_GrossPitaevskii<RealT,OptStdVector<RealT>,OptDualStdVector<RealT> > obj(gnl,V,fd);
    
    // Instantiate normalization constraint
    Normalization_Constraint<RealT,OptStdVector<RealT>,OptDualStdVector<RealT>, 
             ConStdVector<RealT>,ConDualStdVector<RealT> > constr(nx,dx,fd,exactsolve);



    // Define Step
    parlist.set("Nominal SQP Optimality Solver Tolerance", 1.e-4);
    parlist.set("Maximum Number of Krylov Iterations",80);
    parlist.set("Absolute Krylov Tolerance",1e-4);

    CompositeStepSQP<RealT> step(parlist);


    // Define Status Test
    RealT gtol  = 1e-12;  // norm of gradient tolerance
    RealT ctol  = 1e-12;  // norm of constraint tolerance
    RealT stol  = 1e-14;  // norm of step tolerance
    int   maxit = 100;    // maximum number of iterations
    StatusTestSQP<RealT> status(gtol, ctol, stol, maxit);    

    // Define Algorithm
    DefaultAlgorithm<RealT> algo(step,status,false);

    // Run Algorithm
    std::vector<std::string> output = algo.run(psi, g, lam, c, obj, constr, false);
  
    for ( unsigned i = 0; i < output.size(); i++ ) {
      *outStream << output[i];
    }

    if(algo.getState()->gnorm>1e-6) {
        errorFlag += 1; 
    }

    if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
    else
        std::cout << "End Result: TEST PASSED\n";

    delete fd;

    return 0;
}
