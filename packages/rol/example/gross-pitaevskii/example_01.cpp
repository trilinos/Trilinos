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

/** \file   example_01.cpp
    \brief  Minimize the Gross-Pitaevskii functional and demonstrate 
            the effect of choice of function space of the Gradient on
            convergence.

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

    \author Greg von Winckel
    \date   Mon Dec  1 12:55:12 MST 2014
*/

#include "example_01.hpp"

typedef double RealT;

int main(int argc, char **argv) {

    // Set up MPI
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);

    // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
    int iprint     = argc - 1;
    ROL::Ptr<std::ostream> outStream;
    Teuchos::oblackholestream bhs; // outputs nothing
    if (iprint > 0)
        outStream = ROL::makePtrFromRef(std::cout);
    else
        outStream = ROL::makePtrFromRef(bhs);

    int errorFlag = 0;

    Teuchos::ParameterList parlist;
    Teuchos::ParameterList gplist;

    std::string paramfile = "parameters.xml";
    Teuchos::updateParametersFromXmlFile(paramfile,Teuchos::Ptr<Teuchos::ParameterList>(&gplist));
 
    int nx    = gplist.get("Interior Grid Points",100);
    RealT gnl = gplist.get("Nonlinearity Coefficient g",50.0);

    // Grid spacing
    RealT dx = 1.0/(nx+1);

    // Pointer to linspace type vector \f$x_i = \frac{i+1}{n_x+1}\f$ where \f$i=0,\hdots,n_x\f$
    ROL::Ptr<std::vector<RealT> > xi_ptr = ROL::makePtr<std::vector<RealT>>(nx, 0.0);
    
    for(int i=0; i<nx; ++i) {
        (*xi_ptr)[i] = RealT(i+1)/(nx+1);
    }
    
    // Pointer to potential vector (quadratic centered at x=0.5)
    ROL::Ptr<std::vector<RealT> > V_ptr = ROL::makePtr<std::vector<RealT>>(nx, 0.0);
    for(int i=0; i<nx; ++i) {
       (*V_ptr)[i] = 100.0*pow((*xi_ptr)[i]-0.5,2);
    }

    StdVector<RealT> V(V_ptr);
        
    // Iteration Vector (pointer to optimzation vector)
    ROL::Ptr<std::vector<RealT> > psi_ptr = ROL::makePtr<std::vector<RealT>>(nx, 0.0);

       
    // Set Initial Guess (normalized)
    RealT sqrt30 = sqrt(30);

    for (int i=0; i<nx; i++) {
        (*psi_ptr)[i]   = sqrt30*(*xi_ptr)[i]*(1.0-(*xi_ptr)[i]);
    }

    StdVector<RealT> psi(psi_ptr);

    // Constraint value (scalar)  
    ROL::Ptr<std::vector<RealT> > c_ptr = ROL::makePtr<std::vector<RealT>>(1, 0.0);
    StdVector<RealT> c(c_ptr);

    // Lagrange multiplier value (scalar)   
    ROL::Ptr<std::vector<RealT> > lam_ptr = ROL::makePtr<std::vector<RealT>>(1, 0.0);
    StdVector<RealT> lam(lam_ptr);

    // Gradient   
    ROL::Ptr<std::vector<RealT> > g_ptr = ROL::makePtr<std::vector<RealT>>(nx, 0.0);
    StdVector<RealT> g(g_ptr);

    // Instantiate objective function  
    Objective_GrossPitaevskii<RealT> obj(gnl,V);
    
    // Instantiate normalization constraint
    Normalization_Constraint<RealT> constr(nx,dx);

    // Define algorithm.
    std::string stepname = "Composite Step";
    parlist.sublist("Step").sublist(stepname).sublist("Optimality System Solver").set("Nominal Relative Tolerance",1e-4);
    parlist.sublist("Step").sublist(stepname).sublist("Optimality System Solver").set("Fix Tolerance",true);
    parlist.sublist("Step").sublist(stepname).sublist("Tangential Subproblem Solver").set("Iteration Limit",20);
    parlist.sublist("Step").sublist(stepname).sublist("Tangential Subproblem Solver").set("Relative Tolerance",1e-2);
    parlist.sublist("Step").sublist(stepname).set("Output Level",0);
    parlist.sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist.sublist("Status Test").set("Constraint Tolerance",1.e-12);
    parlist.sublist("Status Test").set("Step Tolerance",1.e-14);
    parlist.sublist("Status Test").set("Iteration Limit",100);
    ROL::Algorithm<RealT> algo(stepname, parlist);

    // Run algorithm.
    algo.run(psi, g, lam, c, obj, constr, true, *outStream);


   if(algo.getState()->gnorm>1e-6) {
        errorFlag += 1; 
    }

    if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
    else
        std::cout << "End Result: TEST PASSED\n";



    return 0;

}
