// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "ROL_ConstraintStatusTest.hpp"
#include "ROL_CompositeStep.hpp"

typedef double RealT;

int main(int argc, char **argv) {

    // Set up MPI
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);

    // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
    int iprint     = argc - 1;
    ROL::Ptr<std::ostream> outStream;
    ROL::nullstream bhs; // outputs nothing
    if (iprint > 0)
        outStream = ROL::makePtrFromRef(std::cout);
    else
        outStream = ROL::makePtrFromRef(bhs);

    int errorFlag = 0;

    ROL::ParameterList parlist;

    std::string paramfile = "parameters.xml";
    auto gplist = ROL::getParametersFromXmlFile( paramfile );
 
    int nx    = gplist -> get("Interior Grid Points",100);
    RealT gnl = gplist -> get("Nonlinearity Coefficient g",50.0);

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
    ROL::Ptr<ROL::StatusTest<RealT>>
      status = ROL::makePtr<ROL::ConstraintStatusTest<RealT>>(parlist);
    ROL::Ptr<ROL::Step<RealT>>
      step = ROL::makePtr<ROL::CompositeStep<RealT>>(parlist);
    ROL::Algorithm<RealT> algo(step,status,false);

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
