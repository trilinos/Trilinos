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

/** \file   example_03.cpp
    \brief  Minimize the Gross-Pitaevskii functional and demonstrate 
            the effect of choice of function space of the Gradient on
            convergence. In this version we implement the option to use 
            correct Sobolev inner products and Riesz mapping using nodal 
            Galerkin methods.          
                
    \details Minimize the one-dimensional Gross-Pitaevskii (GP) energy 
             functional
             \f[ J[\psi] = \int \frac{1}{2} |\nabla\psi|^2 + V(x)|\psi|^2 
                           +g|\psi|^4 \,\mathrm{d}x \f]
             Subject to the equality constraint that the particle density be
             normalized. 
             \f[ e(\psi) = \int |\psi|^2\,\mathrm{d}x - 1 = 0 \f]
             For simplicity, we will assume the wavefunction \f$\psi\f$ to 
             be real-valued, the potential function \f$ V(x)\geq 0\f$,
             the computational domain is the interval \f$[-1,1]\f$, and that
             \f$\psi(-1)=\psi(1)=0\f$. We descretize using the nodal Galerkin method
             \f[\psi(x)\approx \sum\limits_{k=1}^{n-1} \psi(x_k) \ell_k(x)\f].
               
    \author Greg von Winckel
    \date   Mon Dec  8 11:01:24 MST 2014
*/


#include "example_03.hpp"


typedef double RealT;

int main(int argc, char* argv[]) { 

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
       
    int    ni         = parlist.get("Interior Grid Points",20);
    double gnl        = parlist.get("Nonlinearity Coefficient g",50.0);
    bool   exactsolve = parlist.get("Solve Exact Augmented System",false);
    bool   useRiesz   = parlist.get("Use Riesz Map",true);

    // Number of quadrature points
    int nq = 2*ni;

    Teuchos::RCP<Teuchos::LAPACK<int,RealT> >lapack = 
        Teuchos::rcp( new Teuchos::LAPACK<int,RealT>() );
    Teuchos::RCP<NodalBasis<RealT> > nb = Teuchos::rcp(new NodalBasis<RealT>(lapack,ni+2,nq));

    // Exclude end point interpolants to impose Dirichlet condition 
    std::vector<RealT> L(ni*nq,0);
    std::copy(nb->L_.begin()+nq,nb->L_.end()-nq,L.begin()); 
    std::vector<RealT> Lp(ni*nq,0);
    std::copy(nb->Lp_.begin()+nq,nb->Lp_.end()-nq,Lp.begin()); 

    // Quadrature points 
    std::vector<RealT> x(nb->xq_); 

    // Quadrature weights
    std::vector<RealT> w(nb->wq_); 

    // Mass matrix
    Teuchos::RCP<InnerProductMatrix<RealT> > mass =    
        Teuchos::rcp( new InnerProductMatrixSolver<RealT>(lapack,L,L,w,1) );
    //Teuchos::RCP<InnerProductMatrix<RealT> > mass =    
    //    Teuchos::null;
/*
    // Kinetic energy matrix
    Teuchos::RCP<InnerProductMatrix<RealT> > kinetic = 
        Teuchos::rcp( new InnerProductMatrixSolver<RealT>(lapack,Lp,Lp,w,1) );

    // Nonlinear (quartic) term. Set \f$\psi\f$ later
    Teuchos::RCP<InnerProductMatrix<RealT> > nonlinear = 
        Teuchos::rcp( new InnerProductMatrix<RealT>(L,L,w,1) ); 

    // Confinement Potential 
    std::vector<RealT> v(nq,0);   
    for(int i=0;i<nq;++i){
        v[i] = 100*x[i]*x[i];
    }

    Teuchos::RCP<InnerProductMatrix<RealT> > potential = 
        Teuchos::rcp( new InnerProductMatrix<RealT>(L,L,w,v) ); 

    // Iteration Vector (pointer to optimzation vector)
    Teuchos::RCP<std::vector<RealT> > psi_rcp = Teuchos::rcp( new std::vector<RealT> (ni, 0.0) );
    
    // Normalized initial guess 
    for(int i=0;i<ni;++i){
        (*psi_rcp)[i]=(1+nb->xi_[i+1])*(1-nb->xi_[i+1])*sqrt(15.0/16.0);
    } 
    OptStdVector<RealT> psi(psi_rcp,useRiesz,kinetic);
  
    // Equality constraint value (scalar)  
    Teuchos::RCP<std::vector<RealT> > c_rcp = Teuchos::rcp( new std::vector<RealT> (1, 0.0) );
    ConStdVector<RealT> c(c_rcp,useRiesz,mass);

    // Lagrange multiplier value (scalar)   
    Teuchos::RCP<std::vector<RealT> > lam_rcp = Teuchos::rcp( new std::vector<RealT> (1, 0.0) );
    ConDualStdVector<RealT> lam(lam_rcp,useRiesz,mass);
*/
    Teuchos::RCP<std::vector<RealT> > aa_rcp = Teuchos::rcp( new std::vector<RealT> (ni, 1.0) );
    OptDualStdVector<RealT> av(aa_rcp,useRiesz,mass);
    Teuchos::RCP<std::vector<RealT> > bb_rcp = Teuchos::rcp( new std::vector<RealT> (ni, 2.0) );
    OptDualStdVector<RealT> bv(bb_rcp,useRiesz,mass);
    Teuchos::RCP<std::vector<RealT> > cc_rcp = Teuchos::rcp( new std::vector<RealT> (ni, 3.0) );
    OptDualStdVector<RealT> cv(cc_rcp,useRiesz,mass);
    av.checkVector(bv,cv);

    Teuchos::RCP<std::vector<RealT> > dd_rcp = Teuchos::rcp( new std::vector<RealT> (1, 1.0) );
    ConDualStdVector<RealT> dv(dd_rcp,useRiesz,mass);
    Teuchos::RCP<std::vector<RealT> > ee_rcp = Teuchos::rcp( new std::vector<RealT> (1, 2.0) );
    ConDualStdVector<RealT> ev(ee_rcp,useRiesz,mass);
    Teuchos::RCP<std::vector<RealT> > ff_rcp = Teuchos::rcp( new std::vector<RealT> (1, 3.0) );
    ConDualStdVector<RealT> fv(ff_rcp,useRiesz,mass);
    dv.checkVector(ev,fv);

return 0;

/*
    // Gradient   
    Teuchos::RCP<std::vector<RealT> > g_rcp = Teuchos::rcp( new std::vector<RealT> (ni, 0.0) );
    OptDualStdVector<RealT> g(g_rcp,useRiesz,kinetic);

    // Instantiate objective function  
    Objective_GrossPitaevskii<RealT,OptStdVector<RealT>,OptDualStdVector<RealT> > 
        obj(ni,gnl,nb,kinetic,potential,nonlinear);

    // Instantiate normalization constraint
    Normalization_Constraint<RealT,OptStdVector<RealT>,OptDualStdVector<RealT>, 
             ConStdVector<RealT>,ConDualStdVector<RealT> > constr(mass,exactsolve);

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

    return 0;
*/
}
