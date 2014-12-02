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

/*! \file   example_01.cpp
    \brief  Minimize the Gross-Pitaevskii functional and demonstrate 
            the effect of choice of function space of the Gradient on
            convergence.

    \details Minimize the one-dimensional Gross-Pitaevskii (GP) energy 
             functional
             \f[ J[\psi] = \int \frac{1}{2} |\nabla\psi|^2 + V(x)|\psi|^2 
                           +\frac{g}{2}|\psi|^4 \,\mathrm{d}x \f]
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

#include <iostream>

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_EqualityConstraint.hpp"
#include "ROL_CompositeStepSQP.hpp"
#include "ROL_Algorithm.hpp"


using namespace ROL;



/* Define Objective Function class */  
template<class Real>
class Objective_GrossPitaevskii : public Objective<Real> {
 
    // STL Vector
    typedef std::vector<Real> svec;

    // ROL StdVector
    typedef StdVector<Real> rvec;

    // Pointer to const STL vector
    typedef Teuchos::RCP<const svec> pcsv;

    // Pointer to STL vector
    typedef Teuchos::RCP<svec> psv;
  

    private:
        Real g_;     // Coefficient \f$g\f$ appearing before quartic term in GP functional     
        int  nx_;    // Number of interior nodes
        Real dx_;    // Mesh spacing \f$ \Delta x = \frac{1}{n_x+1} \f$ 
        
        pcsv Vp_;    // Pointer to potential vector

        void applyK(const Vector<Real> &v, Vector<Real> &kv) {
            /** Pointer to direction vector */
            pcsv vp = (Teuchos::dyn_cast<rvec>(const_cast<Vector<Real> &>(v))).getVector();

            /** Pointer to action of Hessian on direction vector */
            psv kvp = Teuchos::rcp_const_cast<svec>((Teuchos::dyn_cast<rvec>(kv)).getVector());

            Real dx2 = dx_*dx_;

            (*kvp)[0] = (2.0*(*vp)[0]-(*vp)[1])/dx2;
  
            for(int i=1;i<nx_-1;++i) {
                (*kvp)[i] = (2.0*(*vp)[i]-(*vp)[i-1]-(*vp)[i+1])/dx2;
            } 

            (*kvp)[nx_-1] = (2.0*(*vp)[nx_-1]-(*vp)[nx_-2])/dx2;

        } 

    public: 

        Objective_GrossPitaevskii(const Real &g, const Vector<Real> &V) : g_(g),  
            Vp_((Teuchos::dyn_cast<rvec>(const_cast<Vector<Real> &>(V))).getVector())  {

            nx_ = Vp_->size(); 
            dx_ = (1.0/(1.0+nx_));
        }
           

    /** Evaluate cost functional */
    Real value( const Vector<Real> &psi, Real &tol ) {
        
        /** Pointer to opt vector */
        pcsv psip = (Teuchos::dyn_cast<rvec>(const_cast<Vector<Real> &>(psi))).getVector();

        /** Pointer to K applied to opt vector */
        psv kpsip = Teuchos::rcp( new svec(nx_, 0.0) );
        rvec kpsi(kpsip);

        Real J = 0;

        this->applyK(psi,kpsi);

        for(int i=0;i<nx_;++i) {
            J += (*psip)[i]*(*kpsip)[i] + (*Vp_)[i]*pow((*psip)[i],2) + g_*pow((*psip)[i],4);
        } 
      
        // Rescale for numerical integration by trapezoidal rule
        J *= 0.5*dx_;

        return J;
    }


    /** Evaluate gradient */
    void gradient( Vector<Real> &g, const Vector<Real> &psi, Real &tol ) {

        /** Pointer to opt vector */
        pcsv psip = (Teuchos::dyn_cast<rvec>(const_cast<Vector<Real> &>(psi))).getVector();

        /** Pointer to gradient vector */
        psv gp = Teuchos::rcp_const_cast<svec>((Teuchos::dyn_cast<rvec>(g)).getVector());

        /** Pointer to K applied to opt vector */
        psv kpsip = Teuchos::rcp( new svec(nx_, 0.0) );
        rvec kpsi(kpsip);

        this->applyK(psi,kpsi);

        for(int i=0;i<nx_;++i) {
            (*gp)[i] = ((*kpsip)[i] + (*Vp_)[i]*(*psip)[i] + 2.0*g_*pow((*psip)[i],3))*dx_;
        } 
      
    }



    /** Compute the action of the Hessian evaluated at psi on a vector v */
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &psi, Real &tol ) {

        /** Pointer to opt vector */
        pcsv psip = (Teuchos::dyn_cast<rvec>(const_cast<Vector<Real> &>(psi))).getVector();

        /** Pointer to direction vector */
        pcsv vp = (Teuchos::dyn_cast<rvec>(const_cast<Vector<Real> &>(v))).getVector();

        /** Pointer to action of Hessian on direction vector */
        psv hvp = Teuchos::rcp_const_cast<svec>((Teuchos::dyn_cast<rvec>(hv)).getVector());

        this->applyK(v,hv);
 
        for(int i=0;i<nx_;++i) {
            (*hvp)[i] *= dx_;
            (*hvp)[i] += ( (*Vp_)[i] + 6.0*g_*pow((*psip)[i],2) )*(*vp)[i]*dx_;
        } 

   }

};


/** Define class to enforce normalization condition*/  
template<class Real>
class Normalization_Constraint : public EqualityConstraint<Real> {

    // STL Vector
    typedef std::vector<Real> svec;

    // ROL StdVector
    typedef StdVector<Real> rvec;

    // Pointer to const STL vector
    typedef Teuchos::RCP<const svec> pcsv;

    // Pointer to STL vector
    typedef Teuchos::RCP<svec> psv;
 
    private:     
    int nx_;
    Real dx_;

    public:
    Normalization_Constraint(int n, Real dx) : nx_(n), dx_(dx) {}          


    void value(Vector<Real> &c, const Vector<Real> &psi, Real &tol){

        /** Pointer to constraint vector (only one element)*/
        psv cp = Teuchos::rcp_const_cast<svec>((Teuchos::dyn_cast<rvec>(c)).getVector());

        /** Pointer to optimization vector */    
        pcsv psip = (Teuchos::dyn_cast<rvec>(const_cast<Vector<Real> &>(psi))).getVector();

        (*cp)[0] = -1.0;
        for(int i=0;i<nx_;++i) {
            (*cp)[0] += dx_*pow((*psip)[i],2);
        } 
    }


    void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &psi, Real &tol){

        /** Pointer to action of Jacobian of constraint on direction vector (yields scalar)*/
        psv jvp = Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(jv)).getVector());

        /** Pointer to direction vector */    
        pcsv vp = (Teuchos::dyn_cast<rvec>(const_cast<Vector<Real> &>(v))).getVector();

        /** Pointer to optimization vector */    
        pcsv psip = (Teuchos::dyn_cast<rvec>(const_cast<Vector<Real> &>(psi))).getVector();
      
        (*jvp)[0] = 0;
        for(int i=0;i<nx_;++i) {
            (*jvp)[0] += 2.0*dx_*(*psip)[i]*(*vp)[i];
        }
    }


    void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &psi, Real &tol){

        /** Pointer to action of adjoint of Jacobian of constraint on direction vector (yields vector)*/
        psv ajvp = Teuchos::rcp_const_cast<svec>((Teuchos::dyn_cast<rvec>(ajv)).getVector());

        /** Pointer to direction vector */    
        pcsv vp = (Teuchos::dyn_cast<rvec>(const_cast<Vector<Real> &>(v))).getVector();

        /** Pointer to optimization vector */    
        pcsv psip = (Teuchos::dyn_cast<rvec>(const_cast<Vector<Real> &>(psi))).getVector();

        for(int i=0;i<nx_;++i) {
            (*ajvp)[i] = 2.0*dx_*(*psip)[i]*(*vp)[0];
        }
    }

    void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, 
                             const Vector<Real> &psi, Real &tol){

        /** The pointer to action of constraint Hessian in u,v inner product */
        psv ahuvp = Teuchos::rcp_const_cast<svec>((Teuchos::dyn_cast<rvec>(ahuv)).getVector());

        /** Pointer to direction vector u */    
        pcsv up = (Teuchos::dyn_cast<rvec>(const_cast<Vector<Real> &>(u))).getVector();

        /** Pointer to direction vector v */    
        pcsv vp = (Teuchos::dyn_cast<rvec>(const_cast<Vector<Real> &>(v))).getVector();

        /** Pointer to optimization vector */    
        pcsv psip = (Teuchos::dyn_cast<rvec>(const_cast<Vector<Real> &>(psi))).getVector();

        
        for(int i=0;i<nx_;++i) {
            (*ahuvp)[i] = 2.0*dx_*(*vp)[i]*(*up)[0];        
        }  
    }
};


typedef double RealT;

int main(int argc, char **argv) {

    // Output file
    std::ofstream outfile ("output.csv");

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

    // Number of interior grid points \f$n_x\f$
    int nx = 10;
    if(argc > 1) {
        nx = atoi( argv[1] );
    }
 
    // Coefficient of quartic term
    RealT gnl = 10.0;
    if(argc > 2) {
        gnl = atof( argv[2] );
    }
    
    // Grid spacing
    RealT dx = 1.0/(nx+1);

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

       
    /** Set Initial Guess (normalized) */
    double sqrt30 = sqrt(30);

    for (int i=0; i<nx; i++) {
        (*psi_rcp)[i]   = sqrt30*(*xi_rcp)[i]*(1.0-(*xi_rcp)[i]);
    }

    StdVector<RealT> psi(psi_rcp);

    /** Equality constraint value (scalar) */  
    Teuchos::RCP<std::vector<RealT> > c_rcp = Teuchos::rcp( new std::vector<RealT> (1, 0.0) );
    StdVector<RealT> c(c_rcp);

    /** Lagrange multiplier value (scalar) */  
    Teuchos::RCP<std::vector<RealT> > lam_rcp = Teuchos::rcp( new std::vector<RealT> (1, 0.0) );
    StdVector<RealT> lam(lam_rcp);

    /** Gradient */  
    Teuchos::RCP<std::vector<RealT> > g_rcp = Teuchos::rcp( new std::vector<RealT> (nx, 0.0) );
    StdVector<RealT> g(g_rcp);

    /** Instantiate objective function */ 
    Objective_GrossPitaevskii<RealT> obj(gnl,V);
    
    /** Instantiate normalization constraint */
    Normalization_Constraint<RealT> constr(nx,dx);

    Teuchos::ParameterList parlist;

    // Define Step
    parlist.set("Nominal SQP Optimality Solver Tolerance", 1.e-4);
    parlist.set("Maximum Number of Krylov Iterations",80);
    parlist.set("Absolute Krylov Tolerance",1e-4);

    ROL::CompositeStepSQP<RealT> step(parlist);


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

    // Write x,psi(x) to file
    outfile << 0 << "," << 0 << std::endl;

    for(int i=0;i<nx;++i) {
        outfile << (*xi_rcp)[i] << "," << (*psi_rcp)[i] << std::endl ;
    }

    outfile << 1 << "," << 0 << std::endl;

    return 0;

}
