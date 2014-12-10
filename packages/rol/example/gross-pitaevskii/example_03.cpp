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
  
    bool useRiesz = true;

    // Number of (interior) interpolation points 
    int ni = 5;
    // Number of quadrature points
    int nq = 7;

    const Teuchos::LAPACK<int,RealT> * const lapack = new Teuchos::LAPACK<int,RealT>();
    NodalBasis<RealT> nb(lapack,ni+2,nq);
 
    std::vector<RealT> L(ni*nq,0);
    std::copy(nb.L_.begin()+nq,nb.L_.end()-nq,L.begin()); 
    std::vector<RealT> Lp(ni*nq,0);
    std::copy(nb.Lp_.begin()+nq,nb.Lp_.end()-nq,Lp.begin()); 

    // Quadrature points 
    std::vector<RealT> x(nb.xq_); 

    // Quadrature weights
    std::vector<RealT> w(nb.wq_); 

    // Confinement Potential 
    std::vector<RealT> v(nq,0);   
    
    for(int i=0;i<nq;++i){
        v[i] = 100*x[i]*x[i];
    }

    // Mass matrix
    InnerProductMatrix<RealT> *mass      = new InnerProductMatrixSolver<RealT>(lapack,L,L,w,1);

    // Kinetic energy matrix
    InnerProductMatrix<RealT> *kinetic   = new InnerProductMatrixSolver<RealT>(lapack,Lp,Lp,w,1);

    // Iteration Vector (pointer to optimzation vector)
    Teuchos::RCP<std::vector<RealT> > psi_rcp = Teuchos::rcp( new std::vector<RealT> (ni+2, 0.0) );
    Teuchos::RCP<std::vector<RealT> > psii_rcp = Teuchos::rcp( new std::vector<RealT> (ni, 0.0) );
    
    // Normalized initial guess 
    for(int i=0;i<ni+2;++i){
        (*psi_rcp)[i]=(1+nb.xi_[i])*(1-nb.xi_[i])*sqrt(15.0/16.0);
    } 

    // Interpolate onto the quadrature points
    Teuchos::RCP<std::vector<RealT> > psiq2_rcp = Teuchos::rcp( new std::vector<RealT> (nq, 0.0) );
    nb.lagrange_->interp(*psii_rcp,*psiq2_rcp); 
        
    // Square it
    for(int j=0;j<nq;++j) {
        (*psiq2_rcp)[j] *= (*psiq2_rcp)[j];
    }  
    
    InnerProductMatrix<RealT> *nonlinear = new InnerProductMatrix<RealT>(L,L,w,*psiq2_rcp); 
    InnerProductMatrix<RealT> *potential = new InnerProductMatrix<RealT>(L,L,w,v); 

//    std::cout << S << std::endl;
//    std::cout << nonlinear->inner(psii_rcp,psii_rcp) << std::endl;

    delete mass;
    delete kinetic;
    delete potential; 
    delete nonlinear; 
    return 0;
}
