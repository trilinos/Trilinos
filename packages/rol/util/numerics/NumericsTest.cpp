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

//! \brief Verify the numerics are working correctly 

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "OrthogonalPolynomials.hpp"
#include "LinearAlgebra.hpp"

#include <iostream>
#include <iomanip>


typedef double RealT;

int main(int argc, char *argv[]) {

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

    int N = atoi(argv[1]);
    const RealT pi  = 3.14159265359; 
    const RealT tol = 1e-7;

    Teuchos::LAPACK<int,RealT> lapack;

    Teuchos::RCP<std::vector<RealT> > ap  = Teuchos::rcp( new std::vector<RealT>(N,0) );
    Teuchos::RCP<std::vector<RealT> > bp  = Teuchos::rcp( new std::vector<RealT>(N,0) );
    Teuchos::RCP<std::vector<RealT> > xgp = Teuchos::rcp( new std::vector<RealT>(N,0) );
    Teuchos::RCP<std::vector<RealT> > wgp = Teuchos::rcp( new std::vector<RealT>(N,0) );
    Teuchos::RCP<std::vector<RealT> > xlp = Teuchos::rcp( new std::vector<RealT>(N,0) );
    Teuchos::RCP<std::vector<RealT> > wlp = Teuchos::rcp( new std::vector<RealT>(N,0) );

    ROL::StdVector<RealT> a(ap);
    ROL::StdVector<RealT> b(bp);
    ROL::StdVector<RealT> xg(xgp);
    ROL::StdVector<RealT> wg(wgp);
    ROL::StdVector<RealT> xl(xlp);
    ROL::StdVector<RealT> wl(wlp);

    // Compute Legendre Recursion coefficients 
    rec_jacobi(0,0,a,b);

    // Compute Gauss nodes and weights
    gauss(lapack,a,b,xg,wg);
 
    // Modify the recursion coefficients for Lobatto
    rec_lobatto(lapack,-1.0,1.0,a,b); 

    // Compute Gauss nodes and weights
    gauss(lapack,a,b,xl,wl);
 
    // Numerically integrate k*sech(k*x) for somewhat large k
    RealT fg = 0; 
    RealT fl = 0;
    RealT k = 20;
 
    for(int i=0;i<N;++i) {
        fg += k*(*wgp)[i]/cosh(k*(*xgp)[i]);
        fl += k*(*wlp)[i]/cosh(k*(*xlp)[i]);
    }
   
    if(fabs(fg-pi)>tol) {
        ++errorFlag;
    } 
     if(fabs(fl-pi)>tol) {
        ++errorFlag;
    } 
 
    if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
    else
        std::cout << "End Result: TEST PASSED\n";




}
