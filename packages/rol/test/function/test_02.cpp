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

/*! \file  test_02.cpp
    \brief Test derivative checks for log barrier objectives.
*/

#include "ROL_RandomVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_LogBarrierObjective.hpp"
#include "ROL_ObjectiveFromBoundConstraint.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"


typedef double RealT; 

int main(int argc, char *argv[]) {

  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef std::vector<RealT>    vector;
  typedef ROL::Vector<RealT>    V;
  typedef ROL::StdVector<RealT> SV;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  RealT errtol = std::sqrt(ROL::ROL_THRESHOLD);

  int errorFlag  = 0;

  // Specify interval on which to generate uniform random numbers.
  RealT left = 0.01, right = 0.9;

  // *** Test body.

  try {

    int dim = 1;
    RCP<vector>  x_rcp = rcp( new vector(dim,0.0) );
    RCP<vector>  l_rcp = rcp( new vector(dim,0.0) );  // Lower bound
    RCP<vector>  u_rcp = rcp( new vector(dim,1.0) );  // Upper bound
    RCP<vector>  y_rcp = rcp( new vector(dim,0.0) );
    RCP<vector>  v_rcp = rcp( new vector(dim,0.0) );
    RCP<vector>  d_rcp = rcp( new vector(dim,0.0) );
    RCP<vector> gx_rcp = rcp( new vector(dim,0.0) );    
    RCP<vector> gy_rcp = rcp( new vector(dim,0.0) );
    RCP<vector> hv_rcp = rcp( new vector(dim,0.0) );

    SV x( x_rcp);  
    SV y( y_rcp);
    SV v( v_rcp);
    SV d( d_rcp);
    SV gx(gx_rcp);
    SV gy(gy_rcp); 
    SV hv(hv_rcp);

    RandomizeVector(x,left,right);
    RandomizeVector(v,left,right);
    RandomizeVector(d,left,right);

    RCP<V> l = rcp( new SV(l_rcp) );
    RCP<V> u = rcp( new SV(u_rcp) );

    ROL::BoundConstraint<RealT> bc(l,u);

    ROL::ObjectiveFromBoundConstraint<RealT> bc_obj(bc);

    // Fixed difference step size
    RealT delta = 1.e-8; 

    y.set(x);         // y = x
    y.axpy(delta,d);  // y = x+delta*d

    ROL::LogBarrierObjective<RealT> log_obj;
 
    // Do step size sweep
    *outStream << "Test of single logarithmic penalty objective" << std::endl;
    log_obj.checkGradient(x, d, true, *outStream);          *outStream << "\n"; 
    log_obj.checkHessVec(x, v, true, *outStream);           *outStream << "\n";

    *outStream << "Test of bound constraint as logarithmic penalty objective" << std::endl;
    bc_obj.checkGradient(x, d, true, *outStream);          *outStream << "\n"; 
    bc_obj.checkHessVec(x, v, true, *outStream);           *outStream << "\n";

    


    RealT tol = 0;

    // Compute objective at x and y
    RealT fx = log_obj.value(x,tol);
    RealT fy = log_obj.value(y,tol);
    
    // Compute gradient at x and y
    log_obj.gradient(gx,x,tol);
    log_obj.gradient(gy,y,tol);

    // Compute action of Hessian on v at x
    log_obj.hessVec(hv,v,x,tol);

    // FD gradient error 
    RealT graderr = (fy - fx)/delta - gx.dot(d);

    // FD Hessian error
    RCP<V> dg = gx.clone();
    dg->set(gy);
    dg->axpy(-1.0,gx);
    
    RealT hesserr = ( dg->dot(v) )/delta - hv.dot(d);

    

  


    if( std::abs(graderr) > errtol ) {
      ++errorFlag;
    }

    if( std::abs(hesserr) > errtol ) {
      ++errorFlag;
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

