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
#include "ROL_ObjectiveFromBoundConstraint.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"


typedef double RealT; 

int main(int argc, char *argv[]) {

  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef std::vector<RealT>    vector;
  typedef ROL::Vector<RealT>    V;
  typedef ROL::StdVector<RealT> SV;

  typedef typename vector::size_type uint;

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

//  RealT errtol = std::sqrt(ROL::ROL_THRESHOLD<RealT>());

  int errorFlag  = 0;

  // *** Test body.

  try {

    uint dim = 30;
    RealT xmin = 0.5;
    RealT xmax = 2.5;

    RCP<vector>  x_rcp  = rcp( new vector(dim,0.0) );
    RCP<vector>  g_rcp  = rcp( new vector(dim,0.0) );
    RCP<vector>  v_rcp  = rcp( new vector(dim,0.0) );
    RCP<vector>  hv_rcp = rcp( new vector(dim,0.0) );

    RCP<vector>  l_rcp = rcp( new vector(dim,1.0) );
    RCP<vector>  u_rcp = rcp( new vector(dim,2.0) );

    RCP<vector>  xlog0_rcp = rcp( new vector(dim,0.0) );
    RCP<vector>  xlog1_rcp = rcp( new vector(dim,0.0) );
    RCP<vector>  xlog2_rcp = rcp( new vector(dim,0.0) );

    RCP<vector>  xquad0_rcp = rcp( new vector(dim,0.0) );
    RCP<vector>  xquad1_rcp = rcp( new vector(dim,0.0) );
    RCP<vector>  xquad2_rcp = rcp( new vector(dim,0.0) );

    RCP<vector>  xdwell0_rcp = rcp( new vector(dim,0.0) );
    RCP<vector>  xdwell1_rcp = rcp( new vector(dim,0.0) );
    RCP<vector>  xdwell2_rcp = rcp( new vector(dim,0.0) );



    SV x(x_rcp); 
    SV g(g_rcp);
    SV v(v_rcp);
    SV hv(hv_rcp);

    RCP<SV> xlog0 = rcp( new SV(xlog0_rcp) );
    RCP<SV> xlog1 = rcp( new SV(xlog1_rcp) );
    RCP<SV> xlog2 = rcp( new SV(xlog2_rcp) );

    RCP<SV> xquad0 = rcp( new SV(xquad0_rcp) );
    RCP<SV> xquad1 = rcp( new SV(xquad1_rcp) );
    RCP<SV> xquad2 = rcp( new SV(xquad2_rcp) );

    RCP<SV> xdwell0 = rcp( new SV(xdwell0_rcp) );
    RCP<SV> xdwell1 = rcp( new SV(xdwell1_rcp) );
    RCP<SV> xdwell2 = rcp( new SV(xdwell2_rcp) );

    RCP<V> lo = rcp(new SV(l_rcp) );
    RCP<V> up = rcp(new SV(u_rcp) );  

    for(uint i=0; i<dim; ++i) {
      RealT t = static_cast<RealT>(i)/static_cast<RealT>(dim-1);
      (*x_rcp)[i] = xmin*(1-t) + xmax*t;
    }    

    // Create bound constraint
    ROL::BoundConstraint<RealT>  bc(lo,up);

    Teuchos::ParameterList logList;
    Teuchos::ParameterList quadList;
    Teuchos::ParameterList dwellList;

    logList.sublist("Barrier Function").set("Type","Logarithmic");
    quadList.sublist("Barrier Function").set("Type","Quadratic");
    dwellList.sublist("Barrier Function").set("Type","Double Well");

    ROL::ObjectiveFromBoundConstraint<RealT> logObj(bc,logList);
    ROL::ObjectiveFromBoundConstraint<RealT> quadObj(bc,quadList);
    ROL::ObjectiveFromBoundConstraint<RealT> dwellObj(bc,dwellList);

    RealT tol = 0.0;


    logObj.value(x,tol);
    xlog0->set(*Teuchos::rcp_static_cast<SV>(logObj.getBarrierVector()));

    logObj.gradient(g,x,tol);
    xlog1->set(*Teuchos::rcp_static_cast<SV>(logObj.getBarrierVector()));

    logObj.hessVec(hv,v,x,tol);
    xlog2->set(*Teuchos::rcp_static_cast<SV>(logObj.getBarrierVector()));


    quadObj.value(x,tol);
    xquad0->set(*Teuchos::rcp_static_cast<SV>(quadObj.getBarrierVector()));

    quadObj.gradient(g,x,tol);
    xquad1->set(*Teuchos::rcp_static_cast<SV>(quadObj.getBarrierVector()));

    quadObj.hessVec(hv,v,x,tol);
    xquad2->set(*Teuchos::rcp_static_cast<SV>(quadObj.getBarrierVector()));


    dwellObj.value(x,tol);
    xdwell0->set(*Teuchos::rcp_static_cast<SV>(dwellObj.getBarrierVector()));

    dwellObj.gradient(g,x,tol);
    xdwell1->set(*Teuchos::rcp_static_cast<SV>(dwellObj.getBarrierVector()));

    dwellObj.hessVec(hv,v,x,tol);
    xdwell2->set(*Teuchos::rcp_static_cast<SV>(dwellObj.getBarrierVector()));


    *outStream   << std::setw(14) << "x" 
                 << std::setw(14) << "log" 
                 << std::setw(14) << "D(log)" 
                 << std::setw(14) << "D2(log)" 
                 << std::setw(14) << "quad" 
                 << std::setw(14) << "D(quad)" 
                 << std::setw(14) << "D2(quad)" 
                 << std::setw(14) << "dwell" 
                 << std::setw(14) << "D(dwell)" 
                 << std::setw(14) << "D2(dwell)" 
                 << std::endl;
    *outStream   << std::string(140,'-') << std::endl;

    for(uint i=0; i<dim; ++i) {
      *outStream << std::setw(14) << (*x_rcp)[i] 
                 << std::setw(14) << (*xlog0_rcp)[i] 
                 << std::setw(14) << (*xlog1_rcp)[i]
                 << std::setw(14) << (*xlog2_rcp)[i] 
                 << std::setw(14) << (*xquad0_rcp)[i] 
                 << std::setw(14) << (*xquad1_rcp)[i] 
                 << std::setw(14) << (*xquad2_rcp)[i] 
                 << std::setw(14) << (*xdwell0_rcp)[i] 
                 << std::setw(14) << (*xdwell1_rcp)[i] 
                 << std::setw(14) << (*xdwell2_rcp)[i] 
                 << std::endl;
    }    

  
    ROL::RandomizeVector( x,  1.2, 1.8 );
    ROL::RandomizeVector( v, -0.1, 0.1 );

    *outStream << "\n\n";
    *outStream << "Test of logarithmic penalty objective" << std::endl;
    logObj.checkGradient(x,v,true,*outStream);    *outStream << std::endl;
    logObj.checkHessVec(x,v,true,*outStream);     *outStream << std::endl;

    ROL::RandomizeVector( x, -1.0, 1.0 );
    ROL::RandomizeVector( v, -1.0, 1.0 );

    *outStream << "\n\n";
    *outStream << "Test of piecewise quadratic penalty objective" << std::endl;
    quadObj.checkGradient(x,v,true,*outStream);    *outStream << std::endl;
    quadObj.checkHessVec(x,v,true,*outStream);     *outStream << std::endl;


    *outStream << "\n\n";
    *outStream << "Test of double well penalty objective" << std::endl;
    dwellObj.checkGradient(x,v,true,*outStream);    *outStream << std::endl;
    dwellObj.checkHessVec(x,v,true,*outStream);     *outStream << std::endl;
    




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

