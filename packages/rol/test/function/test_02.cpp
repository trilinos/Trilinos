// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_02.cpp
    \brief Test derivative checks for log barrier objectives.
*/

#include "ROL_RandomVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_ObjectiveFromBoundConstraint.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "ROL_ParameterList.hpp"


typedef double RealT; 

int main(int argc, char *argv[]) {

  typedef std::vector<RealT>    vector;
  typedef ROL::Vector<RealT>    V;
  typedef ROL::StdVector<RealT> SV;

  typedef typename vector::size_type uint;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  // Save the format state of the original std::cout.
  ROL::nullstream oldFormatState;
  oldFormatState.copyfmt(std::cout);

//  RealT errtol = std::sqrt(ROL::ROL_THRESHOLD<RealT>());

  int errorFlag  = 0;

  // *** Test body.

  try {

    uint dim = 30;
    RealT xmin = 0.5;
    RealT xmax = 2.5;

    ROL::Ptr<vector>  x_ptr  = ROL::makePtr<vector>(dim,0.0);
    ROL::Ptr<vector>  g_ptr  = ROL::makePtr<vector>(dim,0.0);
    ROL::Ptr<vector>  v_ptr  = ROL::makePtr<vector>(dim,0.0);
    ROL::Ptr<vector>  hv_ptr = ROL::makePtr<vector>(dim,0.0);

    ROL::Ptr<vector>  l_ptr = ROL::makePtr<vector>(dim,1.0);
    ROL::Ptr<vector>  u_ptr = ROL::makePtr<vector>(dim,2.0);

    ROL::Ptr<vector>  xlog0_ptr = ROL::makePtr<vector>(dim,0.0);
    ROL::Ptr<vector>  xlog1_ptr = ROL::makePtr<vector>(dim,0.0);
    ROL::Ptr<vector>  xlog2_ptr = ROL::makePtr<vector>(dim,0.0);

    ROL::Ptr<vector>  xquad0_ptr = ROL::makePtr<vector>(dim,0.0);
    ROL::Ptr<vector>  xquad1_ptr = ROL::makePtr<vector>(dim,0.0);
    ROL::Ptr<vector>  xquad2_ptr = ROL::makePtr<vector>(dim,0.0);

    ROL::Ptr<vector>  xdwell0_ptr = ROL::makePtr<vector>(dim,0.0);
    ROL::Ptr<vector>  xdwell1_ptr = ROL::makePtr<vector>(dim,0.0);
    ROL::Ptr<vector>  xdwell2_ptr = ROL::makePtr<vector>(dim,0.0);



    SV x(x_ptr); 
    SV g(g_ptr);
    SV v(v_ptr);
    SV hv(hv_ptr);

    ROL::Ptr<SV> xlog0 = ROL::makePtr<SV>(xlog0_ptr);
    ROL::Ptr<SV> xlog1 = ROL::makePtr<SV>(xlog1_ptr);
    ROL::Ptr<SV> xlog2 = ROL::makePtr<SV>(xlog2_ptr);

    ROL::Ptr<SV> xquad0 = ROL::makePtr<SV>(xquad0_ptr);
    ROL::Ptr<SV> xquad1 = ROL::makePtr<SV>(xquad1_ptr);
    ROL::Ptr<SV> xquad2 = ROL::makePtr<SV>(xquad2_ptr);

    ROL::Ptr<SV> xdwell0 = ROL::makePtr<SV>(xdwell0_ptr);
    ROL::Ptr<SV> xdwell1 = ROL::makePtr<SV>(xdwell1_ptr);
    ROL::Ptr<SV> xdwell2 = ROL::makePtr<SV>(xdwell2_ptr);

    ROL::Ptr<V> lo = ROL::makePtr<SV>(l_ptr);
    ROL::Ptr<V> up = ROL::makePtr<SV>(u_ptr);  

    for(uint i=0; i<dim; ++i) {
      RealT t = static_cast<RealT>(i)/static_cast<RealT>(dim-1);
      (*x_ptr)[i] = xmin*(1-t) + xmax*t;
    }    

    // Create bound constraint
    ROL::Bounds<RealT>  bc(lo,up);

    ROL::ParameterList logList;
    ROL::ParameterList quadList;
    ROL::ParameterList dwellList;

    logList.sublist("Barrier Function").set("Type","Logarithmic");
    quadList.sublist("Barrier Function").set("Type","Quadratic");
    dwellList.sublist("Barrier Function").set("Type","Double Well");

    ROL::ObjectiveFromBoundConstraint<RealT> logObj(bc,logList);
    ROL::ObjectiveFromBoundConstraint<RealT> quadObj(bc,quadList);
    ROL::ObjectiveFromBoundConstraint<RealT> dwellObj(bc,dwellList);

    RealT tol = 0.0;


    logObj.value(x,tol);
    xlog0->set(*ROL::staticPtrCast<SV>(logObj.getBarrierVector()));

    logObj.gradient(g,x,tol);
    xlog1->set(*ROL::staticPtrCast<SV>(logObj.getBarrierVector()));

    logObj.hessVec(hv,v,x,tol);
    xlog2->set(*ROL::staticPtrCast<SV>(logObj.getBarrierVector()));


    quadObj.value(x,tol);
    xquad0->set(*ROL::staticPtrCast<SV>(quadObj.getBarrierVector()));

    quadObj.gradient(g,x,tol);
    xquad1->set(*ROL::staticPtrCast<SV>(quadObj.getBarrierVector()));

    quadObj.hessVec(hv,v,x,tol);
    xquad2->set(*ROL::staticPtrCast<SV>(quadObj.getBarrierVector()));


    dwellObj.value(x,tol);
    xdwell0->set(*ROL::staticPtrCast<SV>(dwellObj.getBarrierVector()));

    dwellObj.gradient(g,x,tol);
    xdwell1->set(*ROL::staticPtrCast<SV>(dwellObj.getBarrierVector()));

    dwellObj.hessVec(hv,v,x,tol);
    xdwell2->set(*ROL::staticPtrCast<SV>(dwellObj.getBarrierVector()));


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
      *outStream << std::setw(14) << (*x_ptr)[i] 
                 << std::setw(14) << (*xlog0_ptr)[i] 
                 << std::setw(14) << (*xlog1_ptr)[i]
                 << std::setw(14) << (*xlog2_ptr)[i] 
                 << std::setw(14) << (*xquad0_ptr)[i] 
                 << std::setw(14) << (*xquad1_ptr)[i] 
                 << std::setw(14) << (*xquad2_ptr)[i] 
                 << std::setw(14) << (*xdwell0_ptr)[i] 
                 << std::setw(14) << (*xdwell1_ptr)[i] 
                 << std::setw(14) << (*xdwell2_ptr)[i] 
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
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;


}

