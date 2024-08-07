// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file test_04.cpp
 *  \brief Test barrier functions and factory
 */


#include "ROL_StdVector.hpp"
#include "ROL_BarrierFunctions.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "ROL_ParameterList.hpp"

int main( int argc, char *argv[] ) {

  typedef double RealT;

  typedef std::vector<RealT>         vector;
  typedef ROL::Vector<RealT>         V;
  typedef ROL::StdVector<RealT>      SV;

  typedef typename vector::size_type uint;

    

  RealT xmin = -1.0;
  RealT xmax =  2.0;

  uint N = 30;

  ROL::Ptr<vector> x_ptr = ROL::makePtr<vector>(N,0.0);
  ROL::Ptr<vector> y_ptr = ROL::makePtr<vector>(N,0.0);
  ROL::Ptr<vector> z_ptr = ROL::makePtr<vector>(N,0.0);

  ROL::Ptr<vector> u_ptr = ROL::makePtr<vector>(N,1.0);
  ROL::Ptr<vector> l_ptr = ROL::makePtr<vector>(N,0.0);

  SV x(x_ptr);
  SV y(y_ptr);
  SV z(z_ptr);

  SV u(u_ptr);
  SV l(l_ptr);

  ROL::Ptr<V> x_minus_l = x.clone();
  ROL::Ptr<V> u_minus_x = x.clone();

  

/*
  ROL::ParameterList logList;
  ROL::ParameterList quadList;

  logList.sublist("Barrier Function").set("Type","Logarithmic");
  quadList.sublist("Barrier Function").set("Type","Quadratic");

  ROL::Elementwise::BarrierFunctionFactory<RealT> logFactory( logList );
  ROL::Elementwise::BarrierFunctionFactory<RealT> quadFactory( quadList );

  ROL::Ptr<ROL::Elementwise::BarrierFunction<RealT> > logFunction = logFactory.getBarrierFunction();
  ROL::Ptr<ROL::Elementwise::BarrierFunction<RealT> > quadFunction = quadFactory.getBarrierFunction();
*/
  
  ROL::Ptr<ROL::Elementwise::BinaryFunction<RealT> > lesser = ROL::makePtr<ROL::Elementwise::Lesser<RealT>>();
  ROL::Ptr<ROL::Elementwise::BinaryFunction<RealT> > greater = ROL::makePtr<ROL::Elementwise::Greater<RealT>>();
   

  for( uint i=0; i<N; ++i ) {
    RealT t = static_cast<RealT>(i)/(N-1);
    (*x_ptr)[i] = xmin*(1-t) + xmax*t;
  }


  typedef ROL::Elementwise::Axpy<RealT>              Axpy;
  typedef ROL::Elementwise::Scale<RealT>             Scale;
  typedef ROL::Elementwise::Logarithm<RealT>         Logarithm;
  
  



  x_minus_l->set(l);
  x_minus_l->applyBinary(Axpy(-1.0),x);  
  x_minus_l->applyUnary(Logarithm());
  x_minus_l->applyUnary(Scale(-1.0));

  u_minus_x->set(x);
  u_minus_x->applyBinary(Axpy(-1.0),u);
  u_minus_x->applyUnary(Logarithm());
  u_minus_x->applyUnary(Scale(-1.0));


  y.set(x);
  y.applyBinary(*lesser,u);
  y.applyBinary(*greater,l);
  
  ROL::Ptr<vector> xml = ROL::staticPtrCast<SV>(x_minus_l)->getVector();
  ROL::Ptr<vector> umx = ROL::staticPtrCast<SV>(u_minus_x)->getVector();

  for(uint i=0; i<N; ++i ) {
    std::cout << std::setw(16) << (*xml)[i] 
              << std::setw(16) << (*umx)[i] << std::endl;
  } 

 


  return 0;
}


