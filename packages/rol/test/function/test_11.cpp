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

/*! \file  test_11.cpp
    \brief Verify that the implementation of the Coleman-Li Trust-Region
           model passes derivative checks 
*/

#include "ROL_ColemanLiModel.hpp"
#include "ROL_HS2.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_RandomVector.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
//#include <fenv.h>

typedef double RealT;

int main(int argc, char *argv[]) {
//  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  typedef ROL::Vector<RealT>          V;
  typedef ROL::Objective<RealT>       OBJ;
  typedef ROL::BoundConstraint<RealT> CON; 
  

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  RealT zero(0);

  ROL::Ptr<V>   x0;
  ROL::Ptr<V>   x;
  ROL::Ptr<V>   g;
  ROL::Ptr<OBJ> obj;
  ROL::Ptr<CON> con;
  ROL::Ptr<OBJ> model;  

  ROL::ZOO::getHS2<RealT> HS2;
  obj = HS2.getObjective();
  con = HS2.getBoundConstraint();
  x0  = HS2.getInitialGuess();
  x   = HS2.getSolution();

  g = x->dual().clone();

  // Need to evaluate the gradient to construct the model
  obj->gradient(*g,*x,zero);

  model = ROL::makePtr<ROL::ColemanLiModel<RealT>>(*obj,*con,*x,*g);

  ROL::Ptr<V> s = x->clone();
  ROL::Ptr<V> v = x->clone();
  ROL::Ptr<V> u = x->clone();

  ROL::RandomizeVector(*s,-1.0,1.0);
  ROL::RandomizeVector(*u,-1.0,1.0);
  ROL::RandomizeVector(*v,-1.0,1.0);

  model->checkGradient(*s,*v);
  model->checkHessVec(*s,*v);
  model->checkHessSym(*s,*u,*v);

  return 0;
}


