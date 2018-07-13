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

/*! \file  test_03.cpp
 *  \brief Test the elementwise implementation of bound constraints by 
 *         comparison with the StdBoundConstraint
 */


#include "ROL_Bounds.hpp"
#include "ROL_StdBoundConstraint.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>
#include <iomanip>



int main(int argc, char *argv[]) {

  using namespace Teuchos;

  typedef double RealT;

  typedef std::vector<RealT>      vec;
  typedef ROL::Vector<RealT>      V;
  typedef ROL::StdVector<RealT>   SV;

  GlobalMPISession mpiSession(&argc, &argv);
 
  int iprint = argc - 1;
  ROL::Ptr<std::ostream>  outStream;
  ROL::nullstream bhs;
  if( iprint > 0 )
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;

  RealT errtol = ROL::ROL_THRESHOLD<RealT>();

  try {

    int dim = 5;

    // Upper bound vector
    ROL::Ptr<vec> u_ptr = ROL::makePtr<vec>(dim,0.0);

    // Lower bound vector
    ROL::Ptr<vec> l_ptr = ROL::makePtr<vec>(dim,0.0);

    // Gradient vectors
    ROL::Ptr<vec> g1_ptr = ROL::makePtr<vec>(dim,1.0);
    ROL::Ptr<vec> g2_ptr = ROL::makePtr<vec>(dim,-1.0);

    // Test vectors
    ROL::Ptr<vec> y_ptr  = ROL::makePtr<vec>(dim,0.0);
    ROL::Ptr<vec> v_ptr  = ROL::makePtr<vec>(dim,0.0);
    ROL::Ptr<vec> vs_ptr = ROL::makePtr<vec>(dim,0.0);

    (*y_ptr)[0] =  2.0;
    (*y_ptr)[1] = -4.0;
    (*y_ptr)[2] =  1.0;
    (*y_ptr)[3] = -1.0;
    (*y_ptr)[4] =  8.0; 

    for(int i=0;i<dim;++i) {
      (*u_ptr)[i] =   i+1.0;
      (*l_ptr)[i] = -(i+1.0);
    }
 
    ROL::Ptr<V> lp =  ROL::makePtr<SV>(l_ptr);
    ROL::Ptr<V> up  = ROL::makePtr<SV>(u_ptr);
    ROL::Ptr<V> yp  = ROL::makePtr<SV>(y_ptr);
    ROL::Ptr<V> vp  = ROL::makePtr<SV>(v_ptr);
    ROL::Ptr<V> vsp = ROL::makePtr<SV>(vs_ptr);
    ROL::Ptr<V> g1p = ROL::makePtr<SV>(g1_ptr);
    ROL::Ptr<V> g2p = ROL::makePtr<SV>(g2_ptr);


    ROL::Bounds<RealT>     bc(lp,up);
    ROL::StdBoundConstraint<RealT>  sbc(*l_ptr,*u_ptr);

    // Test project
    vp->set(*yp);
    vsp->set(*yp);

    bc.project(*vp); 
    sbc.project(*vsp);  

    vp->axpy(-1.0,*vsp);
 
    *outStream << "Testing BoundConstraint::project() . . . ";

    if(vp->norm() > errtol) {
      *outStream << "failed." << std::endl;
      ++errorFlag; 
    }
    else *outStream << "passed." << std::endl;

    // Test pruneUpperActive (without gradient)
    vp->set(*yp);
    vsp->set(*yp);

    bc.pruneUpperActive(*vp,*yp,0);      
    sbc.pruneUpperActive(*vsp,*yp,0);      
    vp->axpy(-1.0,*vsp);

    *outStream << "Testing BoundConstraint::pruneUpperActive() without gradient . . . ";

    if(vp->norm() > errtol) {
      *outStream << "failed." << std::endl;
      ++errorFlag; 
    }
    else *outStream << "passed." << std::endl;

    // Test pruneLowerActive (without gradient)
    vp->set(*yp);
    vsp->set(*yp);

    bc.pruneLowerActive(*vp,*yp,0);      
    sbc.pruneLowerActive(*vsp,*yp,0);      

    vp->axpy(-1.0,*vsp);

    *outStream << "Testing BoundConstraint::pruneLowerActive() without gradient . . . ";

    if(vp->norm() > errtol) {
      *outStream << "failed." << std::endl;
      ++errorFlag; 
    }
    else *outStream << "passed." << std::endl;


    // Test pruneUpperActive (with positive gradient)
    vp->set(*yp);
    vsp->set(*yp);

    bc.pruneUpperActive(*vp,*yp,*g1p,0);      
    sbc.pruneUpperActive(*vsp,*yp,*g1p,0);      

    vp->axpy(-1.0,*vsp);

    *outStream << "Testing BoundConstraint::pruneUpperActive() with positive gradient . . . ";

    if(vp->norm() > errtol) {
      *outStream << "failed." << std::endl;
      ++errorFlag; 
    }
    else *outStream << "passed." << std::endl;
     
    // Test pruneLowerActive (with positive gradient)
    vp->set(*yp);
    vsp->set(*yp);

    bc.pruneLowerActive(*vp,*yp,*g1p,0);      
    sbc.pruneLowerActive(*vsp,*yp,*g1p,0);      

    vp->axpy(-1.0,*vsp);

    *outStream << "Testing BoundConstraint::pruneLowerActive() with positive gradient . . . ";

    if(vp->norm() > errtol) {
      *outStream << "failed." << std::endl;
      ++errorFlag; 
    }
    else *outStream << "passed." << std::endl;

   
    // Test pruneUpperActive (with negative gradient)
    vp->set(*yp);
    vsp->set(*yp);

    bc.pruneUpperActive(*vp,*yp,*g2p,0);      
    sbc.pruneUpperActive(*vsp,*yp,*g2p,0);      

    vp->axpy(-1.0,*vsp);

    *outStream << "Testing BoundConstraint::pruneLowerActive() with negative gradient . . . ";

    if(vp->norm() > errtol) {
      *outStream << "failed." << std::endl;
      ++errorFlag; 
    }
    else *outStream << "passed." << std::endl;
     
    // Test pruneLowerActive (with negative gradient)
    vp->set(*yp);
    vsp->set(*yp);

    bc.pruneLowerActive(*vp,*yp,*g2p,0);      
    sbc.pruneLowerActive(*vsp,*yp,*g2p,0);      

    vp->axpy(-1.0,*vsp);

    *outStream << "Testing BoundConstraint::pruneLowerActive() with negative gradient . . . ";

    if(vp->norm() > errtol) {
      *outStream << "failed." << std::endl;
      ++errorFlag; 
    }
    else *outStream << "passed." << std::endl;

    vp->set(*yp);
    vsp->set(*yp);

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

