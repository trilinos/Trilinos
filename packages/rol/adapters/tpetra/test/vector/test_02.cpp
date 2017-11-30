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
 *  \brief Test SimulatedVector functionality with TpetraTeuchosBatchManager
 */

#include "ROL_SimulatedVector.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Types.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Comm.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

typedef double RealT;

template<class Real>
void print_vector( const ROL::Vector<Real> &x ) {

  typedef ROL::Vector<Real>            V;
  typedef ROL::StdVector<Real>         SV;
  typedef ROL::SimulatedVector<Real>   PV;
  typedef typename PV::size_type       size_type;

  const PV eb = dynamic_cast<const PV&>(x);
  size_type n = eb.numVectors();
    
  for(size_type k=0; k<n; ++k) {
    std::cout << "[subvector " << k << "]" << std::endl;
    ROL::Ptr<const V> vec = eb.get(k);
    ROL::Ptr<const std::vector<Real> > vp = 
      dynamic_cast<const SV&>(*vec).getVector();  
   for(size_type i=0;i<vp->size();++i) {
      std::cout << (*vp)[i] << std::endl;
    }  
  }
}


int main(int argc, char *argv[]) {

  using namespace Teuchos;

  typedef ROL::Vector<RealT>            V;
  typedef ROL::StdVector<RealT>         SV;
  typedef ROL::SimulatedVector<RealT>   PV;

  GlobalMPISession mpiSession(&argc, &argv);
  ROL::Ptr<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  ROL::Ptr<ROL::TpetraTeuchosBatchManager<RealT> > bman = ROL::makePtr<ROL::TpetraTeuchosBatchManager<RealT>>(comm);

  int iprint = argc - 1;

  ROL::Ptr<std::ostream> outStream;
  oblackholestream bhs; // no output
 
  if( iprint>0 ) 
   outStream = ROL::makePtrFromRef(std::cout);
  else
   outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;

  RealT errtol = ROL::ROL_THRESHOLD<RealT>();

  try {

    int stoch_dim = 4;
    int nSamp = 100;
    std::vector<RealT> tmp(2,0.0); tmp[0] = -1.0; tmp[1] = 1.0;
    std::vector<std::vector<RealT> > bounds(stoch_dim,tmp);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nSamp,bounds,bman);

    int batchID = bman->batchID();
    int nvecloc = sampler->numMySamples();

    *outStream << "Proc " << batchID << ", number of vectors = " << nvecloc << std::endl;

    int dim = 5;

    int total_dim = 0;
     
    RealT left = -1e0, right = 1e0;

    std::vector<ROL::Ptr<V> > x_ptr;
    std::vector<ROL::Ptr<V> > y_ptr;
    std::vector<ROL::Ptr<V> > z_ptr;

    for( int k=0; k<nvecloc; ++k ) {
       
      ROL::Ptr<std::vector<RealT> > xk_ptr = ROL::makePtr<std::vector<RealT>>(dim);
      ROL::Ptr<std::vector<RealT> > yk_ptr = ROL::makePtr<std::vector<RealT>>(dim);
      ROL::Ptr<std::vector<RealT> > zk_ptr = ROL::makePtr<std::vector<RealT>>(dim);

      ROL::Ptr<V> xk = ROL::makePtr<SV>( xk_ptr );
      ROL::Ptr<V> yk = ROL::makePtr<SV>( yk_ptr );
      ROL::Ptr<V> zk = ROL::makePtr<SV>( zk_ptr );

      for( int i=0; i<dim; ++i ) {
        (*xk_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
        (*yk_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
        (*zk_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      }
   
      x_ptr.push_back(xk);
      y_ptr.push_back(yk);
      z_ptr.push_back(zk);
      
      total_dim += dim;
    }

    *outStream << "Proc " << batchID << ", total dimension = " << total_dim << std::endl;

    PV x(x_ptr, bman);
    PV y(y_ptr, bman);
    PV z(z_ptr, bman);

    // Standard tests.
    std::vector<RealT> consistency = x.checkVector(y, z, true, *outStream);
    ROL::StdVector<RealT> checkvec(ROL::makePtrFromRef(consistency));
    if (checkvec.norm() > std::sqrt(errtol)) {
      errorFlag++;
    }

    // Repeat the checkVector tests with a zero vector.
    x.scale(0.0);
    consistency = x.checkVector(x, x, true, *outStream);
    if (checkvec.norm() > 0.0) {
      errorFlag++;
    }

    // Reset x and y vectors to all 1 and 2, resp., and compute dot product.
/*
    x_ptr.resize(0);
    y_ptr.resize(0);
    for( int k=0; k<nvecloc; ++k ) {
      ROL::Ptr<std::vector<RealT> > xk_ptr = ROL::makePtr<std::vector<RealT>>(dim);
      ROL::Ptr<std::vector<RealT> > yk_ptr = ROL::makePtr<std::vector<RealT>>(dim);
      ROL::Ptr<V> xk = ROL::makePtr<SV>( xk_ptr );
      ROL::Ptr<V> yk = ROL::makePtr<SV>( yk_ptr );
      for( int i=0; i<dim; ++i ) {
        (*xk_ptr)[i] = 1.0;
        (*yk_ptr)[i] = 2.0;
      }
      x_ptr.push_back(xk);
      y_ptr.push_back(yk);
    }
    *outStream << "x.dot(y) = " << x.dot(y) << std::endl;
    if (std::abs(x.dot(y) - nSamp*dim*2) > std::sqrt(errtol)) {
      errorFlag++;
    }
*/

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
