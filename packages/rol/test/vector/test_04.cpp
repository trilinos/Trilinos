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

/*! \file  test_04.cpp
 *  \brief Test PartitionedVector functionality
 */

#include "ROL_PartitionedVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Types.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

typedef double RealT;


template<class Real>
void print_vector( const ROL::Vector<Real> &x ) {

  typedef ROL::Vector<Real>            V;
  typedef ROL::StdVector<Real>         SV;
  typedef ROL::PartitionedVector<Real> PV;
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
  typedef ROL::PartitionedVector<RealT> PV;

  GlobalMPISession mpiSession(&argc, &argv);

  int iprint = argc - 1;

  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // no output
 
  if( iprint>0 ) 
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;

  RealT errtol = ROL::ROL_THRESHOLD<RealT>();

  try {

    PV::size_type nvec = 3;    
    std::vector<int> dim;

    dim.push_back(4);
    dim.push_back(3);
    dim.push_back(5);

    int total_dim = 0;
     
    RealT left = -1e0, right = 1e0;

    std::vector<ROL::Ptr<V> > x_ptr;
    std::vector<ROL::Ptr<V> > y_ptr;
    std::vector<ROL::Ptr<V> > z_ptr;

    for( PV::size_type k=0; k<nvec; ++k ) {
      ROL::Ptr<std::vector<RealT> > xk_ptr = ROL::makePtr<std::vector<RealT>>(dim[k]);
      ROL::Ptr<std::vector<RealT> > yk_ptr = ROL::makePtr<std::vector<RealT>>(dim[k]);
      ROL::Ptr<std::vector<RealT> > zk_ptr = ROL::makePtr<std::vector<RealT>>(dim[k]);
       
      for( int i=0; i<dim[k]; ++i ) {
        (*xk_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
        (*yk_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
        (*zk_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      }
   
      ROL::Ptr<V> xk = ROL::makePtr<SV>( xk_ptr );
      ROL::Ptr<V> yk = ROL::makePtr<SV>( yk_ptr );
      ROL::Ptr<V> zk = ROL::makePtr<SV>( zk_ptr );

      x_ptr.push_back(xk);
      y_ptr.push_back(yk);
      z_ptr.push_back(zk);
      
      total_dim += dim[k];
    }

    PV x(x_ptr);
    ROL::Ptr<V> y = ROL::CreatePartitionedVector<RealT>(y_ptr[0],y_ptr[1],y_ptr[2]);
    PV z(z_ptr);

    // Standard tests.
    std::vector<RealT> consistency = x.checkVector(*y, z, true, *outStream);
    ROL::StdVector<RealT> checkvec( ROL::makePtrFromRef(consistency) );
    if (checkvec.norm() > std::sqrt(errtol)) {
      errorFlag++;
    }

    // Basis tests.
    // set x to first basis vector
    ROL::Ptr<ROL::Vector<RealT> > zp = x.clone();
    zp = x.basis(0);
    RealT znorm = zp->norm();
    *outStream << "Norm of ROL::Vector z (first basis vector): " << znorm << "\n";
    if ( std::abs(znorm-1.0) > errtol ) {
      *outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    };
    // set x to middle basis vector
    zp = x.basis(total_dim/2);
    znorm = zp->norm();
    *outStream << "\nNorm of ROL::Vector z ('middle' basis vector): " << znorm << "\n";
    if ( std::abs(znorm-1.0) > errtol ) {
      *outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    };
    // set x to last basis vector
    zp = x.basis(total_dim-1);
    znorm = zp->norm();
    *outStream << "\nNorm of ROL::Vector z (last basis vector): " << znorm << "\n";
    if ( std::abs(znorm-1.0) > errtol ) {
      *outStream << "---> POSSIBLE ERROR ABOVE!\n";
      errorFlag++;
    };

    // Repeat the checkVector tests with a zero vector.
    x.scale(0.0);
    consistency = x.checkVector(x, x, true, *outStream);
    if (checkvec.norm() > 0.0) {
      errorFlag++;
    }

    if(argc>1) {
      int m = atoi(argv[1]);
      print_vector(*(x.basis(m)));

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
