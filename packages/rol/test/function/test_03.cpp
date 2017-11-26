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
    \brief Test action of a BlockOperator on a PartitionedVector 
    \details Apply a \f$ 2\times 2\f$ block operator \f$H\$ to a partitioned vector 
             \f$ x = \begin{pmatrix} x_1 & x_2 \end{pmatrix} \f$ where
             \f[ H=\begin{pmatrix} A & B \\ 0 & C \f] 

*/

#include "ROL_NullOperator.hpp"
#include "ROL_DyadicOperator.hpp"
#include "ROL_BlockOperator2.hpp"
#include "ROL_DiagonalOperator.hpp"
#include "ROL_Elementwise_Function.hpp"
#include "ROL_PartitionedVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Types.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

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
    ROL::SharedPointer<const V> vec = eb.get(k);
    ROL::SharedPointer<const std::vector<Real> > vp = 
      dynamic_cast<SV>(const_cast<V&&>(*vec)).getVector();  
   for(size_type i=0;i<vp->size();++i) {
      std::cout << (*vp)[i] << std::endl;
    }  
  }
}

typedef double RealT;

int main(int argc, char *argv[]) {

  // Define vectors
  typedef std::vector<RealT>            vector;
  typedef ROL::Vector<RealT>            V;
  typedef ROL::StdVector<RealT>         SV;

  // Define operators
  typedef ROL::LinearOperator<RealT>    LinOp;
  typedef ROL::DiagonalOperator<RealT>  DiagOp;
  typedef ROL::DyadicOperator<RealT>    DyadOp;
  typedef ROL::NullOperator<RealT>      NullOp;

  typedef typename vector::size_type    uint;

  using namespace Teuchos;

  GlobalMPISession mpiSession(&argc, &argv);

  int iprint = argc - 1;

  ROL::SharedPointer<std::ostream> outStream;
  oblackholestream bhs; // no output
 
  if( iprint>0 ) 
    outStream = ROL::makeSharedFromRef(std::cout);
  else
    outStream = ROL::makeSharedFromRef(bhs);

  int errorFlag = 0;

  RealT errtol = ROL::ROL_THRESHOLD<RealT>();

  try {

    uint dim   = 3;  // Number of elements in each subvector (could be different)
 
    ROL::SharedPointer<vector> x1_ptr = ROL::makeShared<vector>(dim,1.0);
    ROL::SharedPointer<vector> x2_ptr = ROL::makeShared<vector>(dim,2.0);
 
    ROL::SharedPointer<vector> y1_ptr = ROL::makeShared<vector>(dim,0.0);
    ROL::SharedPointer<vector> y2_ptr = ROL::makeShared<vector>(dim,0.0);

    ROL::SharedPointer<vector> z1_ptr = ROL::makeShared<vector>(dim,0.0);
    ROL::SharedPointer<vector> z2_ptr = ROL::makeShared<vector>(dim,0.0);

    ROL::SharedPointer<V> x1 = ROL::makeShared<SV>( x1_ptr);
    ROL::SharedPointer<V> x2 = ROL::makeShared<SV>( x2_ptr);

    ROL::SharedPointer<V> y1 = ROL::makeShared<SV>( y1_ptr);
    ROL::SharedPointer<V> y2 = ROL::makeShared<SV>( y2_ptr);

    ROL::SharedPointer<V> z1 = ROL::makeShared<SV>( z1_ptr);
    ROL::SharedPointer<V> z2 = ROL::makeShared<SV>( z2_ptr); 

    ROL::SharedPointer<V> x = ROL::CreatePartitionedVector( x1, x2 );
    ROL::SharedPointer<V> y = ROL::CreatePartitionedVector( y1, y2 );
    ROL::SharedPointer<V> z = ROL::CreatePartitionedVector( z1, z2 );

    // Operator diagonals
    ROL::SharedPointer<vector> d1_ptr = ROL::makeShared<vector>(dim,0.0);
    ROL::SharedPointer<vector> d2_ptr = ROL::makeShared<vector>(dim,0.0);

    // Dyadic components
    ROL::SharedPointer<vector> u_ptr = ROL::makeShared<vector>(dim,0.0);
    ROL::SharedPointer<vector> v_ptr = ROL::makeShared<vector>(dim,1.0);

    (*d1_ptr)[0] = 6.0;   (*d2_ptr)[0] = 3.0;
    (*d1_ptr)[1] = 5.0;   (*d2_ptr)[1] = 2.0;
    (*d1_ptr)[2] = 4.0;   (*d2_ptr)[2] = 1.0;

    (*z1_ptr)[0] = 6.0;   (*z2_ptr)[0] = 6.0;    
    (*z1_ptr)[1] = 11.0;  (*z2_ptr)[1] = 4.0;
    (*z1_ptr)[2] = 4.0;   (*z2_ptr)[2] = 2.0;

    (*u_ptr)[1] = 1.0;    

    ROL::SharedPointer<V> d1 = ROL::makeShared<SV>(d1_ptr);
    ROL::SharedPointer<V> d2 = ROL::makeShared<SV>(d2_ptr);
    ROL::SharedPointer<V> u  = ROL::makeShared<SV>(u_ptr);
    ROL::SharedPointer<V> v  = ROL::makeShared<SV>(v_ptr);
    
    ROL::SharedPointer<LinOp> D1 = ROL::makeShared<DiagOp>(*d1);
    ROL::SharedPointer<LinOp> NO = ROL::makeShared<NullOp>();
    ROL::SharedPointer<LinOp> UV = ROL::makeShared<DyadOp>(u,v);
    ROL::SharedPointer<LinOp> D2 = ROL::makeShared<DiagOp>(*d2);

   
    RealT tol = 0.0;

    D1->apply(*x1,*x1,tol);
    D1->applyInverse(*x1,*x1,tol);

    ROL::BlockOperator2<RealT> bkop(D1,NO,UV,D2);


    bkop.apply(*y,*x,tol);  

    z->axpy(-1.0,*y);

    errorFlag += static_cast<int>(z->norm()>errtol);

  }

  catch (std::logic_error err) {


  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";
  
  return 0; 
}

