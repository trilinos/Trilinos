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

#include "ROL_BlockOperator.hpp"
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

  const PV eb = Teuchos::dyn_cast<const PV>(x);
  size_type n = eb.numVectors();
    
  for(size_type k=0; k<n; ++k) {
    std::cout << "[subvector " << k << "]" << std::endl;
    Teuchos::RCP<const V> vec = eb.get(k);
    Teuchos::RCP<const std::vector<Real> > vp = 
      Teuchos::dyn_cast<SV>(const_cast<V&>(*vec)).getVector();  
   for(size_type i=0;i<vp->size();++i) {
      std::cout << (*vp)[i] << std::endl;
    }  
  }
}

// Implemantation of a Dyadic operator x*y'
template<class Real> 
class DyadicOperator : public ROL::LinearOperator<Real> {
  
  typedef ROL::Vector<Real> V;

private:

  const Teuchos::RCP<const V> x_;
  const Teuchos::RCP<const V> y_;

public:
  
  DyadicOperator( const Teuchos::RCP<const V> &x,
                  const Teuchos::RCP<const V> &y ) : x_(x), y_(y) {}

  void apply( V &Hv, const V &v, Real &tol ) const {
    Hv.set(*x_);
    Hv.scale(v.dot(*y_));  
  }
    
};

// Implementation of a Null operator
template<class Real> 
class NullOperator : public ROL::LinearOperator<Real> {

  typedef ROL::Vector<Real> V;

public:

  void apply( V &Hv, const V &v, Real &tol ) const {
    Hv.zero();
  }
};




typedef double RealT;

int main(int argc, char *argv[]) {

  // Define vectors
  typedef std::vector<RealT>            vector;
  typedef ROL::Vector<RealT>            V;
  typedef ROL::StdVector<RealT>         SV;

  // Define operators
  typedef ROL::LinearOperator<RealT>    LinOp;
  typedef ROL::DiagonalOperator<RealT>  DiagOp;
  typedef DyadicOperator<RealT>         DyadOp;
  typedef NullOperator<RealT>           NullOp;

  typedef typename vector::size_type    uint;

  using namespace Teuchos;

  GlobalMPISession mpiSession(&argc, &argv);

  int iprint = argc - 1;

  RCP<std::ostream> outStream;
  oblackholestream bhs; // no output
 
  if( iprint>0 ) 
    outStream = rcp(&std::cout,false);
  else
    outStream = rcp(&bhs,false);

  int errorFlag = 0;

  RealT errtol = ROL::ROL_THRESHOLD<RealT>();

  try {

    uint dim   = 3;  // Number of elements in each subvector (could be different)
 
    RCP<vector> x1_rcp = rcp( new vector(dim,1.0) );
    RCP<vector> x2_rcp = rcp( new vector(dim,2.0) );
 
    RCP<vector> y1_rcp = rcp( new vector(dim,0.0) );
    RCP<vector> y2_rcp = rcp( new vector(dim,0.0) );

    RCP<vector> z1_rcp = rcp( new vector(dim,0.0) );
    RCP<vector> z2_rcp = rcp( new vector(dim,0.0) );

    RCP<V> x1 = rcp( new SV( x1_rcp) );
    RCP<V> x2 = rcp( new SV( x2_rcp) );

    RCP<V> y1 = rcp( new SV( y1_rcp) );
    RCP<V> y2 = rcp( new SV( y2_rcp) );

    RCP<V> z1 = rcp( new SV( z1_rcp) );
    RCP<V> z2 = rcp( new SV( z2_rcp) ); 

    RCP<V> x = ROL::CreatePartitionedVector( x1, x2 );
    RCP<V> y = ROL::CreatePartitionedVector( y1, y2 );
    RCP<V> z = ROL::CreatePartitionedVector( z1, z2 );

    // Operator diagonals
    RCP<vector> d1_rcp = rcp( new vector(dim,0.0) );
    RCP<vector> d2_rcp = rcp( new vector(dim,0.0) );

    // Dyadic components
    RCP<vector> u_rcp = rcp( new vector(dim,0.0) );
    RCP<vector> v_rcp = rcp( new vector(dim,1.0) );

    (*d1_rcp)[0] = 6.0;   (*d2_rcp)[0] = 3.0;
    (*d1_rcp)[1] = 5.0;   (*d2_rcp)[1] = 2.0;
    (*d1_rcp)[2] = 4.0;   (*d2_rcp)[2] = 1.0;

    (*z1_rcp)[0] = 6.0;   (*z2_rcp)[0] = 6.0;    
    (*z1_rcp)[1] = 11.0;  (*z2_rcp)[1] = 4.0;
    (*z1_rcp)[2] = 4.0;   (*z2_rcp)[2] = 2.0;

    (*u_rcp)[1] = 1.0;    

    RCP<V> d1 = rcp( new SV(d1_rcp) );
    RCP<V> d2 = rcp( new SV(d2_rcp) );
    RCP<V> u  = rcp( new SV(u_rcp) );
    RCP<V> v  = rcp( new SV(v_rcp) );
    
    RCP<LinOp> D1 = rcp( new DiagOp(*d1) );
    RCP<LinOp> NO = rcp( new NullOp() );
    RCP<LinOp> UV = rcp( new DyadOp(u,v) );
    RCP<LinOp> D2 = rcp( new DiagOp(*d2) );

   
    RealT tol = 0.0;

    D1->apply(*x1,*x1,tol);
    D1->applyInverse(*x1,*x1,tol);

    ROL::BlockOperator2<RealT> bkop(D1,UV,NO,D2);


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

