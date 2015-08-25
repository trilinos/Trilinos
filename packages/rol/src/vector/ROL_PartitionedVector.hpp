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

#include "ROL_Vector.hpp"

#ifndef ROL_PARTITIONED_VECTOR_H
#define ROL_PARTITIONED_VECTOR_H

/** @ingroup la_group
 *  \class ROL::PartitionedVector
 *  \brief Defines the linear algebra of vector space on a generic partitioned vector
 */
 

namespace ROL {

template<class Real>
class PartitionedVector : public Vector<Real> {

  typedef Vector<Real>                  V;
  typedef Teuchos::RCP<V>               RCPV;
  typedef PartitionedVector<Real>       PV;
  typedef std::vector<Real>::size_type  size_type;

private:
  std::vector<RCPV> vecs_;
  mutable std::vector<RCPV> dual_vecs_;
  mutable Teuchos::RCP<PartitionedVector<Real> > dual_pvec_;

public:
  PartitionedVector( const std::vector<RCPV> &vecs ) : vecs_(vecs) {
   
    for( size_type i=0; i<vecs.size(); ++i ) {
      dual_vecs[i] = (vecs[i]->dual()).clone();  
    }
  }

  void plus( const V &x ) {
    using Teuchos::dyn_cast;
    const PV &xs = dyn_cast<const PV>(dyn_cast<const V>(x));
    
    for( size_type i=0; i<vecs.size(); ++i ) { 
      vecs_[i]->plus(*xs.get(i));
    }
  }

  void scale( const Real alpha ) {
    for( size_type i=0; i<vecs.size(); ++i ) { 
      vecs_[i]->scale(alpha);
    }
  }

  void axpy( const Real alpha, const V &x ) {
    using Teuchos::dyn_cast;
    const PV &xs = dyn_cast<const PV>(dyn_cast<const V>(x));
    for( size_type i=0; i<vecs.size(); ++i ) { 
      vecs_[i]->axpy(alpha,xs.get(i));
    }
  } 
 
  Real dot( const V &x ) const {
    const PV &xs = dyn_cast<const PV>(dyn_cast<const V>(x));
    Real result = 0; 
      for( size_type i=0; i<vecs.size(); ++i ) { 
        result += vecs_[i]->dot(xs.get(i));
      }
    return result;  
  }
 
  Real norm() const {
    Real result = 0; 
      for( size_type i=0; i<vecs_.size(); ++i ) {   
        result += vecs_[i]->(vecs_[i]);
      }
    return result;  
  }

  RCPV clone() const {
    std::vector<RCPV> clonevec;
    
    for( size_type i=0; i<vecs_.size(); ++i ) {   
      clonevec.push_back(vecs_[i]->clone());
    }
    return Teuchos::rcp( new PV(clonevec) );
  }

  const V& dual(void) const {
    for( size_type i=0; i<vecs_.size(); ++i ) {  
      dual_vecs_[i]->set(vecs_[i]); 
    }
    return Teuchos::rcp( new PV(dual_vecs_) );
  }

  RCPV basis( const int i ) const {
    std::vector<RCPV> basisvecs;
    int begin = 0;   
    int end = 0;
    for( size_type j=0; j<vecs_.size(); ++j ) { 
      end += vecs_[j].dimension();
        
      if( begin<= i && i<end ) {
        RCPV e = vecs_[j]->basis(i-begin);    
        basisvecs.push_back(e);  
      }
      else {
        RCPV e = vecs_[j]->clone();
        e->zero();
        basisvecs.push_back(e);     
      }
      begin = end+1;      
    } 
  }

  int dimension() const {
    int total_dim = 0;
    for( size_type j=0; j<vecs_.size(); ++j ) { 
      total_dim += vecs_[j].dimension(); 
    }
    return total_dim;
  }

  RCPV get(size_type i) const {
    return vecs_[i];
  }

  void set(size_type i, const V &x) {
    vecs_[i]->set(x); 
  }

};

} // namespace ROL

#endif // ROL_PARTITIONED_VECTOR_H

