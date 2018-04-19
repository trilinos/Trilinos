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

#pragma once
#ifndef ROL_VECTORWORKSPACE_HPP
#define ROL_VECTORWORKSPACE_HPP

#include "ROL_Vector.hpp"
#include <map>
#include <utility>

/** @ingroup la_group
    \class ROL::VectorWorkspace
    \brief Provides a "smart" cloning manager to be used a member variable in
           a class and called in the member function of the same class. 
           
           VectorWorkspace::clone( const Vector& x )   and
           VectorWorkspace::clone( const Ptr<Vector>& x )    

           Will allocate new memory of a clone of x *if needed* and return 
           a pointer to the clone. A new clone is considered to be needed 
           only if these is not a previously allocated compatible vector 
           stored in the VectorWorkspace. 

           Compatibility is determined by derived type (typeid::hash_code)
           and vector dimension. Together these form a VectorKey. 
           When cloning a vector inside a member function, VectorWorkspace 
           will identify it's VectorKey type. If such a type exists in the
           database, with will then refer to the associated VectorStack that
           is specific to the VectorKey type. 

           The VectorStack will be searched for the first available dynamically
           allocated vector which has no external references to it and return
           a pointer to it. If no such vector exists, a new one will be 
           cloned and added to the stack. When the local pointers to the
           VectorStack elements go out of scope at the end of the member
           function, the reference counts are decremented and the vectors 
           become available for use again. 

           It should also be possible use a single VectorWorkspace in
           multiple levels of nested objects. 

           The concept could probably also generalize to the MPI setting 
           where node locality is an additional distinguishing qualifier
           for cloning.    

*/

namespace ROL {

namespace details {

using namespace std;


template<typename Real>
class VectorWorkspace {

  using V = ROL::Vector<Real>;
  using size_type = typename vector<Real>::type;

private:

  struct VectorKey {
    friend class VectorWorkspace<Real>;
    size_t hash_code;
    int    dimension;

    VectorKey( const V& x ) :
      hash_code(typeid(x).hash_code()), 
      dimension( x.dimension() ) {}
    VectorKey( const Ptr<V>& x ) : 
      VectorKey( *x ) {}
  }; // class VectorKey

  template<typename Real> 
  struct VectorStack {
  
    friend class VectorWorkspace<Real>;
    vector<Ptr<V>> vectors_;
   
    void initialize( const V& x ) { vectors_.push_back( x.clone() ); }
    void initialize( const Ptr<V>& x ) { vectors_.push_back( x->clone(); }

    /** If no next element exists, clone it, increment the index, and
        return a the clone by pointer
    */
    Ptr<V> clone( void ) {
      for( auto e : *vectors_ ) { // Return first unreferenced vector
        if( getCount(e) > 1 ) {   
          return e;
        }
      }   
      // If no unreferenced vectors exist, add a new one
      auto v = prototype_.clone();
      vectors_.push_back( v );
      return v;
    }

    // For testing purposes
    vector<size_type> getRefCounts( void ) const {
      vector<size_type> counts;
      for( auto e: *vectors_ ) counts.push_back( getCount(e) );
      return counts;
    }

  }; // VectorStack

  map<VectorKey,VectorStack<Real>> workspace_;
    
public:

  Ptr<V> clone( const V& x ) {  

    VectorKey key(x);
    auto search = workspace_.find( key );

    if( search != workspace_.end() ) {
      auto entry = workspace_[key];
      return entry.clone();
    }
    else { // First instance of this kind of vector
      auto entry = workspace_[key];
      entry.initialize(x);
      return entry.clone();
    }
  }
 
  Ptr<V> clone( const Ptr<V>& x ) { return clone(*x); }

}; // VectorWorkspace

} // namespace details 

using details::VectorWorkspace;

} // namespace ROL


#endif 
