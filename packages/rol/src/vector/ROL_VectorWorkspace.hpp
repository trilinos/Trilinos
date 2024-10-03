// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_VECTORWORKSPACE_HPP
#define ROL_VECTORWORKSPACE_HPP

#include "ROL_Vector.hpp"
#include <iostream>
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

           NOTE: Stored clones will have a reference count of 2 when there
           are no external pointers to the same object.

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
  using size_type = typename vector<Real>::size_type;

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

    static string to_string( const VectorKey& key ) {
      stringstream ss;
      ss << "VectorKey(" << hex << key.hash_code << ","
                         << dec << key.dimension << ")";
      return ss.str();
    }

    bool operator < ( const VectorKey& x ) const {
      return ( hash_code < x.hash_code ) && ( dimension < x.dimension );
    }

    bool operator == ( const VectorKey& x ) const {
      return ( hash_code == x.hash_code ) && ( dimension == x.dimension );
    }

//    bool operator != ( const VectorKey& x ) const {
//      return ( hash_code != x.hash_code ) || ( dimension != x.dimension );
//    }

  }; // class VectorKey

  struct VectorStack {
  
    friend class VectorWorkspace<Real>;
    vector<Ptr<V>> vectors_;
    VectorKey key_;

    VectorStack( const V& x ) : vectors_( 1, x.clone() ),
      key_(VectorKey(x)) {}

    const VectorKey& getKey() const { return key_; }

    size_type size() const { return vectors_.size(); }

    size_type number_assigned() const {
      size_type count = 0;
      for( auto v : vectors_ ) count += ( getCount(v) > 3 );
      return count;
    }

    /** If no next element exists, clone it, increment the index, and
        return a the clone by pointer
    */
    Ptr<V> clone( const V& x ) {
      VectorKey x_key(x);
      
      ROL_TEST_FOR_EXCEPTION( key_.hash_code != x_key.hash_code, logic_error,
        "VectorWorkspace::VectorStack tried to clone a vector of type "     <<
        hex << key_.hash_code << ", but it can only clone vectors of type " <<
        hex << x_key.hash_code );

      ROL_TEST_FOR_EXCEPTION( key_.dimension != x_key.dimension, logic_error, 
        "VectorWorkspace::VectorStack tried to clone a vector of dimension "     <<
        hex << key_.dimension << ", but it can only clone vectors of dimension " <<
        hex << x_key.dimension );

      for( auto e : vectors_ ) { // Return first unreferenced vector
        if( getCount(e) <= 2 ) { // Storing pointers in vector increments count  
          return e;
        }
      }   
      // If no unreferenced vectors exist, add a new one
      auto v = x.clone();
      vectors_.push_back( v );
      return v;
    }

    // For testing purposes
    vector<size_type> getRefCounts( void ) const {
      vector<size_type> counts;
      for( auto e: vectors_ ) counts.push_back( getCount(e) );
      return counts;
    }

  }; // VectorStack

  map<VectorKey,Ptr<VectorStack>> workspace_;
    
public:

  Ptr<V> clone( const V& x ) {  
   
    VectorKey        key(x);
    size_type        key_count{0};
    Ptr<VectorStack> vstack{nullPtr};

    for( auto e : workspace_ ) key_count += (key == e.first);

    if( key_count == 0 ) { // New key
      vstack = makePtr<VectorStack>(x);
      workspace_.insert( make_pair(key,vstack) );
    }
    else vstack = workspace_[key];

    return vstack->clone(x);
  }
 
  Ptr<V> clone( const Ptr<const V>& x ) { return clone(*x); }

  // Deep copy
  Ptr<V> copy( const V& x ) { 
    auto xc = clone(x);
    xc->set(x);
    return xc;
  }

  Ptr<V> copy( const Ptr<const V>& x ) { return copy(*x); }
  
  void status( ostream& os ) const {
    os << "\n\n" << string(80,'-') << std::endl;
    os << "VectorWorkspace contains the following VectorStack(hash_code,dim) entries:\n\n";
    for( auto entry : workspace_ ) {
      os << "  VectorStack(" << hex << entry.first.hash_code << ","
                               << dec << entry.first.dimension << ")";
      os << "\n  Reference Counts per element" << std::endl;
      for( auto e : entry.second->vectors_ ) {
        os << "        " << getCount( e ) << std::endl;
      }
    }
    os << string(80,'-') << std::endl;
  }


}; // VectorWorkspace

} // namespace details 

using details::VectorWorkspace;

} // namespace ROL


#endif 
