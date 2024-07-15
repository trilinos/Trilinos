#pragma once
#ifndef ROL_VECTORCLONE_HPP
#define ROL_VECTORCLONE_HPP

// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ROL_Vector.hpp"
#include <exception>
#include <typeinfo>
#include <utility>
#include <map>

/** @ingroup la_group
    \class ROL::VectorClone
    \brief Container for wrapping a reusable cloned vector. Declaring an
           object of this type as a class member variable will decrease
           the number of clones needed as memory need only be allocated once
           in the lifetime of the host object. Verifies that member and argument
           types and dimensions agree when called.
*/

namespace ROL {

namespace details {

using namespace std;

template<typename Real>
class VectorClone {
private:

  Ptr<Vector<Real>> vec_;
  bool is_allocated_;

public:

  VectorClone() : vec_(nullPtr), is_allocated_(false) {}

  Ptr<Vector<Real>> operator() ( const Vector<Real>& x ) {
    if( is_allocated_ ) {
      if( typeid(x) != typeid(*vec_) )
        throw logic_error("Argument and member vector types are different!");
      if( x.dimension() != vec_->dimension() )
        throw logic_error("Argument and member vector types have different dimensions!");
    }
    else {
      vec_ = x.clone();
      is_allocated_ = true;
    }
    return vec_;
  }

  Ptr<Vector<Real>> operator() ( const Ptr<const Vector<Real>>& x ) {
    if( is_allocated_ ) {
      if( typeid(*x) != typeid(*vec_) )
        throw logic_error("Argument and member vector types are different!");
      if( x->dimension() != vec_->dimension() )
        throw logic_error("Argument and member vector types have different dimensions!");
    }
    else {
      vec_ = x->clone();
      is_allocated_ = true;
    }
    return vec_;
  }
}; // VectorClone



/** @ingroup la_group
    \class ROL::VectorCloneMap
    \brief Container for wrapping a collection of uniquely-named reusable cloned vectors,
           which in are stored in a map. Uses string-valued ids for keys by default. 
*/


template<typename Real, typename KeyType=const char*>
class VectorCloneMap {
private:
  map<KeyType, VectorClone<Real>> clones_;
 
  template<typename First, typename...Rest>
  void Constructor_Impl( First first, Rest... rest ) {
    clones_[static_cast<KeyType>(first)] = VectorClone<Real>();
    Constructor_Impl( rest... );
  }

  template<typename First>
  void Constructor_Impl( First first ) {
    clones_[static_cast<KeyType>(first)] = VectorClone<Real>();
  }

public:

  /** \brief Preallocate keys if desired */
  template<typename... Keys>
  VectorCloneMap( Keys&&...keys ) {
    Constructor_Impl( forward<Keys>(keys)... );
  }

  Ptr<Vector<Real>> operator() ( const Vector<Real>& x, KeyType key ) {
    return clones_[key](x);
  }

   Ptr<Vector<Real>> operator() ( const Ptr<const Vector<Real>>& x, KeyType key ) {
    return clones_[key](x);
  }
}; // VectorCloneMap





} // namespace details

using details::VectorClone;
using details::VectorCloneMap;

} // namespace ROL


#endif // ROL_VECTORCLONE_HPP

