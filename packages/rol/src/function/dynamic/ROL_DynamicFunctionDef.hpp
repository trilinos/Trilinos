// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_DYNAMICFUNCTIONDEF_HPP
#define ROL_DYNAMICFUNCTIONDEF_HPP


namespace ROL {

template<typename Real>
VectorWorkspace<Real>&
DynamicFunction<Real>::getVectorWorkspace() const { return workspace_; }

template<typename Real>
PartitionedVector<Real>& 
DynamicFunction<Real>::partition( Vector<Real>& x ) const { 
  return static_cast<PartitionedVector<Real>&>(x);
}

template<typename Real>
const PartitionedVector<Real>& 
DynamicFunction<Real>::partition( const Vector<Real>& x ) const { 
  return static_cast<const PartitionedVector<Real>&>(x);
}

template<typename Real>
Vector<Real>& 
DynamicFunction<Real>::getNew( Vector<Real>& x ) const { 
  return *(partition(x).get(1));
}

template<typename Real>
const Vector<Real>& 
DynamicFunction<Real>::getNew( const Vector<Real>& x ) const { 
  return *(partition(x).get(1));
}

template<typename Real>
Vector<Real>& 
DynamicFunction<Real>::getOld( Vector<Real>& x ) const { 
  return *(partition(x).get(0));
}

template<typename Real>
const Vector<Real>& 
DynamicFunction<Real>::getOld( const Vector<Real>& x ) const { 
  return *(partition(x).get(0));
}

} // namespace ROL


#endif // ROL_DYNAMICFUNCTIONDEF_HPP

