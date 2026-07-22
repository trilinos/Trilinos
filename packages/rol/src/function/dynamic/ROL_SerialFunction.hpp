// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_SERIALFUNCION_HPP
#define ROL_SERIALFUNCION_HPP

#include "ROL_PartitionedVector.hpp"
#include "ROL_VectorWorkspace.hpp"
#include "ROL_TimeStamp.hpp"

/** @ingroup func_group
    \class ROL::SerialFunction
    \brief Provides behavior common to SerialObjective as SerialConstaint
*/

namespace ROL {

template<typename Real>
class SerialFunction {
public:
  using size_type = typename std::vector<Real>::size_type;

private:
  using PV = PartitionedVector<Real>;

  Ptr<Vector<Real>>             u_initial_;                // Initial condition 
  Ptr<Vector<Real>>             u_zero_;                   // Zero state vector
  TimeStampsPtr<Real>           timeStampsPtr_;            // Set of all time stamps
  mutable VectorWorkspace<Real> workspace_;                // Memory bank for vector clones
  size_type                     Nt_;                       // Number of type steps
  bool                          skipInitialCond_ = false;

protected:

  const TimeStamp<Real>& ts( size_type i ) const { return timeStampsPtr_->at(i); }
  Ptr<Vector<Real>> clone( const Vector<Real>& x ) { return workspace_.clone(x); }

public:

  SerialFunction( const Vector<Real>& u_initial,
                  const TimeStampsPtr<Real>& timeStampsPtr ) : 
    u_initial_(u_initial.clone()), 
    u_zero_(u_initial.clone()), 
    timeStampsPtr_(timeStampsPtr),
    Nt_(timeStampsPtr->size()) {
    u_initial_->set(u_initial);
    u_zero_->zero();
  }

  size_type numTimeSteps() const { return Nt_; }

  const Vector<Real>& getInitialCondition() const { return *u_initial_; }
  void setInitialCondition( const Vector<Real>& u_initial ) { u_initial_->set(u_initial); }

  const Vector<Real>& getZeroState() const { return *u_zero_; }

  bool getSkipInitialCondition() const { return skipInitialCond_; }
  void setSkipInitialCondition( bool skip ) { skipInitialCond_ = skip; }

  // Get and set methods for all TimeStamps
  TimeStampsPtr<Real> getTimeStampsPtr() const { return timeStampsPtr_; }

  void setTimeStampsPtr( const TimeStampsPtr<Real>& timeStampsPtr ) {
    timeStampsPtr_ = timeStampsPtr; 
    Nt_ = timeStampsPtr_->size(); 
  }

  // Get and set methods for individual TimeStamp objects 
  TimeStamp<Real>& getTimeStamp( size_type i ) { return timeStampsPtr_->at(i); }
  const TimeStamp<Real>& getTimeStamp( size_type i ) const { return timeStampsPtr_->at(i); }

  void setTimeStamp( size_type i, const TimeStamp<Real>& timeStamp ) {
    timeStampsPtr_->at(i) = timeStamp;
  }


}; // class SerialFunction


} // namespace ROL

#endif  // ROL_SERIALFUNCION_HPP
