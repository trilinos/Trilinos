// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_DYNAMICOBJECTIVE_HPP
#define ROL_DYNAMICOBJECTIVE_HPP

#include "ROL_DynamicFunction.hpp"


/** @ingroup func_group
    \class ROL::DynamicObjective
    \brief Defines the time-dependent objective function interface for simulation-based optimization.
           Computes time-local contributions of value, gradient, Hessian-vector product etc to a 
           larger composite objective defined over the simulation time. In contrast to other 
           objective classes Objective_TimeSimOpt has a default implementation of value which 
           returns zero, as time-dependent simulation based optimization problems may have an
           objective value which depends only on the final state of the system. 

           Often in applications, the objective is separable into a state objective and a
           control objective so the state-control cross terms in the Hessian are zero.

           A typical use case involves quadratic terms for the state and control, 
           such that all of the diagonal blocks of the Hessian are nonzero 
           and the off-diagonal blocks are zero.
 
*/

namespace ROL {

template<typename Real> 
class DynamicObjective : public DynamicFunction<Real> {
public:

  using V  = Vector<Real>;
  using TS = TimeStamp<Real>;


  DynamicObjective( std::initializer_list<std::string> zero_deriv_terms={} ) :
    DynamicFunction<Real>( zero_deriv_terms ) {}
  
  virtual ~DynamicObjective() {}

  virtual void update( const V& uo, const V& un, const V& z, const TS& timeStamp ) {
    update_uo( uo, timeStamp );
    update_un( un, timeStamp );
    update_z( z, timeStamp );
  }

  using DynamicFunction<Real>::update_uo;
  using DynamicFunction<Real>::update_un;
  using DynamicFunction<Real>::update_z; 

  virtual Real value( const V& uo, const V& un, 
                      const V& z, const TS& timeStamp ) const = 0;

  //----------------------------------------------------------------------------
  // Gradient Terms
  virtual void gradient_uo( V& g, const V& uo, const V& un, 
                            const V& z, const TS& timeStamp ) const {}

  virtual void gradient_un( V& g, const V& uo, const V& un, 
                            const V& z, const TS& timeStamp ) const {}

  virtual void gradient_z( V& g, const V& uo, const V& un, 
                           const V& z, const TS& timeStamp ) const {}

  //----------------------------------------------------------------------------
  // Hessian-Vector product terms
  virtual void hessVec_uo_uo( V& hv, const V& v, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const {}

  virtual void hessVec_uo_un( V& hv, const V& v, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const {}

  virtual void hessVec_uo_z( V& hv, const V& v, const V& uo, const V& un, 
                             const V& z, const TS& timeStamp ) const {}

  
  virtual void hessVec_un_uo( V& hv, const V& v, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const {}

  virtual void hessVec_un_un( V& hv, const V& v, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const {}

  virtual void hessVec_un_z( V& hv, const V& v, const V& uo, const V& un, 
                             const V& z, const TS& timeStamp ) const {}


  virtual void hessVec_z_uo( V& hv, const V& v, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const {}

  virtual void hessVec_z_un( V& hv, const V& v, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const {}

  virtual void hessVec_z_z( V& hv, const V& v, const V& uo, const V& un, 
                             const V& z, const TS& timeStamp ) const {}
};

} // namespace ROL


#endif // ROL_DYNAMICOBJECTIVE_HPP

