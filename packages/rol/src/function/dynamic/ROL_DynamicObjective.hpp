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
#ifndef ROL_DYNAMICOBJECTIVE_HPP
#define ROL_DYNAMICOBJECTIVE_HPP

#include "ROL_DynamicFunction.hpp"
#include "ROL_TimeStamp.hpp"


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
class DynamicObjective {
public:

  using V  = Vector<Real>;
  using PV = PartitionedVector<Real>;
  using TS = TimeStamp<Real>;

  virtual ~DynamicObjective() {}

  virtual void update( const V& uo, const V& un, const V& z ) {
    update_uo( uo );
    update_un( un );
    update_z( z );
  }

  virtual void update_uo( const V& uo ) { }
  virtual void update_un( const V& un ) { }
  virtual void update_z( const V& z ) { }


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
  virtual void hessVec_uo_uo( V& hv, const V& vo, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const {}

  virtual void hessVec_uo_un( V& hv, const V& vo, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const {}

  virtual void hessVec_uo_z( V& hv, const V& vo, const V& uo, const V& un, 
                             const V& z, const TS& timeStamp ) const {}

  
  virtual void hessVec_un_uo( V& hv, const V& vo, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const {}

  virtual void hessVec_un_un( V& hv, const V& vo, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const {}

  virtual void hessVec_un_z( V& hv, const V& vo, const V& uo, const V& un, 
                             const V& z, const TS& timeStamp ) const {}

  virtual void hessVec_z_uo( V& hv, const V& vo, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const {}

  virtual void hessVec_z_un( V& hv, const V& vo, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const {}

  virtual void hessVec_z_z( V& hv, const V& vo, const V& uo, const V& un, 
                             const V& z, const TS& timeStamp ) const {}
};

} // namespace ROL


#endif // ROL_DYNAMICOBJECTIVE_HPP

