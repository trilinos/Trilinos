// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_FINITEDIFFERENCE_HPP
#define ROL_FINITEDIFFERENCE_HPP


/** \class ROL::FiniteDifference
    \brief General interface for numerical approximation of Constraint 
           or Objective methods using the finite difference method. 

           TODO: Should optionally store values for reuse as needed
*/

#include "ROL_FunctionBindings.hpp"

#include "ROL_Vector.hpp"
#include "ROL_VectorWorkspace.hpp"

namespace ROL {

namespace details {

using namespace std;

template<typename Real>
class FiniteDifference {
public:
 
  using V = ROL::Vector<Real>;

  FiniteDifference( const int order = 1 );
  FiniteDifference( const int order, const Ptr<VectorWorkspace<Real>>& workspace );

  virtual ~FiniteDifference(){}

  /** Approximately compute the derivative of f(x) in the direction v 
      using the step size h */
  virtual Real operator()( f_scalar_t<Real>& f_value,
                           f_update_t<Real>& f_update,
                           const V& v, 
                           const V& x,
                           const Real h ) const;

  /** Approximately compute the Jacobian of f(x) applied to the direction v
      using the step size h */
  virtual void operator()( f_vector_t<Real>& f_value, 
                           f_update_t<Real>& f_update,
                           V& Jv, 
                           const V& v, 
                           const V& x,
                           const Real h ) const;             
private:

  const int order_; // Finite difference order (1,2,3, or 4)
  
  mutable Ptr<VectorWorkspace<Real>> workspace_;

}; // class FiniteDifference

} // namespace details

using details::FiniteDifference;


} // namespace ROL

#include "ROL_FiniteDifferenceDef.hpp"

#endif // ROL_FINITEDIFFERENCE_HPP

