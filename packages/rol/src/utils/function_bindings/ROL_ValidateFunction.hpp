// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_VALIDATEFUNCTION_HPP
#define ROL_VALIDATEFUNCTION_HPP

/** \file  ROL_ValidateFunction.hpp
    \brief Provides a set of tools for validating the behavior of several
           function types that are commonly used in ROL.
 
           - Finite difference check of derivatives
           - Symmetry check of linear operators
           - Adjoint consistency check for linear operators
           - Inverse identity check for linear operators
*/

#include "ROL_FiniteDifference.hpp"

namespace ROL {

namespace details {

template<typename Real>
class ValidateFunction {
public:

  using V = Vector<Real>;

  ValidateFunction( const int order = 1,
                    const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                    const int width = 20,
                    const int precision = 11,
                    const bool printToStream = true,
                    ostream& os = cout );

  virtual ~ValidateFunction(){}

  virtual vector<vector<Real>> derivative_check( f_scalar_t<Real> f_value,
                                                 f_vector_t<Real> f_derivative,
                                                 f_update_t<Real> f_update,
                                                 const V& g,
                                                 const V& v,
                                                 const V& x,
                                                 const string& label ) const;

  virtual vector<vector<Real>> derivative_check( f_vector_t<Real> f_value,
                                                 f_dderiv_t<Real> f_derivative,
                                                 f_update_t<Real> f_update,
                                                 const V& c,
                                                 const V& v,
                                                 const V& x,
                                                 const string& label ) const;

  virtual vector<Real> symmetry_check( f_dderiv_t<Real> A,
                                       f_update_t<Real> A_update,
                                       const V& u,
                                       const V& v,
                                       const V& x,
                                       const string& name="Linear Operator",
                                       const string& symbol="A" ) const;

  virtual vector<Real> adjoint_consistency_check( f_dderiv_t<Real> A,
                                                  f_dderiv_t<Real> A_adj,
                                                  f_update_t<Real> A_update,
                                                  const V& u,
                                                  const V& v,
                                                  const V& x,
                                                  const string& name="Linear Operator",
                                                  const string& symbol="A" ) const;

  virtual vector<Real> inverse_check( f_dderiv_t<Real> A,
                                      f_dderiv_t<Real> A_inv,
                                      f_update_t<Real> A_update,
                                      const V& v,
                                      const V& x,
                                      const string& name="Linear Operator",
                                      const string& symbol="A" ) const;

  virtual vector<Real> solve_check( f_solve_t<Real>  solve,
                                    f_vector_t<Real> value,
                                    f_update_t<Real> update,
                                    const V& c,
                                    const V& x,
                                    const string& name="Function" ) const;

  ostream& getStream() const;

private:

  int                    order_;         // Finite difference order
  int                    numSteps_;      // Number of evalutions of different step sizes
  int                    width_;         // For print formatting
  int                    precision_;     // Number of digits to display
  bool                   printToStream_; // False will suppress output
  vector<Real>           steps_;         // Set of step sizes of FD approximation
  string                 dashline_;      // For print format
  ostream& os_;                          // pointer to Output stream
 
  mutable Ptr<VectorWorkspace<Real>> workspace_;
  FiniteDifference<Real>             fd_;

}; // ValidateFunction

} // namespace details

using details::ValidateFunction;

} // namespace ROL

#include "ROL_ValidateFunctionDef.hpp"

#endif // ROL_CHECKFUNCTION_HPP
