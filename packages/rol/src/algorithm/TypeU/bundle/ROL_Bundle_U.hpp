// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BUNDLE_U_H
#define ROL_BUNDLE_U_H

#include "ROL_Vector.hpp"
#include "ROL_Ptr.hpp"
#include <vector>
#include <set>

/** \class ROL::Bundle_U
    \brief Provides the interface for and implements a bundle.
*/

namespace ROL {

template<typename Real>
class Bundle_U {
/***********************************************************************************************/
/***************** BUNDLE STORAGE **************************************************************/
/***********************************************************************************************/
private: 
  std::vector<Ptr<Vector<Real>>> subgradients_;
  std::vector<Real> linearizationErrors_;
  std::vector<Real> distanceMeasures_;

  std::vector<Real> dualVariables_;

  Ptr<Vector<Real>> tG_;
  Ptr<Vector<Real>> eG_;
  Ptr<Vector<Real>> yG_;
  Ptr<Vector<Real>> gx_;
  Ptr<Vector<Real>> ge_;

  unsigned size_;

  unsigned maxSize_;
  unsigned remSize_;
  Real coeff_;
  Real omega_;

  bool isInitialized_;

  void remove(const std::vector<unsigned> &ind);

  void add(const Vector<Real> &g, const Real le, const Real dm);
  
/***********************************************************************************************/
/***************** BUNDLE MODIFICATION AND ACCESS ROUTINES *************************************/
/***********************************************************************************************/
public:
  virtual ~Bundle_U(void) {}

  Bundle_U(const unsigned maxSize = 10,
           const Real coeff = 0.0,
           const Real omega = 2.0,
           const unsigned remSize = 2);

  virtual void initialize(const Vector<Real> &g);

  virtual unsigned solveDual(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8) = 0;

  const Real linearizationError(const unsigned i) const;

  const Real distanceMeasure(const unsigned i) const;

  const Vector<Real> & subgradient(const unsigned i) const;
  
  const Real getDualVariable(const unsigned i) const;
  
  void setDualVariable(const unsigned i, const Real val);

  void resetDualVariables(void);

  const Real computeAlpha(const Real dm, const Real le) const;

  const Real alpha(const unsigned i) const;

  unsigned size(void) const;

  void aggregate(Vector<Real> &aggSubGrad, Real &aggLinErr, Real &aggDistMeas) const;

  void reset(const Vector<Real> &g, const Real le, const Real dm);

  void update(const bool flag, const Real linErr, const Real distMeas,
              const Vector<Real> &g, const Vector<Real> &s);

protected:
  const Real GiGj(const unsigned i, const unsigned j) const;

  const Real dotGi(const unsigned i, const Vector<Real> &x) const;

  void addGi(const unsigned i, const Real a, Vector<Real> &x) const;

  Real evaluateObjective(std::vector<Real> &g, const std::vector<Real> &x, const Real t) const;

  unsigned solveDual_dim1(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8);

  unsigned solveDual_dim2(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8);

}; // class Bundle_U 

} // namespace ROL

#include "ROL_Bundle_U_Def.hpp"

#endif
