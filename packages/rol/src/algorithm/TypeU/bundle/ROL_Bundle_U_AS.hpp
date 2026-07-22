// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BUNDLE_U_AS_H
#define ROL_BUNDLE_U_AS_H

#include "ROL_Bundle_U.hpp"

/** \class ROL::Bundle_U_AS
    \brief Provides the interface for and implements an active set bundle.
*/

namespace ROL {

template<typename Real>
class Bundle_U_AS : public Bundle_U<Real> {
/***********************************************************************************************/
/***************** BUNDLE STORAGE **************************************************************/
/***********************************************************************************************/
private:

  Ptr<Vector<Real>> tG_;
  Ptr<Vector<Real>> eG_;
  Ptr<Vector<Real>> yG_;
  Ptr<Vector<Real>> gx_;
  Ptr<Vector<Real>> ge_;

  std::set<unsigned> workingSet_;
  std::set<unsigned> nworkingSet_;

  bool isInitialized_;
  
/***********************************************************************************************/
/***************** BUNDLE MODIFICATION AND ACCESS ROUTINES *************************************/
/***********************************************************************************************/
public:
  Bundle_U_AS(const unsigned maxSize = 10,
              const Real coeff = 0.0,
              const Real omega = 2.0,
              const unsigned remSize = 2);

  void initialize(const Vector<Real> &g);

  unsigned solveDual(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8);

/***********************************************************************************************/
/***************** DUAL CUTTING PLANE PROBLEM ROUTINES *****************************************/
/***********************************************************************************************/
private:
  void initializeDualSolver(void);

  void computeLagMult(std::vector<Real> &lam, const Real mu, const std::vector<Real> &g) const;
 
  bool isNonnegative(unsigned &ind, const std::vector<Real> &x) const;

  Real computeStepSize(unsigned &ind, const std::vector<Real> &x, const std::vector<Real> &p) const;

  unsigned solveEQPsubproblem(std::vector<Real> &s, Real &mu,
                        const std::vector<Real> &g, const Real tol) const;

  void applyPreconditioner(std::vector<Real> &Px, const std::vector<Real> &x) const;

  void applyG(std::vector<Real> &Gx, const std::vector<Real> &x) const;

  void applyPreconditioner_Identity(std::vector<Real> &Px, const std::vector<Real> &x) const;

  void applyG_Identity(std::vector<Real> &Gx, const std::vector<Real> &x) const;

  void applyPreconditioner_Jacobi(std::vector<Real> &Px, const std::vector<Real> &x) const;

  void applyG_Jacobi(std::vector<Real> &Gx, const std::vector<Real> &x) const;

  void applyPreconditioner_SymGS(std::vector<Real> &Px, const std::vector<Real> &x) const;

  void applyG_SymGS(std::vector<Real> &Gx, const std::vector<Real> &x) const;

  void computeResidualUpdate(std::vector<Real> &r, std::vector<Real> &g) const;

  void applyFullMatrix(std::vector<Real> &Hx, const std::vector<Real> &x) const;

  void applyMatrix(std::vector<Real> &Hx, const std::vector<Real> &x) const;

  unsigned projectedCG(std::vector<Real> &x, Real &mu, const std::vector<Real> &b, const Real tol) const;

  Real dot(const std::vector<Real> &x, const std::vector<Real> &y) const;

  Real norm(const std::vector<Real> &x) const;

  void axpy(const Real a, const std::vector<Real> &x, std::vector<Real> &y) const;

  void scale(std::vector<Real> &x, const Real a) const;

  void scale(std::vector<Real> &x, const Real a, const std::vector<Real> &y) const;

  unsigned solveDual_arbitrary(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8);

  /************************************************************************/
  /********************** PROJECTION ONTO FEASIBLE SET ********************/
  /************************************************************************/
  void project(std::vector<Real> &x, const std::vector<Real> &v) const;

  Real computeCriticality(const std::vector<Real> &g, const std::vector<Real> &sol) const;
}; // class Bundle_AS

} // namespace ROL

#include "ROL_Bundle_U_AS_Def.hpp"

#endif
