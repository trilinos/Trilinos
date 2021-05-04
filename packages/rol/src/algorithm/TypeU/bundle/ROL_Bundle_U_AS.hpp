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
