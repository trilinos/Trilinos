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

#ifndef ROL_BUNDLE_U_TT_H
#define ROL_BUNDLE_U_TT_H

#include "ROL_Bundle_U.hpp"
#include "ROL_LAPACK.hpp" 
#include "ROL_LinearAlgebra.hpp"
#include <vector>

/** \class ROL::Bundle_U_TT
    \brief Provides the interface for and implements a bundle. The semidefinite
           quadratic subproblem is solved using TT algorithm by Antonio
           Frangioni (1996).
*/

namespace ROL {

template<typename Real>
class Bundle_U_TT : public Bundle_U<Real> {
private: 
  LAPACK<int, Real> lapack_; // TT

  int QPStatus_;           // QP solver status
  int maxind_;             // maximum integer value
  int entering_;           // index of entering item
  int LiMax_;              // index of max element of diag(L)
  int LiMin_;              // index of min element of diag(L)

  unsigned maxSize_;       // maximum bundle size
  unsigned dependent_;     // number of lin. dependent items in base
  unsigned currSize_;      // current size of base

  bool isInitialized_;
  bool optimal_;           // flag for optimality of restricted solution

  Real rho_;
  Real lhNorm;
  Real ljNorm;
  Real lhz1_;
  Real lhz2_;
  Real ljz1_;
  Real kappa_;             // condition number of matrix L ( >= max|L_ii|/min|L_ii| )
  Real objval_;            // value of objective
  Real minobjval_;         // min value of objective (ever reached)
  Real deltaLh_;           // needed in case dependent row becomes independent
  Real deltaLj_;           // needed in case dependent row becomes independent

  std::vector<int> taboo_; // list of "taboo" items
  std::vector<int> base_;  // base

  LA::Matrix<Real> L_;
  LA::Matrix<Real> Id_;
  LA::Vector<Real> tempv_;
  LA::Vector<Real> tempw1_;
  LA::Vector<Real> tempw2_;
  LA::Vector<Real> lh_;
  LA::Vector<Real> lj_;
  LA::Vector<Real> z1_;
  LA::Vector<Real> z2_;

  using Bundle_U<Real>::GiGj;

public:
  Bundle_U_TT(const unsigned maxSize = 10,
              const Real coeff = 0.0,
              const Real omega = 2.0,
              const unsigned remSize = 2);
  
  unsigned solveDual(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8);

/***********************************************************************************************/
/****************** DUAL CUTTING PLANE SUBPROBLEM ROUTINES *************************************/
/***********************************************************************************************/
private:
  Real sgn(const Real x) const;

  void swapRowsL(unsigned ind1, unsigned ind2, bool trans=false);

  void updateK(void);

  void addSubgradToBase(unsigned ind, Real delta);

  void deleteSubgradFromBase(unsigned ind, Real tol);
  
  // TT: solving triangular system for TT algorithm
  void solveSystem(int size, char tran, LA::Matrix<Real> &L, LA::Vector<Real> &v);

  // TT: check that inequality constraints are satisfied for dual variables
  bool isFeasible(LA::Vector<Real> &v, const Real &tol);

  unsigned solveDual_TT(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8);

  unsigned solveDual_arbitrary(const Real t, const unsigned maxit = 1000, const Real tol = 1.e-8);

}; // class Bundle_TT

} // namespace ROL

#include "ROL_Bundle_U_TT_Def.hpp"

#endif

