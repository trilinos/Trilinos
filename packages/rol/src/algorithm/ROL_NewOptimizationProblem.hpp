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

#ifndef ROL_NEWOPTIMIZATIONPROBLEM_HPP
#define ROL_NEWOPTIMIZATIONPROBLEM_HPP

#include <utility>
#include <map>

#include "ROL_Ptr.hpp"
#include "ROL_Types.hpp"
#include "ROL_NewConstraintManager.hpp"
#include "ROL_SlacklessObjective.hpp"
#include "ROL_ReduceLinearConstraint.hpp"
#include "ROL_PolyhedralProjection.hpp"

namespace ROL {

template<typename Real>
class NewOptimizationProblem {
private:
  bool isFinalized_;
  bool hasBounds_;
  bool hasEquality_;
  bool hasInequality_;
  bool hasLinearEquality_;
  bool hasLinearInequality_;

  Ptr<Objective<Real>>                       INPUT_obj_;
  Ptr<Vector<Real>>                          INPUT_xprim_;
  Ptr<Vector<Real>>                          INPUT_xdual_;
  Ptr<BoundConstraint<Real>>                 INPUT_bnd_;
  std::map<std::string,ConstraintData<Real>> INPUT_econ_;
  std::map<std::string,ConstraintData<Real>> INPUT_icon_;
  std::map<std::string,ConstraintData<Real>> INPUT_linear_econ_;
  std::map<std::string,ConstraintData<Real>> INPUT_linear_icon_;

  Ptr<Objective<Real>>            obj_;
  Ptr<Vector<Real>>               xprim_;
  Ptr<Vector<Real>>               xdual_;
  Ptr<BoundConstraint<Real>>      bnd_;
  Ptr<Constraint<Real>>           con_;
  Ptr<Vector<Real>>               mul_;
  Ptr<Vector<Real>>               res_;
  Ptr<PolyhedralProjection<Real>> proj_;

  Ptr<Vector<Real>> xfeas_;
  Ptr<ReduceLinearConstraint<Real>> rlc_;

  EProblem problemType_;

public:
  virtual ~NewOptimizationProblem(void) {}

  // Default constructor
  NewOptimizationProblem(const Ptr<Objective<Real>> &obj,
                         const Ptr<Vector<Real>>    &x,
                         const Ptr<Vector<Real>>    &g = nullPtr);

  /* Set and remove methods for constraints */
  void addBoundConstraint(const Ptr<BoundConstraint<Real>> &bnd);
  void removeBoundConstraint(void);

  void addEqualityConstraint(std::string                  name,
                             const Ptr<Constraint<Real>> &econ,
                             const Ptr<Vector<Real>>     &emul,
                             const Ptr<Vector<Real>>     &eres = nullPtr,
                             bool                         reset = false);
  void removeEqualityConstraint(std::string name);

  void addInequalityConstraint(std::string                       name,
                               const Ptr<Constraint<Real>>      &icon,
                               const Ptr<Vector<Real>>          &imul,
                               const Ptr<BoundConstraint<Real>> &ibnd,
                               const Ptr<Vector<Real>>          &ires = nullPtr,
                               bool                              reset = false);
  void removeInequalityConstraint(std::string name);

  void addLinearEqualityConstraint(std::string                  name,
                                   const Ptr<Constraint<Real>> &linear_econ,
                                   const Ptr<Vector<Real>>     &linear_emul,
                                   const Ptr<Vector<Real>>     &linear_eres = nullPtr,
                                   bool                         reset = false);
  void removeLinearEqualityConstraint(std::string name);

  void addLinearInequalityConstraint(std::string                       name,
                                     const Ptr<Constraint<Real>>      &linear_icon,
                                     const Ptr<Vector<Real>>          &linear_imul,
                                     const Ptr<BoundConstraint<Real>> &linear_ibnd,
                                     const Ptr<Vector<Real>>          &linear_ires = nullPtr,
                                     bool                              reset = false);
  void removeLinearInequalityConstraint(std::string name);

  /* Get methods */
  const Ptr<Objective<Real>>            getObjective(void);
  const Ptr<Vector<Real>>               getPrimalOptimizationVector(void);
  const Ptr<Vector<Real>>               getDualOptimizationVector(void);
  const Ptr<BoundConstraint<Real>>      getBoundConstraint(void);
  const Ptr<Constraint<Real>>           getConstraint(void);
  const Ptr<Vector<Real>>               getMultiplierVector(void);
  const Ptr<Vector<Real>>               getResidualVector(void);
  const Ptr<PolyhedralProjection<Real>> getPolyhedralProjection(void);
  EProblem                              getProblemType(void);

  /* Finalize and edit methods */
  void finalize(bool lumpConstraints = false);
  void edit(void);
  void finalizeIteration(void);

}; // class NewOptimizationProblem

}  // namespace ROL

#include "ROL_NewOptimizationProblem_Def.hpp"

#endif // ROL_NEWOPTIMIZATIONPROBLEM_HPP
