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

#ifndef ROL_CONSTRAINT_ASSEMBLER_H
#define ROL_CONSTRAINT_ASSEMBLER_H

#include "ROL_Constraint_Partitioned.hpp"
#include "ROL_BoundConstraint_Partitioned.hpp"
#include <unordered_map>

/** @ingroup func_group
    \class ROL::ConstraintAssembler
    \brief Provides a wrapper for multiple constraints.

    ---
*/

namespace ROL {

template<typename Real>
struct ConstraintData {
  const Ptr<Constraint<Real>>      constraint;
  const Ptr<Vector<Real>>          multiplier;
  const Ptr<Vector<Real>>          residual;
  const Ptr<BoundConstraint<Real>> bounds;

  ConstraintData(const Ptr<Constraint<Real>>      &con,
                 const Ptr<Vector<Real>>          &mul,
                 const Ptr<Vector<Real>>          &res=nullPtr,
                 const Ptr<BoundConstraint<Real>> &bnd=nullPtr)
    : constraint(con), multiplier(mul),
      residual(res==nullPtr ? mul->dual().clone() : res), bounds(bnd) {}
};

template<typename Real>
class ConstraintAssembler {
private:
  Ptr<Constraint<Real>>      con_;
  Ptr<Vector<Real>>          mul_;
  Ptr<Vector<Real>>          res_;
  Ptr<Constraint<Real>>      linear_con_;
  Ptr<Vector<Real>>          linear_mul_;
  Ptr<Vector<Real>>          linear_res_;
  Ptr<Vector<Real>>          xprim_;
  Ptr<Vector<Real>>          xdual_;
  Ptr<BoundConstraint<Real>> bnd_;

  std::vector<Ptr<Constraint<Real>>>      cvec_;  // General constraints
  std::vector<Ptr<Vector<Real>>>          lvec_;  // General multiplier vector
  std::vector<Ptr<Vector<Real>>>          rvec_;  // General residual vector
  std::vector<Ptr<Constraint<Real>>>      lcvec_; // Linear constraints
  std::vector<Ptr<Vector<Real>>>          llvec_; // Linear multiplier vector
  std::vector<Ptr<Vector<Real>>>          lrvec_; // Linear residual vector
  std::vector<Ptr<Vector<Real>>>          psvec_; // Primal slack vector
  std::vector<Ptr<Vector<Real>>>          dsvec_; // Dual slack vector
  std::vector<Ptr<BoundConstraint<Real>>> sbnd_;  // Slack bound constraint

  std::vector<bool> isInequality_, isLinearInequality_;

  bool isNull_;
  bool hasInequality_;

  void initializeSlackVariable( const Ptr<Constraint<Real>>      &con,
                                const Ptr<BoundConstraint<Real>> &cbnd,
                                const Ptr<Vector<Real>>          &s,
                                const Ptr<Vector<Real>>          &x) const;

  void initialize( const std::unordered_map<std::string,ConstraintData<Real>> &input_con,
                   const Ptr<Vector<Real>>                                    &xprim,
                   const Ptr<Vector<Real>>                                    &xdual,
                   const Ptr<BoundConstraint<Real>>                           &bnd);

  void initialize( const std::unordered_map<std::string,ConstraintData<Real>> &input_con,
                   const std::unordered_map<std::string,ConstraintData<Real>> &input_lcon,
                   const Ptr<Vector<Real>>                                    &xprim,
                   const Ptr<Vector<Real>>                                    &xdual,
                   const Ptr<BoundConstraint<Real>>                           &bnd);

public:
  virtual ~ConstraintAssembler() {}

  ConstraintAssembler( const std::unordered_map<std::string,ConstraintData<Real>> &con,
                     const Ptr<Vector<Real>>                                    &xprim,
                     const Ptr<Vector<Real>>                                    &xdual,
                     const Ptr<BoundConstraint<Real>>                           &bnd = nullPtr);

  ConstraintAssembler( const std::unordered_map<std::string,ConstraintData<Real>> &con,
                     const std::unordered_map<std::string,ConstraintData<Real>> &linear_con,
                     const Ptr<Vector<Real>>                                    &xprim,
                     const Ptr<Vector<Real>>                                    &xdual,
                     const Ptr<BoundConstraint<Real>>                           &bnd = nullPtr);

  const Ptr<Constraint<Real>>&      getConstraint() const;
  const Ptr<Vector<Real>>&          getMultiplier() const;
  const Ptr<Vector<Real>>&          getResidual() const;
  const Ptr<Constraint<Real>>&      getLinearConstraint() const;
  const Ptr<Vector<Real>>&          getLinearMultiplier() const;
  const Ptr<Vector<Real>>&          getLinearResidual() const;
  const Ptr<Vector<Real>>&          getOptVector() const;
  const Ptr<Vector<Real>>&          getDualOptVector() const;
  const Ptr<BoundConstraint<Real>>& getBoundConstraint() const;

  bool isNull() const;
  bool hasInequality() const;

  void resetSlackVariables();

}; // class ConstraintAssembler

} // namespace ROL

#include "ROL_ConstraintAssembler_Def.hpp"

#endif
