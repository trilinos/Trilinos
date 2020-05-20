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

#ifndef ROL_OPTIMIZATIONPROBLEM_PEBBL_H
#define ROL_OPTIMIZATIONPROBLEM_PEBBL_H

#include "ROL_Ptr.hpp"
#include "ROL_NewOptimizationProblem.hpp"
#include "ROL_Transform_PEBBL.hpp"
#include "ROL_TransformedObjective_PEBBL.hpp"
#include "ROL_TransformedConstraint_PEBBL.hpp"


namespace ROL {

template <class Real>
class OptimizationProblem_PEBBL : public NewOptimizationProblem<Real> {
private:
  Ptr<Transform_PEBBL<Real>> trans_;

  Ptr<Objective<Real>>                                 ORIGINAL_obj_;
  std::unordered_map<std::string,ConstraintData<Real>> ORIGINAL_con_;
  std::unordered_map<std::string,ConstraintData<Real>> ORIGINAL_linear_con_;

  using NewOptimizationProblem<Real>::INPUT_xprim_;
  using NewOptimizationProblem<Real>::INPUT_obj_;
  using NewOptimizationProblem<Real>::INPUT_con_;
  using NewOptimizationProblem<Real>::INPUT_linear_con_;
  using NewOptimizationProblem<Real>::isFinalized;

public:
  /** \brief Default constructor for StochasticProblem.

      @param[in] obj  objective function object
      @param[in] x    primal optimization space vector
      @param[in] g    dual optimization space vector
  */
  OptimizationProblem_PEBBL(const Ptr<Objective<Real>> &obj,
                            const Ptr<Vector<Real>>    &x,
                            const Ptr<Vector<Real>>    &g = nullPtr)
    : NewOptimizationProblem<Real>(obj,x,g) {}

  void setTransformation(const Ptr<Transform_PEBBL<Real>> trans) {
    // Throw an exception if problem has been finalized
    ROL_TEST_FOR_EXCEPTION(isFinalized(),std::invalid_argument,
      ">>> ROL::OptimizationProblem_PEBBL::setTransformation: Cannot set transformation after problem has been finalized!");
    trans_ = trans;
  }

  void resetTransformation(void) {
    // Throw an exception if problem has been finalized
    ROL_TEST_FOR_EXCEPTION(isFinalized(),std::invalid_argument,
      ">>> ROL::OptimizationProblem_PEBBL::resetTransformation: Cannot reset transformation after problem has been finalized!");
    ROL_TEST_FOR_EXCEPTION(trans_==nullPtr,std::invalid_argument,
      ">>> ROL::OptimizationProblem_PEBBL::resetTransformation: Transformation has not been set!");
    trans_->reset();
  }

  void updateTransformation(const std::map<int,Real> &in) {
    // Throw an exception if problem has been finalized
    ROL_TEST_FOR_EXCEPTION(isFinalized(),std::invalid_argument,
      ">>> ROL::OptimizationProblem_PEBBL::resetTransformation: Cannot reset transformation after problem has been finalized!");
    // Throw an exception if problem has been finalized
    ROL_TEST_FOR_EXCEPTION(trans_==nullPtr,std::invalid_argument,
      ">>> ROL::OptimizationProblem_PEBBL::updateTransformation: Transformation has not been set!");
    trans_->add(in);
  }

  void updateTransformation(const std::pair<int,Real> &in) {
    // Throw an exception if problem has been finalized
    ROL_TEST_FOR_EXCEPTION(isFinalized(),std::invalid_argument,
      ">>> ROL::OptimizationProblem_PEBBL::resetTransformation: Cannot reset transformation after problem has been finalized!");
    // Throw an exception if problem has been finalized
    ROL_TEST_FOR_EXCEPTION(trans_==nullPtr,std::invalid_argument,
      ">>> ROL::OptimizationProblem_PEBBL::updateTransformation: Transformation has not been set!");
    trans_->add(in);
  }

  void deleteTransformation(void) {
    // Throw an exception if problem has been finalized
    ROL_TEST_FOR_EXCEPTION(isFinalized(),std::invalid_argument,
      ">>> ROL::OptimizationProblem_PEBBL::deleteTransformation: Cannot delete transformation after problem has been finalized!");
    trans_ = nullPtr;
  }

  /***************************************************************************/
  /*** Finalize and edit methods *********************************************/
  /***************************************************************************/

  void finalize(bool lumpConstraints = false, bool printToStream = false,
                std::ostream &outStream = std::cout) override {
    if (!isFinalized()) {
      if (trans_ != nullPtr) {
        ORIGINAL_obj_ = INPUT_obj_;
        INPUT_obj_ = makePtr<TransformedObjective_PEBBL<Real>>(ORIGINAL_obj_, trans_);
        // Apply transformation to optimization vector
        Real tol(std::sqrt(ROL_EPSILON<Real>()));
        trans_->value(*INPUT_xprim_,*INPUT_xprim_,tol);
        // Transform nonlinear constraints
        ORIGINAL_con_.clear();
        ORIGINAL_con_.insert(INPUT_con_.begin(),INPUT_con_.end());
        for (auto it = ORIGINAL_con_.begin(); it != ORIGINAL_con_.end(); ++it) {
          Ptr<Constraint<Real>>      con = makePtr<TransformedConstraint_PEBBL<Real>>(it->second.constraint, trans_);
          Ptr<Vector<Real>>          mul = it->second.multiplier;
          Ptr<Vector<Real>>          res = it->second.residual;
          Ptr<BoundConstraint<Real>> bnd = it->second.bounds;
          NewOptimizationProblem<Real>::removeConstraint(it->first);
          if (bnd != nullPtr) NewOptimizationProblem<Real>::addConstraint(it->first,con,mul,bnd,res);
          else                NewOptimizationProblem<Real>::addConstraint(it->first,con,mul,res);
        }
        // Transform linear constraints
        ORIGINAL_linear_con_.clear();
        ORIGINAL_linear_con_.insert(INPUT_linear_con_.begin(),INPUT_linear_con_.end());
        for (auto it = ORIGINAL_linear_con_.begin(); it != ORIGINAL_linear_con_.end(); ++it) {
          Ptr<Constraint<Real>>      con = makePtr<TransformedConstraint_PEBBL<Real>>(it->second.constraint, trans_);
          Ptr<Vector<Real>>          mul = it->second.multiplier;
          Ptr<Vector<Real>>          res = it->second.residual;
          Ptr<BoundConstraint<Real>> bnd = it->second.bounds;
          NewOptimizationProblem<Real>::removeLinearConstraint(it->first);
          if (bnd != nullPtr) NewOptimizationProblem<Real>::addLinearConstraint(it->first,con,mul,bnd,res);
          else                NewOptimizationProblem<Real>::addLinearConstraint(it->first,con,mul,res);
        }
      }

      // Call default finalize
      NewOptimizationProblem<Real>::finalize(lumpConstraints,printToStream,outStream);
    }
  }

  void edit(void) override {
    NewOptimizationProblem<Real>::edit();
    if (trans_ != nullPtr) {
      INPUT_obj_ = ORIGINAL_obj_;
      ORIGINAL_obj_ = nullPtr;
      // Transform nonlinear constraints
      for (auto it = ORIGINAL_con_.begin(); it != ORIGINAL_con_.end(); ++it) {
        Ptr<Constraint<Real>>      con = it->second.constraint;
        Ptr<Vector<Real>>          mul = it->second.multiplier;
        Ptr<Vector<Real>>          res = it->second.residual;
        Ptr<BoundConstraint<Real>> bnd = it->second.bounds;
        NewOptimizationProblem<Real>::removeConstraint(it->first);
        if (bnd != nullPtr) NewOptimizationProblem<Real>::addConstraint(it->first,con,mul,bnd,res);
        else                NewOptimizationProblem<Real>::addConstraint(it->first,con,mul,res);
      }
      // Transform linear constraints
      for (auto it = ORIGINAL_linear_con_.begin(); it != ORIGINAL_linear_con_.end(); ++it) {
        Ptr<Constraint<Real>>      con = it->second.constraint;
        Ptr<Vector<Real>>          mul = it->second.multiplier;
        Ptr<Vector<Real>>          res = it->second.residual;
        Ptr<BoundConstraint<Real>> bnd = it->second.bounds;
        NewOptimizationProblem<Real>::removeLinearConstraint(it->first);
        if (bnd != nullPtr) NewOptimizationProblem<Real>::addLinearConstraint(it->first,con,mul,bnd,res);
        else                NewOptimizationProblem<Real>::addLinearConstraint(it->first,con,mul,res);
      }
    }
  }

}; // class OptimizationProblemFactory

} // namespace ROL

#endif
