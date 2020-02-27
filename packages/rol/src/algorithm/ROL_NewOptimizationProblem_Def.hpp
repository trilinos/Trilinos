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

#ifndef ROL_NEWOPTIMIZATIONPROBLEM_DEF_HPP
#define ROL_NEWOPTIMIZATIONPROBLEM_DEF_HPP

namespace ROL {

template<typename Real>
NewOptimizationProblem<Real>::NewOptimizationProblem(const Ptr<Objective<Real>> &obj,
                                                     const Ptr<Vector<Real>>    &x,
                                                     const Ptr<Vector<Real>>    &g)
  : isFinalized_(false), hasBounds_(false),
    hasEquality_(false), hasInequality_(false),
    hasLinearEquality_(false), hasLinearInequality_(false),
    obj_(nullPtr), xprim_(nullPtr), xdual_(nullPtr), bnd_(nullPtr),
    con_(nullPtr), mul_(nullPtr), res_(nullPtr), proj_(nullPtr),
    problemType_(TYPE_U) {
  INPUT_obj_ = obj;
  INPUT_xprim_ = x;
  if (g==nullPtr) {
    INPUT_xdual_ = x->dual().clone();
  }
  else {
    INPUT_xdual_ = g;
  }
  INPUT_bnd_ = nullPtr;
  INPUT_econ_.clear();
  INPUT_icon_.clear();
  INPUT_linear_econ_.clear();
  INPUT_linear_icon_.clear();
}

template<typename Real>
void NewOptimizationProblem<Real>::addBoundConstraint(const Ptr<BoundConstraint<Real>> &bnd) {
  if (!isFinalized_) {
    INPUT_bnd_ = bnd;
    hasBounds_ = true;
  }
  else {
    throw Exception::NotImplemented(">>> ROL::NewOptimizationProblem: Cannot add bounds after problem is finalized!");
  }
}

template<typename Real>
void NewOptimizationProblem<Real>::removeBoundConstraint(void) {
  if (!isFinalized_) {
    INPUT_bnd_ = nullPtr;
    hasBounds_ = false;
  }
  else {
    throw Exception::NotImplemented(">>> ROL::NewOptimizationProblem: Cannot remove bounds after problem is finalized!");
  }
}

template<typename Real>
void NewOptimizationProblem<Real>::addEqualityConstraint(std::string                  name,
                                                         const Ptr<Constraint<Real>> &econ,
                                                         const Ptr<Vector<Real>>     &emul,
                                                         const Ptr<Vector<Real>>     &eres,
                                                         bool                         reset) {
  if (!isFinalized_) {
    if (!hasEquality_ || reset) {
      INPUT_econ_.clear();
    }
    INPUT_econ_.insert(std::pair<std::string,ConstraintData<Real>>(name,ConstraintData<Real>(econ,emul,eres)));
    hasEquality_ = true;
  }
  else {
    throw Exception::NotImplemented(">>> ROL::NewOptimizationProblem: Cannot add equality after problem is finalized!");
  }
}

template<typename Real>
void NewOptimizationProblem<Real>::removeEqualityConstraint(std::string name) {
  if (!isFinalized_) {
    auto it = INPUT_econ_.find(name);
    if (it!=INPUT_econ_.end()) {
      INPUT_econ_.erase(it);
    }
    if (INPUT_econ_.empty()) {
      hasEquality_ = false;
    }
  }
  else {
    throw Exception::NotImplemented(">>> ROL::NewOptimizationProblem: Cannot remove equality after problem is finalized!");
  }
}

template<typename Real>
void NewOptimizationProblem<Real>::addInequalityConstraint(std::string                       name,
                                                           const Ptr<Constraint<Real>>      &icon,
                                                           const Ptr<Vector<Real>>          &imul,
                                                           const Ptr<BoundConstraint<Real>> &ibnd,
                                                           const Ptr<Vector<Real>>          &ires,
                                                           bool                              reset) {
  if (!isFinalized_) {
    if (!hasInequality_ || reset) {
      INPUT_icon_.clear();
    }
    INPUT_icon_.insert(std::pair<std::string,ConstraintData<Real>>(name,ConstraintData<Real>(icon,imul,ires,ibnd)));
    hasInequality_ = true;
  }
  else {
    throw Exception::NotImplemented(">>> ROL::NewOptimizationProblem: Cannot add inequality after problem is finalized!");
  }
}

template<typename Real>
void NewOptimizationProblem<Real>::removeInequalityConstraint(std::string name) {
  if (!isFinalized_) {
    auto it = INPUT_icon_.find(name);
    if (it!=INPUT_icon_.end()) {
      INPUT_icon_.erase(it);
    }
    if (INPUT_icon_.empty()) {
      hasInequality_ = false;
    }
  }
  else {
    throw Exception::NotImplemented(">>> ROL::NewOptimizationProblem: Cannot remove inequality after problem is finalized!");
  }
}

template<typename Real>
void NewOptimizationProblem<Real>::addLinearEqualityConstraint(std::string                  name,
                                                               const Ptr<Constraint<Real>> &linear_econ,
                                                               const Ptr<Vector<Real>>     &linear_emul,
                                                               const Ptr<Vector<Real>>     &linear_eres,
                                                               bool                         reset) {
  if (!isFinalized_) {
    if (!hasLinearEquality_ || reset) {
      INPUT_linear_econ_.clear();
    }
    INPUT_linear_econ_.insert(std::pair<std::string,ConstraintData<Real>>(name,ConstraintData<Real>(linear_econ,linear_emul,linear_eres)));
    hasLinearEquality_ = true;
  }
  else {
    throw Exception::NotImplemented(">>> ROL::NewOptimizationProblem: Cannot add linear equality after problem is finalized!");
  }
}

template<typename Real>
void NewOptimizationProblem<Real>::removeLinearEqualityConstraint(std::string name) {
  if (!isFinalized_) {
    auto it = INPUT_linear_econ_.find(name);
    if (it!=INPUT_linear_econ_.end()) {
      INPUT_linear_econ_.erase(it);
    }
    if (INPUT_linear_econ_.empty()) {
      hasLinearEquality_ = false;
    }
  }
  else {
    throw Exception::NotImplemented(">>> ROL::NewOptimizationProblem: Cannot remove linear equality after problem is finalized!");
  }
}

template<typename Real>
void NewOptimizationProblem<Real>::addLinearInequalityConstraint(std::string                       name,
                                                                 const Ptr<Constraint<Real>>      &linear_icon,
                                                                 const Ptr<Vector<Real>>          &linear_imul,
                                                                 const Ptr<BoundConstraint<Real>> &linear_ibnd,
                                                                 const Ptr<Vector<Real>>          &linear_ires,
                                                                 bool                              reset) {
  if (!isFinalized_) {
    if (!hasLinearInequality_ || reset) {
      INPUT_linear_icon_.clear();
    }
    INPUT_linear_icon_.insert(std::pair<std::string,ConstraintData<Real>>(name,ConstraintData<Real>(linear_icon,linear_imul,linear_ires,linear_ibnd)));
    hasLinearInequality_ = true;
  }
  else {
    throw Exception::NotImplemented(">>> ROL::NewOptimizationProblem: Cannot add linear inequality after problem is finalized!");
  }
}

template<typename Real>
void NewOptimizationProblem<Real>::removeLinearInequalityConstraint(std::string name) {
  if (!isFinalized_) {
    auto it = INPUT_linear_icon_.find(name);
    if (it!=INPUT_linear_icon_.end()) {
      INPUT_linear_icon_.erase(it);
    }
    if (INPUT_linear_icon_.empty()) {
      hasLinearInequality_ = false;
    }
  }
  else {
    throw Exception::NotImplemented(">>> ROL::NewOptimizationProblem: Cannot remove linear inequality after problem is finalized!");
  }
}

template<typename Real>
void NewOptimizationProblem<Real>::finalize(bool lumpConstraints) {
  if (!isFinalized_) {
    int cnt = 0, lcnt = 0;
    std::map<std::string,ConstraintData<Real>> con, lcon, icon;
    bool hasEquality         = hasEquality_;
    bool hasLinearEquality   = hasLinearEquality_;
    bool hasInequality       = hasInequality_;
    bool hasLinearInequality = hasLinearInequality_;
    if (hasEquality_) {
      for (auto it = INPUT_econ_.begin(); it != INPUT_econ_.end(); ++it) { 
        con.insert(std::pair<std::string,ConstraintData<Real>>(std::to_string(cnt),it->second));
        cnt++;
      }
    }
    if (hasLinearEquality_) {
      for (auto it = INPUT_linear_econ_.begin(); it != INPUT_linear_econ_.end(); ++it) { 
        if (lumpConstraints) {
          con.insert(std::pair<std::string,ConstraintData<Real>>(std::to_string(cnt),it->second));
          cnt++;
        }
        else {
          lcon.insert(std::pair<std::string,ConstraintData<Real>>(std::to_string(lcnt),it->second));
          lcnt++;
        }
      }
      if (lumpConstraints) {
        hasEquality = true;
        hasLinearEquality = false;
      }
    }
    if (hasInequality_) {
      for (auto it = INPUT_icon_.begin(); it != INPUT_icon_.end(); ++it) { 
        con.insert(std::pair<std::string,ConstraintData<Real>>(std::to_string(cnt),it->second));
        cnt++;
      }
    }
    if (hasLinearInequality_) {
      for (auto it = INPUT_linear_icon_.begin(); it != INPUT_linear_icon_.end(); ++it) { 
        if (lumpConstraints) {
          con.insert(std::pair<std::string,ConstraintData<Real>>(std::to_string(cnt),it->second));
          cnt++;
        }
        else {
          lcon.insert(std::pair<std::string,ConstraintData<Real>>(std::to_string(lcnt),it->second));
          lcnt++;
        }
      }
      if (lumpConstraints) {
        hasInequality = true;
        hasLinearInequality = false;
      }
    }
    // Transform optimization problem
    //std::cout << hasBounds_ << "  " << hasEquality << "  " << hasInequality << "  " << hasLinearEquality << "  " << hasLinearInequality << std::endl;
    if (!hasLinearEquality && !hasLinearInequality) {
      proj_ = nullPtr;
      if (!hasEquality && !hasInequality && !hasBounds_) {
        problemType_ = TYPE_U;
        obj_         = INPUT_obj_;
        xprim_       = INPUT_xprim_;
        xdual_       = INPUT_xdual_;
        bnd_         = nullPtr;
        con_         = nullPtr;
        mul_         = nullPtr;
        res_         = nullPtr;
      }
      else if (!hasEquality && !hasInequality && hasBounds_) {
        problemType_ = TYPE_B;
        obj_         = INPUT_obj_;
        xprim_       = INPUT_xprim_;
        xdual_       = INPUT_xdual_;
        bnd_         = INPUT_bnd_;
        con_         = nullPtr;
        mul_         = nullPtr;
        res_         = nullPtr;
      }
      else if (hasEquality && !hasInequality && !hasBounds_) {
        NewConstraintManager<Real> cm(con,INPUT_xprim_,INPUT_xdual_);
        problemType_ = TYPE_E;
        obj_         = INPUT_obj_;
        xprim_       = INPUT_xprim_;
        xdual_       = INPUT_xdual_;
        bnd_         = nullPtr;
        con_         = cm.getConstraint();
        mul_         = cm.getMultiplier();
        res_         = cm.getResidual();
      }
      else {
        NewConstraintManager<Real> cm(con,INPUT_xprim_,INPUT_xdual_,INPUT_bnd_);
        problemType_ = TYPE_EB;
        obj_         = INPUT_obj_;
        if (cm.hasInequality()) {
          obj_      = makePtr<SlacklessObjective<Real>>(INPUT_obj_);
        }
        xprim_       = cm.getOptVector();
        xdual_       = cm.getDualOptVector();
        bnd_         = cm.getBoundConstraint();
        con_         = cm.getConstraint();
        mul_         = cm.getMultiplier();
        res_         = cm.getResidual();
      }
    }
    else {
      if (!hasBounds_ && !hasLinearInequality) {
        NewConstraintManager<Real> cm(lcon,INPUT_xprim_,INPUT_xdual_);
        xfeas_ = cm.getOptVector()->clone(); xfeas_->set(*cm.getOptVector());
        rlc_   = makePtr<ReduceLinearConstraint<Real>>(cm.getConstraint(),xfeas_,cm.getResidual());
        proj_  = nullPtr;
        if (!hasEquality && !hasInequality) {
          problemType_ = TYPE_U;
          obj_         = rlc_->transform(INPUT_obj_);
          xprim_       = xfeas_->clone(); xprim_->zero();
          xdual_       = cm.getDualOptVector();
          bnd_         = nullPtr;
          con_         = nullPtr;
          mul_         = nullPtr;
          res_         = nullPtr;
        }
        else {
          for (auto it = con.begin(); it != con.end(); ++it) {
            icon.insert(std::pair<std::string,ConstraintData<Real>>(it->first,
              ConstraintData<Real>(rlc_->transform(it->second.constraint),
                it->second.multiplier,it->second.residual,it->second.bounds)));
          }
          Ptr<Vector<Real>> xtmp = xfeas_->clone(); xtmp->zero();
          NewConstraintManager<Real> cm1(icon,xtmp,cm.getDualOptVector());
          xprim_         = cm1.getOptVector();
          xdual_         = cm1.getDualOptVector();
          con_           = cm1.getConstraint();
          mul_           = cm1.getMultiplier();
          res_           = cm1.getResidual();
          if (!hasInequality) {
            problemType_ = TYPE_E;
            obj_         = rlc_->transform(INPUT_obj_);
            bnd_         = nullPtr;
          }
          else {
            problemType_ = TYPE_EB;
            obj_         = makePtr<SlacklessObjective<Real>>(rlc_->transform(INPUT_obj_));
            bnd_         = cm1.getBoundConstraint();
          }
        }
      }
      else if ((hasBounds_ || hasLinearInequality) && !hasEquality && !hasInequality) {
        NewConstraintManager<Real> cm(lcon,INPUT_xprim_,INPUT_xdual_,INPUT_bnd_);
        problemType_ = TYPE_B;
        obj_         = INPUT_obj_;
        if (cm.hasInequality()) {
          obj_       = makePtr<SlacklessObjective<Real>>(INPUT_obj_);
        }
        xprim_       = cm.getOptVector();
        xdual_       = cm.getDualOptVector();
        bnd_         = cm.getBoundConstraint();
        con_         = nullPtr;
        mul_         = nullPtr;
        res_         = nullPtr;
        proj_        = makePtr<PolyhedralProjection<Real>>(*xprim_,*xdual_,bnd_,
                         cm.getConstraint(),*cm.getMultiplier(),*cm.getResidual());
      }
      else {
        NewConstraintManager<Real> cm(con,lcon,INPUT_xprim_,INPUT_xdual_,INPUT_bnd_);
        problemType_ = TYPE_EB;
        obj_         = INPUT_obj_;
        if (cm.hasInequality()) {
          obj_       = makePtr<SlacklessObjective<Real>>(INPUT_obj_);
        }
        xprim_       = cm.getOptVector();
        xdual_       = cm.getDualOptVector();
        con_         = cm.getConstraint();
        mul_         = cm.getMultiplier();
        res_         = cm.getResidual();
        bnd_         = cm.getBoundConstraint();
        proj_        = makePtr<PolyhedralProjection<Real>>(*xprim_,*xdual_,bnd_,
                        cm.getLinearConstraint(),*cm.getLinearMultiplier(),
                        *cm.getLinearResidual());
      }
    }
    isFinalized_ = true;
  }
}

template<typename Real>
const Ptr<Objective<Real>> NewOptimizationProblem<Real>::getObjective(void) {
  finalize();
  return obj_;
}

template<typename Real>
const Ptr<Vector<Real>> NewOptimizationProblem<Real>::getPrimalOptimizationVector(void) {
  finalize();
  return xprim_;
}

template<typename Real>
const Ptr<Vector<Real>> NewOptimizationProblem<Real>::getDualOptimizationVector(void) {
  finalize();
  return xprim_;
}

template<typename Real>
const Ptr<BoundConstraint<Real>> NewOptimizationProblem<Real>::getBoundConstraint(void) {
  finalize();
  return bnd_;
}

template<typename Real>
const Ptr<Constraint<Real>> NewOptimizationProblem<Real>::getConstraint(void) {
  finalize();
  return con_;
}

template<typename Real>
const Ptr<Vector<Real>> NewOptimizationProblem<Real>::getMultiplierVector(void) {
  finalize();
  return mul_;
}

template<typename Real>
const Ptr<Vector<Real>> NewOptimizationProblem<Real>::getResidualVector(void) {
  finalize();
  return res_;
}

template<typename Real>
const Ptr<PolyhedralProjection<Real>> NewOptimizationProblem<Real>::getPolyhedralProjection(void) {
  finalize();
  return proj_;
}

template<typename Real>
EProblem NewOptimizationProblem<Real>::getProblemType(void) {
  finalize();
  return problemType_;
}

template<typename Real>
void NewOptimizationProblem<Real>::edit(void) {
  isFinalized_ = false;
  rlc_  = nullPtr;
  proj_ = nullPtr;
}

template<typename Real>
void NewOptimizationProblem<Real>::finalizeIteration(void) {
  if (rlc_ != nullPtr) {
    if (!hasInequality_) {
      rlc_->project(*INPUT_xprim_,*xprim_);
      INPUT_xprim_->plus(*rlc_->getFeasibleVector());
    }
    else {
      Ptr<Vector<Real>> xprim = dynamic_cast<PartitionedVector<Real>&>(*xprim_).get(0)->clone();
      xprim->set(*dynamic_cast<PartitionedVector<Real>&>(*xprim_).get(0));
      rlc_->project(*INPUT_xprim_,*xprim);
      INPUT_xprim_->plus(*rlc_->getFeasibleVector());
    }
  }
}

}  // namespace ROL

#endif // ROL_NEWOPTIMIZATIONPROBLEM_DEF_HPP
