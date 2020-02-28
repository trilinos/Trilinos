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
    cnt_econ_(0), cnt_icon_(0), cnt_linear_econ_(0), cnt_linear_icon_(0),
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
  INPUT_con_.clear();
  INPUT_linear_con_.clear();
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
void NewOptimizationProblem<Real>::addConstraint(std::string                  name,
                                                 const Ptr<Constraint<Real>> &econ,
                                                 const Ptr<Vector<Real>>     &emul,
                                                 const Ptr<Vector<Real>>     &eres,
                                                 bool                         reset) {
  if (!isFinalized_) {
    if (reset) {
      INPUT_con_.clear();
    }
    INPUT_con_.insert(std::pair<std::string,ConstraintData<Real>>(name,ConstraintData<Real>(econ,emul,eres)));
    hasEquality_ = true;
    cnt_econ_++;
  }
  else {
    throw Exception::NotImplemented(">>> ROL::NewOptimizationProblem: Cannot add constraint after problem is finalized!");
  }
}

template<typename Real>
void NewOptimizationProblem<Real>::addConstraint(std::string                       name,
                                                 const Ptr<Constraint<Real>>      &icon,
                                                 const Ptr<Vector<Real>>          &imul,
                                                 const Ptr<BoundConstraint<Real>> &ibnd,
                                                 const Ptr<Vector<Real>>          &ires,
                                                 bool                              reset) {
  if (!isFinalized_) {
    if (reset) {
      INPUT_con_.clear();
    }
    INPUT_con_.insert(std::pair<std::string,ConstraintData<Real>>(name,ConstraintData<Real>(icon,imul,ires,ibnd)));
    hasInequality_ = true;
    cnt_icon_++;
  }
  else {
    throw Exception::NotImplemented(">>> ROL::NewOptimizationProblem: Cannot add constraint after problem is finalized!");
  }
}

template<typename Real>
void NewOptimizationProblem<Real>::removeConstraint(std::string name) {
  if (!isFinalized_) {
    auto it = INPUT_con_.find(name);
    if (it!=INPUT_con_.end()) {
      if (it->second.bounds==nullPtr) {
        cnt_econ_--;
      }
      else {
        cnt_icon_--;
      }
      INPUT_con_.erase(it);
    }
    if (cnt_econ_==0) {
      hasEquality_ = false;
    }
    if (cnt_icon_==0) {
      hasInequality_ = false;
    }
  }
  else {
    throw Exception::NotImplemented(">>> ROL::NewOptimizationProblem: Cannot remove constraint after problem is finalized!");
  }
}

template<typename Real>
void NewOptimizationProblem<Real>::addLinearConstraint(std::string                  name,
                                                       const Ptr<Constraint<Real>> &linear_econ,
                                                       const Ptr<Vector<Real>>     &linear_emul,
                                                       const Ptr<Vector<Real>>     &linear_eres,
                                                       bool                         reset) {
  if (!isFinalized_) {
    if (reset) {
      INPUT_linear_con_.clear();
    }
    INPUT_linear_con_.insert(std::pair<std::string,ConstraintData<Real>>(name,ConstraintData<Real>(linear_econ,linear_emul,linear_eres)));
    hasLinearEquality_ = true;
    cnt_linear_econ_++;
  }
  else {
    throw Exception::NotImplemented(">>> ROL::NewOptimizationProblem: Cannot add linear constraint after problem is finalized!");
  }
}

template<typename Real>
void NewOptimizationProblem<Real>::addLinearConstraint(std::string                       name,
                                                       const Ptr<Constraint<Real>>      &linear_icon,
                                                       const Ptr<Vector<Real>>          &linear_imul,
                                                       const Ptr<BoundConstraint<Real>> &linear_ibnd,
                                                       const Ptr<Vector<Real>>          &linear_ires,
                                                       bool                              reset) {
  if (!isFinalized_) {
    if (reset) {
      INPUT_linear_con_.clear();
    }
    INPUT_linear_con_.insert(std::pair<std::string,ConstraintData<Real>>(name,ConstraintData<Real>(linear_icon,linear_imul,linear_ires,linear_ibnd)));
    hasLinearInequality_ = true;
    cnt_linear_icon_++;
  }
  else {
    throw Exception::NotImplemented(">>> ROL::NewOptimizationProblem: Cannot add linear constraint after problem is finalized!");
  }
}

template<typename Real>
void NewOptimizationProblem<Real>::removeLinearConstraint(std::string name) {
  if (!isFinalized_) {
    auto it = INPUT_linear_con_.find(name);
    if (it!=INPUT_linear_con_.end()) {
      if (it->second.bounds==nullPtr) {
        cnt_linear_econ_--;
      }
      else {
        cnt_linear_icon_--;
      }
      INPUT_linear_con_.erase(it);
    }
    if (cnt_linear_econ_==0) {
      hasLinearEquality_ = false;
    }
    if (cnt_linear_icon_==0) {
      hasLinearInequality_ = false;
    }
  }
  else {
    throw Exception::NotImplemented(">>> ROL::NewOptimizationProblem: Cannot remove linear inequality after problem is finalized!");
  }
}

template<typename Real>
void NewOptimizationProblem<Real>::finalize(bool lumpConstraints, bool printToStream, std::ostream &outStream) {
  if (!isFinalized_) {
    std::map<std::string,ConstraintData<Real>> con, lcon, icon;
    bool hasEquality         = hasEquality_;
    bool hasLinearEquality   = hasLinearEquality_;
    bool hasInequality       = hasInequality_;
    bool hasLinearInequality = hasLinearInequality_;
    con.insert(INPUT_con_.begin(),INPUT_con_.end());
    if (lumpConstraints) {
      con.insert(INPUT_linear_con_.begin(),INPUT_linear_con_.end());
      hasEquality = (hasEquality || hasLinearEquality);
      hasInequality = (hasInequality || hasLinearInequality);
      hasLinearEquality = false;
      hasLinearInequality = false;
    }
    else {
      lcon.insert(INPUT_linear_con_.begin(),INPUT_linear_con_.end());
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
    if (printToStream) {
      //std::ios_base_fmtflags state(outStream.flags());
      outStream << std::endl;
      outStream << "  ROL::NewOptimizationProblem::finalize" << std::endl;
      outStream << "    Problem Summary:" << std::endl;
      outStream << "      Has Bound Constraint? .............. " << (hasBounds_ ? "yes" : "no") << std::endl;
      outStream << "      Has Equality Constraint? ........... " << (hasEquality ? "yes" : "no") << std::endl; 
      if (hasEquality) {
        int cnt = 0;
	for (auto it = con.begin(); it != con.end(); ++it) {
          if (it->second.bounds==nullPtr) {
            if (cnt==0) {
              outStream << "        Names: ........................... ";
	      cnt++;
	    }
	    else {
              outStream << "                                           ";
	    }
            outStream << it->first << std::endl;
	  }
	}
        outStream << "        Total: ........................... " << cnt_econ_+(lumpConstraints ? cnt_linear_econ_ : 0) << std::endl; 
      }
      outStream << "      Has Inequality Constraint? ......... " << (hasInequality ? "yes" : "no") << std::endl; 
      if (hasInequality) {
        int cnt = 0;
	for (auto it = con.begin(); it != con.end(); ++it) {
          if (it->second.bounds!=nullPtr) {
            if (cnt==0) {
              outStream << "        Names: ........................... ";
	      cnt++;
	    }
	    else {
              outStream << "                                           ";
	    }
            outStream << it->first << std::endl;
	  }
	}
        outStream << "        Total: ........................... " << cnt_icon_+(lumpConstraints ? cnt_linear_icon_ : 0) << std::endl; 
      }
      if (!lumpConstraints) {
        outStream << "      Has Linear Equality Constraint? .... " << (hasLinearEquality ? "yes" : "no") << std::endl;
        if (hasLinearEquality) {
          int cnt = 0;
	  for (auto it = lcon.begin(); it != lcon.end(); ++it) {
            if (it->second.bounds==nullPtr) {
              if (cnt==0) {
                outStream << "        Names: ........................... ";
		cnt++;
	      }
	      else {
                outStream << "                                           ";
	      }
              outStream << it->first << std::endl;
	    }
	  }
          outStream << "        Total: ........................... " << cnt_linear_econ_ << std::endl; 
        }
        outStream << "      Has Linear Inequality Constraint? .. " << (hasLinearInequality ? "yes" : "no") << std::endl;
        if (hasLinearInequality) {
          int cnt = 0;
	  for (auto it = lcon.begin(); it != lcon.end(); ++it) {
            if (it->second.bounds!=nullPtr) {
              if (cnt==0) {
                outStream << "        Names: ........................... ";
		cnt++;
	      }
	      else {
                outStream << "                                           ";
	      }
              outStream << it->first << std::endl;
	    }
	  }
          outStream << "        Total: ........................... " << cnt_linear_icon_ << std::endl; 
        }
      }
      outStream << std::endl;
      //outStream.flags(state);
    }
  }
  else {
    if (printToStream) {
      outStream << std::endl;
      outStream << "  ROL::NewOptimizationProblem::finalize" << std::endl;
      outStream << "    Problem already finalized!" << std::endl;
      outStream << std::endl;
    }
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
