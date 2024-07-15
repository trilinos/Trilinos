// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PROBLEM_DEF_HPP
#define ROL_PROBLEM_DEF_HPP

#include <iostream>

namespace ROL {

template<typename Real>
Problem<Real>::Problem( const Ptr<Objective<Real>> &obj,
                        const Ptr<Vector<Real>>    &x,
                        const Ptr<Vector<Real>>    &g)
  : isFinalized_(false), hasBounds_(false),
    hasEquality_(false), hasInequality_(false),
    hasLinearEquality_(false), hasLinearInequality_(false),
    cnt_econ_(0), cnt_icon_(0), cnt_linear_econ_(0), cnt_linear_icon_(0),
    obj_(nullPtr), xprim_(nullPtr), xdual_(nullPtr), bnd_(nullPtr),
    con_(nullPtr), mul_(nullPtr), res_(nullPtr), proj_(nullPtr),
    problemType_(TYPE_U) {
  INPUT_obj_   = obj;
  INPUT_xprim_ = x;
  INPUT_bnd_   = nullPtr;
  INPUT_con_.clear();
  INPUT_linear_con_.clear();
  if (g==nullPtr) INPUT_xdual_ = x->dual().clone();
  else            INPUT_xdual_ = g;
}

template<typename Real>
void Problem<Real>::addBoundConstraint(const Ptr<BoundConstraint<Real>> &bnd) {
  ROL_TEST_FOR_EXCEPTION(isFinalized_,std::invalid_argument,
    ">>> ROL::Problem: Cannot add bounds after problem is finalized!");

  INPUT_bnd_ = bnd;
  hasBounds_ = true;
}

template<typename Real>
void Problem<Real>::removeBoundConstraint() {
  ROL_TEST_FOR_EXCEPTION(isFinalized_,std::invalid_argument,
    ">>> ROL::Problem: Cannot remove bounds after problem is finalized!");

  INPUT_bnd_ = nullPtr;
  hasBounds_ = false;
}

template<typename Real>
void Problem<Real>::addProxObjective(const Ptr<ProxObjective<Real>> &prox) {
  ROL_TEST_FOR_EXCEPTION(isFinalized_,std::invalid_argument,
    ">>> ROL::Problem: Cannot add prox objective after problem is finalized!");

  INPUT_prox_ = prox;
  hasProx_    = true;
}

template<typename Real>
void Problem<Real>::removeProxObjective() {
  ROL_TEST_FOR_EXCEPTION(isFinalized_,std::invalid_argument,
    ">>> ROL::Problem: Cannot remove prox objective after problem is finalized!");

  INPUT_prox_ = nullPtr;
  hasProx_    = false;
}

template<typename Real>
void Problem<Real>::addConstraint( std::string                  name,
                                   const Ptr<Constraint<Real>> &econ,
                                   const Ptr<Vector<Real>>     &emul,
                                   const Ptr<Vector<Real>>     &eres,
                                   bool                         reset) {
  ROL_TEST_FOR_EXCEPTION(isFinalized_,std::invalid_argument,
    ">>> ROL::Problem: Cannot add constraint after problem is finalized!");

  if (reset) INPUT_con_.clear();

  auto it = INPUT_con_.find(name);
  ROL_TEST_FOR_EXCEPTION(it != INPUT_con_.end(),std::invalid_argument,
    ">>> ROL::Problem: Constraint names must be distinct!");

  INPUT_con_.insert({name,ConstraintData<Real>(econ,emul,eres)});
  hasEquality_ = true;
  cnt_econ_++;
}

template<typename Real>
void Problem<Real>::addConstraint( std::string                       name,
                                   const Ptr<Constraint<Real>>      &icon,
                                   const Ptr<Vector<Real>>          &imul,
                                   const Ptr<BoundConstraint<Real>> &ibnd,
                                   const Ptr<Vector<Real>>          &ires,
                                   bool                              reset) {
  ROL_TEST_FOR_EXCEPTION(isFinalized_,std::invalid_argument,
    ">>> ROL::Problem: Cannot add constraint after problem is finalized!");

  if (reset) INPUT_con_.clear();

  auto it = INPUT_con_.find(name);
  ROL_TEST_FOR_EXCEPTION(it != INPUT_con_.end(),std::invalid_argument,
    ">>> ROL::Problem: Constraint names must be distinct!");

  INPUT_con_.insert({name,ConstraintData<Real>(icon,imul,ires,ibnd)});
  hasInequality_ = true;
  cnt_icon_++;
}

template<typename Real>
void Problem<Real>::removeConstraint(std::string name) {
  ROL_TEST_FOR_EXCEPTION(isFinalized_,std::invalid_argument,
    ">>> ROL::Problem: Cannot remove constraint after problem is finalized!");

  auto it = INPUT_con_.find(name);
  if (it!=INPUT_con_.end()) {
    if (it->second.bounds==nullPtr) cnt_econ_--;
    else                            cnt_icon_--;
    INPUT_con_.erase(it);
  }
  if (cnt_econ_==0) hasEquality_   = false;
  if (cnt_icon_==0) hasInequality_ = false;
}

template<typename Real>
void Problem<Real>::addLinearConstraint( std::string                  name,
                                         const Ptr<Constraint<Real>> &linear_econ,
                                         const Ptr<Vector<Real>>     &linear_emul,
                                         const Ptr<Vector<Real>>     &linear_eres,
                                         bool                         reset) {
  ROL_TEST_FOR_EXCEPTION(isFinalized_,std::invalid_argument,
    ">>> ROL::Problem: Cannot add linear constraint after problem is finalized!");

  if (reset) INPUT_linear_con_.clear();

  auto it = INPUT_linear_con_.find(name);
  ROL_TEST_FOR_EXCEPTION(it != INPUT_linear_con_.end(),std::invalid_argument,
    ">>> ROL::Problem: Linear constraint names must be distinct!");

  INPUT_linear_con_.insert({name,ConstraintData<Real>(linear_econ,linear_emul,linear_eres)});
  hasLinearEquality_ = true;
  cnt_linear_econ_++;
}

template<typename Real>
void Problem<Real>::addLinearConstraint( std::string                       name,
                                         const Ptr<Constraint<Real>>      &linear_icon,
                                         const Ptr<Vector<Real>>          &linear_imul,
                                         const Ptr<BoundConstraint<Real>> &linear_ibnd,
                                         const Ptr<Vector<Real>>          &linear_ires,
                                         bool                              reset) {
  ROL_TEST_FOR_EXCEPTION(isFinalized_,std::invalid_argument,
    ">>> ROL::Problem: Cannot add linear constraint after problem is finalized!");

  if (reset) INPUT_linear_con_.clear();

  auto it = INPUT_linear_con_.find(name);
  ROL_TEST_FOR_EXCEPTION(it != INPUT_linear_con_.end(),std::invalid_argument,
    ">>> ROL::Problem: Linear constraint names must be distinct!");

  INPUT_linear_con_.insert({name,ConstraintData<Real>(linear_icon,linear_imul,linear_ires,linear_ibnd)});
  hasLinearInequality_ = true;
  cnt_linear_icon_++;
}

template<typename Real>
void Problem<Real>::removeLinearConstraint(std::string name) {
  ROL_TEST_FOR_EXCEPTION(isFinalized_,std::invalid_argument,
    ">>> ROL::Problem: Cannot remove linear inequality after problem is finalized!");

  auto it = INPUT_linear_con_.find(name);
  if (it!=INPUT_linear_con_.end()) {
    if (it->second.bounds==nullPtr) cnt_linear_econ_--;
    else                            cnt_linear_icon_--;
    INPUT_linear_con_.erase(it);
  }
  if (cnt_linear_econ_==0) hasLinearEquality_ = false;
  if (cnt_linear_icon_==0) hasLinearInequality_ = false;
}

template<typename Real>
void Problem<Real>::setProjectionAlgorithm(ParameterList &list) {
  ROL_TEST_FOR_EXCEPTION(isFinalized_,std::invalid_argument,
    ">>> ROL::Problem: Cannot set polyhedral projection algorithm after problem is finalized!");

  ppa_list_ = list;
}

template<typename Real>
void Problem<Real>::finalize(bool lumpConstraints, bool printToStream, std::ostream &outStream) {
  if (!isFinalized_) {
    std::unordered_map<std::string,ConstraintData<Real>> con, lcon, icon;
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
    bool proxCompatible = false;
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
        if (hasProx_) { prox_ = INPUT_prox_; }
        else          { prox_ = nullPtr;     }
        proxCompatible = true;
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
        ConstraintAssembler<Real> cm(con,INPUT_xprim_,INPUT_xdual_);
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
        ConstraintAssembler<Real> cm(con,INPUT_xprim_,INPUT_xdual_,INPUT_bnd_);
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
        ConstraintAssembler<Real> cm(lcon,INPUT_xprim_,INPUT_xdual_);
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
          ConstraintAssembler<Real> cm1(icon,xtmp,cm.getDualOptVector());
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
        ConstraintAssembler<Real> cm(lcon,INPUT_xprim_,INPUT_xdual_,INPUT_bnd_);
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
        proj_        = PolyhedralProjectionFactory<Real>(*xprim_,*xdual_,bnd_,
                         cm.getConstraint(),*cm.getMultiplier(),*cm.getResidual(),ppa_list_);
      }
      else {
        ConstraintAssembler<Real> cm(con,lcon,INPUT_xprim_,INPUT_xdual_,INPUT_bnd_);
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
        proj_        = PolyhedralProjectionFactory<Real>(*xprim_,*xdual_,bnd_,
                        cm.getLinearConstraint(),*cm.getLinearMultiplier(),
                        *cm.getLinearResidual(),ppa_list_);
      }
    }

    std::stringstream out;
    out << ">>> ROL::Problem: Cannot solve ";
    if (problemType_ == TYPE_B)       out << "TypeB";
    else if (problemType_ == TYPE_E)  out << "TypeE";
    else if (problemType_ == TYPE_EB) out << "TypeG";
    else                              out << "TypeU with linear constraints";
    out << " problems with a prox objective!";
    ROL_TEST_FOR_EXCEPTION(hasProx_&&!proxCompatible,std::invalid_argument,out.str());

    isFinalized_ = true;
    if (printToStream) {
      outStream << std::endl;
      outStream << "  ROL::Problem::finalize" << std::endl;
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
    }
  }
  else {
    if (printToStream) {
      outStream << std::endl;
      outStream << "  ROL::Problem::finalize" << std::endl;
      outStream << "    Problem already finalized!" << std::endl;
      outStream << std::endl;
    }
  }
}

template<typename Real>
const Ptr<Objective<Real>>& Problem<Real>::getObjective() {
  finalize();
  return obj_;
}

template<typename Real>
const Ptr<Vector<Real>>& Problem<Real>::getPrimalOptimizationVector() {
  finalize();
  return xprim_;
}

template<typename Real>
const Ptr<Vector<Real>>& Problem<Real>::getDualOptimizationVector() {
  finalize();
  return xdual_;
}

template<typename Real>
const Ptr<BoundConstraint<Real>>& Problem<Real>::getBoundConstraint() {
  finalize();
  return bnd_;
}

template<typename Real>
const Ptr<Constraint<Real>>& Problem<Real>::getConstraint() {
  finalize();
  return con_;
}

template<typename Real>
const Ptr<Vector<Real>>& Problem<Real>::getMultiplierVector() {
  finalize();
  return mul_;
}

template<typename Real>
const Ptr<Vector<Real>>& Problem<Real>::getResidualVector() {
  finalize();
  return res_;
}

template<typename Real>
const Ptr<PolyhedralProjection<Real>>& Problem<Real>::getPolyhedralProjection() {
  finalize();
  return proj_;
}

template<typename Real>
EProblem Problem<Real>::getProblemType() {
  finalize();
  return problemType_;
}

template<typename Real>
Real Problem<Real>::checkLinearity(bool printToStream, std::ostream &outStream) const {
  std::ios_base::fmtflags state(outStream.flags());
  if (printToStream) {
    outStream << std::setprecision(8) << std::scientific;
    outStream << std::endl;
    outStream << "  ROL::Problem::checkLinearity" << std::endl;
  }
  const Real one(1), two(2), eps(1e-2*std::sqrt(ROL_EPSILON<Real>()));
  Real tol(std::sqrt(ROL_EPSILON<Real>())), cnorm(0), err(0), maxerr(0);
  Ptr<Vector<Real>> x  = INPUT_xprim_->clone(); x->randomize(-one,one);
  Ptr<Vector<Real>> y  = INPUT_xprim_->clone(); y->randomize(-one,one);
  Ptr<Vector<Real>> z  = INPUT_xprim_->clone(); z->zero();
  Ptr<Vector<Real>> xy = INPUT_xprim_->clone();
  Real alpha = two*static_cast<Real>(rand())/static_cast<Real>(RAND_MAX)-one;
  xy->set(*x); xy->axpy(alpha,*y);
  Ptr<Vector<Real>> c1, c2;
  for (auto it = INPUT_linear_con_.begin(); it != INPUT_linear_con_.end(); ++it) {
    c1 = it->second.residual->clone();
    c2 = it->second.residual->clone();
    it->second.constraint->update(*xy,UpdateType::Temp);
    it->second.constraint->value(*c1,*xy,tol);
    cnorm = c1->norm();
    it->second.constraint->update(*x,UpdateType::Temp);
    it->second.constraint->value(*c2,*x,tol);
    c1->axpy(-one,*c2);
    it->second.constraint->update(*y,UpdateType::Temp);
    it->second.constraint->value(*c2,*y,tol);
    c1->axpy(-alpha,*c2);
    it->second.constraint->update(*z,UpdateType::Temp);
    it->second.constraint->value(*c2,*z,tol);
    c1->axpy(alpha,*c2);
    err = c1->norm();
    maxerr = std::max(err,maxerr);
    if (printToStream) {
      outStream << "    Constraint " << it->first;
      outStream << ":  ||c(x+alpha*y) - (c(x)+alpha*(c(y)-c(0)))|| = " << err << std::endl;
      if (err > eps*cnorm) {
        outStream << "      Constraint " << it->first << " may not be linear!" << std::endl;
      }
    }
  }
  if (printToStream) {
    outStream << std::endl;
  }
  outStream.flags(state);
  return maxerr;
}

template<typename Real>
void Problem<Real>::checkVectors(bool printToStream, std::ostream &outStream) const {
  const Real one(1);
  Ptr<Vector<Real>> x, y;
  // Primal optimization space vector
  x = INPUT_xprim_->clone(); x->randomize(-one,one);
  y = INPUT_xprim_->clone(); y->randomize(-one,one);
  if (printToStream) {
    outStream << std::endl << "  Check primal optimization space vector" << std::endl;
  }
  INPUT_xprim_->checkVector(*x,*y,printToStream,outStream);

  // Dual optimization space vector
  x = INPUT_xdual_->clone(); x->randomize(-one,one);
  y = INPUT_xdual_->clone(); y->randomize(-one,one);
  if (printToStream) {
    outStream << std::endl << "  Check dual optimization space vector" << std::endl;
  }
  INPUT_xdual_->checkVector(*x,*y,printToStream,outStream);
  
  // Check constraint space vectors
  for (auto it = INPUT_con_.begin(); it != INPUT_con_.end(); ++it) {
    // Primal constraint space vector
    x = it->second.residual->clone(); x->randomize(-one,one);
    y = it->second.residual->clone(); y->randomize(-one,one);
    if (printToStream) {
      outStream << std::endl << "  " << it->first << ": Check primal constraint space vector" << std::endl;
    }
    it->second.residual->checkVector(*x,*y,printToStream,outStream);

    // Dual optimization space vector
    x = it->second.multiplier->clone(); x->randomize(-one,one);
    y = it->second.multiplier->clone(); y->randomize(-one,one);
    if (printToStream) {
      outStream << std::endl << "  " << it->first << ": Check dual constraint space vector" << std::endl;
    }
    it->second.multiplier->checkVector(*x,*y,printToStream,outStream);
  }
  
  // Check constraint space vectors
  for (auto it = INPUT_linear_con_.begin(); it != INPUT_linear_con_.end(); ++it) {
    // Primal constraint space vector
    x = it->second.residual->clone(); x->randomize(-one,one);
    y = it->second.residual->clone(); y->randomize(-one,one);
    if (printToStream) {
      outStream << std::endl << "  " << it->first << ": Check primal linear constraint space vector" << std::endl;
    }
    it->second.residual->checkVector(*x,*y,printToStream,outStream);

    // Dual optimization space vector
    x = it->second.multiplier->clone(); x->randomize(-one,one);
    y = it->second.multiplier->clone(); y->randomize(-one,one);
    if (printToStream) {
      outStream << std::endl << "  " << it->first << ": Check dual linear constraint space vector" << std::endl;
    }
    it->second.multiplier->checkVector(*x,*y,printToStream,outStream);
  }
}

template<typename Real>
void Problem<Real>::checkDerivatives(bool printToStream, std::ostream &outStream) const {
  const Real one(1);
  Ptr<Vector<Real>> x, d, v, g, c, w;
  // Objective check
  x = INPUT_xprim_->clone(); x->randomize(-one,one);
  d = INPUT_xprim_->clone(); d->randomize(-one,one);
  v = INPUT_xprim_->clone(); v->randomize(-one,one);
  g = INPUT_xdual_->clone(); g->randomize(-one,one);
  if (printToStream) {
    outStream << std::endl << "  Check objective function" << std::endl << std::endl;
  }
  INPUT_obj_->checkGradient(*x,*g,*d,printToStream,outStream);
  INPUT_obj_->checkHessVec(*x,*g,*d,printToStream,outStream);
  INPUT_obj_->checkHessSym(*x,*g,*d,*v,printToStream,outStream);
  
  // Constraint check
  for (auto it = INPUT_con_.begin(); it != INPUT_con_.end(); ++it) {
    c = it->second.residual->clone();   c->randomize(-one,one);
    w = it->second.multiplier->clone(); w->randomize(-one,one);
    if (printToStream) {
      outStream << std::endl << "  " << it->first << ": Check constraint function" << std::endl << std::endl;
    }
    it->second.constraint->checkApplyJacobian(*x,*v,*c,printToStream,outStream);
    it->second.constraint->checkAdjointConsistencyJacobian(*w,*v,*x,printToStream,outStream);
    it->second.constraint->checkApplyAdjointHessian(*x,*w,*v,*g,printToStream,outStream);
  }
  
  // Linear constraint check
  for (auto it = INPUT_linear_con_.begin(); it != INPUT_linear_con_.end(); ++it) {
    c = it->second.residual->clone();   c->randomize(-one,one);
    w = it->second.multiplier->clone(); w->randomize(-one,one);
    if (printToStream) {
      outStream << std::endl << "  " << it->first << ": Check constraint function" << std::endl << std::endl;
    }
    it->second.constraint->checkApplyJacobian(*x,*v,*c,printToStream,outStream);
    it->second.constraint->checkAdjointConsistencyJacobian(*w,*v,*x,printToStream,outStream);
    it->second.constraint->checkApplyAdjointHessian(*x,*w,*v,*g,printToStream,outStream);
  }
}

template<typename Real>
void Problem<Real>::check(bool printToStream, std::ostream &outStream) const {
  checkVectors(printToStream,outStream);
  if (hasLinearEquality_ || hasLinearInequality_) {
    checkLinearity(printToStream,outStream);
  }
  checkDerivatives(printToStream,outStream);
}

template<typename Real>
bool Problem<Real>::isFinalized() const {
  return isFinalized_;
}

template<typename Real>
void Problem<Real>::edit() {
  isFinalized_ = false;
  rlc_  = nullPtr;
  proj_ = nullPtr;
}

template<typename Real>
void Problem<Real>::finalizeIteration() {
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
