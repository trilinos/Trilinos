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

#ifndef ROL_OPTIMIZATIONPROBLEM_HPP
#define ROL_OPTIMIZATIONPROBLEM_HPP

#include "ROL_ConstraintManager.hpp"
#include "ROL_SlacklessObjective.hpp"
#include "ROL_RandomVector.hpp"

// Stochastic Includes
#include "ROL_SampleGenerator.hpp"
#include "ROL_RiskVector.hpp"
#include "ROL_RiskBoundConstraint.hpp"
// Objective includes
#include "ROL_RiskNeutralObjective.hpp" 
#include "ROL_StochasticObjective.hpp"
#include "ROL_RiskLessObjective.hpp"
// Constraint includes
#include "ROL_RiskNeutralConstraint.hpp" 
#include "ROL_StochasticConstraint.hpp" 
#include "ROL_RiskLessConstraint.hpp"
// Almost sure constraint includes
#include "ROL_AlmostSureConstraint.hpp"
#include "ROL_SimulatedBoundConstraint.hpp"
#include "ROL_SimulatedVector.hpp"

namespace ROL {

/* Represents optimization problems in Type-EB form 
 */

template<class Real>
class OptimizationProblem {
private:
  ROL::Ptr<Objective<Real> >                     INPUT_obj_;
  ROL::Ptr<Vector<Real> >                        INPUT_sol_;
  ROL::Ptr<BoundConstraint<Real> >               INPUT_bnd_;
  std::vector<ROL::Ptr<Constraint<Real> > >      INPUT_econ_;
  std::vector<ROL::Ptr<Vector<Real> > >          INPUT_emul_;
  std::vector<ROL::Ptr<Constraint<Real> > >      INPUT_icon_;
  std::vector<ROL::Ptr<Vector<Real> > >          INPUT_imul_;
  std::vector<ROL::Ptr<BoundConstraint<Real> > > INPUT_ibnd_;

  ROL::Ptr<Objective<Real> >                     INTERMEDIATE_obj_;
  ROL::Ptr<Vector<Real> >                        INTERMEDIATE_sol_;
  ROL::Ptr<BoundConstraint<Real> >               INTERMEDIATE_bnd_;
  std::vector<ROL::Ptr<Constraint<Real> > >      INTERMEDIATE_econ_;
  std::vector<ROL::Ptr<Vector<Real> > >          INTERMEDIATE_emul_;
  std::vector<ROL::Ptr<Constraint<Real> > >      INTERMEDIATE_icon_;
  std::vector<ROL::Ptr<Vector<Real> > >          INTERMEDIATE_imul_;
  std::vector<ROL::Ptr<BoundConstraint<Real> > > INTERMEDIATE_ibnd_;

  ROL::Ptr<SampleGenerator<Real> >               vsampler_;
  ROL::Ptr<SampleGenerator<Real> >               gsampler_;
  ROL::Ptr<SampleGenerator<Real> >               hsampler_;
  std::vector<ROL::Ptr<SampleGenerator<Real> > > exsampler_;
  std::vector<ROL::Ptr<BatchManager<Real> > >    ecbman_;
  std::vector<ROL::Ptr<SampleGenerator<Real> > > ixsampler_;
  std::vector<ROL::Ptr<BatchManager<Real> > >    icbman_;

  Teuchos::RCP<Teuchos::ParameterList>               parlistObj_;
  std::vector<Teuchos::RCP<Teuchos::ParameterList> > parlistCon_;

  ROL::Ptr<Objective<Real> >       obj_;
  ROL::Ptr<Vector<Real> >          sol_;
  ROL::Ptr<BoundConstraint<Real> > bnd_;
  ROL::Ptr<Constraint<Real> >      con_;
  ROL::Ptr<Vector<Real> >          mul_;

  ROL::Ptr<ConstraintManager<Real> > conManager_;

  EProblem problemType_;

  bool isInitialized_;

  bool              needRiskLessObj_;
  std::vector<bool> needRiskLessEcon_;
  std::vector<bool> needRiskLessIcon_;
  bool              isStochastic_;

  void initialize( const ROL::Ptr<Objective<Real> >                     &obj,
                   const ROL::Ptr<Vector<Real> >                        &x,
                   const ROL::Ptr<BoundConstraint<Real> >               &bnd,
                   const std::vector<ROL::Ptr<Constraint<Real> > >      &econ,
                   const std::vector<ROL::Ptr<Vector<Real> > >          &emul,
                   const std::vector<ROL::Ptr<Constraint<Real> > >      &icon,
                   const std::vector<ROL::Ptr<Vector<Real> > >          &imul,
                   const std::vector<ROL::Ptr<BoundConstraint<Real> > > &ibnd ) {
    if (!isInitialized_) {
      int esize = static_cast<int>(econ.size());
      int isize = static_cast<int>(icon.size());
      std::vector<ROL::Ptr<Constraint<Real> > >      cvec;
      std::vector<ROL::Ptr<Vector<Real> > >          lvec;
      std::vector<ROL::Ptr<BoundConstraint<Real> > > bvec;
      for (int i = 0; i < esize; ++i) {
        if ( econ[i] != ROL::nullPtr ) {
          if (isStochastic_) {
            cvec.push_back(setRiskLessCon(econ[i],needRiskLessEcon_[i]));
          }
          else {
            cvec.push_back(econ[i]);
          }
          lvec.push_back(emul[i]);
          bvec.push_back(ROL::nullPtr);
        }
      }
      for (int i = 0; i < isize; ++i) {
        if ( icon[i] != ROL::nullPtr ) {
          if (isStochastic_) {
            cvec.push_back(setRiskLessCon(icon[i],needRiskLessIcon_[i]));
          }
          else {
            cvec.push_back(icon[i]);
          }
          lvec.push_back(imul[i]);
          bvec.push_back(ibnd[i]);
        }
      }
 
      conManager_ = ROL::makePtr<ConstraintManager<Real>>(cvec,lvec,bvec,x,bnd);
      con_        = conManager_->getConstraint();
      mul_        = conManager_->getMultiplier();
      sol_        = conManager_->getOptVector();
      bnd_        = conManager_->getBoundConstraint();
      ROL::Ptr<Objective<Real> > obj0;
      if (isStochastic_) {
        obj0 = setRiskLessObj(obj,needRiskLessObj_);
      }
      else {
        obj0 = obj;
      }
      if ( conManager_->hasInequality() ) {
        obj_      = ROL::makePtr<SlacklessObjective<Real>>(obj0);
      }
      else {
        obj_      = obj0;
      }

      if ( conManager_->isNull() ) {
        if( bnd_ == ROL::nullPtr || !bnd_->isActivated() ) {  // Type-U
          problemType_ = TYPE_U;        
        }
        else { // Type-B
          problemType_ = TYPE_B; 
        }
      }
      else {
        if( bnd_ == ROL::nullPtr || !bnd_->isActivated() ) { // Type-E
          problemType_ = TYPE_E;     
        }
        else { // Type-EB
          problemType_ = TYPE_EB; 
        }
      }
      isInitialized_ = true;
    }
  }

  const ROL::Ptr<Constraint<Real> > setRiskLessCon(const ROL::Ptr<Constraint<Real> > &con, const bool needRiskLess) const {
    if (needRiskLess) {
      return ROL::makePtr<RiskLessConstraint<Real>>(con);
    }
    else {
      return con;
    }
  }

  const ROL::Ptr<Objective<Real> > setRiskLessObj(const ROL::Ptr<Objective<Real> > &obj, const bool needRiskLess) const {
    if (needRiskLess) {
      return ROL::makePtr<RiskLessObjective<Real>>(obj);
    }
    else {
      return obj;
    }
  }

  std::vector<Real> computeSampleMean(const ROL::Ptr<SampleGenerator<Real> > &sampler) const {
    // Compute mean value of inputs and set parameter in objective
    int dim = sampler->getMyPoint(0).size(), nsamp = sampler->numMySamples();
    std::vector<Real> loc(dim), mean(dim), pt(dim);
    Real wt(0);
    for (int i = 0; i < nsamp; i++) {
      pt = sampler->getMyPoint(i);
      wt = sampler->getMyWeight(i);
      for (int j = 0; j < dim; j++) {
        loc[j] += wt*pt[j];
      }
    }
    sampler->sumAll(&loc[0],&mean[0],dim);
    return mean;
  }

  void initStochastic(void) {
    if (!isStochastic_) {
      int econSize = INPUT_econ_.size();
      int iconSize = INPUT_icon_.size();
      needRiskLessObj_ = true;
      needRiskLessEcon_.clear(); needRiskLessEcon_.resize(econSize,true);
      needRiskLessIcon_.clear(); needRiskLessIcon_.resize(iconSize,true);
      parlistObj_ = ROL::nullPtr;
      parlistCon_.clear(); parlistCon_.resize(iconSize,ROL::nullPtr);

      exsampler_.clear(); exsampler_.resize(econSize,ROL::nullPtr);
      ecbman_.clear();    ecbman_.resize(econSize,ROL::nullPtr);

      ixsampler_.clear(); ixsampler_.resize(iconSize,ROL::nullPtr);
      icbman_.clear();    icbman_.resize(iconSize,ROL::nullPtr);

      isStochastic_ = true;
    }
  }

  void buildRiskVec(ROL::Ptr<Vector<Real> > &x) {
    // Build risk vector and risk bound constraint
    INTERMEDIATE_sol_
      = ROL::makePtr<RiskVector<Real>>(parlistObj_,parlistCon_,x);
    if (parlistObj_ != ROL::nullPtr) {
      Real statObj = parlistObj_->sublist("SOL").get("Initial Statistic",1.0);
      ROL::dynamicPtrCast<RiskVector<Real> >(INTERMEDIATE_sol_)->setStatistic(statObj,0);
    }
    int nc = INPUT_icon_.size();
    for (int i = 0; i < nc; ++i) {
      if (parlistCon_[i] != ROL::nullPtr) {
        Real statCon = parlistCon_[i]->sublist("SOL").get("Initial Statistic",1.0);
        ROL::dynamicPtrCast<RiskVector<Real> >(INTERMEDIATE_sol_)->setStatistic(statCon,1,i);
      }
    }
  }

  void buildRiskBnd(ROL::Ptr<BoundConstraint<Real> > &bnd) {
    if ( INPUT_bnd_ != ROL::nullPtr ) {
      INTERMEDIATE_bnd_
        = ROL::makePtr<RiskBoundConstraint<Real>>(parlistObj_,parlistCon_,bnd);
    }
  }

public:
  virtual ~OptimizationProblem(void) {}

  // Default constructor [0]
  OptimizationProblem(void)
   : isInitialized_(false), isStochastic_(false) {}

  // Complete option constructor [1]
  OptimizationProblem( const ROL::Ptr<Objective<Real> >                     &obj,
                       const ROL::Ptr<Vector<Real> >                        &x,
                       const ROL::Ptr<BoundConstraint<Real> >               &bnd,
                       const std::vector<ROL::Ptr<Constraint<Real> > >      &econ,
                       const std::vector<ROL::Ptr<Vector<Real> > >          &emul,
                       const std::vector<ROL::Ptr<Constraint<Real> > >      &icon,
                       const std::vector<ROL::Ptr<Vector<Real> > >          &imul,
                       const std::vector<ROL::Ptr<BoundConstraint<Real> > > &ibnd )
    : INPUT_obj_(obj), INPUT_sol_(x), INPUT_bnd_(bnd),
      INPUT_econ_(econ), INPUT_emul_(emul),
      INPUT_icon_(icon), INPUT_imul_(imul), INPUT_ibnd_(ibnd),
      INTERMEDIATE_obj_(obj), INTERMEDIATE_sol_(x), INTERMEDIATE_bnd_(bnd),
      INTERMEDIATE_econ_(econ), INTERMEDIATE_emul_(emul),
      INTERMEDIATE_icon_(icon), INTERMEDIATE_imul_(imul), INTERMEDIATE_ibnd_(ibnd),
      isInitialized_(false), isStochastic_(false) {}

  OptimizationProblem( const ROL::Ptr<Objective<Real> >                     &obj,
                       const ROL::Ptr<Vector<Real> >                        &x,
                       const ROL::Ptr<BoundConstraint<Real> >               &bnd,
                       const ROL::Ptr<Constraint<Real> >                    &econ,
                       const ROL::Ptr<Vector<Real> >                        &emul,
                       const std::vector<ROL::Ptr<Constraint<Real> > >      &icon,
                       const std::vector<ROL::Ptr<Vector<Real> > >          &imul,
                       const std::vector<ROL::Ptr<BoundConstraint<Real> > > &ibnd )
    : INPUT_obj_(obj), INPUT_sol_(x), INPUT_bnd_(bnd),
      INPUT_icon_(icon), INPUT_imul_(imul), INPUT_ibnd_(ibnd),
      INTERMEDIATE_obj_(obj), INTERMEDIATE_sol_(x), INTERMEDIATE_bnd_(bnd),
      INTERMEDIATE_icon_(icon), INTERMEDIATE_imul_(imul), INTERMEDIATE_ibnd_(ibnd),
      isInitialized_(false), isStochastic_(false) {
    std::vector<ROL::Ptr<Constraint<Real> > > econ0(1,econ);
    std::vector<ROL::Ptr<Vector<Real> > >     emul0(1,emul);
    INPUT_econ_ = econ0;
    INPUT_emul_ = emul0;
    INTERMEDIATE_econ_ = econ0;
    INTERMEDIATE_emul_ = emul0;
  }

  OptimizationProblem( const ROL::Ptr<Objective<Real> >                     &obj,
                       const ROL::Ptr<Vector<Real> >                        &x,
                       const ROL::Ptr<BoundConstraint<Real> >               &bnd,
                       const std::vector<ROL::Ptr<Constraint<Real> > >      &econ,
                       const std::vector<ROL::Ptr<Vector<Real> > >          &emul,
                       const ROL::Ptr<Constraint<Real> >                    &icon,
                       const ROL::Ptr<Vector<Real> >                        &imul,
                       const ROL::Ptr<BoundConstraint<Real> >               &ibnd )
    : INPUT_obj_(obj), INPUT_sol_(x), INPUT_bnd_(bnd),
      INPUT_econ_(econ), INPUT_emul_(emul),
      INTERMEDIATE_obj_(obj), INTERMEDIATE_sol_(x), INTERMEDIATE_bnd_(bnd),
      INTERMEDIATE_econ_(econ), INTERMEDIATE_emul_(emul),
      isInitialized_(false), isStochastic_(false) {
    std::vector<ROL::Ptr<Constraint<Real> > >      icon0(1,icon);
    std::vector<ROL::Ptr<Vector<Real> > >          imul0(1,imul);
    std::vector<ROL::Ptr<BoundConstraint<Real> > > ibnd0(1,ibnd);
    INPUT_icon_ = icon0;
    INPUT_imul_ = imul0;
    INPUT_ibnd_ = ibnd0;
    INTERMEDIATE_icon_ = icon0;
    INTERMEDIATE_imul_ = imul0;
    INTERMEDIATE_ibnd_ = ibnd0;
  }

  OptimizationProblem( const ROL::Ptr<Objective<Real> >                     &obj,
                       const ROL::Ptr<Vector<Real> >                        &x,
                       const ROL::Ptr<BoundConstraint<Real> >               &bnd,
                       const ROL::Ptr<Constraint<Real> >                    &econ,
                       const ROL::Ptr<Vector<Real> >                        &emul,
                       const ROL::Ptr<Constraint<Real> >                    &icon,
                       const ROL::Ptr<Vector<Real> >                        &imul,
                       const ROL::Ptr<BoundConstraint<Real> >               &ibnd )
    : INPUT_obj_(obj), INPUT_sol_(x), INPUT_bnd_(bnd),
      INTERMEDIATE_obj_(obj), INTERMEDIATE_sol_(x), INTERMEDIATE_bnd_(bnd),
      isInitialized_(false), isStochastic_(false) {
    std::vector<ROL::Ptr<Constraint<Real> > >      econ0(1,econ);
    std::vector<ROL::Ptr<Vector<Real> > >          emul0(1,emul);
    std::vector<ROL::Ptr<Constraint<Real> > >      icon0(1,icon);
    std::vector<ROL::Ptr<Vector<Real> > >          imul0(1,imul);
    std::vector<ROL::Ptr<BoundConstraint<Real> > > ibnd0(1,ibnd);
    INPUT_econ_ = econ0;
    INPUT_emul_ = emul0;
    INPUT_icon_ = icon0;
    INPUT_imul_ = imul0;
    INPUT_ibnd_ = ibnd0;
    INTERMEDIATE_econ_ = econ0;
    INTERMEDIATE_emul_ = emul0;
    INTERMEDIATE_icon_ = icon0;
    INTERMEDIATE_imul_ = imul0;
    INTERMEDIATE_ibnd_ = ibnd0;
  }

  // No bound constuctor [2]
  OptimizationProblem( const ROL::Ptr<Objective<Real> >                     &obj,
                       const ROL::Ptr<Vector<Real> >                        &x,
                       const std::vector<ROL::Ptr<Constraint<Real> > >      &econ,
                       const std::vector<ROL::Ptr<Vector<Real> > >          &emul,
                       const std::vector<ROL::Ptr<Constraint<Real> > >      &icon,
                       const std::vector<ROL::Ptr<Vector<Real> > >          &imul,
                       const std::vector<ROL::Ptr<BoundConstraint<Real> > > &ibnd )
    : OptimizationProblem( obj, x, ROL::nullPtr, econ, emul, icon, imul, ibnd ) {}

  OptimizationProblem( const ROL::Ptr<Objective<Real> >                     &obj,
                       const ROL::Ptr<Vector<Real> >                        &x,
                       const ROL::Ptr<Constraint<Real> >                    &econ,
                       const ROL::Ptr<Vector<Real> >                        &emul,
                       const std::vector<ROL::Ptr<Constraint<Real> > >      &icon,
                       const std::vector<ROL::Ptr<Vector<Real> > >          &imul,
                       const std::vector<ROL::Ptr<BoundConstraint<Real> > > &ibnd )
    : OptimizationProblem( obj, x, ROL::nullPtr, econ, emul, icon, imul, ibnd) {}

  OptimizationProblem( const ROL::Ptr<Objective<Real> >                     &obj,
                       const ROL::Ptr<Vector<Real> >                        &x,
                       const std::vector<ROL::Ptr<Constraint<Real> > >      &econ,
                       const std::vector<ROL::Ptr<Vector<Real> > >          &emul,
                       const ROL::Ptr<Constraint<Real> >                    &icon,
                       const ROL::Ptr<Vector<Real> >                        &imul,
                       const ROL::Ptr<BoundConstraint<Real> >               &ibnd )
    : OptimizationProblem( obj, x, ROL::nullPtr, econ, emul, icon, imul, ibnd) {}

  OptimizationProblem( const ROL::Ptr<Objective<Real> >                     &obj,
                       const ROL::Ptr<Vector<Real> >                        &x,
                       const ROL::Ptr<Constraint<Real> >                    &econ,
                       const ROL::Ptr<Vector<Real> >                        &emul,
                       const ROL::Ptr<Constraint<Real> >                    &icon,
                       const ROL::Ptr<Vector<Real> >                        &imul,
                       const ROL::Ptr<BoundConstraint<Real> >               &ibnd )
    : OptimizationProblem( obj, x, ROL::nullPtr, econ, emul, icon, imul, ibnd) {}

  // No inequality constraint [3]
  OptimizationProblem( const ROL::Ptr<Objective<Real> >                     &obj,
                       const ROL::Ptr<Vector<Real> >                        &x,
                       const ROL::Ptr<BoundConstraint<Real> >               &bnd,
                       const std::vector<ROL::Ptr<Constraint<Real> > >      &econ,
                       const std::vector<ROL::Ptr<Vector<Real> > >          &emul )
    : OptimizationProblem( obj, x, bnd, econ, emul, ROL::nullPtr, ROL::nullPtr, ROL::nullPtr ) {}

  OptimizationProblem( const ROL::Ptr<Objective<Real> >                     &obj,
                       const ROL::Ptr<Vector<Real> >                        &x,
                       const ROL::Ptr<BoundConstraint<Real> >               &bnd,
                       const ROL::Ptr<Constraint<Real> >                    &econ,
                       const ROL::Ptr<Vector<Real> >                        &emul )
    : OptimizationProblem( obj, x, bnd, econ, emul, ROL::nullPtr, ROL::nullPtr, ROL::nullPtr) {}

  // No equality constraint [4]
  OptimizationProblem( const ROL::Ptr<Objective<Real> >                     &obj,
                       const ROL::Ptr<Vector<Real> >                        &x,
                       const ROL::Ptr<BoundConstraint<Real> >               &bnd,
                       const std::vector<ROL::Ptr<Constraint<Real> > >      &icon,
                       const std::vector<ROL::Ptr<Vector<Real> > >          &imul,
                       const std::vector<ROL::Ptr<BoundConstraint<Real> > > &ibnd )
    : OptimizationProblem( obj, x, bnd, ROL::nullPtr, ROL::nullPtr, icon, imul, ibnd ) {}

  OptimizationProblem( const ROL::Ptr<Objective<Real> >                     &obj,
                       const ROL::Ptr<Vector<Real> >                        &x,
                       const ROL::Ptr<BoundConstraint<Real> >               &bnd,
                       const ROL::Ptr<Constraint<Real> >                    &icon,
                       const ROL::Ptr<Vector<Real> >                        &imul,
                       const ROL::Ptr<BoundConstraint<Real> >               &ibnd )
    : OptimizationProblem( obj, x, bnd, ROL::nullPtr, ROL::nullPtr, icon, imul, ibnd) {}

  // No inequality or bound constraint [5]
  OptimizationProblem( const ROL::Ptr<Objective<Real> >                     &obj,
                       const ROL::Ptr<Vector<Real> >                        &x,
                       const std::vector<ROL::Ptr<Constraint<Real> > >      &econ,
                       const std::vector<ROL::Ptr<Vector<Real> > >          &emul )
    : OptimizationProblem( obj, x, ROL::nullPtr, econ, emul, ROL::nullPtr, ROL::nullPtr, ROL::nullPtr ) {}

  OptimizationProblem( const ROL::Ptr<Objective<Real> >                     &obj,
                       const ROL::Ptr<Vector<Real> >                        &x,
                       const ROL::Ptr<Constraint<Real> >                    &econ,
                       const ROL::Ptr<Vector<Real> >                        &emul )
    : OptimizationProblem( obj, x, ROL::nullPtr, econ, emul, ROL::nullPtr, ROL::nullPtr, ROL::nullPtr) {}

  // No equality or bound constraint [6]
  OptimizationProblem( const ROL::Ptr<Objective<Real> >                     &obj,
                       const ROL::Ptr<Vector<Real> >                        &x,
                       const std::vector<ROL::Ptr<Constraint<Real> > >      &icon,
                       const std::vector<ROL::Ptr<Vector<Real> > >          &imul,
                       const std::vector<ROL::Ptr<BoundConstraint<Real> > > &ibnd )
    : OptimizationProblem( obj, x, ROL::nullPtr, ROL::nullPtr, ROL::nullPtr, icon, imul, ibnd ) {}

  OptimizationProblem( const ROL::Ptr<Objective<Real> >                     &obj,
                       const ROL::Ptr<Vector<Real> >                        &x,
                       const ROL::Ptr<Constraint<Real> >                    &icon,
                       const ROL::Ptr<Vector<Real> >                        &imul,
                       const ROL::Ptr<BoundConstraint<Real> >               &ibnd )
    : OptimizationProblem( obj, x, ROL::nullPtr, ROL::nullPtr, ROL::nullPtr, icon, imul, ibnd) {}

  // Bound constrained problem [7]
  OptimizationProblem( const ROL::Ptr<Objective<Real> >                     &obj,
                       const ROL::Ptr<Vector<Real> >                        &x,
                       const ROL::Ptr<BoundConstraint<Real> >               &bnd )
    : OptimizationProblem( obj, x, bnd, ROL::nullPtr, ROL::nullPtr, ROL::nullPtr, ROL::nullPtr, ROL::nullPtr ) {}

  // Unconstrained problem [8]
  OptimizationProblem( const ROL::Ptr<Objective<Real> > &obj,
                       const ROL::Ptr<Vector<Real> >    &x ) :
     OptimizationProblem( obj, x, ROL::nullPtr, ROL::nullPtr, ROL::nullPtr, ROL::nullPtr, ROL::nullPtr, ROL::nullPtr ) {} 

  /* Get methods */

  virtual ROL::Ptr<Objective<Real> > getObjective(void) {
    if ( INTERMEDIATE_obj_ == ROL::nullPtr ) {
      throw Exception::NotImplemented(">>> ROL::OptimizationProblem::getObjective: No objective inputed!");
    }
    initialize(INTERMEDIATE_obj_,INTERMEDIATE_sol_,INTERMEDIATE_bnd_,
               INTERMEDIATE_econ_,INTERMEDIATE_emul_,
               INTERMEDIATE_icon_,INTERMEDIATE_imul_,INTERMEDIATE_ibnd_);
    return obj_;
  }

  virtual ROL::Ptr<Vector<Real> > getSolutionVector(void) {
    if ( INTERMEDIATE_sol_ == ROL::nullPtr ) {
      throw Exception::NotImplemented(">>> ROL::OptimizationProblem::getSolutionVector: No solution vector inputed!");
    }
    initialize(INTERMEDIATE_obj_,INTERMEDIATE_sol_,INTERMEDIATE_bnd_,
               INTERMEDIATE_econ_,INTERMEDIATE_emul_,
               INTERMEDIATE_icon_,INTERMEDIATE_imul_,INTERMEDIATE_ibnd_);
    return sol_;
  }

  virtual ROL::Ptr<BoundConstraint<Real> > getBoundConstraint(void) {
    initialize(INTERMEDIATE_obj_,INTERMEDIATE_sol_,INTERMEDIATE_bnd_,
               INTERMEDIATE_econ_,INTERMEDIATE_emul_,
               INTERMEDIATE_icon_,INTERMEDIATE_imul_,INTERMEDIATE_ibnd_);
    return bnd_;
  }

  virtual ROL::Ptr<Constraint<Real> > getConstraint(void) {
    initialize(INTERMEDIATE_obj_,INTERMEDIATE_sol_,INTERMEDIATE_bnd_,
               INTERMEDIATE_econ_,INTERMEDIATE_emul_,
               INTERMEDIATE_icon_,INTERMEDIATE_imul_,INTERMEDIATE_ibnd_);
    return con_;
  }

  virtual ROL::Ptr<Vector<Real> > getMultiplierVector(void) {
    initialize(INTERMEDIATE_obj_,INTERMEDIATE_sol_,INTERMEDIATE_bnd_,
               INTERMEDIATE_econ_,INTERMEDIATE_emul_,
               INTERMEDIATE_icon_,INTERMEDIATE_imul_,INTERMEDIATE_ibnd_);
    return mul_;
  }

  EProblem getProblemType(void) {
    initialize(INTERMEDIATE_obj_,INTERMEDIATE_sol_,INTERMEDIATE_bnd_,
               INTERMEDIATE_econ_,INTERMEDIATE_emul_,
               INTERMEDIATE_icon_,INTERMEDIATE_imul_,INTERMEDIATE_ibnd_);
    return problemType_;
  }

  /* Set Stochastic Methods */

  /* Objective function */
  /** \brief Set objective function to mean value objective.

      We assume the objective function is parametrized by an additional
      variable (other than the optimization variable).  This variable could,
      e.g., be random.  The mean value objective function evaluates the
      the parametrized objective function at the sample average of the
      auxiliary variable.

      @param[in]    sampler  is the SampleGenerator defining the distribution of the auxiliary parameter
  */
  void setMeanValueObjective(const ROL::Ptr<SampleGenerator<Real> > &sampler) {
    initStochastic();
    // Set objective function samplers
    vsampler_ = sampler;
    gsampler_ = sampler;
    hsampler_ = sampler;
    // Construct risk-averse/probabilistic objective function
    if ( vsampler_ == ROL::nullPtr ) {
      throw Exception::NotImplemented(">>> ROL::OptimizationProblem::setMeanValueObjective: Objective function value sampler is null!");
    }
    else {
      std::vector<Real> mean = computeSampleMean(vsampler_);
      INTERMEDIATE_obj_ = INPUT_obj_;
      INTERMEDIATE_obj_->setParameter(mean);
    }
    // Set vector and bound constraint
    buildRiskVec(INPUT_sol_);
    buildRiskBnd(INPUT_bnd_);

    isInitialized_ = false;
  }

  /** \brief Set objective function to risk neutral objective.

      We assume the objective function is parametrized by an additional
      variable (other than the optimization variable).  This variable could,
      e.g., be random.  The risk neutral objective function evaluates the
      the average of parametrized objective function.

      @param[in]    vsampler  is the SampleGenerator defining the distribution of the auxiliary parameter for the value
      @param[in]    gsampler  is the SampleGenerator defining the distribution of the auxiliary parameter for the gradient
      @param[in]    hsampler  is the SampleGenerator defining the distribution of the auxiliary parameter for the Hessian
      @param[in]    storage   whether or not to store the sampled value and gradient
  */
  void setRiskNeutralObjective(const ROL::Ptr<SampleGenerator<Real> > &vsampler,
                               const ROL::Ptr<SampleGenerator<Real> > &gsampler = ROL::nullPtr,
                               const ROL::Ptr<SampleGenerator<Real> > &hsampler = ROL::nullPtr,
                               const bool storage = true) {
    initStochastic();
    // Set objective function samplers
    vsampler_ = vsampler;
    gsampler_ = gsampler;
    hsampler_ = hsampler;
    if ( gsampler == ROL::nullPtr ) {
      gsampler_ = vsampler_;
    }
    if ( hsampler == ROL::nullPtr ) {
      hsampler_ = gsampler_;
    }
    // Construct risk-averse/probabilistic objective function
    if ( vsampler_ == ROL::nullPtr ) {
      throw Exception::NotImplemented(">>> ROL::OptimizationProblem::setRiskNeutralObjective: Objective function value sampler is null!");
    }
    else {
      INTERMEDIATE_obj_
        = ROL::makePtr<RiskNeutralObjective<Real>>(INPUT_obj_,vsampler_,gsampler_,hsampler_,storage);
    }
    // Set vector and bound constraint
    buildRiskVec(INPUT_sol_);
    buildRiskBnd(INPUT_bnd_);

    isInitialized_ = false;
  }

  /** \brief Set objective function to risk averse objective.

      We assume the objective function is parametrized by an additional
      variable (other than the optimization variable).  This variable could,
      e.g., be random.  The risk averse objective function evaluates the
      the ``risk'' of parametrized objective function.

      @param[in]    parlist   contains the information defining the risk measure
      @param[in]    vsampler  is the SampleGenerator defining the distribution of the auxiliary parameter for the value
      @param[in]    gsampler  is the SampleGenerator defining the distribution of the auxiliary parameter for the gradient
      @param[in]    hsampler  is the SampleGenerator defining the distribution of the auxiliary parameter for the Hessian
  */
  void setRiskAverseObjective(Teuchos::ParameterList &parlist,
                              const ROL::Ptr<SampleGenerator<Real> > &vsampler,
                              const ROL::Ptr<SampleGenerator<Real> > &gsampler = ROL::nullPtr,
                              const ROL::Ptr<SampleGenerator<Real> > &hsampler = ROL::nullPtr) {
    initStochastic();
    // Set objective function samplers
    vsampler_ = vsampler;
    gsampler_ = gsampler;
    hsampler_ = hsampler;
    if ( gsampler == ROL::nullPtr ) {
      gsampler_ = vsampler_;
    }
    if ( hsampler == ROL::nullPtr ) {
      hsampler_ = gsampler_;
    }
    // Construct risk-averse/probabilistic objective function
    if ( vsampler_ == ROL::nullPtr ) {
      throw Exception::NotImplemented(">>> ROL::OptimizationProblem::setRiskAverseObjective: Objective function value sampler is null!");
    }
    else {
      needRiskLessObj_ = false;
      parlistObj_      = ROL::makePtrFromRef(parlist);
      INTERMEDIATE_obj_
        = ROL::makePtr<StochasticObjective<Real>>(INPUT_obj_,parlist,vsampler_,gsampler_,hsampler_);
    }
    // Set vector and bound constraint
    buildRiskVec(INPUT_sol_);
    buildRiskBnd(INPUT_bnd_);

    isInitialized_ = false;
  }

  void setStochasticObjective(Teuchos::ParameterList &parlist,
                              const ROL::Ptr<SampleGenerator<Real> > &vsampler,
                              const ROL::Ptr<SampleGenerator<Real> > &gsampler = ROL::nullPtr,
                              const ROL::Ptr<SampleGenerator<Real> > &hsampler = ROL::nullPtr) {
    // Determine Stochastic Objective Type
    std::string type = parlist.sublist("SOL").get("Stochastic Component Type","Risk Neutral");
    if ( type == "Risk Neutral" ) {
      bool storage = parlist.sublist("SOL").get("Store Sampled Value and Gradient",true);
      setRiskNeutralObjective(vsampler,gsampler,hsampler,storage);
    }
    else if ( type == "Risk Averse" ||
              type == "Deviation"   ||
              type == "Error"       ||
              type == "Regret"      ||
              type == "Probability" ) {
      setRiskAverseObjective(parlist,vsampler,gsampler,hsampler);
    }
    else if ( type == "Mean Value" ) {
      setMeanValueObjective(vsampler);
    }
    else {
      throw Exception::NotImplemented(">>> ROL::OptimizationProblem::setStochasticObjective: Invalid stochastic optimization type!");
    }
    // Set vector and bound constraint
    buildRiskVec(INPUT_sol_);
    buildRiskBnd(INPUT_bnd_);

    isInitialized_ = false;
  }

  /* Equality Constraint */
  void setMeanValueEquality(const ROL::Ptr<SampleGenerator<Real> > &sampler, const int index = 0) {
    initStochastic();
    exsampler_[index] = sampler;
    if ( INPUT_econ_[index] != ROL::nullPtr && sampler != ROL::nullPtr ) {
      std::vector<Real> mean = computeSampleMean(sampler);
      INTERMEDIATE_econ_[index] = INPUT_econ_[index];
      INTERMEDIATE_econ_[index]->setParameter(mean);
      INTERMEDIATE_emul_[index] = INPUT_emul_[index];
    }
    else {
      throw Exception::NotImplemented(">>> ROL::OptimizationProblem::setMeanValueEquality: Either SampleGenerator or Constraint is NULL!");
    }
    // Set vector and bound constraint
    buildRiskVec(INPUT_sol_);
    buildRiskBnd(INPUT_bnd_);

    isInitialized_ = false;
  }

  void setRiskNeutralEquality(const ROL::Ptr<SampleGenerator<Real> > &xsampler,
                              const ROL::Ptr<BatchManager<Real> >    &cbman,
                              const int index = 0) {
    initStochastic();
    exsampler_[index] = xsampler;
    ecbman_[index]    = cbman;
    if ( INPUT_econ_[index] != ROL::nullPtr
         &&        xsampler != ROL::nullPtr
         &&           cbman != ROL::nullPtr ) {
      INTERMEDIATE_econ_[index]
        = ROL::makePtr<RiskNeutralConstraint<Real>>(INPUT_econ_[index],xsampler,cbman);
      INTERMEDIATE_emul_[index] = INPUT_emul_[index];
    }
    else {
      throw Exception::NotImplemented(">>> ROL::OptimizationProblem::SetRiskNeutralEquality: Either SampleGenerator, Constraint or BatchManager is NULL!");
    }
    // Set vector and bound constraint
    buildRiskVec(INPUT_sol_);
    buildRiskBnd(INPUT_bnd_);

    isInitialized_ = false;
  }

  void setAlmostSureEquality(const ROL::Ptr<SampleGenerator<Real> > &sampler,
                             const int index = 0) {
    initStochastic();
    exsampler_[index] = sampler;
    if ( INPUT_econ_[index] != ROL::nullPtr && sampler != ROL::nullPtr ) {
      int nsamp = sampler->numMySamples();
      INTERMEDIATE_econ_[index]
        = ROL::makePtr<AlmostSureConstraint<Real>>(sampler,INPUT_econ_[index]);
      std::vector<ROL::Ptr<Vector<Real> > > emul(nsamp,ROL::nullPtr);
      for (int j = 0; j < nsamp; ++j) {
        emul[j] = INPUT_emul_[index]->clone();
        emul[j]->set(*INPUT_emul_[index]);
      }
      INTERMEDIATE_emul_[index]
        = ROL::makePtr<DualSimulatedVector<Real>>(emul, sampler->getBatchManager(), sampler);
    }
    else {
      throw Exception::NotImplemented(">>> ROL::OptimizationProblem::SetAlmostSureEquality: Either SampleGenerator or Constraint is NULL!");
    }
    // Set vector and bound constraint
    buildRiskVec(INPUT_sol_);
    buildRiskBnd(INPUT_bnd_);

    isInitialized_ = false;
  }


  void setStochasticEquality(std::vector<Teuchos::ParameterList> &parlist,
                             const std::vector<ROL::Ptr<SampleGenerator<Real> > > &xsampler,
                             const std::vector<ROL::Ptr<BatchManager<Real> > > &cbman) {
    initStochastic();
    int nc = static_cast<int>(INPUT_econ_.size());
    if ( nc != static_cast<int>(xsampler.size()) || nc != static_cast<int>(cbman.size()) ) {
      throw Exception::NotImplemented(">>> ROL::OptimizationProblem::setStochasticEquality: Constraint vector and SampleGenerator vector are not the same size!");
    }
    for (int i = 0; i < nc; ++i) {
      if (xsampler[i] != ROL::nullPtr) {
        std::string type = parlist[i].sublist("SOL").get("Stochastic Component Type","Risk Neutral");
        if ( type == "Risk Neutral" ) {
          setRiskNeutralEquality(xsampler[i],cbman[i],i);
        }
        else if ( type == "Almost Sure" ) {
          setAlmostSureEquality(xsampler[i],i);
        }
        else if ( type == "Mean Value" ) {
          setMeanValueEquality(xsampler[i],i);
        }
        else {
          throw Exception::NotImplemented(">>> ROL::OptimizationProblem::SetStochasticEquality: Invalid stochastic constraint type!");
        }
      }
      else {
        INTERMEDIATE_econ_[i] = INPUT_econ_[i];
        INTERMEDIATE_emul_[i] = INPUT_emul_[i];
      }
    }
    // Set vector and bound constraint
    buildRiskVec(INPUT_sol_);
    buildRiskBnd(INPUT_bnd_);

    isInitialized_ = false;
  }

  void setStochasticEquality(Teuchos::ParameterList &parlist,
                             const ROL::Ptr<SampleGenerator<Real> > &xsampler,
                             const ROL::Ptr<BatchManager<Real> > &cbman) {
    std::vector<Teuchos::ParameterList> cparlist(1,parlist);
    std::vector<ROL::Ptr<SampleGenerator<Real> > > cxsampler(1,xsampler);
    std::vector<ROL::Ptr<BatchManager<Real> > > ccbman(1,cbman);
    setStochasticEquality(cparlist,cxsampler,ccbman);
  }

  /* Inequality constraint */
  void setMeanValueInequality(const ROL::Ptr<SampleGenerator<Real> > &sampler,
                              const int index = 0) {
    initStochastic();
    ixsampler_[index] = sampler;
    if ( INPUT_icon_[index] != ROL::nullPtr && sampler != ROL::nullPtr ) {
      std::vector<Real> mean = computeSampleMean(sampler);
      INTERMEDIATE_icon_[index] = INPUT_icon_[index];
      INTERMEDIATE_icon_[index]->setParameter(mean);
      INTERMEDIATE_ibnd_[index] = INPUT_ibnd_[index];
      INTERMEDIATE_imul_[index] = INPUT_imul_[index];
    }
    else {
      throw Exception::NotImplemented(">>> ROL::OptimizationProblem::SetMeanValueInequality: Either Constraint or SampleGenerator is NULL!");
    }
    // Set vector and bound constraint
    buildRiskVec(INPUT_sol_);
    buildRiskBnd(INPUT_bnd_);

    isInitialized_ = false;
  }

  void setRiskNeutralInequality(const ROL::Ptr<SampleGenerator<Real> > &xsampler,
                                const ROL::Ptr<BatchManager<Real> >    &cbman,
                                const int index = 0) {
    initStochastic();
    ixsampler_[index] = xsampler;
    icbman_[index]    = cbman;
    if ( INPUT_icon_[index] != ROL::nullPtr
         &&    xsampler     != ROL::nullPtr
         &&       cbman     != ROL::nullPtr ) {
      INTERMEDIATE_icon_[index]
        = ROL::makePtr<RiskNeutralConstraint<Real>>(INPUT_icon_[index],xsampler,cbman);
      INTERMEDIATE_ibnd_[index] = INPUT_ibnd_[index];
      INTERMEDIATE_imul_[index] = INPUT_imul_[index];
    }
    else {
      throw Exception::NotImplemented(">>> ROL::OptimizationProblem::SetRiskNeutralInequality: Either Constraint, SampleGenerator or BatchManager is NULL!");
    }
    // Set vector and bound constraint
    buildRiskVec(INPUT_sol_);
    buildRiskBnd(INPUT_bnd_);

    isInitialized_ = false;
  }

  void setRiskAverseInequality(Teuchos::ParameterList &parlist,
                               const ROL::Ptr<SampleGenerator<Real> > &sampler,
                               const int index = 0) {
    initStochastic();
    ixsampler_[index] = sampler;
    if ( INPUT_icon_[index] != ROL::nullPtr && sampler != ROL::nullPtr ) {
      needRiskLessIcon_[index] = false;
      parlistCon_[index]       = ROL::makePtrFromRef(parlist);
      INTERMEDIATE_icon_[index]
        = ROL::makePtr<StochasticConstraint<Real>>(INPUT_icon_[index],sampler,parlist,index);
      INTERMEDIATE_ibnd_[index] = INPUT_ibnd_[index];
      INTERMEDIATE_imul_[index] = INPUT_imul_[index];
    }
    else {
      throw Exception::NotImplemented(">>> ROL::OptimizationProblem::SetRiskAverseInequality: Either Constraint or SampleGenerator is NULL!");
    }
    // Set vector and bound constraint
    buildRiskVec(INPUT_sol_);
    buildRiskBnd(INPUT_bnd_);

    isInitialized_ = false;
  }

  void setAlmostSureInequality(const ROL::Ptr<SampleGenerator<Real> > &sampler,
                               const int index = 0) {
    initStochastic();
    ixsampler_[index] = sampler;
    if ( INPUT_icon_[index] != ROL::nullPtr && sampler != ROL::nullPtr ) {
      int nsamp = sampler->numMySamples();
      INTERMEDIATE_icon_[index]
        = ROL::makePtr<AlmostSureConstraint<Real>>(sampler, INPUT_icon_[index]);
      std::vector<ROL::Ptr<Vector<Real> > > imul(nsamp,ROL::nullPtr);
      for (int j = 0; j < nsamp; ++j) {
        imul[j] = INPUT_imul_[index]->clone();
        imul[j]->set(*INPUT_imul_[index]);
      }
      INTERMEDIATE_imul_[index]
        = ROL::makePtr<DualSimulatedVector<Real>>(imul, sampler->getBatchManager(), sampler);
      INTERMEDIATE_ibnd_[index]
        = ROL::makePtr<SimulatedBoundConstraint<Real>>(sampler, INPUT_ibnd_[index]);
    }
    else {
      throw Exception::NotImplemented(">>> ROL::OptimizationProblem::SetAlmostSureInequality: Either Constraint or SampleGenerator is NULL!");
    }
    // Set vector and bound constraint
    buildRiskVec(INPUT_sol_);
    buildRiskBnd(INPUT_bnd_);

    isInitialized_ = false;
  }

  void setStochasticInequality(std::vector<Teuchos::ParameterList> &parlist,
                               const std::vector<ROL::Ptr<SampleGenerator<Real> > > &xsampler,
                               const std::vector<ROL::Ptr<BatchManager<Real> > >    &cbman) {
    initStochastic();
    int nc = static_cast<int>(INPUT_icon_.size());
    if ( nc != static_cast<int>(xsampler.size()) || nc != static_cast<int>(cbman.size()) ) {
      throw Exception::NotImplemented(">>> ROL::OptimizationProblem::setStochasticInequality: Constraint vector and SampleGenerator vector are not the same size!");
    }
    for (int i = 0; i < nc; ++i) {
      if ( xsampler[i] != ROL::nullPtr ) {
        std::string type = parlist[i].sublist("SOL").get("Stochastic Component Type","Risk Neutral");
        if ( type == "Risk Neutral" ) {
          setRiskNeutralInequality(xsampler[i],cbman[i],i);
        }
        else if ( type == "Risk Averse" ||
                  type == "Deviation"   ||
                  type == "Error"       ||
                  type == "Regret"      ||
                  type == "Probability" ) {
          setRiskAverseInequality(parlist[i],xsampler[i],i);
        }
        else if ( type == "Almost Sure" ) {
          setAlmostSureInequality(xsampler[i],i);
        }
        else if ( type == "Mean Value" ) {
          setMeanValueInequality(xsampler[i],i);
        }
        else {
          throw Exception::NotImplemented(">>> ROL::OptimizationProblem::SetStochasticInequality: Invalid stochastic constraint type!");
        }
      }
      else {
        INTERMEDIATE_icon_[i] = INPUT_icon_[i];
        INTERMEDIATE_imul_[i] = INPUT_imul_[i];
        INTERMEDIATE_ibnd_[i] = INPUT_ibnd_[i];
      }
    }
    // Set vector and bound constraint
    buildRiskVec(INPUT_sol_);
    buildRiskBnd(INPUT_bnd_);

    isInitialized_ = false;
  }

  void setStochasticInequality(Teuchos::ParameterList &parlist,
                               const ROL::Ptr<SampleGenerator<Real> > &xsampler,
                               const ROL::Ptr<BatchManager<Real> > &cbman) {
    std::vector<Teuchos::ParameterList> cparlist(1,parlist);
    std::vector<ROL::Ptr<SampleGenerator<Real> > > cxsampler(1,xsampler);
    std::vector<ROL::Ptr<BatchManager<Real> > > ccbman(1,cbman);
    setStochasticInequality(cparlist,cxsampler,ccbman);
  }

  /** \brief Returns the statistic from the soluton vector.

      @param[in]    comp   is the component of the risk vector (0 for objective, 1 for inequality constraint)
      @param[in]    index  is the inequality constraint index
  */
  Real getSolutionStatistic(int comp = 0, int index = 0) {
    Real val(0);
    if (comp == 0) {
      try {
        val = dynamicPtrCast<StochasticObjective<Real>>(INTERMEDIATE_obj_)->computeStatistic(*INTERMEDIATE_sol_);
      }
      catch (std::exception &e) {
        throw Exception::NotImplemented(">>> ROL::OptimizationProblem::getSolutionStatistic: Objective does not have computeStatistic function!");
      }
    }
    else if (comp == 1) {
      int np = INTERMEDIATE_icon_.size();
      if (np <= index || index < 0) {
        throw Exception::NotImplemented(">>> ROL::OptimizationProblem::getSolutionStatistic: Index out of bounds!");
      }
      try {
        val = dynamicPtrCast<StochasticConstraint<Real>>(INTERMEDIATE_icon_[index])->computeStatistic(*INTERMEDIATE_sol_);
      }
      catch (std::exception &e) {
        throw Exception::NotImplemented(">>> ROL::OptimizationProblem::getSolutionStatistic: Constraint does not have computeStatistic function!");
      }
    }
    else {
      throw Exception::NotImplemented(">>> ROL::OptimizationProblem::getSolutionStatistic: Component must be either 0 or 1!");
    }
    return val;
  }

  void reset(void) {
    initialize(INTERMEDIATE_obj_,INTERMEDIATE_sol_,INTERMEDIATE_bnd_,
               INTERMEDIATE_econ_,INTERMEDIATE_emul_,
               INTERMEDIATE_icon_,INTERMEDIATE_imul_,INTERMEDIATE_ibnd_);
    conManager_->resetSlackVariables();
  }

  // Check derivatives, and consistency 
  void checkSolutionVector( Vector<Real> &x, // Optimization space
                            Vector<Real> &y, // Optimization space
                            Vector<Real> &u, // Optimization space
                            std::ostream &outStream = std::cout ) {
    initialize(INTERMEDIATE_obj_,INTERMEDIATE_sol_,INTERMEDIATE_bnd_,
               INTERMEDIATE_econ_,INTERMEDIATE_emul_,
               INTERMEDIATE_icon_,INTERMEDIATE_imul_,INTERMEDIATE_ibnd_);
    if (obj_ != ROL::nullPtr) {
      outStream << "\nPerforming OptimizationProblem diagnostics." << std::endl << std::endl;

      outStream << "Checking vector operations in optimization vector space X." << std::endl;
      x.checkVector(y,u,true,outStream);
    }
  }

  void checkObjective( Vector<Real> &x, // Optimization space
                       Vector<Real> &u, // Optimization space
                       Vector<Real> &v, // Optimization space
                       std::ostream &outStream = std::cout,
                       const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                       const int order = 1 ) {
    initialize(INTERMEDIATE_obj_,INTERMEDIATE_sol_,INTERMEDIATE_bnd_,
               INTERMEDIATE_econ_,INTERMEDIATE_emul_,
               INTERMEDIATE_icon_,INTERMEDIATE_imul_,INTERMEDIATE_ibnd_);
    if (obj_ != ROL::nullPtr) {
      outStream << "\nPerforming OptimizationProblem diagnostics." << std::endl << std::endl;

      outStream << "Checking objective function." << std::endl;
      obj_->checkGradient(x,v,true,outStream,numSteps,order); outStream << std::endl;
      obj_->checkHessVec(x,u,true,outStream,numSteps,order);  outStream << std::endl;
      obj_->checkHessSym(x,u,v,true,outStream);               outStream << std::endl;
    }
  }

  void checkMultiplierVector( Vector<Real> &w, // Dual constraint space
                              Vector<Real> &q, // Dual constraint space
                              Vector<Real> &l, // Dual constraint space
                              std::ostream &outStream = std::cout ) {
    initialize(INTERMEDIATE_obj_,INTERMEDIATE_sol_,INTERMEDIATE_bnd_,
               INTERMEDIATE_econ_,INTERMEDIATE_emul_,
               INTERMEDIATE_icon_,INTERMEDIATE_imul_,INTERMEDIATE_ibnd_);
    if(con_ != ROL::nullPtr) {
      outStream << "\nPerforming OptimizationProblem diagnostics." << std::endl << std::endl;

      outStream << "Checking vector operations in constraint multiplier space C*." << std::endl;
      l.checkVector(q,w,true,outStream);
    }
  }

  void checkConstraint( Vector<Real> &x, // Optimization space
                        Vector<Real> &u, // Optimization space
                        Vector<Real> &v, // Optimization space
                        Vector<Real> &c, // Constraint space
                        Vector<Real> &l, // Dual constraint space
                        std::ostream &outStream = std::cout,
                        const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                        const int order = 1 ) {
    initialize(INTERMEDIATE_obj_,INTERMEDIATE_sol_,INTERMEDIATE_bnd_,
               INTERMEDIATE_econ_,INTERMEDIATE_emul_,
               INTERMEDIATE_icon_,INTERMEDIATE_imul_,INTERMEDIATE_ibnd_);
    if(con_ != ROL::nullPtr) {
      outStream << "\nPerforming OptimizationProblem diagnostics." << std::endl << std::endl;

      outStream << "Checking equality constraint." << std::endl;
      con_->checkApplyJacobian(x,v,c,true,outStream,numSteps,order);         outStream << std::endl;
      con_->checkAdjointConsistencyJacobian(l,u,x,true,outStream);           outStream << std::endl;
      con_->checkApplyAdjointHessian(x,l,v,u,true,outStream,numSteps,order); outStream << std::endl;  
    }
  }

  // Check derivatives, and consistency 
  void check( std::ostream &outStream = std::cout,
              const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
              const int order = 1 ) {
    initialize(INTERMEDIATE_obj_,INTERMEDIATE_sol_,INTERMEDIATE_bnd_,
               INTERMEDIATE_econ_,INTERMEDIATE_emul_,
               INTERMEDIATE_icon_,INTERMEDIATE_imul_,INTERMEDIATE_ibnd_);

    ROL::Ptr<Vector<Real> > x, y, u, v;
    try {
      x = sol_->clone(); RandomizeVector(*x);
      y = sol_->clone(); RandomizeVector(*y);
      u = sol_->clone(); RandomizeVector(*u);
      v = sol_->clone(); RandomizeVector(*v);

      checkSolutionVector(*x,*y,*u,outStream);
      checkObjective(*x,*u,*v,outStream,numSteps,order);
    }
    catch (std::exception &e) {
//      throw Exception::NotImplemented(">>> ROL::OptimizationProblem::check: Elementwise is not implemented for optimization space vectors");
    }

    if(con_ != ROL::nullPtr) {
      ROL::Ptr<Vector<Real> > c, l, w, q;
      try {
        c = mul_->dual().clone(); RandomizeVector(*c);
        l = mul_->clone();        RandomizeVector(*l);
        w = mul_->clone();        RandomizeVector(*w);
        q = mul_->clone();        RandomizeVector(*q);   

        checkMultiplierVector(*w,*q,*l,outStream);
        checkConstraint(*x,*u,*v,*c,*l,outStream,numSteps,order);
      }
      catch (std::exception &e) {
        throw Exception::NotImplemented(">>> ROL::OptimizationProblem::check: Elementwise is not implemented for constraint space vectors");
      }
    }
  }

}; // class OptimizationProblem

}  // namespace ROL

#endif // ROL_OPTIMIZATIONPROBLEM_HPP
