// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_REDUCEDDYNAMICOBJECTIVE_HPP
#define ROL_REDUCEDDYNAMICOBJECTIVE_HPP

#include "ROL_Ptr.hpp"
#include "ROL_Sketch.hpp"
#include "ROL_Objective.hpp"
#include "ROL_DynamicObjective.hpp"
#include "ROL_DynamicConstraint.hpp"

#include <fstream>

/** @ingroup func_group
    \class ROL::ReducedDynamicObjective
    \brief Defines the reduced time-dependent objective function interface
           for simulation-based optimization.

    This objective function implements the implicitly defined objective
    function given by
    \f[
       F(z) := \sum_{n=1}^{N_t} f_n(u_{n-1}(z),u_n(z),z_n)
    \f]
    where \f$f_n:\mathcal{U}\times\mathcal{U}\times\mathcal{Z}\to\mathbb{R}\f$,
    and \f$u_n\in\mathcal{U}\f$ solves the system of equations
    \f[
       c_n(u_{n-1},u_n,z_n) = 0,\quad n=1,\ldots,N_t
    \f]
    with \f$u_0\f$ provided.

    Disclaimer: This is currently only set up for single step time integrators
    and piecewise-constant-in-time controls.
*/


namespace ROL {

template<typename Real> 
class ReducedDynamicObjective : public Objective<Real> {
  using size_type = typename std::vector<Real>::size_type;
private:
  // Problem data.
  const Ptr<DynamicObjective<Real>>  obj_;
  const Ptr<DynamicConstraint<Real>> con_;
  const Ptr<Vector<Real>>            u0_;
  const std::vector<TimeStamp<Real>> timeStamp_;
  const size_type                    Nt_;
  // General sketch information.
  const bool                         useSketch_;
  // State sketch information.
  size_type                          rankState_;
  Ptr<Sketch<Real>>                  stateSketch_;
  Ptr<Sketch<Real>>                  stateSketchCache_;
  // Adjoint sketch information.
  size_type                          rankAdjoint_;
  Ptr<Sketch<Real>>                  adjointSketch_;
  // State sensitivity sketch information.
  size_type                          rankStateSens_;
  Ptr<Sketch<Real>>                  stateSensSketch_;
  // Inexactness information.
  const bool                         useInexact_;
  const size_type                    updateFactor_;
  const size_type                    maxRank_;
  const bool                         syncHessRank_;
  const bool                         sumRankUpdate_;
  const bool                         useDefaultRankUpdate_;
  const Real                         a_, b_;                    
  const Real                         maxTol_;
  // Vector storage for intermediate computations.
  std::vector<Ptr<Vector<Real>>>     uhist_;             // State history.
  std::vector<Ptr<Vector<Real>>>     lhist_;             // Adjoint history.
  std::vector<Ptr<Vector<Real>>>     whist_;             // State sensitivity history.
  std::vector<Ptr<Vector<Real>>>     phist_;             // Adjoint sensitivity history.
  Ptr<Vector<Real>>                  cprimal_;           // Primal constraint vector.
  Ptr<Vector<Real>>                  crhs_;              // Primal constraint vector.
  Ptr<Vector<Real>>                  udual_;             // Dual state vector.
  Ptr<Vector<Real>>                  rhs_;               // Dual state vector.
  Ptr<Vector<Real>>                  zdual_;             // Dual control vector.
  Real                               val_, valCache_;    // Objective function value.
  // Flags.
  bool                               isValueComputed_;   // Whether objective has been evaluated.
  bool                               isStateComputed_;   // Whether state has been solved.
  bool                               isAdjointComputed_; // Whether adjoint has been solved.
  bool                               isValueCached_;     // Whether objective has been evaluated.
  bool                               isStateCached_;     // Whether state has been solved.
  bool                               useHessian_;        // Whether to use Hessian or FD.
  bool                               useSymHess_;        // Whether to use symmetric Hessian approximation.
  // Output information
  Ptr<std::ostream>                  stream_;
  const bool                         print_;
  const int                          freq_;

  PartitionedVector<Real>& partition ( Vector<Real>& x ) const {
    return static_cast<PartitionedVector<Real>&>(x);
  }

  const PartitionedVector<Real>& partition ( const Vector<Real>& x ) const {
    return static_cast<const PartitionedVector<Real>&>(x);
  }

  void throwError(const int err, const std::string &sfunc,
                  const std::string &func, const int line) const {
    if (err != 0) {
      std::stringstream out;
      out << ">>> ROL::ReducedDynamicObjective::" << func << " (Line " << line << "): ";
      if (err == 1)
        out << "Reconstruction has already been called!";
      else if (err == 2)
        out << "Input column index exceeds total number of columns!";
      else if (err == 3)
        out << "Reconstruction failed to compute domain-space basis!";
      else if (err == 4)
        out << "Reconstruction failed to compute range-space basis!";
      else if (err == 5)
        out << "Reconstruction failed to generate core matrix!";
      throw Exception::NotImplemented(out.str());
    }
  }

public:
  ReducedDynamicObjective(const Ptr<DynamicObjective<Real>>  &obj,
                          const Ptr<DynamicConstraint<Real>> &con,
                          const Ptr<Vector<Real>>            &u0,
                          const Ptr<Vector<Real>>            &zvec,
                          const Ptr<Vector<Real>>            &cvec,
                          const std::vector<TimeStamp<Real>> &timeStamp,
                          ROL::ParameterList                 &pl,
                          const Ptr<std::ostream>            &stream = nullPtr)
    : obj_                  ( obj ),                                           // Dynamic objective function.
      con_                  ( con ),                                           // Dynamic constraint function.
      u0_                   ( u0 ),                                            // Initial condition.
      timeStamp_            ( timeStamp ),                                     // Vector of time stamps.
      Nt_                   ( timeStamp.size() ),                              // Number of time intervals.
      useSketch_            ( pl.get("Use Sketching", false) ),                // Use state sketch if true.
      rankState_            ( pl.get("State Rank", 10) ),                      // Rank of state sketch.
      stateSketch_          ( nullPtr ),                                       // State sketch object.
      rankAdjoint_          ( pl.get("Adjoint Rank", 10) ),                    // Rank of adjoint sketch.
      adjointSketch_        ( nullPtr ),                                       // Adjoint sketch object.
      rankStateSens_        ( pl.get("State Sensitivity Rank", 10) ),          // Rank of state sensitivity sketch.
      stateSensSketch_      ( nullPtr ),                                       // State sensitivity sketch object.
      useInexact_           ( pl.get("Adaptive Rank", false) ),                // Update rank adaptively.
      updateFactor_         ( pl.get("Rank Update Factor", 2) ),               // Rank update factor.
      maxRank_              ( pl.get("Maximum Rank", 100) ),                   // Maximum rank.
      syncHessRank_         ( pl.get("Sync Hessian Rank", useInexact_) ),      // Sync rank for Hessian storage.
      sumRankUpdate_        ( pl.get("Additive Rank Update", true) ),          // Use additive rank update, otherwise use multiplicative
      useDefaultRankUpdate_ ( pl.get("Use Basic Rank Update", true) ),         // Use basic additive/multiplicative rank update
      a_                    ( pl.get("Log Rank Update Slope", 1.0) ),          // Slope of semilogy tail energy
      b_                    ( pl.get("Log Rank Update Shift", 1.0) ),          // Shift of semilogy tail energy
      maxTol_               ( pl.get("Maximum Tolerance", ROL_INF<Real>()) ),  // Maximum rank update tolerance
      isValueComputed_      ( false ),                                         // Flag indicating whether value has been computed.
      isStateComputed_      ( false ),                                         // Flag indicating whether state has been computed.
      isAdjointComputed_    ( false ),                                         // Flag indicating whether adjoint has been computed.
      isValueCached_        ( false ),                                         // Flag indicating whether value has been computed.
      isStateCached_        ( false ),                                         // Flag indicating whether state has been computed.
      useHessian_           ( pl.get("Use Hessian", true) ),                   // Flag indicating whether to use the Hessian.
      useSymHess_           ( pl.get("Use Only Sketched Sensitivity", true) ), // Flag indicating whether to use symmetric sketched Hessian.
      stream_               ( stream ),                                        // Output stream to print sketch information.
      print_                ( pl.get("Print Optimization Vector", false) ),    // Print control vector to file
      freq_                 ( pl.get("Output Frequency", 5) ) {                // Print frequency for control vector
    uhist_.clear(); lhist_.clear(); whist_.clear(); phist_.clear();
    if (useSketch_) { // Only maintain a sketch of the state time history
      Real orthTol   = pl.get("Orthogonality Tolerance", 1e2*ROL_EPSILON<Real>());
      int  orthIt    = pl.get("Reorthogonalization Iterations", 5);
      bool trunc     = pl.get("Truncate Approximation", false);
      if (syncHessRank_) {
        rankAdjoint_   = rankState_;
        rankStateSens_ = rankState_;
      }
      unsigned dseed = pl.get("State Domain Seed",0);
      unsigned rseed = pl.get("State Range Seed",0);
      stateSketch_ = makePtr<Sketch<Real>>(*u0_,static_cast<int>(Nt_)-1,
        rankState_,orthTol,orthIt,trunc,dseed,rseed);
      stateSketchCache_ = makePtr<Sketch<Real>>(*u0_,static_cast<int>(Nt_)-1,
        rankState_,orthTol,orthIt,trunc,dseed,rseed);
      uhist_.push_back(u0_->clone());
      uhist_.push_back(u0_->clone());
      lhist_.push_back(cvec->dual().clone());
      dseed = pl.get("Adjoint Domain Seed",0);
      rseed = pl.get("Adjoint Range Seed",0);
      adjointSketch_ = makePtr<Sketch<Real>>(*u0_,static_cast<int>(Nt_)-1,rankAdjoint_,
        orthTol,orthIt,trunc,dseed,rseed);
      if (useHessian_) {
        dseed = pl.get("State Sensitivity Domain Seed",0);
        rseed = pl.get("State Sensitivity Range Seed",0);
        stateSensSketch_ = makePtr<Sketch<Real>>(*u0_,static_cast<int>(Nt_)-1,
          rankStateSens_,orthTol,orthIt,trunc,dseed,rseed);
        whist_.push_back(u0_->clone());
        whist_.push_back(u0_->clone());
        phist_.push_back(cvec->dual().clone());
      }
    }
    else {            // Store entire state time history
      for (size_type i = 0; i < Nt_; ++i) {
        uhist_.push_back(u0_->clone());
        lhist_.push_back(cvec->dual().clone());
        if (useHessian_) {
          whist_.push_back(u0_->clone());
          phist_.push_back(cvec->dual().clone());
        }
      }
    }
    cprimal_ = cvec->clone();
    crhs_    = cvec->clone();
    udual_   = u0_->dual().clone();
    rhs_     = u0_->dual().clone();
    zdual_   = zvec->dual().clone();
  }

  Ptr<Vector<Real>> makeDynamicVector(const Vector<Real> &x) const {
    return ROL::PartitionedVector<Real>::create(x, Nt_);
  }

  size_type getStateRank() const { return rankState_; } 
  size_type getAdjointRank() const { return rankAdjoint_; } 
  size_type getStateSensitivityRank() const { return rankStateSens_; } 

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    if (useSketch_) {
      stateSketch_->reset(true);
      if (flag == true) {
        adjointSketch_->reset(true);
        if (useHessian_) stateSensSketch_->reset(true);
      }
    }
    for (size_type i = 0; i < uhist_.size(); ++i) uhist_[i]->zero();
    val_ = static_cast<Real>(0);
    isValueComputed_ = false;
    isStateComputed_ = false;
    if (flag == true) isAdjointComputed_ = false;

    if (iter >= 0 && iter%freq_==0 && print_) {
      std::stringstream name;
      name << "optvector." << iter << ".txt";
      std::ofstream file;
      file.open(name.str());
      x.print(file);
      file.close();
    }
  }

  void update(const Vector<Real> &x, UpdateType type, int iter = -1) {
    // This just resets all storage independent of update type.
    // If sketching is not in use, we could be smarter with the
    // storage.  When sketching is in use, you may want to reset
    // everything and re-randomize the sketch object.  It is
    // unclear to me (Drew) if there is a benefit to saving the
    // previous computations in this case.
    if (useSketch_) {
      switch(type) {
        case UpdateType::Initial:
        {
          for (size_type i = 0; i < uhist_.size(); ++i) uhist_[i]->zero();
          val_               = static_cast<Real>(0);
          isValueComputed_   = false;
          isStateComputed_   = false;
          isAdjointComputed_ = false;
          stateSketch_->reset(true);
          stateSketchCache_->reset(true);
          adjointSketch_->reset(true);
          if (useHessian_) stateSensSketch_->reset(true);
          break;
        }
        case UpdateType::Trial:
        {
          for (size_type i = 0; i < uhist_.size(); ++i) uhist_[i]->zero();
          valCache_          = val_;
          isValueCached_     = isValueComputed_;
          isStateCached_     = isStateComputed_;
          val_               = static_cast<Real>(0);
          isValueComputed_   = false;
          isStateComputed_   = false;
          auto tmp           = stateSketch_;
          stateSketch_       = stateSketchCache_;
          stateSketchCache_  = tmp;
          stateSketch_->reset(true);
          break;
        }
        case UpdateType::Accept:
        {
          isAdjointComputed_ = false;
          adjointSketch_->reset(true);
          if (useHessian_) stateSensSketch_->reset(true);
          break;
        }
        case UpdateType::Revert:
        {
          for (size_type i = 0; i < uhist_.size(); ++i) uhist_[i]->zero();
          val_               = valCache_;
          isValueComputed_   = isValueCached_;
          isStateComputed_   = isStateCached_;
          auto tmp           = stateSketchCache_;
          stateSketchCache_  = stateSketch_;
          stateSketch_       = tmp;
          if (useHessian_) stateSensSketch_->reset(true);
          break;
        }
        case UpdateType::Temp:
        {
          for (size_type i = 0; i < uhist_.size(); ++i) uhist_[i]->zero();
          valCache_          = val_;
          isValueCached_     = isValueComputed_;
          isStateCached_     = isStateComputed_;
          val_               = static_cast<Real>(0);
          isValueComputed_   = false;
          isStateComputed_   = false;
          isAdjointComputed_ = false;
          auto tmp           = stateSketch_;
          stateSketch_       = stateSketchCache_;
          stateSketchCache_  = tmp;
          stateSketch_->reset(true);
          adjointSketch_->reset(true);
          if (useHessian_) stateSensSketch_->reset(true);
          break;
        }
      }
    }
    else {
      switch(type) {
        case UpdateType::Initial:
        case UpdateType::Accept:
        case UpdateType::Revert:
        case UpdateType::Trial:
        case UpdateType::Temp:
        {
          for (size_type i = 0; i < uhist_.size(); ++i) uhist_[i]->zero();
          val_ = static_cast<Real>(0);
          isValueComputed_   = false;
          isStateComputed_   = false;
          isAdjointComputed_ = false;
          break;
        }
      }
    }

    if (iter >= 0 && iter%freq_==0 && print_) {
      std::stringstream name;
      name << "optvector." << iter << ".txt";
      std::ofstream file;
      file.open(name.str());
      x.print(file);
      file.close();
    }
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    if (!isValueComputed_) {
      int eflag(0);
      const Real one(1);
      const PartitionedVector<Real> &xp = partition(x);
      // Set initial condition
      uhist_[0]->set(*u0_);
      if (useSketch_) {
        stateSketch_->update();
        uhist_[1]->set(*u0_);
      }
      // Run time stepper
      Real valk(0);
      size_type index;
      for (size_type k = 1; k < Nt_; ++k) {
        index = (useSketch_ ? 1 : k);
        if (!useSketch_) uhist_[index]->set(*uhist_[index-1]);
        // Update dynamic constraint
        con_->update_uo(*uhist_[index-1], timeStamp_[k]);
        con_->update_z(*xp.get(k), timeStamp_[k]);
        // Solve state on current time interval
        con_->solve(*cprimal_, *uhist_[index-1], *uhist_[index], *xp.get(k), timeStamp_[k]);
        // Update dynamic objective
        obj_->update(*uhist_[index-1], *uhist_[index], *xp.get(k), timeStamp_[k]);
        // Compute objective function value on current time interval
        valk  = obj_->value(*uhist_[index-1], *uhist_[index], *xp.get(k), timeStamp_[k]);
        // Update total objective function value
        val_ += valk;
        // Sketch state
        if (useSketch_) {
          eflag = stateSketch_->advance(one,*uhist_[1],static_cast<int>(k)-1,one);
          throwError(eflag,"advance","value",315);
          uhist_[0]->set(*uhist_[1]);
        }
      }
      isValueComputed_ = true;
      isStateComputed_ = true;
    }
    return val_;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    int eflag(0);
    PartitionedVector<Real>       &gp = partition(g);
    gp.get(0)->zero(); // zero for the nonexistant zeroth control interval
    const PartitionedVector<Real> &xp = partition(x);
    const Real one(1);
    size_type uindex = (useSketch_ ? 1 : Nt_-1);
    size_type lindex = (useSketch_ ? 0 : Nt_-1);
    // Must first compute the value
    solveState(x);
    if (useSketch_ && useInexact_) {
      if (stream_ != nullPtr) {
        *stream_ << std::string(80,'=') << std::endl; 
        *stream_ << "  ROL::ReducedDynamicObjective::gradient" << std::endl;
      }
      tol = updateSketch(x,tol);
      if (stream_ != nullPtr) {
        *stream_ << "    State Rank for Gradient Computation: " << rankState_ << std::endl;
        *stream_ << "    Residual Norm:                       " << tol << std::endl;
        *stream_ << std::string(80,'=') << std::endl; 
      }
    }
    // Recover state from sketch
    if (useSketch_) {
      uhist_[1]->set(*uhist_[0]);
      eflag = stateSketch_->reconstruct(*uhist_[0],static_cast<int>(Nt_)-3);
      throwError(eflag,"reconstruct","gradient",351);
      if (isAdjointComputed_) {
        eflag = adjointSketch_->reconstruct(*lhist_[0],static_cast<int>(Nt_)-2);
        throwError(eflag,"reconstruct","gradient",354);
      }
    }
    // Update dynamic constraint and objective
    con_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(Nt_-1), timeStamp_[Nt_-1]);
    obj_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(Nt_-1), timeStamp_[Nt_-1]);
    // Solve for terminal condition
    if (!isAdjointComputed_) {
      setTerminalCondition(*lhist_[lindex],
                           *uhist_[uindex-1], *uhist_[uindex],
                           *xp.get(Nt_-1),    timeStamp_[Nt_-1]);
      if (useSketch_) {
        eflag = adjointSketch_->advance(one,*lhist_[0],static_cast<int>(Nt_)-2,one);
        throwError(eflag,"advance","gradient",367);
      }
    }
    // Update gradient on terminal interval
    updateGradient(*gp.get(Nt_-1),    *lhist_[lindex],
                   *uhist_[uindex-1], *uhist_[uindex],
                   *xp.get(Nt_-1),    timeStamp_[Nt_-1]);
    // Run reverse time stepper
    for (size_type k = Nt_-2; k > 0; --k) {
      if (!isAdjointComputed_) {
        // Compute k+1 component of rhs
        computeAdjointRHS(*rhs_,             *lhist_[lindex],
                          *uhist_[uindex-1], *uhist_[uindex],
                          *xp.get(k+1),      timeStamp_[k+1]);
      }
      uindex = (useSketch_ ? 1 : k);
      lindex = (useSketch_ ? 0 : k);
      // Recover state from sketch
      if (useSketch_) {
        uhist_[1]->set(*uhist_[0]);
        if (k==1) {
          uhist_[0]->set(*u0_);
        }
        else {
          eflag = stateSketch_->reconstruct(*uhist_[0],static_cast<int>(k)-2);
          throwError(eflag,"reconstruct","gradient",392);
        }
        if (isAdjointComputed_) {
          eflag = adjointSketch_->reconstruct(*lhist_[0],static_cast<int>(k)-1);
          throwError(eflag,"reconstruct","gradient",396);
        }
      }
      // Update dynamic constraint and objective
      con_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(k), timeStamp_[k]);
      obj_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(k), timeStamp_[k]);
      // Solve for adjoint on interval k
      if (!isAdjointComputed_) {
        advanceAdjoint(*lhist_[lindex],   *rhs_,
                       *uhist_[uindex-1], *uhist_[uindex],
                       *xp.get(k),        timeStamp_[k]);
        if (useSketch_) {
          eflag = adjointSketch_->advance(one,*lhist_[0],static_cast<int>(k)-1,one);
          throwError(eflag,"advance","gradient",367);
        }
      }
      // Update gradient on interval k
      updateGradient(*gp.get(k),        *lhist_[lindex],
                     *uhist_[uindex-1], *uhist_[uindex],
                     *xp.get(k),        timeStamp_[k]);
    }
    isAdjointComputed_ = true;
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    int eflag(0);
    if (useHessian_) {
      // Must first solve the state and adjoint equations
      solveState(x);
      solveAdjoint(x);
      // Now compute Hessian
      const Real one(1);
      const PartitionedVector<Real> &xp  = partition(x);
      const PartitionedVector<Real> &vp  = partition(v);
      PartitionedVector<Real>       &hvp = partition(hv);
      hvp.get(0)->zero(); // zero for the nonexistant zeroth control interval
      // Compute state sensitivity
      whist_[0]->zero();
      if (useSketch_) {
        stateSensSketch_->reset(false);
        //stateSensSketch_->advance(one,*whist_[0],0,one);
        uhist_[0]->set(*u0_);
      }
      size_type uindex, lindex;
      for (size_type k = 1; k < Nt_; ++k) {
        uindex = (useSketch_ ? 1 : k);
        lindex = (useSketch_ ? 0 : k);
        // Reconstruct sketched state
        if (useSketch_) {
          eflag = stateSketch_->reconstruct(*uhist_[1],static_cast<int>(k)-1);
          throwError(eflag,"reconstruct","hessVec",446);
          eflag = adjointSketch_->reconstruct(*lhist_[0],static_cast<int>(k)-1);
          throwError(eflag,"reconstruct","hessVec",448);
        }
        // Update dynamic constraint and objective
        con_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(k), timeStamp_[k]);
        obj_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(k), timeStamp_[k]);
        // Advance state sensitivity on current time interval
        advanceStateSens(*whist_[uindex],
                         *vp.get(k),        *whist_[uindex-1],
                         *uhist_[uindex-1], *uhist_[uindex],
                         *xp.get(k),        timeStamp_[k]);
        // Set control only Hessian
        computeControlHessLag(*hvp.get(k),
                              *vp.get(k),        *lhist_[lindex],
                              *uhist_[uindex-1], *uhist_[uindex],
                              *xp.get(k),        timeStamp_[k]);
        if (!useSymHess_) {
          // Add mixed derivative Hessian
          addMixedHessLag(*hvp.get(k),       *lhist_[lindex],
                          *whist_[uindex-1], *whist_[uindex],
                          *uhist_[uindex-1], *uhist_[uindex],
                          *xp.get(k),        timeStamp_[k]);
        }
        // Sketch state sensitivity
        if (useSketch_) {
          eflag = stateSensSketch_->advance(one,*whist_[1],static_cast<int>(k)-1,one);
          throwError(eflag,"advance","hessVec",473);
          whist_[0]->set(*whist_[1]);
          uhist_[0]->set(*uhist_[1]);
        }
      }

      // Compute adjoint sensitivity
      uindex = (useSketch_ ? 1 : Nt_-1);
      lindex = (useSketch_ ? 0 : Nt_-1);
      if (useSketch_) {
        // Recover terminal state
        uhist_[1]->set(*uhist_[0]);
        eflag = stateSketch_->reconstruct(*uhist_[0],static_cast<int>(Nt_)-3);
        throwError(eflag,"reconstruct","hessVec",486);
        // Recover terminal adjoint
        eflag = adjointSketch_->reconstruct(*lhist_[0],static_cast<int>(Nt_)-2);
        throwError(eflag,"reconstruct","hessVec",489);
        // Recover terminal state sensitivity
        whist_[1]->set(*whist_[0]);
        eflag = stateSensSketch_->reconstruct(*whist_[0],static_cast<int>(Nt_)-3);
        throwError(eflag,"reconstruct","hessVec",493);
      }
      // Update dynamic constraint and objective
      con_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(Nt_-1), timeStamp_[Nt_-1]);
      obj_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(Nt_-1), timeStamp_[Nt_-1]);
      // Solve for terminal condition
      setTerminalConditionHess(*phist_[lindex],
                               *vp.get(Nt_-1),    *lhist_[lindex],
                               *whist_[uindex-1], *whist_[uindex],
                               *uhist_[uindex-1], *uhist_[uindex],
                               *xp.get(Nt_-1),    timeStamp_[Nt_-1]);
      if (useSymHess_) {
        // Add mixed derivative Hessian
        addMixedHessLag(*hvp.get(Nt_-1),   *lhist_[lindex],
                        *whist_[uindex-1], *whist_[uindex],
                        *uhist_[uindex-1], *uhist_[uindex],
                        *xp.get(Nt_-1),    timeStamp_[Nt_-1]);
      }
      // Add adjoint sensitivity to Hessian
      addAdjointSens(*hvp.get(Nt_-1),   *phist_[lindex],
                     *uhist_[uindex-1], *uhist_[uindex],
                     *xp.get(Nt_-1),    timeStamp_[Nt_-1]);
      // Run reverse time stepper
      for (size_type k = Nt_-2; k > 0; --k) {
        // Compute new components of rhs
        computeNewStateHessLag(*rhs_,             *lhist_[lindex],
                               *whist_[uindex-1], *whist_[uindex],
                               *uhist_[uindex-1], *uhist_[uindex],
                               *xp.get(k+1),      timeStamp_[k+1],
                               false);
        computeNewMixedHessLag(*rhs_,
                               *vp.get(k+1),      *lhist_[lindex],
                               *uhist_[uindex-1], *uhist_[uindex],
                               *xp.get(k+1),      timeStamp_[k+1],
                               true);
        computeNewStateJacobian(*rhs_,             *phist_[lindex],
                                *uhist_[uindex-1], *uhist_[uindex],
                                *xp.get(k+1),      timeStamp_[k+1],
                                true);
        // Recover state, adjoint and state sensitivity from sketch
        if (useSketch_) {
          uhist_[1]->set(*uhist_[0]);
          whist_[1]->set(*whist_[0]);
          if (k==1) {
            uhist_[0]->set(*u0_);
            whist_[0]->zero();
          }
          else {
            eflag = stateSketch_->reconstruct(*uhist_[0],static_cast<int>(k)-2);
            throwError(eflag,"reconstruct","hessVec",542);
            eflag = stateSensSketch_->reconstruct(*whist_[0],static_cast<int>(k)-2);
            throwError(eflag,"reconstruct","hessVec",544);
          }
          eflag = adjointSketch_->reconstruct(*lhist_[0],static_cast<int>(k)-1);
          throwError(eflag,"reconstruct","hessVec",547);
        }
        uindex = (useSketch_ ? 1 : k);
        lindex = (useSketch_ ? 0 : k);
        // Update dynamic constraint and objective
        con_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(k), timeStamp_[k]);
        obj_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(k), timeStamp_[k]);
        // Compute old components of rhs
        computeOldStateHessLag(*rhs_,             *lhist_[lindex],
                               *whist_[uindex-1], *whist_[uindex],
                               *uhist_[uindex-1], *uhist_[uindex],
                               *xp.get(k),        timeStamp_[k],
                               true);
        computeOldMixedHessLag(*rhs_,
                               *vp.get(k),        *lhist_[lindex],
                               *uhist_[uindex-1], *uhist_[uindex],
                               *xp.get(k),        timeStamp_[k],
                               true);
        // Solve for adjoint on interval k
        advanceAdjointSens(*phist_[lindex],   *rhs_,
                           *uhist_[uindex-1], *uhist_[uindex],
                           *xp.get(k),        timeStamp_[k]);
        if (useSymHess_) {
          // Add mixed derivative Hessian
          addMixedHessLag(*hvp.get(k),       *lhist_[lindex],
                          *whist_[uindex-1], *whist_[uindex],
                          *uhist_[uindex-1], *uhist_[uindex],
                          *xp.get(k),        timeStamp_[k]);
        }
        // Add adjoint sensitivity to Hessian
        addAdjointSens(*hvp.get(k),       *phist_[lindex],
                       *uhist_[uindex-1], *uhist_[uindex],
                       *xp.get(k),        timeStamp_[k]);
      }
    }
    else {
      Objective<Real>::hessVec(hv,v,x,tol);
    }
  }

private:
  /***************************************************************************/
  /************ Method to solve state equation *******************************/
  /***************************************************************************/
  Real solveState(const Vector<Real> &x) {
    Real cnorm(0);
    if (!isStateComputed_) {
      int eflag(0);
      const Real one(1);
      const PartitionedVector<Real> &xp = partition(x);
      // Set initial condition
      uhist_[0]->set(*u0_);
      if (useSketch_) { // Set initial guess for solve
        stateSketch_->update();
        uhist_[1]->set(*u0_);
      }
      // Run time stepper
      size_type index;
      for (size_type k = 1; k < Nt_; ++k) {
        index = (useSketch_ ? 1 : k);
        // Set initial guess for solve
        if (!useSketch_) uhist_[index]->set(*uhist_[index-1]);
        // Solve state on current time interval
        con_->update_uo(*uhist_[index-1], timeStamp_[k]);
        con_->update_z(*xp.get(k), timeStamp_[k]);
        con_->solve(*cprimal_, *uhist_[index-1], *uhist_[index], *xp.get(k), timeStamp_[k]);
        cnorm = std::max(cnorm,cprimal_->norm()); 
        // Sketch state
        if (useSketch_) {
          eflag = stateSketch_->advance(one, *uhist_[1], static_cast<int>(k)-1, one);
          throwError(eflag,"advance","solveState",574);
          uhist_[0]->set(*uhist_[1]);
        }
      }
      isStateComputed_ = true;
    }
    return cnorm;
  }

  Real updateSketch(const Vector<Real> &x, const Real tol) {
    int eflag(0);
    const PartitionedVector<Real> &xp = partition(x);
    Real err(0), err2(0), serr(0), cdot(0), dt(0); //, cnorm(0)
    Real tol0 = std::min(maxTol_,tol);
    bool flag = true;
    while (flag) {
      err  = static_cast<Real>(0);
      err2 = err;
      uhist_[0]->set(*u0_);
      for (size_type k = 1; k < Nt_; ++k) {
        eflag = stateSketch_->reconstruct(*uhist_[1],static_cast<int>(k)-1);
        throwError(eflag,"reconstruct","updateSketch",592);
        con_->update(*uhist_[0], *uhist_[1], *xp.get(k), timeStamp_[k]);
        con_->value(*cprimal_, *uhist_[0], *uhist_[1], *xp.get(k), timeStamp_[k]);
        /**** Linf norm for residual error ****/
        //cnorm = cprimal_->norm();
        //err   = (cnorm > err ? cnorm : err);
        /**** L2 norm for residual error ****/
        cdot  = cprimal_->dot(*cprimal_);
        dt    = timeStamp_[k].t[1]-timeStamp_[k].t[0];
        err2 += cdot / dt;
        err   = std::sqrt(err2);
        if (err > tol0) break;
        uhist_[0]->set(*uhist_[1]);
      }
      if (stream_ != nullPtr) {
        *stream_ << "      *** State Rank:                    " << rankState_ << std::endl;
        *stream_ << "      *** Required Tolerance:            " << tol0 << std::endl;
        *stream_ << "      *** Residual Norm:                 " << err << std::endl;
      }
      if (err > tol0) {
        rankState_ = (sumRankUpdate_ ? rankState_ + updateFactor_ : rankState_ * updateFactor_); 
        if (!useDefaultRankUpdate_)
          rankState_ = std::max(rankState_,static_cast<size_type>(std::ceil((b_-std::log(tol0))/a_)));
        //Real a(0.1838), b(3.1451); // Navier-Stokes
        //Real a(2.6125), b(2.4841); // Semilinear
        //rankState_  = std::max(rankState_+2,static_cast<size_t>(std::ceil((b-std::log(tol0))/a)));
        //rankState_ *= updateFactor_; // Perhaps there is a better update strategy
        rankState_  = (maxRank_ < rankState_ ? maxRank_ : rankState_);
        stateSketch_->setRank(rankState_);
        if (syncHessRank_) {
          rankAdjoint_   = rankState_;
          rankStateSens_ = rankState_;
          adjointSketch_->setRank(rankAdjoint_);
          stateSensSketch_->setRank(rankStateSens_);
        }
        isStateComputed_   = false;
        isAdjointComputed_ = false;
        serr = solveState(x);
        if (stream_ != nullPtr)
          *stream_ << "      *** Maximum Solver Error:          " << serr << std::endl;
      }
      else {
        flag = false;
        break;
      }
    }
    return err;
  }

  /***************************************************************************/
  /************ Methods to solve adjoint equation ****************************/
  /***************************************************************************/
  void setTerminalCondition(Vector<Real> &l,
                      const Vector<Real> &uold, const Vector<Real> &unew,
                      const Vector<Real> &z,    const TimeStamp<Real> &ts) {
    obj_->gradient_un(*udual_, uold, unew, z, ts);
    con_->applyInverseAdjointJacobian_un(l, *udual_, uold, unew, z, ts);
  }

  void computeAdjointRHS(Vector<Real> &rhs,  const Vector<Real>    &l,
                   const Vector<Real> &uold, const Vector<Real>    &unew,
                   const Vector<Real> &z,    const TimeStamp<Real> &ts) {
    const Real one(1);
    obj_->gradient_uo(rhs, uold, unew, z, ts);
    con_->applyAdjointJacobian_uo(*udual_, l, uold, unew, z, ts);
    rhs.axpy(-one,*udual_);
  }

  void advanceAdjoint(Vector<Real> &l,          Vector<Real>    &rhs,
                const Vector<Real> &uold, const Vector<Real>    &unew,
                const Vector<Real> &z,    const TimeStamp<Real> &ts) {
    obj_->gradient_un(*udual_, uold, unew, z, ts);
    rhs.plus(*udual_);
    con_->applyInverseAdjointJacobian_un(l, rhs, uold, unew, z, ts);
  }

  void solveAdjoint(const Vector<Real> &x) {
    int eflag(0);
    if (!isAdjointComputed_) {
      const PartitionedVector<Real> &xp = partition(x);
      const Real one(1);
      size_type uindex = (useSketch_ ? 1 : Nt_-1);
      size_type lindex = (useSketch_ ? 0 : Nt_-1);
      // Must first compute solve the state equation
      solveState(x);
      // Recover state from sketch
      if (useSketch_) {
        uhist_[1]->set(*uhist_[0]);
        eflag = stateSketch_->reconstruct(*uhist_[0],static_cast<int>(Nt_)-3);
        throwError(eflag,"reconstruct","solveAdjoint",672); 
      }
      // Update dynamic constraint and objective
      con_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(Nt_-1), timeStamp_[Nt_-1]);
      obj_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(Nt_-1), timeStamp_[Nt_-1]);
      // Solve for terminal condition
      setTerminalCondition(*lhist_[lindex],
                           *uhist_[uindex-1], *uhist_[uindex],
                           *xp.get(Nt_-1),    timeStamp_[Nt_-1]);
      if (useSketch_) {
        eflag = adjointSketch_->advance(one,*lhist_[lindex],static_cast<int>(Nt_)-2,one);
        throwError(eflag,"advance","solveAdjoint",683);
      }
      // Run reverse time stepper
      for (size_type k = Nt_-2; k > 0; --k) {
        // Compute k+1 component of rhs
        computeAdjointRHS(*rhs_,             *lhist_[lindex],
                          *uhist_[uindex-1], *uhist_[uindex],
                          *xp.get(k+1),      timeStamp_[k+1]);
        // Recover state from sketch
        if (useSketch_) {
          uhist_[1]->set(*uhist_[0]);
          if (k==1) {
            uhist_[0]->set(*u0_);
          }
          else {
            eflag = stateSketch_->reconstruct(*uhist_[0],static_cast<int>(k)-2);
            throwError(eflag,"reconstruct","solveAdjoint",699); 
          }
        }
        uindex = (useSketch_ ? 1 : k);
        lindex = (useSketch_ ? 0 : k);
        // Update dynamic constraint and objective
        con_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(k), timeStamp_[k]);
        obj_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(k), timeStamp_[k]);
        // Solve for adjoint on interval k
        advanceAdjoint(*lhist_[lindex],   *rhs_,
                       *uhist_[uindex-1], *uhist_[uindex],
                       *xp.get(k),        timeStamp_[k]);
        if (useSketch_) {
          eflag = adjointSketch_->advance(one,*lhist_[lindex],static_cast<int>(k)-1,one);
          throwError(eflag,"advance","solveAdjoint",713);
        }
      }
      isAdjointComputed_ = true;
    }
  }

  /***************************************************************************/
  /************ Method for gradient computation ******************************/
  /***************************************************************************/
  void updateGradient(Vector<Real> &g,    const Vector<Real> &l,
                const Vector<Real> &uold, const Vector<Real> &unew,
                const Vector<Real> &z,    const TimeStamp<Real> &ts) {
    const Real one(1);
    obj_->gradient_z(g, uold, unew, z, ts);
    con_->applyAdjointJacobian_z(*zdual_, l, uold, unew, z, ts);
    g.axpy(-one,*zdual_);
  }

  /***************************************************************************/
  /************ Method to solve state sensitivity equation *******************/
  /***************************************************************************/
  void advanceStateSens(Vector<Real>    &wnew, const Vector<Real> &v,
                  const Vector<Real>    &wold, const Vector<Real> &uold,
                  const Vector<Real>    &unew, const Vector<Real> &z,
                  const TimeStamp<Real> &ts) {
    const Real one(1);
    con_->applyJacobian_z(*crhs_, v, uold, unew, z, ts);
    con_->applyJacobian_uo(*cprimal_, wold, uold, unew, z, ts);
    crhs_->axpy(-one, *cprimal_);
    con_->applyInverseJacobian_un(wnew, *crhs_, uold, unew, z, ts);
  }

  /***************************************************************************/
  /************ Methods to solve adjoint sensitivity equation ****************/
  /***************************************************************************/
  void setTerminalConditionHess(Vector<Real> &p,
                          const Vector<Real> &v,    const Vector<Real>    &l,
                          const Vector<Real> &wold, const Vector<Real>    &wnew,
                          const Vector<Real> &uold, const Vector<Real>    &unew,
                          const Vector<Real> &z,    const TimeStamp<Real> &ts) {
    const Real one(1);
    // Mixed derivative rhs term
    con_->applyAdjointHessian_z_un(*rhs_, l, v, uold, unew, z, ts);
    obj_->hessVec_un_z(*udual_, v, uold, unew, z, ts);
    rhs_->axpy(-one, *udual_);
    // State derivative rhs term
    con_->applyAdjointHessian_un_un(*udual_, l, wnew, uold, unew, z, ts);
    rhs_->axpy(-one, *udual_);
    obj_->hessVec_un_un(*udual_, wnew, uold, unew, z, ts);
    rhs_->plus(*udual_);
    con_->applyAdjointHessian_uo_un(*udual_, l, wold, uold, unew, z, ts);
    rhs_->axpy(-one, *udual_);
    obj_->hessVec_un_uo(*udual_, wold, uold, unew, z, ts);
    rhs_->plus(*udual_);
    // Invert adjoint Jacobian
    con_->applyInverseAdjointJacobian_un(p, *rhs_, uold, unew, z, ts);
  }

  void computeOldStateHessLag(Vector<Real> &Hv,   const Vector<Real>    &l,
                        const Vector<Real> &wold, const Vector<Real>    &wnew,
                        const Vector<Real> &uold, const Vector<Real>    &unew,
                        const Vector<Real> &z,    const TimeStamp<Real> &ts,
                        const bool sumInto = false) {
    const Real one(1);
    // Compute new/old Hessian of Lagrangian
    if (!sumInto) {
      obj_->hessVec_un_uo(Hv, wold, uold, unew, z, ts);
    }
    else {
      obj_->hessVec_un_uo(*udual_, wold, uold, unew, z, ts);
      Hv.plus(*udual_);
    }
    con_->applyAdjointHessian_uo_un(*udual_, l, wold, uold, unew, z, ts);
    Hv.axpy(-one,*udual_);
    // Compute new/new Hessian of Lagrangian
    obj_->hessVec_un_un(*udual_, wnew, uold, unew, z, ts);
    Hv.plus(*udual_);
    con_->applyAdjointHessian_un_un(*udual_, l, wnew, uold, unew, z, ts);
    Hv.axpy(-one,*udual_);
  }

  void computeOldMixedHessLag(Vector<Real> &Hv,
                        const Vector<Real> &v,    const Vector<Real>    &l,
                        const Vector<Real> &uold, const Vector<Real>    &unew,
                        const Vector<Real> &z,    const TimeStamp<Real> &ts,
                        const bool sumInto = false) {
    const Real one(1);
    // Compute new/old Hessian of Lagrangian
    if (!sumInto) {
      con_->applyAdjointHessian_z_un(Hv, l, v, uold, unew, z, ts);
    }
    else {
      con_->applyAdjointHessian_z_un(*udual_, l, v, uold, unew, z, ts);
      Hv.plus(*udual_);
    }
    obj_->hessVec_un_z(*udual_, v, uold, unew, z, ts);
    Hv.axpy(-one, *udual_);
  }

  void computeNewStateJacobian(Vector<Real> &Hv,   const Vector<Real>    &p,
                         const Vector<Real> &uold, const Vector<Real>    &unew,
                         const Vector<Real> &z,    const TimeStamp<Real> &ts,
                         const bool sumInto = false) {
    const Real one(1);
    if (!sumInto) {
      con_->applyAdjointJacobian_uo(Hv, p, uold, unew, z, ts);
      Hv.scale(-one);
    }
    else {
      con_->applyAdjointJacobian_uo(*udual_, p, uold, unew, z, ts);
      Hv.axpy(-one, *udual_);
    }
  }

  void computeNewStateHessLag(Vector<Real> &Hv,   const Vector<Real>    &l,
                        const Vector<Real> &wold, const Vector<Real>    &wnew,
                        const Vector<Real> &uold, const Vector<Real>    &unew,
                        const Vector<Real> &z,    const TimeStamp<Real> &ts,
                        const bool sumInto = false) {
    const Real one(1);
    // Compute old/new Hessian of Lagrangian
    if (!sumInto) {
      obj_->hessVec_uo_un(Hv, wnew, uold, unew, z, ts);
    }
    else {
      obj_->hessVec_uo_un(*udual_, wnew, uold, unew, z, ts);
      Hv.plus(*udual_);
    }
    con_->applyAdjointHessian_un_uo(*udual_, l, wnew, uold, unew, z, ts);
    Hv.axpy(-one,*udual_);
    // Compute old/old Hessian of Lagrangian
    obj_->hessVec_uo_uo(*udual_, wold, uold, unew, z, ts);
    Hv.plus(*udual_);
    con_->applyAdjointHessian_uo_uo(*udual_, l, wold, uold, unew, z, ts);
    Hv.axpy(-one,*udual_);
  }

  void computeNewMixedHessLag(Vector<Real> &Hv,
                        const Vector<Real> &v,    const Vector<Real>    &l,
                        const Vector<Real> &uold, const Vector<Real>    &unew,
                        const Vector<Real> &z,    const TimeStamp<Real> &ts,
                        const bool sumInto = false) {
    const Real one(1);
    // Compute new/old Hessian of Lagrangian
    if (!sumInto) {
      con_->applyAdjointHessian_z_uo(Hv, l, v, uold, unew, z, ts);
    }
    else {
      con_->applyAdjointHessian_z_uo(*udual_, l, v, uold, unew, z, ts);
      Hv.plus(*udual_);
    }
    obj_->hessVec_uo_z(*udual_, v, uold, unew, z, ts);
    Hv.axpy(-one,*udual_);
  }

  void advanceAdjointSens(Vector<Real> &p,          Vector<Real>    &rhs,
                    const Vector<Real> &uold, const Vector<Real>    &unew,
                    const Vector<Real> &z,    const TimeStamp<Real> &ts) {
    // Solve adjoint sensitivity on current time interval
    con_->applyInverseAdjointJacobian_un(p, rhs, uold, unew, z, ts);
  }

  /***************************************************************************/
  /************ Method for Hessian-times-a-vector computation ****************/
  /***************************************************************************/
  void computeControlHessLag(Vector<Real> &Hv,
                       const Vector<Real> &v,    const Vector<Real>    &l,
                       const Vector<Real> &uold, const Vector<Real>    &unew,
                       const Vector<Real> &z,    const TimeStamp<Real> &ts) {
    const Real one(1);
    // Compute Hessian of Lagrangian
    obj_->hessVec_z_z(Hv, v, uold, unew, z, ts);
    con_->applyAdjointHessian_z_z(*zdual_, l, v, uold, unew, z, ts);
    Hv.axpy(-one, *zdual_);
  }

  void addMixedHessLag(Vector<Real> &Hv,   const Vector<Real>    &l,
                 const Vector<Real> &wold, const Vector<Real>    &wnew,
                 const Vector<Real> &uold, const Vector<Real>    &unew,
                 const Vector<Real> &z,    const TimeStamp<Real> &ts) {
    const Real one(1);
    // Compute Hessian of Lagrangian on previous time interval
    obj_->hessVec_z_uo(*zdual_, wold, uold, unew, z, ts);
    Hv.axpy(-one, *zdual_);
    con_->applyAdjointHessian_uo_z(*zdual_, l, wold, uold, unew, z, ts);
    Hv.plus(*zdual_);
    // Compute Hessian of Lagrangian on current time interval
    obj_->hessVec_z_un(*zdual_, wnew, uold, unew, z, ts);
    Hv.axpy(-one, *zdual_);
    con_->applyAdjointHessian_un_z(*zdual_, l, wnew, uold, unew, z, ts);
    Hv.plus(*zdual_);
  }

  void addAdjointSens(Vector<Real> &Hv,   const Vector<Real>    &p,
                const Vector<Real> &uold, const Vector<Real>    &unew,
                const Vector<Real> &z,    const TimeStamp<Real> &ts) {
    con_->applyAdjointJacobian_z(*zdual_, p, uold, unew, z, ts);
    Hv.plus(*zdual_);
  }
};

} // namespace ROL

#endif // ROL_REDUCEDDYNAMICOBJECTIVE_HPP
