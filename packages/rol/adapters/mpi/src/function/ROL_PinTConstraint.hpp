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


#ifndef ROL_PINTCONSTRAINT_HPP
#define ROL_PINTCONSTRAINT_HPP

#include "ROL_TimeStamp.hpp"
#include "ROL_PinTVector.hpp"
#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_SerialConstraint.hpp"

#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StackedTimer.hpp"

namespace ROL {

/** This helper method builds a pint "state" vector for use in the 
    PinTConstraint class. 

    \param[in] communicators Structure with objects required for parallel-in-time communication.
                             (hint: this can be used to properly allocate the state and control vectors)
    \param[in] steps Number of steps 
    \param[in] localVector Spatial vector for a single time step.
  */
template <class Real>
Ptr<PinTVector<Real>>
buildStatePinTVector(const Ptr<const PinTCommunicators> & communicators,
                     int steps,
                     const Ptr<Vector<Real>> & localVector)
{ std::vector<int> stencil = {-1,0};
  return makePtr<PinTVector<Real>>(communicators,localVector,steps,stencil); }

/** This helper method builds a pint "control" vector for use in the 
    PinTConstraint class. 

    \param[in] communicators Structure with objects required for parallel-in-time communication.
                             (hint: this can be used to properly allocate the state and control vectors)
    \param[in] steps Number of steps 
    \param[in] localVector Spatial vector for a single time step.
  */
template <class Real>
Ptr<PinTVector<Real>>
buildControlPinTVector(const Ptr<const PinTCommunicators> & communicators,
                       int steps,
                       const Ptr<Vector<Real>> & localVector)
{ std::vector<int> stencil = {0};
  return makePtr<PinTVector<Real>>(communicators,localVector,steps,stencil); }

template<typename Real> 
class PinTConstraint : public ROL::Constraint_SimOpt<Real> {

  using V  = ROL::Vector<Real>;
  using PV = ROL::PartitionedVector<Real>;

  using size_type = typename std::vector<Real>::size_type;
  template<typename T> using Ptr = ROL::Ptr<T>;

private:

  // internal state members
  bool isInitialized_;

  //! Enumeration to define which components of a a vector are required.
  typedef enum {PAST,CURRENT,FUTURE,ALL} ETimeAccessor;

  Ptr<DynamicConstraint<Real>> stepConstraint_;
  Ptr<SerialConstraint<Real>>  timeDomainConstraint_;
  Ptr<const Vector<Real>>      initialCond_;

  Ptr<const std::vector<TimeStamp<Real>>>  userTimeStamps_;   // these are the untouched ones the user passes in
  Ptr<std::vector<TimeStamp<Real>>>        timeStamps_;       // these are used internally

  // preconditioner settings
  bool applyMultigrid_;         // default to block jacobi
  int maxLevels_;               // must turn on multigrid 
  int numSweeps_;
  Real omega_;
  Real globalScale_;
  std::vector<double> timePerLevel_; // given a MGRIT solve, what is the runtime for each level

  // Get the state vector required by the serial constraint object
  Ptr<Vector<Real>> getStateVector(const PinTVector<Real> & src)
  {
    // notice this ignores the stencil for now
   
    std::vector<Ptr<Vector<Real>>> vecs;

    vecs.push_back(src.getVectorPtr(-1));  // inject the initial condition

    for(int i=0;i<src.numOwnedSteps();i++) { // ignoring the stencil
      vecs.push_back(src.getVectorPtr(i)); 
    }
    return makePtr<PartitionedVector<Real>>(vecs);
  }

  // Get the control vector required by the serial constraint object
  Ptr<Vector<Real>> getControlVector(const PinTVector<Real> & src)
  {
    // notice this ignores the stencil for now

    std::vector<Ptr<Vector<Real>>> vecs;
 
    vecs.push_back(src.getVectorPtr(0)->clone());
      // FIXME-A: this control index is never used. This stems from not wanting to set
      //          an intial condition for the test function. But instead use the 
      //          "setSkipInitialCondition" option. This all needs to be rethought. 
      //          This is linked to the FIXME's in "getStateConstraint"
      //             -- ECC, June 11, 2018

    for(int i=0;i<src.numOwnedSteps();i++) { // ignore the intial step...ignoring the stencil
      vecs.push_back(src.getVectorPtr(i)); 
    }
    return makePtr<PartitionedVector<Real>>(vecs);
  }

public:

  //! Default constructor
  PinTConstraint()
    : Constraint_SimOpt<Real>()
    , isInitialized_(false)
    , applyMultigrid_(false)
    , maxLevels_(-1)
    , numSweeps_(1)
    , omega_(2.0/3.0)
    , globalScale_(0.99e0)
  { }

  /**
   * \brief Constructor
   *
   * Build a parallel-in-time constraint with a specified step constraint. This specifies
   * any communication you might need between processors and steps using "stencils".
   * Multigrid in time preconditioning is disabled by default.
   * 
   * \param[in] stepConstraint Constraint for a single step.
   * \param[in] initialCond Initial condition
   */
  PinTConstraint(const Ptr<DynamicConstraint<Real>> & stepConstraint,
                 const Ptr<const Vector<Real>> & initialCond,
                 const Ptr<std::vector<TimeStamp<Real>>> & timeStamps)
    : Constraint_SimOpt<Real>()
    , isInitialized_(false)
    , applyMultigrid_(false)
    , maxLevels_(-1)
    , numSweeps_(1)
    , omega_(2.0/3.0)
    , globalScale_(0.99e0)
  { 
    initialize(stepConstraint,initialCond,timeStamps);
  }

  /**
   * Set sweeps for jacobi and multigrid
   */
  void setSweeps(int s)
  { numSweeps_ = s; }

  /**
   * Set relaxation parater for jacobi and multigrid
   */
  void setRelaxation(Real o)
  { omega_ = o; }

  /**
   * Set the global scaling for the coupling constraints.
   */
  void setGlobalScale(Real o)
  { globalScale_ = o; }

  /**
   * Get the aggregate time per multigrid level. 
   */
  const std::vector<Real> & getTimePerLevel() const
  { return timePerLevel_; }

  /**
   * Clear out the aggregate MGRIT timers
   */
  void clearTimePerLevel()
  { return timePerLevel_.clear(); }

  /**
   * Turn on multigrid preconditioning in time with a specified number of levels.
   *
   * \param[in] maxLevels Largest number of levels
   */
  void applyMultigrid(int maxLevels)
  {
    applyMultigrid_ = true;
    maxLevels_ = maxLevels;
  }

  /**
   * \brief Turn off multigrid preconditioning in time (default off).
   */
  void disableMultigrid()
  {
    applyMultigrid_ = false;
    maxLevels_ = -1;
  }

  /**
   * \brief Check if multigrid is enabled.
   */
  bool multigridEnabled() const
  { return applyMultigrid_; }
  
  /**
   * \brief Get the time stamps for by level
   *
   * The current implementation is recursive and expensive, but its the easiest
   * way to do this right now. FIXME-B
   */
  Ptr<std::vector<TimeStamp<Real>>> getTimeStampsByLevel(int level) const
  {
    assert(level>=0); // precondition

    // base case
    if(level==0)
      return timeStamps_;

    Ptr<std::vector<TimeStamp<Real>>> higherLevel = getTimeStampsByLevel(level-1);

    // let's start easy!
    assert((higherLevel->size()-1) % 2 == 0);

    Ptr<std::vector<TimeStamp<Real>>> currentLevel 
        = ROL::makePtr<std::vector<ROL::TimeStamp<Real>>>((higherLevel->size()-1)/2+1);

    // build up the current level
    for(size_t k=0;k<currentLevel->size()-1;k++) {
      currentLevel->at(k+1).t.resize(2);
      currentLevel->at(k+1).t.at(0) = higherLevel->at(1+2*k+0).t.at(0);
      currentLevel->at(k+1).t.at(1) = higherLevel->at(1+2*k+1).t.at(1);
    }

    // setup false step (this is to ensure the control is properly handled
    {
      Real ta = currentLevel->at(1).t.at(0);
      Real tb = currentLevel->at(1).t.at(1);
      Real dt =  tb-ta;

      // the serial constraint should never see this!
      currentLevel->at(0).t.resize(2);
      currentLevel->at(0).t.at(0) = ta-dt;
      currentLevel->at(0).t.at(1) = ta;
    }

    return currentLevel;
  }

  /**
   * \brief Get the serial constraint associated with a level and 
   *        with an initial condition
   *
   * \param[in] ui    Initial condition 
   * \param[in] level Multigrid level (used only in preconditioning), level=0 is the finest
   */
  Ptr<SerialConstraint<Real>> 
  getSerialConstraint(const Ptr<ROL::Vector<Real>> & ui,int level)
  {
    auto timeStamps = getTimeStampsByLevel(level);

    if(ROL::is_nullPtr(timeDomainConstraint_)) {
      timeDomainConstraint_ = ROL::makePtr<SerialConstraint<Real>>(stepConstraint_,*ui,timeStamps);
        // FIXME-A: Change the Nt+1 back to Nt...
    }
    else {
      timeDomainConstraint_->setTimeStamp(timeStamps);
        // FIXME-A: Change the Nt+1 back to Nt...
      timeDomainConstraint_->setInitialCondition(*ui);
    }

    timeDomainConstraint_->setSkipInitialCondition(true);
    
    return timeDomainConstraint_;
  }

  /** \brief Initialize this class, setting up parallel distribution.
   
   */
  void initialize(const Ptr<DynamicConstraint<Real>> & stepConstraint,
                  const Ptr<const Vector<Real>> & initialCond,
                 const Ptr<std::vector<TimeStamp<Real>>> & timeStamps)
  {
    // initialize user member variables
    stepConstraint_ = stepConstraint;
    initialCond_ = initialCond;
    userTimeStamps_ = timeStamps;

    // build up the internally used time stamps
    ////////////////////////////////////////////////////////////////
   
    timeStamps_ = makePtr<std::vector<TimeStamp<Real>>>(timeStamps->size()+1);

    // setup false step (this is to ensure the control is properly handled
    {
      Real ta = timeStamps->at(0).t.at(0);
      Real tb = timeStamps->at(0).t.at(1);
      Real dt =  tb-ta;

      // the serial constraint should never see this!
      timeStamps_->at(0).t.resize(2);
      timeStamps_->at(0).t.at(0) = ta-dt;
      timeStamps_->at(0).t.at(1) = ta;
    }

    for(size_type k=0;k<timeStamps->size();k++) 
      timeStamps_->at(k+1).t = timeStamps->at(k).t;
  
    isInitialized_ = true;
  }

  void solve( V& c, V& u, const V& z, Real& tol ) override { 
    // solve is weird because it serializes in time. But we want to get it
    // working so that we can test, but it won't be a performant way to do this. 
    // We could do a parallel in time solve here but thats not really the focus.
    
    int level = 0;
    
    PinTVector<Real>       & pint_c = dynamic_cast<PinTVector<Real>&>(c);
    PinTVector<Real>       & pint_u = dynamic_cast<PinTVector<Real>&>(u);
    const PinTVector<Real> & pint_z = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here
 
    assert(pint_c.numOwnedSteps()==pint_u.numOwnedSteps());

    pint_z.boundaryExchange();
    pint_u.boundaryExchange(PinTVector<Real>::RECV_ONLY); // this is going to block

    // satisfy the time continuity constraint, note that this just using the identity matrix
    const std::vector<int> & stencil = pint_c.stencil();
    int timeRank = pint_u.communicators().getTimeRank();
    for(size_t i=0;i<stencil.size();i++) {
      int offset = stencil[i];
      if(offset<0 and timeRank>0) {

        // update the old time information
        pint_u.getVectorPtr(offset)->set(*pint_u.getRemoteBufferPtr(offset));        // this is the local value

        // this should be zero!
        pint_c.getVectorPtr(offset)->set(*pint_u.getVectorPtr(offset));             // this is the local value
        pint_c.getVectorPtr(offset)->axpy(-1.0,*pint_u.getRemoteBufferPtr(offset)); // this is the remote value
        pint_c.getVectorPtr(offset)->scale(globalScale_);
      }
      else if(offset<0 and timeRank==0) {  
        // this is the intial condition

        // update the old time information
        pint_u.getVectorPtr(offset)->set(*initialCond_);        

        // this should be zero!
        pint_c.getVectorPtr(offset)->set(*pint_u.getVectorPtr(offset));
        pint_c.getVectorPtr(offset)->axpy(-1.0,*initialCond_);
        pint_c.getVectorPtr(offset)->scale(globalScale_);
      }
    }

    auto part_c = getStateVector(pint_c);   
    auto part_u = getStateVector(pint_u);   
    auto part_z = getControlVector(pint_z); 

    // compute the constraint for this subdomain
    auto constraint = getSerialConstraint(pint_u.getVectorPtr(-1),level);

    constraint->solve(*part_c,*part_u,*part_z,tol);

    pint_u.boundaryExchange(PinTVector<Real>::SEND_ONLY);
  }  

  void value( V& c, const V& u, const V& z, Real& tol ) override {

    int level = 0;

    c.zero();   

    PinTVector<Real>       & pint_c = dynamic_cast<PinTVector<Real>&>(c);
    const PinTVector<Real> & pint_u = dynamic_cast<const PinTVector<Real>&>(u);
    const PinTVector<Real> & pint_z = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here
       
    assert(pint_c.numOwnedSteps()==pint_u.numOwnedSteps());

    // communicate neighbors, these are block calls
    pint_u.boundaryExchange();
    pint_z.boundaryExchange();

    int timeRank = pint_u.communicators().getTimeRank();

    // build in the time continuity constraint, note that this is just using the identity matrix
    const std::vector<int> & stencil = pint_c.stencil();
    for(size_t i=0;i<stencil.size();i++) {
      int offset = stencil[i];
      if(offset<0 and timeRank>0) {
        pint_c.getVectorPtr(offset)->set(*pint_u.getVectorPtr(offset));             // this is the local value
        pint_c.getVectorPtr(offset)->axpy(-1.0,*pint_u.getRemoteBufferPtr(offset)); // this is the remote value
        pint_c.getVectorPtr(offset)->scale(globalScale_);
      }
      else if(offset<0 and timeRank==0) {  
        // set the initial condition
        pint_c.getVectorPtr(offset)->set(*pint_u.getVectorPtr(offset));            
        pint_c.getVectorPtr(offset)->axpy(-1.0,*initialCond_); 
        pint_c.getVectorPtr(offset)->scale(globalScale_);
      }
    }

    auto part_c = getStateVector(pint_c);   // strip out initial condition/constraint
    auto part_u = getStateVector(pint_u);   // strip out initial condition/constraint
    auto part_z = getControlVector(pint_z); 

    // compute the constraint for this subdomain
    auto constraint = getSerialConstraint(pint_u.getVectorPtr(-1),level);

    constraint->value(*part_c,*part_u,*part_z,tol);
  } 

  void applyJacobian_1( V& jv, const V& v, const V& u, 
                       const V& z, Real& tol ) override {
    int level = 0;
    applyJacobian_1_leveled(jv,v,u,z,tol,level);
  }

  void applyJacobian_1_leveled( V& jv, const V& v, const V& u, 
                       const V& z, Real& tol,int level) {
    jv.zero();

    PinTVector<Real>       & pint_jv = dynamic_cast<PinTVector<Real>&>(jv);
    const PinTVector<Real> & pint_v  = dynamic_cast<const PinTVector<Real>&>(v);
    const PinTVector<Real> & pint_u  = dynamic_cast<const PinTVector<Real>&>(u);
    const PinTVector<Real> & pint_z  = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here
 
    assert(pint_jv.numOwnedSteps()==pint_u.numOwnedSteps());
    assert(pint_jv.numOwnedSteps()==pint_v.numOwnedSteps());

    // communicate neighbors, these are block calls
    pint_v.boundaryExchange();
    pint_u.boundaryExchange();
    pint_z.boundaryExchange();

    // differentiate the time continuity constraint, note that this just using the identity matrix
    int timeRank = pint_u.communicators().getTimeRank();
    const std::vector<int> & stencil = pint_jv.stencil();
    for(size_t i=0;i<stencil.size();i++) {
      int offset = stencil[i];
      if(offset<0) {
        pint_jv.getVectorPtr(offset)->set(*pint_v.getVectorPtr(offset));             // this is the local value

        // this is a hack to make sure that sensitivities with respect to the initial condition
        // don't get accidentally included
        if(timeRank>0) 
          pint_jv.getVectorPtr(offset)->axpy(-1.0,*pint_v.getRemoteBufferPtr(offset)); // this is the remote value

        pint_jv.getVectorPtr(offset)->scale(globalScale_);
      }
    }

    auto part_jv = getStateVector(pint_jv);   // strip out initial condition/constraint
    auto part_v  = getStateVector(pint_v);    // strip out initial condition/constraint
    auto part_u  = getStateVector(pint_u);    // strip out initial condition/constraint
    auto part_z  = getControlVector(pint_z); 

    // compute the constraint for this subdomain
    auto constraint = getSerialConstraint(pint_u.getVectorPtr(-1),level);

    constraint->applyJacobian_1(*part_jv,*part_v,*part_u,*part_z,tol);
  }

  void applyJacobian_1_leveled_approx( V& jv, const V& v, const V& u, 
                       const V& z, Real& tol,int level) {
    jv.zero();

    PinTVector<Real>       & pint_jv = dynamic_cast<PinTVector<Real>&>(jv);
    const PinTVector<Real> & pint_v  = dynamic_cast<const PinTVector<Real>&>(v);
    const PinTVector<Real> & pint_u  = dynamic_cast<const PinTVector<Real>&>(u);
    const PinTVector<Real> & pint_z  = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here
 
    assert(pint_jv.numOwnedSteps()==pint_u.numOwnedSteps());
    assert(pint_jv.numOwnedSteps()==pint_v.numOwnedSteps());

    // communicate neighbors, these are block calls
    pint_v.boundaryExchange();
    pint_u.boundaryExchange();
    pint_z.boundaryExchange();

    // differentiate the time continuity constraint, note that this just using the identity matrix
    int timeRank = pint_u.communicators().getTimeRank();
    const std::vector<int> & stencil = pint_jv.stencil();
    for(size_t i=0;i<stencil.size();i++) {
      int offset = stencil[i];
      if(offset<0) {
        pint_jv.getVectorPtr(offset)->set(*pint_v.getVectorPtr(offset));             // this is the local value

        // this is a hack to make sure that sensitivities with respect to the initial condition
        // don't get accidentally included
        // if(timeRank>0) 
        //   pint_jv.getVectorPtr(offset)->axpy(-1.0,*pint_v.getRemoteBufferPtr(offset)); // this is the remote value

        pint_jv.getVectorPtr(offset)->scale(globalScale_);
      }
    }

    auto part_jv = getStateVector(pint_jv);   // strip out initial condition/constraint
    auto part_v  = getStateVector(pint_v);    // strip out initial condition/constraint
    auto part_u  = getStateVector(pint_u);    // strip out initial condition/constraint
    auto part_z  = getControlVector(pint_z); 

    // compute the constraint for this subdomain
    auto constraint = getSerialConstraint(pint_u.getVectorPtr(-1),level);

    constraint->applyJacobian_1(*part_jv,*part_v,*part_u,*part_z,tol);
  }

  void applyJacobian_2( V& jv, const V& v, const V& u,
                        const V &z, Real &tol ) override { 
    int level = 0;
    applyJacobian_2_leveled(jv,v,u,z,tol,level);
  }

  void applyJacobian_2_leveled( V& jv, const V& v, const V& u,
                        const V &z, Real &tol,int level ) { 

    jv.zero();

    PinTVector<Real>       & pint_jv = dynamic_cast<PinTVector<Real>&>(jv);
    const PinTVector<Real> & pint_v  = dynamic_cast<const PinTVector<Real>&>(v);
    const PinTVector<Real> & pint_u  = dynamic_cast<const PinTVector<Real>&>(u);
    const PinTVector<Real> & pint_z  = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here
 
    assert(pint_jv.numOwnedSteps()==pint_u.numOwnedSteps());
    assert(pint_jv.numOwnedSteps()==pint_v.numOwnedSteps());

    // communicate neighbors, these are block calls
    pint_v.boundaryExchange();
    pint_u.boundaryExchange();
    pint_z.boundaryExchange();

    // differentiate the time continuity constraint, note that this just using the identity matrix
    const std::vector<int> & stencil = pint_jv.stencil();
    for(size_t i=0;i<stencil.size();i++) {
      int offset = stencil[i];
      if(offset<0) {
        pint_jv.getVectorPtr(offset)->zero();
      }
    }

    auto part_jv = getStateVector(pint_jv);   // strip out initial condition/constraint
    auto part_v  = getControlVector(pint_v);   
    auto part_u  = getStateVector(pint_u);    // strip out initial condition/constraint
    auto part_z  = getControlVector(pint_z); 

    // compute the constraint for this subdomain
    auto constraint = getSerialConstraint(pint_u.getVectorPtr(-1),level);

    constraint->applyJacobian_2(*part_jv,*part_v,*part_u,*part_z,tol);
   }

   void applyInverseJacobian_1( V& ijv, const V& v, const V& u,
                                const V& z, Real& tol) override {
     int level = 0;
     applyInverseJacobian_1_leveled(ijv,v,u,z,tol,level);
   }

   void applyInverseJacobian_1_leveled(V& ijv, 
                                       const V& v, 
                                       const V& u,
                                       const V& z, 
                                       Real& tol,
                                       int level) {
    // applyInverseJacobian_1 is weird because it serializes in time. But we want to get it
    // working so that we can test, but it won't be a performant way to do this. 
    // We could do a parallel in time solve here but thats not really the focus.
    
    PinTVector<Real>       & pint_ijv = dynamic_cast<PinTVector<Real>&>(ijv);
    const PinTVector<Real> & pint_v   = dynamic_cast<const PinTVector<Real>&>(v);
    const PinTVector<Real> & pint_u   = dynamic_cast<const PinTVector<Real>&>(u);
    const PinTVector<Real> & pint_z   = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here
 
    assert(pint_ijv.numOwnedSteps()==pint_u.numOwnedSteps());
    assert(pint_ijv.numOwnedSteps()==pint_v.numOwnedSteps());

    pint_z.boundaryExchange();
    pint_v.boundaryExchange();
    pint_u.boundaryExchange();
    pint_ijv.boundaryExchange(PinTVector<Real>::RECV_ONLY); // this is going to block

    // satisfy the time continuity constraint, note that this just using the identity matrix
    const std::vector<int> & stencil = pint_ijv.stencil();
    int timeRank = pint_u.communicators().getTimeRank();
    for(size_t i=0;i<stencil.size();i++) {
      int offset = stencil[i];
      if(offset<0 and timeRank>0) {

        // update the old time information
        pint_ijv.getVectorPtr(offset)->set(*pint_ijv.getRemoteBufferPtr(offset));    
        pint_ijv.getVectorPtr(offset)->axpy(1.0/globalScale_,*pint_v.getVectorPtr(offset));        
      }
      else if(offset<0 and timeRank==0) {  
        // this is the intial condition

        // update the old time information
        pint_ijv.getVectorPtr(offset)->set(*pint_v.getVectorPtr(offset));
        pint_ijv.getVectorPtr(offset)->scale(1.0/globalScale_);
      }
    }

    auto part_ijv = getStateVector(pint_ijv);   
    auto part_v   = getStateVector(pint_v);   
    auto part_u   = getStateVector(pint_u);   
    auto part_z   = getControlVector(pint_z); 

    // compute the constraint for this subdomain
    auto constraint = getSerialConstraint(pint_u.getVectorPtr(-1),level);

    constraint->applyInverseJacobian_1(*part_ijv,*part_v,*part_u,*part_z,tol);

    pint_ijv.boundaryExchange(PinTVector<Real>::SEND_ONLY);
   }


   void applyAdjointJacobian_1( V& ajv, 
                                const V& v, 
                                const V& u,
                                const V& z, 
                                Real& tol) override {
     int level = 0;
     applyAdjointJacobian_1_leveled(ajv,v,u,z,tol,level);
   }

   /**
    * This is a convenience function for multi-grid that gives access
    * to a "leveled" version of the adjoint jacobian.
    */
   void applyAdjointJacobian_1_leveled( V& ajv, 
                                        const V& v, 
                                        const V& u,
                                        const V& z, 
                                        Real& tol,
                                        int level) {
    PinTVector<Real>       & pint_ajv = dynamic_cast<PinTVector<Real>&>(ajv);
    const PinTVector<Real> & pint_v   = dynamic_cast<const PinTVector<Real>&>(v);
    const PinTVector<Real> & pint_u   = dynamic_cast<const PinTVector<Real>&>(u);
    const PinTVector<Real> & pint_z   = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here
      
    assert(pint_v.numOwnedSteps()==pint_u.numOwnedSteps());
    assert(pint_ajv.numOwnedSteps()==pint_u.numOwnedSteps());

    // we need to make sure this has all zeros to begin with (this includes boundary exchange components)
    pint_ajv.zeroAll();

    // communicate neighbors, these are block calls
    pint_v.boundaryExchange();
    pint_u.boundaryExchange();
    pint_z.boundaryExchange();

    auto part_ajv = getStateVector(pint_ajv);   // strip out initial condition/constraint
    auto part_v   = getStateVector(pint_v);   
    auto part_u   = getStateVector(pint_u);    // strip out initial condition/constraint
    auto part_z   = getControlVector(pint_z); 

    // compute the constraint for this subdomain
    auto constraint = getSerialConstraint(pint_u.getVectorPtr(-1),level);

    constraint->applyAdjointJacobian_1(*part_ajv,*part_v,*part_u,*part_z,*part_v,tol);

    // handle the constraint adjoint
    const std::vector<int> & stencil = pint_u.stencil();
    for(size_t i=0;i<stencil.size();i++) {
      int offset = stencil[i];
      if(offset<0) {
        pint_ajv.getVectorPtr(offset)->axpy(globalScale_,*pint_v.getVectorPtr(offset));
        pint_ajv.getRemoteBufferPtr(offset)->set(*pint_v.getVectorPtr(offset)); // this will be sent to the remote processor
        pint_ajv.getRemoteBufferPtr(offset)->scale(-globalScale_);
      }
    }

    // this sums from the remote buffer into the local buffer which is part of the vector
    pint_ajv.boundaryExchangeSumInto();
   }

   /**
    * This is a convenience function for multi-grid that gives access
    * to a "leveled" version of the adjoint jacobian.
    */
   void applyAdjointJacobian_1_leveled_approx( V& ajv, 
                                        const V& v, 
                                        const V& u,
                                        const V& z, 
                                        Real& tol,
                                        int level) {
    PinTVector<Real>       & pint_ajv = dynamic_cast<PinTVector<Real>&>(ajv);
    const PinTVector<Real> & pint_v   = dynamic_cast<const PinTVector<Real>&>(v);
    const PinTVector<Real> & pint_u   = dynamic_cast<const PinTVector<Real>&>(u);
    const PinTVector<Real> & pint_z   = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here
      
    assert(pint_v.numOwnedSteps()==pint_u.numOwnedSteps());
    assert(pint_ajv.numOwnedSteps()==pint_u.numOwnedSteps());

    // we need to make sure this has all zeros to begin with (this includes boundary exchange components)
    pint_ajv.zeroAll();

    // communicate neighbors, these are block calls
    pint_v.boundaryExchange();
    pint_u.boundaryExchange();
    pint_z.boundaryExchange();

    auto part_ajv = getStateVector(pint_ajv);   // strip out initial condition/constraint
    auto part_v   = getStateVector(pint_v);   
    auto part_u   = getStateVector(pint_u);    // strip out initial condition/constraint
    auto part_z   = getControlVector(pint_z); 

    // compute the constraint for this subdomain
    auto constraint = getSerialConstraint(pint_u.getVectorPtr(-1),level);

    constraint->applyAdjointJacobian_1(*part_ajv,*part_v,*part_u,*part_z,*part_v,tol);

    // handle the constraint adjoint
    const std::vector<int> & stencil = pint_u.stencil();
    for(size_t i=0;i<stencil.size();i++) {
      int offset = stencil[i];
      if(offset<0) {
        pint_ajv.getVectorPtr(offset)->axpy(globalScale_,*pint_v.getVectorPtr(offset));
        pint_ajv.getRemoteBufferPtr(offset)->set(*pint_v.getVectorPtr(offset)); // this will be sent to the remote processor
        pint_ajv.getRemoteBufferPtr(offset)->scale(-globalScale_);
        pint_ajv.getRemoteBufferPtr(offset)->zero();
      }
    }

    // this sums from the remote buffer into the local buffer which is part of the vector
    pint_ajv.boundaryExchangeSumInto();
   }

   void applyAdjointJacobian_2( V& ajv,  const V& v, const V& u,
                                const V& z, Real& tol ) override {
     int level = 0;
     applyAdjointJacobian_2_leveled(ajv,v,u,z,tol,level);
   }

   void applyAdjointJacobian_2_leveled( V& ajv,  
                                        const V& v, 
                                        const V& u,
                                        const V& z, 
                                        Real& tol,
                                        int level ) {

    PinTVector<Real>       & pint_ajv = dynamic_cast<PinTVector<Real>&>(ajv);
    const PinTVector<Real> & pint_v   = dynamic_cast<const PinTVector<Real>&>(v);
    const PinTVector<Real> & pint_u   = dynamic_cast<const PinTVector<Real>&>(u);
    const PinTVector<Real> & pint_z   = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here
      
    assert(pint_v.numOwnedSteps()==pint_u.numOwnedSteps());
    assert(pint_ajv.numOwnedSteps()==pint_u.numOwnedSteps());

    // we need to make sure this has all zeros to begin with (this includes boundary exchange components)
    pint_ajv.zeroAll();

    // communicate neighbors, these are block calls
    pint_v.boundaryExchange();
    pint_u.boundaryExchange();
    pint_z.boundaryExchange();

    auto part_ajv = getControlVector(pint_ajv);   // strip out initial condition/constraint
    auto part_v   = getStateVector(pint_v);   
    auto part_u   = getStateVector(pint_u);    // strip out initial condition/constraint
    auto part_z   = getControlVector(pint_z); 

    // compute the constraint for this subdomain
    auto constraint = getSerialConstraint(pint_u.getVectorPtr(-1),level);

    constraint->applyAdjointJacobian_2(*part_ajv,*part_v,*part_u,*part_z,tol);

    // no constraint for controls hanlding because that doesn't work yet!
   }

   void applyInverseAdjointJacobian_1( V& iajv, const V& v, const V& u,
                                       const V& z, Real& tol) override {
     int level = 0;
     applyInverseAdjointJacobian_1_leveled(iajv,v,u,z,tol,level);
   }

   void applyInverseAdjointJacobian_1_leveled( V& iajv, 
                                               const V& v, 
                                               const V& u,
                                               const V& z, 
                                               Real& tol,
                                               int level) {

    // applyInverseAdjointJacobian_1 is weird because it serializes in time. But we want to get it
    // working so that we can test, but it won't be a performant way to do this. 
    // We could do a parallel in time solve here but thats not really the focus.
    
    // this is an inefficient hack (see line below where pint_v is modified!!!!)
    ///////////////////////////////////
    auto v_copy = v.clone();
    v_copy->set(v);
    ///////////////////////////////////
    
    PinTVector<Real>       & pint_iajv = dynamic_cast<PinTVector<Real>&>(iajv);
    const PinTVector<Real> & pint_v    = dynamic_cast<const PinTVector<Real>&>(*v_copy);
    const PinTVector<Real> & pint_u    = dynamic_cast<const PinTVector<Real>&>(u);
    const PinTVector<Real> & pint_z    = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here
       //
    int timeRank = pint_u.communicators().getTimeRank();
    int timeSize = pint_u.communicators().getTimeSize();
 
    assert(pint_iajv.numOwnedSteps()==pint_u.numOwnedSteps());
    assert(pint_iajv.numOwnedSteps()==pint_v.numOwnedSteps());

    pint_iajv.zero();

    pint_z.boundaryExchange();
    pint_v.boundaryExchange();
    pint_u.boundaryExchange();
    pint_iajv.boundaryExchangeSumInto(PinTVector<Real>::RECV_ONLY); // this is going to block

    auto part_iajv = getStateVector(pint_iajv);   
    auto part_v    = getStateVector(pint_v);   
    auto part_u    = getStateVector(pint_u);   
    auto part_z    = getControlVector(pint_z); 

    const std::vector<int> & stencil = pint_iajv.stencil();

    //////////////////////////////////////////////////////////////////////////////////////

    // we need to modify the RHS vector if this is an interior boundary
    if(timeRank+1 < timeSize) {
      int owned_steps = pint_v.numOwnedSteps();
      for(size_t i=0;i<stencil.size();i++) {
        if(stencil[i]<0) {
          // this is why the hack above is necessary
          pint_v.getVectorPtr(owned_steps+stencil[i])->axpy(globalScale_,*pint_iajv.getVectorPtr(owned_steps+stencil[i]));
        }
      }
    }

    //////////////////////////////////////////////////////////////////////////////////////
    
    // compute the constraint for this subdomain
    auto constraint = getSerialConstraint(pint_u.getVectorPtr(-1),level);

    constraint->applyInverseAdjointJacobian_1(*part_iajv,*part_v,*part_u,*part_z,tol);

    //////////////////////////////////////////////////////////////////////////////////////
    
    // satisfy the time continuity constraint, note that this just using the identity matrix
    for(size_t i=0;i<stencil.size();i++) {
      int offset = stencil[i];
      if(offset<0) {
        // you would think that was the weird term here. But there isn't as its automatically
        // included by the applyInverseAdjointJacobian_1 call with the "skipInitialCondition"
        // flag on.

        // sending this to the left
        pint_iajv.getVectorPtr(offset)->scale(1.0/globalScale_);
        pint_iajv.getRemoteBufferPtr(offset)->set(*pint_iajv.getVectorPtr(offset));        
      }
    }

    pint_iajv.boundaryExchangeSumInto(PinTVector<Real>::SEND_ONLY);

   }

   // Done in parallel, no blocking, solve a linear system on this processor
   void invertTimeStepJacobian(Vector<Real>       & pv,
                               const Vector<Real> & v,
                               const Vector<Real> & u,
                               const Vector<Real> & z,
                               Real &tol,
                               int level) 
   {
     PinTVector<Real> & pint_pv = dynamic_cast<PinTVector<Real>&>(pv);
     const PinTVector<Real> & pint_v = dynamic_cast<const PinTVector<Real>&>(v);
     const PinTVector<Real> & pint_u = dynamic_cast<const PinTVector<Real>&>(u);
     const PinTVector<Real> & pint_z = dynamic_cast<const PinTVector<Real>&>(z);

     invertTimeStepJacobian(pint_pv,pint_v,pint_u,pint_z,tol,level);
   }
   void invertAdjointTimeStepJacobian(Vector<Real>       & pv,
                               const Vector<Real> & v,
                               const Vector<Real> & u,
                               const Vector<Real> & z,
                               Real &tol,
                               int level) 
   {
     PinTVector<Real> & pint_pv = dynamic_cast<PinTVector<Real>&>(pv);
     const PinTVector<Real> & pint_v = dynamic_cast<const PinTVector<Real>&>(v);
     const PinTVector<Real> & pint_u = dynamic_cast<const PinTVector<Real>&>(u);
     const PinTVector<Real> & pint_z = dynamic_cast<const PinTVector<Real>&>(z);

     invertAdjointTimeStepJacobian(pint_pv,pint_v,pint_u,pint_z,tol,level);
   }

   // Done in parallel, no blocking, solve a linear system on this processor
   void invertTimeStepJacobian(PinTVector<Real>       & pint_pv,
                               const PinTVector<Real> & pint_v,
                               const PinTVector<Real> & pint_u,
                               const PinTVector<Real> & pint_z,
                               Real &tol,
                               int level) 
   {
     // apply the time continuity constraint, note that this just using the identity matrix
     const std::vector<int> & stencil = pint_pv.stencil();
     for(size_t i=0;i<stencil.size();i++) {
       int offset = stencil[i];
       if(offset<0) {
         pint_pv.getVectorPtr(offset)->set(*pint_v.getVectorPtr(offset));             // this is the local value
         pint_pv.getVectorPtr(offset)->scale(1.0/globalScale_);
       }
     }

     auto part_pv  = getStateVector(pint_pv);   // strip out initial condition/constraint
     auto part_v   = getStateVector(pint_v);   
     auto part_u   = getStateVector(pint_u);    // strip out initial condition/constraint
     auto part_z   = getControlVector(pint_z); 

     // compute the constraint for this subdomain
     auto constraint = getSerialConstraint(pint_u.getVectorPtr(-1),level);

     constraint->applyInverseJacobian_1(*part_pv,*part_v,*part_u,*part_z,tol);
   }
 
   // Done in parallel, no blocking, solve a linear system on this processor
   void invertAdjointTimeStepJacobian(PinTVector<Real>       & pint_pv,
                                      const PinTVector<Real> & pint_v,
                                      const PinTVector<Real> & pint_u,
                                      const PinTVector<Real> & pint_z,
                                      Real &tol,
                                      int level)
   {
     auto part_pv  = getStateVector(pint_pv);   // strip out initial condition/constraint
     auto part_v   = getStateVector(pint_v);   
     auto part_u   = getStateVector(pint_u);    // strip out initial condition/constraint
     auto part_z   = getControlVector(pint_z); 

     // compute the constraint for this subdomain
     auto constraint = getSerialConstraint(pint_u.getVectorPtr(-1),level);

     constraint->applyInverseAdjointJacobian_1(*part_pv,*part_v,*part_u,*part_z,tol);
    
     // apply the time continuity constraint, note that this just using the identity matrix
     const std::vector<int> & stencil = pint_pv.stencil();
     for(size_t i=0;i<stencil.size();i++) {
       int offset = stencil[i];
       if(offset<0) {
         // pint_pv.getVectorPtr(offset)->set(*temp);
           // NOT good enough for multi step!!!!
         pint_pv.getVectorPtr(offset)->scale(1.0/globalScale_);
       }
     }
   }
 
   void blockJacobiSweep(PinTVector<Real>       & pint_pv,
                         const PinTVector<Real> & pint_rhs,
                         const PinTVector<Real> & pint_u,
                         const PinTVector<Real> & pint_z,
                         Real &tol,
                         int level)
   {
     pint_rhs.boundaryExchange();
 
     Ptr<Vector<Real>> scratch = pint_rhs.clone();
     PinTVector<Real> pint_scratch  = dynamic_cast<const PinTVector<Real>&>(*scratch);

     // do block Jacobi smoothing
     invertTimeStepJacobian(pint_scratch, pint_rhs, pint_u,pint_z,tol,level);
     // applyInverseJacobian_1_leveled(pint_scratch, pint_rhs, pint_u,pint_z,tol,level);

     invertAdjointTimeStepJacobian(pint_pv, pint_scratch, pint_u,pint_z,tol,level);
     // applyInverseAdjointJacobian_1_leveled(pint_pv, pint_scratch, pint_u,pint_z,tol,level);

     pint_pv.scale(-1.0);

     int timeRank = pint_u.communicators().getTimeRank();
     if(false)
     {
       std::stringstream ss;
       ss << timeRank << " AFTER SOLVE = " << std::endl;
       for(int i=-1;i<pint_u.numOwnedSteps();i++) 
         ss << "   " << timeRank << " - " << i << ": "<< pint_pv.getVectorPtr(i)->norm()  << std::endl;
       std::cout << ss.str() << std::endl;; 
     }
   }
 
   void schurComplementAction(Vector<Real>       & sc_v,
                              const Vector<Real> & v,
                              const Vector<Real> & u,
                              const Vector<Real> & z,
                              Real &tol,
                              int level)
   {
     auto scratch_u = u.clone(); 
     auto scratch_z = z.clone(); 

     scratch_u->zero();
     scratch_z->zero();
     sc_v.zero();
 
     // do boundary exchange
     {
       const PinTVector<Real>       & pint_v = dynamic_cast<const PinTVector<Real>&>(v);
       pint_v.boundaryExchange();
     }
 
     // compute J_1 * J_1^T
     applyAdjointJacobian_1_leveled(*scratch_u,v,u,z,tol,level)  ;
     applyJacobian_1_leveled(sc_v,*scratch_u,u,z,tol,level);
 
     // compute J_2 * J_2^T
     applyAdjointJacobian_2_leveled(*scratch_z,v,u,z,tol,level);
     applyJacobian_2_leveled(*scratch_u,*scratch_z,u,z,tol,level);
     
     // compute J_1 * J_1^T + J_2 * J_2^T
     // sc_v.axpy(1.0,*scratch_u);
     sc_v.scale(-1.0);
   }
 
   virtual void weightedJacobiRelax(PinTVector<Real> &pint_pv,
                                    const PinTVector<Real> & pint_v,
                                    const PinTVector<Real> & pint_u,
                                    const PinTVector<Real> & pint_z,
                                    Real omega,
                                    int numSweeps,
                                    Real &tol,
                                    int level) 
   {
     auto residual = pint_v.clone();
     auto dx = pint_pv.clone();

     int timeRank = pint_u.communicators().getTimeRank();
     if(false)
     {
       std::stringstream ss;
       ss << timeRank << " INITIAL = " << std::endl;
       for(int i=-1;i<pint_u.numOwnedSteps();i++) 
         ss << "   " << timeRank << " - " << i << ": "<< pint_pv.getVectorPtr(i)->norm()  << std::endl;
       std::cout << ss.str() << std::endl;; 
     }
     if(false)
     {
       std::stringstream ss;
       ss << timeRank << " RHS = " << std::endl;
       for(int i=-1;i<pint_u.numOwnedSteps();i++) 
         ss << "   " << timeRank << " - " << i << ": "<< pint_v.getVectorPtr(i)->norm()  << std::endl;
       std::cout << ss.str() << std::endl;; 
     }
     

     // force intial guess to zero
     dx->scale(0.0);
     residual->set(pint_v);           
 
     PinTVector<Real> & pint_residual = dynamic_cast<PinTVector<Real>&>(*residual);
     PinTVector<Real> & pint_dx       = dynamic_cast<PinTVector<Real>&>(*dx);

     for(int s=0;s<numSweeps;s++) {

       // update the residual
       schurComplementAction(pint_residual,pint_pv,pint_u,pint_z,tol,level); // A*x
 
       pint_residual.scale(-1.0);                            // -A*x
       pint_residual.plus(pint_v);                           // b - A*x

       // compute and apply the correction
       blockJacobiSweep(pint_dx,pint_residual,pint_u,pint_z,tol,level);
 
       pint_pv.axpy(omega,pint_dx);
     }
   }

   virtual void multigridInTime(Vector<Real> & pv,
                                const Vector<Real> & v,
                                const Vector<Real> & u,
                                const Vector<Real> & z,
                                Real omega,
                                int numSweeps,
                                Real &tol,
                                int level) 
   {
     PinTVector<Real> & pint_pv = dynamic_cast<PinTVector<Real>&>(pv);
     const PinTVector<Real> & pint_v = dynamic_cast<const PinTVector<Real>&>(v);
     const PinTVector<Real> & pint_u = dynamic_cast<const PinTVector<Real>&>(u);
     const PinTVector<Real> & pint_z = dynamic_cast<const PinTVector<Real>&>(z);

     multigridInTime(pint_pv,pint_v,pint_u,pint_z,omega,numSweeps,tol,level);
   }

   virtual void multigridInTime(PinTVector<Real> &pint_x,
                                const PinTVector<Real> & pint_b,
                                const PinTVector<Real> & pint_u,
                                const PinTVector<Real> & pint_z,
                                Real omega,
                                int numSweeps,
                                Real &tol,
                                int level) 
   {
     // sanity check the input (preconditions)
     assert(multigridEnabled());
     assert(maxLevels_>=0);

     auto residual = pint_b.clone();
     auto correction = pint_x.clone();
     PinTVector<Real> & pint_residual = dynamic_cast<PinTVector<Real>&>(*residual);

     pint_b.boundaryExchange();
     pint_u.boundaryExchange();
     pint_z.boundaryExchange();

     if(level==maxLevels_) {
       // compute fine residual
       schurComplementAction(pint_residual,pint_x,pint_u,pint_z,tol,level); // A*x
 
       pint_residual.scale(-1.0);                       // -A*x 
       pint_residual.plus(pint_b);                      // b - A*x

       // base case !!!
       auto scratch = pint_x.clone();
       scratch->zero();
     
       // compute -J_1^-T * J_1^-1
       applyInverseJacobian_1_leveled(*scratch,pint_residual,pint_u,pint_z,tol,level);
       applyInverseAdjointJacobian_1_leveled(*correction,*scratch,pint_u,pint_z,tol,level);

       pint_x.axpy(-1.0,*correction);

       return;
     }

     // pre-smooth
     /////////////////////////////////////////////////////////////////////////////////////
     weightedJacobiRelax(pint_x,pint_b,pint_u,pint_z,omega,numSweeps,tol,level);

     // compute fine residual
     schurComplementAction(pint_residual,pint_x,pint_u,pint_z,tol,level); // A*x
 
     pint_residual.scale(-1.0);                       // -A*x 
     pint_residual.plus(pint_b);                      // b - A*x

     int timeRank = pint_residual.communicators().getTimeRank();

     // coarse solve
     /////////////////////////////////////////////////////////////////////////////////////
     if(true)
     {
       auto crs_u          = allocateSimVector(pint_x,level+1);
       auto crs_z          = allocateOptVector(pint_z,level+1);
       auto crs_residual   = allocateSimVector(pint_x,level+1);
       auto crs_correction = allocateSimVector(pint_x,level+1);

       crs_correction->zero();

       restrictSimVector(pint_u,*crs_u);               // restrict the control to the coarse level
       restrictOptVector(pint_z,*crs_z);               // restrict the state to the coarse level
       restrictSimVector(pint_residual,*crs_residual);  

       multigridInTime(*crs_correction,*crs_residual,*crs_u,*crs_z,omega,numSweeps,tol,level+1);

       prolongSimVector(*crs_correction,*correction);  // prolongate the correction 

       pint_x.plus(*correction);                          // x+c
     }

     // post-smooth
     /////////////////////////////////////////////////////////////////////////////////////
     
     weightedJacobiRelax(pint_x,pint_b,pint_u,pint_z,omega,numSweeps,tol,level);
   }
 
 
   virtual void applyPreconditioner(Vector<Real> &pv,
                                    const Vector<Real> &v,
                                    const Vector<Real> &x,
                                    const Vector<Real> &g,
                                    Real &tol) 
   {
     const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(x);
     auto u = xs.get_1(); // sim
     auto z = xs.get_2(); // opt
 
     // update the solution vectors
     this->update(x);
 
     // apply inverse of state jacobian on a single time step
     
     PinTVector<Real>       & pint_pv = dynamic_cast<PinTVector<Real>&>(pv);
     const PinTVector<Real> & pint_v  = dynamic_cast<const PinTVector<Real>&>(v);
     const PinTVector<Real> & pint_u  = dynamic_cast<const PinTVector<Real>&>(*u);
     const PinTVector<Real> & pint_z  = dynamic_cast<const PinTVector<Real>&>(*z);
       // its possible we won't always want to cast to a PinT vector here
  
     assert(pint_pv.numOwnedSteps()==pint_u.numOwnedSteps());
     assert(pint_pv.numOwnedSteps()==pint_v.numOwnedSteps());
 
     // these don't change, do boundary exchange
     pint_u.boundaryExchange();
     pint_z.boundaryExchange();
 
     // now we are going to do a relaxation sweep: x = xhat+omega * Minv * (b-A*xhat)
     /////////////////////////////////////////////////////////////////////////////////////////
     
     int numSweeps = numSweeps_;
     Real omega = omega_;
     
     int level = 0;
     if(not multigridEnabled()) {
       pint_pv.zero();
       weightedJacobiRelax(pint_pv,pint_v,pint_u,pint_z,omega,numSweeps,tol,level);
     }
     else if(multigridEnabled()) {
       pint_pv.zero();
       multigridInTime(pint_pv,pint_v,pint_u,pint_z,omega,numSweeps,tol,level);
     }
     else {
       // unreachable :(
       pv.set(v.dual());
     } 
   }

   // restriction and prolongation functions
   ///////////////////////////////////////////////////////////////////////////////////////////
   
   /**
    * \brief Allocate a simulation space vector at a particular multigrid level.
    */
   Ptr<Vector<Real>> allocateSimVector(const Vector<Real> & level_0_ref,int level) const
   {
     const PinTVector<Real> & pint_ref  = dynamic_cast<const PinTVector<Real>&>(level_0_ref);
     auto comm = pint_ref.communicatorsPtr();
     
     Ptr<std::vector<TimeStamp<Real>>> stamps  = getTimeStampsByLevel(level);
     return buildStatePinTVector(comm,comm->getTimeSize()*(stamps->size()-1),pint_ref.getVectorPtr(0)->clone());
   }

   /**
    * \brief Allocate a simulation space vector at a particular multigrid level.
    */
   Ptr<Vector<Real>> allocateOptVector(const Vector<Real> & level_0_ref,int level) const
   {
     const PinTVector<Real> & pint_ref  = dynamic_cast<const PinTVector<Real>&>(level_0_ref);
     auto comm = pint_ref.communicatorsPtr();
     
     Ptr<std::vector<TimeStamp<Real>>> stamps  = getTimeStampsByLevel(level);
     return buildControlPinTVector(comm,comm->getTimeSize()*(stamps->size()-1),pint_ref.getVectorPtr(0)->clone());
   }
   
   /**
    * \brief Restrict a simulation space vector
    *
    * Currently doing the following. Let the coarse distribution be defined locally
    * to contain Np steps. The fine discretization is denoted by a X, while the coarse
    * discretization is denoted by a x. A superscript indicates the relative processor index.
    * If no index is included then it is assumed to be 0 (denoting no relative change in
    * processor index). The restriction for the sim vector thus computes
    *
    *   x_i = X_i                                  i = -1           (this is the "virtual variable")
    *   x_i = X_{2*i} + X_{2*i+1} + X_{2*i+2}      0 <= i < Np - 1
    *   x_i = X_{2*i} + X_{2*i+1} + X_0^1          i = Np-1         (this is the end boundary on this processor)
    *
    * This ensures that for all interior points the coarse grid is the sum of the three neighbors.
    * So now we need to divide the coarse x vector by 1/3 except at the end point at the final
    * time step. This needs to be scaled by 1/2. Note that in the code we scale each component as needed
    * rather than scaling the global vector and the fixing up only the final time step. There is a (very) modest
    * efficiency gain here.
    *
    * \param[in] inputLevel Multigrid level of the input vector (not currently used)
    * \param[in] input The input vector to be "restricted" 
    * \param[in] output The output vector resulting from restriction, must be at level inputLevel+1
    */
   void restrictSimVector(const Vector<Real> & input,Vector<Real> & output)
   {
     const PinTVector<Real> & pint_input  = dynamic_cast<const PinTVector<Real>&>(input);
     PinTVector<Real>       & pint_output = dynamic_cast<PinTVector<Real>&>(output);

     int Np = pint_output.numOwnedSteps();

     // this is the virtual variable (set it to -1)
     //     x_i = X_i                                  i = -1           (this is the "virtual variable")
     ////////////////////////////////////////////////////////////////////////////////////////////////////////
     {
       ROL::Vector<Real> & out = *pint_output.getVectorPtr(-1);
       out.set(*pint_input.getVectorPtr(-1));
     }

     // handle interior points
     //     x_i = X_{2*i} + X_{2*i+1} + X_{2*i+2}      0 <= i < Np - 1
     ////////////////////////////////////////////////////////////////////////////////////////////////////////
     for(int k=0;k<Np-1;k++) {
       ROL::Vector<Real> & out = *pint_output.getVectorPtr(k);

       out.zero();
       out.axpy(1.0,*pint_input.getVectorPtr(2*k+0)); 
       out.axpy(1.0,*pint_input.getVectorPtr(2*k+1)); 
       out.axpy(1.0,*pint_input.getVectorPtr(2*k+2)); 
       
       out.scale(1.0/3.0);
     }

     // handle the end point on this processor
     //   x_i = X_{2*i} + X_{2*i+1} + X_0^1          i = Np-1         (this is the end boundary on this processor)
     ////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     {
       int k = Np-1;
       ROL::Vector<Real> & out = *pint_output.getVectorPtr(k);
       out.zero();
       out.axpy(1.0,*pint_input.getVectorPtr(2*k+0)); 
       out.axpy(1.0,*pint_input.getVectorPtr(2*k+1)); 
     
       pint_output.getRemoteBufferPtr(-1)->set(*pint_input.getVectorPtr(0));
       pint_output.boundaryExchangeSumInto();   // this adds the X_0^1 into the previous rank

       // scale end boundary point by 1/2, or the other end points by 1/3
       int timeRank = pint_input.communicators().getTimeRank();
       int timeSize = pint_input.communicators().getTimeSize();

       if(timeRank+1==timeSize)
         out.scale(1.0/2.0);
       else
         out.scale(1.0/3.0);
     }
   }
   
   /**
    * \brief Restrict a control space vector
    *
    * Currently assumes a piecewise constant control
    *
    * \param[in] inputLevel Multigrid level of the input vector (not currently used)
    * \param[in] input The input vector to be "restricted" 
    * \param[in] output The output vector resulting from restriction, must be at level inputLevel+1
    */
   void restrictOptVector(const Vector<Real> & input,Vector<Real> & output)
   {
     const PinTVector<Real> & pint_input  = dynamic_cast<const PinTVector<Real>&>(input);
     PinTVector<Real>       & pint_output = dynamic_cast<PinTVector<Real>&>(output);

     // handle interior
     for(int k=0;k<pint_output.numOwnedSteps();k++) {
       ROL::Vector<Real> & out = *pint_output.getVectorPtr(k);

       out.axpy(1.0/2.0,*pint_input.getVectorPtr(2*k+0)); 
       out.axpy(1.0/2.0,*pint_input.getVectorPtr(2*k+1)); 
     }
   }

   /**
    * \brief Prolong a simulation space vector
    *
    * Currently doing the following. Let the coarse distribution be defined locally
    * to contain Np steps. The fine discretization is denoted by a X, while the coarse
    * discretization is denoted by a x. A superscript indicates the relative processor index.
    * If no index is included then it is assumed to be 0 (denoting no relative change in
    * processor index). The restriction for the sim vector thus computes
    *
    *   X_0 = (x_0 + x_{Np-1}^{-1})/2                               (initial condition on time domain)
    *   X_{2*i+1} = x_i                           -1 <= i < Np      (injection variables, including virtual variable)
    *   X_{2*i+2} = (x_i + x_{i+1})/2              0 <= i < Np - 1  (averaged variables, created at fine level)
    *
    * Note: Currently if the timeRank==0, the initial condition is assumed to be zero. 
    *
    * \param[in] inputLevel Multigrid level of the input vector (not currently used)
    * \param[in] input The input vector to be "restricted" 
    * \param[in] output The output vector resulting from restriction, must be at level inputLevel+1
    */
   void prolongSimVector(const Vector<Real> & input,Vector<Real> & output)
   {
     const PinTVector<Real> & pint_input  = dynamic_cast<const PinTVector<Real>&>(input); // coarse vector
     PinTVector<Real>       & pint_output = dynamic_cast<PinTVector<Real>&>(output);

     int Np = pint_input.numOwnedSteps();

     //   X_0 = (x_0 + x_{Np-1}^{-1})/2                               (initial condition on time domain)
     ////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     {
       int timeRank = pint_input.communicators().getTimeRank();

       ROL::Vector<Real> & out = *pint_output.getVectorPtr(0);

       pint_input.boundaryExchange(); // fill remote buffer

       out.zero(); 

       if(timeRank!=0)
         out.axpy(1.0/2.0,*pint_input.getRemoteBufferPtr(-1)); 
       out.axpy(1.0/2.0,*pint_input.getVectorPtr(0)); 
     } 

     //   X_{2*i+1} = x_i                            -1 <= i < Np          (includes "virtual variable")
     ////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     for(int k=-1;k<Np;k++) {
       ROL::Vector<Real> & out = *pint_output.getVectorPtr(2*k+1);
       out.set(*pint_input.getVectorPtr(k));
     }

     //   X_{2*i+2} = (x_i + x_{i+1})/2              0 <= i < Np - 1  (averaged variables, created at fine level)
     ////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     for(int k=0;k<Np-1;k++) {
       ROL::Vector<Real> & out = *pint_output.getVectorPtr(2*k+2);

       out.zero(); 
       out.axpy(1.0/2.0,*pint_input.getVectorPtr(k+0)); 
       out.axpy(1.0/2.0,*pint_input.getVectorPtr(k+1)); 
     }
   }

   /**
    * \brief Prolong a control space vector
    *
    * Currently assumes a piecewise constant control
    *
    * \param[in] inputLevel Multigrid level of the input vector (not currently used)
    * \param[in] input The input vector to be "restricted" 
    * \param[in] output The output vector resulting from restriction, must be at level inputLevel+1
    */
   void prolongOptVector(const Vector<Real> & input,Vector<Real> & output)
   {
     const PinTVector<Real> & pint_input  = dynamic_cast<const PinTVector<Real>&>(input);
     PinTVector<Real>       & pint_output = dynamic_cast<PinTVector<Real>&>(output);

     // handle interior
     for(int k=0;k<pint_input.numOwnedSteps();k++) {

       pint_output.getVectorPtr(2*k+0)->set(*pint_input.getVectorPtr(k)); 
       pint_output.getVectorPtr(2*k+1)->set(*pint_input.getVectorPtr(k)); 
     }
   }

   // KKT multigrid preconditioner
   /////////////////////////////////////////////
   
   /**
    * \brief Apply the augmented KKT operator
    *
    * \param[out] output Action of the KKT operator 
    * \param[in] input Input vector to the operator
    * \param[in] u State vector to linearize around
    * \param[in] z Control vector to linearize around
    * \param[in] tol Tolerance
    * \param[in] level Level of the multigrid hierarchy, default is 0 (fine level)        
    */
   void applyAugmentedKKT(Vector<Real> & output, 
                          const Vector<Real> & input,
                          const Vector<Real> & u, 
                          const Vector<Real> & z,
                          Real & tol,
                          int level=0) 
   {
     using PartitionedVector = PartitionedVector<Real>;

     auto part_output = dynamic_cast<PartitionedVector&>(output);
     auto part_input = dynamic_cast<const PartitionedVector&>(input);

     part_output.zero();

     auto output_u = part_output.get(0);
     auto output_z = part_output.get(1);
     auto output_v = part_output.get(2);
     auto output_v_tmp = output_v->clone();

     auto input_u  = part_input.get(0); // state
     auto input_z  = part_input.get(1); // control
     auto input_v  = part_input.get(2); // lagrange multiplier

     // objective
     applyAdjointJacobian_1_leveled(*output_u,*input_v,u,z,tol,level);
     output_u->axpy(1.0,*input_u);

     applyAdjointJacobian_2_leveled(*output_z,*input_v,u,z,tol,level);
     output_z->axpy(1.0,*input_z);

     // constraint
     applyJacobian_1_leveled(*output_v_tmp,*input_u,u,z,tol,level);
     applyJacobian_2_leveled(*output_v,*input_z,u,z,tol,level);

     output_v->axpy(1.0,*output_v_tmp);
   }

   /**
    * \brief Apply an inverse of the KKT system, or an approximation.
    *
    * Compute the inverse of the KKT system. This uses a block factorization
    * and the exact Jacobian. Currently the implementation ignores the constraint
    * Jacobian and looks more like a Wathen style preconditioner.
    *
    * Additionally, if an approximate version is desired then a block jacobi preconditioner
    * split at processor boundaries is used.
    *
    * \param[out] output Action of the inverse KKT matrix
    * \param[in] input Input vector to the inverse operator
    * \param[in] u State vector to linearize around
    * \param[in] z Control vector to linearize around
    * \param[in] tol Tolerance
    * \param[in] approx Apply the parallel block Jacobi inverse if true,
    *                   otherwise do the serial Wathen preconditioner, default false
    * \param[in] level Level of the multigrid hierarchy, default is 0 (fine level)        
    */
   void applyAugmentedInverseKKT(Vector<Real> & output, 
                                 const Vector<Real> & input,
                                 const Vector<Real> & u, 
                                 const Vector<Real> & z,
                                 Real & tol,
                                 bool approx=false,
                                 int level=0) 
   {
     using PartitionedVector = PartitionedVector<Real>;

     auto timer = Teuchos::TimeMonitor::getStackedTimer();

     std::string levelStr = "";
     {
       std::stringstream ss;
       ss << "-" << level;
       levelStr = ss.str();
     }

     timer->start("applyAugmentedInverseKKT"+levelStr); 

     auto part_output = dynamic_cast<PartitionedVector&>(output);
     auto part_input = dynamic_cast<const PartitionedVector&>(input);
 
     part_output.zero();
 
     auto output_u = part_output.get(0);
     auto output_z = part_output.get(1);
     auto output_v = part_output.get(2);
     auto temp_u = output_u->clone();
     auto temp_z = output_z->clone();
     auto temp_v = output_v->clone();
 
     auto input_u  = part_input.get(0);
     auto input_z  = part_input.get(1);
     auto input_v  = part_input.get(2);
 
     temp_u->zero();
     temp_z->zero();
     temp_v->zero();
 
     // [ I   0  J' * inv(J*J') ] [  I        ]
     // [     I  K' * inv(J*J') ] [  0  I     ]
     // [            -inv(J*J') ] [ -J -K  I  ]
    
     // L Factor
     /////////////////////
     temp_u->axpy(1.0,*input_u);
     temp_z->axpy(1.0,*input_z);
 
     // apply -J
     if(not approx)
       applyJacobian_1_leveled(*temp_v,*input_u,u,z,tol,level);
     else
       applyJacobian_1_leveled_approx(*temp_v,*input_u,u,z,tol,level);

     // apply -K ???? (not yet)
 
     temp_v->scale(-1.0);
     temp_v->axpy(1.0,*input_v);
 
     // U Factor
     /////////////////////
     
     // schur complement (Wathen style)
     {
       auto temp = output_v->clone();
       temp->zero();
 
       if(not approx) {
         applyInverseJacobian_1_leveled(*temp, *temp_v, u,z,tol,level);
         applyInverseAdjointJacobian_1_leveled(*output_v,*temp,u,z,tol,level);
       }
       else {
         invertTimeStepJacobian(*temp, *temp_v, u,z,tol,level);
         invertAdjointTimeStepJacobian(*output_v,*temp,u,z,tol,level);
       }
       output_v->scale(-1.0);
     }

     output_z->set(*temp_z);
 
     if(not approx)
       applyAdjointJacobian_1_leveled(*output_u,*output_v,u,z,tol,level); 
     else 
       applyAdjointJacobian_1_leveled_approx(*output_u,*output_v,u,z,tol,level); 
 
     output_u->scale(-1.0);
     output_u->axpy(1.0,*temp_u);

     timer->stop("applyAugmentedInverseKKT"+levelStr); 
   }

   /**
    * \brief Apply a multigrid in time approximate inverse operator to the KKT matrix
    *
    * This uses a multigrid in time algorithm with a parallel domain decomposition block
    * Jacobi smoother for each time sub domain. 
    *
    * \param[out] x Action of the multigrid operator 
    * \param[in] b Input vector to the operator
    * \param[in] u State vector to linearize around
    * \param[in] z Control vector to linearize around
    * \param[in] tol Tolerance
    * \param[in] level Level of the multigrid hierarchy, default is 0 (fine level)        
    */
   void applyMultigridAugmentedKKT(Vector<Real> & x, 
                                const Vector<Real> & b,
                                const Vector<Real> & u, 
                                const Vector<Real> & z,
                                Real & tol,
                                int level=0) 
   {
     using PartitionedVector = PartitionedVector<Real>;

     auto timer = Teuchos::TimeMonitor::getStackedTimer();

     std::string levelStr = "";
     {
       std::stringstream ss;
       ss << "-" << level;
       levelStr = ss.str();
     }

     timer->start("applyMGAugmentedKKT"+levelStr);

     // base case: solve the KKT system directly
     if(level+1==maxLevels_) {
       bool approxSmoother = false;
  
       applyAugmentedInverseKKT(x,b,u,z,tol,approxSmoother,level);

       timer->stop("applyMGAugmentedKKT"+levelStr);
       return;
     }

     bool approxSmoother = true;

     timer->start("applyMGAugmentedKKT-preSmooth");

     auto dx = x.clone();
     auto residual = b.clone();
     residual->set(b);

     // apply one smoother sweep
     for(int i=0;i<numSweeps_;i++) {
       applyAugmentedInverseKKT(*dx,*residual,u,z,tol,approxSmoother,level); 
       x.axpy(omega_,*dx);

       // compute the residual
       applyAugmentedKKT(*residual,x,u,z,tol,level);
       residual->scale(-1.0);
       residual->axpy(1.0,b);
     }
     timer->stop("applyMGAugmentedKKT-preSmooth");

     // solve the coarse system
     timer->start("applyMGAugmentedKKT-coarse");
     {
       auto pint_u = dynamic_cast<const PinTVector<Real>&>(u);
       auto pint_z = dynamic_cast<const PinTVector<Real>&>(z);

       auto dx_u = dynamic_cast<PartitionedVector&>(*dx).get(0);
       auto dx_z = dynamic_cast<PartitionedVector&>(*dx).get(1);
       auto dx_v = dynamic_cast<PartitionedVector&>(*dx).get(2);

       auto residual_u = dynamic_cast<PartitionedVector&>(*residual).get(0);
       auto residual_z = dynamic_cast<PartitionedVector&>(*residual).get(1);
       auto residual_v = dynamic_cast<PartitionedVector&>(*residual).get(2);

       auto crs_u            = allocateSimVector(pint_u,level+1);
       auto crs_z            = allocateOptVector(pint_z,level+1);

       auto crs_residual_u   = allocateSimVector(pint_u,level+1);
       auto crs_residual_z   = allocateOptVector(pint_z,level+1);
       auto crs_residual_v   = allocateSimVector(pint_u,level+1);

       auto crs_correction_u = allocateSimVector(pint_u,level+1);
       auto crs_correction_z = allocateOptVector(pint_z,level+1);
       auto crs_correction_v = allocateSimVector(pint_u,level+1);

       restrictSimVector(*residual_u,*crs_residual_u);
       restrictOptVector(*residual_z,*crs_residual_z);
       restrictSimVector(*residual_v,*crs_residual_v);

       restrictSimVector(u,*crs_u);               // restrict the state to the coarse level
       restrictOptVector(z,*crs_z);               // restrict the control to the coarse level

       auto crs_correction = makePtr<PartitionedVector>({crs_correction_u,crs_correction_z,crs_correction_v});
       auto crs_residual   = makePtr<PartitionedVector>({  crs_residual_u,  crs_residual_z,  crs_residual_v});

       applyMultigridAugmentedKKT(*crs_correction,*crs_residual,*crs_u,*crs_z,tol,level+1);

       prolongSimVector(*crs_correction_u,*dx_u);
       prolongOptVector(*crs_correction_z,*dx_z);
       prolongSimVector(*crs_correction_v,*dx_v);

       x.axpy(1.0,*dx);
     }
     timer->stop("applyMGAugmentedKKT-coarse");

     // apply one smoother sweep
     timer->start("applyMGAugmentedKKT-postSmooth");
     for(int i=0;i<numSweeps_;i++) {
       // compute the residual
       applyAugmentedKKT(*residual,x,u,z,tol,level);
       residual->scale(-1.0);
       residual->axpy(1.0,b);

       applyAugmentedInverseKKT(*dx,*residual,u,z,tol,approxSmoother,level); 
       x.axpy(omega_,*dx);
     }
     timer->stop("applyMGAugmentedKKT-postSmooth");

     timer->stop("applyMGAugmentedKKT"+levelStr);
   }

}; // ROL::PinTConstraint


} // namespace ROL 

#endif // ROL_PINTCONSTRAINT_HPP

