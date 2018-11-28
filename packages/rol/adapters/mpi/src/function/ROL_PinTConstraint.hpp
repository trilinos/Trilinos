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
#include "ROL_SerialConstraint.hpp"
#include "ROL_PinTCommunicationUtilities.hpp"

#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StackedTimer.hpp"

namespace ROL {

/** This helper method builds a pint "state" vector for use in the 
    PinTConstraint class. 

    \param[in] communicators Structure with objects required for parallel-in-time communication.
                             (hint: this can be used to properly allocate the state and control vectors)
    \param[in] vectorComm Communication mechansim for sharing vectors between processors.
    \param[in] steps Number of steps 
    \param[in] localVector Spatial vector for a single time step.
  */
template <class Real>
Ptr<PinTVector<Real>>
buildStatePinTVector(const Ptr<const PinTCommunicators> & communicators,
                     const Ptr<const PinTVectorCommunication<Real>> & vectorComm,
                     int steps,
                     const Ptr<Vector<Real>> & localVector)
{ std::vector<int> stencil = {-1,0};
  const int replicate = 2; // use virtual variables
  return makePtr<PinTVector<Real>>(communicators,vectorComm,localVector,steps,stencil,replicate); }

/** This helper method builds a pint "control" vector for use in the 
    PinTConstraint class. 

    \param[in] communicators Structure with objects required for parallel-in-time communication.
                             (hint: this can be used to properly allocate the state and control vectors)
    \param[in] vectorComm Communication mechansim for sharing vectors between processors.
    \param[in] steps Number of steps 
    \param[in] localVector Spatial vector for a single time step.
  */
template <class Real>
Ptr<PinTVector<Real>>
buildControlPinTVector(const Ptr<const PinTCommunicators> & communicators,
                       const Ptr<const PinTVectorCommunication<Real>> & vectorComm,
                       int steps,
                       const Ptr<Vector<Real>> & localVector)
{ std::vector<int> stencil = {0};
  return makePtr<PinTVector<Real>>(communicators,vectorComm,localVector,steps,stencil); }

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

  Ptr<const std::vector<TimeStamp<Real>>>         userTimeStamps_;    // these are the untouched ones the user passes in
  Ptr<std::vector<TimeStamp<Real>>>               timeStamps_;        // these are used internally
  std::vector<Ptr<std::vector<TimeStamp<Real>>>>  stamps_;            // these are used internally on each level

  // preconditioner settings
  bool applyMultigrid_;         // default to block jacobi
  int maxLevels_;               // must turn on multigrid 

  int numSweeps_;
  int numCoarseSweeps_;
  Real omega_;
  Real omegaCoarse_;

  Real globalScale_;

  std::vector<Ptr<const PinTCommunicators>> communicators_;
  Ptr<const PinTVectorCommunication<Real>> vectorComm_; // object for sending whole vectors between processors

  bool useRepart_;

  // Get the state vector required by the serial constraint object
  Ptr<Vector<Real>> getStateVector(const PinTVector<Real> & src,int s=0)
  {
    // notice this ignores the stencil for now
   
    /*
    std::vector<Ptr<Vector<Real>>> vecs;

    vecs.push_back(src.getVectorPtr(-1));  // inject the initial condition

    for(int i=0;i<src.numOwnedSteps();i+=2) { // ignoring the stencil
      vecs.push_back(src.getVectorPtr(i+1)); 
    }

    return makePtr<PartitionedVector<Real>>(vecs);
    */

    std::vector<Ptr<Vector<Real>>> vecs;
    vecs.push_back(src.getVectorPtr(2*s-1)); // virtual value lives on the odds (and negatives)
    vecs.push_back(src.getVectorPtr(2*s));   // true value lives on the evens (the vector ends on the evens)

    return makePtr<PartitionedVector<Real>>(vecs);
  }

  // Get the control vector required by the serial constraint object
  Ptr<Vector<Real>> getControlVector(const PinTVector<Real> & src,int s=0)
  {
    // notice this ignores the stencil for now
 
    /*
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
    */

    std::vector<Ptr<Vector<Real>>> vecs;
    vecs.push_back(src.getVectorPtr(0)->clone());
      // FIXME-A: this control index is never used. This stems from not wanting to set
      //          an intial condition for the test function. But instead use the 
      //          "setSkipInitialCondition" option. This all needs to be rethought. 
      //          This is linked to the FIXME's in "getStateConstraint"
      //             -- ECC, June 11, 2018
    vecs.push_back(src.getVectorPtr(s));
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
    , numCoarseSweeps_(1)
    , omega_(2.0/3.0)
    , omegaCoarse_(1.0)
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
    , useRepart_(false)
  { 
    initialize(stepConstraint,initialCond,timeStamps);
  }

  /**
   * Set sweeps for jacobi and multigrid
   */
  void setSweeps(int s)
  { numSweeps_ = s; }

  /**
   * Set sweeps for coarse level solver
   */
  void setCoarseSweeps(int s)
  { numCoarseSweeps_ = s; }

  /**
   * Set relaxation parater for jacobi and multigrid
   */
  void setRelaxation(Real o)
  { omega_ = o; }

  /**
   * Set relaxation parater for jacobi and multigrid
   */
  void setCoarseRelaxation(Real o)
  { omegaCoarse_ = o; }

  /**
   * Set the global scaling for the coupling constraints.
   */
  void setGlobalScale(Real o)
  { globalScale_ = o; }

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
  getSerialConstraint(const Ptr<ROL::Vector<Real>> & ui,int level,int step=0)
  {
    typedef Ptr<std::vector<TimeStamp<Real>>> VecPtr;

    // should you access this by repartioned scheme?
    VecPtr timeStamps;
    if(useRepart_)
      timeStamps = getTimeStampsByLevel_repart(level);
    else
      timeStamps = getTimeStampsByLevel(level);

    // build a singl epoint time stamp
    VecPtr localTimeStamps = makePtr<std::vector<TimeStamp<Real>>>(2);
    localTimeStamps->at(0) = timeStamps->at(step);
    localTimeStamps->at(1) = timeStamps->at(step+1);

    if(ROL::is_nullPtr(timeDomainConstraint_)) {

      timeDomainConstraint_ = ROL::makePtr<SerialConstraint<Real>>(stepConstraint_,*ui,localTimeStamps);
        // FIXME-A: Change the Nt+1 back to Nt...
    }
    else {
      timeDomainConstraint_->setTimeStampsPtr(localTimeStamps);
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

    int timeRank = pint_u.communicators().getTimeRank();
    const std::vector<int> & stencil = pint_c.stencil();

    // assert a forward stencil
    assert(stencil.size()==2);
    assert(stencil[0]==-1);
    assert(stencil[1]== 0);

    // v0 = u0
    if(timeRank==0)
      pint_u.getVectorPtr(-1)->set(*initialCond_);
    else
      pint_u.getVectorPtr(-1)->set(*pint_u.getRemoteBufferPtr(-1));

    size_t numSteps = getTimeStampsByLevel(0)->size()-1;
    for(size_t s=0;s<numSteps;s++) { // num time steps == num time stamps-1 

      auto part_c = getStateVector(pint_c,s);   
      auto part_u = getStateVector(pint_u,s);   
      auto part_z = getControlVector(pint_z,s); 

      // compute the constraint for this subdomain
      auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s-1),level,s);

      constraint->solve(*part_c,*part_u,*part_z,tol);

      // satisfy the time continuity constraint, note that this just using the identity matrix
      if(s<numSteps-1) {
        pint_u.getVectorPtr(2*s+1)->set(*pint_u.getVectorPtr(2*s));     // assign true value to virtual value
      }
    }

    pint_u.boundaryExchange(PinTVector<Real>::SEND_ONLY);

    // c.scale(-10000.0);

    // call the original implementation
    value(c,u,z,tol);

    // std::cout << "U NORM = " << u.norm() << std::endl;
    // std::cout << "Z NORM = " << z.norm() << std::endl;
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
    const std::vector<int> & stencil = pint_c.stencil();

    // assert a forward stencil
    assert(stencil.size()==2);
    assert(stencil[0]==-1);
    assert(stencil[1]== 0);

    // v0 - u0
    pint_c.getVectorPtr(-1)->set(*pint_u.getVectorPtr(-1)); 
    if(timeRank==0)
      pint_c.getVectorPtr(-1)->axpy(-1,*initialCond_);
    else
      pint_c.getVectorPtr(-1)->axpy(-1,*pint_u.getRemoteBufferPtr(-1));
    pint_c.getVectorPtr(-1)->scale(globalScale_);

    size_t numSteps = getTimeStampsByLevel(0)->size()-1;
    for(size_t s=0;s<numSteps;s++) { // num time steps == num time stamps-1 

      auto part_c = getStateVector(pint_c,s);   // strip out initial condition/constraint
      auto part_u = getStateVector(pint_u,s);   // strip out initial condition/constraint
      auto part_z = getControlVector(pint_z,s); 
  
      // compute the constraint for this subdomain
      //    u_s = F(v_{s-1},z_s)
      auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s-1),level,s);
  
      constraint->value(*part_c,*part_u,*part_z,tol);
  
      // build in the time continuity constraint, note that this is just using the identity matrix
      //    v_s = u_s
      if(s<numSteps-1) {
        pint_c.getVectorPtr(2*s+1)->set(*pint_u.getVectorPtr(2*s+1));     // this is the virtual value
        pint_c.getVectorPtr(2*s+1)->axpy(-1.0,*pint_u.getVectorPtr(2*s)); // this is the u value
        pint_c.getVectorPtr(2*s+1)->scale(globalScale_);
      }
    }
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

    // communicate neighbors, these are blocking calls
    pint_v.boundaryExchange();
    pint_u.boundaryExchange();
    pint_z.boundaryExchange();

    int timeRank = pint_u.communicators().getTimeRank();
    const std::vector<int> & stencil = pint_jv.stencil();

    // assert a forward stencil
    assert(stencil.size()==2);
    assert(stencil[0]==-1);
    assert(stencil[1]== 0);

    // v0 - u0
    pint_jv.getVectorPtr(-1)->set(*pint_v.getVectorPtr(-1)); 
    if(timeRank!=0)
      pint_jv.getVectorPtr(-1)->axpy(-1,*pint_v.getRemoteBufferPtr(-1));
    pint_jv.getVectorPtr(-1)->scale(globalScale_);

    size_t numSteps = getTimeStampsByLevel(level)->size()-1;
    for(size_t s=0;s<numSteps;s++) { // num time steps == num time stamps-1 

      auto part_jv = getStateVector(pint_jv,s); 
      auto part_v  = getStateVector(pint_v,s);    
      auto part_u  = getStateVector(pint_u,s);   
      auto part_z  = getControlVector(pint_z,s); 
  
      // compute the constraint for this subdomain
      //    u_s = F(v_{s-1},z_s)
      auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s-1),level,s);
  
      constraint->applyJacobian_1(*part_jv,*part_v,*part_u,*part_z,tol);
  
      // build in the time continuity constraint, note that this is just using the identity matrix
      //    v_s = u_s
      if(s<numSteps-1) {
        pint_jv.getVectorPtr(2*s+1)->set(*pint_v.getVectorPtr(2*s+1));     // this is the virtual value
        pint_jv.getVectorPtr(2*s+1)->axpy(-1.0,*pint_v.getVectorPtr(2*s)); // this is the u value
        pint_jv.getVectorPtr(2*s+1)->scale(globalScale_);
      }
    }
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

    // communicate neighbors, these are blocking calls
    pint_v.boundaryExchange();
    pint_u.boundaryExchange();
    pint_z.boundaryExchange();

    // int timeRank = pint_u.communicators().getTimeRank(); // Unused
    const std::vector<int> & stencil = pint_jv.stencil();

    // assert a forward stencil
    assert(stencil.size()==2);
    assert(stencil[0]==-1);
    assert(stencil[1]== 0);

    // v0 - u0
    pint_jv.getVectorPtr(-1)->set(*pint_v.getVectorPtr(-1)); 
    pint_jv.getVectorPtr(-1)->scale(globalScale_);

    size_t numSteps = getTimeStampsByLevel(level)->size()-1;
    for(size_t s=0;s<numSteps;s++) { // num time steps == num time stamps-1 

      auto part_jv = getStateVector(pint_jv,s); 
      auto part_v  = getStateVector(pint_v,s);    
      auto part_u  = getStateVector(pint_u,s);   
      auto part_z  = getControlVector(pint_z,s); 
  
      // compute the constraint for this subdomain
      //    u_s = F(v_{s-1},z_s)
      auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s-1),level,s);
  
      constraint->applyJacobian_1(*part_jv,*part_v,*part_u,*part_z,tol);
  
      // build in the time continuity constraint, note that this is just using the identity matrix
      //    v_s = u_s
      if(s<numSteps-1) {
        pint_jv.getVectorPtr(2*s+1)->set(*pint_v.getVectorPtr(2*s+1));     // this is the virtual value
        pint_jv.getVectorPtr(2*s+1)->axpy(-1.0,*pint_v.getVectorPtr(2*s)); // this is the u value
        pint_jv.getVectorPtr(2*s+1)->scale(globalScale_);
      }
    }
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
    assert(pint_z.numOwnedSteps()==pint_v.numOwnedSteps());

    // communicate neighbors, these are block calls
    pint_v.boundaryExchange();
    pint_u.boundaryExchange();
    pint_z.boundaryExchange();
//  Unused
//    int timeRank = pint_u.communicators().getTimeRank();
    const std::vector<int> & stencil = pint_jv.stencil();

    // assert a forward stencil
    assert(stencil.size()==2);
    assert(stencil[0]==-1);
    assert(stencil[1]== 0);

    // v0 - u0
    pint_jv.getVectorPtr(-1)->zero();

    size_t numSteps = getTimeStampsByLevel(level)->size()-1;
    for(size_t s=0;s<numSteps;s++) { // num time steps == num time stamps-1 

      auto part_jv = getStateVector(pint_jv,s); 
      auto part_v  = getControlVector(pint_v,s);    
      auto part_u  = getStateVector(pint_u,s);   
      auto part_z  = getControlVector(pint_z,s); 
  
      // compute the constraint for this subdomain
      //    u_s = F(v_{s-1},z_s)
      auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s-1),level,s);
  
      constraint->applyJacobian_2(*part_jv,*part_v,*part_u,*part_z,tol);
  
      // build in the time continuity constraint, note that this is just using the identity matrix
      //    v_s = u_s
      if(s<numSteps-1) {
        pint_jv.getVectorPtr(2*s+1)->zero();
      }
    }
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

    int timeRank = pint_u.communicators().getTimeRank();
    const std::vector<int> & stencil = pint_ijv.stencil();

    // assert a forward stencil
    assert(stencil.size()==2);
    assert(stencil[0]==-1);
    assert(stencil[1]== 0);

    // initial constraint
    if(timeRank==0) {
      // update the old time information
      pint_ijv.getVectorPtr(-1)->set(*pint_v.getVectorPtr(-1));
      pint_ijv.getVectorPtr(-1)->scale(1.0/globalScale_);
    }
    else {
      pint_ijv.getVectorPtr(-1)->set(*pint_ijv.getRemoteBufferPtr(-1));    
      pint_ijv.getVectorPtr(-1)->axpy(1.0/globalScale_,*pint_v.getVectorPtr(-1));        
    }

    size_t numSteps = getTimeStampsByLevel(level)->size()-1;
    for(size_t s=0;s<numSteps;s++) { // num time steps == num time stamps-1 

      auto part_ijv = getStateVector(pint_ijv,s);   
      auto part_v   = getStateVector(pint_v,s);   
      auto part_u   = getStateVector(pint_u,s);   
      auto part_z   = getControlVector(pint_z,s); 

      // compute the constraint for this subdomain
      auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s-1),level,s);

      constraint->applyInverseJacobian_1(*part_ijv,*part_v,*part_u,*part_z,tol);

      // satisfy the time continuity constraint, note that this just using the identity matrix
      if(s<numSteps-1) {
        pint_ijv.getVectorPtr(2*s+1)->set(*pint_ijv.getVectorPtr(2*s));    
        pint_ijv.getVectorPtr(2*s+1)->axpy(1.0/globalScale_,*pint_v.getVectorPtr(2*s+1));        
      }
    }

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

     // int timeRank = pint_u.communicators().getTimeRank();
     const std::vector<int> & stencil = pint_ajv.stencil();

     // assert a forward stencil
     assert(stencil.size()==2);
     assert(stencil[0]==-1);
     assert(stencil[1]== 0);

     int numSteps = Teuchos::as<int>(getTimeStampsByLevel(level)->size())-1;
     for(int s=numSteps-1;s>=0;s--) { // num time steps == num time stamps-1 
       auto part_ajv = getStateVector(pint_ajv,s);
       auto part_v   = getStateVector(pint_v,s);   
       auto part_u   = getStateVector(pint_u,s);   
       auto part_z   = getControlVector(pint_z,s); 

       // handle the constraint adjoint: Required for correctly applying step adjoint
       if(s<numSteps-1)
         pint_ajv.getVectorPtr(2*s+1)->axpy(globalScale_,*pint_v.getVectorPtr(2*s+1));

       // compute the constraint for this subdomain
       auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s-1),level,s);

       constraint->applyAdjointJacobian_1(*part_ajv,*part_v,*part_u,*part_z,*part_v,tol);

       // this is the remainder of the constraint application
       if(s<numSteps-1)
         pint_ajv.getVectorPtr(2*s)->axpy(-globalScale_,*pint_v.getVectorPtr(2*s+1));
     }

     pint_ajv.getVectorPtr(-1)->axpy(globalScale_,*pint_v.getVectorPtr(-1));
     pint_ajv.getRemoteBufferPtr(-1)->set(*pint_v.getVectorPtr(-1)); // this will be sent to the remote processor
     pint_ajv.getRemoteBufferPtr(-1)->scale(-globalScale_);

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

     // int timeRank = pint_u.communicators().getTimeRank();
     const std::vector<int> & stencil = pint_ajv.stencil();

     // assert a forward stencil
     assert(stencil.size()==2);
     assert(stencil[0]==-1);
     assert(stencil[1]== 0);

     int numSteps = Teuchos::as<int>(getTimeStampsByLevel(level)->size())-1;
     for(int s=numSteps-1;s>=0;s--) { // num time steps == num time stamps-1 
       auto part_ajv = getStateVector(pint_ajv,s);
       auto part_v   = getStateVector(pint_v,s);   
       auto part_u   = getStateVector(pint_u,s);   
       auto part_z   = getControlVector(pint_z,s); 

       // handle the constraint adjoint: Required for correctly applying step adjoint
       if(s<numSteps-1)
         pint_ajv.getVectorPtr(2*s+1)->axpy(globalScale_,*pint_v.getVectorPtr(2*s+1));

       // compute the constraint for this subdomain
       auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s-1),level,s);

       constraint->applyAdjointJacobian_1(*part_ajv,*part_v,*part_u,*part_z,*part_v,tol);

       // this is the remainder of the constraint application
       if(s<numSteps-1)
         pint_ajv.getVectorPtr(2*s)->axpy(-globalScale_,*pint_v.getVectorPtr(2*s+1));
     }

     pint_ajv.getVectorPtr(-1)->axpy(globalScale_,*pint_v.getVectorPtr(-1));
     // pint_ajv.getRemoteBufferPtr(-1)->set(*pint_v.getVectorPtr(-1)); // this will be sent to the remote processor
     // pint_ajv.getRemoteBufferPtr(-1)->scale(-globalScale_);
     pint_ajv.getRemoteBufferPtr(-1)->zero();

     // this sums from the remote buffer into the local buffer which is part of the vector
     // pint_ajv.boundaryExchangeSumInto();
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
    assert(pint_ajv.numOwnedSteps()==pint_z.numOwnedSteps());

    const std::vector<int> & stencil = pint_u.stencil();

    // assert a forward stencil
    assert(stencil.size()==2);
    assert(stencil[0]==-1);
    assert(stencil[1]== 0);

    // we need to make sure this has all zeros to begin with (this includes boundary exchange components)
    pint_ajv.zeroAll();

    // communicate neighbors, these are block calls
    pint_v.boundaryExchange();
    pint_u.boundaryExchange();
    pint_z.boundaryExchange();

    int numSteps = Teuchos::as<int>(getTimeStampsByLevel(level)->size())-1;
    for(int s=numSteps-1;s>=0;s--) { // num time steps == num time stamps-1 
  
      auto part_ajv = getControlVector(pint_ajv,s);
      auto part_v   = getStateVector(pint_v,s);   
      auto part_u   = getStateVector(pint_u,s); 
      auto part_z   = getControlVector(pint_z,s); 

      // compute the constraint for this subdomain
      auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s-1),level,s);

      constraint->applyAdjointJacobian_2(*part_ajv,*part_v,*part_u,*part_z,tol);
    }

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
     PinTVector<Real> & pint_v          = dynamic_cast<PinTVector<Real>&>(*v_copy);
     const PinTVector<Real> & pint_u    = dynamic_cast<const PinTVector<Real>&>(u);
     const PinTVector<Real> & pint_z    = dynamic_cast<const PinTVector<Real>&>(z);
     // its possible we won't always want to cast to a PinT vector here
     //

     int timeRank = pint_u.communicators().getTimeRank();
     int timeSize = pint_u.communicators().getTimeSize();
     const std::vector<int> & stencil = pint_iajv.stencil();

     // assert a forward stencil
     assert(stencil.size()==2);
     assert(stencil[0]==-1);
     assert(stencil[1]== 0);

     assert(pint_iajv.numOwnedSteps()==pint_u.numOwnedSteps());
     assert(pint_iajv.numOwnedSteps()==pint_v.numOwnedSteps());

     pint_iajv.zero();

     pint_z.boundaryExchange();
     pint_v.boundaryExchange();
     pint_u.boundaryExchange();
     pint_iajv.boundaryExchangeSumInto(PinTVector<Real>::RECV_ONLY); // this is going to block

     int numSteps = Teuchos::as<int>(getTimeStampsByLevel(level)->size())-1;

     // on the last processor (with the last time step) this doesn't happen
     if(timeRank+1 < timeSize) {
       int s = numSteps-1;
       pint_v.getVectorPtr(2*s)->axpy(1.0,*pint_iajv.getVectorPtr(2*s));
         // here we are abusing the iajv vector and storing the RHS to pass
         // between processors, this will be overwritten below
     }

     for(int s=numSteps-1;s>=0;s--) { // num time steps == num time stamps-1 
       auto part_iajv = getStateVector(pint_iajv,s);   
       auto part_v    = getStateVector(pint_v,s);   
       auto part_u    = getStateVector(pint_u,s);   
       auto part_z    = getControlVector(pint_z,s); 

       if(s<numSteps-1) {
         pint_iajv.getVectorPtr(2*s+1)->scale(1.0/globalScale_);
         pint_v.getVectorPtr(2*s)->axpy(globalScale_,*pint_iajv.getVectorPtr(2*s+1));
       }

       // compute the constraint for this subdomain
       auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s-1),level,s);

       constraint->applyInverseAdjointJacobian_1(*part_iajv,*part_v,*part_u,*part_z,tol);
     }

     pint_iajv.getRemoteBufferPtr(-1)->set(*pint_iajv.getVectorPtr(-1));
         // here we are abusing the iajv vector and storing the RHS to pass between processors
     pint_iajv.getVectorPtr(-1)->scale(1.0/globalScale_);

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
   void invertTimeStepJacobian(PinTVector<Real>       & pint_ijv,
                               const PinTVector<Real> & pint_v,
                               const PinTVector<Real> & pint_u,
                               const PinTVector<Real> & pint_z,
                               Real &tol,
                               int level) 
   {
//  Unused
//     int timeRank = pint_u.communicators().getTimeRank();

     // update the old time information
     pint_ijv.getVectorPtr(-1)->set(*pint_v.getVectorPtr(-1));
     pint_ijv.getVectorPtr(-1)->scale(1.0/globalScale_);

     size_t numSteps = getTimeStampsByLevel(level)->size()-1;
     for(size_t s=0;s<numSteps;s++) { // num time steps == num time stamps-1 

       auto part_ijv = getStateVector(pint_ijv,s);   
       auto part_v   = getStateVector(pint_v,s);   
       auto part_u   = getStateVector(pint_u,s);   
       auto part_z   = getControlVector(pint_z,s); 

       // compute the constraint for this subdomain
       auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s-1),level,s);

       constraint->applyInverseJacobian_1(*part_ijv,*part_v,*part_u,*part_z,tol);

       // satisfy the time continuity constraint, note that this just using the identity matrix
       if(s<numSteps-1) {
         pint_ijv.getVectorPtr(2*s+1)->set(*pint_ijv.getVectorPtr(2*s));    
         pint_ijv.getVectorPtr(2*s+1)->axpy(1.0/globalScale_,*pint_v.getVectorPtr(2*s+1));        
       }
     }
   }
 
   // Done in parallel, no blocking, solve a linear system on this processor
   void invertAdjointTimeStepJacobian(PinTVector<Real>       & pint_iajv,
                                      const PinTVector<Real> & pint_v_src,
                                      const PinTVector<Real> & pint_u,
                                      const PinTVector<Real> & pint_z,
                                      Real &tol,
                                      int level)
   {
     // this is an inefficient hack (see line below where pint_v is modified!!!!)
     ///////////////////////////////////
     auto v_copy = pint_v_src.clone();
     v_copy->set(pint_v_src);
     PinTVector<Real> & pint_v = dynamic_cast<PinTVector<Real>&>(*v_copy);
     ///////////////////////////////////
     int numSteps = Teuchos::as<int>(getTimeStampsByLevel(level)->size())-1;

     for(int s=numSteps-1;s>=0;s--) { // num time steps == num time stamps-1 
       auto part_iajv = getStateVector(pint_iajv,s);   
       auto part_v    = getStateVector(pint_v,s);   
       auto part_u    = getStateVector(pint_u,s);   
       auto part_z    = getControlVector(pint_z,s); 

       if(s<numSteps-1) {
         pint_iajv.getVectorPtr(2*s+1)->scale(1.0/globalScale_);
         pint_v.getVectorPtr(2*s)->axpy(globalScale_,*pint_iajv.getVectorPtr(2*s+1));
       }

       // compute the constraint for this subdomain
       auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s-1),level,s);

       constraint->applyInverseAdjointJacobian_1(*part_iajv,*part_v,*part_u,*part_z,tol);
     }

     pint_iajv.getVectorPtr(-1)->scale(1.0/globalScale_);
   }
 
   virtual void applyPreconditioner(Vector<Real> &pv,
                                    const Vector<Real> &v,
                                    const Vector<Real> &x,
                                    const Vector<Real> &g,
                                    Real &tol) override
   {
     assert(false);
   }

   // restriction and prolongation functions
   ///////////////////////////////////////////////////////////////////////////////////////////
   
   ROL::Ptr<const PinTCommunicators> getLevelCommunicators(int level) const
   { return communicators_[level]; }

   Ptr<std::vector<TimeStamp<Real>>> getTimeStampsByLevel_repart(int level) const
   {
     return stamps_[level];
   }
   
   /** 
    * \brief Builds up communicators and time stamps for each level
    *
    * This currently doubles the size of the time step as it coarsens. This calls
    * the recursive helper function buildLevelDataStructures to do most of the work
    *
    * \param[in] level_0_ref Reference vector that is used to pull out all the communication
    *                        devices.
    */
   void buildLevels(const Vector<Real> & level_0_ref)
   {
     const PinTVector<Real> & pint_ref  = dynamic_cast<const PinTVector<Real>&>(level_0_ref);

     // set the vector communicator
     vectorComm_ = pint_ref.vectorCommunicationPtr();

     // we allocate the communicators with null
     communicators_.resize(maxLevels_);
     communicators_[0] = pint_ref.communicatorsPtr();

     // we allocate the time stamps vectors with null
     stamps_.resize(maxLevels_);
     stamps_[0] = timeStamps_;

     buildLevelDataStructures(1);

     useRepart_ = true;
   }

   //! Recursive function that builds up the communicators and time stamps on coarse grids
   void buildLevelDataStructures(int level)
   {
     // protect sanity
     assert(level>0);

     // don't go too deep (base case)
     if(level>=maxLevels_)
       return;

     // this processor is no longer participating (base case)
     if(ROL::is_nullPtr(communicators_[level-1]))
       return;

     // this will subdivide the communicators by two
     communicators_[level] = communicators_[level-1]->buildCoarseCommunicators();

     // we now build the satmps to be distributed to other processors, we are halfing the size
     auto & fineStamps = *stamps_[level-1];
     std::vector<ROL::TimeStamp<Real>> sourceStamps((fineStamps.size()-1)/2); 
     for(size_type i=1;i<fineStamps.size();i+=2) {
       auto & stamp = sourceStamps[(i-1)/2];

       // buld a stamp skipping one step
       stamp.t.resize(2);
       stamp.t[0] = fineStamps[i].t[0];
       stamp.t[1] = fineStamps[i+1].t[1];
     }

     // export to coarse distribution always requires a reference to the stamps, however that is only filled
     // if its on the right processor. The builtStamps states if that vector was filled.
     ROL::Ptr<std::vector<ROL::TimeStamp<Real>>> stamps = ROL::makePtr<std::vector<ROL::TimeStamp<Real>>>();
     bool builtStamps = PinT::exportToCoarseDistribution_TimeStamps(sourceStamps,*stamps,*communicators_[level-1],0);

     // setup false step (this is to ensure the control is properly handled
     if(builtStamps) {
       stamps_[level] = stamps;

       Real ta = stamps_[level]->at(0).t.at(0);
       Real tb = stamps_[level]->at(0).t.at(1);
       Real dt =  tb-ta;
 
       // the serial constraint should never see this!
       ROL::TimeStamp<Real> stamp; 
       stamp.t.resize(2);
       stamp.t.at(0) = ta-dt;
       stamp.t.at(1) = ta;

       stamps_[level]->insert(stamps_[level]->begin(),stamp);
     }

     buildLevelDataStructures(level+1);
   }

   /**
    * \brief Allocate a simulation space vector at a particular multigrid level using repartioning
    *         
    * Here repartioning means that the coarse levels are on fewer processors.
    */
   Ptr<Vector<Real>> allocateSimVector_repart(const Vector<Real> & level_0_ref,int level) const
   {
     const PinTVector<Real> & pint_ref  = dynamic_cast<const PinTVector<Real>&>(level_0_ref);
     auto comm = communicators_[level];

     if(ROL::is_nullPtr(comm)) {
       // no vector exists on this processor at this level
       return ROL::nullPtr;
     } else {
       Ptr<std::vector<TimeStamp<Real>>> stamps  = getTimeStampsByLevel_repart(level);
       return buildStatePinTVector(comm,vectorComm_,comm->getTimeSize()*(stamps->size()-1),pint_ref.getVectorPtr(0)->clone());
     }
   }

   Ptr<Vector<Real>> allocateSimVector_local(const Vector<Real> & level_0_ref,int level) const
   {
     const PinTVector<Real> & pint_ref  = dynamic_cast<const PinTVector<Real>&>(level_0_ref);
     auto comm = communicators_[level-1];

     if(ROL::is_nullPtr(comm)) {
       // no vector exists on this processor at this level
       return ROL::nullPtr;
     } else {
       Ptr<std::vector<TimeStamp<Real>>> stamps  = getTimeStampsByLevel_repart(level-1);
       return buildStatePinTVector(comm,vectorComm_,comm->getTimeSize()*(stamps->size()-1)/2,pint_ref.getVectorPtr(0)->clone());
     }
   }

   Ptr<Vector<Real>> allocateSimVector_local_fine(const Vector<Real> & level_0_ref,int level) const
   {
     const PinTVector<Real> & pint_ref  = dynamic_cast<const PinTVector<Real>&>(level_0_ref);
     auto comm = communicators_[level];

     if(ROL::is_nullPtr(comm)) {
       // no vector exists on this processor at this level
       return ROL::nullPtr;
     } else {
       Ptr<std::vector<TimeStamp<Real>>> stamps  = getTimeStampsByLevel_repart(level);
       return buildStatePinTVector(comm,vectorComm_,comm->getTimeSize()*(stamps->size()-1)*2,pint_ref.getVectorPtr(0)->clone());
     }
   }

   /**
    * \brief Allocate a simulation space vector at a particular multigrid level using repartitioning.
    *         
    * Here repartioning means that the coarse levels are on fewer processors.
    */
   Ptr<Vector<Real>> allocateOptVector_repart(const Vector<Real> & level_0_ref,int level) const
   {
     const PinTVector<Real> & pint_ref  = dynamic_cast<const PinTVector<Real>&>(level_0_ref);
     auto comm = communicators_[level];
     
     if(ROL::is_nullPtr(comm)) {
       // no vector exists on this processor at this level
       return ROL::nullPtr;
     } else {
       // no vector exists on this processor at this level
       Ptr<std::vector<TimeStamp<Real>>> stamps  = getTimeStampsByLevel_repart(level);
       return buildControlPinTVector(comm,vectorComm_,comm->getTimeSize()*(stamps->size()-1),pint_ref.getVectorPtr(0)->clone());
     }
   }

   /**
    * \brief Allocate a simulation space vector at a particular multigrid level using repartitioning.
    *         
    * Here repartioning means that the coarse levels are on fewer processors.
    */
   Ptr<Vector<Real>> allocateOptVector_local(const Vector<Real> & level_0_ref,int level) const
   {
     const PinTVector<Real> & pint_ref  = dynamic_cast<const PinTVector<Real>&>(level_0_ref);
     auto comm = communicators_[level-1];
     
     if(ROL::is_nullPtr(comm)) {
       // no vector exists on this processor at this level
       return ROL::nullPtr;
     } else {
       // no vector exists on this processor at this level
       Ptr<std::vector<TimeStamp<Real>>> stamps  = getTimeStampsByLevel_repart(level-1);
       return buildControlPinTVector(comm,vectorComm_,comm->getTimeSize()*(stamps->size()-1)/2,pint_ref.getVectorPtr(0)->clone());
     }
   }

   /**
    * \brief Allocate a simulation space vector at a particular multigrid level using repartitioning.
    *         
    * Here repartioning means that the coarse levels are on fewer processors.
    */
   Ptr<Vector<Real>> allocateOptVector_fine(const Vector<Real> & level_0_ref,int level) const
   {
     const PinTVector<Real> & pint_ref  = dynamic_cast<const PinTVector<Real>&>(level_0_ref);
     auto comm = communicators_[level];
     
     if(ROL::is_nullPtr(comm)) {
       // no vector exists on this processor at this level
       return ROL::nullPtr;
     } else {
       // no vector exists on this processor at this level
       Ptr<std::vector<TimeStamp<Real>>> stamps  = getTimeStampsByLevel_repart(level);
       return buildControlPinTVector(comm,vectorComm_,comm->getTimeSize()*(stamps->size()-1)*2,pint_ref.getVectorPtr(0)->clone());
     }
   }
   
   /**
    * \brief Allocate a simulation space vector at a particular multigrid level.
    */
   Ptr<Vector<Real>> allocateSimVector(const Vector<Real> & level_0_ref,int level) const
   {
     // this is to unify the code
     if(useRepart_) 
       return allocateSimVector_repart(level_0_ref,level);

     const PinTVector<Real> & pint_ref  = dynamic_cast<const PinTVector<Real>&>(level_0_ref);
     auto comm = pint_ref.communicatorsPtr();
     auto vectorComm = pint_ref.vectorCommunicationPtr();
     
     Ptr<std::vector<TimeStamp<Real>>> stamps  = getTimeStampsByLevel(level);
     return buildStatePinTVector(comm,vectorComm,comm->getTimeSize()*(stamps->size()-1),pint_ref.getVectorPtr(0)->clone());
   }

   /**
    * \brief Allocate a simulation space vector at a particular multigrid level.
    */
   Ptr<Vector<Real>> allocateOptVector(const Vector<Real> & level_0_ref,int level) const
   {
     if(useRepart_) 
       return allocateOptVector_repart(level_0_ref,level);

     const PinTVector<Real> & pint_ref  = dynamic_cast<const PinTVector<Real>&>(level_0_ref);
     auto comm = pint_ref.communicatorsPtr();
     auto vectorComm = pint_ref.vectorCommunicationPtr();
     
     Ptr<std::vector<TimeStamp<Real>>> stamps  = getTimeStampsByLevel(level);
     return buildControlPinTVector(comm,vectorComm,comm->getTimeSize()*(stamps->size()-1),pint_ref.getVectorPtr(0)->clone());
   }
   
   /**
    * \brief Restrict a simulation space vector
    *
    * Currently doing the following. Let the coarse distribution be defined locally
    * to contain Np steps. The fine discretization virtual variable is denoted by a V, while the coarse
    * discretization is denoted by a v (similarly U and u for the primary variable). A superscript indicates the relative processor index.
    * If no index is included then it is assumed to be 0 (denoting no relative change in
    * processor index). The restriction for the sim vector thus computes
    *
    *   v_i = V_{2*i} + V_{2*i+1} + V_{2*i+2}      0 <= i < Np - 1
    *   u_i = U_{2*i} + U_{2*i+1} + U_{2*i+2}      0 <= i < Np - 1
    *
    * This ensures that for all interior points the coarse grid is the sum of the three neighbors. These
    * averages are dividecd by three.
    * At the boundaries of processors the same thing holds, but with the appropriate sharing. Finally,
    * at the origin for v and the terminus for u only the relavant two points are taken and the sum is divided
    * by two.
    *
    * \param[in] inputLevel Multigrid level of the input vector (not currently used)
    * \param[in] input The input vector to be "restricted" 
    * \param[in] output The output vector resulting from restriction, must be at level inputLevel+1
    */
   void restrictSimVector(const Vector<Real> & input,Vector<Real> & output)
   {
     const PinTVector<Real> & pint_input  = dynamic_cast<const PinTVector<Real>&>(input);
     PinTVector<Real>       & pint_output = dynamic_cast<PinTVector<Real>&>(output);

     int Np = (pint_output.numOwnedVectors()-2)/2 + 1;

     // handle the left end point on this processor
     ////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     {     
       int si = -1;                  // -1           
       int ns = si+1;                // 0
       int v_crs_index = 2*(ns-1)+1; // -1
       int v_index = 2*(2*ns-1)+1;   // -1

       ROL::Vector<Real> & out_v = *pint_output.getVectorPtr(v_crs_index);
       out_v.zero();
       out_v.axpy(1.0,*pint_input.getVectorPtr(v_index)); 
       out_v.axpy(1.0,*pint_input.getVectorPtr(v_index+2)); 
     
       int timeRank = pint_input.communicators().getTimeRank();
       {
         int crs_final_index = pint_output.numOwnedVectors()-2;
         int fine_final_index = pint_input.numOwnedVectors()-2;

         // copy final index into the buffer to be sent, this will be overwritten
         pint_output.getVectorPtr(crs_final_index)->set(*pint_input.getVectorPtr(fine_final_index-1));
       }

       pint_output.boundaryExchange();
       out_v.axpy(1.0,*pint_output.getRemoteBufferPtr(-1));

       if(timeRank==0)
         out_v.scale(1.0/2.0);
       else
         out_v.scale(1.0/3.0);
     }

     // handle interior points
     ////////////////////////////////////////////////////////////////////////////////////////////////////////
     for(int si=0;si<Np-1;si++) {

       // si = the step index, ns = the number of steps = si + 1
       int ns = si+1;

       // coarse indices
       int u_crs_index = 2*(ns-1);
       int v_crs_index = 2*(ns-1)+1;

       // fine indices
       int u_index = 2*(2*ns-1);   
       int v_index = 2*(2*ns-1)+1;

       ROL::Vector<Real> & out_u = *pint_output.getVectorPtr(u_crs_index);
       ROL::Vector<Real> & out_v = *pint_output.getVectorPtr(v_crs_index);

       out_u.zero();
       out_u.axpy(1.0,*pint_input.getVectorPtr(u_index-2)); 
       out_u.axpy(1.0,*pint_input.getVectorPtr(u_index+0)); 
       out_u.axpy(1.0,*pint_input.getVectorPtr(u_index+2)); 

       out_u.scale(1.0/3.0);

       out_v.zero();
       out_v.axpy(1.0,*pint_input.getVectorPtr(v_index-2)); 
       out_v.axpy(1.0,*pint_input.getVectorPtr(v_index+0)); 
       out_v.axpy(1.0,*pint_input.getVectorPtr(v_index+2)); 

       out_v.scale(1.0/3.0);
     }

     // handle the right end point on this processor
     ////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     {
       int si = Np-1;                   
       int ns = si+1;                   
       int u_crs_index = 2*(ns-1);      
       int u_index = 2*(2*ns-1);       

       ROL::Vector<Real> & out_u = *pint_output.getVectorPtr(u_crs_index);
       out_u.zero();
       out_u.axpy(1.0,*pint_input.getVectorPtr(u_index-2)); 
       out_u.axpy(1.0,*pint_input.getVectorPtr(u_index+0)); 
     
       pint_output.getRemoteBufferPtr(-1)->set(*pint_input.getVectorPtr(0));
       pint_output.boundaryExchangeSumInto();   // this adds the X_0^1 into the previous rank

       // scale end boundary point by 1/2, or the other end points by 1/3
       int timeRank = pint_input.communicators().getTimeRank();
       int timeSize = pint_input.communicators().getTimeSize();

       if(timeRank+1==timeSize)
         out_u.scale(1.0/2.0);
       else
         out_u.scale(1.0/3.0);
     }
   }

   void restrictSimVector(const Vector<Real> & input,ROL::Ptr<Vector<Real>> & output,int level)
   {
     // if you don't want to repartition then just restrict the vector
     if(not useRepart_) {
       restrictSimVector(input,*output);

       return;
     }

     // get the level communicators
     auto comm    = getLevelCommunicators(level);
     auto crsComm = getLevelCommunicators(level+1); // this will be null for processors that don't need it

     // now we need to repartion, allocate a local output vector
     auto output_local = allocateSimVector_local(input,level+1);

     // compute restricted vector
     restrictSimVector(input,*output_local);

     const PinTVector<Real> & pint_local  = dynamic_cast<const PinTVector<Real>&>(*output_local);

     // build a STL vector with the fine vectors
     int startIndex = (comm->getTimeRank() % 2==0) ? -1 : 0;
     std::vector<ROL::Ptr<ROL::Vector<Real>>> localVectors(pint_local.numOwnedSteps()-startIndex);
     bool is_null = false; // Unused
     for(int i=startIndex;i<pint_local.numOwnedSteps();i++) {
       localVectors[i-startIndex] = pint_local.getVectorPtr(i); 
       is_null =  ROL::is_nullPtr(localVectors[i-startIndex]);
//       if(ROL::is_nullPtr(localVectors[i-startIndex]))
//         is_null = true; 
     }

     // now send the coarse vectors this to the coarse data structure
     ROL::PinT::sendToCoarseDistribution_Vector(localVectors,*vectorComm_,*comm);

     // recv data
     if(not ROL::is_nullPtr(crsComm)) {
       PinTVector<Real>       & pint_remote = dynamic_cast<PinTVector<Real>&>(*output);

       // build a STL vector with the fine vectors
       int localSize = pint_local.numOwnedSteps();
       std::vector<ROL::Ptr<ROL::Vector<Real>>> remoteVectors(2*localSize+1);
       assert(pint_remote.numOwnedSteps()==2*localSize);

       for(int i=-1;i<pint_remote.numOwnedSteps();i++)
         remoteVectors[i+1] = pint_remote.getVectorPtr(i); 
   
       ROL::PinT::recvFromFineDistribution_Vector(remoteVectors,*vectorComm_,*comm,localSize+1,localSize);
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

   void restrictOptVector(const Vector<Real> & input,ROL::Ptr<Vector<Real>> & output,int level)
   {
     // if you don't want to repartition then just restrict the vector
     if(not useRepart_) {
       restrictOptVector(input,*output);

       return;
     }

     // get the level communicators
     auto comm    = getLevelCommunicators(level);
     auto crsComm = getLevelCommunicators(level+1); // this will be null for processors that don't need it

     // now we need to repartion, allocate a local output vector
     auto output_local = allocateOptVector_local(input,level+1);

     // compute restricted vector
     restrictOptVector(input,*output_local);

     const PinTVector<Real> & pint_local  = dynamic_cast<const PinTVector<Real>&>(*output_local);

     // build a STL vector with the fine vectors
     std::vector<ROL::Ptr<ROL::Vector<Real>>> localVectors(pint_local.numOwnedSteps());
     for(int i=0;i<pint_local.numOwnedSteps();i++)
       localVectors[i] = pint_local.getVectorPtr(i); 

     // now send the coarse vectors this to the coarse data structure
     ROL::PinT::sendToCoarseDistribution_Vector(localVectors,*vectorComm_,*comm);

     // recv data
     if(not ROL::is_nullPtr(crsComm)) {
       PinTVector<Real>       & pint_remote = dynamic_cast<PinTVector<Real>&>(*output);

       // build a STL vector with the fine vectors
       int localSize = pint_local.numOwnedSteps();
       std::vector<ROL::Ptr<ROL::Vector<Real>>> remoteVectors(2*localSize);
       for(int i=0;i<pint_remote.numOwnedSteps();i++)
         remoteVectors[i] = pint_remote.getVectorPtr(i); 
   
       ROL::PinT::recvFromFineDistribution_Vector(remoteVectors,*vectorComm_,*comm,localSize,localSize);
     }
   }

   /**
    * \brief Prolong a simulation space vector
    *
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

     int Np = (pint_output.numOwnedVectors()-2)/2 + 1;

     // handle left exterior points
     ////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     pint_input.getRemoteBufferPtr(-1)->setScalar(0.0);
     pint_input.boundaryExchange();
     
     {
       // injection of the virtual point
       pint_output.getVectorPtr(-1)->set(*pint_input.getVectorPtr(-1));

       // averaging of the state points
       pint_output.getVectorPtr( 0)->set(*pint_input.getRemoteBufferPtr(-1));
       pint_output.getVectorPtr( 0)->axpy(1.0,*pint_input.getVectorPtr(0));
       pint_output.getVectorPtr( 0)->scale(0.5);

       // averaging of the virtual points
       pint_output.getVectorPtr( 1)->set(*pint_input.getVectorPtr(-1));
       pint_output.getVectorPtr( 1)->axpy(1.0,*pint_input.getVectorPtr(1));
       pint_output.getVectorPtr( 1)->scale(0.5);
     }

     // handle interior points
     ////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     for(int si=1;si<Np-2;si++) {

       // fine indices
       int u_index = 2*si;   
       int v_index = 2*si+1;   

       if(si % 2 == 1) {
         int coarse_si = (si-1)/2;
         int u_crs_index = 2*coarse_si;
         int v_crs_index = 2*coarse_si+1;

         pint_output.getVectorPtr(u_index)->set(*pint_input.getVectorPtr(u_crs_index));
         pint_output.getVectorPtr(v_index)->set(*pint_input.getVectorPtr(v_crs_index));
       }
       else {
         int coarse_si_0 = si/2-1;
         int coarse_si_1 = si/2;
         int u_crs_index = -1; 
         int v_crs_index = -1;
 
         u_crs_index = 2*coarse_si_0;
         v_crs_index = 2*coarse_si_0+1;

         pint_output.getVectorPtr(u_index)->set(*pint_input.getVectorPtr(u_crs_index));
         pint_output.getVectorPtr(v_index)->set(*pint_input.getVectorPtr(v_crs_index));

         u_crs_index = 2*coarse_si_1;
         v_crs_index = 2*coarse_si_1+1;

         pint_output.getVectorPtr(u_index)->axpy(1.0,*pint_input.getVectorPtr(u_crs_index));
         pint_output.getVectorPtr(v_index)->axpy(1.0,*pint_input.getVectorPtr(v_crs_index));

         pint_output.getVectorPtr(u_index)->scale(0.5);
         pint_output.getVectorPtr(v_index)->scale(0.5);
       }
     }

     // handle left exterior points
     ////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     {
       int fine_last = pint_output.numOwnedVectors();
       int crs_last = pint_input.numOwnedVectors();

       pint_output.getRemoteBufferPtr(-1)->set(*pint_input.getVectorPtr(-1)); // copy the virtual variable
       pint_output.getVectorPtr(fine_last-2)->zero();  // make space
       pint_output.boundaryExchangeSumInto();

       // averaging of the virtual points
       pint_output.getVectorPtr(fine_last-3)->set(*pint_output.getVectorPtr(fine_last-2));
       pint_output.getVectorPtr(fine_last-3)->axpy(1.0,*pint_input.getVectorPtr(crs_last-3));
       pint_output.getVectorPtr(fine_last-3)->scale(0.5);

       // averaging of the virtual points
       pint_output.getVectorPtr(fine_last-4)->set(*pint_input.getVectorPtr(crs_last-2));
       pint_output.getVectorPtr(fine_last-4)->axpy(1.0,*pint_input.getVectorPtr(crs_last-4));
       pint_output.getVectorPtr(fine_last-4)->scale(0.5);

       // injection of the state point
       pint_output.getVectorPtr(fine_last-2)->set(*pint_input.getVectorPtr(crs_last-2));

     }
   }

   void prolongSimVector(const ROL::Ptr<const Vector<Real>> & input,Vector<Real> & output,int level)
   {
     // if you don't want to repartition then just restrict the vector
     if(not useRepart_) {
       prolongSimVector(*input,output);

       return;
     }

     // get the level communicators
     auto comm    = getLevelCommunicators(level);

     // do this only on the coarse grid
     if(not ROL::is_nullPtr(input)) {
       // now we need to repartion, allocate a local output vector
       auto output_local = allocateSimVector_local_fine(*input,level);

       prolongSimVector(*input,*output_local);

       const PinTVector<Real> & pint_local  = dynamic_cast<const PinTVector<Real>&>(*output_local);

       int localSize = pint_local.numOwnedSteps()/2;

       // build a STL vector with the fine vectors
       std::vector<ROL::Ptr<ROL::Vector<Real>>> localVectors_a(localSize+1,ROL::nullPtr);
       std::vector<ROL::Ptr<ROL::Vector<Real>>> localVectors_b(localSize+1,ROL::nullPtr);
       for(int i=-1;i<localSize;i++)
         localVectors_a[i+1] = pint_local.getVectorPtr(i); 

       for(int i=localSize-1;i<pint_local.numOwnedSteps();i++)
         localVectors_b[i-localSize+1] = pint_local.getVectorPtr(i); 

       // send data
       ROL::PinT::sendToFineDistribution_Vector(localVectors_a,localVectors_b,*vectorComm_,*comm);
     }

     PinTVector<Real> & pint_remote = dynamic_cast<PinTVector<Real>&>(output);

     // build a STL vector with the fine vectors
     std::vector<ROL::Ptr<ROL::Vector<Real>>> remoteVectors(pint_remote.numOwnedSteps()+1);
     for(int i=-1;i<pint_remote.numOwnedSteps();i++)
       remoteVectors[i+1] = pint_remote.getVectorPtr(i); 

     ROL::PinT::recvFromCoarseDistribution_Vector(remoteVectors,*vectorComm_,*comm);
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

   void prolongOptVector(const ROL::Ptr<const Vector<Real>> & input,Vector<Real> & output,int level)
   {
     // if you don't want to repartition then just restrict the vector
     if(not useRepart_) {
       prolongOptVector(*input,output);

       return;
     }

     // get the level communicators
     auto comm    = getLevelCommunicators(level);

     std::string levelRankStr;
     {
       std::stringstream ss;
       ss << "-" << level << "-";
       if(ROL::is_nullPtr(comm))
         ss << "null";
       else
         ss << comm->getTimeRank();
       levelRankStr = ss.str();
     }


     // do this only on the coarse grid
     if(not ROL::is_nullPtr(input)) {
       // now we need to repartion, allocate a local output vector
       auto output_local = allocateOptVector_fine(*input,level);

       prolongOptVector(*input,*output_local);

       const PinTVector<Real> & pint_local  = dynamic_cast<const PinTVector<Real>&>(*output_local);

       int localSize = pint_local.numOwnedSteps()/2;

       // build a STL vector with the fine vectors
       std::vector<ROL::Ptr<ROL::Vector<Real>>> localVectors_a(localSize,ROL::nullPtr);
       std::vector<ROL::Ptr<ROL::Vector<Real>>> localVectors_b(localSize,ROL::nullPtr);
       for(int i=0;i<localSize;i++)
         localVectors_a[i] = pint_local.getVectorPtr(i); 

       for(int i=localSize;i<pint_local.numOwnedSteps();i++)
         localVectors_b[i-localSize] = pint_local.getVectorPtr(i); 

       // send data
       ROL::PinT::sendToFineDistribution_Vector(localVectors_a,localVectors_b,*vectorComm_,*comm);
     }

     PinTVector<Real> & pint_remote = dynamic_cast<PinTVector<Real>&>(output);

     // build a STL vector with the fine vectors
     std::vector<ROL::Ptr<ROL::Vector<Real>>> remoteVectors(pint_remote.numOwnedSteps());
     for(int i=0;i<pint_remote.numOwnedSteps();i++)
       remoteVectors[i] = pint_remote.getVectorPtr(i); 

     ROL::PinT::recvFromCoarseDistribution_Vector(remoteVectors,*vectorComm_,*comm);
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
     applyJacobian_2_leveled(*output_v,*input_z,u,z,tol,level);   // BAD ???

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

     bool printResidualInfo = true;

     auto timer = Teuchos::TimeMonitor::getStackedTimer();

     const PinTVector<Real>       & pint_u = dynamic_cast<const PinTVector<Real>&>(u);

     // std::cout << "SOLVING ON LEVEL " << level << std::endl;

     std::string levelStr = "";
     std::string levelRankStr = "";
     {
       std::stringstream ss;
       ss << "-" << level;
       levelStr = ss.str();
     }

     if(useRepart_) {
       // get the level communicators
       auto comm = getLevelCommunicators(level);

       std::stringstream ss;
       ss << levelStr << "-";
       if(ROL::is_nullPtr(comm))
         ss << "null";
       else
         ss << comm->getTimeRank();
       levelRankStr = ss.str();
     }

     timer->start("applyMGAugmentedKKT"+levelStr);

     double b_norm = b.norm();

     std::string levelIndent = "";
     if(printResidualInfo) {
       int timeRank = pint_u.communicators().getTimeRank();
       for(int i=0;i<level+1;i++)
         levelIndent += "   ";

       if(timeRank==0 && level==0) {
         std::cout << "Start new solve: " << std::endl;
       }
       else if(timeRank==0) {
         std::cout << levelIndent << "Next level " << level << std::endl;
       }
     }

     // base case: solve the KKT system directly
     if(level+1==maxLevels_) {
       bool approxSmoother = false;

       auto dx = x.clone();
       auto residual = b.clone();
       residual->set(b);

       double relax = omegaCoarse_;

       for(int i=0;i<numCoarseSweeps_;i++) {
         // compute the residual
         applyAugmentedKKT(*residual,x,u,z,tol,level);
         residual->scale(-1.0);
         residual->axpy(1.0,b);

         applyAugmentedInverseKKT(*dx,*residual,u,z,tol,approxSmoother,level); 
         x.axpy(relax,*dx);
       }

       if(printResidualInfo) {
         applyAugmentedKKT(*residual,x,u,z,tol,level);
         residual->scale(-1.0);
         residual->axpy(1.0,b);
         int timeRank = pint_u.communicators().getTimeRank();
         double res_norm = residual->norm();
         if(timeRank == 0)
           std::cout << levelIndent << "(level " << level << ") coarse reduction = " << res_norm/b_norm << std::endl;
       }

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

     if(printResidualInfo) {
       int timeRank = pint_u.communicators().getTimeRank();
       double res_norm = residual->norm();
       if(timeRank == 0)
         std::cout << levelIndent << "(level " << level << ") pre reduction = " << res_norm/b_norm << std::endl;
     }

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

       restrictSimVector(*residual_u,crs_residual_u,level);
       restrictOptVector(*residual_z,crs_residual_z,level);
       restrictSimVector(*residual_v,crs_residual_v,level);

       restrictSimVector(u,crs_u,level);               // restrict the state to the coarse level
       restrictOptVector(z,crs_z,level);               // restrict the control to the coarse level

       // if(not ROL::is_nullPtr(getLevelCommunicators(level+1))) 
       {
         auto crs_correction = makePtr<PartitionedVector>({crs_correction_u,crs_correction_z,crs_correction_v});
         auto crs_residual   = makePtr<PartitionedVector>({  crs_residual_u,  crs_residual_z,  crs_residual_v});

         applyMultigridAugmentedKKT(*crs_correction,*crs_residual,*crs_u,*crs_z,tol,level+1);
       }

       prolongSimVector(crs_correction_u,*dx_u,level);
       prolongOptVector(crs_correction_z,*dx_z,level);
       prolongSimVector(crs_correction_v,*dx_v,level);

       x.axpy(1.0,*dx);
     }
     timer->stop("applyMGAugmentedKKT-coarse");

     if(printResidualInfo) {
       applyAugmentedKKT(*residual,x,u,z,tol,level);
       residual->scale(-1.0);
       residual->axpy(1.0,b);

       int timeRank = pint_u.communicators().getTimeRank();
       double res_norm = residual->norm();
       if(timeRank == 0)
         std::cout << levelIndent << "(level " << level << ") coarse solve reduction = " << res_norm/b_norm << std::endl;
     }

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

     if(printResidualInfo) {
       applyAugmentedKKT(*residual,x,u,z,tol,level);
       residual->scale(-1.0);
       residual->axpy(1.0,b);

       int timeRank = pint_u.communicators().getTimeRank();
       double res_norm = residual->norm();
       if(timeRank == 0)
         std::cout << levelIndent << "(level " << level << ") post reduction = " << res_norm/b_norm << std::endl;
     }

     timer->stop("applyMGAugmentedKKT"+levelStr);
   }

}; // ROL::PinTConstraint


} // namespace ROL 

#endif // ROL_PINTCONSTRAINT_HPP

