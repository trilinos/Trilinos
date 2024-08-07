// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PINTCONSTRAINT_HPP
#define ROL_PINTCONSTRAINT_HPP

#include "ROL_TimeStamp.hpp"
#include "ROL_PinTVector.hpp"
#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_SerialConstraint.hpp"
#include "ROL_SerialConstraint.hpp"
#include "ROL_PinTCommunicationUtilities.hpp"
#include "ROL_PinTHierarchy.hpp"

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
  if((steps - 1) % 2!=0) // steps must be 2^n+1
    throw std::logic_error("Wrong number of steps, must be 2^n+1");
  return makePtr<PinTVector<Real>>(communicators,vectorComm,localVector,steps,-1,1,replicate); }

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
{ const int replicate = 1; // no virtual variables
  if((steps - 1) % 2!=0) // steps must be 2^n+1
    throw std::logic_error("Wrong number of steps, must be 2^n+1");
  return makePtr<PinTVector<Real>>(communicators,vectorComm,localVector,steps,-1,1 /* buffer size of 1*/,replicate); }

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
  Ptr<Vector<Real>>      initialCond_;

  Ptr<const std::vector<TimeStamp<Real>>>         userTimeStamps_;    // these are the untouched ones the user passes in

  ROL::PinTHierarchy<Real> hierarchy_;

  // preconditioner settings
  bool applyMultigrid_;         // default to block jacobi
  int maxLevels_;               // must turn on multigrid 

  int numSweeps_;
  int numCoarseSweeps_;
  Real omega_;
  Real omegaCoarse_;

  int numCGIter_; // number of CG interations for the smoother

  Real globalScale_;

  Real controlRegParam_;

  // clone vector storage
  struct WathenInverseStorage {
    Ptr<Vector<Real>> temp_u;
    Ptr<Vector<Real>> temp_z;
    Ptr<Vector<Real>> temp_v;
    Ptr<Vector<Real>> temp_schur;
  };

  // clone vector storage
  struct MGAugmentedKKTStorage {
    Ptr<Vector<Real>> dx;
    Ptr<Vector<Real>> residual;
  };

  // CG Reduced Hessian Solver
  struct CGReducedHessStorage {
    Ptr<Vector<Real>> temp_r;
    Ptr<Vector<Real>> temp_Ap;
    Ptr<Vector<Real>> temp_p;
  };

  std::map<int,WathenInverseStorage> inverseKKTStorage_; 

  std::map<int,MGAugmentedKKTStorage> mgAugmentedKKTStorage_;

  std::map<int,CGReducedHessStorage> cgReducedHessStorage_;

  std::map<int,Ptr<Vector<Real>>> inverseAdjointStorage_;

  bool recordResidualReductions_;
  std::map<int,std::vector<double>> preSmoothResidualReduction_;
  std::map<int,std::vector<double>> postSmoothResidualReduction_;
  std::vector<double> coarseResidualReduction_;

  /** Get the state vector required by the serial constraint object
    *
    * \param[in] src Vector to build the state vector from
    * \param[in] s  Time step index
    * \param[in] useInitialCond At the boundary, use the initial condition if true, otherwise
    *                           use z zero vector.
    */
  Ptr<PartitionedVector<Real>> getStateVector(const PinTVector<Real> & src,int s,bool useInitialCond)
  {
    int timeRank = src.communicators().getTimeRank();
    int oldTimeIndex = 2*s-1;

    // the state PinTVector is assuemd to have two vectors per time step,
    // one for the "primary" variable (first, even index), and the "virtual" variable (second, odd index)

    std::vector<Ptr<Vector<Real>>> vecs;
    if(oldTimeIndex==-1 && timeRank==0) {
      auto val = initialCond_->clone();
      if(useInitialCond) 
        val->set(*initialCond_); 
      else 
        val->scale(0.0);
      vecs.push_back(val);
    }
    else if(oldTimeIndex==-1) 
      vecs.push_back(src.getRemoteBufferPtr(1));  // this index is 1 to get the the virtual variable from the left processor
    else if(oldTimeIndex>=0) 
      vecs.push_back(src.getVectorPtr(oldTimeIndex));   
    else {
      std::stringstream ss; 
      ss << "ROL::PinTConstraint::getStateVector: unxpected step index\n"
         << "   (dev info: s=" << s 
         << ", timeRank=" << timeRank 
         << ", oldTimeIndex=" << oldTimeIndex 
         << ")" << std::endl;
  
      throw std::logic_error(ss.str());
    }

    vecs.push_back(src.getVectorPtr(oldTimeIndex+1));   // true value lives on the evens

    return makePtr<PartitionedVector<Real>>(vecs);
  }

  // Get the control vector required by the serial constraint object
  Ptr<Vector<Real>> getControlVector(const PinTVector<Real> & src,int s=0)
  {
    // notice this ignores the stencil for now

    std::vector<Ptr<Vector<Real>>> vecs;
    vecs.push_back(src.getVectorPtr(0)->clone());
      // FIXME-A: this control index is never used. This stems from not wanting to set
      //          an intial condition for the test function. But instead use the 
      //          "setSkipInitialCondition" option. This all needs to be rethought. 
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
    , numCGIter_(2)
    , globalScale_(0.99e0)
    , recordResidualReductions_(false)
  { }

  /**
   * \brief Constructor
   *
   * Build a parallel-in-time constraint with a specified step constraint. This specifies
   * any communication you might need between processors and steps using "stencils".
   * Multigrid in time preconditioning is disabled by default.
   * 
   * \param[in] stepConstraint Constraint for a single step.
   * \param[in] initialCond Initial condition, this will be replicated with a deep copy
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
    , numCGIter_(2)
    , globalScale_(0.99e0)
    , controlRegParam_(1.0)
    , recordResidualReductions_(false)
  { 
    initialize(stepConstraint,initialCond,timeStamps);
  }

  /**
   * Set sweeps for jacobi and multigrid
   */
  void setSweeps(int s)
  { numSweeps_ = s; }

  void setCGIterations(int iter)
  { numCGIter_ = iter; }

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
   * Set the scaling for the control regularization parameter,
   * this defaults to 1.0. Only relevant for the KKT system and
   * multigrid solvers.
   */
  void setControlRegParam(Real o)
  { controlRegParam_ = o; }

  /**
   *  For analysis purposes, record all the residual reductions.
   */
  void setRecordResidualReductions(bool record)
  { recordResidualReductions_ = record; } 

  const std::map<int,std::vector<double>> & getPreSmoothResidualReductions() const
  { return preSmoothResidualReduction_; }

  const std::map<int,std::vector<double>> & getPostSmoothResidualReduction() const
  { return postSmoothResidualReduction_; }

  const std::vector<double> & getCoarseResidualReduction() const
  { return coarseResidualReduction_; }

  void clearResidualReduction() 
  { 
    preSmoothResidualReduction_.clear();
    postSmoothResidualReduction_.clear();
    coarseResidualReduction_.clear();
  }

  /**
   * Turn on multigrid preconditioning in time with a specified number of levels.
   *
   * \param[in] maxLevels Largest number of levels
   */
  void applyMultigrid(int maxLevels,
                      const ROL::Ptr<const ROL::PinTCommunicators> & pintComm,
                      const ROL::Ptr<const ROL::PinTVectorCommunication<Real>> & vectorComm,
                      bool rebalance=false)
  {
    applyMultigrid_ = true;
    maxLevels_ = maxLevels;

    hierarchy_.setMaxLevels(maxLevels_);
    hierarchy_.buildLevels(pintComm,vectorComm,rebalance);
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
   */
  Ptr<std::vector<TimeStamp<Real>>> getTimeStampsByLevel(int level) const
  {
    return hierarchy_.getTimeStampsByLevel(level);
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
    timeStamps = hierarchy_.getTimeStampsByLevel(level);

    // build a single point time stamp

    VecPtr localTimeStamps = makePtr<std::vector<TimeStamp<Real>>>(2);
    if(step==0) {
      localTimeStamps->at(0) = timeStamps->at(step);

      // shift it backwards in time one step... this is a hack
      Real dt = localTimeStamps->at(0).t[1] - localTimeStamps->at(0).t[0];
      localTimeStamps->at(0).t[0] -= dt; 
      localTimeStamps->at(0).t[1] -= dt;
 
      localTimeStamps->at(0).k--;
    }
    else {
      localTimeStamps->at(0) = timeStamps->at(step-1);
    }

    localTimeStamps->at(1) = timeStamps->at(step);

    if(ROL::is_nullPtr(timeDomainConstraint_)) {

      timeDomainConstraint_ = ROL::makePtr<SerialConstraint<Real>>(stepConstraint_,*ui,localTimeStamps);
    }
    else {
      timeDomainConstraint_->setTimeStampsPtr(localTimeStamps);
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
    initialCond_ = initialCond->clone();
    initialCond_->set(*initialCond);
    userTimeStamps_ = timeStamps;

    hierarchy_.initialize(timeStamps);

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

    pint_z.boundaryExchangeLeftToRight();
    pint_u.boundaryExchangeLeftToRight(PinTVector<Real>::RECV_ONLY); // this is going to block

    size_t numSteps = hierarchy_.getTimeStampsByLevel(0)->size();
    for(size_t s=0;s<numSteps;s++) { // num time steps == num time stamps-1 

      auto part_c = getStateVector(pint_c,s,false);   
      auto part_u = getStateVector(pint_u,s,true);   
      auto part_z = getControlVector(pint_z,s); 

      // compute the constraint for this subdomain
      //    u_s = F(v_{s-1},z_s)
      auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s),level,s);

      constraint->solve(*part_c,*part_u,*part_z,tol);

      // satisfy the time continuity constraint, note that this just using the identity matrix
      pint_u.getVectorPtr(2*s+1)->set(*pint_u.getVectorPtr(2*s));
    }

    pint_u.boundaryExchangeLeftToRight(PinTVector<Real>::SEND_ONLY);

    // call the original implementation
    value(c,u,z,tol);
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
    pint_u.boundaryExchangeLeftToRight();
    pint_z.boundaryExchangeLeftToRight();

    size_t numSteps = hierarchy_.getTimeStampsByLevel(0)->size();
    for(size_t s=0;s<numSteps;s++) { // num time steps == num time stamps-1 

      auto part_c = getStateVector(pint_c,s,false);   
      auto part_u = getStateVector(pint_u,s,true);   
      auto part_z = getControlVector(pint_z,s); 

      // compute the constraint for this subdomain
      //    u_s = F(v_{s-1},z_s)
      auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s),level,s);

      constraint->value(*part_c,*part_u,*part_z,tol);

      // build in the time continuity constraint, note that this is just using the identity matrix
      //    v_s = u_s
      pint_c.getVectorPtr(2*s+1)->set(*pint_u.getVectorPtr(2*s+1));     // this is the virtual value
      pint_c.getVectorPtr(2*s+1)->axpy(-1.0,*pint_u.getVectorPtr(2*s)); // this is the u value
      pint_c.getVectorPtr(2*s+1)->scale(globalScale_);
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
    pint_v.boundaryExchangeLeftToRight();
    pint_u.boundaryExchangeLeftToRight();
    pint_z.boundaryExchangeLeftToRight();

    size_t numSteps = hierarchy_.getTimeStampsByLevel(level)->size();
    for(size_t s=0;s<numSteps;s++) { // num time steps == num time stamps-1 

      auto part_jv = getStateVector(pint_jv,s,false); 
      auto part_v  = getStateVector(pint_v,s,false);    
      auto part_u  = getStateVector(pint_u,s,true);   
      auto part_z  = getControlVector(pint_z,s); 

      // compute the constraint for this subdomain
      //    u_s = F(v_{s-1},z_s)
      auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s),level,s);

      constraint->applyJacobian_1(*part_jv,*part_v,*part_u,*part_z,tol);

      // build in the time continuity constraint, note that this is just using the identity matrix
      //    v_s = u_s
      pint_jv.getVectorPtr(2*s+1)->set(*pint_v.getVectorPtr(2*s+1));     // this is the virtual value
      pint_jv.getVectorPtr(2*s+1)->axpy(-1.0,*pint_v.getVectorPtr(2*s)); // this is the u value
      pint_jv.getVectorPtr(2*s+1)->scale(globalScale_);
    }
  }

  void applyJacobian_1_leveled_approx( V& jv, const V& v, const V& u, 
                       const V& z, Real& tol,int level) {
    PinTVector<Real>       & pint_jv = dynamic_cast<PinTVector<Real>&>(jv);
    const PinTVector<Real> & pint_v  = dynamic_cast<const PinTVector<Real>&>(v);
    const PinTVector<Real> & pint_u  = dynamic_cast<const PinTVector<Real>&>(u);
    const PinTVector<Real> & pint_z  = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here
      
    int timeRank = pint_u.communicators().getTimeRank();
    int timeSize = pint_u.communicators().getTimeSize();
    bool lastRank = (timeRank+1==timeSize); // do something special on the last rank
 
    assert(pint_jv.numOwnedSteps()==pint_u.numOwnedSteps());
    assert(pint_jv.numOwnedSteps()==pint_v.numOwnedSteps());

    // communicate neighbors, these are blocking calls
    pint_v.boundaryExchangeLeftToRight();
    // pint_u.boundaryExchangeLeftToRight();
    // pint_z.boundaryExchangeLeftToRight();

    size_t numSteps = hierarchy_.getTimeStampsByLevel(level)->size();
    for(size_t s=0;s<numSteps;s++) { // num time steps == num time stamps-1 

      auto part_jv = getStateVector(pint_jv,s,false); 
      auto part_v  = getStateVector(pint_v,s,false);    
      auto part_u  = getStateVector(pint_u,s,true);   
      auto part_z  = getControlVector(pint_z,s); 
  
      // compute the constraint for this subdomain
      //    u_s = F(v_{s-1},z_s)
      auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s),level,s);
  
      constraint->applyJacobian_1(*part_jv,*part_v,*part_u,*part_z,tol);
  
      // build in the time continuity constraint, note that this is just using the identity matrix
      //    v_s = u_s
      pint_jv.getVectorPtr(2*s+1)->set(*pint_v.getVectorPtr(2*s+1));       // this is the virtual value
      if(s+1<numSteps || lastRank)
        pint_jv.getVectorPtr(2*s+1)->axpy(-1.0,*pint_v.getVectorPtr(2*s)); // this is the u value
      pint_jv.getVectorPtr(2*s+1)->scale(globalScale_);
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
    pint_v.boundaryExchangeLeftToRight();
    pint_u.boundaryExchangeLeftToRight();
    pint_z.boundaryExchangeLeftToRight();

    size_t numSteps = hierarchy_.getTimeStampsByLevel(level)->size();
    for(size_t s=0;s<numSteps;s++) { // num time steps == num time stamps-1 

      auto part_jv = getStateVector(pint_jv,s,false); 
      auto part_v  = getControlVector(pint_v,s);    
      auto part_u  = getStateVector(pint_u,s,true);   
      auto part_z  = getControlVector(pint_z,s); 
  
      // compute the constraint for this subdomain
      //    u_s = F(v_{s-1},z_s)
      auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s),level,s);
  
      constraint->applyJacobian_2(*part_jv,*part_v,*part_u,*part_z,tol);
  
      // build in the time continuity constraint, note that this is just using the identity matrix
      //    v_s = u_s
      pint_jv.getVectorPtr(2*s+1)->zero();
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
     
     // ijv.zero();

     PinTVector<Real>       & pint_ijv = dynamic_cast<PinTVector<Real>&>(ijv);
     const PinTVector<Real> & pint_v   = dynamic_cast<const PinTVector<Real>&>(v);
     const PinTVector<Real> & pint_u   = dynamic_cast<const PinTVector<Real>&>(u);
     const PinTVector<Real> & pint_z   = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here

     assert(pint_ijv.numOwnedSteps()==pint_u.numOwnedSteps());
     assert(pint_ijv.numOwnedSteps()==pint_v.numOwnedSteps());

     pint_z.boundaryExchangeLeftToRight();
     pint_v.boundaryExchangeLeftToRight();
     pint_u.boundaryExchangeLeftToRight();
     pint_ijv.boundaryExchangeLeftToRight(PinTVector<Real>::RECV_ONLY); // this is going to block

     size_t numSteps = hierarchy_.getTimeStampsByLevel(level)->size();
     for(size_t s=0;s<numSteps;s++) {

       auto part_ijv = getStateVector(pint_ijv,s,false);   
       auto part_v   = getStateVector(pint_v,s,false);   
       auto part_u   = getStateVector(pint_u,s,true);   
       auto part_z   = getControlVector(pint_z,s); 

       // compute the constraint for this subdomain
       auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s),level,s);

       constraint->applyInverseJacobian_1(*part_ijv,*part_v,*part_u,*part_z,tol);

       // satisfy the time continuity constraint, note that this just using the identity matrix
       // g * v - g * u = r
       // v = u + r/g;

       pint_ijv.getVectorPtr(2*s+1)->set(*pint_ijv.getVectorPtr(2*s));    
       pint_ijv.getVectorPtr(2*s+1)->axpy(1.0/globalScale_,*pint_v.getVectorPtr(2*s+1));        
     }

     pint_ijv.boundaryExchangeLeftToRight(PinTVector<Real>::SEND_ONLY); // this releases the block
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
     pint_v.boundaryExchangeLeftToRight();
     pint_u.boundaryExchangeLeftToRight();
     pint_z.boundaryExchangeLeftToRight();

     Ptr<Vector<Real>> virt;
 
     int numSteps = Teuchos::as<int>(hierarchy_.getTimeStampsByLevel(level)->size());
     for(int s=numSteps-1;s>=0;s--) { 
       auto part_ajv = getStateVector(pint_ajv,s,false);
       auto part_v   = getStateVector(pint_v,s,false);   
       auto part_u   = getStateVector(pint_u,s,true);   
       auto part_z   = getControlVector(pint_z,s); 

       // handle the constraint adjoint: Required for correctly applying step adjoint
       pint_ajv.getVectorPtr(2*s+1)->axpy(globalScale_,*pint_v.getVectorPtr(2*s+1));

       // compute the constraint for this subdomain
       auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s),level,s);

       constraint->applyAdjointJacobian_1(*part_ajv,*part_v,*part_u,*part_z,*part_v,tol);

       // this is the remainder of the constraint application
       pint_ajv.getVectorPtr(2*s)->axpy(-globalScale_,*pint_v.getVectorPtr(2*s+1));

       virt = part_ajv->get(0);
     }

     //
     //  [ I       ]'   [ I A'    ]
     //  [ A  B    ]  = [   B' -I ]
     //  [   -I  I ]    [       I ]
     //
     
     assert(not is_null(virt));
     
     pint_ajv.boundaryExchangeRightToLeft({virt});

     // this sums from the remote buffer into the local buffer which is part of the vector
     int timeRank = pint_ajv.communicators().getTimeRank();
     int timeSize = pint_ajv.communicators().getTimeSize();
     if(timeRank+1<timeSize) {
       pint_ajv.getVectorPtr(2*(numSteps-1)+1)->axpy(1.0,*pint_ajv.getRemoteBufferPtr(0));
     }
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
     
     int timeRank = pint_u.communicators().getTimeRank();
     int timeSize = pint_u.communicators().getTimeSize();
     bool lastRank = (timeRank+1==timeSize); // do something special on the last rank

     assert(pint_v.numOwnedSteps()==pint_u.numOwnedSteps());
     assert(pint_ajv.numOwnedSteps()==pint_u.numOwnedSteps());

     // we need to make sure this has all zeros to begin with (this includes boundary exchange components)
     pint_ajv.zeroAll();

     // communicate neighbors, these are block calls
     pint_v.boundaryExchangeLeftToRight();
     // pint_u.boundaryExchangeLeftToRight();
     // pint_z.boundaryExchangeLeftToRight();

     std::vector<Ptr<Vector<Real>>> sendBuffer(1);
 
     int numSteps = Teuchos::as<int>(hierarchy_.getTimeStampsByLevel(level)->size());
     for(int s=numSteps-1;s>=0;s--) { 
       auto part_ajv = getStateVector(pint_ajv,s,false);
       auto part_v   = getStateVector(pint_v,s,false);   
       auto part_u   = getStateVector(pint_u,s,true);   
       auto part_z   = getControlVector(pint_z,s); 

       // handle the constraint adjoint: Required for correctly applying step adjoint
       pint_ajv.getVectorPtr(2*s+1)->axpy(globalScale_,*pint_v.getVectorPtr(2*s+1));

       // compute the constraint for this subdomain
       auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s),level,s);

       constraint->applyAdjointJacobian_1(*part_ajv,*part_v,*part_u,*part_z,*part_v,tol);

       // this is the remainder of the constraint application
       if(s+1<numSteps || lastRank) 
         pint_ajv.getVectorPtr(2*s)->axpy(-globalScale_,*pint_v.getVectorPtr(2*s+1));

       sendBuffer[0] = part_ajv->get(0);
     }

     pint_ajv.boundaryExchangeRightToLeft(sendBuffer);

     // this sums from the remote buffer into the local buffer which is part of the vector
     if(not lastRank) {
       pint_ajv.getVectorPtr(2*(numSteps-1)+1)->axpy(1.0,*pint_ajv.getRemoteBufferPtr(0));
     }
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

     // communicate neighbors, these are block calls
     pint_v.boundaryExchangeLeftToRight();
     pint_u.boundaryExchangeLeftToRight();
     pint_z.boundaryExchangeLeftToRight();

     int numSteps = Teuchos::as<int>(hierarchy_.getTimeStampsByLevel(level)->size());
     for(int s=numSteps-1;s>=0;s--) { 

       auto part_ajv = getControlVector(pint_ajv,s);
       auto part_v   = getStateVector(pint_v,s,false);   
       auto part_u   = getStateVector(pint_u,s,true);   
       auto part_z   = getControlVector(pint_z,s); 

       // compute the constraint for this subdomain
       auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s),level,s);

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

     iajv.zero();

     // applyInverseAdjointJacobian_1 is weird because it serializes in time. But we want to get it
     // working so that we can test, but it won't be a performant way to do this. 
     // We could do a parallel in time solve here but thats not really the focus.

     // this is an inefficient hack (see line below where pint_v is modified!!!!)
     ///////////////////////////////////
     Ptr<Vector<Real>> v_copy;

     // handle the lazy construction by level
     {
       auto store = inverseAdjointStorage_.find(level);
       if(store==inverseAdjointStorage_.end()) {
         v_copy = v.clone();
         inverseAdjointStorage_[level] = v_copy;
       }
       else {
         v_copy = store->second;
       }
     }

     v_copy->set(v);
     ///////////////////////////////////

     PinTVector<Real>       & pint_iajv = dynamic_cast<PinTVector<Real>&>(iajv);
     PinTVector<Real> & pint_v          = dynamic_cast<PinTVector<Real>&>(*v_copy);
     const PinTVector<Real> & pint_u    = dynamic_cast<const PinTVector<Real>&>(u);
     const PinTVector<Real> & pint_z    = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here

     assert(pint_iajv.numOwnedSteps()==pint_u.numOwnedSteps());
     assert(pint_iajv.numOwnedSteps()==pint_v.numOwnedSteps());

     std::vector<Ptr<Vector<Real>>> sendBuffer(1);

     pint_z.boundaryExchangeLeftToRight();
     pint_v.boundaryExchangeLeftToRight();
     pint_u.boundaryExchangeLeftToRight();
     pint_iajv.boundaryExchangeRightToLeft(sendBuffer,PinTVector<Real>::RECV_ONLY); // this is going to block

     int timeRank = pint_iajv.communicators().getTimeRank();
     int timeSize = pint_iajv.communicators().getTimeSize();

     int numSteps = Teuchos::as<int>(hierarchy_.getTimeStampsByLevel(level)->size());
     if(timeRank+1<timeSize) {
       // the serial constraint inverse adjoint jacobian computes this term fully
       // if the pint_v exachange left to right is called
       pint_iajv.getVectorPtr(2*(numSteps-1)+1)->set(*pint_iajv.getRemoteBufferPtr(0));
     }
     else {
       // this is the last rank case
       pint_iajv.getVectorPtr(2*(numSteps-1)+1)->set(*pint_v.getVectorPtr(2*(numSteps-1)+1));
     }

     for(int s=numSteps-1;s>=0;s--) {
       auto part_iajv = getStateVector(pint_iajv,s,false);   
       auto part_v    = getStateVector(pint_v,s,false);   
       auto part_u    = getStateVector(pint_u,s,true);   
       auto part_z    = getControlVector(pint_z,s); 

       pint_iajv.getVectorPtr(2*s+1)->scale(1.0/globalScale_);
       pint_v.getVectorPtr(2*s)->axpy(globalScale_,*pint_iajv.getVectorPtr(2*s+1));

       // compute the constraint for this subdomain
       auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s),level,s);

       constraint->applyInverseAdjointJacobian_1(*part_iajv,*part_v,*part_u,*part_z,tol);

       sendBuffer[0] = part_iajv->get(0);
     }

     pint_iajv.boundaryExchangeRightToLeft(sendBuffer,PinTVector<Real>::SEND_ONLY);
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

   // Solve a linear system on this processor
   void invertTimeStepJacobian(PinTVector<Real>       & pint_ijv,
                               const PinTVector<Real> & pint_v,
                               const PinTVector<Real> & pint_u,
                               const PinTVector<Real> & pint_z,
                               Real &tol,
                               int level) 
   {
     int timeRank = pint_u.communicators().getTimeRank();
     int timeSize = pint_u.communicators().getTimeSize();
     bool lastRank = (timeRank+1==timeSize); // do something special on the last rank

     pint_v.boundaryExchangeLeftToRight();
     // pint_z.boundaryExchangeLeftToRight();
     // pint_u.boundaryExchangeLeftToRight();

     // fix up old data with previous time step information: This is the match to *** below
     pint_ijv.getRemoteBufferPtr(1)->set(*pint_v.getRemoteBufferPtr(1));
     pint_ijv.getRemoteBufferPtr(1)->scale(1.0/globalScale_);

     size_t numSteps = hierarchy_.getTimeStampsByLevel(level)->size();
     for(size_t s=0;s<numSteps;s++) { 

       auto part_ijv = getStateVector(pint_ijv,s,false);   
       auto part_v   = getStateVector(pint_v,s,false);   
       auto part_u   = getStateVector(pint_u,s,true);   
       auto part_z   = getControlVector(pint_z,s); 

       // compute the constraint for this subdomain
       auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s),level,s);

       constraint->applyInverseJacobian_1(*part_ijv,*part_v,*part_u,*part_z,tol);

       // satisfy the time continuity constraint, note that this just using the identity matrix
       if(s+1<numSteps || lastRank) {
         pint_ijv.getVectorPtr(2*s+1)->set(*pint_ijv.getVectorPtr(2*s));    
         pint_ijv.getVectorPtr(2*s+1)->axpy(1.0/globalScale_,*pint_v.getVectorPtr(2*s+1));        
       }
       else {
         // ***
         pint_ijv.getVectorPtr(2*s+1)->set(*pint_v.getVectorPtr(2*s+1));      
         pint_ijv.getVectorPtr(2*s+1)->scale(1.0/globalScale_);
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
     int timeRank = pint_u.communicators().getTimeRank();
     int timeSize = pint_u.communicators().getTimeSize();
     bool lastRank = (timeRank+1==timeSize); // do something special on the last rank

     // this is an inefficient hack (see line below where pint_v is modified!!!!)
     ///////////////////////////////////
     Ptr<Vector<Real>> v_copy;

     // handle the lazy construction by level (this is shared with the inverse adjoint call)
     {
       auto store = inverseAdjointStorage_.find(level);
       if(store==inverseAdjointStorage_.end()) {
         v_copy = pint_v_src.clone();
         inverseAdjointStorage_[level] = v_copy;
       }
       else {
         v_copy = store->second;
       }
     }
     v_copy->set(pint_v_src);

     PinTVector<Real> & pint_v = dynamic_cast<PinTVector<Real>&>(*v_copy);
     ///////////////////////////////////

     assert(pint_iajv.numOwnedSteps()==pint_u.numOwnedSteps());
     assert(pint_iajv.numOwnedSteps()==pint_v.numOwnedSteps());

     // pint_z.boundaryExchangeLeftToRight();
     // pint_u.boundaryExchangeLeftToRight();
     pint_v.boundaryExchangeLeftToRight();

     int numSteps = Teuchos::as<int>(hierarchy_.getTimeStampsByLevel(level)->size());
     if(lastRank) {
       // this is the last rank case
       pint_iajv.getVectorPtr(2*(numSteps-1)+1)->set(*pint_v.getVectorPtr(2*(numSteps-1)+1));
     }
     else {
       pint_iajv.getVectorPtr(2*(numSteps-1)+1)->zero();
     }

     std::vector<Ptr<Vector<Real>>> sendBuffer(1);
     for(int s=numSteps-1;s>=0;s--) {
       auto part_iajv = getStateVector(pint_iajv,s,false);   
       auto part_v    = getStateVector(pint_v,s,false);   
       auto part_u    = getStateVector(pint_u,s,true);   
       auto part_z    = getControlVector(pint_z,s); 

       pint_iajv.getVectorPtr(2*s+1)->scale(1.0/globalScale_);
       pint_v.getVectorPtr(2*s)->axpy(globalScale_,*pint_iajv.getVectorPtr(2*s+1));

       // compute the constraint for this subdomain
       auto constraint = getSerialConstraint(pint_u.getVectorPtr(2*s),level,s);

       constraint->applyInverseAdjointJacobian_1(*part_iajv,*part_v,*part_u,*part_z,tol);

       sendBuffer[0] = part_iajv->get(0);
     }

     pint_iajv.boundaryExchangeRightToLeft(sendBuffer);

     if(not lastRank) {
       // grab the pre-computed solution for the right hand side
       pint_iajv.getVectorPtr(2*(numSteps-1)+1)->set(*pint_iajv.getRemoteBufferPtr(0));
       pint_iajv.getVectorPtr(2*(numSteps-1)+1)->scale(1.0/globalScale_);
     }
   }
 
   virtual void applyPreconditioner(Vector<Real> &pv,
                                    const Vector<Real> &v,
                                    const Vector<Real> &x,
                                    const Vector<Real> &g,
                                    Real &tol) override
   {
     assert(false);
   }

   /////////////////////////////////////////////////////////////////////////////////////////
   // KKT multigrid preconditioner
   /////////////////////////////////////////////////////////////////////////////////////////
   
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

     auto timer = Teuchos::TimeMonitor::getStackedTimer();

     std::string levelStr = "";
     {
       std::stringstream ss;
       ss << "-" << level;
       levelStr = ss.str();
     }

     timer->start("applyAugmentedKKT"+levelStr); 

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
     output_u->axpy(1.0,*input_u); // TODO: Add Riesz scaling

     applyAdjointJacobian_2_leveled(*output_z,*input_v,u,z,tol,level);
     output_z->axpy(controlRegParam_,*input_z); // multiply by \alpha * I
                                   // TODO: Add Riesz scaling

     // constraint
     applyJacobian_1_leveled(*output_v_tmp,*input_u,u,z,tol,level);
     applyJacobian_2_leveled(*output_v,*input_z,u,z,tol,level);   // BAD ???

     output_v->axpy(1.0,*output_v_tmp);

     timer->stop("applyAugmentedKKT"+levelStr); 
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
   void applyWathenInverse(Vector<Real> & output, 
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

     timer->start("applyWathenInverse"+levelStr); 

     auto part_output = dynamic_cast<PartitionedVector&>(output);
     auto part_input = dynamic_cast<const PartitionedVector&>(input);
 
     part_output.zero();
 
     auto output_u = part_output.get(0);
     auto output_z = part_output.get(1);
     auto output_v = part_output.get(2);

     auto input_u  = part_input.get(0);
     auto input_z  = part_input.get(1);
     auto input_v  = part_input.get(2);

     Ptr<Vector<Real>> temp_u, temp_z, temp_v, temp_schur;

     auto store = inverseKKTStorage_.find(level);
     if(store==inverseKKTStorage_.end()) {
       temp_u = output_u->clone();
       temp_z = output_z->clone();
       temp_v = output_v->clone();
       temp_schur = output_v->clone();
       
       WathenInverseStorage data;
       data.temp_u = temp_u;
       data.temp_z = temp_z;
       data.temp_v = temp_v;
       data.temp_schur = temp_schur;

       inverseKKTStorage_[level] = data;
     }
     else {
       temp_u = store->second.temp_u;
       temp_z = store->second.temp_z;
       temp_v = store->second.temp_v;
       temp_schur = store->second.temp_schur;
     }
 
     temp_u->zero();
     temp_z->zero();
     temp_v->zero();
 
     // ==============================
     // [ I         0  J' * inv(J*J') ] [  I              ]
     // [     I/alpha  K' * inv(J*J') ] [  0     I        ]  // why is there no alpha parameter on the K' term?????
     // [                  -inv(J*J') ] [ -J -K/alpha  I  ]
    
     // with riesz maps (TODO)
     // ==============================
     // [ iMu          0        iMu * J' * inv(J*iMu*J') ] [  I                      ]
     // [      iMz/alpha  iMz * K' * inv(J*iMu*J')/alpha ] [  0           I          ]
     // [                                 -inv(J*iMu*J') ] [ -J*iMu -K/alpha*iMz  I  ]
     //
     // iMu - state Riesz map
     // iMz - control Riesz map
     // ==============================
    
     // L Factor
     /////////////////////
     temp_u->axpy(1.0,*input_u);
     temp_z->axpy(1.0/controlRegParam_,*input_z); // controlRegParam seems like a bug...
                                                  // its not! temp_z is not used or modified
 
     // apply -J   TODO: Apply inverse Riesz map to temp_v
     if(not approx)
       applyJacobian_1_leveled(*temp_v,*input_u,u,z,tol,level);
     else
       applyJacobian_1_leveled_approx(*temp_v,*input_u,u,z,tol,level);

     // apply -K (not yet, implies (2,1) block of L is not applied)
 
     temp_v->scale(-1.0);
     temp_v->axpy(1.0,*input_v);
 
     // U Factor
     /////////////////////
     
     // schur complement (Wathen style)
     {
       temp_schur->zero();
 
       if(not approx) {
         applyInverseJacobian_1_leveled(*temp_schur, *temp_v, u,z,tol,level);
         // TODO: Apply inverse Riesz map
         applyInverseAdjointJacobian_1_leveled(*output_v,*temp_schur,u,z,tol,level);
       }
       else {
         invertTimeStepJacobian(*temp_schur, *temp_v, u,z,tol,level);
         // TODO: Apply inverse Riesz map
         invertAdjointTimeStepJacobian(*output_v,*temp_schur,u,z,tol,level);
       }
       output_v->scale(-1.0);
     }

     if(not approx)
       applyAdjointJacobian_1_leveled(*output_u,*output_v,u,z,tol,level); 
     else 
       applyAdjointJacobian_1_leveled_approx(*output_u,*output_v,u,z,tol,level); 

     output_z->set(*temp_z); // above temp_z scaling by control parameter is included
     // TODO apply inverse Riesz map
 
     output_u->scale(-1.0);
     output_u->axpy(1.0,*temp_u);
     // TODO apply inverse Riesz map

     timer->stop("applyWathenInverse"+levelStr); 
   }


   void computeInvP(Vector<Real> & output, 
                    const Vector<Real> & input,
                    const Vector<Real> & u, 
                    const Vector<Real> & z,
                    Real & tol,
                    int level) 
   {
     auto temp_schur = output.clone();
     temp_schur->zero();

     invertTimeStepJacobian(*temp_schur, input, u,z,tol,level);
     invertAdjointTimeStepJacobian(output,*temp_schur,u,z,tol,level);

     output.scale(-1.0);
   }

   void applyLocalInverse(Vector<Real> & output, 
                          const Vector<Real> & input,
                          const Vector<Real> & u, 
                          const Vector<Real> & z,
                          Real & tol,
                          int level=0,
                          bool approx=true)
   {
     // TODO: Protect for use of Riesz maps

     using PartitionedVector = PartitionedVector<Real>;

     auto timer = Teuchos::TimeMonitor::getStackedTimer();

     std::string levelStr = "";
     {
       std::stringstream ss;
       ss << "-" << level;
       levelStr = ss.str();
     }

     timer->start("applyLocalInverse"+levelStr); 

     auto part_output = dynamic_cast<PartitionedVector&>(output);
     auto part_input = dynamic_cast<const PartitionedVector&>(input);
 
     part_output.zero();
 
     auto output_u = part_output.get(0);
     auto output_z = part_output.get(1);
     auto output_v = part_output.get(2);

     auto input_u  = part_input.get(0);
     auto input_z  = part_input.get(1);
     auto input_v  = part_input.get(2);

     Ptr<Vector<Real>> temp_u, temp_z, temp_v, temp_schur;

     auto store = inverseKKTStorage_.find(level);
     if(store==inverseKKTStorage_.end()) {
       temp_u = output_u->clone();
       temp_v = output_v->clone();
       temp_z = output_z->clone();
       temp_schur = output_v->clone();
       
       WathenInverseStorage data;
       data.temp_u = temp_u;
       data.temp_z = temp_z;
       data.temp_v = temp_v;
       data.temp_schur = temp_schur;

       inverseKKTStorage_[level] = data;
     }
     else {
       temp_u = store->second.temp_u;
       temp_v = store->second.temp_v;
       temp_z = store->second.temp_z;
       temp_schur = store->second.temp_schur;
     }
 
     temp_u->zero();
     temp_v->zero();
     temp_z->zero();

     //
     // [ I  J'         ]   [ I               ]   [ I  J'   ]
     // [ J           K ] = [ J     I         ] * [    P  K ]
     // [    K' alpha*I ]   [    K'*inv(P)  I ]   [       S ]
     //
     //   P = -J*J',   S = alpha*I - K'*inv(P)*K
     // 
   
     // L Factor
     /////////////////////

     // t_u = rhs_u
     temp_u->axpy(1.0,*input_u);
 
     // t_v = -J * t_u + f_v
     if(not approx)
       applyJacobian_1_leveled(*temp_v,*temp_u,u,z,tol,level);
     else
       applyJacobian_1_leveled_approx(*temp_v,*temp_u,u,z,tol,level);

     temp_v->scale(-1.0);
     temp_v->axpy(1.0,*input_v);

     // t_z = - K' * inv(P)*t_v + f_w
     {
       auto scratch = temp_schur;
       scratch->zero();

       computeInvP(*scratch,*temp_v,u,z,tol,level);
       applyAdjointJacobian_2_leveled(*temp_z,*scratch,u,z,tol,level);

       temp_z->scale(-1.0);
       temp_z->axpy(1.0,*input_z);
     }
 
     // U Factor
     /////////////////////

     // o_z = inv(S)*t_z
     output_z->zero();
     applyLocalReducedInverseHessian(*output_z,*temp_z,u,z,tol,level,approx);
     
     // o_v = inv(P)*(t_v-K*o_z)
     {
       auto scratch = temp_schur;
       scratch->zero();

       applyJacobian_2_leveled(*scratch,*output_z,u,z,tol,level);
       scratch->scale(-1.0);
       scratch->axpy(1.0,*temp_v);
       
       computeInvP(*output_v,*scratch,u,z,tol,level);
     }

     // o_u = t_u - J'*o_v
     if(not approx)
       applyAdjointJacobian_1_leveled(*output_u,*output_v,u,z,tol,level); 
     else
       applyAdjointJacobian_1_leveled_approx(*output_u,*output_v,u,z,tol,level); 
     output_u->scale(-1.0);
     output_u->axpy(1.0,*temp_u);

     timer->stop("applyLocalInverse"+levelStr); 
   }

   /**
    * The reduced Hessian operator for the augmented KKT sytem looks
    * like I-K'*inv(-J*J')*K. Set up this operator's action. 
    */
   void applyLocalReducedHessian(Vector<Real> & y, 
                                 const Vector<Real> & x,
                                 const Vector<Real> & u, 
                                 const Vector<Real> & z,
                                 Real & tol,
                                 int level,
                                 bool approx) 
   {
     auto temp_k = u.clone();
     auto temp_schur = u.clone();

     // K*x
     applyJacobian_2_leveled(*temp_k,x,u,z,tol,level);

     if(not approx) {
       // inv(J)*K*x
       applyInverseJacobian_1_leveled(*temp_schur,*temp_k, u,z,tol,level);

       // inv(J')*inv(J)*K*x
       applyInverseAdjointJacobian_1_leveled(*temp_k,*temp_schur,u,z,tol,level);
     }
     else {
       // inv(J)*K*x
       invertTimeStepJacobian(*temp_schur,*temp_k, u,z,tol,level);

       // inv(J')*inv(J)*K*x
       invertAdjointTimeStepJacobian(*temp_k,*temp_schur,u,z,tol,level);
     }

     // K'*inv(J')*inv(J)*K*x
     applyAdjointJacobian_2_leveled(y,*temp_k,u,z,tol,level);

     // (I+K'*inv(J')*inv(J)*K)*x
     y.axpy(controlRegParam_,x);
   }

   /**
    * The reduced Hessian operator for the augmented KKT sytem looks
    * like I-K'*inv(-J*J')*K. Compute its inverse using a hand rolled CG.
    */
   void applyLocalReducedInverseHessian(Vector<Real> & x, 
                                        const Vector<Real> & b,
                                        const Vector<Real> & u, 
                                        const Vector<Real> & z,
                                        Real & tol,
                                        int level,
                                        bool skipInitialMatVec=true,  // if x==0, no need
                                        bool approx=true) 
   {
     const PinTVector<Real> & pint_u = dynamic_cast<const PinTVector<Real>&>(u);
     int timeRank = pint_u.communicators().getTimeRank();

     Ptr<Vector<Real>> r;
     Ptr<Vector<Real>> Ap;
     Ptr<Vector<Real>> p;

     auto store = cgReducedHessStorage_.find(level);
     if(store==cgReducedHessStorage_.end()) {
       r  = z.clone();   
       Ap = z.clone();   
       p  = z.clone();   

       CGReducedHessStorage data;
       data.temp_r  = r;
       data.temp_Ap = Ap;
       data.temp_p  = p;

       cgReducedHessStorage_[level] = data;
     }
     else {
       r  = store->second.temp_r;
       Ap = store->second.temp_Ap;
       p  = store->second.temp_p;
     }
     
     // r = b-A*x
     if(skipInitialMatVec) {
       // this performance optimization is really important as the
       // applyLocal is super expensive
       r->set(b);
     }
     else {
       applyLocalReducedHessian(*r,x,u,z,tol,level,approx);
       r->axpy(-1.0,b);
       r->scale(-1.0);
     }

     // rold = r.r
     Real rold = r->dot(*r);
     Real rorg = rold;

     if(timeRank==0)
       std::cout << "START CG = " << std::sqrt(rorg) << std::endl;

     p->set(*r);
     
     int iters = numCGIter_;
     for(int i=0;i<iters;i++) {
       applyLocalReducedHessian(*Ap,*p,u,z,tol,level,approx);

       Real pAp = Ap->dot(*p);
       Real alpha = rold / pAp;

       x.axpy(alpha,*p);

       // exist early
       // if(i==iters-1)
       //   break;

       r->axpy(-alpha,*Ap);

       Real rnew = r->dot(*r);

       if(timeRank==0)
         std::cout << "CG Residual " << i << " = " << std::sqrt(rnew / rorg)  << std::endl;

       p->scale(rnew/rold);
       p->axpy(1.0,*r);

       rold = rnew;
     }

     if(timeRank==0)
       std::cout << "CG Residual Reduction  = " << std::sqrt(rold / rorg)  << std::endl;
   }
   

   /**
    * Apply the specified relaxation given an initial condition (output value)
    * x and the righthand side b. Returns, if requested, the residual reductions
    * are maintened.
    *
    * \param [in,out] x                    Solution vector
    * \param [in]     b                    RHS vector
    * \param [in]     u                    State vector
    * \param [in]     z                    Control vector
    * \param [in]     tol                  ROL tolerance
    * \param [in]     level                Multigrid level (level=0 is the fine level)
    * \param [in]     assumeZeroInitial    Assume the zero initial, saves a mat-vec
    * \param [in]     computeFinalResidual Compute the final residual before exiting
    */ 
   double applySmoother(Vector<Real> & x, 
                        const Vector<Real> & b,
                        const Vector<Real> & u, 
                        const Vector<Real> & z,
                        Real & tol,
                        int level,
                        bool assumeZeroInitial,
                        bool computeFinalResidual,
                        int numSweeps) 
   {
     bool approxSmoother = true;

     auto store = mgAugmentedKKTStorage_.find(level);
     Ptr<Vector<Real>> dx       = store->second.dx;
     Ptr<Vector<Real>> residual = store->second.residual;

     double res_norm = 0.0;
     if(assumeZeroInitial) {
       residual->set(b);
     }
     else {
       applyAugmentedKKT(*residual,x,u,z,tol,level);
       residual->scale(-1.0);
       residual->axpy(1.0,b);
     }

     if(recordResidualReductions_)
       res_norm = residual->norm();

     // apply one smoother sweep
     for(int i=0;i<numSweeps;i++) {

       if(numCGIter_==0)
         applyWathenInverse(*dx,*residual,u,z,tol,approxSmoother,level); 
       else
         applyLocalInverse(*dx,*residual,u,z,tol,level); 

       x.axpy(omega_,*dx);

       // compute the residual
       applyAugmentedKKT(*residual,x,u,z,tol,level);
 
       if(i<numSweeps-1 || computeFinalResidual || recordResidualReductions_) {
         residual->scale(-1.0);
         residual->axpy(1.0,b);
       }
     }

     if(recordResidualReductions_) {
       std::cout << "RES-Norm " <<  residual->norm()/res_norm << std::endl;
       return residual->norm()/res_norm;
     }

     return -1.0;
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

     const PinTVector<Real>       & pint_u = dynamic_cast<const PinTVector<Real>&>(u);
     const PinTVector<Real>       & pint_z = dynamic_cast<const PinTVector<Real>&>(z);

     pint_u.boundaryExchangeLeftToRight();
     pint_z.boundaryExchangeLeftToRight();

     std::string levelStr = "";
     std::string levelRankStr = "";
     {
       std::stringstream ss;
       ss << "-" << level;
       levelStr = ss.str();
     }

     timer->start("applyMGAugmentedKKT"+levelStr);

     Ptr<Vector<Real>> dx;
     Ptr<Vector<Real>> residual;

     auto store = mgAugmentedKKTStorage_.find(level);
     if(store==mgAugmentedKKTStorage_.end()) {
       dx       = x.clone();
       residual = b.clone();

       MGAugmentedKKTStorage data;
       data.dx       = dx;
       data.residual = residual;

       mgAugmentedKKTStorage_[level] = data;
     }
     else {
       dx       = store->second.dx;
       residual = store->second.residual;
     }

     std::string levelIndent = "";

     // base case: solve the KKT system directly
     if(level+1==maxLevels_) {
       bool approxSmoother = false;

       residual->set(b);

       double relax = omegaCoarse_;

       double res_norm = -1.0;
       if(recordResidualReductions_) 
         res_norm = residual->norm();

       for(int i=0;i<numCoarseSweeps_;i++) {
         // compute the residual
         applyAugmentedKKT(*residual,x,u,z,tol,level);
         residual->scale(-1.0);
         residual->axpy(1.0,b);

         if(numCGIter_==0)
           applyWathenInverse(*dx,*residual,u,z,tol,approxSmoother,level); 
         else
           applyLocalInverse(*dx,*residual,u,z,tol,level,approxSmoother); 
         x.axpy(relax,*dx);
       }

       if(recordResidualReductions_) {
         // compute the residual
         applyAugmentedKKT(*residual,x,u,z,tol,level);
         residual->scale(-1.0);
         residual->axpy(1.0,b);

         coarseResidualReduction_.push_back(residual->norm()/res_norm);
       }

       timer->stop("applyMGAugmentedKKT"+levelStr);

       return;
     }

     timer->start("applyMGAugmentedKKT-preSmooth");

     // pre-smooth
     /////////////////////////////////////////////////////////////////////////////////
     {
       bool finalResidualRequired  = true;
       bool assumeZeroInitialGuess = true;
       double relativeResidual = applySmoother(x,b,u,z,tol,level,assumeZeroInitialGuess,finalResidualRequired,numSweeps_);
       preSmoothResidualReduction_[level].push_back(relativeResidual);
     }

     timer->stop("applyMGAugmentedKKT-preSmooth");

     // solve the coarse system
     timer->start("applyMGAugmentedKKT-coarse");
     {
       auto pint_u = ROL::makePtrFromRef(dynamic_cast<const PinTVector<Real>&>(u));
       auto pint_z = ROL::makePtrFromRef(dynamic_cast<const PinTVector<Real>&>(z));

       auto dx_u = dynamicPtrCast<PinTVector<Real>>(dynamic_cast<PartitionedVector&>(*dx).get(0));
       auto dx_z = dynamicPtrCast<PinTVector<Real>>(dynamic_cast<PartitionedVector&>(*dx).get(1));
       auto dx_v = dynamicPtrCast<PinTVector<Real>>(dynamic_cast<PartitionedVector&>(*dx).get(2));

       auto residual_u = dynamicPtrCast<PinTVector<Real>>(dynamic_cast<PartitionedVector&>(*residual).get(0));
       auto residual_z = dynamicPtrCast<PinTVector<Real>>(dynamic_cast<PartitionedVector&>(*residual).get(1));
       auto residual_v = dynamicPtrCast<PinTVector<Real>>(dynamic_cast<PartitionedVector&>(*residual).get(2));

       auto crs_u            = hierarchy_.allocateSimVector(*pint_u,level+1);

       auto crs_residual_u   = hierarchy_.allocateSimVector(*pint_u,level+1);
       auto crs_residual_v   = hierarchy_.allocateSimVector(*pint_u,level+1);
       auto crs_correction_u = hierarchy_.allocateSimVector(*pint_u,level+1);
       auto crs_correction_v = hierarchy_.allocateSimVector(*pint_u,level+1);

       auto crs_z            = hierarchy_.allocateOptVector(*pint_z,level+1);
       auto crs_residual_z   = hierarchy_.allocateOptVector(*pint_z,level+1);
       auto crs_correction_z = hierarchy_.allocateOptVector(*pint_z,level+1);

       hierarchy_.restrictSimVector(residual_u,crs_residual_u,level);
       hierarchy_.restrictOptVector(residual_z,crs_residual_z,level);
       hierarchy_.restrictSimVector(residual_v,crs_residual_v,level);

       hierarchy_.restrictSimVector(pint_u,crs_u,level);               // restrict the state to the coarse level
       hierarchy_.restrictOptVector(pint_z,crs_z,level);               // restrict the control to the coarse level

       if(hierarchy_.levelIsActiveOnMyRank(level+1)) {
         typedef std::vector<ROL::Ptr<ROL::Vector<Real>>> vector;

         auto crs_correction = makePtr<PartitionedVector>(vector({crs_correction_u,crs_correction_z,crs_correction_v}));
         auto crs_residual   = makePtr<PartitionedVector>(vector({  crs_residual_u,  crs_residual_z,  crs_residual_v}));

         applyMultigridAugmentedKKT(*crs_correction,*crs_residual,*crs_u,*crs_z,tol,level+1);
       }

       hierarchy_.prolongSimVector(crs_correction_u,dx_u,level+1);
       hierarchy_.prolongOptVector(crs_correction_z,dx_z,level+1);
       hierarchy_.prolongSimVector(crs_correction_v,dx_v,level+1);

       x.axpy(1.0,*dx);
     }
     timer->stop("applyMGAugmentedKKT-coarse");

     // apply one smoother sweep
     timer->start("applyMGAugmentedKKT-postSmooth");

     // post-smooth
     /////////////////////////////////////////////////////////////////////////////////
     {
       bool finalResidualRequired = false;
       bool assumeZeroInitialGuess     = false;
       double relativeResidual = applySmoother(x,b,u,z,tol,level,assumeZeroInitialGuess,finalResidualRequired,numSweeps_);
       postSmoothResidualReduction_[level].push_back(relativeResidual);
     }

     timer->stop("applyMGAugmentedKKT-postSmooth");

     timer->stop("applyMGAugmentedKKT"+levelStr);
   }

}; // ROL::PinTConstraint


} // namespace ROL 

#endif // ROL_PINTCONSTRAINT_HPP

