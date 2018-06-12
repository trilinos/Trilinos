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

  // Get the state vector required by the serial constraint object
  Ptr<SerialConstraint<Real>> 
  getSerialConstraint(const Ptr<ROL::Vector<Real>> & ui,int Nt)
  {
    if(ROL::is_nullPtr(timeDomainConstraint_)) {
      timeDomainConstraint_ = ROL::makePtr<SerialConstraint<Real>>(stepConstraint_,*ui,timeStamps_);
        // FIXME-A: Change the Nt+1 back to Nt...
    }
    else {
      timeDomainConstraint_->setTimeStamp(timeStamps_);
        // FIXME-A: Change the Nt+1 back to Nt...
      timeDomainConstraint_->setInitialCondition(*ui);
    }

    timeDomainConstraint_->setSkipInitialCondition(true);
    
    return timeDomainConstraint_;
  }

public:

  //! Default constructor
  PinTConstraint()
    : Constraint_SimOpt<Real>()
    , isInitialized_(false)
  { }

  /**
   * \brief Constructor
   *
   * Build a parallel-in-time constraint with a specified step constraint. This specifies
   * any communication you might need between processors and steps using "stencils".
   * 
   * \param[in] stepConstraint Constraint for a single step.
   * \param[in] initialCond Initial condition
   */
  PinTConstraint(const Ptr<DynamicConstraint<Real>> & stepConstraint,
                 const Ptr<const Vector<Real>> & initialCond,
                 const Ptr<std::vector<TimeStamp<Real>>> & timeStamps)
    : Constraint_SimOpt<Real>()
    , isInitialized_(false)
  { 
    initialize(stepConstraint,initialCond,timeStamps);
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
    
    PinTVector<Real>       & pint_c = dynamic_cast<PinTVector<Real>&>(c);
    PinTVector<Real>       & pint_u = dynamic_cast<PinTVector<Real>&>(u);
    const PinTVector<Real> & pint_z = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here
 
    TEUCHOS_ASSERT(pint_c.numOwnedSteps()==pint_u.numOwnedSteps());

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
      }
      else if(offset<0 and timeRank==0) {  
        // this is the intial condition

        // update the old time information
        pint_u.getVectorPtr(offset)->set(*initialCond_);        

        // this should be zero!
        pint_c.getVectorPtr(offset)->set(*pint_u.getVectorPtr(offset));
        pint_c.getVectorPtr(offset)->axpy(-1.0,*initialCond_);
      }
    }

    auto part_c = getStateVector(pint_c);   
    auto part_u = getStateVector(pint_u);   
    auto part_z = getControlVector(pint_z); 

    // compute the constraint for this subdomain
    auto constraint = getSerialConstraint(pint_u.getVectorPtr(-1),pint_u.numOwnedSteps());

    constraint->solve(*part_c,*part_u,*part_z,tol);

    pint_u.boundaryExchange(PinTVector<Real>::SEND_ONLY);
  }  

  void value( V& c, const V& u, const V& z, Real& tol ) override {

    c.zero();   

    PinTVector<Real>       & pint_c = dynamic_cast<PinTVector<Real>&>(c);
    const PinTVector<Real> & pint_u = dynamic_cast<const PinTVector<Real>&>(u);
    const PinTVector<Real> & pint_z = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here
       
    TEUCHOS_ASSERT(pint_c.numOwnedSteps()==pint_u.numOwnedSteps());

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
      }
      else if(offset<0 and timeRank==0) {  
        // set the initial condition
        pint_c.getVectorPtr(offset)->set(*pint_u.getVectorPtr(offset));            
        pint_c.getVectorPtr(offset)->axpy(-1.0,*initialCond_); 
      }
    }

    auto part_c = getStateVector(pint_c);   // strip out initial condition/constraint
    auto part_u = getStateVector(pint_u);   // strip out initial condition/constraint
    auto part_z = getControlVector(pint_z); 

    // compute the constraint for this subdomain
    auto constraint = getSerialConstraint(pint_u.getVectorPtr(-1),pint_u.numOwnedSteps());

    constraint->value(*part_c,*part_u,*part_z,tol);
  } 

  void applyJacobian_1( V& jv, const V& v, const V& u, 
                       const V& z, Real& tol ) override {

    jv.zero();

    PinTVector<Real>       & pint_jv = dynamic_cast<PinTVector<Real>&>(jv);
    const PinTVector<Real> & pint_v  = dynamic_cast<const PinTVector<Real>&>(v);
    const PinTVector<Real> & pint_u  = dynamic_cast<const PinTVector<Real>&>(u);
    const PinTVector<Real> & pint_z  = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here
 
    TEUCHOS_ASSERT(pint_jv.numOwnedSteps()==pint_u.numOwnedSteps());
    TEUCHOS_ASSERT(pint_jv.numOwnedSteps()==pint_v.numOwnedSteps());

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
      }
    }

    auto part_jv = getStateVector(pint_jv);   // strip out initial condition/constraint
    auto part_v  = getStateVector(pint_v);    // strip out initial condition/constraint
    auto part_u  = getStateVector(pint_u);    // strip out initial condition/constraint
    auto part_z  = getControlVector(pint_z); 

    // compute the constraint for this subdomain
    auto constraint = getSerialConstraint(pint_u.getVectorPtr(-1),pint_u.numOwnedSteps());

    constraint->applyJacobian_1(*part_jv,*part_v,*part_u,*part_z,tol);
  }

  void applyJacobian_2( V& jv, const V& v, const V& u,
                        const V &z, Real &tol ) override { 

    jv.zero();

    PinTVector<Real>       & pint_jv = dynamic_cast<PinTVector<Real>&>(jv);
    const PinTVector<Real> & pint_v  = dynamic_cast<const PinTVector<Real>&>(v);
    const PinTVector<Real> & pint_u  = dynamic_cast<const PinTVector<Real>&>(u);
    const PinTVector<Real> & pint_z  = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here
 
    TEUCHOS_ASSERT(pint_jv.numOwnedSteps()==pint_u.numOwnedSteps());
    TEUCHOS_ASSERT(pint_jv.numOwnedSteps()==pint_v.numOwnedSteps());

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
    auto constraint = getSerialConstraint(pint_u.getVectorPtr(-1),pint_jv.numOwnedSteps());

    constraint->applyJacobian_2(*part_jv,*part_v,*part_u,*part_z,tol);
   }


   void applyInverseJacobian_1( V& ijv, const V& v, const V& u,
                                const V& z, Real& tol) override {

   }


   void applyAdjointJacobian_1( V& ajv, 
                                const V& v, 
                                const V& u,
                                const V& z, 
                                Real& tol) override {
    PinTVector<Real>       & pint_ajv = dynamic_cast<PinTVector<Real>&>(ajv);
    const PinTVector<Real> & pint_v   = dynamic_cast<const PinTVector<Real>&>(v);
    const PinTVector<Real> & pint_u   = dynamic_cast<const PinTVector<Real>&>(u);
    const PinTVector<Real> & pint_z   = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here
      
    TEUCHOS_ASSERT(pint_v.numOwnedSteps()==pint_u.numOwnedSteps());
    TEUCHOS_ASSERT(pint_ajv.numOwnedSteps()==pint_u.numOwnedSteps());

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
    auto constraint = getSerialConstraint(pint_u.getVectorPtr(-1),pint_u.numOwnedSteps());

    constraint->applyAdjointJacobian_1(*part_ajv,*part_v,*part_u,*part_z,*part_v,tol);

    // handle the constraint adjoint
    const std::vector<int> & stencil = pint_u.stencil();
    for(size_t i=0;i<stencil.size();i++) {
      int offset = stencil[i];
      if(offset<0) {
        pint_ajv.getVectorPtr(offset)->axpy(1.0,*pint_v.getVectorPtr(offset));
        pint_ajv.getRemoteBufferPtr(offset)->set(*pint_v.getVectorPtr(offset)); // this will be sent to the remote processor
        pint_ajv.getRemoteBufferPtr(offset)->scale(-1.0);
      }
    }

    // this sums from the remote buffer into the local buffer which is part of the vector
    pint_ajv.boundaryExchangeSumInto();
   }

   void applyAdjointJacobian_2( V& ajv,  const V& v, const V& u,
                                const V& z, Real& tol ) override {

    PinTVector<Real>       & pint_ajv = dynamic_cast<PinTVector<Real>&>(ajv);
    const PinTVector<Real> & pint_v   = dynamic_cast<const PinTVector<Real>&>(v);
    const PinTVector<Real> & pint_u   = dynamic_cast<const PinTVector<Real>&>(u);
    const PinTVector<Real> & pint_z   = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here
      
    TEUCHOS_ASSERT(pint_v.numOwnedSteps()==pint_u.numOwnedSteps());
    TEUCHOS_ASSERT(pint_ajv.numOwnedSteps()==pint_u.numOwnedSteps());

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
    auto constraint = getSerialConstraint(pint_u.getVectorPtr(-1),pint_u.numOwnedSteps());

    constraint->applyAdjointJacobian_2(*part_ajv,*part_v,*part_u,*part_z,tol);

    // no constraint for controls hanlding because that doesn't work yet!
   }

   void applyInverseAdjointJacobian_1( V& iajv, const V& v, const V& u,
                                       const V& z, Real& tol) override {

   }

   // Done in parallel, no blocking, solve a linear system on this processor
   void invertTimeStepJacobian(PinTVector<Real>       & pint_pv,
                               const PinTVector<Real> & pint_v,
                               const PinTVector<Real> & pint_u,
                               const PinTVector<Real> & pint_z,
                               Real &tol) 
   {
     // apply the time continuity constraint, note that this just using the identity matrix
     const std::vector<int> & stencil = pint_pv.stencil();
     for(size_t i=0;i<stencil.size();i++) {
       int offset = stencil[i];
       if(offset<0) {
         pint_pv.getVectorPtr(offset)->set(*pint_v.getVectorPtr(offset));             // this is the local value
       }
     }

     auto part_pv  = getStateVector(pint_pv);   // strip out initial condition/constraint
     auto part_v   = getStateVector(pint_v);   
     auto part_u   = getStateVector(pint_u);    // strip out initial condition/constraint
     auto part_z   = getControlVector(pint_z); 

     // compute the constraint for this subdomain
     auto constraint = getSerialConstraint(pint_u.getVectorPtr(-1),pint_u.numOwnedSteps());

     constraint->applyInverseJacobian_1(*part_pv,*part_v,*part_u,*part_z,tol);
   }
 
   // Done in parallel, no blocking, solve a linear system on this processor
   void invertAdjointTimeStepJacobian(PinTVector<Real>       & pint_pv,
                                      const PinTVector<Real> & pint_v,
                                      const PinTVector<Real> & pint_u,
                                      const PinTVector<Real> & pint_z,
                                      Real &tol) 
   {
     auto part_pv  = getStateVector(pint_pv);   // strip out initial condition/constraint
     auto part_v   = getStateVector(pint_v);   
     auto part_u   = getStateVector(pint_u);    // strip out initial condition/constraint
     auto part_z   = getControlVector(pint_z); 

     // compute the constraint for this subdomain
     auto constraint = getSerialConstraint(pint_u.getVectorPtr(-1),pint_u.numOwnedSteps());

     constraint->applyInverseAdjointJacobian_1(*part_pv,*part_v,*part_u,*part_z,tol);
    
     // apply the time continuity constraint, note that this just using the identity matrix
     const std::vector<int> & stencil = pint_pv.stencil();
     for(size_t i=0;i<stencil.size();i++) {
       int offset = stencil[i];
       if(offset<0) {
         // pint_pv.getVectorPtr(offset)->set(*temp);
           // NOT good enough for multi step!!!!
       }
     }
   }
 
   void blockJacobiSweep(PinTVector<Real>       & pint_pv,
                         const PinTVector<Real> & pint_rhs,
                         const PinTVector<Real> & pint_u,
                         const PinTVector<Real> & pint_z,
                         Real &tol) 
   {
     pint_rhs.boundaryExchange();
 
     Ptr<Vector<Real>> scratch = pint_rhs.clone();
     PinTVector<Real> pint_scratch  = dynamic_cast<const PinTVector<Real>&>(*scratch);
 
     // do block Jacobi smoothing
     invertTimeStepJacobian(pint_scratch, pint_rhs, pint_u,pint_z,tol);
 
     invertAdjointTimeStepJacobian(pint_pv, pint_scratch, pint_u,pint_z,tol);
   }
 
   void schurComplementAction(Vector<Real>       & sc_v,
                              const Vector<Real> & v,
                              const Vector<Real> & u,
                              const Vector<Real> & z,
                              Real &tol) 
   {
     auto scratch_u = u.clone(); 
     auto scratch_z = z.clone(); 
 
     // do boundary exchange
     {
       const PinTVector<Real>       & pint_v = dynamic_cast<const PinTVector<Real>&>(v);
       pint_v.boundaryExchange();
     }
 
     // compute J_1 * J_1^T
     applyAdjointJacobian_1(*scratch_u,v,u,z,tol);
     applyJacobian_1(sc_v,*scratch_u,u,z,tol);
 
     // compute J_2 * J_2^T
     applyAdjointJacobian_2(*scratch_z,v,u,z,tol);
     applyJacobian_2(*scratch_u,*scratch_z,u,z,tol);
     
     // compute J_1 * J_1^T + J_2 * J_2^T
     sc_v.axpy(1.0,*scratch_u);
   }
 
   virtual void weightedJacobiRelax(PinTVector<Real> &pint_pv,
                                    const PinTVector<Real> & pint_v,
                                    const PinTVector<Real> & pint_u,
                                    const PinTVector<Real> & pint_z,
                                    Real omega,
                                    int numSweeps,
                                    Real &tol) 
   {
     auto residual = pint_v.clone();
     auto dx = pint_pv.clone();
 
     // force intial guess to zero
     pint_pv.scale(0.0);
     dx->scale(0.0);
     residual->set(pint_v);           
 
     PinTVector<Real> & pint_residual = dynamic_cast<PinTVector<Real>&>(*residual);
     PinTVector<Real> & pint_dx       = dynamic_cast<PinTVector<Real>&>(*dx);
 
     for(int s=0;s<numSweeps;s++) {
       blockJacobiSweep(pint_dx,pint_residual,pint_u,pint_z,tol);
 
       pint_pv.axpy(omega,pint_dx);
 
       // now update the residual
       if(s!=numSweeps-1) {
         schurComplementAction(pint_residual,pint_pv,pint_u,pint_z,tol); // A*x
 
         pint_residual.scale(-1.0);                            // -A*x
         pint_residual.axpy(-1.0,pint_v);                      // b - A*x
       }
     }
 
     pint_pv.boundaryExchange();
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
  
     TEUCHOS_ASSERT(pint_pv.numOwnedSteps()==pint_u.numOwnedSteps());
     TEUCHOS_ASSERT(pint_pv.numOwnedSteps()==pint_v.numOwnedSteps());
 
     // these don't change, do boundary exchange
     pint_u.boundaryExchange();
     pint_z.boundaryExchange();
 
     // now we are going to do a relaxation sweep: x = xhat+omega * Minv * (b-A*xhat)
     /////////////////////////////////////////////////////////////////////////////////////////
     
     int numSweeps = 4;
     Real omega = 2.0/3.0;
     
     weightedJacobiRelax(pint_pv,pint_v,pint_u,pint_z,omega,numSweeps,tol);
     
     // pv.set(v.dual());
   }


}; // ROL::PinTConstraint


} // namespace ROL 

#endif // ROL_PINTCONSTRAINT_HPP

