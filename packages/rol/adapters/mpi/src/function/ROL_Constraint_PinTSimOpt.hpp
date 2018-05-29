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

#ifndef ROL_CONSTRAINT_PINTSIMOPT_H
#define ROL_CONSTRAINT_PINTSIMOPT_H

#include <cassert>

#include "ROL_PinTVector.hpp"
#include "ROL_Constraint_SimOpt.hpp"

/** @ingroup func_group
    \class ROL::Constraint_TimeSimOpt
    \brief Defines the time dependent constraint operator interface for simulation-based optimization.

    This constraint interface inherits from ROL_Constraint_SimOpt. Though the interface
    takes two simulation space vectors from spaces
    \f$\mathcal{U_o}\times\mathcal{U_n}\f$. The space \f$\mathcal{U_o}\f$ is ``old'' information
    that accounts for the initial condition on the time interval. The space \f$\mathcal{U_n}\f$ is the
    ``new'' variables that can be determined by satisfying constraints in the form
    \f[
      c(u_o,u_n,z) = 0 \,.
    \f]
    where \f$u_0 \in \mathcal{U_o},\; u_n\in\mathcal{U_n},\f$ and \f$z\in\mathcal{Z}\f$. In this way
    this constraint defines a sequence of state variables.
    The basic operator interface, to be implemented by the user, requires:
    \li #value -- constraint evaluation.
    \li #applyJacobian_1_old         -- action of the partial constraint Jacobian --derivatives are 
                                        with respect to the first component \f$\mathcal{U_o}\f$;
    \li #applyJacobian_1_new         -- action of the partial constraint Jacobian --derivatives are 
                                        with respect to the first component \f$\mathcal{U_n}\f$;
    \li #applyJacobian_2             -- action of the partial constraint Jacobian --derivatives are 
                                        with respect to the second component \f$\mathcal{Z}\f$;
    \li #applyAdjointJacobian_1_old  -- action of the adjoint of the partial constraint Jacobian --derivatives are 
                                        with respect to the first component \f$\mathcal{U_o}\f$;
    \li #applyAdjointJacobian_1_new  -- action of the adjoint of the partial constraint Jacobian --derivatives are 
                                        with respect to the first component \f$\mathcal{U_n}\f$;
    \li #applyAdjointJacobian_2_time -- action of the adjoint of the partial constraint Jacobian --derivatives are 
                                        with respect to the second component \f$\mathcal{Z}\f$; (note the time
                                        suffix here is to prevent collisions with the parent class)

    The user may also overload:
    \li #applyAdjointHessian_11  -- action of the adjoint of the partial constraint Hessian --derivatives 
                                    are with respect to the first component only;
    \li #applyAdjointHessian_12  -- action of the adjoint of the partial constraint Hessian --derivatives 
                                    are with respect to the first and second components;
    \li #applyAdjointHessian_21  -- action of the adjoint of the partial constraint Hessian --derivatives 
                                    are with respect to the second and first components;
    \li #applyAdjointHessian_22  -- action of the adjoint of the partial constraint Hessian --derivatives 
                                    are with respect to the second component only;
    \li #solveAugmentedSystem -- solution of the augmented system --the default is an iterative
                                 scheme based on the action of the Jacobian and its adjoint.
    \li #applyPreconditioner  -- action of a constraint preconditioner --the default is null-op.

    ---
*/


namespace ROL {


/** This helper method builds a pint vector for use in the 
    Constraint_PinTSimOpt class. Notice no difference is made between
    the state and contraint vectors.

    \param[in] communicators Structure with objects required for parallel-in-time communication.
                             (hint: this can be used to properly allocate the state and control vectors)
    \param[in] steps Number of steps 
    \param[in] stencil Stencil for each time step.
    \param[in] localVector Spatial vector for a single time step.
  */
template <class Real>
Ptr<PinTVector<Real>>
buildPinTVector(const Ptr<const PinTCommunicators> & communicators,
                int steps,
                const std::vector<int> & stencil,
                const Ptr<Vector<Real>> & localVector)
{ return makePtr<PinTVector<Real>>(communicators,localVector,steps,stencil); }

template <class Real>
class Constraint_PinTSimOpt : public Constraint_SimOpt<Real> {
private:
  Ptr<Constraint_TimeSimOpt<Real>> stepConstraint_;

  // internal state members
  
  bool isInitialized_;

  //! Enumeration to define which components of a a vector are required.
  typedef enum {PAST,CURRENT,FUTURE,ALL} ETimeAccessor;

  /** \brief Build a vector composed of all vectors required by accessor
   
      Build a vector composed of all vectors require by accessor. If this
      is only one vector, that vector is returned directly. If there is more
      than one vector than a partioned vector is returned with each entry.
   */
  Ptr<const Vector<Real>> getVector(const PinTVector<double> & src,
                                    int step,
                                    ETimeAccessor accessor)
  {
    return getNonconstVector(src,step,accessor);
  }

  /** \brief Build a vector composed of all vectors required by accessor
   
      Build a vector composed of all vectors require by accessor. If this
      is only one vector, that vector is returned directly. If there is more
      than one vector than a partioned vector is returned with each entry.
   */
  Ptr<Vector<Real>> getNonconstVector(const PinTVector<double> & src,
                                      int step,
                                      ETimeAccessor accessor)
  {
    const std::vector<int> stencil = src.stencil();

    std::vector<Ptr<Vector<Real>>> vecs;
    for(std::size_t i=0;i<stencil.size();i++) {

      // make sure this aligns with the accesor
      bool valid = false;
      if(accessor==PAST    && stencil[i]<0)       
        valid = true;
      if(accessor==CURRENT && stencil[i]==0)
        valid = true;
      if(accessor==FUTURE  && stencil[i]>0)
        valid = true;
      if(accessor==ALL)
        valid = true;
        
      // stencil entry matches with the accessor
      if(valid)
        vecs.push_back(src.getVectorPtr(step+stencil[i]));
    }

    TEUCHOS_ASSERT(vecs.size()>0);

    if(vecs.size()==1) 
      return vecs[0];

    return makePtr<PartitionedVector<Real>>(vecs);
  }

  Ptr<Vector<Real>> getStateVector(const PinTVector<double> & src,int step)
  {
      Ptr<Vector<Real>> u_old = getNonconstVector(src,step,PAST);
      Ptr<Vector<Real>> u_now = getNonconstVector(src,step,CURRENT);
    
      // this is specified by the Constraint_TimeSimOpt class
      std::vector<Ptr<Vector<Real>>> vecs;
      vecs.push_back(u_old);
      vecs.push_back(u_now);
      return makePtr<PartitionedVector<Real>>(vecs);
  }

public:

  //! Default constructor
  Constraint_PinTSimOpt()
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
   * \param[in] steps The number of steps to take.
   * \param[in] stepState State vector step value, also the initial condition
   * \param[in] stateStencil Stencil for parallel stepping.
   * \param[in] controlState Control vector step value, also the initial condition
   * \param[in] controlStencil Stencil for control parallel stepping.
   */
  Constraint_PinTSimOpt(const Ptr<Constraint_TimeSimOpt<Real>> & stepConstraint)
    : Constraint_SimOpt<Real>()
    , isInitialized_(false)
  { 
    initialize(stepConstraint);
  }

  /** \brief Initialize this class, setting up parallel distribution.
   
   */
  void initialize(const Ptr<Constraint_TimeSimOpt<Real>> & stepConstraint)
  {
    // initialize user member variables
    stepConstraint_ = stepConstraint;

    isInitialized_ = true;
  }

  /** \brief Update constraint functions with respect to Opt variable.
                u is the state variable, 
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update_1( const Vector<Real> &u, bool flag = true, int iter = -1 ) override {}

  /** \brief Update constraint functions with respect to Opt variable.
                z is the control variable, 
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update_2( const Vector<Real> &z, bool flag = true, int iter = -1 ) override {}

  virtual void value(Vector<Real> &c,
                     const Vector<Real> &u,
                     const Vector<Real> &z,
                     Real &tol) override {

    PinTVector<Real>       & pint_c = dynamic_cast<PinTVector<Real>&>(c);
    const PinTVector<Real> & pint_u = dynamic_cast<const PinTVector<Real>&>(u);
    const PinTVector<Real> & pint_z = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here
 
    TEUCHOS_ASSERT(pint_c.numOwnedSteps()==pint_u.numOwnedSteps());

    // communicate neighbors, these are block calls
    pint_u.boundaryExchange();
    pint_z.boundaryExchange();

    // build in the time continuity constraint, note that this just using the identity matrix
    const std::vector<int> & stencil = pint_c.stencil();
    for(size_t i=0;i<stencil.size();i++) {
      int offset = stencil[i];
      if(offset<0) {
        pint_c.getVectorPtr(offset)->set(*pint_u.getVectorPtr(offset));             // this is the local value
        pint_c.getVectorPtr(offset)->axpy(-1.0,*pint_u.getRemoteBufferPtr(offset)); // this is the remote value
      }
    }

    // now setup the residual for the time evolution
    for(int i=0;i<pint_c.numOwnedSteps();i++) {
      Ptr<const Vector<Real>> u_state = getStateVector(pint_u,i);
      Ptr<const Vector<Real>> z_all = getVector(pint_z,i,ALL);

      Ptr<Vector<Real>> c_now = getNonconstVector(pint_c,i,CURRENT);

      stepConstraint_->value(*c_now,*u_state,*z_all,tol);
    }
  } 

  virtual void solve(Vector<Real> &c,
                     Vector<Real> &u, 
                     const Vector<Real> &z,
                     Real &tol) override {
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
    for(size_t i=0;i<stencil.size();i++) {
      int offset = stencil[i];
      if(offset<0) {

        // update the old time information
        pint_u.getVectorPtr(offset)->set(*pint_u.getRemoteBufferPtr(offset));        // this is the local value

        // this should be zero!
        pint_c.getVectorPtr(offset)->set(*pint_u.getVectorPtr(offset));             // this is the local value
        pint_c.getVectorPtr(offset)->axpy(-1.0,*pint_u.getRemoteBufferPtr(offset)); // this is the remote value
      }
    }

    // solve the time steps
    for(int i=0;i<pint_c.numOwnedSteps();i++) {
      Ptr<const Vector<Real>> z_all = getVector(pint_z,i,ALL);

      Ptr<Vector<Real>> u_old = getNonconstVector(pint_u,i,PAST);
      Ptr<Vector<Real>> u_now = getNonconstVector(pint_u,i,CURRENT);
      Ptr<Vector<Real>> c_now = getNonconstVector(pint_c,i,CURRENT);

      // we solve for u_now
      stepConstraint_->solve(*c_now,*u_old,*u_now,*z_all,tol);
    }

    pint_u.boundaryExchange(PinTVector<Real>::SEND_ONLY);
  }

  virtual void applyJacobian_1(Vector<Real> &jv,
                               const Vector<Real> &v,
                               const Vector<Real> &u,
                               const Vector<Real> &z,
                               Real &tol) override {
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

    for(int i=0;i<pint_jv.numOwnedSteps();i++) {
      // pull out all the time steps required
      Ptr<const Vector<Real>> v_state = getStateVector(pint_v,i);
      Ptr<const Vector<Real>> u_state = getStateVector(pint_u,i);
      Ptr<const Vector<Real>> z_all = getVector(pint_z,i,ALL);

      Ptr<Vector<Real>> jv_now = getNonconstVector(pint_jv,i,CURRENT);

      stepConstraint_->applyJacobian_1(*jv_now,*v_state,*u_state,*z_all,tol);
    }
  }

  virtual void applyJacobian_2(Vector<Real> &jv,
                               const Vector<Real> &v,
                               const Vector<Real> &u,
                               const Vector<Real> &z,
                               Real &tol) override {
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

    const std::vector<int> & stencil = pint_jv.stencil();
    for(size_t i=0;i<stencil.size();i++) {
      int offset = stencil[i];
      if(offset<0) {
        pint_jv.getVectorPtr(offset)->zero();
      }
    }

    for(int i=0;i<pint_jv.numOwnedSteps();i++) {
      // pull out all the time steps required
      Ptr<const Vector<Real>> u_state   = getStateVector(pint_u,i);
      Ptr<const Vector<Real>> v_control = getVector(pint_v,i,ALL);
      Ptr<const Vector<Real>> z_all     = getVector(pint_z,i,ALL);

      Ptr<Vector<Real>> jv_now = getNonconstVector(pint_jv,i,CURRENT);

      stepConstraint_->applyJacobian_2(*jv_now,*v_control,*u_state,*z_all,tol);
    }
  }

  virtual void applyInverseJacobian_1(Vector<Real> &ijv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) override final {
    ROL_TEST_FOR_EXCEPTION(true, std::logic_error,
      "The method applyInverseJacobian_1 is used but not implemented!\n");
  }

  virtual void applyAdjointJacobian_1(Vector<Real> &ajv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) override {

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

    Ptr<Vector<Real>> scratch;
    for(int i=0;i<pint_ajv.numOwnedSteps();i++) {
      // pull out all the time steps required
      Ptr<const Vector<Real>> v_dual = getVector(pint_v,i,CURRENT);
      Ptr<const Vector<Real>> u_state = getStateVector(pint_u,i);
      Ptr<const Vector<Real>> z_all = getVector(pint_z,i,ALL);

      Ptr<Vector<Real>> ajv_now = getStateVector(pint_ajv,i);

      // this will lead to problems if problem size changes in time
      if(scratch==ROL::nullPtr) {
        scratch = ajv_now->clone();
      }

      stepConstraint_->applyAdjointJacobian_1(*scratch,*v_dual,*u_state,*z_all,tol);

      ajv_now->axpy(1.0,*scratch);   
    }

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

  virtual void applyAdjointJacobian_2(Vector<Real> &ajv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) override {

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

    Ptr<Vector<Real>> scratch;
    for(int i=0;i<pint_ajv.numOwnedSteps();i++) {
      // pull out all the time steps required
      Ptr<const Vector<Real>> v_dual  = getVector(pint_v,i,CURRENT);
      Ptr<const Vector<Real>> u_state = getStateVector(pint_u,i);
      Ptr<const Vector<Real>> z_all   = getVector(pint_z,i,ALL);

      Ptr<Vector<Real>> ajv_now     = getNonconstVector(pint_ajv,i,ALL);

      stepConstraint_->applyAdjointJacobian_2(*ajv_now,*v_dual,*u_state,*z_all,tol);
    }

    // no constraint for controls hanlding because that doesn't work yet!
  }

  virtual void applyInverseAdjointJacobian_1(Vector<Real> &iajv,
                                             const Vector<Real> &v,
                                             const Vector<Real> &u,
                                             const Vector<Real> &z,
                                             Real &tol) override final {
    ROL_TEST_FOR_EXCEPTION(true, std::logic_error,
      "The method applyInverseAdjointJacobian_1 is used but not implemented!\n");
  };

  /** \brief Apply the simulation-space derivative of the adjoint of the constraint
             simulation-space Jacobian at \f$(u,z)\f$ to the vector \f$w\f$ in the
             direction \f$v\f$, according to \f$v\mapsto c_{uu}(u,z)(v,\cdot)^*w\f$.

             @param[out]      ahwv is the result of applying the simulation-space derivative of the adjoint of the constraint simulation-space Jacobian at @b \f$(u,z)\f$ to the vector @b \f$w\f$ in direction @b \f$w\f$; a dual simulation-space vector
             @param[in]       w    is the direction vector; a dual constraint-space vector
             @param[in]       v    is a simulation-space vector
             @param[in]       u    is the constraint argument; a simulation-space vector
             @param[in]       z    is the constraint argument; an optimization-space vector
             @param[in,out]   tol  is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ahwv} = c_{uu}(u,z)(v,\cdot)^*w\f$, where
             \f$w \in \mathcal{C}^*\f$, \f$v \in \mathcal{U}\f$, and
             \f$\mathsf{ahwv} \in \mathcal{U}^*\f$.

             ---
  */
  virtual void applyAdjointHessian_11(Vector<Real> &ahwv,
                                      const Vector<Real> &w,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) override
  {
    PinTVector<Real>       & pint_ahwv = dynamic_cast<PinTVector<Real>&>(ahwv);
    const PinTVector<Real> & pint_w    = dynamic_cast<const PinTVector<Real>&>(w);
    const PinTVector<Real> & pint_v    = dynamic_cast<const PinTVector<Real>&>(v);
    const PinTVector<Real> & pint_u    = dynamic_cast<const PinTVector<Real>&>(u);
    const PinTVector<Real> & pint_z    = dynamic_cast<const PinTVector<Real>&>(z);
       // its possible we won't always want to cast to a PinT vector here
      
    TEUCHOS_ASSERT(pint_ahwv.numOwnedSteps()==pint_u.numOwnedSteps());
    TEUCHOS_ASSERT(pint_ahwv.numOwnedSteps()==pint_v.numOwnedSteps());

    // communicate neighbors, these are block calls
    pint_w.boundaryExchange();
    pint_v.boundaryExchange();
    pint_u.boundaryExchange();
    pint_z.boundaryExchange();

    // we need to make sure this has all zeros to begin with (this includes boundary exchange components)
    pint_ahwv.zeroAll();

    Ptr<Vector<Real>> scratch;
    for(int i=0;i<pint_ahwv.numOwnedSteps();i++) {
      // pull out all the time steps required
      Ptr<const Vector<Real>> w_dual = getVector(pint_w,i,CURRENT);
      Ptr<const Vector<Real>> v_dual = getStateVector(pint_v,i);
      Ptr<const Vector<Real>> u_state = getStateVector(pint_u,i);
      Ptr<const Vector<Real>> z_all = getVector(pint_z,i,ALL);

      Ptr<Vector<Real>> ahwv_now = getStateVector(pint_ahwv,i);

      // this will lead to problems if problem size changes in time
      if(scratch==ROL::nullPtr) {
        scratch = ahwv_now->clone();
      }

      stepConstraint_->applyAdjointHessian_11(*scratch,*w_dual,*v_dual,*u_state,*z_all,tol);

      ahwv_now->axpy(1.0,*scratch);   
    }

    // handle the constraint adjoint
    //   nothing to do here, as the constraint is linear, 
  }

  /** \brief Apply the optimization-space derivative of the adjoint of the constraint
             simulation-space Jacobian at \f$(u,z)\f$ to the vector \f$w\f$ in the
             direction \f$v\f$, according to \f$v\mapsto c_{uz}(u,z)(v,\cdot)^*w\f$.

             @param[out]      ahwv is the result of applying the optimization-space derivative of the adjoint of the constraint simulation-space Jacobian at @b \f$(u,z)\f$ to the vector @b \f$w\f$ in direction @b \f$w\f$; a dual optimization-space vector
             @param[in]       w    is the direction vector; a dual constraint-space vector
             @param[in]       v    is a simulation-space vector
             @param[in]       u    is the constraint argument; a simulation-space vector
             @param[in]       z    is the constraint argument; an optimization-space vector
             @param[in,out]   tol  is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ahwv} = c_{uz}(u,z)(v,\cdot)^*w\f$, where
             \f$w \in \mathcal{C}^*\f$, \f$v \in \mathcal{U}\f$, and
             \f$\mathsf{ahwv} \in \mathcal{Z}^*\f$.

             ---
  */
  virtual void applyAdjointHessian_12(Vector<Real> &ahwv,
                                      const Vector<Real> &w,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) override  
  {
    ahwv.zero();
  }

  /** \brief Apply the simulation-space derivative of the adjoint of the constraint 
             optimization-space Jacobian at \f$(u,z)\f$ to the vector \f$w\f$ in the 
             direction \f$v\f$, according to \f$v\mapsto c_{zu}(u,z)(v,\cdot)^*w\f$. 
 
             @param[out]      ahwv is the result of applying the simulation-space derivative of the adjoint of the constraint optimization-space Jacobian at @b \f$(u,z)\f$ to the vector @b \f$w\f$ in direction @b \f$w\f$; a dual simulation-space vector 
             @param[in]       w    is the direction vector; a dual constraint-space vector 
             @param[in]       v    is a optimization-space vector 
             @param[in]       u    is the constraint argument; a simulation-space vector 
             @param[in]       z    is the constraint argument; an optimization-space vector 
             @param[in,out]   tol  is a tolerance for inexact evaluations; currently unused 
 
             On return, \f$\mathsf{ahwv} = c_{zu}(u,z)(v,\cdot)^*w\f$, where 
             \f$w \in \mathcal{C}^*\f$, \f$v \in \mathcal{Z}\f$, and 
             \f$\mathsf{ahwv} \in \mathcal{U}^*\f$. 
 
             --- 
  */ 
  virtual void applyAdjointHessian_21(Vector<Real> &ahwv, 
                                      const Vector<Real> &w, 
                                      const Vector<Real> &v, 
                                      const Vector<Real> &u, 
                                      const Vector<Real> &z, 
                                      Real &tol) override { 
    ahwv.zero();
  }

  /** \brief Apply the optimization-space derivative of the adjoint of the constraint
             optimization-space Jacobian at \f$(u,z)\f$ to the vector \f$w\f$ in the
             direction \f$v\f$, according to \f$v\mapsto c_{zz}(u,z)(v,\cdot)^*w\f$.

             @param[out]      ahwv is the result of applying the optimization-space derivative of the adjoint of the constraint optimization-space Jacobian at @b \f$(u,z)\f$ to the vector @b \f$w\f$ in direction @b \f$w\f$; a dual optimization-space vector
             @param[in]       w    is the direction vector; a dual constraint-space vector
             @param[in]       v    is a optimization-space vector
             @param[in]       u    is the constraint argument; a simulation-space vector
             @param[in]       z    is the constraint argument; an optimization-space vector
             @param[in,out]   tol  is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ahwv} = c_{zz}(u,z)(v,\cdot)^*w\f$, where
             \f$w \in \mathcal{C}^*\f$, \f$v \in \mathcal{Z}\f$, and
             \f$\mathsf{ahwv} \in \mathcal{Z}^*\f$.

             ---
  */
  virtual void applyAdjointHessian_22(Vector<Real> &ahwv,
                                      const Vector<Real> &w,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) override {
    ahwv.zero();
  }

}; // class Constraint_SimOpt

} // namespace ROL

#endif
