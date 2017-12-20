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

template <class Real>
class Constraint_PinTSimOpt : public Constraint_SimOpt<Real> {
private:
  // user input
  MPI_Comm comm_;
 
  Ptr<ROL::Constraint_TimeSimOpt<Real>> stepConstraint_;
  int steps_;

  Ptr<ROL::Vector<Real>> stepState_;
  Ptr<ROL::Vector<Real>> stepControl_;

  std::vector<int> stateStencil_;
  std::vector<int> controlStencil_;

  // internal state members
  
  bool isInitialized_;
  
  int stepStart_;
  int stepEnd_;

public:

  //! Default constructor
  Constraint_PinTSimOpt()
    : Constraint_SimOpt<Real>()
    , isInitialized_(false)
    , stepStart_(-1)
    , stepEnd_(-1)
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
  Constraint_PinTSimOpt(MPI_Comm comm,
                        const Ptr<ROL::Constraint_TimeSimOpt<Real>> & stepConstraint,
                        int steps,
                        const Ptr<ROL::Vector<Real>> & stepState,
                        const std::vector<int> & stateStencil,
                        const Ptr<ROL::Vector<Real>> & stepControl,
                        const std::vector<int> & controlStencil)
    : Constraint_SimOpt<Real>()
    , isInitialized_(false)
    , stepStart_(-1)
    , stepEnd_(-1)
  { 
    initialize(comm,stepConstraint,steps,stepState,stateStencil,stepControl,controlStencil);
  }

  /** What step to start at on this processor (inclusive).
   */
  int stepStart() const { return stepStart_; }

  /** What step to stop before on this processor (non inclusive).
   */
  int stepEnd() const { return stepEnd_; }

  /** \brief Initialize this class, setting up parallel distribution.
   
   */
  void initialize(MPI_Comm comm,
                  const Ptr<ROL::Constraint_TimeSimOpt<Real>> & stepConstraint,
                  int steps,
                  const Ptr<ROL::Vector<Real>> & stepState,
                  const std::vector<int> & stateStencil,
                  const Ptr<ROL::Vector<Real>> & stepControl,
                  const std::vector<int> & controlStencil) 
  {
    // initialize user member variables
    comm_           = comm;
    steps_          = steps;
    stepConstraint_ = stepConstraint;
    stepState_      = stepState;
    stepControl_    = stepControl;
    stateStencil_   = stateStencil;
    controlStencil_ = controlStencil;

    // initialize the parallel variables
    
    int numRanks = -1;
    int myRank = -1;

    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // determine which steps are owned by this processor
    {
      int stepsPerRank = steps_ / numRanks;
      int remainder    = steps_ % numRanks; 

      stepStart_ = 0;

      if(myRank<remainder) {
        stepStart_ = myRank*(stepsPerRank+1);
        stepEnd_   = (myRank+1)*(stepsPerRank+1);
      }
      else if(myRank==remainder) {
        stepStart_ = myRank*(stepsPerRank+1);
        stepEnd_   = (myRank+1)*stepsPerRank + myRank;
      }
      else if(myRank>remainder) {
        stepStart_ = myRank*stepsPerRank + remainder;
        stepEnd_   = (myRank+1)*stepsPerRank + remainder;
      }
    }

    isInitialized_ = true;
  }

  /*
  ROL::Ptr<ROL::Vector<Real>> buildStateVector() const
  {
    assert(isInitialized_);

    return ROL::makePtr<PinTVector<Real>>(stepStart_,stepStop_,stateStencil_);
  }

  ROL::Ptr<ROL::Vector<Real>> buildControlVector() const
  {
    assert(isInitialized_);

    return ROL::makePtr<PinTVector<Real>>(stepStart_,stepStop_,controlStencil_);
  }
  */

  /** \brief Update constraint functions with respect to Opt variable.
                u is the state variable, 
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update_1( const Vector<Real> &z, bool flag = true, int iter = -1 ) {}

  /** \brief Update constraint functions with respect to Opt variable.
                z is the control variable, 
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update_2( const Vector<Real> &z, bool flag = true, int iter = -1 ) {}

  virtual void value(Vector<Real> &c,
                     const Vector<Real> &u,
                     const Vector<Real> &z,
                     Real &tol) override {
    TEUCHOS_ASSERT(false);
  } 

  virtual void solve(Vector<Real> &c,
                     Vector<Real> &u, 
                     const Vector<Real> &z,
                     Real &tol) override {
    TEUCHOS_ASSERT(false);
  }

  virtual void applyJacobian_1(Vector<Real> &jv,
                               const Vector<Real> &v,
                               const Vector<Real> &u,
                               const Vector<Real> &z,
                               Real &tol) override {
    TEUCHOS_ASSERT(false);
  }

  virtual void applyInverseJacobian_1(Vector<Real> &ijv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) override final {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
      "The method applyInverseJacobian_1 is used but not implemented!\n");
  }

  virtual void applyAdjointJacobian_1(Vector<Real> &ajv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      const Vector<Real> &dualv,
                                      Real &tol) override {
    TEUCHOS_ASSERT(false);
  }

  virtual void applyAdjointJacobian_2(Vector<Real> &ajv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      const Vector<Real> &dualv,
                                      Real &tol) override {
    TEUCHOS_ASSERT(false);
  }

  virtual void applyInverseAdjointJacobian_1(Vector<Real> &iajv,
                                             const Vector<Real> &v,
                                             const Vector<Real> &u,
                                             const Vector<Real> &z,
                                             Real &tol) override final {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
      "The method applyInverseAdjointJacobian_1 is used but not implemented!\n");
  };

}; // class Constraint_SimOpt

} // namespace ROL

#endif
