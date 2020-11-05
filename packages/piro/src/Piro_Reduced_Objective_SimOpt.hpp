// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#ifndef PIRO_REDUCED_OBJECTIVE_SIMOPT_H
#define PIRO_REDUCED_OBJECTIVE_SIMOPT_H

#include "ROL_Reduced_Objective_SimOpt.hpp"

namespace Piro {

template <class Real>
class Reduced_Objective_SimOpt : public ROL::Reduced_Objective_SimOpt<Real> {
   public:
    /** \brief Constructor.

        @param[in] obj          is a pointer to a SimOpt objective function.
        @param[in] con          is a pointer to a SimOpt equality constraint.
        @param[in] state        is a pointer to a state space vector, \f$\mathcal{U}\f$.
        @param[in] control      is a pointer to a optimization space vector, \f$\mathcal{Z}\f$.
        @param[in] adjoint      is a pointer to a dual constraint space vector, \f$\mathcal{C}^*\f$.
        @param[in] storage      is a flag whether or not to store computed states and adjoints.
        @param[in] useFDhessVec is a flag whether or not to use a finite-difference Hessian approximation.
    */
    Reduced_Objective_SimOpt(
        const ROL::Ptr<ROL::Objective_SimOpt<Real> > &obj,
        const ROL::Ptr<ROL::Constraint_SimOpt<Real> > &con,
        const ROL::Ptr<ROL::Vector<Real> > &state,
        const ROL::Ptr<ROL::Vector<Real> > &control,
        const ROL::Ptr<ROL::Vector<Real> > &adjoint,
        const bool storage = true,
        const bool useFDhessVec = false) : ROL::Reduced_Objective_SimOpt<Real>(obj, con, state, control, adjoint, storage, useFDhessVec) {}

    /** \brief Secondary, general constructor for use with dual optimization vector spaces where the user does not define the dual() method.

        @param[in] obj          is a pointer to a SimOpt objective function.
        @param[in] con          is a pointer to a SimOpt equality constraint.
        @param[in] state        is a pointer to a state space vector, \f$\mathcal{U}\f$.
        @param[in] control      is a pointer to a optimization space vector, \f$\mathcal{Z}\f$.
        @param[in] adjoint      is a pointer to a dual constraint space vector, \f$\mathcal{C}^*\f$.
        @param[in] dualstate    is a pointer to a dual state space vector, \f$\mathcal{U}^*\f$.
        @param[in] dualadjoint  is a pointer to a constraint space vector, \f$\mathcal{C}\f$.
        @param[in] storage      is a flag whether or not to store computed states and adjoints.
        @param[in] useFDhessVec is a flag whether or not to use a finite-difference Hessian approximation.
    */
    Reduced_Objective_SimOpt(
        const ROL::Ptr<ROL::Objective_SimOpt<Real> > &obj,
        const ROL::Ptr<ROL::Constraint_SimOpt<Real> > &con,
        const ROL::Ptr<ROL::Vector<Real> > &state,
        const ROL::Ptr<ROL::Vector<Real> > &control,
        const ROL::Ptr<ROL::Vector<Real> > &adjoint,
        const ROL::Ptr<ROL::Vector<Real> > &dualstate,
        const ROL::Ptr<ROL::Vector<Real> > &dualcontrol,
        const ROL::Ptr<ROL::Vector<Real> > &dualadjoint,
        const bool storage = true,
        const bool useFDhessVec = false) : ROL::Reduced_Objective_SimOpt<Real>(obj, con, state, control, adjoint, dualstate, dualcontrol, dualadjoint, storage, useFDhessVec) {}

    /** \brief Constructor.

        @param[in] obj          is a pointer to a SimOpt objective function.
        @param[in] con          is a pointer to a SimOpt equality constraint.
        @param[in] stateStore   is a pointer to a SimController object.
        @param[in] state        is a pointer to a state space vector, \f$\mathcal{U}\f$.
        @param[in] control      is a pointer to a optimization space vector, \f$\mathcal{Z}\f$.
        @param[in] adjoint      is a pointer to a dual constraint space vector, \f$\mathcal{C}^*\f$.
        @param[in] storage      is a flag whether or not to store computed states and adjoints.
        @param[in] useFDhessVec is a flag whether or not to use a finite-difference Hessian approximation.
    */
    Reduced_Objective_SimOpt(
        const ROL::Ptr<ROL::Objective_SimOpt<Real> > &obj,
        const ROL::Ptr<ROL::Constraint_SimOpt<Real> > &con,
        const ROL::Ptr<ROL::SimController<Real> > &stateStore,
        const ROL::Ptr<ROL::Vector<Real> > &state,
        const ROL::Ptr<ROL::Vector<Real> > &control,
        const ROL::Ptr<ROL::Vector<Real> > &adjoint,
        const bool storage = true,
        const bool useFDhessVec = false) : ROL::Reduced_Objective_SimOpt<Real>(obj, con, stateStore, state, control, adjoint, storage, useFDhessVec) {}

    /** \brief Secondary, general constructor for use with dual optimization vector spaces
        where the user does not define the dual() method.

        @param[in] obj          is a pointer to a SimOpt objective function.
        @param[in] con          is a pointer to a SimOpt equality constraint.
        @param[in] stateStore   is a pointer to a SimController object.
        @param[in] state        is a pointer to a state space vector, \f$\mathcal{U}\f$.
        @param[in] control      is a pointer to a optimization space vector, \f$\mathcal{Z}\f$.
        @param[in] adjoint      is a pointer to a dual constraint space vector, \f$\mathcal{C}^*\f$.
        @param[in] dualstate    is a pointer to a dual state space vector, \f$\mathcal{U}^*\f$.
        @param[in] dualadjoint  is a pointer to a constraint space vector, \f$\mathcal{C}\f$.
        @param[in] storage      is a flag whether or not to store computed states and adjoints.
        @param[in] useFDhessVec is a flag whether or not to use a finite-difference Hessian approximation.
    */
    Reduced_Objective_SimOpt(
        const ROL::Ptr<ROL::Objective_SimOpt<Real> > &obj,
        const ROL::Ptr<ROL::Constraint_SimOpt<Real> > &con,
        const ROL::Ptr<ROL::SimController<Real> > &stateStore,
        const ROL::Ptr<ROL::Vector<Real> > &state,
        const ROL::Ptr<ROL::Vector<Real> > &control,
        const ROL::Ptr<ROL::Vector<Real> > &adjoint,
        const ROL::Ptr<ROL::Vector<Real> > &dualstate,
        const ROL::Ptr<ROL::Vector<Real> > &dualcontrol,
        const ROL::Ptr<ROL::Vector<Real> > &dualadjoint,
        const bool storage = true,
        const bool useFDhessVec = false) : ROL::Reduced_Objective_SimOpt<Real>(obj, con, stateStore, state, control, adjoint, dualstate, dualcontrol, dualadjoint, storage, useFDhessVec) {}

    /** \brief Given \f$z\in\mathcal{Z}\f$, store a precomputed state
        \f$u=u(z)\in\mathcal{U}\f$ which solves \f$e(u,z) = 0\f$ and
        update the SimOpt objective function and equality constraint
        accordingly.
    */
    void set_precomputed_state(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z) {
        this->con_->update(u, z, this->updateFlag_, this->updateIter_);
        this->obj_->update(u, z, this->updateFlag_, this->updateIter_);
        if (this->storage_) {
            this->stateStore_->set(u, ROL::Objective<Real>::getParameter());
        }
    }
};  // class Reduced_Objective_SimOpt

}  // namespace Piro

#endif
