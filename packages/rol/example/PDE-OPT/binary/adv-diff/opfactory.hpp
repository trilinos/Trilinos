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

#ifndef OPFACTORY_BINARY_ADVDIFF_HPP
#define OPFACTORY_BINARY_ADVDIFF_HPP

#include "ROL_Bounds.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_OptimizationProblemFactory.hpp"

#include "../../TOOLS/linearpdeconstraint.hpp"
#include "../../TOOLS/pdeconstraint.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "../../TOOLS/integralconstraint.hpp"
#include "pde_adv_diff.hpp"
#include "qoi_adv_diff.hpp"
#include "mesh_adv_diff.hpp"
#include "extractQP.hpp"

template<class Real>
class BinaryAdvDiffFactory : public ROL::OptimizationProblemFactory<Real> {
private:
  int dim_;

  ROL::ParameterList pl_;
  ROL::Ptr<const Teuchos::Comm<int>> comm_;
  ROL::Ptr<std::ostream> os_;

  ROL::Ptr<MeshManager_adv_diff<Real>>  mesh_;
  ROL::Ptr<PDE_adv_diff<Real>>          pde_;
  ROL::Ptr<Linear_PDE_Constraint<Real>> con_;
  ROL::Ptr<Assembler<Real>>             assembler_;
  ROL::Ptr<ROL::Vector<Real>>           u_, z_, p_;

public:
  BinaryAdvDiffFactory(ROL::ParameterList                 &pl,
                 const ROL::Ptr<const Teuchos::Comm<int>> &comm,
                 const ROL::Ptr<std::ostream>             &os)
    : pl_(pl), comm_(comm), os_(os) {
    int nx = pl.sublist("Problem").get("Number Controls - X", 3);
    int ny = pl.sublist("Problem").get("Number Controls - Y", 3);
    dim_ = nx*ny;
  }

  void update(void) {
    mesh_      = ROL::makePtr<MeshManager_adv_diff<Real>>(pl_);
    pde_       = ROL::makePtr<PDE_adv_diff<Real>>(pl_);
    con_       = ROL::makePtr<Linear_PDE_Constraint<Real>>(pde_,mesh_,comm_,pl_,*os_);
    assembler_ = con_->getAssembler();

    ROL::Ptr<Tpetra::MultiVector<>> u_ptr, p_ptr, r_ptr;
    u_ptr = assembler_->createStateVector();
    p_ptr = assembler_->createStateVector();
    r_ptr = assembler_->createResidualVector();
    ROL::Ptr<ROL::Vector<Real>> up, pp, rp, zp;
    u_ = ROL::makePtr<PDE_PrimalSimVector<Real>>(u_ptr,pde_,assembler_,pl_);
    p_ = ROL::makePtr<PDE_PrimalSimVector<Real>>(p_ptr,pde_,assembler_,pl_);
    z_ = ROL::makePtr<PDE_OptVector<Real>>(ROL::makePtr<ROL::StdVector<Real>>(dim_));
  }

  ROL::Ptr<ROL::Objective<Real>> buildObjective(void) {
    std::vector<ROL::Ptr<QoI<Real>>> qoi_vec(1,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_State_Cost_adv_diff<Real>>(pde_->getFE());
    Real stateCost = pl_.sublist("Problem").get("State Cost",4.0);
    std::vector<Real> wts = {stateCost};
    ROL::Ptr<ROL::Objective_SimOpt<Real>> obj
      = ROL::makePtr<PDE_Objective<Real>>(qoi_vec,wts,assembler_);
    bool storage = pl_.sublist("Problem").get("Use Storage",true);
    return ROL::makePtr<ROL::Reduced_Objective_SimOpt<Real>>(obj, con_, u_, z_, p_, storage, false);
  }

  ROL::Ptr<ROL::Vector<Real>> buildSolutionVector(void) {
    return z_;
  }

  ROL::Ptr<ROL::BoundConstraint<Real>> buildBoundConstraint(void) {
    ROL::Ptr<ROL::Vector<Real>> zlop
      = ROL::makePtr<PDE_OptVector<Real>>(ROL::makePtr<ROL::StdVector<Real>>(dim_,0.0));
    ROL::Ptr<ROL::Vector<Real>> zhip
      = ROL::makePtr<PDE_OptVector<Real>>(ROL::makePtr<ROL::StdVector<Real>>(dim_,1.0));
    return ROL::makePtr<ROL::Bounds<Real>>(zlop,zhip);
  }

  ROL::Ptr<ROL::Constraint<Real>> buildEqualityConstraint(void) {
    bool useIneq  = pl_.sublist("Problem").get("Use Inequality", false);
    if (!useIneq) {
      Real budget = pl_.sublist("Problem").get("Control Cost", 4.0);
      ROL::Ptr<QoI<Real>> qoi
        = ROL::makePtr<QoI_Control_Cost_adv_diff<Real>>(budget);
      return ROL::makePtr<IntegralOptConstraint<Real>>(qoi,assembler_);
    }
    return ROL::nullPtr;
  }

  ROL::Ptr<ROL::Vector<Real>> buildEqualityMultiplier(void) {
    bool useIneq  = pl_.sublist("Problem").get("Use Inequality", false);
    if (!useIneq) {
      return ROL::makePtr<ROL::StdVector<Real>>(1,0.0);
    }
    return ROL::nullPtr;
  }

  ROL::Ptr<ROL::Constraint<Real>> buildInequalityConstraint(void) {
    bool useIneq  = pl_.sublist("Problem").get("Use Inequality", false);
    if (useIneq) {
      Real budget   = static_cast<Real>(0);
      ROL::Ptr<QoI<Real>> qoi
        = ROL::makePtr<QoI_Control_Cost_adv_diff<Real>>(budget);
      return ROL::makePtr<IntegralOptConstraint<Real>>(qoi,assembler_);
    }
    return ROL::nullPtr;
  }

  ROL::Ptr<ROL::Vector<Real>> buildInequalityMultiplier(void) {
    bool useIneq  = pl_.sublist("Problem").get("Use Inequality", false);
    if (useIneq) {
      return ROL::makePtr<ROL::StdVector<Real>>(1,0.0);
    }
    return ROL::nullPtr;
  }

  ROL::Ptr<ROL::BoundConstraint<Real>> buildInequalityBoundConstraint(void) {
    bool useIneq  = pl_.sublist("Problem").get("Use Inequality", false);
    if (useIneq) {
      Real budget = pl_.sublist("Problem").get("Control Cost", 4.0);
      ROL::Ptr<ROL::Vector<Real>> klop, khip;
      klop = ROL::makePtr<ROL::StdVector<Real>>(1,static_cast<Real>(0));
      khip = ROL::makePtr<ROL::StdVector<Real>>(1,budget);
      return ROL::makePtr<ROL::Bounds<Real>>(klop,khip);
    }
    return ROL::nullPtr;
  }
};

template<class Real>
class BinaryAdvDiffQPFactory : public ROL::OptimizationProblemFactory<Real> {
private:
  ROL::ParameterList pl_;
  ROL::Ptr<const Teuchos::Comm<int>> comm_;
  ROL::Ptr<std::ostream> os_;

  ROL::Ptr<BinaryAdvDiffFactory<Real>> factory_;
  ROL::Ptr<extractQP<Real>> extract_;
  ROL::Ptr<ROL::OptimizationProblem<Real>> problem_;

public:
  BinaryAdvDiffQPFactory(ROL::ParameterList                 &pl,
                   const ROL::Ptr<const Teuchos::Comm<int>> &comm,
                   const ROL::Ptr<std::ostream>             &os)
    : pl_(pl), comm_(comm), os_(os) {}

  void update(void) {
    factory_ = ROL::makePtr<BinaryAdvDiffFactory<Real>>(pl_,comm_,os_);
    factory_->update();
    bool useIneq  = pl_.sublist("Problem").get("Use Inequality", false);
    if (useIneq) {
      extract_ = ROL::makePtr<extractQP<Real>>(factory_->buildObjective(),
                                               factory_->buildSolutionVector(),
                                               factory_->buildBoundConstraint(),
                                               factory_->buildInequalityConstraint(),
                                               factory_->buildInequalityMultiplier(),
                                               factory_->buildInequalityBoundConstraint());
    }
    else {
      extract_ = ROL::makePtr<extractQP<Real>>(factory_->buildObjective(),
                                               factory_->buildSolutionVector(),
                                               factory_->buildBoundConstraint(),
                                               factory_->buildEqualityConstraint(),
                                               factory_->buildEqualityMultiplier());
    }
    problem_ = (*extract_)();
  }

  ROL::Ptr<ROL::Objective<Real>> buildObjective(void) {
    return problem_->getObjective();
  }

  ROL::Ptr<ROL::Vector<Real>> buildSolutionVector(void) {
    return problem_->getSolutionVector();
  }

  ROL::Ptr<ROL::BoundConstraint<Real>> buildBoundConstraint(void) {
    return problem_->getBoundConstraint();
  }

  ROL::Ptr<ROL::Constraint<Real>> buildEqualityConstraint(void) {
    return problem_->getConstraint();
  }

  ROL::Ptr<ROL::Vector<Real>> buildEqualityMultiplier(void) {
    return problem_->getMultiplierVector();
  }
};
  

#endif
