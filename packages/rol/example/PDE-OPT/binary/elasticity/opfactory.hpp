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

#ifndef OPFACTORY_BINARY_ELASTICITY_HPP
#define OPFACTORY_BINARY_ELASTICITY_HPP

#include "ROL_Bounds.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_OptimizationProblemFactory.hpp"

#include "../../TOOLS/pdeconstraint.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "../../TOOLS/meshreader.hpp"
#include "src/obj_compliance.hpp"
#include "src/con_volume.hpp"
#include "src/mesh_topo-opt.hpp"
#include "src/pde_elasticity.hpp"
#include "src/bilinearpdeconstraint.hpp"
#include "src/femdata.hpp"

template<class Real>
class ElasticityFactory : public ROL::OptimizationProblemFactory<Real> {
private:
  ROL::ParameterList                 pl_;
  ROL::Ptr<const Teuchos::Comm<int>> comm_;
  ROL::Ptr<std::ostream>             os_;

  ROL::Ptr<MeshManager<Real>>             mesh_;
  ROL::Ptr<PDE_Elasticity<Real>>          pde_;
  ROL::Ptr<Bilinear_PDE_Constraint<Real>> con_;
  ROL::Ptr<Assembler<Real>>               assembler_;
  ROL::Ptr<FEM_Data<Real>>                fem_;
  ROL::Ptr<ROL::Vector<Real>>             u_, z_, p_, r_;

  Real cmpScaling_;
  bool storage_;
  int M_, N_, T_, dim_, probDim_;

public:
  ElasticityFactory(ROL::ParameterList                 &pl,
              const ROL::Ptr<const Teuchos::Comm<int>> &comm,
              const ROL::Ptr<std::ostream>             &os)
    : pl_(pl), comm_(comm), os_(os) {
    std::vector<Real> ym = ROL::getArrayFromStringParameter<Real>(pl.sublist("Problem"), "Young's Modulus");
    probDim_    = pl.sublist("Problem").get("Problem Dimension", 2);
    cmpScaling_ = pl.sublist("Problem").get("Compliance Scaling", 1e0);
    storage_    = pl.sublist("Problem").get("Use Storage",true);
    M_          = pl.sublist("Problem").get("Number of Horizontal Cells",10);
    N_          = pl.sublist("Problem").get("Number of Vertical Cells",20);
    T_          = ym.size();
    dim_        = M_*N_*T_;

    if (probDim_ == 2) {
      mesh_ = ROL::makePtr<MeshManager_TopoOpt<Real>>(pl_);
    } else if (probDim_ == 3) {
      mesh_ = ROL::makePtr<MeshReader<Real>>(pl_);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/binary/elasticity/example_01.cpp: Problem dim is not 2 or 3!");
    }
    pde_       = ROL::makePtr<PDE_Elasticity<Real>>(pl_);
    fem_       = ROL::makePtr<FEM_Data<Real>>(pde_,mesh_,comm_,pl_,*os_);
    assembler_ = fem_->getAssembler();
  }

  void check(void) {
    if (con_ == ROL::nullPtr) {
      update();
    }
    ROL::Ptr<ROL::Vector<Real>> ru = u_->clone(); ru->randomize();
    ROL::Ptr<ROL::Vector<Real>> rz = z_->clone(); rz->randomize();
    ROL::Ptr<ROL::Vector<Real>> rp = p_->clone(); rp->randomize();
    ROL::Ptr<ROL::Vector<Real>> rv = u_->clone(); rv->randomize();
    ROL::Ptr<ROL::Vector<Real>> rr = r_->clone(); rr->randomize();

    *os_ << std::endl;
    rr = r_->clone(); rr->randomize();
    con_->checkSolve(*ru,*rz,*rr,true,*os_);

    *os_ << std::endl;
    rv = u_->clone(); rv->randomize();
    rr = r_->clone(); rr->randomize();
    con_->checkApplyJacobian_1(*ru,*rz,*rv,*rr,true,*os_);

    *os_ << std::endl;
    rv = u_->clone(); rv->randomize();
    rr = r_->clone(); rr->randomize();
    con_->checkInverseJacobian_1(*rr,*rv,*ru,*rz,true,*os_);

    *os_ << std::endl;
    rv = u_->clone(); rv->randomize();
    rr = r_->clone(); rr->randomize();
    con_->checkInverseAdjointJacobian_1(*rr,*rv,*ru,*rz,true,*os_);

    *os_ << std::endl;
    rv = z_->clone(); rv->randomize();
    rr = r_->clone(); rr->randomize();
    con_->checkApplyJacobian_2(*ru,*rz,*rv,*rr,true,*os_);

    *os_ << std::endl;
    rv = u_->clone(); rv->randomize();
    rr = r_->clone(); rr->randomize();
    con_->checkApplyAdjointHessian_11(*ru,*rz,*rp,*rv,*rr,true,*os_);
    
    *os_ << std::endl;
    rv = z_->clone(); rv->randomize();
    rr = u_->clone(); rr->randomize();
    con_->checkApplyAdjointHessian_21(*ru,*rz,*rp,*rv,*rr,true,*os_);

    *os_ << std::endl;
    rv = u_->clone(); rv->randomize();
    rr = z_->clone(); rr->randomize();
    con_->checkApplyAdjointHessian_12(*ru,*rz,*rp,*rv,*rr,true,*os_);

    *os_ << std::endl;
    rv = z_->clone(); rv->randomize();
    rr = z_->clone(); rr->randomize();
    con_->checkApplyAdjointHessian_22(*ru,*rz,*rp,*rv,*rr,true,*os_);
  }

  void update(void) {
    //con_ = ROL::makePtr<PDE_Constraint<Real>>(pde_,mesh_,comm_,pl_,*os_);
    con_ = ROL::makePtr<Bilinear_PDE_Constraint<Real>>(fem_,pl_,*os_);
    con_->setSolveParameters(pl_);

    ROL::Ptr<Tpetra::MultiVector<>> u_ptr, p_ptr, r_ptr;
    u_ptr = assembler_->createStateVector();
    p_ptr = assembler_->createStateVector();
    r_ptr = assembler_->createResidualVector();
    u_ = ROL::makePtr<PDE_PrimalSimVector<Real>>(u_ptr,pde_,assembler_,pl_);
    p_ = ROL::makePtr<PDE_PrimalSimVector<Real>>(p_ptr,pde_,assembler_,pl_);
    r_ = ROL::makePtr<PDE_DualSimVector<Real>>(r_ptr,pde_,assembler_,pl_);
    z_ = ROL::makePtr<ROL::StdVector<Real>>(dim_);
  }

  ROL::Ptr<ROL::Objective<Real>> buildObjective(void) {
    ROL::Ptr<QoI<Real>> qoi
      = ROL::makePtr<QoI_Compliance_TopoOpt<Real>>(pde_->getFE(),
                                                   pde_->getLoad(),
                                                   pde_->getBdryFE(),
                                                   pde_->getBdryCellLocIds(),
                                                   pde_->getTraction(),
                                                   pde_->getFieldHelper(),
                                                   cmpScaling_);
    ROL::Ptr<ROL::Objective_SimOpt<Real>> obj
      = ROL::makePtr<PDE_Objective<Real>>(qoi,assembler_);
    return ROL::makePtr<ROL::Reduced_Objective_SimOpt<Real>>(obj, con_, u_, z_, p_, storage_, false);
  }

  ROL::Ptr<ROL::Vector<Real>> buildSolutionVector(void) {
    return z_;
  }

  ROL::Ptr<ROL::BoundConstraint<Real>> buildBoundConstraint(void) {
    ROL::Ptr<ROL::Vector<Real>> zlop, zhip;
    zlop = ROL::makePtr<ROL::StdVector<Real>>(dim_,static_cast<Real>(0));
    zhip = ROL::makePtr<ROL::StdVector<Real>>(dim_,static_cast<Real>(1));
    return ROL::makePtr<ROL::Bounds<Real>>(zlop,zhip);
  }

  ROL::Ptr<ROL::Constraint<Real>> buildEqualityConstraint(void) {
    bool useIneq = pl_.sublist("Problem").get("Use Inequality", false);
    if (!useIneq) {
      return ROL::makePtr<Selection_TopoOpt<Real>>(pl_);
    }
    return ROL::nullPtr;
  }

  ROL::Ptr<ROL::Vector<Real>> buildEqualityMultiplier(void) {
    bool useIneq = pl_.sublist("Problem").get("Use Inequality", false);
    if (!useIneq) {
      return ROL::makePtr<ROL::StdVector<Real>>(M_*N_,static_cast<Real>(0));
    }
    return ROL::nullPtr;
  }

  ROL::Ptr<ROL::Constraint<Real>> buildInequalityConstraint(void) {
    //return ROL::nullPtr;
    bool useIneq = pl_.sublist("Problem").get("Use Inequality", false);
    if (!useIneq) {
      return ROL::makePtr<Volume_TopoOpt<Real>>(pl_);
    }
    std::vector<ROL::Ptr<ROL::Constraint<Real>>> icon(2);
    icon[0] = ROL::makePtr<Volume_TopoOpt<Real>>(pl_);
    icon[1] = ROL::makePtr<Selection_TopoOpt<Real>>(pl_);
    return ROL::makePtr<ROL::Constraint_Partitioned<Real>>(icon,true);
  }

  ROL::Ptr<ROL::Vector<Real>> buildInequalityMultiplier(void) {
    //return ROL::nullPtr;
    bool useIneq = pl_.sublist("Problem").get("Use Inequality", false);
    if (!useIneq) {
      return ROL::makePtr<ROL::StdVector<Real>>(dim_,static_cast<Real>(0));
    }
    std::vector<ROL::Ptr<ROL::Vector<Real>>> imul(2);
    imul[0] = ROL::makePtr<ROL::StdVector<Real>>(dim_,static_cast<Real>(0));
    imul[1] = ROL::makePtr<ROL::StdVector<Real>>(M_*N_,static_cast<Real>(0));
    return ROL::makePtr<ROL::PartitionedVector<Real>>(imul);
  }

  ROL::Ptr<ROL::BoundConstraint<Real>> buildInequalityBoundConstraint(void) {
    //return ROL::nullPtr;
    bool useIneq = pl_.sublist("Problem").get("Use Inequality", false);
    if (!useIneq) {
      ROL::Ptr<ROL::Vector<Real>> iup;
      iup = ROL::makePtr<ROL::StdVector<Real>>(dim_,static_cast<Real>(0));
      return ROL::makePtr<ROL::Bounds<Real>>(*iup,false);
    }
    std::vector<ROL::Ptr<ROL::Vector<Real>>> ivec(2);
    std::vector<ROL::Ptr<ROL::BoundConstraint<Real>>> ibnd(2);
    ivec[0] = ROL::makePtr<ROL::StdVector<Real>>(dim_,static_cast<Real>(0));
    ivec[1] = ROL::makePtr<ROL::StdVector<Real>>(M_*N_,static_cast<Real>(0));
    ibnd[0] = ROL::makePtr<ROL::Bounds<Real>>(*ivec[0],false);
    ibnd[1] = ROL::makePtr<ROL::Bounds<Real>>(*ivec[1],false);
    return ROL::makePtr<ROL::BoundConstraint_Partitioned<Real>>(ibnd,ivec);
  }

  void print(void) const {
    Real tol(1e-8);
    con_->printMeshData(*os_);
    //con_->solve(*r_,*u_,*z_,tol);
    //con_->outputTpetraVector(ROL::staticPtrCast<ROL::
    std::vector<Real> &z = *ROL::staticPtrCast<ROL::StdVector<Real>>(z_)->getVector();
    std::ofstream file;
    file.open("design.txt");
    for (int i = 0; i < M_; ++i) {
      for (int j = 0; j < N_; ++j) {
        for (int k = 0; k < T_; ++k) {
          file << i << "  " << j << "  " << k << "  " << z[i+M_*(j+N_*k)] << std::endl;
        }
      }
    }
    file.close();
  }
};

#endif
