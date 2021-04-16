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

#ifndef OPFACTORY_TOPOOPT_MULTIMAT_HPP
#define OPFACTORY_TOPOOPT_MULTIMAT_HPP

#include "ROL_Bounds.hpp"
#include "ROL_Problem.hpp"

#include "../../TOOLS/pdevector.hpp"
#include "../../TOOLS/meshreader.hpp"
#include "src/mesh_topo-opt.hpp"
#include "src/pde_elasticity.hpp"
#include "src/pde_filter.hpp"
#include "src/pde_selection.hpp"
#include "src/obj_volume.hpp"
#include "src/compliance_robj.hpp"
#include "src/filtered_compliance_robj.hpp"
#include "src/con_selection.hpp"
#include "src/volume_con.hpp"

template<typename Real>
class ElasticityFactory {
private:
  ROL::ParameterList                 pl_;
  ROL::Ptr<const Teuchos::Comm<int>> comm_;
  ROL::Ptr<std::ostream>             os_;

  ROL::Ptr<MeshManager<Real>>             mesh_;
  ROL::Ptr<PDE_MultiMat_Elasticity<Real>> pde_;
  ROL::Ptr<MultiMat_PDE_Filter<Real>>     filter_;
  bool useFilter_;

  ROL::Ptr<PDE_MultiMat_Selection<Real>> spde_;
  ROL::Ptr<ROL::Constraint<Real>>        scon_;
  ROL::Ptr<Assembler<Real>>              sassembler_;
  ROL::Ptr<ROL::Vector<Real>>            smul_;
  ROL::Ptr<ROL::Vector<Real>>            slo_, sup_;
  ROL::Ptr<ROL::BoundConstraint<Real>>   sbnd_;

  ROL::Ptr<QoI_MultiMat_Weight<Real>>  vqoi_;
  ROL::Ptr<ROL::Constraint<Real>>      vcon_;
  ROL::Ptr<ROL::Vector<Real>>          vmul_;
  ROL::Ptr<ROL::Vector<Real>>          vlo_, vup_;
  ROL::Ptr<ROL::BoundConstraint<Real>> vbnd_;

  bool useIneq_;

  ROL::Ptr<ROL::Vector<Real>> z_, zlo_, zup_;
  ROL::Ptr<ROL::BoundConstraint<Real>> bnd_;

  ROL::Ptr<ROL::Objective<Real>> obj_;

public:
  ElasticityFactory(ROL::ParameterList                 &pl,
              const ROL::Ptr<const Teuchos::Comm<int>> &comm,
              const ROL::Ptr<std::ostream>             &os)
    : pl_(pl), comm_(comm), os_(os) {
    // Elasticity PDE
    int probDim = pl_.sublist("Problem").get("Problem Dimension",2);
    if (probDim == 2) {
      mesh_ = ROL::makePtr<MeshManager_MultiMat<Real>>(pl_);
    } else if (probDim == 3) {
      mesh_ = ROL::makePtr<MeshReader<Real>>(pl_);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/topo-opt/multimat/example_01.cpp: Problem dim is not 2 or 3!");
    }
    pde_ = ROL::makePtr<PDE_MultiMat_Elasticity<Real>>(pl_);

    // Selection constraint
    spde_ = ROL::makePtr<PDE_MultiMat_Selection<Real>>(pl_);
    scon_ = ROL::makePtr<MultiMat_Selection_Constraint<Real>>(spde_, mesh_, comm_, pl_, *os_);
    sassembler_ = ROL::dynamicPtrCast<MultiMat_Selection_Constraint<Real>>(scon_)->getAssembler();
    smul_ = ROL::makePtr<PDE_DualSimVector<Real>>(sassembler_->createStateVector(),spde_,sassembler_,pl_);
    slo_  = smul_->dual().clone(); slo_->setScalar(static_cast<Real>(-1));
    sup_  = smul_->dual().clone(); sup_->setScalar(static_cast<Real>(0));
    sbnd_ = ROL::makePtr<ROL::Bounds<Real>>(slo_,sup_);

    // Create density vector
    z_ = ROL::makePtr<PDE_PrimalOptVector<Real>>(sassembler_->createControlVector(),spde_,sassembler_,pl_);
    //z_->setScalar(static_cast<Real>(1.0/3.0));
    //z_->zero();
    z_->setScalar(static_cast<Real>(1));

    // Volume constraint
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    ROL::Ptr<ROL::Vector<Real>> z0 = z_->clone(); z0->zero();
    vqoi_ = ROL::makePtr<QoI_MultiMat_Weight<Real>>(spde_->getDensityFE(), spde_->getDensityFieldInfo(), pl_);
    vcon_ = ROL::makePtr<MultiMat_Volume_Constraint<Real>>(vqoi_,sassembler_,z_);
    vmul_ = ROL::makePtr<ROL::SingletonVector<Real>>();
    vlo_  = vmul_->dual().clone(); vcon_->value(*vlo_,*z0,tol);
    vup_  = vmul_->dual().clone(); vup_->setScalar(static_cast<Real>(0));
    vbnd_ = ROL::makePtr<ROL::Bounds<Real>>(vlo_,vup_);

    useIneq_ = pl_.sublist("Problem").get("Use Inequality", false);

    // Density bound vectors
    zlo_ = z_->clone(); zlo_->setScalar(static_cast<Real>(0));
    zup_ = z_->clone(); zup_->setScalar(static_cast<Real>(1));
    bnd_ = ROL::makePtr<ROL::Bounds<Real>>(zlo_,zup_);

    // Objective function
    useFilter_ = pl_.sublist("Problem").get("Use Filter", false);
    if (useFilter_) { 
      filter_ = ROL::makePtr<MultiMat_PDE_Filter<Real>>(pl_);
      pde_->setDensityFields(filter_->getFields());
      obj_ = ROL::makePtr<MultiMat_Filtered_Compliance_Objective<Real>>(pde_, filter_, mesh_, comm_, pl_, *os_);
    }
    else {
      obj_ = ROL::makePtr<MultiMat_Compliance_Objective<Real>>(pde_, mesh_, comm_, pl_, *os_);
    }
  }

  void check(void) {
    build()->check(true,*os_);
  }

  ROL::Ptr<ROL::Problem<Real>> build(void) {
    ROL::Ptr<ROL::Problem<Real>> prob
      = ROL::makePtr<ROL::Problem<Real>>(obj_,z_);
    prob->addBoundConstraint(bnd_);
    prob->addLinearConstraint("Weight",vcon_,vmul_,vbnd_);
    if (!useIneq_) {
      prob->addLinearConstraint("Selection",scon_,smul_);
    }
    else {
      prob->addLinearConstraint("Selection",scon_,smul_,sbnd_);
    }
    prob->setProjectionAlgorithm(pl_);
    prob->finalize(false,true,*os_);
    return prob;
  }

  void print(const ROL::Vector<Real> &z) {
    if (useFilter_)
      ROL::dynamicPtrCast<MultiMat_Filtered_Compliance_Objective<Real>>(obj_)->printToFile(
        *dynamic_cast<const ROL::PartitionedVector<Real>&>(z).get(0),*os_);
    else
      ROL::dynamicPtrCast<MultiMat_Compliance_Objective<Real>>(obj_)->printToFile(
        *dynamic_cast<const ROL::PartitionedVector<Real>&>(z).get(0),*os_);
  }
};

#endif
