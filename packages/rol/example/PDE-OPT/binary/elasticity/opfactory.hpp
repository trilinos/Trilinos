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
#include "ROL_PEBBL_IntegerProblemFactory.hpp"

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

template<class Real>
class ElasticityFactory : public ROL::PEBBL::IntegerProblemFactory<Real> {
private:
  ROL::ParameterList                 pl_;
  ROL::Ptr<const Teuchos::Comm<int>> comm_;
  ROL::Ptr<std::ostream>             os_;

  ROL::Ptr<MeshManager<Real>>         mesh_;
  ROL::Ptr<PDE_Elasticity<Real>>      pde_;
  ROL::Ptr<MultiMat_PDE_Filter<Real>> pde_filter_;
  ROL::Ptr<Filter<Real>>              filter_;
  ROL::Ptr<Assembler<Real>>           assembler_;
  bool useFilter_;

  ROL::Ptr<PDE_Selection<Real>>        spde_;
  ROL::Ptr<Selection_Constraint<Real>> scon_;
  ROL::Ptr<Assembler<Real>>            sassembler_;
  ROL::Ptr<ROL::Vector<Real>>          smul_;
  ROL::Ptr<ROL::Vector<Real>>          slo_, sup_;
  ROL::Ptr<ROL::BoundConstraint<Real>> sbnd_;

  ROL::Ptr<QoI_Weight<Real>>           vqoi_;
  ROL::Ptr<Volume_Constraint<Real>>    vcon_;
  ROL::Ptr<ROL::Vector<Real>>          vmul_;
  ROL::Ptr<ROL::Vector<Real>>          vlo_, vup_;
  ROL::Ptr<ROL::BoundConstraint<Real>> vbnd_;

  bool useIneq_;
  bool useLinCon_;
  size_t numMat_;

  ROL::Ptr<ROL::Vector<Real>> z_, zlo_, zup_;
  ROL::Ptr<ROL::Bounds<Real>> bnd_;

  bool serial_;
  ROL::Ptr<ROL::Objective<Real>> obj_;

public:
  ElasticityFactory(ROL::ParameterList                 &pl,
              const ROL::Ptr<const Teuchos::Comm<int>> &comm,
              const ROL::Ptr<std::ostream>             &os)
    : pl_(pl), comm_(comm), os_(os) {
    // Elasticity PDE
    int probDim = pl_.sublist("Problem").get("Problem Dimension",2);
    TEUCHOS_TEST_FOR_EXCEPTION(probDim<2||probDim>3, std::invalid_argument,
      ">>> PDE-OPT/binary/elasticity/example_01.cpp: Problem dim is not 2 or 3!");
    if (probDim == 2)      mesh_ = ROL::makePtr<MeshManager_TopoOpt<Real>>(pl_);
    else if (probDim == 3) mesh_ = ROL::makePtr<MeshReader<Real>>(pl_);
    pde_ = ROL::makePtr<PDE_Elasticity<Real>>(pl_);

    useFilter_ = pl.sublist("Problem").get("Use Filter",true);
    if (useFilter_) {
      pde_filter_ = ROL::makePtr<MultiMat_PDE_Filter<Real>>(pl_);
      pde_->setDensityFields(pde_filter_->getFields());
      filter_ = ROL::makePtr<Filter<Real>>(pde_filter_, mesh_, comm_, pl_, *os_);
    }
    assembler_ = ROL::makePtr<Assembler<Real>>(pde_->getFields(),
                                               pde_->getFields2(),
                                               mesh_, comm_, pl_, *os_);
    assembler_->setCellNodes(*pde_);

    // Selection constraint
    spde_ = ROL::makePtr<PDE_Selection<Real>>(pl_);
    scon_ = ROL::makePtr<Selection_Constraint<Real>>(spde_, mesh_, comm_, pl_, *os_);
    sassembler_ = scon_->getAssembler();
    smul_ = ROL::makePtr<PDE_DualSimVector<Real>>(sassembler_->createStateVector(),spde_,sassembler_,pl_);
    slo_  = smul_->dual().clone(); slo_->setScalar(static_cast<Real>(-1));
    sup_  = smul_->dual().clone(); sup_->setScalar(static_cast<Real>(0));
    sbnd_ = ROL::makePtr<ROL::Bounds<Real>>(slo_,sup_);
    useIneq_ = pl_.sublist("Problem").get("Use Inequality", false);
    useLinCon_ = pl_.sublist("Problem").get("Use Linear Constraints", true);

    // Create density vector
    z_ = ROL::makePtr<PDE_PrimalOptVector<Real>>(sassembler_->createControlVector(),spde_,sassembler_,pl_);
    z_->setScalar(static_cast<Real>(1));
    bool inDens = pl_.sublist("Problem").get("Input Density",false);
    if (inDens) {
      std::string inDensName = pl_.sublist("Problem").get("Input Density Name","density.txt");
      ROL::Ptr<Tpetra::MultiVector<>> zdata = ROL::dynamicPtrCast<ROL::TpetraMultiVector<Real>>(z_)->getVector();
      sassembler_->inputTpetraVector(zdata, inDensName);
//      ROL::Ptr<ROL::Vector<Real>> zr = z_->clone(); zr->randomize(-0.1,0.1);
//      z_->plus(*zr);
    }

    // Volume constraint
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    ROL::Ptr<ROL::Vector<Real>> z0 = z_->clone(); z0->zero();
    vqoi_ = ROL::makePtr<QoI_Weight<Real>>(spde_->getDensityFE(), spde_->getDensityFieldInfo(), pl_);
    vcon_ = ROL::makePtr<Volume_Constraint<Real>>(vqoi_,sassembler_,z_);
    vmul_ = ROL::makePtr<ROL::SingletonVector<Real>>();
    vlo_  = vmul_->dual().clone(); vcon_->value(*vlo_,*z0,tol);
    vup_  = vmul_->dual().clone(); vup_->setScalar(static_cast<Real>(0));
    vbnd_ = ROL::makePtr<ROL::Bounds<Real>>(vlo_,vup_);

    // Density bound vectors
    zlo_ = z_->clone(); zlo_->setScalar(static_cast<Real>(0));
    zup_ = z_->clone(); zup_->setScalar(static_cast<Real>(1));
    bnd_ = ROL::makePtr<ROL::Bounds<Real>>(zlo_,zup_);

    // Get number of materials
    std::vector<Real> ym = ROL::getArrayFromStringParameter<Real>(pl_.sublist("Problem"), "Young's Modulus");
    numMat_ = ym.size();

    // Build objective if running in serial
    serial_ = pl_.sublist("Problem").get("Serial",true);
    if (serial_) {
      if (useFilter_)
        obj_ = ROL::makePtr<Filtered_Compliance_Objective<Real>>(filter_, pde_, assembler_, pl_);
      else
        obj_ = ROL::makePtr<Compliance_Objective<Real>>(pde_, assembler_, pl_);
    }
  }

  ROL::Ptr<ROL::PEBBL::IntegerProblem<Real>> build(void) override {
    ROL::Ptr<ROL::PEBBL::IntegerProblem<Real>>
      problem = ROL::makePtr<ROL::PEBBL::IntegerProblem<Real>>(buildObjective(),buildSolutionVector());
    problem->addBoundConstraint(bnd_);
    ROL::Ptr<ROL::Vector<Real>> smul = smul_->clone();
    ROL::Ptr<ROL::Vector<Real>> vmul = vmul_->clone();
    if (useLinCon_) {
      if (numMat_ > size_t(1)) {
        if (useIneq_) problem->addLinearConstraint("Selection",scon_,smul,sbnd_);
        else          problem->addLinearConstraint("Selection",scon_,smul);
      }
      problem->addLinearConstraint("Weight",vcon_,vmul,vbnd_);
    }
    else {
      if (numMat_ > size_t(1)) {
        if (useIneq_) problem->addConstraint("Selection",scon_,smul,sbnd_);
        else          problem->addConstraint("Selection",scon_,smul);
      }
      problem->addConstraint("Weight",vcon_,vmul,vbnd_);
    }
    problem->setProjectionAlgorithm(pl_);
    return problem;
  }

  void check(void) {
    build()->check(true,*os_);
  }

  ROL::Ptr<ROL::Objective<Real>> buildObjective(void) {
    if (serial_) {
      return obj_;
    }
    else {
      if (useFilter_)
        return ROL::makePtr<Filtered_Compliance_Objective<Real>>(filter_, pde_, assembler_, pl_);
      else
        return ROL::makePtr<Compliance_Objective<Real>>(pde_, assembler_, pl_);
    }
  }

  ROL::Ptr<ROL::Vector<Real>> buildSolutionVector(void) {
    ROL::Ptr<ROL::Vector<Real>> z = z_->clone();
    z->set(*z_);
    return z;
  }

  void print(const ROL::Vector<Real> &z) {
    if (useFilter_) {
      ROL::Ptr<Filtered_Compliance_Objective<Real>>
        obj = ROL::makePtr<Filtered_Compliance_Objective<Real>>(filter_, pde_, mesh_, comm_, pl_, *os_);
      obj->printToFile(*dynamic_cast<const ROL::PartitionedVector<Real>&>(z).get(0),*os_);
    }
    else {
      ROL::Ptr<Compliance_Objective<Real>>
        obj = ROL::makePtr<Compliance_Objective<Real>>(pde_, mesh_, comm_, pl_, *os_);
      obj->printToFile(*dynamic_cast<const ROL::PartitionedVector<Real>&>(z).get(0),*os_);
    }
  }
};

#endif
