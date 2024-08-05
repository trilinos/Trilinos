// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef POTENTIALFLOW_HPP
#define POTENTIALFLOW_HPP

#include "../../TOOLS/fe.hpp"
#include "../../TOOLS/fieldhelper.hpp"
#include "Tpetra_MultiVector.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Ptr.hpp"

template <class Real>
class PotentialFlow {
private:
  typedef Tpetra::Map<>::global_ordinal_type GO;  

  const ROL::Ptr<FE<Real>> feVel_, fePrs_;
  const ROL::Ptr<Intrepid::FieldContainer<Real>> cellNodes_;
  const ROL::Ptr<Intrepid::FieldContainer<int>> cellDofs_;
  const Teuchos::Array<GO> cellIds_;
  const ROL::Ptr<FieldHelper<Real>> fieldHelper_;
  Real cx_, cy_, r_;

public:
  PotentialFlow(const ROL::Ptr<FE<Real>> &feVel,
                const ROL::Ptr<FE<Real>> &fePrs,
                const ROL::Ptr<Intrepid::FieldContainer<Real>> &cellNodes,
                const ROL::Ptr<Intrepid::FieldContainer<int>> &cellDofs,
                const Teuchos::Array<GO> &cellIds,
                const ROL::Ptr<FieldHelper<Real>> &fieldHelper,
                ROL::ParameterList &parlist)
    : feVel_(feVel), fePrs_(fePrs),
      cellNodes_(cellNodes), cellDofs_(cellDofs), cellIds_(cellIds),
      fieldHelper_(fieldHelper) {
    // Get parameters from parameter list
    cx_ = parlist.sublist("Problem").get("Cylinder Center X", 0.0);
    cy_ = parlist.sublist("Problem").get("Cylinder Center Y", 0.0);
    r_  = parlist.sublist("Problem").get("Cylinder Radius",   0.5);
  }

  Real evalVelocity(const std::vector<Real> &x, const int dir) const {
    const Real one(1);
    Real px  = x[0] - cx_;
    Real py  = x[1] - cy_;
    Real rad = std::sqrt(std::pow(px,2) + std::pow(py,2));
    Real arg = std::atan2(py,px);
    Real rdt =  (one - std::pow(r_/rad,2))*std::cos(arg);
    Real adt = -(one + std::pow(r_/rad,2))*std::sin(arg);
    return (dir == 0 ? rdt*std::cos(arg) - adt*std::sin(arg)
                     : rdt*std::sin(arg) + adt*std::cos(arg));
  }

  Real evalPressure(const std::vector<Real> &x) const {
    const Real half(0.5), one(1);
    Real vx = evalVelocity(x,0);
    Real vy = evalVelocity(x,1);
    return half*(one-(std::pow(vx,2)+std::pow(vy,2)));
  }

  void build(const ROL::Ptr<Tpetra::MultiVector<>> &vec) const {
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int d  = feVel_->gradN()->dimension(3);
    // Grab DOF coordinates
    ROL::Ptr<Intrepid::FieldContainer<Real>> dofPointsVel, dofPointsPrs, u;
    dofPointsVel = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,fv,d);
    dofPointsPrs = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,fp,d);
    feVel_->computeDofCoords(dofPointsVel, cellNodes_);
    fePrs_->computeDofCoords(dofPointsPrs, cellNodes_);
    // Initialize velocity and pressure vectors
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U(d+1);
    for (int i = 0; i < d; ++i) {
      U[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,fv);
    }
    U[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,fp);
    // Evaluate velocity and pressure
    std::vector<Real> coord(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < fv; ++j) {
        for (int k = 0; k < d; ++k) {
          coord[k] = (*dofPointsVel)(i,j,k);
        }
        for (int k = 0; k < d; ++k) {
          (*U[k])(i,j) = evalVelocity(coord,k);
        }
      }
      for (int j = 0; j < fp; ++j) {
        for (int k = 0; k < d; ++k) {
          coord[k] = (*dofPointsPrs)(i,j,k);
        }
        (*U[d])(i,j) = evalPressure(coord);
      }
    }
    // Combine velocity and pressure into a single vector
    fieldHelper_->combineFieldCoeff(u,U);
    // Add velocity and pressure values to input vector
    int numDofs = u->dimension(1);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < numDofs; ++j) {
        vec->replaceGlobalValue((*cellDofs_)(cellIds_[i],j),0,(*u)(i,j));
      }
    }
  }
};

#endif
