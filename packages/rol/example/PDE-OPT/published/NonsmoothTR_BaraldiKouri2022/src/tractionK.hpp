// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  traction.hpp
    \brief Implements traction loads for the structural topology
           optimization problem.
*/

#ifndef ROL_PDEOPT_ELASTICITY_TRACTIONK_HPP
#define ROL_PDEOPT_ELASTICITY_TRACTIONK_HPP

#include <vector>
#include "ROL_Ptr.hpp"
#include "ROL_ParameterList.hpp"
#include "../../../TOOLS/feK.hpp"

template<class Real, class DeviceType>
class Traction {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  std::vector<int> sidesets_, sideids_;
  std::vector<Real> loadMagnitude_;
  std::vector<Real> loadPolar_, loadAzimuth_;
  std::vector<std::vector<scalar_view>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  int dim_, offset_;

protected:
  Real DegreesToRadians(Real deg) const {
    return deg * static_cast<Real>(M_PI) / static_cast<Real>(180);
  }

public:
  virtual ~Traction() {}

  Traction(void) {};

  Traction(int dim) : dim_(dim) {
    TEUCHOS_TEST_FOR_EXCEPTION(dim > 3 || dim < 2, std::invalid_argument,
      ">>> PDE-OPT/published/NonsmoothTR_BaraldiKouri2022/src/traction.hpp: Problem dimension is not 2 or 3!");
    if (dim==2) {
      offset_ = 0;
      sidesets_.push_back(2);
      sideids_.push_back(0);
      loadMagnitude_.push_back(static_cast<Real>(1));
      loadPolar_.push_back(static_cast<Real>(270));
      loadAzimuth_.push_back(static_cast<Real>(0));
    }
    else {
      offset_ = 0;
      sidesets_.push_back(2);
      sideids_.push_back(1);
      loadMagnitude_.push_back(static_cast<Real>(1));
      loadPolar_.push_back(static_cast<Real>(270));
      loadAzimuth_.push_back(static_cast<Real>(90));
    }

    int polarSize = loadPolar_.size();
    for (int i=0; i<polarSize; ++i)
      loadPolar_[i] = DegreesToRadians(loadPolar_[i]);

    int azimuthSize = loadAzimuth_.size();
    for (int i=0; i<azimuthSize; ++i)
      loadAzimuth_[i] = DegreesToRadians(loadAzimuth_[i]);
  }

  bool isNull(void) const {
    return (sidesets_.size()==0);
  }

  void setCellNodes(const std::vector<std::vector<scalar_view>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds) {
    bdryCellNodes_  = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
  }

  virtual Real evaluate(const int                dir,
                        const int                sideset,
                        const std::vector<Real> &param) const {
//    const int numSideSets = sidesets_.size();
    const int paramSize   = param.size();
    Real val(0);
    Real loadMagNoise(0), loadAngNoise0(0);
    if (paramSize > 0 + dim_*sideset)
      loadMagNoise  = param[offset_ + 0 + dim_*sideset];
    if (paramSize > 1 + dim_*sideset)
      loadAngNoise0 = DegreesToRadians(param[offset_ + 1 + dim_*sideset]);
    const Real loadMagnitude = loadMagnitude_[sideset] + loadMagNoise;
    const Real loadAngle0    = loadPolar_[sideset]     + loadAngNoise0;

    if (dim_==2) {
      if (dir==0)
        val = loadMagnitude*std::cos(loadAngle0);
      if (dir==1)
        val = loadMagnitude*std::sin(loadAngle0);
    }
    if (dim_==3) {
      Real loadAngNoise1(0);
      if (paramSize > 2 + dim_*sideset) {
        loadAngNoise1 = DegreesToRadians(param[offset_ + 2 + dim_*sideset]);
      }
      const Real loadAngle1 = loadAzimuth_[sideset] + loadAngNoise1;

      if (dir==0)
        val = loadMagnitude*std::sin(loadAngle0)*std::cos(loadAngle1);
      if (dir==1)
        val = loadMagnitude*std::sin(loadAngle0)*std::sin(loadAngle1);
      if (dir==2)
        val = loadMagnitude*std::cos(loadAngle0);
    }
    return val;
  }

  void compute(std::vector<std::vector<std::vector<scalar_view>>> &traction,
               const std::vector<std::vector<ROL::Ptr<fe_type>>> &fe,
               const std::vector<Real> &param,
               Real scale = Real(1)) const {
    traction.clear();
    traction.resize(bdryCellLocIds_.size());
    const int numSideSets = sidesets_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        const int numLocalSideIds = bdryCellLocIds_[sidesets_[i]].size();
        traction[sidesets_[i]].resize(numLocalSideIds);
        for (int j = 0; j < numLocalSideIds; ++j) {
          const int numCellsSide = bdryCellLocIds_[sidesets_[i]][j].size();
          if (numCellsSide > 0) {
            const int numCubPerSide = fe[sidesets_[i]][j]->cubPts().extent_int(1);
            traction[sidesets_[i]][j].resize(dim_);
            for (int k = 0; k < dim_; ++k) {
              traction[sidesets_[i]][j][k] = scalar_view("traction", numCellsSide, numCubPerSide);
              for (int c = 0; c < numCellsSide; ++c) {
                for (int p = 0; p < numCubPerSide; ++p)
                  (traction[sidesets_[i]][j][k])(c,p) = scale*evaluate(k,i,param);
              }
            }
          }
        }
      }
    }
  }

  void apply(std::vector<scalar_view> &R,
             const std::vector<std::vector<ROL::Ptr<fe_type>>> &fe,
             const std::vector<Real> &param,
             Real scale = Real(1)) const {
    const int numSideSets = sidesets_.size();
    const int nf = R[0].extent_int(1);
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        const int numLocalSideIds = bdryCellLocIds_[sidesets_[i]].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          const int numCellsSide = bdryCellLocIds_[sidesets_[i]][j].size();
          if (numCellsSide > 0) {
            const int numCubPerSide = fe[sidesets_[i]][j]->cubPts().extent_int(1);
            for (int k = 0; k < dim_; ++k) {
              scalar_view traction("traction", numCellsSide, numCubPerSide);
              for (int c = 0; c < numCellsSide; ++c) {
                for (int p = 0; p < numCubPerSide; ++p)
                  traction(c,p) = scale*evaluate(k,i,param);
              }
              scalar_view tractionResidual("tractionResidual", numCellsSide, nf);
              fst::integrate(tractionResidual,traction,fe[sidesets_[i]][j]->NdetJ(),false);
              // Add traction residual to volume residual
              for (int l = 0; l < numCellsSide; ++l) {
                int cidx = bdryCellLocIds_[sidesets_[i]][j][l];
                for (int m = 0; m < nf; ++m)
                  (R[k])(cidx,m) += tractionResidual(l,m);
              }
            }
          }
        }
      }
    }
  }
};

#endif
