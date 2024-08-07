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

#ifndef ROL_PDEOPT_ELASTICITY_TRACTION_HPP
#define ROL_PDEOPT_ELASTICITY_TRACTION_HPP

#include "ROL_Ptr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Intrepid_FieldContainer.hpp"

#include "../../../TOOLS/fe.hpp"

#include <vector>

template<class Real>
class Traction {
private:
  std::vector<int> sidesets_, sideids_;
  std::vector<Real> loadMagnitude_;
  std::vector<Real> loadPolar_, loadAzimuth_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  int dim_, offset_;

protected:
  Real DegreesToRadians(const Real deg) const {
    return deg * static_cast<Real>(M_PI) / static_cast<Real>(180);
  }

public:
  virtual ~Traction() {}

  Traction(void) {};

  Traction(int dim) : dim_(dim) {
    assert(dim > 3 || dim < 2);
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

  void setCellNodes(const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &bdryCellNodes,
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
    if (paramSize > 0 + dim_*sideset) {
      loadMagNoise  = param[offset_ + 0 + dim_*sideset];
    }
    if (paramSize > 1 + dim_*sideset) {
      loadAngNoise0 = DegreesToRadians(param[offset_ + 1 + dim_*sideset]);
    }
    const Real loadMagnitude = loadMagnitude_[sideset] + loadMagNoise;
    const Real loadAngle0    = loadPolar_[sideset]     + loadAngNoise0;

    if (dim_==2) {
      if (dir==0) {
        val = loadMagnitude*std::cos(loadAngle0);
      }
      if (dir==1) {
        val = loadMagnitude*std::sin(loadAngle0);
      }
    }
    if (dim_==3) {
      Real loadAngNoise1(0);
      if (paramSize > 2 + dim_*sideset) {
        loadAngNoise1 = DegreesToRadians(param[offset_ + 2 + dim_*sideset]);
      }
      const Real loadAngle1 = loadAzimuth_[sideset] + loadAngNoise1;

      if (dir==0) {
        val = loadMagnitude*std::sin(loadAngle0)*std::cos(loadAngle1);
      }
      if (dir==1) {
        val = loadMagnitude*std::sin(loadAngle0)*std::sin(loadAngle1);
      }
      if (dir==2) {
        val = loadMagnitude*std::cos(loadAngle0);
      }
    }
    return val;
  }

  void compute(std::vector<std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>>> &traction,
               const std::vector<std::vector<ROL::Ptr<FE<Real>>>>                              &fe,
               const std::vector<Real>                                                         &param,
               const Real                                                                       scale = Real(1)) const {
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
            const int numCubPerSide = fe[sidesets_[i]][j]->cubPts()->dimension(1);
            traction[sidesets_[i]][j].resize(dim_);
            for (int k = 0; k < dim_; ++k) {
              traction[sidesets_[i]][j][k]
                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
              for (int c = 0; c < numCellsSide; ++c) {
                for (int p = 0; p < numCubPerSide; ++p) {
                  (*traction[sidesets_[i]][j][k])(c,p) = scale*evaluate(k,i,param);
                }
              }
            }
          }
        }
      }
    }
  }

  void apply(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &R,
             const std::vector<std::vector<ROL::Ptr<FE<Real>>>>    &fe,
             const std::vector<Real>                               &param,
             const Real                                             scale = Real(1)) const {
    const int numSideSets = sidesets_.size();
    const int nf          = R[0]->dimension(1);
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        const int numLocalSideIds = bdryCellLocIds_[sidesets_[i]].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          const int numCellsSide = bdryCellLocIds_[sidesets_[i]][j].size();
          if (numCellsSide > 0) {
            const int numCubPerSide = fe[sidesets_[i]][j]->cubPts()->dimension(1);
            for (int k = 0; k < dim_; ++k) {
              Intrepid::FieldContainer<Real> traction(numCellsSide, numCubPerSide);
              for (int c = 0; c < numCellsSide; ++c) {
                for (int p = 0; p < numCubPerSide; ++p) {
                  traction(c,p) = scale*evaluate(k,i,param);
                }
              }
              Intrepid::FieldContainer<Real> tractionResidual(numCellsSide, nf);
              Intrepid::FunctionSpaceTools::integrate<Real>(tractionResidual,
                                                            traction,
                                                            *fe[sidesets_[i]][j]->NdetJ(),
                                                            Intrepid::COMP_CPP, false);
              // Add traction residual to volume residual
              for (int l = 0; l < numCellsSide; ++l) {
                int cidx = bdryCellLocIds_[sidesets_[i]][j][l];
                for (int m = 0; m < nf; ++m) {
                  (*R[k])(cidx,m) += tractionResidual(l,m);
                }
              }
            }
          }
        }
      }
    }
  }
};

#endif
