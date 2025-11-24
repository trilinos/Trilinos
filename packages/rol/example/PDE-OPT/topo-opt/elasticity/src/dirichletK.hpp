// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  dirichlet.hpp
    \brief Implements Dirichlet boundary conditions for the structural
           topology optimization problem.
*/

#ifndef ROL_PDEOPT_ELASTICITY_DIRICHLETK_HPP
#define ROL_PDEOPT_ELASTICITY_DIRICHLETK_HPP

#include "ROL_Ptr.hpp"
#include "ROL_ParameterList.hpp"
#include "../../../TOOLS/feK.hpp"
#include <vector>

template<class Real, class DeviceType>
class Dirichlet {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;

private:
  std::vector<int> sidesets_, types_;
  std::vector<std::vector<scalar_view>> bdryCellNodes_;
  std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  std::vector<std::vector<int>> fidx_;

  std::vector<int> getConstrainedDimensions(int type, int dim) const {
    std::vector<int> conDim;
    if ( type == 1 ) {
      conDim.push_back(0);
    }
    else if ( type == 2 ) {
      conDim.push_back(1);
    }
    else if ( type == 3 && dim > 2 ) {
      conDim.push_back(2);
    }
    else if ( type == 4 && dim > 2 ) {
      conDim.push_back(0);
      conDim.push_back(1);
    }
    else if ( type == 5 && dim > 2 ) {
      conDim.push_back(0);
      conDim.push_back(2);
    }
    else if ( type == 6 && dim > 2 ) {
      conDim.push_back(1);
      conDim.push_back(2);
    }
    else {
      for (int i = 0; i < dim; ++i) {
        conDim.push_back(i);
      }
    }
    return conDim;
  }

public:
  Dirichlet(ROL::ParameterList &parlist, std::string ex = "Default") {
    if (ex == "2D Cantilever with 1 Load"  ||
        ex == "3D Cantilever"              ||
        ex == "2D Truss"                   ||
        ex == "2D Beams") {
      sidesets_.push_back(3);
      types_.push_back(0);
    }
    else if (ex == "2D Cantilever with 3 Loads") {
      sidesets_.push_back(8);
      types_.push_back(0);
    }
    else if (ex == "2D Wheel") {
      sidesets_.push_back(0);
      sidesets_.push_back(4);
      types_.push_back(0);
      types_.push_back(2);
    }
    else if (ex == "2D Carrier Plate") {
      sidesets_.push_back(0);
      types_.push_back(0);
    }
    else if (ex == "3D L Beam") {
      sidesets_.push_back(0);
      types_.push_back(0);
    }
    else {
      // Grab sidesets
      sidesets_ = ROL::getArrayFromStringParameter<int>(parlist.sublist("Dirichlet"), "Sidesets");
      // Grab type
      types_ = ROL::getArrayFromStringParameter<int>(parlist.sublist("Dirichlet"), "Types");
    }
  }

  void setCellNodes(const std::vector<std::vector<scalar_view>> &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds,
                    const std::vector<std::vector<int>> &fidx) {
    bdryCellNodes_  = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    fidx_           = fidx;
  }

  void applyResidual(std::vector<scalar_view> &R,
                     const std::vector<scalar_view> &U) const {
    const int d = R.size();
    const int numSideSets = sidesets_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        const std::vector<int> conDim = getConstrainedDimensions(types_[i],d);
        const int numConDim = conDim.size();
        const int numLocalSideIds = bdryCellLocIds_[sidesets_[i]].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          const int numCellsSide = bdryCellLocIds_[sidesets_[i]][j].size();
          const int numBdryDofs = fidx_[j].size();
          for (int k = 0; k < numCellsSide; ++k) {
            const int cidx = bdryCellLocIds_[sidesets_[i]][j][k];
            for (int l = 0; l < numBdryDofs; ++l) {
              for (int m=0; m < numConDim; ++m) {
                (R[conDim[m]])(cidx,fidx_[j][l]) = (U[conDim[m]])(cidx,fidx_[j][l]);
              }
            }
          }
        }
      }
    }
  }

  void applyJacobian1(std::vector<std::vector<scalar_view>> &J) const {
    const int d = J.size();
    const int f = J[0][0].extent_int(1);
    const int numSideSets = sidesets_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        const std::vector<int> conDim = getConstrainedDimensions(types_[i],d);
        const int numConDim = conDim.size();
        const int numLocalSideIds = bdryCellLocIds_[sidesets_[i]].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          const int numCellsSide = bdryCellLocIds_[sidesets_[i]][j].size();
          const int numBdryDofs = fidx_[j].size();
          for (int k = 0; k < numCellsSide; ++k) {
            const int cidx = bdryCellLocIds_[sidesets_[i]][j][k];
            for (int l = 0; l < numBdryDofs; ++l) {
              for (int m=0; m < f; ++m) {
                for (int n=0; n < numConDim; ++n) {
                  for (int p=0; p < d; ++p) {
                    (J[conDim[n]][p])(cidx,fidx_[j][l],m) = static_cast<Real>(0);
                  }
                  (J[conDim[n]][conDim[n]])(cidx,fidx_[j][l],fidx_[j][l]) = static_cast<Real>(1);
                }
              }
            }
          }
        }
      }
    }
  }

  void applyJacobian2(std::vector<std::vector<scalar_view>> &J) const {
    const int d = J.size();
    const int f = J[0][0].extent_int(2);
    const int numSideSets = sidesets_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        const std::vector<int> conDim = getConstrainedDimensions(types_[i],d);
        const int numConDim = conDim.size();
        const int numLocalSideIds = bdryCellLocIds_[sidesets_[i]].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          const int numCellsSide = bdryCellLocIds_[sidesets_[i]][j].size();
          const int numBdryDofs = fidx_[j].size();
          for (int k = 0; k < numCellsSide; ++k) {
            const int cidx = bdryCellLocIds_[sidesets_[i]][j][k];
            for (int l = 0; l < numBdryDofs; ++l) {
              for (int m=0; m < f; ++m) {
                for (int n=0; n < numConDim; ++n) {
                  (J[conDim[n]][0])(cidx,fidx_[j][l],m) = static_cast<Real>(0);
                }
              }
            }
          }
        }
      }
    }
  }

  void applyMultiplier(std::vector<scalar_view> &L) const {
    const int d = L.size();
    const int numSideSets = sidesets_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        const std::vector<int> conDim = getConstrainedDimensions(types_[i],d);
        const int numConDim = conDim.size();
        const int numLocalSideIds = bdryCellLocIds_[sidesets_[i]].size();
        for (int j = 0; j < numLocalSideIds; ++j) {
          const int numCellsSide = bdryCellLocIds_[sidesets_[i]][j].size();
          const int numBdryDofs = fidx_[j].size();
          for (int k = 0; k < numCellsSide; ++k) {
            const int cidx = bdryCellLocIds_[sidesets_[i]][j][k];
            for (int l = 0; l < numBdryDofs; ++l) {
              for (int m=0; m < numConDim; ++m) {
                (L[conDim[m]])(cidx,fidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
      }
    }
  }
};

#endif
