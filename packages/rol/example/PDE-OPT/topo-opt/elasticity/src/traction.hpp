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
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
  std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  int dim_, offset_;

protected:
  Real DegreesToRadians(const Real deg) const {
    return deg * static_cast<Real>(M_PI) / static_cast<Real>(180);
  }

public:
  virtual ~Traction() {}

  Traction(void) {};

  Traction(Teuchos::ParameterList & parlist, std::string ex = "Default") {
    if (ex == "2D Wheel") {
      dim_ = 2;
      offset_ = 0;
      sidesets_.push_back(2);
      sideids_.push_back(0);
      loadMagnitude_.push_back(static_cast<Real>(1));
      loadPolar_.push_back(static_cast<Real>(270));
      loadAzimuth_.push_back(static_cast<Real>(0));
    }
    else if (ex == "2D Cantilever with 1 Load" ||
             ex == "2D Truss") {
      dim_ = 2;
      offset_ = 2;
      sidesets_.clear();
      sideids_.clear();
      loadMagnitude_.clear();
      loadPolar_.clear();
      loadAzimuth_.clear();
    }
    else if (ex == "2D Cantilever with 3 Loads") {
      dim_ = 2;
      offset_ = 0;
      sidesets_.push_back(1);
      sidesets_.push_back(3);
      sidesets_.push_back(5);
      sideids_.push_back(0);
      sideids_.push_back(0);
      sideids_.push_back(0);
      loadMagnitude_.push_back(static_cast<Real>(1));
      loadMagnitude_.push_back(static_cast<Real>(1));
      loadMagnitude_.push_back(static_cast<Real>(1));
      loadPolar_.push_back(static_cast<Real>(270));
      loadPolar_.push_back(static_cast<Real>(270));
      loadPolar_.push_back(static_cast<Real>(270));
      loadAzimuth_.push_back(static_cast<Real>(0));
      loadAzimuth_.push_back(static_cast<Real>(0));
      loadAzimuth_.push_back(static_cast<Real>(0));
    }
    else if (ex == "2D Beams") {
      dim_ = 2;
      offset_ = 4;
      sidesets_.clear();
      sideids_.clear();
      loadMagnitude_.clear();
      loadPolar_.clear();
      loadAzimuth_.clear();
    }
    else if (ex == "2D Carrier Plate") {
      dim_ = 2;
      offset_ = 0;
      sidesets_.push_back(2);
      sidesets_.push_back(4);
      sideids_.push_back(2);
      sideids_.push_back(2);
      loadMagnitude_.push_back(static_cast<Real>(1.0));
      loadMagnitude_.push_back(static_cast<Real>(0.5));
      loadPolar_.push_back(static_cast<Real>(270));
      loadPolar_.push_back(static_cast<Real>(270));
      loadAzimuth_.push_back(static_cast<Real>(0));
      loadAzimuth_.push_back(static_cast<Real>(0));
    }
    else if (ex == "3D Cantilever") {
      dim_ = 3;
      offset_ = 3;
      sidesets_.clear();
      sideids_.clear();
      loadMagnitude_.clear();
      loadPolar_.clear();
      loadAzimuth_.clear();
    }
    else {
      if (parlist.isSublist("Traction")) {
        // Grab problem dimension
        dim_ = parlist.get("Problem Dimension",2);

        offset_ = 0;
        if (parlist.isSublist("Load")) {
          Teuchos::Array<Real> loadMag
            = Teuchos::getArrayFromStringParameter<double>(parlist.sublist("Load"), "Magnitude");
          offset_ = dim_*static_cast<int>(loadMag.size());
        }

        // Grab Sidesets
        Teuchos::Array<int> sidesets
          = Teuchos::getArrayFromStringParameter<int>(parlist.sublist("Traction"), "Sidesets");
        sidesets_ = sidesets.toVector();

        // Grab Side IDs
        Teuchos::Array<int> sideids
          = Teuchos::getArrayFromStringParameter<int>(parlist.sublist("Traction"), "Side IDs");
        sideids_ = sideids.toVector();

        // Grab Magnitudes
        Teuchos::Array<Real> magnitude
          = Teuchos::getArrayFromStringParameter<double>(parlist.sublist("Traction"), "Magnitude");
        loadMagnitude_ = magnitude.toVector();

        // Grab Polar Angle
        Teuchos::Array<Real> polar
          = Teuchos::getArrayFromStringParameter<double>(parlist.sublist("Traction"), "Polar Angle");
        loadPolar_ = polar.toVector();

        if (dim_ == 3) {
          // Grab Azimuth Angle
          Teuchos::Array<Real> azimuth
            = Teuchos::getArrayFromStringParameter<double>(parlist.sublist("Traction"), "Azimuth Angle");
          loadAzimuth_ = azimuth.toVector();
        }
      }
    }

    int polarSize = loadPolar_.size();
    for (int i=0; i<polarSize; ++i) {
      loadPolar_[i] = DegreesToRadians(loadPolar_[i]);
    }

    int azimuthSize = loadAzimuth_.size();
    for (int i=0; i<azimuthSize; ++i) {
      loadAzimuth_[i] = DegreesToRadians(loadAzimuth_[i]);
    }
  }

  bool isNull(void) const {
    return (sidesets_.size()==0);
  }

  void setCellNodes(const std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds) {
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

  void compute(std::vector<std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > > &traction,
               const std::vector<std::vector<ROL::Ptr<FE<Real> > > >                               &fe,
               const std::vector<Real>                                                                 &param,
               const Real                                                                               scale = 1) const {
    traction.clear();
    traction.resize(bdryCellLocIds_.size());
    const int numSideSets = sidesets_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        traction[sidesets_[i]].resize(bdryCellLocIds_[sidesets_[i]].size());
        const int numCellsSide = bdryCellLocIds_[sidesets_[i]][sideids_[i]].size();
        if (numCellsSide > 0) {
          traction[sidesets_[i]][sideids_[i]].resize(dim_);
          const int numCubPerSide = fe[0][0]->gradN()->dimension(2);
          for (int k = 0; k < dim_; ++k) {
            traction[sidesets_[i]][sideids_[i]][k]
              = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
            for (int c = 0; c < numCellsSide; ++c) {
              for (int p = 0; p < numCubPerSide; ++p) {
                (*traction[sidesets_[i]][sideids_[i]][k])(c,p) = scale*evaluate(k,i,param);
              }
            }
          }
        }
      }
    }
  }

  void apply(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > &R,
             const std::vector<std::vector<ROL::Ptr<FE<Real> > > >   &fe,
             const std::vector<Real>                                     &param,
             const Real                                                   scale = 1) const {
    const int numSideSets = sidesets_.size();
    const int nf          = R[0]->dimension(1);
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        const int numCellsSide = bdryCellLocIds_[sidesets_[i]][sideids_[i]].size();
        if (numCellsSide > 0) {
          const int numCubPerSide = fe[0][0]->gradN()->dimension(2);
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
                                                          *(fe[sidesets_[i]][sideids_[i]]->NdetJ()),
                                                          Intrepid::COMP_CPP, false);
            // Add Robin residual to volume residual
            for (int l = 0; l < numCellsSide; ++l) {
              int cidx = bdryCellLocIds_[sidesets_[i]][sideids_[i]][l];
              for (int m = 0; m < nf; ++m) { 
                (*R[k])(cidx,m) += tractionResidual(l,m);
              }
            }
          }
        }
      }
    }
  }
};

#endif
