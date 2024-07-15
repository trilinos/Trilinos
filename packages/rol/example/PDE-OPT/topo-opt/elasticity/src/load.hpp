// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  load.hpp
    \brief Implements points loads for the structural topology
           optimization problem.
*/

#ifndef ROL_PDEOPT_ELASTICITY_LOAD_HPP
#define ROL_PDEOPT_ELASTICITY_LOAD_HPP

#include "ROL_Ptr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Intrepid_FieldContainer.hpp"

#include "../../../TOOLS/fe.hpp"

#include <vector>

template<class Real>
class Load {
private:
  std::vector<Real> loadMagnitude_;
  std::vector<Real> loadPolar_, loadAzimuth_;
  std::vector<Real> loadX_, loadY_, loadZ_;
  std::vector<Real> loadWidthX_, loadWidthY_, loadWidthZ_;

protected:
  Real DegreesToRadians(const Real deg) const {
    return deg * static_cast<Real>(M_PI) / static_cast<Real>(180);
  }

public:
  virtual ~Load() {}

  Load(void) {}

  Load(Teuchos::ParameterList & parlist, std::string ex = "Default") {
    if (ex == "2D Cantilever with 1 Load") {
      loadMagnitude_.push_back(static_cast<Real>(1));
      loadPolar_.push_back(static_cast<Real>(270));
      loadAzimuth_.push_back(static_cast<Real>(0));
      loadX_.push_back(static_cast<Real>(1));
      loadY_.push_back(static_cast<Real>(0));
      loadZ_.push_back(static_cast<Real>(0));
      loadWidthX_.push_back(static_cast<Real>(0.125));
      loadWidthY_.push_back(static_cast<Real>(0.125));
      loadWidthZ_.push_back(static_cast<Real>(0.125));
    }
    if (ex == "3D Cantilever") {
      loadMagnitude_.push_back(static_cast<Real>(1));
      loadPolar_.push_back(static_cast<Real>(270));
      loadAzimuth_.push_back(static_cast<Real>(90));
      loadX_.push_back(static_cast<Real>(1));
      loadY_.push_back(static_cast<Real>(0));
      loadZ_.push_back(static_cast<Real>(0));
      loadWidthX_.push_back(static_cast<Real>(0.125));
      loadWidthY_.push_back(static_cast<Real>(0.125));
      loadWidthZ_.push_back(static_cast<Real>(0.125));
    }
    else if (ex == "2D Truss") {
      loadMagnitude_.push_back(static_cast<Real>(1));
      loadPolar_.push_back(static_cast<Real>(270));
      loadAzimuth_.push_back(static_cast<Real>(0));
      loadX_.push_back(static_cast<Real>(2));
      loadY_.push_back(static_cast<Real>(0.5));
      loadZ_.push_back(static_cast<Real>(0));
      loadWidthX_.push_back(static_cast<Real>(1e-2));
      loadWidthY_.push_back(static_cast<Real>(1e-2));
      loadWidthZ_.push_back(static_cast<Real>(1e-2));
    }
    else if (ex == "2D Beams") {
      loadMagnitude_.push_back(static_cast<Real>(1));
      loadMagnitude_.push_back(static_cast<Real>(1));
      loadPolar_.push_back(static_cast<Real>(180));
      loadPolar_.push_back(static_cast<Real>(180));
      loadAzimuth_.push_back(static_cast<Real>(0));
      loadAzimuth_.push_back(static_cast<Real>(0));
      loadX_.push_back(static_cast<Real>(2));
      loadX_.push_back(static_cast<Real>(2));
      loadY_.push_back(static_cast<Real>(0));
      loadY_.push_back(static_cast<Real>(1));
      loadZ_.push_back(static_cast<Real>(0));
      loadZ_.push_back(static_cast<Real>(0));
      loadWidthX_.push_back(static_cast<Real>(1e-2));
      loadWidthX_.push_back(static_cast<Real>(1e-2));
      loadWidthY_.push_back(static_cast<Real>(1e-2));
      loadWidthY_.push_back(static_cast<Real>(1e-2));
      loadWidthZ_.push_back(static_cast<Real>(1e-2));
      loadWidthZ_.push_back(static_cast<Real>(1e-2));
    }
    else if (ex == "2D Cantilever with 3 Loads" ||
             ex == "2D Carrier Plate"           ||
             ex == "2D Wheel"                   ||
             ex == "3D L Beam") {
      loadMagnitude_.clear();
      loadPolar_.clear();
      loadAzimuth_.clear();
      loadX_.clear();
      loadY_.clear();
      loadZ_.clear();
      loadWidthX_.clear();
      loadWidthY_.clear();
      loadWidthZ_.clear();
    }
    else {
      if (parlist.isSublist("Load")) {
        // Grab Magnitudes
        Teuchos::Array<Real> magnitude
          = Teuchos::getArrayFromStringParameter<double>(parlist.sublist("Load"), "Magnitude");
        loadMagnitude_ = magnitude.toVector();
        // Grab Polar Angle
        Teuchos::Array<Real> polar
          = Teuchos::getArrayFromStringParameter<double>(parlist.sublist("Load"), "Polar Angle");
        loadPolar_ = polar.toVector();
        // Grab Azimuth Angle
        Teuchos::Array<Real> azimuth
          = Teuchos::getArrayFromStringParameter<double>(parlist.sublist("Load"), "Azimuth Angle");
        loadAzimuth_ = azimuth.toVector();
        // Grab X Location
        Teuchos::Array<Real> Xlocation
          = Teuchos::getArrayFromStringParameter<double>(parlist.sublist("Load"), "X Location");
        loadX_ = Xlocation.toVector();
        // Grab Y Location
        Teuchos::Array<Real> Ylocation
          = Teuchos::getArrayFromStringParameter<double>(parlist.sublist("Load"), "Y Location");
        loadY_ = Ylocation.toVector();
        // Grab Z Location
        Teuchos::Array<Real> Zlocation
          = Teuchos::getArrayFromStringParameter<double>(parlist.sublist("Load"), "Z Location");
        loadZ_ = Zlocation.toVector();
        // Grab X Width
        Teuchos::Array<Real> Xwidth
          = Teuchos::getArrayFromStringParameter<double>(parlist.sublist("Load"), "X Width");
        loadWidthX_ = Xwidth.toVector();
        // Grab Y Width
        Teuchos::Array<Real> Ywidth
          = Teuchos::getArrayFromStringParameter<double>(parlist.sublist("Load"), "Y Width");
        loadWidthY_ = Ywidth.toVector();
        // Grab Z Width
        Teuchos::Array<Real> Zwidth
          = Teuchos::getArrayFromStringParameter<double>(parlist.sublist("Load"), "Z Width");
        loadWidthZ_ = Zwidth.toVector();
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
    return (loadMagnitude_.size()==0);
  }

  virtual Real loadFunc(const std::vector<Real> &coords,
                        const int                dir,
                        const std::vector<Real> &param) const {
    Real val(0);
    const Real half(0.5);
    const int numLoads  = loadMagnitude_.size();
    const int paramSize = param.size();
    const int d         = coords.size();
    for (int i = 0; i < numLoads; ++i) {
      Real loadMagNoise(0), loadAngNoise0(0);
      if (paramSize > 0 + d*i) {
        loadMagNoise  = param[0 + d*i];
      }
      if (paramSize > 1 + d*i) {
        loadAngNoise0 = DegreesToRadians(param[1 + d*i]);
      }
      const Real loadMagnitude = loadMagnitude_[i] + loadMagNoise;
      const Real loadAngle0    = loadPolar_[i] + loadAngNoise0;
      const Real Gx = std::exp(-half*std::pow(coords[0]-loadX_[i],2)/std::pow(loadWidthX_[i],2));
      const Real Gy = std::exp(-half*std::pow(coords[1]-loadY_[i],2)/std::pow(loadWidthY_[i],2));

      if (d==2) {
        if (dir==0) {
          val += loadMagnitude*std::cos(loadAngle0)*Gx*Gy;
        }
        if (dir==1) {
          val += loadMagnitude*std::sin(loadAngle0)*Gx*Gy;
        }
      }
      if (d==3) {
        Real loadAngNoise1(0);
        if (paramSize > 2 + d*i) {
          loadAngNoise1 = DegreesToRadians(param[2 + d*i]);
        }
        const Real loadAngle1 = loadAzimuth_[i] + loadAngNoise1;
        const Real Gz = std::exp(-half*std::pow(coords[2]-loadZ_[i],2)/std::pow(loadWidthZ_[i],2));

        if (dir==0) {
          val += loadMagnitude*std::sin(loadAngle0)*std::cos(loadAngle1)*Gx*Gy*Gz;
        }
        if (dir==1) {
          val += loadMagnitude*std::sin(loadAngle0)*std::sin(loadAngle1)*Gx*Gy*Gz;
        }
        if (dir==2) {
          val += loadMagnitude*std::cos(loadAngle0)*Gx*Gy*Gz;
        }
      }
    }
    return val;
  }

  void compute(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > &load,
               const ROL::Ptr<FE<Real> >                               &fe,
               const std::vector<Real>                                     &param,
               const Real                                                   scale = 1) const {
    // Retrieve dimensions.
    int c = fe->gradN()->dimension(0);
    int p = fe->gradN()->dimension(2);
    int d = fe->gradN()->dimension(3);
    std::vector<Real> coord(d);

    for (int i=0; i<c; ++i) {
      for (int j=0; j<p; ++j) {
        for (int k=0; k<d; ++k) {
          coord[k] = (*fe->cubPts())(i,j,k);
        }
        for (int k=0; k<d; ++k) {
          (*load[k])(i,j) = loadFunc(coord, k, param)*scale;
        }
      }
    }
  }
};

#endif
