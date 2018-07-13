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

/*! \file  sampler.hpp
    \brief Generates sampler for stochastic structural
           topology optimization problem.
*/

#ifndef ROL_PDEOPT_ELASTICITY_SAMPLER_HPP
#define ROL_PDEOPT_ELASTICITY_SAMPLER_HPP

#include "ROL_Ptr.hpp"
#include "ROL_DistributionFactory.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "Teuchos_ParameterList.hpp"

template<class Real>
class BuildSampler {
private:
  int stochDim_;
  std::vector<ROL::Ptr<ROL::Distribution<Real>>> distVec_;

public:
  BuildSampler(Teuchos::ParameterList &parlist, int probDim, std::string ex = "Default") : stochDim_(0) {
    if (ex == "2D Wheel") {
      stochDim_ = probDim;
      // Magnitude distribution.
      Teuchos::ParameterList magList;
      magList.sublist("Distribution").set("Name","Uniform");
      magList.sublist("Distribution").sublist("Uniform").set("Lower Bound",-0.5);
      magList.sublist("Distribution").sublist("Uniform").set("Upper Bound", 0.5);
      distVec_.push_back(ROL::DistributionFactory<Real>(magList));
      // Polar angle distribution.
      Teuchos::ParameterList polList;
      polList.sublist("Distribution").set("Name","Truncated Gaussian");
      polList.sublist("Distribution").sublist("Truncated Gaussian").set("Mean",                0.0);
      polList.sublist("Distribution").sublist("Truncated Gaussian").set("Standard Deviation", 17.0);
      polList.sublist("Distribution").sublist("Truncated Gaussian").set("Lower Bound",       -90.0);
      polList.sublist("Distribution").sublist("Truncated Gaussian").set("Upper Bound",        90.0);
      distVec_.push_back(ROL::DistributionFactory<Real>(polList));
    }
    else if (ex == "2D Truss" || ex == "3D Cantilever") {
      stochDim_ = probDim;
      // Magnitude distribution.
      Teuchos::ParameterList magList;
      magList.sublist("Distribution").set("Name","Uniform");
      magList.sublist("Distribution").sublist("Uniform").set("Lower Bound",-0.5);
      magList.sublist("Distribution").sublist("Uniform").set("Upper Bound", 0.5);
      distVec_.push_back(ROL::DistributionFactory<Real>(magList));
      // Polar angle distribution.
      Teuchos::ParameterList polList;
      polList.sublist("Distribution").set("Name","Uniform");
      polList.sublist("Distribution").sublist("Uniform").set("Lower Bound",-45.0);
      polList.sublist("Distribution").sublist("Uniform").set("Upper Bound", 45.0);
      distVec_.push_back(ROL::DistributionFactory<Real>(polList));
      if (probDim == 3) {
        Teuchos::ParameterList aziList;
        aziList.sublist("Distribution").set("Name","Uniform");
        aziList.sublist("Distribution").sublist("Uniform").set("Lower Bound",-45.0);
        aziList.sublist("Distribution").sublist("Uniform").set("Upper Bound", 45.0);
        distVec_.push_back(ROL::DistributionFactory<Real>(aziList));
      }
    }
    else if (ex == "2D Beams") {
      stochDim_ = 2*probDim;
      // Load 1: Magnitude distribution.
      Teuchos::ParameterList magList1;
      magList1.sublist("Distribution").set("Name","Uniform");
      magList1.sublist("Distribution").sublist("Uniform").set("Lower Bound",-0.25);
      magList1.sublist("Distribution").sublist("Uniform").set("Upper Bound", 0.25);
      distVec_.push_back(ROL::DistributionFactory<Real>(magList1));
      // Load 1: Polar angle distribution.
      Teuchos::ParameterList polList1;
      polList1.sublist("Distribution").set("Name","Uniform");
      polList1.sublist("Distribution").sublist("Uniform").set("Lower Bound",-25.0);
      polList1.sublist("Distribution").sublist("Uniform").set("Upper Bound", 25.0);
      distVec_.push_back(ROL::DistributionFactory<Real>(polList1));
      // Load 2: Magnitude distribution.
      Teuchos::ParameterList magList2;
      magList2.sublist("Distribution").set("Name","Uniform");
      magList2.sublist("Distribution").sublist("Uniform").set("Lower Bound",-0.75);
      magList2.sublist("Distribution").sublist("Uniform").set("Upper Bound", 0.75);
      distVec_.push_back(ROL::DistributionFactory<Real>(magList2));
      // Load 2: Polar angle distribution.
      Teuchos::ParameterList polList2;
      polList2.sublist("Distribution").set("Name","Uniform");
      polList2.sublist("Distribution").sublist("Uniform").set("Lower Bound",-45.0);
      polList2.sublist("Distribution").sublist("Uniform").set("Upper Bound", 45.0);
      distVec_.push_back(ROL::DistributionFactory<Real>(polList2));
    }
    else {
      if (parlist.isSublist("Load")) {
        Teuchos::Array<Real> loadMag
          = Teuchos::getArrayFromStringParameter<double>(parlist.sublist("Load"), "Magnitude");
        int nLoads = loadMag.size();
        stochDim_  = probDim*nLoads;
        for (int i = 0; i < nLoads; ++i) {
          std::stringstream sli;
          sli << "Stochastic Load " << i;
          // Magnitude distribution.
          Teuchos::ParameterList magList;
          magList.sublist("Distribution") = parlist.sublist(sli.str()).sublist("Magnitude");
          distVec_.push_back(ROL::DistributionFactory<Real>(magList));
          // Polar angle distribution.
          Teuchos::ParameterList polList;
          polList.sublist("Distribution") = parlist.sublist(sli.str()).sublist("Polar Angle");
          distVec_.push_back(ROL::DistributionFactory<Real>(polList));
          // Azimuth angle distribution.
          if (probDim == 3) {
            Teuchos::ParameterList aziList;
            aziList.sublist("Distribution") = parlist.sublist(sli.str()).sublist("Azimuth Angle");
            distVec_.push_back(ROL::DistributionFactory<Real>(aziList));
          }
        }
      }
      if (parlist.isSublist("Traction")) {
        Teuchos::Array<Real> tractionMag
          = Teuchos::getArrayFromStringParameter<double>(parlist.sublist("Traction"), "Magnitude");
        int nTractions = tractionMag.size();
        stochDim_     += probDim*nTractions;
        for (int i = 0; i < nTractions; ++i) {
          std::stringstream sli;
          sli << "Stochastic Traction " << i;
          // Magnitude distribution.
          Teuchos::ParameterList magList;
          magList.sublist("Distribution") = parlist.sublist(sli.str()).sublist("Magnitude");
          distVec_.push_back(ROL::DistributionFactory<Real>(magList));
          // Polar angle distribution.
          Teuchos::ParameterList polList;
          polList.sublist("Distribution") = parlist.sublist(sli.str()).sublist("Polar Angle");
          distVec_.push_back(ROL::DistributionFactory<Real>(polList));
          // Azimuth angle distribution.
          if (probDim == 3) {
            Teuchos::ParameterList aziList;
            aziList.sublist("Distribution") = parlist.sublist(sli.str()).sublist("Azimuth Angle");
            distVec_.push_back(ROL::DistributionFactory<Real>(aziList));
          }
        }
      }
    }
  }

  ROL::Ptr<ROL::SampleGenerator<Real>> get(const int nsamp,
                                           const ROL::Ptr<ROL::BatchManager<Real>> &bman) const {
    return ROL::makePtr<ROL::MonteCarloGenerator<Real>>(nsamp,distVec_,bman);
  }

  int getDimension(void) const {
    return stochDim_;
  }

};

#endif
