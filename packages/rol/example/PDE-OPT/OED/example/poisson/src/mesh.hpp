// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "../../../../TOOLS/meshmanager.hpp"

template <class Real>
class MeshManager_OED_Poisson : public MeshManager_Rectangle<Real> {

private:

  int nx_;
  int ny_;
  ROL::Ptr<std::vector<std::vector<std::vector<int>>>>  meshSideSets_;

public: 

  MeshManager_OED_Poisson(Teuchos::ParameterList &parlist) : MeshManager_Rectangle<Real>(parlist)
  {
    nx_ = parlist.sublist("Geometry").get("NX", 3);
    ny_ = parlist.sublist("Geometry").get("NY", 3);
    computeSideSets();
  }


  void computeSideSets() {

    int numSideSets = 4;
    meshSideSets_ = ROL::makePtr<std::vector<std::vector<std::vector<int> > >>(numSideSets);

    // Neumann 0
    (*meshSideSets_)[0].resize(4);
    (*meshSideSets_)[0][0].resize(nx_);
    (*meshSideSets_)[0][1].resize(0);
    (*meshSideSets_)[0][2].resize(0);
    (*meshSideSets_)[0][3].resize(0);
    // Neumann 1
    (*meshSideSets_)[1].resize(4);
    (*meshSideSets_)[1][0].resize(0);
    (*meshSideSets_)[1][1].resize(ny_);
    (*meshSideSets_)[1][2].resize(0);
    (*meshSideSets_)[1][3].resize(0);
    // Neumann 2
    (*meshSideSets_)[2].resize(4);
    (*meshSideSets_)[2][0].resize(0);
    (*meshSideSets_)[2][1].resize(0);
    (*meshSideSets_)[2][2].resize(nx_);
    (*meshSideSets_)[2][3].resize(0);
    // Neumann 3
    (*meshSideSets_)[3].resize(4);
    (*meshSideSets_)[3][0].resize(0);
    (*meshSideSets_)[3][1].resize(0);
    (*meshSideSets_)[3][2].resize(0);
    (*meshSideSets_)[3][3].resize(ny_);
    
    for (int i=0; i<nx_; ++i) {
      (*meshSideSets_)[0][0][i] = i;
    }
    for (int i=0; i<ny_; ++i) {
      (*meshSideSets_)[1][1][i] = (i+1)*nx_-1;
    }
    for (int i=0; i<nx_; ++i) {
      (*meshSideSets_)[2][2][i] = i + nx_*(ny_-1);
    }
    for (int i=0; i<ny_; ++i) {
      (*meshSideSets_)[3][3][i] = i*nx_;
    }

  } // computeSideSets

  ROL::Ptr<std::vector<std::vector<std::vector<int>>>> getSideSets(
              const bool verbose = false,
              std::ostream & outStream = std::cout) const { 
    if (verbose) {
      outStream << "Mesh_stoch_adv_diff: getSideSets called" << std::endl;
      outStream << "Mesh_stoch_adv_diff: numSideSets = "     << meshSideSets_->size() << std::endl;
    }
    return meshSideSets_;
  }
  
}; // MeshManager_OED_Poisson
