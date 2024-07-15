// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "../TOOLS/meshmanager.hpp"

template <class Real>
class MeshManager_ThermalFluids : public MeshManager_Rectangle<Real> {
private:
  int nx_;
  int ny_;

  ROL::Ptr<std::vector<std::vector<std::vector<int> > > >  meshSideSets_;

public:

  MeshManager_ThermalFluids(Teuchos::ParameterList &parlist) : MeshManager_Rectangle<Real>(parlist) {
    nx_ = parlist.sublist("Geometry").get("NX", 3);
    ny_ = parlist.sublist("Geometry").get("NY", 1);
    computeSideSets();
  }

  void computeSideSets() {

    int numSideSets = 8;
    meshSideSets_ = ROL::makePtr<std::vector<std::vector<std::vector<int> > >>(numSideSets);

    Real patchFrac = static_cast<Real>(1)/static_cast<Real>(3);
    int np1 = static_cast<int>(patchFrac * static_cast<Real>(nx_));
    int np2 = static_cast<int>(patchFrac * static_cast<Real>(nx_));
    int np3 = nx_-(np1+np2);

    // Bottom
    (*meshSideSets_)[0].resize(4);
    (*meshSideSets_)[0][0].resize(nx_);
    (*meshSideSets_)[0][1].resize(0);
    (*meshSideSets_)[0][2].resize(0);
    (*meshSideSets_)[0][3].resize(0);
    // Right
    (*meshSideSets_)[1].resize(4);
    (*meshSideSets_)[1][0].resize(0);
    (*meshSideSets_)[1][1].resize(ny_);
    (*meshSideSets_)[1][2].resize(0);
    (*meshSideSets_)[1][3].resize(0);
    // Left
    (*meshSideSets_)[2].resize(4);
    (*meshSideSets_)[2][0].resize(0);
    (*meshSideSets_)[2][1].resize(0);
    (*meshSideSets_)[2][2].resize(0);
    (*meshSideSets_)[2][3].resize(ny_);
    // Top
    (*meshSideSets_)[3].resize(4);
    (*meshSideSets_)[3][0].resize(0);
    (*meshSideSets_)[3][1].resize(0);
    (*meshSideSets_)[3][2].resize(nx_);
    (*meshSideSets_)[3][3].resize(0);
    // Thermal boundaries
    // Top Right
    (*meshSideSets_)[4].resize(4);
    (*meshSideSets_)[4][0].resize(0);
    (*meshSideSets_)[4][1].resize(0);
    (*meshSideSets_)[4][2].resize(np3);
    (*meshSideSets_)[4][3].resize(0);
    // Top Center
    (*meshSideSets_)[5].resize(4);
    (*meshSideSets_)[5][0].resize(0);
    (*meshSideSets_)[5][1].resize(0);
    (*meshSideSets_)[5][2].resize(np2);
    (*meshSideSets_)[5][3].resize(0);
    // Top Left
    (*meshSideSets_)[6].resize(4);
    (*meshSideSets_)[6][0].resize(0);
    (*meshSideSets_)[6][1].resize(0);
    (*meshSideSets_)[6][2].resize(np1);
    (*meshSideSets_)[6][3].resize(0);
    // Pressure Pinning
    (*meshSideSets_)[7].resize(4);
    //(*meshSideSets_)[7][0].resize(1);
    //(*meshSideSets_)[7][1].resize(1);
    //(*meshSideSets_)[7][2].resize(1);
    //(*meshSideSets_)[7][3].resize(1);
    (*meshSideSets_)[7][0].resize(1);
    (*meshSideSets_)[7][1].resize(0);
    (*meshSideSets_)[7][2].resize(0);
    (*meshSideSets_)[7][3].resize(0);
    
    // Bottom
    for (int i=0; i<nx_; ++i) {
      (*meshSideSets_)[0][0][i] = i;
    }
    // Right
    for (int i=0; i<ny_; ++i) {
      (*meshSideSets_)[1][1][i] = (i+1)*nx_-1;
    }
    // Left
    for (int i=0; i<ny_; ++i) {
      (*meshSideSets_)[2][3][i] = i*nx_;
    }
    // Top
    for (int i=0; i<nx_; ++i) {
      (*meshSideSets_)[3][2][i] = i + nx_*(ny_-1);
    }
    // Top Left
    for (int i=0; i<np1; ++i) {
      (*meshSideSets_)[6][2][i] = i + nx_*(ny_-1);
    }
    // Top Center
    for (int i=0; i<np2; ++i) {
      (*meshSideSets_)[5][2][i] = i + np1 + nx_*(ny_-1);
    }
    // Top Right
    for (int i=0; i<np3; ++i) {
      (*meshSideSets_)[4][2][i] = i + (np1+np2) + nx_*(ny_-1);
    }
    // Pressure Pinning
    //int offset = nx_ * (ny_/2) + (nx_-1)/2;
    //(*meshSideSets_)[7][0][0] = offset + nx_ + 1;
    //(*meshSideSets_)[7][1][0] = offset + nx_;
    //(*meshSideSets_)[7][2][0] = offset;
    //(*meshSideSets_)[7][3][0] = offset + 1;
    (*meshSideSets_)[7][0][0] = 0;

  } // computeSideSets

  ROL::Ptr<std::vector<std::vector<std::vector<int> > > > getSideSets(
              const bool verbose = false,
              std::ostream & outStream = std::cout) const { 
    return meshSideSets_;
  }

}; // MeshManager_ThermalFluids
