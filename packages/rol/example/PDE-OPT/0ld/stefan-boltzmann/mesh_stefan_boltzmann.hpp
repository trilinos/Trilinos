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
class MeshManager_Stefan_Boltzmann : public MeshManager_Rectangle<Real> {

private:

  int nx_;
  int ny_;
  ROL::Ptr<std::vector<std::vector<Intrepid::FieldContainer<int> > > >  meshSideSets_;

public: 

  MeshManager_Stefan_Boltzmann(Teuchos::ParameterList &parlist) : MeshManager_Rectangle<Real>(parlist)
  {
    nx_ = parlist.sublist("Geometry").get("NX", 3);
    ny_ = parlist.sublist("Geometry").get("NY", 3);
    computeSideSets();
  }


  void computeSideSets() {

    int numSideSets = 3;
    meshSideSets_ = ROL::makePtr<std::vector<std::vector<Intrepid::FieldContainer<int> > >>(numSideSets);

    // Dirichlet
    (*meshSideSets_)[0].resize(2);
    (*meshSideSets_)[0][0].resize(nx_);
    (*meshSideSets_)[0][1].resize(ny_);
    // Robin
    (*meshSideSets_)[1].resize(1);
    (*meshSideSets_)[1][0].resize(ny_);
    // Neumann
    (*meshSideSets_)[2].resize(1);
    (*meshSideSets_)[2][0].resize(nx_);
    
    for (int i=0; i<nx_; ++i) {
      (*meshSideSets_)[0][0](i) = i;
    }
    for (int i=0; i<ny_; ++i) {
      (*meshSideSets_)[1][0](i) = (i+1)*nx_-1;
    }
    for (int i=0; i<nx_; ++i) {
      (*meshSideSets_)[2][0](i) = i + nx_*(ny_-1);
    }
    for (int i=0; i<ny_; ++i) {
      (*meshSideSets_)[0][1](i) = i*nx_;
    }

  } // computeSideSets

  ROL::Ptr<std::vector<std::vector<Intrepid::FieldContainer<int> > > > getSideSets(
      std::ostream & outStream = std::cout,
      const bool verbose = false) const {
    if ( verbose ) {
      outStream << "Mesh_SB: getSideSets called" << std::endl;
      outStream << "Mesh_SB: numSideSets = " << meshSideSets_->size() << std::endl;
    }
    return meshSideSets_;
  }
  
}; // MeshManager_Stefan_Boltzmann
