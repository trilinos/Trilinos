// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BINARY_ELASTICITY_MESHMANAGER_HPP
#define BINARY_ELASTICITY_MESHMANAGER_HPP

#include "../../../TOOLS/meshmanager.hpp"

template <class Real>
class MeshManager_TopoOpt : public MeshManager_Rectangle<Real> {

private:

  int nx_;
  int ny_;
  Real height_;
  Real width_;
  std::string name_;
  ROL::Ptr<std::vector<std::vector<std::vector<int> > > >  meshSideSets_;

public: 

  MeshManager_TopoOpt(Teuchos::ParameterList &parlist) : MeshManager_Rectangle<Real>(parlist)
  {
    nx_     = parlist.sublist("Geometry").get("NX", 3);
    ny_     = parlist.sublist("Geometry").get("NY", 3);
    height_ = parlist.sublist("Geometry").get("Height", 1.0);
    width_  = parlist.sublist("Geometry").get("Width", 2.0);
    name_   = parlist.sublist("Problem").get("Example","Default");
    computeSideSets();
  }

  void computeSideSets() {
    if (name_ == "2D Cantilever with 3 Loads") {
      cantilever_computeSideSets();
    }
    else if (name_ == "2D Carrier Plate") {
      carrierplate_computeSideSets();
    }
    else if (name_ == "2D Wheel") {
      wheel_computeSideSets();
    }
    else {
      default_computeSideSets();
    }
  }

private:

  void cantilever_computeSideSets() {
    int numSideSets = 9;
    meshSideSets_ = ROL::makePtr<std::vector<std::vector<std::vector<int>>>>(numSideSets);

    Real pf(1.0/6.0);
    int np = static_cast<int>(pf * static_cast<Real>(nx_));
    int nf = nx_-5*np;

    // Bottom
    (*meshSideSets_)[0].resize(4);
    (*meshSideSets_)[0][0].resize(nf);
    (*meshSideSets_)[0][1].resize(0);
    (*meshSideSets_)[0][2].resize(0);
    (*meshSideSets_)[0][3].resize(0);
    (*meshSideSets_)[1].resize(4);
    (*meshSideSets_)[1][0].resize(np);
    (*meshSideSets_)[1][1].resize(0);
    (*meshSideSets_)[1][2].resize(0);
    (*meshSideSets_)[1][3].resize(0);
    (*meshSideSets_)[2].resize(4);
    (*meshSideSets_)[2][0].resize(np);
    (*meshSideSets_)[2][1].resize(0);
    (*meshSideSets_)[2][2].resize(0);
    (*meshSideSets_)[2][3].resize(0);
    (*meshSideSets_)[3].resize(4);
    (*meshSideSets_)[3][0].resize(np);
    (*meshSideSets_)[3][1].resize(0);
    (*meshSideSets_)[3][2].resize(0);
    (*meshSideSets_)[3][3].resize(0);
    (*meshSideSets_)[4].resize(4);
    (*meshSideSets_)[4][0].resize(np);
    (*meshSideSets_)[4][1].resize(0);
    (*meshSideSets_)[4][2].resize(0);
    (*meshSideSets_)[4][3].resize(0);
    (*meshSideSets_)[5].resize(4);
    (*meshSideSets_)[5][0].resize(np);
    (*meshSideSets_)[5][1].resize(0);
    (*meshSideSets_)[5][2].resize(0);
    (*meshSideSets_)[5][3].resize(0);
    // Right
    (*meshSideSets_)[6].resize(4);
    (*meshSideSets_)[6][0].resize(0);
    (*meshSideSets_)[6][1].resize(ny_);
    (*meshSideSets_)[6][2].resize(0);
    (*meshSideSets_)[6][3].resize(0);
    // Top
    (*meshSideSets_)[7].resize(4);
    (*meshSideSets_)[7][0].resize(0);
    (*meshSideSets_)[7][1].resize(0);
    (*meshSideSets_)[7][2].resize(nx_);
    (*meshSideSets_)[7][3].resize(0);
    // Left
    (*meshSideSets_)[8].resize(4);
    (*meshSideSets_)[8][0].resize(0);
    (*meshSideSets_)[8][1].resize(0);
    (*meshSideSets_)[8][2].resize(0);
    (*meshSideSets_)[8][3].resize(ny_);
    
    for (int i=0; i<nf; ++i) {
      (*meshSideSets_)[0][0][i] = i;
    }
    for (int i=0; i<np; ++i) {
      (*meshSideSets_)[1][0][i] = i + nf;
    }
    for (int i=0; i<np; ++i) {
      (*meshSideSets_)[2][0][i] = i + (nf+np);
    }
    for (int i=0; i<np; ++i) {
      (*meshSideSets_)[3][0][i] = i + (nf+2*np);
    }
    for (int i=0; i<np; ++i) {
      (*meshSideSets_)[4][0][i] = i + (nf+3*np);
    }
    for (int i=0; i<np; ++i) {
      (*meshSideSets_)[5][0][i] = i + (nf+4*np);
    }
    for (int i=0; i<ny_; ++i) {
      (*meshSideSets_)[6][1][i] = (i+1)*nx_-1;
    }
    for (int i=0; i<nx_; ++i) {
      (*meshSideSets_)[7][2][i] = i + nx_*(ny_-1);
    }
    for (int i=0; i<ny_; ++i) {
      (*meshSideSets_)[8][3][i] = i*nx_;
    }
  } // cantilever_computeSideSets

  void carrierplate_computeSideSets() {
    int numSideSets = 6;
    meshSideSets_ = ROL::makePtr<std::vector<std::vector<std::vector<int>>>>(numSideSets);

    Real pf(0.25);
    int np = static_cast<int>(pf * static_cast<Real>(nx_));
    int nc = nx_-2*np;

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
    // Top
    (*meshSideSets_)[2].resize(4);
    (*meshSideSets_)[2][0].resize(0);
    (*meshSideSets_)[2][1].resize(0);
    (*meshSideSets_)[2][2].resize(np);
    (*meshSideSets_)[2][3].resize(0);
    (*meshSideSets_)[3].resize(4);
    (*meshSideSets_)[3][0].resize(0);
    (*meshSideSets_)[3][1].resize(0);
    (*meshSideSets_)[3][2].resize(nc);
    (*meshSideSets_)[3][3].resize(0);
    (*meshSideSets_)[4].resize(4);
    (*meshSideSets_)[4][0].resize(0);
    (*meshSideSets_)[4][1].resize(0);
    (*meshSideSets_)[4][2].resize(np);
    (*meshSideSets_)[4][3].resize(0);
    // Left
    (*meshSideSets_)[5].resize(4);
    (*meshSideSets_)[5][0].resize(0);
    (*meshSideSets_)[5][1].resize(0);
    (*meshSideSets_)[5][2].resize(0);
    (*meshSideSets_)[5][3].resize(ny_);
    
    for (int i=0; i<nx_; ++i) {
      (*meshSideSets_)[0][0][i] = i;
    }
    for (int i=0; i<ny_; ++i) {
      (*meshSideSets_)[1][1][i] = (i+1)*nx_-1;
    }
    for (int i=0; i<np; ++i) {
      (*meshSideSets_)[2][2][i] = i + nx_*(ny_-1);
    }
    for (int i=0; i<nc; ++i) {
      (*meshSideSets_)[3][2][i] = i + nx_*(ny_-1) + np;
    }
    for (int i=0; i<np; ++i) {
      (*meshSideSets_)[4][2][i] = i + nx_*(ny_-1) + (np+nc);
    }
    for (int i=0; i<ny_; ++i) {
      (*meshSideSets_)[5][3][i] = i*nx_;
    }
  } // carrierplate_computeSideSets

  void wheel_computeSideSets() {
    int numSideSets = 8;
    meshSideSets_ = ROL::makePtr<std::vector<std::vector<std::vector<int>>>>(numSideSets);

    Real pf(0.125);
    Real patchFrac = (pf < width_) ? pf/width_ : pf;
    int np = static_cast<int>(patchFrac * static_cast<Real>(nx_));
    int nz = static_cast<int>(0.5*static_cast<Real>(nx_-4*np));
    int nc = nx_-2*(np+nz);

    // Bottom
    (*meshSideSets_)[0].resize(4);
    (*meshSideSets_)[0][0].resize(np);
    (*meshSideSets_)[0][1].resize(0);
    (*meshSideSets_)[0][2].resize(0);
    (*meshSideSets_)[0][3].resize(0);
    (*meshSideSets_)[1].resize(4);
    (*meshSideSets_)[1][0].resize(nz);
    (*meshSideSets_)[1][1].resize(0);
    (*meshSideSets_)[1][2].resize(0);
    (*meshSideSets_)[1][3].resize(0);
    (*meshSideSets_)[2].resize(4);
    (*meshSideSets_)[2][0].resize(nc);
    (*meshSideSets_)[2][1].resize(0);
    (*meshSideSets_)[2][2].resize(0);
    (*meshSideSets_)[2][3].resize(0);
    (*meshSideSets_)[3].resize(4);
    (*meshSideSets_)[3][0].resize(nz);
    (*meshSideSets_)[3][1].resize(0);
    (*meshSideSets_)[3][2].resize(0);
    (*meshSideSets_)[3][3].resize(0);
    (*meshSideSets_)[4].resize(4);
    (*meshSideSets_)[4][0].resize(np);
    (*meshSideSets_)[4][1].resize(0);
    (*meshSideSets_)[4][2].resize(0);
    (*meshSideSets_)[4][3].resize(0);
    // Right
    (*meshSideSets_)[5].resize(4);
    (*meshSideSets_)[5][0].resize(0);
    (*meshSideSets_)[5][1].resize(ny_);
    (*meshSideSets_)[5][2].resize(0);
    (*meshSideSets_)[5][3].resize(0);
    // Top
    (*meshSideSets_)[6].resize(4);
    (*meshSideSets_)[6][0].resize(0);
    (*meshSideSets_)[6][1].resize(0);
    (*meshSideSets_)[6][2].resize(nx_);
    (*meshSideSets_)[6][3].resize(0);
    // Left
    (*meshSideSets_)[7].resize(4);
    (*meshSideSets_)[7][0].resize(0);
    (*meshSideSets_)[7][1].resize(0);
    (*meshSideSets_)[7][2].resize(0);
    (*meshSideSets_)[7][3].resize(ny_);
    
    for (int i=0; i<np; ++i) {
      (*meshSideSets_)[0][0][i] = i;
    }
    for (int i=0; i<nz; ++i) {
      (*meshSideSets_)[1][0][i] = i + np;
    }
    for (int i=0; i<nc; ++i) {
      (*meshSideSets_)[2][0][i] = i + (np+nz);
    }
    for (int i=0; i<nz; ++i) {
      (*meshSideSets_)[3][0][i] = i + (np+nz+nc);
    }
    for (int i=0; i<np; ++i) {
      (*meshSideSets_)[4][0][i] = i + (np+nz+nc+nz);
    }
    for (int i=0; i<ny_; ++i) {
      (*meshSideSets_)[5][1][i] = (i+1)*nx_-1;
    }
    for (int i=0; i<nx_; ++i) {
      (*meshSideSets_)[6][2][i] = i + nx_*(ny_-1);
    }
    for (int i=0; i<ny_; ++i) {
      (*meshSideSets_)[7][3][i] = i*nx_;
    }
  } // wheel_computeSideSets

  void default_computeSideSets() {

    int numSideSets = 4;
    meshSideSets_ = ROL::makePtr<std::vector<std::vector<std::vector<int> > >>(numSideSets);

    //
    (*meshSideSets_)[0].resize(4);
    (*meshSideSets_)[0][0].resize(nx_);
    (*meshSideSets_)[0][1].resize(0);
    (*meshSideSets_)[0][2].resize(0);
    (*meshSideSets_)[0][3].resize(0);
    //
    (*meshSideSets_)[1].resize(4);
    (*meshSideSets_)[1][0].resize(0);
    (*meshSideSets_)[1][1].resize(ny_);
    (*meshSideSets_)[1][2].resize(0);
    (*meshSideSets_)[1][3].resize(0);
    //
    (*meshSideSets_)[2].resize(4);
    (*meshSideSets_)[2][0].resize(0);
    (*meshSideSets_)[2][1].resize(0);
    (*meshSideSets_)[2][2].resize(nx_);
    (*meshSideSets_)[2][3].resize(0);
    //
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
  } // default_computeSideSets

  ROL::Ptr<std::vector<std::vector<std::vector<int>>>> getSideSets(
              const bool verbose = false,
              std::ostream & outStream = std::cout) const { 
    if ( verbose ) {
      outStream << "Mesh_TopoOpt: getSideSets called" << std::endl;
      outStream << "Mesh_TopoOpt: numSideSets = " << meshSideSets_->size() << std::endl;
    }
    return meshSideSets_;
  }
  
}; // MeshManager_TopoOpt

#endif
