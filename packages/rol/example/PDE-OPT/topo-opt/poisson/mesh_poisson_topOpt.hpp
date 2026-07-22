// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "../../TOOLS/meshmanager.hpp"

template <class Real>
class MeshManager_Poisson_TopOpt : public MeshManager_Rectangle<Real> {

private:

  int nx_;
  int ny_;
  ROL::Ptr<std::vector<std::vector<std::vector<int> > > >  meshSideSets_;

public: 

  MeshManager_Poisson_TopOpt(Teuchos::ParameterList &parlist) : MeshManager_Rectangle<Real>(parlist)
  {
    nx_ = parlist.sublist("Geometry").get("NX", 3);
    ny_ = parlist.sublist("Geometry").get("NY", 3);
    computeSideSets();
  }


  void computeSideSets() {

    int numSideSets = 4;
    meshSideSets_ = ROL::makePtr<std::vector<std::vector<std::vector<int> > >>(numSideSets);

    int ny3 = ny_ / 3;
    int nyr = ny_ % 3;

    int n1 = (nyr == 2) ? ny3+1 : ny3;
    int n2 = (nyr == 1) ? ny3+1 : ny3;
    int n3 = (nyr == 2) ? ny3+1 : ny3;

    // Neumann sides
    (*meshSideSets_)[0].resize(4);
    (*meshSideSets_)[0][0].resize(nx_); // x = (0,1), y = {0}
    (*meshSideSets_)[0][1].resize(ny_); // x = {1}, y = (0,1)
    (*meshSideSets_)[0][2].resize(nx_); // x = (0,1), y = {1}
    (*meshSideSets_)[0][3].resize(0);
    // Neumann side
    (*meshSideSets_)[1].resize(4);
    (*meshSideSets_)[1][0].resize(0);
    (*meshSideSets_)[1][1].resize(0);
    (*meshSideSets_)[1][2].resize(0);
    (*meshSideSets_)[1][3].resize(n1);
    // Dirichlet side
    (*meshSideSets_)[2].resize(4);
    (*meshSideSets_)[2][0].resize(0);
    (*meshSideSets_)[2][1].resize(0);
    (*meshSideSets_)[2][2].resize(0);
    (*meshSideSets_)[2][3].resize(n2);
    // Neuman side
    (*meshSideSets_)[3].resize(4);
    (*meshSideSets_)[3][0].resize(0);
    (*meshSideSets_)[3][1].resize(0);
    (*meshSideSets_)[3][2].resize(0);
    (*meshSideSets_)[3][3].resize(n3);
    
    for (int i=0; i<nx_; ++i) {
      (*meshSideSets_)[0][0][i] = i;
    }
    for (int i=0; i<ny_; ++i) {
      (*meshSideSets_)[0][1][i] = (i+1)*nx_-1;
    }
    for (int i=0; i<nx_; ++i) {
      (*meshSideSets_)[0][2][i] = i + nx_*(ny_-1);
    }
    for (int i=0; i<ny_; ++i) {
      if ( i < n1 ) {
        (*meshSideSets_)[1][3][i] = i*nx_;
      }
      if ( i >= n1 && i < n1+n2 ) {
        (*meshSideSets_)[2][3][i-n1] = i*nx_;
      }
      if ( i >= n1+n2 ) {
        (*meshSideSets_)[3][3][i-n1-n2] = i*nx_;
      }
    }

  } // computeSideSets

  ROL::Ptr<std::vector<std::vector<std::vector<int> > > > getSideSets(
              const bool verbose = false,
              std::ostream & outStream = std::cout) const { 
    if ( verbose ) {
      outStream << "Mesh_Poisson_TopOpt: getSideSets called" << std::endl;
      outStream << "Mesh_Poisson_TopOpt: numSideSets = " << meshSideSets_->size() << std::endl;
    }
    return meshSideSets_;
  }
  
}; // MeshManager_TopoOpt
