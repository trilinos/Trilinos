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

#include "../TOOLS/meshmanager.hpp"

template <class Real>
class MeshManager_Example02 : public MeshManager_Rectangle<Real> {

private:

  int nx_;
  int ny_;
  Real width_;
  ROL::Ptr<std::vector<std::vector<std::vector<int> > > >  meshSideSets_;

public: 

  MeshManager_Example02(Teuchos::ParameterList &parlist) : MeshManager_Rectangle<Real>(parlist)
  {
    nx_ = parlist.sublist("Geometry").get("NX", 3);
    ny_ = parlist.sublist("Geometry").get("NY", 3);
    width_ = parlist.sublist("Geometry").get("Width", 1.0);
    computeSideSets();
  }


  void computeSideSets() {

    int numSideSets = 4;
    meshSideSets_ = ROL::makePtr<std::vector<std::vector<std::vector<int> > >>(numSideSets);

    Real patchFrac = static_cast<Real>(1)/static_cast<Real>(6);
    int np1 = static_cast<int>(patchFrac * static_cast<Real>(nx_));   // Source
    int np2 = static_cast<int>(patchFrac * static_cast<Real>(nx_));
    int np3 = 2*static_cast<int>(patchFrac * static_cast<Real>(nx_)); // Gate
    int np4 = static_cast<int>(patchFrac * static_cast<Real>(nx_));
    int np5 = nx_-(np1+np2+np3+np4);                                  // Drain

    // Neumann
    (*meshSideSets_)[0].resize(4);
    (*meshSideSets_)[0][0].resize(nx_);
    (*meshSideSets_)[0][1].resize(ny_);
    (*meshSideSets_)[0][2].resize(np2+np4);
    (*meshSideSets_)[0][3].resize(ny_);
    // Source
    (*meshSideSets_)[1].resize(4);
    (*meshSideSets_)[1][0].resize(0);
    (*meshSideSets_)[1][1].resize(0);
    (*meshSideSets_)[1][2].resize(np1);
    (*meshSideSets_)[1][3].resize(0);
    // Gate
    (*meshSideSets_)[2].resize(4);
    (*meshSideSets_)[2][0].resize(0);
    (*meshSideSets_)[2][1].resize(0);
    (*meshSideSets_)[2][2].resize(np3);
    (*meshSideSets_)[2][3].resize(0);
    // Drain
    (*meshSideSets_)[3].resize(4);
    (*meshSideSets_)[3][0].resize(0);
    (*meshSideSets_)[3][1].resize(0);
    (*meshSideSets_)[3][2].resize(np5);
    (*meshSideSets_)[3][3].resize(0);
    
    for (int i=0; i<nx_; ++i) {
      (*meshSideSets_)[0][0][i] = i;
    }
    for (int i=0; i<ny_; ++i) {
      (*meshSideSets_)[0][1][i] = (i+1)*nx_-1;
    }
    for (int i=0; i<ny_; ++i) {
      (*meshSideSets_)[0][3][i] = i*nx_;
    }
    for (int i=0; i<np1; ++i) {
      (*meshSideSets_)[1][2][i] = i + nx_*(ny_-1);
    }
    for (int i=0; i<np2; ++i) {
      (*meshSideSets_)[0][2][i] = i + np1 + nx_*(ny_-1);
    }
    for (int i=0; i<np3; ++i) {
      (*meshSideSets_)[2][2][i] = i + (np1+np2) + nx_*(ny_-1);
    }
    for (int i=0; i<np4; ++i) {
      (*meshSideSets_)[0][2][i] = i + (np1+np2+np3) + nx_*(ny_-1);
    }
    for (int i=0; i<np5; ++i) {
      (*meshSideSets_)[3][2][i] = i + (np1+np2+np3+np4) + nx_*(ny_-1);
    }

  } // computeSideSets

  ROL::Ptr<std::vector<std::vector<std::vector<int> > > > getSideSets(
              const bool verbose = false,
              std::ostream & outStream = std::cout) const { 
    if ( verbose ) {
      outStream << "Mesh_TopoOpt: getSideSets called" << std::endl;
      outStream << "Mesh_TopoOpt: numSideSets = " << meshSideSets_->size() << std::endl;
    }
    return meshSideSets_;
  }
  
}; // MeshManager_TopoOpt
