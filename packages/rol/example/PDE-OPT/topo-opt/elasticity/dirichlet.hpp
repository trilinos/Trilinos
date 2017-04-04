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

/*! \file  dirichlet.hpp
    \brief Implements Dirichlet boundary conditions for the structural
           topology optimization problem.
*/

#ifndef ROL_PDEOPT_ELASTICITY_DIRICHLET_HPP
#define ROL_PDEOPT_ELASTICITY_DIRICHLET_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Intrepid_FieldContainer.hpp"
#include <vector>

template<class Real>
class Dirichlet {
private:
  std::vector<int> sidesets_, types_;
  std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
  std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  std::vector<std::vector<int> > fidx_;

  std::vector<int> getConstrainedDimensions(const int type, const int dim) const {
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
  Dirichlet(Teuchos::ParameterList &parlist) {
    // Grab sidesets
    Teuchos::Array<int> sidesets
      = Teuchos::getArrayFromStringParameter<int>(parlist.sublist("Dirichlet"), "Sidesets");
    sidesets_ = sidesets.toVector();
    // Grab type
    Teuchos::Array<int> types
      = Teuchos::getArrayFromStringParameter<int>(parlist.sublist("Dirichlet"), "Types");
    types_ = types.toVector();
  }

  void setCellNodes(const std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds,
                    const std::vector<std::vector<int> > fidx) {
    bdryCellNodes_  = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    fidx_           = fidx;
  }

  void applyResidual(std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > &R,
                     const std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > &U) const {
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
                (*R[conDim[m]])(cidx,fidx_[j][l]) = (*U[conDim[m]])(cidx,fidx_[j][l]);
              }
            }
          }
        }
      }
    }
  }

  void applyJacobian1(std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > &J) const {
    const int d = J.size();
    const int f = J[0][0]->dimension(1);
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
                    (*J[conDim[n]][p])(cidx,fidx_[j][l],m) = static_cast<Real>(0);
                  }
                  (*J[conDim[n]][conDim[n]])(cidx,fidx_[j][l],fidx_[j][l]) = static_cast<Real>(1);
                }
              }
            }
          }
        }
      }
    }
  }

  void applyJacobian2(std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > &J) const {
    const int d = J.size();
    const int f = J[0][0]->dimension(1);
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
                    (*J[conDim[n]][p])(cidx,fidx_[j][l],m) = static_cast<Real>(0);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  void applyMultiplier(std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > &L) const {
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
                (*L[conDim[m]])(cidx,fidx_[j][l]) = static_cast<Real>(0);
              }
            }
          }
        }
      }
    }
  }
};

#endif

//  Real dirichletFunc(const std::vector<Real> & coords, const int sideset, const int locSideId, const int dir) const {
//    Real val(0);
//    return val;
//  }
//
//  void computeDirichlet(void) {
//    // Compute Dirichlet values at DOFs.
//    int d = basisPtr_->getBaseCellTopology().getDimension();
//    int numSidesets = bdryCellLocIds_.size();
//    bdryCellDofValues_.resize(numSidesets);
//    for (int i=0; i<numSidesets; ++i) {
//      int numLocSides = bdryCellLocIds_[i].size();
//      bdryCellDofValues_[i].resize(numLocSides);
//      for (int j=0; j<numLocSides; ++j) {
//        int c = bdryCellLocIds_[i][j].size();
//        int f = basisPtr_->getCardinality();
//        bdryCellDofValues_[i][j] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, d));
//        Teuchos::RCP<Intrepid::FieldContainer<Real> > coords =
//          Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, d));
//        if (c > 0) {
//          fe_->computeDofCoords(coords, bdryCellNodes_[i][j]);
//        }
//        for (int k=0; k<c; ++k) {
//          for (int l=0; l<f; ++l) {
//            std::vector<Real> dofpoint(d);
//            //std::cout << "Sideset " << i << " LocalSide " << j << "  Cell " << k << "  Field " << l << "  Coord ";
//            for (int m=0; m<d; ++m) {
//              dofpoint[m] = (*coords)(k, l, m);
//              //std::cout << dofpoint[m] << "  ";
//            }
//
//            for (int m=0; m<d; ++m) {
//              (*bdryCellDofValues_[i][j])(k, l, m) = dirichletFunc(dofpoint, i, j, m);
//              //std::cout << "  " << m << "-Value " << dirichletFunc(dofpoint, i, j, m);
//            }
//            //std::cout << std::endl;
//          }
//        }
//      }
//    }
//  }
