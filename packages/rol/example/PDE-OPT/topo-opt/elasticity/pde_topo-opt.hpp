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

/*! \file  pde_topo-opt.hpp
    \brief Implements the local PDE interface for the structural topology
           optimization problem.
*/

#ifndef PDE_TOPOOPT_HPP
#define PDE_TOPOOPT_HPP

#include "../../TOOLS/pde.hpp"
#include "../../TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "Teuchos_RCP.hpp"



template<class Real>
class Load {
private:
  Real loadMagnitude_;
  std::vector<Real> loadAngle_;
  std::vector<Real> loadLocation_;
  std::vector<Real> loadWidth_;

  Real DegreesToRadians(const Real deg) const {
    return deg * static_cast<Real>(M_PI) / static_cast<Real>(180);
  }

public:
  Load(Teuchos::ParameterList & parlist) {
    loadMagnitude_ = parlist.get("Load Magnitude", 1.0);
    Teuchos::Array<Real> angle
      = Teuchos::getArrayFromStringParameter<double>(parlist, "Load Angle");
    loadAngle_ = angle.toVector();
    Teuchos::Array<Real> location
      = Teuchos::getArrayFromStringParameter<double>(parlist, "Load Location");
    loadLocation_ = location.toVector();
    Teuchos::Array<Real> width
      = Teuchos::getArrayFromStringParameter<double>(parlist, "Load Width");
    loadWidth_ = width.toVector();

    int size = loadAngle_.size();
    for (int i=0; i<size; ++i) {
      loadAngle_[i] = DegreesToRadians(loadAngle_[i]);
    }
  }

  Real loadFunc(const std::vector<Real> & coords, const int dir, const std::vector<Real> & param) const {
    Real loadMagNoise(0), loadAngNoise0(0), loadAngNoise1(0);
    if (param.size() > 0) {
      loadMagNoise  = param[0];
    }
    if (param.size() > 1) {
      loadAngNoise0 = DegreesToRadians(param[1]);
    }
    if (param.size() > 2) {
      loadAngNoise1 = DegreesToRadians(param[2]);
    }
    const Real half(0.5);
    const Real loadMagnitude = loadMagnitude_ + loadMagNoise;
    const Real loadAngle0    = loadAngle_[0] + loadAngNoise0;
    const Real Gx = std::exp(-half*std::pow(coords[0]-loadLocation_[0],2)/std::pow(loadWidth_[0],2));
    const Real Gy = std::exp(-half*std::pow(coords[1]-loadLocation_[1],2)/std::pow(loadWidth_[1],2));

    Real val=0;
    int d = coords.size();
    if (d==2) {
      if (dir==0) {
        val = loadMagnitude*std::cos(loadAngle0)*Gx*Gy;
      }
      if (dir==1) {
        val = loadMagnitude*std::sin(loadAngle0)*Gx*Gy;
      }
    }
    if (d==3) {
      const Real loadAngle1 = loadAngle_[1] + loadAngNoise1;
      const Real Gz = std::exp(-half*std::pow(coords[2]-loadLocation_[2],2)/std::pow(loadWidth_[2],2));
      if (dir==0) {
        val = loadMagnitude*std::sin(loadAngle0)*std::cos(loadAngle1)*Gx*Gy*Gz;
      }
      if (dir==1) {
        val = loadMagnitude*std::sin(loadAngle0)*std::sin(loadAngle1)*Gx*Gy*Gz;
      }
      if (dir==2) {
        val = loadMagnitude*std::cos(loadAngle0)*Gx*Gy*Gz;
      }
    }
    return val;
  }

  void compute(std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > &load, const Teuchos::RCP<FE<Real> > & fe, const std::vector<Real> & param, const Real scale=1) const {
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


template<class Real>
class FieldHelper {
  private:
  const int numFields_, numDofs_;
  const std::vector<int> numFieldDofs_;
  const std::vector<std::vector<int> > fieldPattern_;

  public:
  FieldHelper(const int numFields, const int numDofs,
              const std::vector<int> &numFieldDofs,
              const std::vector<std::vector<int> > &fieldPattern)
    : numFields_(numFields), numDofs_(numDofs),
      numFieldDofs_(numFieldDofs), fieldPattern_(fieldPattern) {}

  void splitFieldCoeff(std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > & U,
                       const Teuchos::RCP<const Intrepid::FieldContainer<Real> >   & u_coeff) const {
    U.resize(numFields_);
    int  c = u_coeff->dimension(0);
    for (int i=0; i<numFields_; ++i) {
      U[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,numFieldDofs_[i]));
      for (int j=0; j<c; ++j) {
        for (int k=0; k<numFieldDofs_[i]; ++k) {
          //U[i](j,k) = u_coeff(j,offset[i]+k);
          (*U[i])(j,k) = (*u_coeff)(j,fieldPattern_[i][k]);
        }
      }
    }
  }

  void combineFieldCoeff(Teuchos::RCP<Intrepid::FieldContainer<Real> >   & res,
                         const std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > & R) const {
    int c = R[0]->dimension(0);  // number of cells
    res = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, numDofs_));
    for (int i=0; i<numFields_; ++i) {
      for (int j=0; j<c; ++j) {
        for (int k=0; k<numFieldDofs_[i]; ++k) {
          (*res)(j,fieldPattern_[i][k]) = (*R[i])(j,k);
        }
      }
    }
  }

  void combineFieldCoeff(Teuchos::RCP<Intrepid::FieldContainer<Real> >   & jac,
                         const std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > & J) const {
    int c = J[0][0]->dimension(0);  // number of cells
    jac = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, numDofs_, numDofs_));        
    for (int i=0; i<numFields_; ++i) {
      for (int j=0; j<numFields_; ++j) {
        for (int k=0; k<c; ++k) {
          for (int l=0; l<numFieldDofs_[i]; ++l) {
            for (int m=0; m<numFieldDofs_[j]; ++m) {
              (*jac)(k,fieldPattern_[i][l],fieldPattern_[j][m]) = (*J[i][j])(k,l,m);
            }
          }
        }
      }
    }
  }

  int numFields(void) const {
    return numFields_;
  }

};

template <class Real>
class PDE_TopoOpt : public PDE<Real> {
private:
  // Finite element basis information
  Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtr_;
  std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
  // Cell cubature information
  Teuchos::RCP<Intrepid::Cubature<Real> > cellCub_;
  Teuchos::RCP<Intrepid::Cubature<Real> > bdryCub_;
  // Cell node information
  Teuchos::RCP<Intrepid::FieldContainer<Real> > volCellNodes_;
  std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;
  std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  // Finite element definition
  Teuchos::RCP<FE<Real> > fe_;
  std::vector<Teuchos::RCP<FE<Real> > > feBdry_;
  // Local degrees of freedom on boundary, for each side of the reference cell (first index).
  std::vector<std::vector<int> > fidx_;
  // Coordinates of degrees freedom on boundary cells.
  // Indexing:  [sideset number][local side id](cell number, value at dof)
  std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > bdryCellDofValues_;
  // Field pattern, offsets, etc.
  std::vector<std::vector<int> > fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                                // number of fields (equations in the PDE)
  int numDofs_;                                  // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                      // for each field, a counting offset
  std::vector<int> numFieldDofs_;                // for each field, number of degrees of freedom
  
  // Problem parameters.
  Real youngsModulus_;
  Real poissonRatio_;
  bool isPlainStress_;
  Real minDensity_;
  Real maxDensity_;
  Real powerSIMP_;
  Teuchos::RCP<Load<Real> > load_; 
  std::vector<std::vector<Real> > materialMat_;

  Teuchos::RCP<FieldHelper<Real> > fieldHelper_;

  // Precomputed quantities.
  std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > BMat_;
  std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > BdetJMat_;
  std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > NMat_;
  std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > NdetJMat_;
  std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > CBdetJMat_;

  Real dirichletFunc(const std::vector<Real> & coords, const int sideset, const int locSideId, const int dir) const {
    Real val(0);
    return val;
  }

  void computeDirichlet(void) {
    // Compute Dirichlet values at DOFs.
    int d = basisPtr_->getBaseCellTopology().getDimension();
    int numSidesets = bdryCellLocIds_.size();
    bdryCellDofValues_.resize(numSidesets);
    for (int i=0; i<numSidesets; ++i) {
      int numLocSides = bdryCellLocIds_[i].size();
      bdryCellDofValues_[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        int c = bdryCellLocIds_[i][j].size();
        int f = basisPtr_->getCardinality();
        bdryCellDofValues_[i][j] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, d));
        Teuchos::RCP<Intrepid::FieldContainer<Real> > coords =
          Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, d));
        if (c > 0) {
          fe_->computeDofCoords(coords, bdryCellNodes_[i][j]);
        }
        for (int k=0; k<c; ++k) {
          for (int l=0; l<f; ++l) {
            std::vector<Real> dofpoint(d);
            //std::cout << "Sideset " << i << " LocalSide " << j << "  Cell " << k << "  Field " << l << "  Coord ";
            for (int m=0; m<d; ++m) {
              dofpoint[m] = (*coords)(k, l, m);
              //std::cout << dofpoint[m] << "  ";
            }

            for (int m=0; m<d; ++m) {
              (*bdryCellDofValues_[i][j])(k, l, m) = dirichletFunc(dofpoint, i, j, m);
              //std::cout << "  " << m << "-Value " << dirichletFunc(dofpoint, i, j, m);
            }
            //std::cout << std::endl;
          }
        }
      }
    }
  }

  void computeMaterialTensor(const int d) {
    int matd = (d*(d+1))/2;
    std::vector<Real> tmpvec(matd);
    materialMat_.resize(matd, tmpvec);
    if ((d==2) && (isPlainStress_)) {
      Real one(1), half(0.5);
      Real factor1 = youngsModulus_ /(one-poissonRatio_*poissonRatio_);
      materialMat_[0][0] = factor1;
      materialMat_[0][1] = factor1 * poissonRatio_;
      materialMat_[1][0] = factor1 * poissonRatio_;
      materialMat_[1][1] = factor1;
      materialMat_[2][2] = factor1 * half * (one-poissonRatio_);
    }		
    else if ((d==2) && (!isPlainStress_)) {
      Real one(1), two(2), half(0.5);
      Real factor2 = youngsModulus_ /(one+poissonRatio_)/(one-two*poissonRatio_);
      materialMat_[0][0] = factor2 * (one-poissonRatio_);
      materialMat_[0][1] = factor2 * poissonRatio_;
      materialMat_[1][0] = factor2 * poissonRatio_;
      materialMat_[1][1] = factor2 * (one-poissonRatio_);
      materialMat_[2][2] = factor2 * half * (one-two*poissonRatio_);
    }
    else {
      Real one(1), two(2), half(0.5);
      Real lam = youngsModulus_*poissonRatio_/(one+poissonRatio_)/(one-two*poissonRatio_);
      Real mu = half*youngsModulus_/(one+poissonRatio_);
      materialMat_[0][0] = lam + two*mu;
      materialMat_[0][1] = lam;
      materialMat_[0][2] = lam;
      materialMat_[1][0] = lam;
      materialMat_[1][1] = lam + two*mu;
      materialMat_[1][2] = lam;
      materialMat_[2][0] = lam;
      materialMat_[2][1] = lam;
      materialMat_[2][2] = lam + two*mu;
      materialMat_[3][3] = mu;
      materialMat_[4][4] = mu;
      materialMat_[5][5] = mu;
    }				
  }

  void computeNBmats(void) {
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = materialMat_.size();

    for (int i=0; i<d; ++i) {
      BMat_.push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,p,matd)));
      BdetJMat_.push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,p,matd)));
      CBdetJMat_.push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,p,matd)));
      NMat_.push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,p,d)));
      NdetJMat_.push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,p,d)));
    }

    if (d==2) {
      for (int i=0; i<c; ++i) {
        for (int j=0; j<f; ++j) {
      	  for (int k=0; k<p; ++k) {
            (*NMat_[0])(i, j, k, 0) = (*fe_->N())(i, j, k);
            (*NMat_[1])(i, j, k, 1) = (*fe_->N())(i, j, k);
            (*NdetJMat_[0])(i, j, k, 0) = (*fe_->NdetJ())(i, j, k);
            (*NdetJMat_[1])(i, j, k, 1) = (*fe_->NdetJ())(i, j, k);
                
            (*BMat_[0])(i, j, k, 0) = (*fe_->gradN())(i, j, k, 0);
            (*BMat_[1])(i, j, k, 1) = (*fe_->gradN())(i, j, k, 1);
            (*BMat_[0])(i, j, k, 2) = (*fe_->gradN())(i, j, k, 1);
            (*BMat_[1])(i, j, k, 2) = (*fe_->gradN())(i, j, k, 0);

            (*BdetJMat_[0])(i, j, k, 0) = (*fe_->gradNdetJ())(i, j, k, 0);
            (*BdetJMat_[1])(i, j, k, 1) = (*fe_->gradNdetJ())(i, j, k, 1);
            (*BdetJMat_[0])(i, j, k, 2) = (*fe_->gradNdetJ())(i, j, k, 1);
            (*BdetJMat_[1])(i, j, k, 2) = (*fe_->gradNdetJ())(i, j, k, 0);
          }
        }
      }
    }

    if(d==3) {
      for (int i=0; i<c; ++i) {
        for (int j=0; j<f; ++j) {
          for (int k=0; k<p; ++k) {
            (*NMat_[0])(i, j, k, 0) = (*fe_->N())(i, j, k);
            (*NMat_[1])(i, j, k, 1) = (*fe_->N())(i, j, k);
            (*NMat_[2])(i, j, k, 2) = (*fe_->N())(i, j, k);
            (*NdetJMat_[0])(i, j, k, 0) = (*fe_->NdetJ())(i, j, k);
            (*NdetJMat_[1])(i, j, k, 1) = (*fe_->NdetJ())(i, j, k);
            (*NdetJMat_[2])(i, j, k, 2) = (*fe_->NdetJ())(i, j, k);
            
            (*BMat_[0])(i, j, k, 0) = (*fe_->gradN())(i, j, k, 0);
            (*BMat_[1])(i, j, k, 1) = (*fe_->gradN())(i, j, k, 1);
            (*BMat_[2])(i, j, k, 2) = (*fe_->gradN())(i, j, k, 2);
            (*BMat_[1])(i, j, k, 3) = (*fe_->gradN())(i, j, k, 2);
            (*BMat_[2])(i, j, k, 3) = (*fe_->gradN())(i, j, k, 1);
            (*BMat_[0])(i, j, k, 4) = (*fe_->gradN())(i, j, k, 2);
            (*BMat_[2])(i, j, k, 4) = (*fe_->gradN())(i, j, k, 0);
            (*BMat_[0])(i, j, k, 5) = (*fe_->gradN())(i, j, k, 1);
            (*BMat_[1])(i, j, k, 5) = (*fe_->gradN())(i, j, k, 0);
                
            (*BdetJMat_[0])(i, j, k, 0) = (*fe_->gradNdetJ())(i, j, k, 0);
            (*BdetJMat_[1])(i, j, k, 1) = (*fe_->gradNdetJ())(i, j, k, 1);
            (*BdetJMat_[2])(i, j, k, 2) = (*fe_->gradNdetJ())(i, j, k, 2);
            (*BdetJMat_[1])(i, j, k, 3) = (*fe_->gradNdetJ())(i, j, k, 2);
            (*BdetJMat_[2])(i, j, k, 3) = (*fe_->gradNdetJ())(i, j, k, 1);
            (*BdetJMat_[0])(i, j, k, 4) = (*fe_->gradNdetJ())(i, j, k, 2);
            (*BdetJMat_[2])(i, j, k, 4) = (*fe_->gradNdetJ())(i, j, k, 0);
            (*BdetJMat_[0])(i, j, k, 5) = (*fe_->gradNdetJ())(i, j, k, 1);
            (*BdetJMat_[1])(i, j, k, 5) = (*fe_->gradNdetJ())(i, j, k, 0);
          }
        }
      }
    } 

    for (int i=0; i<c; ++i) {
      for (int j=0; j<f; ++j) {
        for (int k=0; k<p; ++k) {
          for (int l=0; l<d; ++l) {
            for (int m=0; m<matd; ++m) {
              for (int n=0; n<matd; ++n) {
                (*CBdetJMat_[l])(i,j,k,m) += materialMat_[m][n] * (*BdetJMat_[l])(i,j,k,n);
              }
            }
          }
        }
      }
    }

  }

  void applyTensor(Teuchos::RCP<Intrepid::FieldContainer<Real> > & out, const Teuchos::RCP<Intrepid::FieldContainer<Real> > & inData) const {
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int matd = materialMat_.size();

    for (int i=0; i<c; ++i) {
      for (int j=0; j<p; ++j) {
        for (int k=0; k<matd; ++k) {
          for (int l=0; l<matd; ++l) {
            (*out)(i,j,k) += materialMat_[k][l] * (*inData)(i,j,l);
          }
        }
      }
    }
  }

  void computeUmat(Teuchos::RCP<Intrepid::FieldContainer<Real> > & UMat, std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > & gradU) const {
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = materialMat_.size();

    UMat = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,p,matd));

    if (d==2) {
      for (int i=0; i<c; ++i) {
        for (int j=0; j<p; ++j) {
          (*UMat)(i, j, 0) = (*gradU[0])(i, j, 0);
          (*UMat)(i, j, 1) = (*gradU[1])(i, j, 1);
          (*UMat)(i, j, 2) = (*gradU[0])(i, j, 1) + (*gradU[1])(i, j, 0);
        }
      }
    }

    if(d==3) {
      for (int i=0; i<c; ++i) {
        for (int j=0; j<p; ++j) {
          (*UMat)(i, j, 0) = (*gradU[0])(i, j, 0);
          (*UMat)(i, j, 1) = (*gradU[1])(i, j, 1);
          (*UMat)(i, j, 2) = (*gradU[2])(i, j, 2);
          (*UMat)(i, j, 3) = (*gradU[1])(i, j, 2) + (*gradU[2])(i, j, 1);
          (*UMat)(i, j, 4) = (*gradU[0])(i, j, 2) + (*gradU[2])(i, j, 0);
          (*UMat)(i, j, 5) = (*gradU[0])(i, j, 1) + (*gradU[1])(i, j, 0);
        }
      }
    } 

  }

  void computeDensity(Teuchos::RCP<Intrepid::FieldContainer<Real> > & rho,
                      const Teuchos::RCP<Intrepid::FieldContainer<Real> > & Z,
                      const int deriv = 0) const {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);

    for (int i=0; i<c; ++i) {
      for (int j=0; j<p; ++j) {
        if (deriv==0) {
          (*rho)(i,j) = minDensity_ + (maxDensity_-minDensity_)*std::pow((*Z)(i,j), powerSIMP_);
        }
        else if (deriv==1) { 
          (*rho)(i,j) = powerSIMP_*(maxDensity_-minDensity_)*std::pow((*Z)(i,j), powerSIMP_-1);
        }
        else if	(deriv==2) {
          (*rho)(i,j) = powerSIMP_*(powerSIMP_-1)*(maxDensity_-minDensity_)*std::pow((*Z)(i,j), powerSIMP_-2);
        } 
      }
    }
    
  }

public:
  PDE_TopoOpt(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    basisPtr_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();            // get the cell type from the basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                           // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 2);        // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                           // create default cubature
    int d = cellType.getDimension();

    basisPtrs_.clear();
    for (int i=0; i<d; ++i) {
      basisPtrs_.push_back(basisPtr_);  // Displacement component
    }

    shards::CellTopology bdryCellType = cellType.getCellTopologyData(d-1, 0);
    int bdryCubDegree = parlist.sublist("Problem").get("Boundary Cubature Degree",2); // set cubature degree, e.g., 2
    bdryCub_ = cubFactory.create(bdryCellType, bdryCubDegree);

    // Other problem parameters.
    youngsModulus_  = parlist.sublist("Problem").get("Young's Modulus",     1.0);
    poissonRatio_   = parlist.sublist("Problem").get("Poisson Ratio",       0.3);
    isPlainStress_  = parlist.sublist("Problem").get("Use Plain Stress",    true);
    minDensity_     = parlist.sublist("Problem").get("Minimum Density",     1e-4);
    maxDensity_     = parlist.sublist("Problem").get("Maximum Density",     1.0);
    powerSIMP_      = parlist.sublist("Problem").get("SIMP Power",          3.0);
    computeMaterialTensor(d);

    load_  = Teuchos::rcp(new Load<Real>(parlist.sublist("Problem")));

    numDofs_ = 0;
    numFields_ = basisPtrs_.size();
    offset_.resize(numFields_);
    numFieldDofs_.resize(numFields_);
    for (int i=0; i<numFields_; ++i) {
      if (i==0) {
        offset_[i]  = 0;
      }
      else {
        offset_[i]  = offset_[i-1] + basisPtrs_[i-1]->getCardinality();
      }
      numFieldDofs_[i] = basisPtrs_[i]->getCardinality();
      numDofs_ += numFieldDofs_[i];
    }
  }

  void residual(Teuchos::RCP<Intrepid::FieldContainer<Real> > & res,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = materialMat_.size();
 
    // Initialize residuals.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > R(d);
    for (int i=0; i<d; ++i) {
      R[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f));
    }

    // Split u_coeff into components.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U;
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Evaluate/interpolate finite element fields on cells.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > rho =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > UMat;
    Teuchos::RCP<Intrepid::FieldContainer<Real> > rhoUMat =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, matd));
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > gradDisp_eval(d);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > load(d);
    for (int i=0; i<d; ++i) {
      gradDisp_eval[i] =  Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
      load[i] =  Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    }
    for (int i=0; i<d; ++i) {
      fe_->evaluateGradient(gradDisp_eval[i], U[i]);
    }
    fe_->evaluateValue(valZ_eval, Z[0]);
    computeUmat(UMat, gradDisp_eval);
    computeDensity(rho, valZ_eval);
    load_->compute(load, fe_, PDE<Real>::getParameter(), -1.0);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoUMat, *rho, *UMat);

    /*** Evaluate weak form of the residual. ***/
    for (int i=0; i<d; ++i) {
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *rhoUMat,           // rho B U
                                                    *CBdetJMat_[i],     // B' C
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *load[i],           // F
                                                    *fe_->NdetJ(),      // N
                                                    Intrepid::COMP_CPP,
                                                    true);
    }

    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      computeDirichlet();
      for (int i = 0; i < numSideSets; ++i) {
        if ((i==3)) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                //std::cout << "\n  i=" << i << "   cidx=" << cidx << "   j=" << j << "  l=" << l << "  " << fidx_[j][l] << " " << (*bdryCellDofValues_[i][j])(k,fidx_[j][l],0);
                //std::cout << "\n  i=" << i << "   cidx=" << cidx << "   j=" << j << "  l=" << l << "  " << fidx_[j][l] << " " << (*bdryCellDofValues_[i][j])(k,fidx_[j][l],1);
                for (int m=0; m < d; ++m) {
                  (*R[m])(cidx,fidx_[j][l]) = (*U[m])(cidx,fidx_[j][l]) - (*bdryCellDofValues_[i][j])(k,fidx_[j][l],m);
                }
              }
            }
          }
        }
      }
    }

    // Combine the residuals.
    fieldHelper_->combineFieldCoeff(res, R);
  }

  void Jacobian_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = materialMat_.size();
 
    // Initialize Jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j) {
        J[i].push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,f)));
      }
    }

    // Split z_coeff into components.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    Teuchos::RCP<Intrepid::FieldContainer<Real> > rho =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > rhoBMat(d);
    for (int i=0; i<d; ++i) {
      rhoBMat[i] =  Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, p, matd));
    }
    fe_->evaluateValue(valZ_eval, Z[0]);
    computeDensity(rho, valZ_eval);
    for (int i=0; i<d; ++i) {
      Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*rhoBMat[i], *rho, *BMat_[i]);
    }

    /*** Evaluate weak form of the Jacobian. ***/
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j) {
        Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][j],
                                                      *rhoBMat[i],          // rho B
                                                      *CBdetJMat_[j],       // B' C
                                                      Intrepid::COMP_CPP,
                                                      false);
      }
    }

    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if ((i==3)) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < f; ++m) {
                  for (int n=0; n < d; ++n) {
                    for (int p=0; p < d; ++p) {
                      (*J[n][p])(cidx,fidx_[j][l],m) = static_cast<Real>(0);
                    }
                    (*J[n][n])(cidx,fidx_[j][l],fidx_[j][l]) = static_cast<Real>(1);
                  }
                }
              }
            }
          }
        }
      }
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);

  }


  void Jacobian_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = materialMat_.size();

    // Initialize Jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j) {
        J[i].push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,f)));
      }
    }

    // Split u_coeff into components.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U;
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Evaluate/interpolate finite element fields on cells.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > rho =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > UMat;
    Teuchos::RCP<Intrepid::FieldContainer<Real> > rhoUMat =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, matd));
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > gradDisp_eval(d);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > CBrhoUMat(d);
    for (int i=0; i<d; ++i) {
      gradDisp_eval[i] =  Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
      CBrhoUMat[i] =  Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, p));
    }
    for (int i=0; i<d; ++i) {
      fe_->evaluateGradient(gradDisp_eval[i], U[i]);
    }
    fe_->evaluateValue(valZ_eval, Z[0]);
    computeUmat(UMat, gradDisp_eval);
    computeDensity(rho, valZ_eval, 1);  // first derivative
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoUMat, *rho, *UMat);

    for (int i=0; i<d; ++i) {
      Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*CBrhoUMat[i], *rhoUMat, *CBdetJMat_[i]);
    }

    /*** Evaluate weak form of the residual. ***/
    for (int i=0; i<d; ++i) {
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][0],
                                                    *CBrhoUMat[i],      // B' C drho B U
                                                    *fe_->N(),          // N
                                                    Intrepid::COMP_CPP,
                                                    false);
    }

    // APPLY DIRICHLET CONDITIONS
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if ((i==3)) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                for (int m=0; m < f; ++m) {
                  for (int n=0; n < d; ++n) {
                    for (int p=0; p < d; ++p) {
                      (*J[n][p])(cidx,fidx_[j][l],m) = static_cast<Real>(0);
                    }
                  }
                }
              }
            }
          }
	}
      }
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);

  }

  void Hessian_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_TopoOpt::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = materialMat_.size();

    // Initialize Hessians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j) {
        J[i].push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,f)));
      }
    }

    // Split u_coeff into components.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > L;
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(L, l_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Apply Dirichlet conditions to the multipliers.
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if ((i==3)) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
                for (int m=0; m < d; ++m) {
                  (*L[m])(cidx,fidx_[j][l]) = static_cast<Real>(0);
                }
              }
            }
          }
        }
      }
    }

    // Evaluate/interpolate finite element fields on cells.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > rho =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > LMat;
    Teuchos::RCP<Intrepid::FieldContainer<Real> > rhoLMat =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, matd));
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > gradDisp_eval(d);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > CBrhoLMat(d);
    for (int i=0; i<d; ++i) {
      gradDisp_eval[i] =  Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
      CBrhoLMat[i] =  Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, p));
    }
    for (int i=0; i<d; ++i) {
      fe_->evaluateGradient(gradDisp_eval[i], L[i]);
    }
    fe_->evaluateValue(valZ_eval, Z[0]);
    computeUmat(LMat, gradDisp_eval);
    computeDensity(rho, valZ_eval, 1);  // first derivative
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoLMat, *rho, *LMat);

    for (int i=0; i<d; ++i) {
      Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*CBrhoLMat[i], *rhoLMat, *CBdetJMat_[i]);
    }

    /*** Evaluate weak form of the residual. ***/
    for (int i=0; i<d; ++i) {
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][i],
                                                    *fe_->N(),          // N
                                                    *CBrhoLMat[i],      // B' C drho B U
                                                    Intrepid::COMP_CPP,
                                                    false);
    }

    // Combine the Hessians.
    fieldHelper_->combineFieldCoeff(hess, J);

  }

  void Hessian_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = materialMat_.size();

    // Initialize Hessians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j) {
        J[i].push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,f)));
      }
    }

    // Split u_coeff into components.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > L;
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(L, l_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Apply Dirichlet conditions to the multipliers.
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if ((i==3)) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
                for (int m=0; m < d; ++m) {
                  (*L[m])(cidx,fidx_[j][l]) = static_cast<Real>(0);
                }
              }
            }
          }
        }
      }
    }

    // Evaluate/interpolate finite element fields on cells.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > rho =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > LMat;
    Teuchos::RCP<Intrepid::FieldContainer<Real> > rhoLMat =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, matd));
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > gradDisp_eval(d);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > CBrhoLMat(d);
    for (int i=0; i<d; ++i) {
      gradDisp_eval[i] =  Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
      CBrhoLMat[i] =  Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, p));
    }
    for (int i=0; i<d; ++i) {
      fe_->evaluateGradient(gradDisp_eval[i], L[i]);
    }
    fe_->evaluateValue(valZ_eval, Z[0]);
    computeUmat(LMat, gradDisp_eval);
    computeDensity(rho, valZ_eval, 1);  // first derivative
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoLMat, *rho, *LMat);

    for (int i=0; i<d; ++i) {
      Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*CBrhoLMat[i], *rhoLMat, *CBdetJMat_[i]);
    }

    /*** Evaluate weak form of the residual. ***/
    for (int i=0; i<d; ++i) {
      Intrepid::FunctionSpaceTools::integrate<Real>(*J[i][0],
                                                    *CBrhoLMat[i],      // B' C drho B U
                                                    *fe_->N(),          // N
                                                    Intrepid::COMP_CPP,
                                                    false);
    }

    // Combine the Hessians.
    fieldHelper_->combineFieldCoeff(hess, J);

  }

  void Hessian_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    int matd = materialMat_.size();

    // Initialize Hessians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j) {
        J[i].push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,f)));
      }
    }

    // Split u_coeff into components.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U;
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > L;
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(L, l_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Apply Dirichlet conditions to the multipliers.
    int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if ((i==3)) {
          int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            int numCellsSide = bdryCellLocIds_[i][j].size();
            int numBdryDofs = fidx_[j].size();
            for (int k = 0; k < numCellsSide; ++k) {
              int cidx = bdryCellLocIds_[i][j][k];
              for (int l = 0; l < numBdryDofs; ++l) {
                //std::cout << "\n   j=" << j << "  l=" << l << "  " << fidx[j][l];
                for (int m=0; m < d; ++m) {
                  (*L[m])(cidx,fidx_[j][l]) = static_cast<Real>(0);
                }
              }
            }
          }
        }
      }
    }

    // Evaluate/interpolate finite element fields on cells.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > rho =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > UMat;
    Teuchos::RCP<Intrepid::FieldContainer<Real> > LMat;
    Teuchos::RCP<Intrepid::FieldContainer<Real> > rhoLMat =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, matd));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > CUMat =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, matd));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > CUrhoLMat =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > NCUrhoLMat =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f, p));
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > gradDispU_eval(d);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > gradDispL_eval(d);
    for (int i=0; i<d; ++i) {
      gradDispU_eval[i] =  Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
      gradDispL_eval[i] =  Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    }
    for (int i=0; i<d; ++i) {
      fe_->evaluateGradient(gradDispU_eval[i], U[i]);
      fe_->evaluateGradient(gradDispL_eval[i], L[i]);
    }
    fe_->evaluateValue(valZ_eval, Z[0]);
    computeUmat(UMat, gradDispU_eval);
    computeUmat(LMat, gradDispL_eval);
    computeDensity(rho, valZ_eval, 2);  // second derivative
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*rhoLMat, *rho, *LMat);
    applyTensor(CUMat, UMat);

    Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*CUrhoLMat, *rhoLMat, *CUMat);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*NCUrhoLMat, *CUrhoLMat, *fe_->N());

    /*** Evaluate weak form of the residual. ***/
    Intrepid::FunctionSpaceTools::integrate<Real>(*J[0][0],
                                                  *NCUrhoLMat,        // N L' B' C ddrho B U
                                                  *fe_->NdetJ(),      // N
                                                  Intrepid::COMP_CPP,
                                                  false);

    // Combine the Hessians.
    fieldHelper_->combineFieldCoeff(hess, J);
  }

  void RieszMap_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_1): Not implemented.");

    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int d = fe_->gradN()->dimension(3);
 
    // Initialize Jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j) {
        J[i].push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,f)));
      }
    }

    for (int i=0; i<d; ++i) {
      *(J[i][i]) = *(fe_->stiffMat());
      Intrepid::RealSpaceTools<Real>::add(*(J[i][i]),*(fe_->massMat()));
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  void RieszMap_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & riesz) {
    //throw Exception::NotImplemented(">>> (PDE_TopoOpt::RieszMap_2): Not implemented.");

    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int d = fe_->gradN()->dimension(3);
 
    // Initialize Jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j) {
        J[i].push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,f)));
      }
    }

    for (int i=0; i<d; ++i) {
      *(J[i][i]) = *(fe_->massMat());
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > getFields() {
    return basisPtrs_;
  }

  void setCellNodes(const Teuchos::RCP<Intrepid::FieldContainer<Real> > &volCellNodes,
                    const std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds) {
    volCellNodes_ = volCellNodes;
    bdryCellNodes_ = bdryCellNodes;
    bdryCellLocIds_ = bdryCellLocIds;
    // Finite element definition.
    fe_ = Teuchos::rcp(new FE<Real>(volCellNodes_,basisPtr_,cellCub_));
    fidx_ = fe_->getBoundaryDofs();
    // Construct boundary FE
    /*int sideset = 6;
    int numLocSides = bdryCellNodes[sideset].size();
    feBdry_.resize(numLocSides);
    for (int j = 0; j < numLocSides; ++j) {
      if (bdryCellNodes[sideset][j] != Teuchos::null) {
        feBdry_[j] = Teuchos::rcp(new FE<Real>(bdryCellNodes[sideset][j],basisPtr_,bdryCub_,j));
      }
    }*/
    computeNBmats();
  }

  void setFieldPattern(const std::vector<std::vector<int> > & fieldPattern) {
    fieldPattern_ = fieldPattern;
    fieldHelper_ = Teuchos::rcp(new FieldHelper<Real>(numFields_, numDofs_, numFieldDofs_, fieldPattern_));
  }

  const Teuchos::RCP<FE<Real> > getFE(void) const {
    return fe_;
  }

  const std::vector<Teuchos::RCP<FE<Real> > > getBdryFE(void) const {
    return feBdry_;
  }

  const std::vector<std::vector<int> > getBdryCellLocIds(const int sideset = 6) const {
    return bdryCellLocIds_[sideset];
  }

  const Teuchos::RCP<FieldHelper<Real> > getFieldHelper(void) const {
    return fieldHelper_;
  }

  const Teuchos::RCP<Load<Real> > getLoad(void) const {
    return load_;
  }

}; // PDE_TopoOpt



template <class Real>
class PDE_Filter : public PDE<Real> {
private:
  // Finite element basis information
  Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtr_;
  std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
  // Cell cubature information
  Teuchos::RCP<Intrepid::Cubature<Real> > cellCub_;
  // Cell node information
  Teuchos::RCP<Intrepid::FieldContainer<Real> > volCellNodes_;
  // Finite element definition
  Teuchos::RCP<FE<Real> > fe_;
  // Field pattern, offsets, etc.
  std::vector<std::vector<int> > fieldPattern_;  // local Field/DOF pattern; set from DOF manager 
  int numFields_;                                // number of fields (equations in the PDE)
  int numDofs_;                                  // total number of degrees of freedom for all (local) fields
  std::vector<int> offset_;                      // for each field, a counting offset
  std::vector<int> numFieldDofs_;                // for each field, number of degrees of freedom
  
  // Problem parameters.
  Real lengthScale_;

  Teuchos::RCP<FieldHelper<Real> > fieldHelper_;

public:
  PDE_Filter(Teuchos::ParameterList &parlist) {
    // Finite element fields.
    basisPtr_ = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >);
    // Quadrature rules.
    shards::CellTopology cellType = basisPtr_->getBaseCellTopology();            // get the cell type from the basis
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                           // create cubature factory
    int cubDegree = parlist.sublist("Problem").get("Cubature Degree", 2);        // set cubature degree, e.g., 2
    cellCub_ = cubFactory.create(cellType, cubDegree);                           // create default cubature
    int d = cellType.getDimension();

    basisPtrs_.clear();
    for (int i=0; i<d; ++i) {
      basisPtrs_.push_back(basisPtr_);  // Filter components; there is only one, but we need d because of the infrastructure.
    }

    // Other problem parameters.
    Real filterRadius = parlist.sublist("Problem").get("Filter Radius",  0.1);
    lengthScale_ = std::pow(filterRadius, 2)/12.0;

    numDofs_ = 0;
    numFields_ = basisPtrs_.size();
    offset_.resize(numFields_);
    numFieldDofs_.resize(numFields_);
    for (int i=0; i<numFields_; ++i) {
      if (i==0) {
        offset_[i]  = 0;
      }
      else {
        offset_[i]  = offset_[i-1] + basisPtrs_[i-1]->getCardinality();
      }
      numFieldDofs_[i] = basisPtrs_[i]->getCardinality();
      numDofs_ += numFieldDofs_[i];
    }
  }

  void residual(Teuchos::RCP<Intrepid::FieldContainer<Real> > & res,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
 
    // Initialize residuals.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > R(d);
    for (int i=0; i<d; ++i) {
      R[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f));
    }

    // Split u_coeff and z_coeff into components.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U;
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(Z, z_coeff);

    // Evaluate/interpolate finite element fields on cells.
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > valU_eval(d);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > gradU_eval(d);
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > valZ_eval(d);
    for (int i=0; i<d; ++i) {
      valU_eval[i]  =  Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      valZ_eval[i]  =  Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      gradU_eval[i] =  Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p, d));
    }
    for (int i=0; i<d; ++i) {
      fe_->evaluateValue(valU_eval[i], U[i]);
      fe_->evaluateValue(valZ_eval[i], Z[i]);
      fe_->evaluateGradient(gradU_eval[i], U[i]);
    }

    for (int i=0; i<d; ++i) {
      Intrepid::RealSpaceTools<Real>::scale(*gradU_eval[i], lengthScale_);
      Intrepid::RealSpaceTools<Real>::scale(*valZ_eval[i],  static_cast<Real>(-1));
    }

    /*** Evaluate weak form of the residual. ***/
    for (int i=0; i<d; ++i) {
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *gradU_eval[i],        // R*gradU
                                                    *fe_->gradNdetJ(),     // gradN
                                                    Intrepid::COMP_CPP,
                                                    false);
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *valU_eval[i],         // U
                                                    *fe_->NdetJ(),         // N
                                                    Intrepid::COMP_CPP,
                                                    true);
      Intrepid::FunctionSpaceTools::integrate<Real>(*R[i],
                                                    *valZ_eval[i],         // -Z
                                                    *fe_->NdetJ(),         // N
                                                    Intrepid::COMP_CPP,
                                                    true);
    }

    // Combine the residuals.
    fieldHelper_->combineFieldCoeff(res, R);
  }

  void Jacobian_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int d = fe_->gradN()->dimension(3);
 
    // Initialize Jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j) {
        J[i].push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,f)));
      }
    }

    /*** Evaluate weak form of the Jacobian. ***/
    for (int i=0; i<d; ++i) {
      *(J[i][i]) = *(fe_->stiffMat());
      Intrepid::RealSpaceTools<Real>::scale(*(J[i][i]), lengthScale_);    // ls*gradN1 . gradN2
      Intrepid::RealSpaceTools<Real>::add(*(J[i][i]),*(fe_->massMat()));  // + N1 * N2
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }


  void Jacobian_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int d = fe_->gradN()->dimension(3);
 
    // Initialize Jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j) {
        J[i].push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,f)));
      }
    }

    /*** Evaluate weak form of the Jacobian. ***/
    for (int i=0; i<d; ++i) {
      *(J[i][i]) = *(fe_->massMat());
      Intrepid::RealSpaceTools<Real>::scale(*J[i][i], static_cast<Real>(-1));  // -N1 * N2
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(jac, J);
  }

  void Hessian_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Filter::Hessian_11): Hessian is zero.");
  }

  void Hessian_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Filter::Hessian_12): Hessian is zero.");
  }

  void Hessian_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Filter::Hessian_21): Hessian is zero.");
  }

  void Hessian_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & l_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> (PDE_Filter::Hessian_22): Hessian is zero.");
  }

  void RieszMap_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & riesz) {
    throw Exception::NotImplemented(">>> (PDE_Filter::RieszMap_1): Not implemented.");
    // Retrieve dimensions.
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int d = fe_->gradN()->dimension(3);
 
    // Initialize Jacobians.
    std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > J(d);
    for (int i=0; i<d; ++i) {
      for (int j=0; j<d; ++j) {
        J[i].push_back(Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,f)));
      }
    }

    for (int i=0; i<d; ++i) {
      *(J[i][i]) = *(fe_->stiffMat());
      Intrepid::RealSpaceTools<Real>::add(*(J[i][i]),*(fe_->massMat()));
    }

    // Combine the jacobians.
    fieldHelper_->combineFieldCoeff(riesz, J);
  }

  void RieszMap_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & riesz) {
    throw Exception::NotImplemented(">>> (PDE_Filter::RieszMap_2): Not implemented.");
  }

  std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > getFields() {
    return basisPtrs_;
  }

  void setCellNodes(const Teuchos::RCP<Intrepid::FieldContainer<Real> > &volCellNodes,
                    const std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > &bdryCellNodes,
                    const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds) {
    volCellNodes_ = volCellNodes;
    // Finite element definition.
    fe_ = Teuchos::rcp(new FE<Real>(volCellNodes_,basisPtr_,cellCub_));
  }

  void setFieldPattern(const std::vector<std::vector<int> > & fieldPattern) {
    fieldPattern_ = fieldPattern;
    fieldHelper_ = Teuchos::rcp(new FieldHelper<Real>(numFields_, numDofs_, numFieldDofs_, fieldPattern_));
  }

}; // PDE_Filter

#endif
