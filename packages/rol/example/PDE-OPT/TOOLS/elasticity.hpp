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

#ifndef ROL_PDEOPT_ELASTICITY_H
#define ROL_PDEOPT_ELASTICITY_H

#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "../TOOLS/PDE_FEM.hpp"
#include "../TOOLS/materials.hpp"

template<class Real>
class Elasticity : public PDE_FEM <Real> {
protected:

  Real E_;
  Real poissonR_;
  bool planeStrain_;
  
  std::vector<Teuchos::RCP<Material<Real> > > material_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > BMat_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > BMatWeighted_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > CBMat_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > NMat_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > NMatWeighted_;
  int materialTensorDim_;

public:

  Elasticity() {}
  virtual ~Elasticity() {}

  virtual void Initialize(const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                          const Teuchos::RCP<Teuchos::ParameterList> &parlist,
                          const Teuchos::RCP<std::ostream> &outStream) {
    /****************************************************************************/
    /*** Initialize the base PDE_FEM class. *************************************/
    /****************************************************************************/
    PDE_FEM<Real>::Initialize(comm,parlist,outStream);

    /****************************************************************************/
    /*** Grab the elasticity information. ***************************************/
    /****************************************************************************/
    planeStrain_      = parlist->sublist("Elasticity").get("Plane Strain", false);
    E_ 	              = parlist->sublist("Elasticity").get("Young's Modulus", 1.0);
    poissonR_         = parlist->sublist("Elasticity").get("Poisson Ratio", 0.3);
     
    /****************************************************************************/
    /*** Initialize mesh / finite element fields / degree-of-freedom manager. ***/
    /****************************************************************************/
    int basisOrder = parlist->sublist("PDE FEM").get("Order of FE Discretization", 1);
    Teuchos::RCP<MeshManager<Real> > meshMgr = Teuchos::rcp(new MeshManager_Rectangle<Real>(*parlist));
    int spaceDim = 2;
    std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs(spaceDim,Teuchos::null);
    for (int k = 0; k < spaceDim; ++k) {
      if (basisOrder == 1) {
        basisPtrs[k] = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >);
      }
      else if (basisOrder == 2) {
        basisPtrs[k] = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >);
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          ">>> (Elasticity::Initialize): Basis Order is out of bounds!");
      }
    }
    PDE_FEM<Real>::SetDiscretization(meshMgr,basisPtrs);
    PDE_FEM<Real>::printMeshData(*outStream);
  }

  virtual void CreateMaterial(void) {
    int numCells = PDE_FEM<Real>::GetNumCells();
    int spaceDim = PDE_FEM<Real>::GetSpaceDim();
    for(int i = 0; i < numCells; ++i) {
      Teuchos::RCP<Material<Real> > CellMaterial = Teuchos::rcp(new Material<Real>());
      CellMaterial-> InitializeMaterial(spaceDim, planeStrain_, E_, poissonR_);
      materialTensorDim_ = CellMaterial->GetMaterialTensorDim();
      material_.push_back(CellMaterial);
    }
  }

  virtual void ComputeLocalSystemMats(void) {
    int spaceDim = PDE_FEM<Real>::GetSpaceDim();
    int numCells = PDE_FEM<Real>::GetNumCells();
    int numCubPoints = PDE_FEM<Real>::GetNumCubPoints();
    int lfs = PDE_FEM<Real>::GetLocalFieldSize();
    int full_lfs = lfs * spaceDim;
/* 
    int numLocalDofs = PDE_FEM<Real>::GetNumLocalDofs();
    if(numLocalDofs != full_lfs) {
      std::cout << "numLocalDofs DOES NOT match full_lfs, numLocalDofs = " << numLocalDofs
                << ", full_lfs = " << full_lfs << std::endl;
    }
    else {
      std::cout << "numLocalDofs DOES match full_lfs, numLocalDofs = " << numLocalDofs
                << ", full_lfs = " << full_lfs << std::endl;
    }
*/	
    this->gradgradMats_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells, full_lfs, full_lfs));
    this->valvalMats_   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells, full_lfs, full_lfs));
    
    // construct mats
    CreateMaterial();
    
    BMat_         = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells, full_lfs, numCubPoints, materialTensorDim_));
    BMatWeighted_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells, full_lfs, numCubPoints, materialTensorDim_));
    CBMat_        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells, full_lfs, numCubPoints, materialTensorDim_));
    NMat_         = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells, full_lfs, numCubPoints, spaceDim));
    NMatWeighted_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells, full_lfs, numCubPoints, spaceDim));
    Construct_N_B_mats();
    Construct_CBmats();

    Intrepid::FunctionSpaceTools::integrate<Real>(*this->gradgradMats_, // compute local grad.grad (stiffness) matrices
                                                  *CBMat_,
                                                  *BMatWeighted_,
                                                  Intrepid::COMP_CPP);
    Intrepid::FunctionSpaceTools::integrate<Real>(*this->valvalMats_,   // compute local val.val (mass) matrices
                                                  *NMat_,
                                                  *NMatWeighted_,
                                                  Intrepid::COMP_CPP);
  } 
   
  virtual void ComputeLocalForceVec(void) {
    int spaceDim = PDE_FEM<Real>::GetSpaceDim();
    int numCells = PDE_FEM<Real>::GetNumCells();
    int numCubPoints = PDE_FEM<Real>::GetNumCubPoints();
    int lfs = PDE_FEM<Real>::GetLocalFieldSize();
    int full_lfs = lfs * spaceDim;
    this->dataF_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells, numCubPoints, spaceDim));
    this->datavalVecF_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells, full_lfs));

    std::vector<Real> x(spaceDim), F(spaceDim);
    for (int i = 0; i < numCells; ++i) { // evaluate functions at these points
      for (int j = 0; j < numCubPoints; ++j) {
        for (int k = 0; k < spaceDim; ++k) {
          x[k] = (*this->cubPointsPhysical_)(i,j,k);
	}
        funcRHS(F,x);
        for (int k = 0; k < spaceDim; ++k) {
          (*this->dataF_)(i,j,k) = F[k];
        }
      }
    }
    Intrepid::FunctionSpaceTools::integrate<Real>(*this->datavalVecF_, // compute local data.val vectors for RHS F
						  *this->dataF_,
                                                  *NMatWeighted_,
                                                  Intrepid::COMP_CPP);
  }

  void Construct_N_B_mats(void) {
    //std::cout<<"Computing N and B mats."<<std::endl;
    BMat_->initialize(0.0);
    NMat_->initialize(0.0);
    BMatWeighted_->initialize(0.0);
    NMatWeighted_->initialize(0.0);
    
    if(this->spaceDim_==2) {
      for (int i=0; i<this->numCells_; ++i) { // evaluate functions at these points
        for (int j=0; j<this->numCubPoints_; ++j) {
      	  for (int k=0; k<this->lfs_; ++k) {
            (*NMat_)(i, k*this->spaceDim_+0, j, 0) = (*this->valPhysical_)(i, k, j);
            (*NMat_)(i, k*this->spaceDim_+1, j, 1) = (*this->valPhysical_)(i, k, j);
            (*NMatWeighted_)(i, k*this->spaceDim_+0, j, 0) = (*this->valPhysicalWeighted_)(i, k, j);
            (*NMatWeighted_)(i, k*this->spaceDim_+1, j, 1) = (*this->valPhysicalWeighted_)(i, k, j);
                
            (*BMat_)(i, k*this->spaceDim_+0, j, 0) = (*this->gradPhysical_)(i, k, j, 0);
            (*BMat_)(i, k*this->spaceDim_+1, j, 1) = (*this->gradPhysical_)(i, k, j, 1);
            (*BMat_)(i, k*this->spaceDim_+0, j, 2) = (*this->gradPhysical_)(i, k, j, 1);
            (*BMat_)(i, k*this->spaceDim_+1, j, 2) = (*this->gradPhysical_)(i, k, j, 0);
            
            (*BMatWeighted_)(i, k*this->spaceDim_+0, j, 0) = (*this->gradPhysicalWeighted_)(i, k, j, 0);
            (*BMatWeighted_)(i, k*this->spaceDim_+1, j, 1) = (*this->gradPhysicalWeighted_)(i, k, j, 1);
            (*BMatWeighted_)(i, k*this->spaceDim_+0, j, 2) = (*this->gradPhysicalWeighted_)(i, k, j, 1);
            (*BMatWeighted_)(i, k*this->spaceDim_+1, j, 2) = (*this->gradPhysicalWeighted_)(i, k, j, 0);
	  }
        }
      }
    }

    if(this->spaceDim_==3) {
      for (int i=0; i<this->numCells_; ++i) { // evaluate functions at these points
        for (int j=0; j<this->numCubPoints_; ++j) {
          for (int k=0; k<this->lfs_; ++k) {
            (*NMat_)(i, k*this->spaceDim_+0, j, 0) = (*this->valPhysical_)(i, k, j);
            (*NMat_)(i, k*this->spaceDim_+1, j, 1) = (*this->valPhysical_)(i, k, j);
            (*NMat_)(i, k*this->spaceDim_+2, j, 2) = (*this->valPhysical_)(i, k, j);
            (*NMatWeighted_)(i, k*this->spaceDim_+0, j, 0) = (*this->valPhysicalWeighted_)(i, k, j);
            (*NMatWeighted_)(i, k*this->spaceDim_+1, j, 1) = (*this->valPhysicalWeighted_)(i, k, j);
            (*NMatWeighted_)(i, k*this->spaceDim_+2, j, 2) = (*this->valPhysicalWeighted_)(i, k, j);
            
            (*BMat_)(i, k*this->spaceDim_+0, j, 0) = (*this->gradPhysical_)(i, k, j, 0);
            (*BMat_)(i, k*this->spaceDim_+1, j, 1) = (*this->gradPhysical_)(i, k, j, 1);
            (*BMat_)(i, k*this->spaceDim_+2, j, 2) = (*this->gradPhysical_)(i, k, j, 2);
            (*BMat_)(i, k*this->spaceDim_+1, j, 3) = (*this->gradPhysical_)(i, k, j, 2);
            (*BMat_)(i, k*this->spaceDim_+2, j, 3) = (*this->gradPhysical_)(i, k, j, 1);
            (*BMat_)(i, k*this->spaceDim_+0, j, 4) = (*this->gradPhysical_)(i, k, j, 2);
            (*BMat_)(i, k*this->spaceDim_+2, j, 4) = (*this->gradPhysical_)(i, k, j, 0);
            (*BMat_)(i, k*this->spaceDim_+0, j, 5) = (*this->gradPhysical_)(i, k, j, 1);
            (*BMat_)(i, k*this->spaceDim_+1, j, 5) = (*this->gradPhysical_)(i, k, j, 0);
                
            (*BMatWeighted_)(i, k*this->spaceDim_+0, j, 0) = (*this->gradPhysicalWeighted_)(i, k, j, 0);
            (*BMatWeighted_)(i, k*this->spaceDim_+1, j, 1) = (*this->gradPhysicalWeighted_)(i, k, j, 1);
            (*BMatWeighted_)(i, k*this->spaceDim_+2, j, 2) = (*this->gradPhysicalWeighted_)(i, k, j, 2);
            (*BMatWeighted_)(i, k*this->spaceDim_+1, j, 3) = (*this->gradPhysicalWeighted_)(i, k, j, 2);
            (*BMatWeighted_)(i, k*this->spaceDim_+2, j, 3) = (*this->gradPhysicalWeighted_)(i, k, j, 1);
            (*BMatWeighted_)(i, k*this->spaceDim_+0, j, 4) = (*this->gradPhysicalWeighted_)(i, k, j, 2);
            (*BMatWeighted_)(i, k*this->spaceDim_+2, j, 4) = (*this->gradPhysicalWeighted_)(i, k, j, 0);
            (*BMatWeighted_)(i, k*this->spaceDim_+0, j, 5) = (*this->gradPhysicalWeighted_)(i, k, j, 1);
            (*BMatWeighted_)(i, k*this->spaceDim_+1, j, 5) = (*this->gradPhysicalWeighted_)(i, k, j, 0);
          }
        }
      }
    } 
  }

  // test matrices
  virtual void test_mats(void) {
    Teuchos::RCP<Intrepid::FieldContainer<Real> > test_Jaco_Mat;
    Teuchos::RCP<Intrepid::FieldContainer<Real> > test_Grad_Mat;
    Teuchos::RCP<Intrepid::FieldContainer<Real> > test_N_Mat;
    Teuchos::RCP<Intrepid::FieldContainer<Real> > test_B_Mat;
    Teuchos::RCP<Intrepid::FieldContainer<Real> > test_K_Mat;
    Teuchos::RCP<Intrepid::FieldContainer<Real> > test_M_Mat;
    Teuchos::RCP<Intrepid::FieldContainer<Real> > test_F_Vec;
    
    test_Jaco_Mat = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->spaceDim_, this->spaceDim_));
    test_Grad_Mat = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->lfs_, this->spaceDim_));
    test_N_Mat = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numLocalDofs_, this->spaceDim_));
    test_B_Mat = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numLocalDofs_, materialTensorDim_));
    test_K_Mat = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numLocalDofs_, this->numLocalDofs_));
    test_M_Mat = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numLocalDofs_, this->numLocalDofs_));
    test_F_Vec = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numLocalDofs_, 1));
    
    for(int i=0; i<this->spaceDim_; i++) {
       for(int j=0; j<this->spaceDim_; j++) {
         (*test_Jaco_Mat)(i, j) = (*this->cellJac_)(0, 0, i, j);
       }
    }
    for(int i=0; i<this->numLocalDofs_; i++) {
      for(int j=0; j<this->spaceDim_; j++) {
        if(i<this->lfs_) {
          (*test_Grad_Mat)(i, j) = (*this->gradReference_)(i, 0, j);
        }
        (*test_N_Mat)(i, j) = (*NMat_)(0, i, 0, j);	
      }
      for(int j=0; j<materialTensorDim_; j++) {
        (*test_B_Mat)(i, j) = (*BMat_)(0, i, 0, j);	
      }
      for(int j=0; j<this->numLocalDofs_; j++) {
        (*test_K_Mat)(i, j) = (*this->gradgradMats_)(0, i, j);
        (*test_M_Mat)(i, j) =   (*this->valvalMats_)(0, i, j);
      }
      (*test_F_Vec)(i, 0) = (*this->datavalVecF_)(0, i);
    }
    std::cout<<*test_Jaco_Mat<<std::endl;
    std::cout<<*test_Grad_Mat<<std::endl;
    std::cout<<*test_N_Mat<<std::endl;
    std::cout<<*test_B_Mat<<std::endl;
    std::cout<<*test_K_Mat<<std::endl;
    std::cout<<*test_M_Mat<<std::endl;
    std::cout<<*test_F_Vec<<std::endl;
  }


  virtual void Construct_CBmats(void) {
    //std::cout<<"Computing CB mats."<<std::endl;
    CBMat_->initialize(0.0);
    for (int i=0; i<this->numCells_; ++i) {
      Teuchos::RCP<Intrepid::FieldContainer<Real> > materialMat = material_[i]->GetMaterialTensor();
      for (int j=0; j<this->numCubPoints_; ++j) {
        for (int m=0; m<this->lfs_*this->spaceDim_; m++) {
          for (int n=0; n<materialTensorDim_; n++) {
            for (int k=0; k<materialTensorDim_; k++) {
              (*CBMat_)(i, m, j, n) += (*BMat_)(i, m, j, k) * (*materialMat)(k, n);
            }
          }
        }
      }
    }
  }

  virtual void updateF(const std::vector<Real> &param) {
    int spaceDim = PDE_FEM<Real>::GetSpaceDim();
    int numCells = PDE_FEM<Real>::GetNumCells();
    int numCubPoints = PDE_FEM<Real>::GetNumCubPoints();
    int lfs = PDE_FEM<Real>::GetLocalFieldSize();
    int full_lfs = lfs * spaceDim;
    this->dataF_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells, numCubPoints, spaceDim));
    this->datavalVecF_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells, full_lfs));

    std::vector<Real> x(spaceDim), F(spaceDim);
    for (int i = 0; i < numCells; ++i) { // evaluate functions at these points
      for (int j = 0; j < numCubPoints; ++j) {
        for (int k = 0; k < spaceDim; ++k) {
          x[k] = (*this->cubPointsPhysical_)(i,j,k);
	}
        funcRHS(F,x,param);
        for (int k = 0; k < spaceDim; ++k) {
          (*this->dataF_)(i,j,k) = F[k];
        }
      }
    }
    Intrepid::FunctionSpaceTools::integrate<Real>(*this->datavalVecF_, // compute local data.val vectors for RHS F
						  *this->dataF_,
                                                  *NMatWeighted_,
                                                  Intrepid::COMP_CPP);
    // vecF_ requires assembly using vecF_overlap_ and redistribution
    this->vecF_ = Tpetra::rcp(new Tpetra::MultiVector<>(this->matA_->getRangeMap(), 1, true));
    this->vecF_overlap_ = Tpetra::rcp(new Tpetra::MultiVector<>(this->myOverlapMap_, 1, true));
    for (int i = 0; i < numCells; ++i) {
      for (int j = 0; j < this->numLocalDofs_; ++j) {
        this->vecF_overlap_->sumIntoGlobalValue(this->cellDofs_(this->myCellIds_[i],j),0,
                                                (*this->datavalVecF_)[i*(this->numLocalDofs_)+j]);
      }
    }
    Tpetra::Export<> exporter(this->vecF_overlap_->getMap(), this->vecF_->getMap());
    this->vecF_->doExport(*this->vecF_overlap_, exporter, Tpetra::ADD);

    Tpetra::deep_copy(*this->vecF_dirichlet_, *this->vecF_);
    int gDof = 0;
    for (int i = 0; i < numCells; ++i) {
      for (int j = 0; j < this->numLocalDofs_; ++j) {
        gDof = this->cellDofs_(this->myCellIds_[i], j);
        if (this->myUniqueMap_->isNodeGlobalElement(gDof)
            && this->check_myGlobalDof_on_boundary(gDof)) {
          this->vecF_dirichlet_->replaceGlobalValue(gDof, 0, 0);
        }
      }
    }
  }

  virtual void funcRHS(std::vector<Real> &F,
                 const std::vector<Real> &x) const {
    Real zero(0);
    F.clear();
    F.resize(this->spaceDim_,zero);
  }

  virtual void funcRHS(std::vector<Real> &F,
                 const std::vector<Real> &x,
                 const std::vector<Real> &param) const {
    Real zero(0);
    F.clear();
    F.resize(this->spaceDim_,zero);
  }

}; // class Elasticity

#endif
