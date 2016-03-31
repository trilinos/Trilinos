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
    /************************************/
    /*** Retrieve communication data. ***/
    /************************************/
    this->commPtr_   = comm;
    this->parlist_   = parlist;
    this->outStream_ = outStream;
    this->myRank_    = comm->getRank();
    this->numProcs_  = comm->getSize();
    *outStream << "Total number of processors: " << this->numProcs_ << std::endl;
    /****************************************************************************/
    /*** Initialize mesh / finite element fields / degree-of-freedom manager. ***/
    /****************************************************************************/
    this->basisOrder_ = parlist->sublist("Elasticity").get("Order of FE discretization", 1);
    planeStrain_      = parlist->sublist("Elasticity").get("Plane Strain", false);
    E_ 	              = parlist->sublist("Elasticity").get("Young's Modulus", 1.0);
    poissonR_         = parlist->sublist("Elasticity").get("Poisson Ratio", 0.3);
     
    // Create Mesh manager.
    this->meshMgr_ = Teuchos::rcp(new MeshManager_Rectangle<Real>(*parlist));
    PDE_FEM<Real>::printMeshData(*outStream);

    // Finite element fields.
    Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtr;
    if (this->basisOrder_ == 1) {
      basisPtr = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >);
    }
    else if (this->basisOrder_ == 2) {
      basisPtr = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >);
    }
    
    // Retrieve some basic cell information.
    this->cellType_        = basisPtr -> getBaseCellTopology(); // get the cell type from any basis
    this->spaceDim_        = this->cellType_.getDimension();    // retrieve spatial dimension
    this->numNodesPerCell_ = this->cellType_.getNodeCount();    // retrieve number of nodes per cell
   
    // Create basis  
    this->basisPtrs_.resize(this->spaceDim_, Teuchos::null);
    for(int i=0; i<this->spaceDim_; i++) {
      this->basisPtrs_[i] = basisPtr;
    }
	
    // DOF coordinate interface.
    this->coord_iface_ = Teuchos::rcp_dynamic_cast<Intrepid::DofCoordsInterface<Intrepid::FieldContainer<Real> > >(this->basisPtrs_[0]);
    // Degree-of-freedom manager.
    this->dofMgr_ = Teuchos::rcp(new DofManager<Real>(this->meshMgr_, this->basisPtrs_));
    // Retrieve total number of cells in the mesh.
    this->totalNumCells_ = this->meshMgr_->getNumCells();
    // Retrieve total number of degrees of freedom in the mesh.
    this->totalNumDofs_ = this->dofMgr_->getNumDofs();
  }


  virtual void CreateMaterial(void) {
    for(int i=0; i<this->numCells_; i++) {
      Teuchos::RCP<Material<Real> > CellMaterial = Teuchos::rcp(new Material<Real>());
      CellMaterial-> InitializeMaterial(this->spaceDim_, planeStrain_, E_, poissonR_);
      materialTensorDim_ = CellMaterial->GetMaterialTensorDim();
      material_.push_back(CellMaterial);
    }
  }

  virtual void ComputeLocalSystemMats(void) {
    int full_lfs = this->lfs_ * this->spaceDim_;
/*  if(this->numLocalDofs_ != full_lfs)
	std::cout<<"numLocalDofs_ DOES NOT match full_lfs, numLocalDofs_="<<this->numLocalDofs_<<", full_lfs="<<full_lfs<<std::endl;
    else
	std::cout<<"numLocalDofs_ DOES match full_lfs, numLocalDofs_="<<this->numLocalDofs_<<", full_lfs="<<full_lfs<<std::endl;
*/	
    this->gradgradMats_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numCells_, full_lfs, full_lfs));
    this->valvalMats_   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numCells_, full_lfs, full_lfs));
    
    // construct mats
    CreateMaterial();
    
    BMat_         = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numCells_, full_lfs, this->numCubPoints_, materialTensorDim_));
    BMatWeighted_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numCells_, full_lfs, this->numCubPoints_, materialTensorDim_));
    CBMat_        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numCells_, full_lfs, this->numCubPoints_, materialTensorDim_));
    NMat_         = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numCells_, full_lfs, this->numCubPoints_, this->spaceDim_));
    NMatWeighted_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numCells_, full_lfs, this->numCubPoints_, this->spaceDim_));
    Construct_N_B_mats();
    Construct_CBmats();

    Intrepid::FunctionSpaceTools::integrate<Real>(*this->gradgradMats_,                            // compute local grad.grad (stiffness) matrices
                                                  *CBMat_,
                                                  *BMatWeighted_,
                                                  Intrepid::COMP_CPP);
    Intrepid::FunctionSpaceTools::integrate<Real>(*this->valvalMats_,                              // compute local val.val (mass) matrices
                                                  *NMat_,
                                                  *NMatWeighted_,
                                                  Intrepid::COMP_CPP);
  } 
   
  virtual void ComputeLocalForceVec(void) {
    int full_lfs = this->lfs_ * this->spaceDim_;
    this->dataF_                = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numCells_, this->numCubPoints_, this->spaceDim_));
    this->datavalVecF_          = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numCells_, full_lfs));

    for (int i=0; i<this->numCells_; ++i) {                                                         // evaluate functions at these points
      for (int j=0; j<this->numCubPoints_; ++j) {
	for (int k=0; k<this->spaceDim_; k++) {
		if(this->spaceDim_==2)
        	{	(*this->dataF_)(i, j, k) = funcRHS_2D((*this->cubPointsPhysical_)(i, j, 0), (*this->cubPointsPhysical_)(i, j, 1), k);
			//std::cout<<(*this->dataF_)(i, j, k)<<std::endl;
		}
		else if(this->spaceDim_==3)
        		(*this->dataF_)(i, j, k) = funcRHS_3D((*this->cubPointsPhysical_)(i, j, 0), (*this->cubPointsPhysical_)(i, j, 1), (*this->cubPointsPhysical_)(i, j, 2), k);
		else
			std::cout<<"1D RHS not implemented"<<std::endl;
	}
      }
    }
    Intrepid::FunctionSpaceTools::integrate<Real>(*this->datavalVecF_,                              // compute local data.val vectors for RHS F
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

  virtual Real funcRHS_2D(const Real &x1, const Real &x2, const int k) {
    Real zero(0);
    return zero; 
  }

  virtual Real funcRHS_3D(const Real &x1, const Real &x2, const Real &x3, const int k) {
    Real zero(0);
    return zero; 
  }

}; // class Elasticity

#endif
