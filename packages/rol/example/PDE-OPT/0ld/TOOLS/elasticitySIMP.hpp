// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PDEOPT_ELASTICITYSIMP_H
#define ROL_PDEOPT_ELASTICITYSIMP_H

#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "../TOOLS/elasticity.hpp"
#include "../TOOLS/materials.hpp"

template<class Real>
class ElasticitySIMP : public Elasticity <Real> {
protected:
  Real initDensity_;
  Real minDensity_;
  int powerP_;
  
  std::vector<ROL::Ptr<Material_SIMP<Real> > >SIMPmaterial_;
  ROL::Ptr<Tpetra::MultiVector<> >  myDensity_;
  ROL::Ptr<Tpetra::MultiVector<> >  myCellArea_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > CBMat0_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > gradgradMats0_;

public:
  ElasticitySIMP() {}
  ~ElasticitySIMP() {}

  virtual void ElasticitySIMP_Initialize(const ROL::Ptr<const Teuchos::Comm<int> > &comm,
                                         const Teuchos::RCP<Teuchos::ParameterList> &parlist,
                                         const ROL::Ptr<std::ostream> &outStream) {
    this->Initialize(comm, parlist, outStream);
    // new material parameters
    Teuchos::ParameterList &list = this->parlist_->sublist("ElasticitySIMP");
    powerP_      = list.get<int>("SIMP Power");
    initDensity_ = list.get<Real>("Initial Density");
    minDensity_  = list.get<Real>("Minimum Density");
  }


  virtual void SetSIMPParallelStructure() { 
    PDE_FEM<Real>::SetParallelStructure();
    this->myCellMap_ = ROL::makePtr<Tpetra::Map<>>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                              this->myCellIds_, 0, this->commPtr_);
  }

  ROL::Ptr<const Tpetra::Map<> > getDomainMapA() const {
    return this->matA_->getDomainMap();
  }

  ROL::Ptr<const Tpetra::Map<> > getCellMap() const {
    return this->myCellMap_;
  }
  
  ROL::Ptr<Tpetra::MultiVector<> > getMaterialDensity() const {
    return myDensity_;
  }
 
// return cell measures 
  ROL::Ptr<Tpetra::MultiVector<> > getCellAreas() {
    myCellArea_ = ROL::makePtr<Tpetra::MultiVector<>>(this->myCellMap_, 1, true);
    for (int i=0; i<this->numCells_; i++){
    	myCellArea_ -> replaceGlobalValue(this->myCellIds_[i], 0, this->myCellMeasure_[i]);
    }
    return myCellArea_;
  }
//
  void ApplyBCToVec (const ROL::Ptr<Tpetra::MultiVector<> > & u) {
    Real zero(0.0);
    // u is myOverlapMap_
    for (int i=0; i<this->myDirichletDofs_.size(); ++i) {
      if (this->myOverlapMap_->isNodeGlobalElement(this->myDirichletDofs_[i]))
	u->replaceGlobalValue(this->myDirichletDofs_[i], 0, zero);
    }
  }

  void resetMaterialDensity (const Real val) {
    myDensity_ = ROL::makePtr<Tpetra::MultiVector<>>(this->myCellMap_, 1, true);
    myDensity_->putScalar(val);
    renewMaterialVector();
  }

  void updateMaterialDensity (const ROL::Ptr<const Tpetra::MultiVector<> > & newDensity) {
    myDensity_ = ROL::makePtr<Tpetra::MultiVector<>>(this->myCellMap_, 1, true);
    Tpetra::deep_copy(*myDensity_, *newDensity);	
    renewMaterialVector();
  }

  void renewMaterialVector () {
    Teuchos::ArrayRCP<const Real> densData = myDensity_->get1dView();
    for(int i=0; i<this->numCells_; i++) {
      Real dens = densData[myDensity_->getMap()->getLocalElement(this->myCellIds_[i])];
      SIMPmaterial_[i]->setDensity(dens);
    }
  }

  virtual void CreateMaterial() {
    for(int i=0; i<this->numCells_; i++) {
      ROL::Ptr<Material_SIMP<Real> > CellMaterial = ROL::makePtr<Material_SIMP<Real>>();
      CellMaterial->InitializeSIMP(this->spaceDim_, this->planeStrain_, this->E_,
                                   this->poissonR_, initDensity_, powerP_, minDensity_);
      this->materialTensorDim_ = CellMaterial->GetMaterialTensorDim();
      SIMPmaterial_.push_back(CellMaterial);
    }
    resetMaterialDensity (initDensity_);
  }


  virtual void ComputeLocalSystemMats(bool ifInitial) {
    int full_lfs = this->lfs_ * this->spaceDim_;
    if(!ifInitial) {
      #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*LocalAssemblyTime_example_PDEOPT_TOOLS_PDEFEM_GLOB);
      #endif
      renewMaterialVector();	
      //this->gradgradMats_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numCells_, full_lfs, full_lfs);
      Construct_CBmats(ifInitial);
      Intrepid::FunctionSpaceTools::integrate<Real>(*this->gradgradMats_, // compute local grad.grad (stiffness) matrices
                                                    *this->CBMat_,
                                                    *this->BMatWeighted_,
                                                    Intrepid::COMP_CPP);
      return;
    }	
    CreateMaterial();
    this->gradgradMats0_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numCells_, full_lfs, full_lfs);
    this->gradgradMats_  = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numCells_, full_lfs, full_lfs);
    this->valvalMats_    = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numCells_, full_lfs, full_lfs);
    
    this->BMat_		 = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numCells_, full_lfs, this->numCubPoints_, this->materialTensorDim_);
    this->BMatWeighted_	 = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numCells_, full_lfs, this->numCubPoints_, this->materialTensorDim_);
    CBMat0_ 		 = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numCells_, full_lfs, this->numCubPoints_, this->materialTensorDim_);
    this->CBMat_         = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numCells_, full_lfs, this->numCubPoints_, this->materialTensorDim_);
    this->NMat_		 = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numCells_, full_lfs, this->numCubPoints_, this->spaceDim_);
    this->NMatWeighted_  = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numCells_, full_lfs, this->numCubPoints_, this->spaceDim_);
    this->Construct_N_B_mats();
    Construct_CBmats(ifInitial);
    Intrepid::FunctionSpaceTools::integrate<Real>(*this->gradgradMats0_, // compute local grad.grad (stiffness) matrices
                                                  *CBMat0_,
                                                  *this->BMatWeighted_,
                                                  Intrepid::COMP_CPP);
    Intrepid::FunctionSpaceTools::integrate<Real>(*this->gradgradMats_,  // compute local grad.grad (stiffness) matrices
                                                  *this->CBMat_,
                                                  *this->BMatWeighted_,
                                                  Intrepid::COMP_CPP);
    Intrepid::FunctionSpaceTools::integrate<Real>(*this->valvalMats_,    // compute local val.val (mass) matrices
                                                  *this->NMat_,
                                                  *this->NMatWeighted_,
                                                  Intrepid::COMP_CPP);
  }

  // test matrices
  virtual void test_mats() {
    ROL::Ptr<Intrepid::FieldContainer<Real> > test_Jaco_Mat;
    ROL::Ptr<Intrepid::FieldContainer<Real> > test_Grad_Mat;
    ROL::Ptr<Intrepid::FieldContainer<Real> > test_N_Mat;
    ROL::Ptr<Intrepid::FieldContainer<Real> > test_B_Mat;
    ROL::Ptr<Intrepid::FieldContainer<Real> > test_K0_Mat;
    ROL::Ptr<Intrepid::FieldContainer<Real> > test_K_Mat;
    ROL::Ptr<Intrepid::FieldContainer<Real> > test_M_Mat;
    ROL::Ptr<Intrepid::FieldContainer<Real> > test_F_Vec;
    
    test_Jaco_Mat = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->spaceDim_, this->spaceDim_);
    test_Grad_Mat = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->lfs_, this->spaceDim_);
    test_N_Mat = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numLocalDofs_, this->spaceDim_);
    test_B_Mat = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numLocalDofs_, this->materialTensorDim_);
    test_K0_Mat = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numLocalDofs_, this->numLocalDofs_);
    test_K_Mat = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numLocalDofs_, this->numLocalDofs_);
    test_M_Mat = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numLocalDofs_, this->numLocalDofs_);
    test_F_Vec = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numLocalDofs_, 1);
    
    for(int i=0; i<this->spaceDim_; i++) {
       for(int j=0; j<this->spaceDim_; j++) {
          (*test_Jaco_Mat)(i, j) = (*this->cellJac_)(0, 0, i, j);
       }
    }
    for(int i=0; i<this->numLocalDofs_; i++) {
       for(int j=0; j<this->spaceDim_; j++) {
    	if(i<this->lfs_)
    		(*test_Grad_Mat)(i, j) = (*this->gradReference_)(i, 0, j);
    	
    	(*test_N_Mat)(i, j) = (*this->NMat_)(0, i, 0, j);	
       }
       for(int j=0; j<this->materialTensorDim_; j++) {
    	(*test_B_Mat)(i, j) = (*this->BMat_)(0, i, 0, j);	
       }
       for(int j=0; j<this->numLocalDofs_; j++) {
    	(*test_K0_Mat)(i, j) = (*this->gradgradMats0_)(0, i, j);
    	(*test_K_Mat)(i, j) = (*this->gradgradMats_)(0, i, j);
    	(*test_M_Mat)(i, j) =   (*this->valvalMats_)(0, i, j);
       }
       (*test_F_Vec)(i, 0) = (*this->datavalVecF_)(0, i);
    }
    std::cout<<*(SIMPmaterial_[0]->GetMaterialTensor())<<std::endl;
    std::cout<<*test_Jaco_Mat<<std::endl;
    std::cout<<*test_Grad_Mat<<std::endl;
    std::cout<<*test_N_Mat<<std::endl;
    std::cout<<*test_B_Mat<<std::endl;
    std::cout<<*test_M_Mat<<std::endl;
    std::cout<<*test_F_Vec<<std::endl;
    std::cout<<*test_K0_Mat<<std::endl;
    std::cout<<*test_K_Mat<<std::endl;
  }


  void Construct_CBmats(const bool ifInitial) {
    if(ifInitial)
      CBMat0_->initialize(0.0);

     this->CBMat_->initialize(0.0);
     Real SIMPScale;
     for (int i=0; i<this->numCells_; ++i) {
      
      SIMPScale = SIMPmaterial_[i]->getSIMPScaleFactor();
      ROL::Ptr<Intrepid::FieldContainer<Real> > materialMat = SIMPmaterial_[i]->GetMaterialTensor();
      for (int j=0; j<this->numCubPoints_; ++j) {
        for (int m=0; m<this->lfs_*this->spaceDim_; m++) {
          for (int n=0; n<this->materialTensorDim_; n++) {
            for (int k=0; k<this->materialTensorDim_; k++) {
              if(ifInitial)
                (*CBMat0_)(i, m, j, n) += (*this->BMat_)(i, m, j, k) * (*materialMat)(k, n);
              (*this->CBMat_)(i, m, j, n) += SIMPScale * (*this->BMat_)(i, m, j, k) * (*materialMat)(k, n);
            }
          }
        }
      }
    }
  }


}; // class ElasticitySIMPData

#endif
