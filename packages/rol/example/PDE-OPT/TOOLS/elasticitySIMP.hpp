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
  Real xmax_;
  Real ymax_;
  Real cx_;
  Real cy_;

  int loadCase_;
  std::vector<Real> param_;
  
  std::vector<Teuchos::RCP<Material_SIMP<Real> > >SIMPmaterial_;
  Teuchos::RCP<const Tpetra::Map<> >    myCellMap_;
  Teuchos::RCP<Tpetra::MultiVector<> >  myDensity_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > CBMat0_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradgradMats0_;

public:
  ElasticitySIMP() {}
  ~ElasticitySIMP() {}

  virtual void ElasticitySIMP_Initialize(const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                                         const Teuchos::RCP<Teuchos::ParameterList> &parlist,
                                         const Teuchos::RCP<std::ostream> &outStream) {
    Real zero(0), one(1), pi(M_PI);
    this->Initialize(comm, parlist, outStream);

    // new material parameters
    powerP_      = this->parlist_->sublist("ElasticitySIMP").get("SIMP Power", 3);
    initDensity_ = this->parlist_->sublist("ElasticitySIMP").get("Initial Density", one);
    minDensity_  = this->parlist_->sublist("ElasticitySIMP").get("Minimum Density", one);

    // Loading magnitude and angles
    loadCase_    = this->parlist_->sublist("ElasticitySIMP").get("Load Case", 0);
    param_.clear(); param_.resize(3);
    param_[0] = this->parlist_->sublist("ElasticitySIMP").get("Load Magnitude", one);
    param_[1] = this->parlist_->sublist("ElasticitySIMP").get("Load Polar Angle", -pi);
    param_[2] = this->parlist_->sublist("ElasticitySIMP").get("Load Azimuth Angle", zero);

    // Domain width and height
    xmax_        = this->parlist_->sublist("Geometry").get("Width", one);
    ymax_        = this->parlist_->sublist("Geometry").get("Height", one);
    int NX       = this->parlist_->sublist("Geometry").get("NX", 1);
    int NY       = this->parlist_->sublist("Geometry").get("NY", 1);
    cx_          = one/static_cast<Real>(NX);
    cy_          = one/static_cast<Real>(NY);
  }


  virtual void SetSIMPParallelStructure() { 
    PDE_FEM<Real>::SetParallelStructure();
    myCellMap_ = Teuchos::rcp(new Tpetra::Map<>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                              this->myCellIds_, 0, this->commPtr_));
  }

  Teuchos::RCP<const Tpetra::Map<> > getDomainMapA() const {
    return this->matA_->getDomainMap();
  }

  Teuchos::RCP<const Tpetra::Map<> > getCellMap() const {
    return myCellMap_;
  }
  
  Teuchos::RCP<Tpetra::MultiVector<> > getMaterialDensity() const {
    return myDensity_;
  }

  void ApplyBCToVec (const Teuchos::RCP<Tpetra::MultiVector<> > & u) {
    Real zero(0);
    // u is myOverlapMap_
    for (int i=0; i<this->myDirichletDofs_.size(); ++i) {
      if (this->myOverlapMap_->isNodeGlobalElement(this->myDirichletDofs_[i]))
	u->replaceGlobalValue(this->myDirichletDofs_[i], 0, zero);
    }
  }

  void resetMaterialDensity (const Real val) {
    myDensity_ = Teuchos::rcp(new Tpetra::MultiVector<>(myCellMap_, 1, true));
    myDensity_->putScalar(val);
    renewMaterialVector ();
  }

  void updateMaterialDensity (const Teuchos::RCP<const Tpetra::MultiVector<> > & newDensity) {
    myDensity_ = Teuchos::rcp(new Tpetra::MultiVector<>(myCellMap_, 1, true));
    Tpetra::deep_copy(*myDensity_, *newDensity);	
    renewMaterialVector ();
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
      Teuchos::RCP<Material_SIMP<Real> > CellMaterial = Teuchos::rcp(new Material_SIMP<Real>());
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
      renewMaterialVector();	
      this->gradgradMats_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numCells_, full_lfs, full_lfs));
      Construct_CBmats(ifInitial);
      Intrepid::FunctionSpaceTools::integrate<Real>(*this->gradgradMats_, // compute local grad.grad (stiffness) matrices
                                                    *this->CBMat_,
                                                    *this->BMatWeighted_,
                                                    Intrepid::COMP_CPP);
      return;
    }	
    CreateMaterial();
    this->gradgradMats0_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numCells_, full_lfs, full_lfs));
    this->gradgradMats_  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numCells_, full_lfs, full_lfs));
    this->valvalMats_    = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numCells_, full_lfs, full_lfs));
    
    this->BMat_		 = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numCells_, full_lfs, this->numCubPoints_, this->materialTensorDim_));
    this->BMatWeighted_	 = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numCells_, full_lfs, this->numCubPoints_, this->materialTensorDim_));
    CBMat0_ 		 = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numCells_, full_lfs, this->numCubPoints_, this->materialTensorDim_));
    this->CBMat_         = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numCells_, full_lfs, this->numCubPoints_, this->materialTensorDim_));
    this->NMat_		 = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numCells_, full_lfs, this->numCubPoints_, this->spaceDim_));
    this->NMatWeighted_   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numCells_, full_lfs, this->numCubPoints_, this->spaceDim_));
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
    Teuchos::RCP<Intrepid::FieldContainer<Real> > test_Jaco_Mat;
    Teuchos::RCP<Intrepid::FieldContainer<Real> > test_Grad_Mat;
    Teuchos::RCP<Intrepid::FieldContainer<Real> > test_N_Mat;
    Teuchos::RCP<Intrepid::FieldContainer<Real> > test_B_Mat;
    Teuchos::RCP<Intrepid::FieldContainer<Real> > test_K0_Mat;
    Teuchos::RCP<Intrepid::FieldContainer<Real> > test_K_Mat;
    Teuchos::RCP<Intrepid::FieldContainer<Real> > test_M_Mat;
    Teuchos::RCP<Intrepid::FieldContainer<Real> > test_F_Vec;
    
    test_Jaco_Mat = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->spaceDim_, this->spaceDim_));
    test_Grad_Mat = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->lfs_, this->spaceDim_));
    test_N_Mat = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numLocalDofs_, this->spaceDim_));
    test_B_Mat = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numLocalDofs_, this->materialTensorDim_));
    test_K0_Mat = Teuchos::rcp(new Intrepid::FieldContainer<Real>(this->numLocalDofs_, this->numLocalDofs_));
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
      Teuchos::RCP<Intrepid::FieldContainer<Real> > materialMat = SIMPmaterial_[i]->GetMaterialTensor();
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


  /*virtual Real funcRHS_2D(const Real &x1, const Real &x2, const int k) {
       if(loadCase_==1 && k==0 && x2 > ymax_-0.05 && x1 > xmax_/2-0.05 && x1 < xmax_/2+0.05)
               return -10.0;
       else if(loadCase_==2 && k==0 && x2 > ymax_-0.05 && x1 < 0.05)
               return -10.0;
       else if(loadCase_==3 && k==1 && x2 > ymax_-0.05 && x1 > xmax_/2-0.05 && x1 < xmax_/2+0.05)
               return 10.0;
       else if(loadCase_==4 && k==1 && x2 > ymax_-0.05 && x1 < 0.05)
               return 10.0;
       else if(loadCase_==5 && k==0 && x2 > ymax_-0.05)
               return -2.0;
       else if(loadCase_==6 && k==1 && x2 > ymax_-0.05)
               return 2.0;
       else if(loadCase_==7 && k==0 && x1 > xmax_-0.05)
               return -2.0;
       else if(loadCase_==8 && k==1 && x1 > xmax_-0.05)
               return 2.0;
       else 
               return 0.0; 
  }*/


  virtual void funcRHS(std::vector<Real> &F,
                 const std::vector<Real> &x) const {
    Real zero(0);
    Real loadMag = param_[0], loadAngle1 = param_[1], loadAngle2 = param_[2];
    F.clear(); F.resize(this->spaceDim_,zero);
    if (this->spaceDim_ == 2) {
      //Real eps(std::sqrt(ROL::ROL_EPSILON<Real>()));
      //Real cx(2*cx_+eps), cy(2*cy_+eps), half(0.5);
      Real cx(0.1), cy(0.05), half(0.5);
      Real x1 = x[0], x2 = x[1];
      switch(loadCase_) {
        case 0: { // Force applied to center of top boundary
          if ( x1 > half*(xmax_-cx) && x1 < half*(xmax_+cx) &&
               x2 > ymax_-cy ) {
            F[0] = loadMag*std::sin(loadAngle1);
            F[1] = loadMag*std::cos(loadAngle1);
          }
          break;
        }
        case 1: { // Force applied to top right corner
          if ( (x1 > xmax_-cx) &&
               (x2 > ymax_-cy) ) {
            F[0] = loadMag*std::sin(loadAngle1);
            F[1] = loadMag*std::cos(loadAngle1);
          }
          break;
        }
        case 2: { // Force applied to center of right boundary
          if ( (x1 > xmax_-cx) &&
               (x2 > half*(ymax_-cy) && x2 < half*(ymax_+cy)) ) {
            F[0] = loadMag*std::sin(loadAngle1);
            F[1] = loadMag*std::cos(loadAngle1);
          }
          break;
        }
        case 3: { // Force applied to lower right corner
          if ( (x1 > xmax_-cx) &&
               (x2 < cy) ) {
            F[0] = loadMag*std::sin(loadAngle1);
            F[1] = loadMag*std::cos(loadAngle1);
          }
          break;
        }
        case 4: { // Force applied to center of bottom boundary
          if ( (x1 > half*(xmax_-cx) && x1 < half*(xmax_+cx)) &&
               (x2 < cy) ) {
            F[0] = loadMag*std::sin(loadAngle1);
            F[1] = loadMag*std::cos(loadAngle1);
          }
          break;
        }
      }
    }
  }

  virtual void funcRHS(std::vector<Real> &F,
                 const std::vector<Real> &x,
                 const std::vector<Real> &param) const {
    Real zero(0), half(0.5), six(6), ten(10);
    Real loadMag = std::pow(ten,param[0]), loadAngle1 = M_PI*(param[1]/six - half); // loadAngle2 = param[2];
    F.clear(); F.resize(this->spaceDim_,zero);
    if (this->spaceDim_ == 2) {
      //Real eps(std::sqrt(ROL::ROL_EPSILON<Real>()));
      //Real cx(2*cx_+eps), cy(2*cy_+eps), half(0.5);
      Real cx(0.1), cy(0.05), half(0.5);
      Real x1 = x[0], x2 = x[1];
      switch(loadCase_) {
        case 0: { // Force applied to center of top boundary
          if ( x1 > half*(xmax_-cx) && x1 < half*(xmax_+cx) &&
               x2 > ymax_-cy ) {
            F[0] = loadMag*std::sin(loadAngle1);
            F[1] = loadMag*std::cos(loadAngle1);
          }
          break;
        }
        case 1: { // Force applied to top right corner
          if ( (x1 > xmax_-cx) &&
               (x2 > ymax_-cy) ) {
            F[0] = loadMag*std::sin(loadAngle1);
            F[1] = loadMag*std::cos(loadAngle1);
          }
          break;
        }
        case 2: { // Force applied to center of right boundary
          if ( (x1 > xmax_-cx) &&
               (x2 > half*(ymax_-cy) && x2 < half*(ymax_+cy)) ) {
            F[0] = loadMag*std::sin(loadAngle1);
            F[1] = loadMag*std::cos(loadAngle1);
          }
          break;
        }
        case 3: { // Force applied to lower right corner
          if ( (x1 > xmax_-cx) &&
               (x2 < cy) ) {
            F[0] = loadMag*std::sin(loadAngle1);
            F[1] = loadMag*std::cos(loadAngle1);
          }
          break;
        }
        case 4: { // Force applied to center of bottom boundary
          if ( (x1 > half*(xmax_-cx) && x1 < half*(xmax_+cx)) &&
               (x2 < cy) ) {
            F[0] = loadMag*std::sin(loadAngle1);
            F[1] = loadMag*std::cos(loadAngle1);
          }
          break;
        }
      }
    }
  }

}; // class ElasticitySIMPData

#endif
