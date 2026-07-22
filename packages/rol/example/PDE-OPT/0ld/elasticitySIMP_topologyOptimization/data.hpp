// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PDEOPT_ELASTICITYSIMP_OPERATORS_H
#define ROL_PDEOPT_ELASTICITYSIMP_OPERATORS_H

#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "../TOOLS/elasticitySIMP.hpp"

template<class Real>
class ElasticitySIMPOperators : public ElasticitySIMP <Real> {
private:
  using GO = typename Tpetra::Map<>::global_ordinal_type;

  Real volFrac_;

public:

  ElasticitySIMPOperators(const ROL::Ptr<const Teuchos::Comm<int> > &comm,
                          const Teuchos::RCP<Teuchos::ParameterList> &parlist,
                          const ROL::Ptr<std::ostream> &outStream) {
    ElasticitySIMPOperators_Initialize(comm, parlist, outStream);
    this->SetSIMPParallelStructure();
    this->SetUpLocalIntrepidArrays();
    this->ComputeLocalSystemMats(true);
//Setup DBC information, do not specify any bc sides, use coordinates to determine the BC instead
    std::vector<GO> dbc_side {};
    this->SetUpMyDBCInfo(true, dbc_side);
    this->process_loading_information(parlist);
//With new modification on boundary traction, ComputeLocalForceVec should go after SetUpMyBCInfo and process_loading_information
    //this->AssembleSystemGraph();
    this->AssembleSystemGraph();
    this->AssembleSystemMats();
    this->ComputeLocalForceVec();
    this->AssembleRHSVector();
    this->EnforceDBC();
    this->ConstructSolvers();
  }

  void ElasticitySIMPOperators_Initialize(const ROL::Ptr<const Teuchos::Comm<int> > &comm,
                                          const Teuchos::RCP<Teuchos::ParameterList> &parlist,
                                          const ROL::Ptr<std::ostream> &outStream) {
    ElasticitySIMP<Real>::ElasticitySIMP_Initialize(comm, parlist, outStream);
    volFrac_           = parlist->sublist("ElasticityTopoOpt").get<Real>("Volume Fraction");
    this->initDensity_ = volFrac_;
  }

  // construct solvers with new material
  void constructSolverWithNewMaterial(void) {
    bool ifInit = false;
    {
    #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*SolverUpdateTime_example_PDEOPT_TOOLS_PDEFEM_GLOB);
    #endif
    this->ComputeLocalSystemMats(ifInit);
    this->AssembleSystemMats();
    this->MatrixRemoveDBC();
    }
    this->ConstructSolvers();
  }
  //

  ROL::Ptr<Intrepid::FieldContainer<int> > getPosCell(void) const {
    return this->meshMgr_->getPosCell();
  }

  void ApplyMatAToVec(ROL::Ptr<Tpetra::MultiVector<> > &Ju,
                const ROL::Ptr<const Tpetra::MultiVector<> > &u) {
    this->matA_dirichlet_->apply(*u, *Ju);
  }

  void ApplyJacobian1ToVec(ROL::Ptr<Tpetra::MultiVector<> > &Ju,
                     const ROL::Ptr<const Tpetra::MultiVector<> > &u) {
    #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ConstraintDerivativeTime_example_PDEOPT_TOOLS_PDEFEM_GLOB);
    #endif
    this->matA_dirichlet_->apply(*u, *Ju);
    // Jv is myUniqueMap_
    // u is myUniqueMap_
    // v should be myCellMap_
    //
    /*
    ROL::Ptr<Tpetra::MultiVector<> > Ju_overlap = ROL::makePtr<Tpetra::MultiVector<>>(this->myOverlapMap_, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > u_overlap = ROL::makePtr<Tpetra::MultiVector<>>(this->myOverlapMap_, 1, true);
    Tpetra::Import<> importer(u->getMap(), u_overlap->getMap());
    u_overlap->doImport(*u, importer, Tpetra::REPLACE);
    //simply apply matrix to a vector, no need to apply EBC to the vec
    //this-> ApplyBCToVec (u_overlap);
    Teuchos::ArrayRCP<const Real> uData = u_overlap->get1dView();
    
    Real SIMPScale(0), sum(0), zero(0);
    Intrepid::FieldContainer<Real> localDisp(this->numLocalDofs_);
    
    for(int i=0; i<this->numCells_; i++) {
      SIMPScale = this->SIMPmaterial_[i]->getSIMPScaleFactor();
      localDisp.initialize(zero);
      for(int j=0; j<this->numLocalDofs_; j++) {
        localDisp(j) = uData[u_overlap->getMap()->getLocalElement(this->cellDofs_(this->myCellIds_[i], j))];
      }
      for(int j=0; j<this->numLocalDofs_; j++) {	
        if(this->check_myGlobalDof_on_boundary(this->cellDofs_(this->myCellIds_[i],j))) {
          Ju_overlap -> replaceGlobalValue(this->cellDofs_(this->myCellIds_[i], j), 0, localDisp(j));
        }
        else {
          sum = zero;
          for(int k=0; k<this->numLocalDofs_; k++) {
            if( !this->check_myGlobalDof_on_boundary(this->cellDofs_(this->myCellIds_[i],k)) ) {
              sum += SIMPScale * (*this->gradgradMats0_)(i, j, k) * localDisp(k);
            }
          }
          Ju_overlap -> sumIntoGlobalValue(this->cellDofs_(this->myCellIds_[i], j), 0, sum);	
        }
      }
    }
    
    Tpetra::Export<> exporter(Ju_overlap->getMap(), Ju->getMap());
    //Ju->doExport(*Ju_overlap, exporter, Tpetra::REPLACE);
    Ju->doExport(*Ju_overlap, exporter, Tpetra::ADD);
    */
  }

  void ApplyInverseJacobian1ToVec(ROL::Ptr<Tpetra::MultiVector<> > &InvJu,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const bool ifTrans) {
    //Ju is matA->getDomainMap();
    #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*LUSubstitutionTime_example_PDEOPT_TOOLS_PDEFEM_GLOB);
    #endif
    this->getSolver(ifTrans)->setB(u);
    this->getSolver(ifTrans)->setX(InvJu);
    this->getSolver(ifTrans)->solve();
  }

  void ApplyJacobian2ToVec(ROL::Ptr<Tpetra::MultiVector<> > &Jv,
                     const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                     const ROL::Ptr<const Tpetra::MultiVector<> > &v) {
    // Jv is myUniqueMap_
    // u is myUniqueMap_
    // v should be myCellMap_
    //
    #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ConstraintDerivativeTime_example_PDEOPT_TOOLS_PDEFEM_GLOB);
    #endif
    ROL::Ptr<Tpetra::MultiVector<> > Jv_overlap = ROL::makePtr<Tpetra::MultiVector<>>(this->myOverlapMap_, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > u_overlap = ROL::makePtr<Tpetra::MultiVector<>>(this->myOverlapMap_, 1, true);
    Tpetra::Import<> importer(u->getMap(), u_overlap->getMap());
    u_overlap->doImport(*u, importer, Tpetra::REPLACE);
    //apply BC here, KU
    this->ApplyBCToVec (u_overlap);
    Teuchos::ArrayRCP<const Real> uData = u_overlap->get1dView();

    ROL::Ptr<Tpetra::MultiVector<> > v_local = ROL::makePtr<Tpetra::MultiVector<>>(this->myCellMap_, 1, true);
    Tpetra::Export<> exporter1(v->getMap(), this->myCellMap_);
    v_local->doExport(*v, exporter1, Tpetra::REPLACE);
    Teuchos::ArrayRCP<const Real> vData = v_local->get1dView();

    Real SIMPDScale(0), localV(0), sum(0), zero(0);
    Intrepid::FieldContainer<Real> localDisp(this->numLocalDofs_);

    for(int i=0; i<this->numCells_; i++) {
      SIMPDScale = this->SIMPmaterial_[i]->getSIMPFirstDerivativeScaleFactor();
      localV = vData[v_local->getMap()->getLocalElement(this->myCellIds_[i])];

      localDisp.initialize(zero);
      for(int j=0; j<this->numLocalDofs_; j++) {
        localDisp(j) = uData[u_overlap->getMap()->getLocalElement(this->cellDofs_(this->myCellIds_[i], j))];
      }

      for(int j=0; j<this->numLocalDofs_; j++) {
        if(this->check_myGlobalDof_on_boundary(this->cellDofs_(this->myCellIds_[i],j))) {
          Jv_overlap -> replaceGlobalValue(this->cellDofs_(this->myCellIds_[i], j), 0, localDisp(j));
        }
        else {
          sum = zero;
          for(int k=0; k<this->numLocalDofs_; k++) {
            if( !this->check_myGlobalDof_on_boundary(this->cellDofs_(this->myCellIds_[i],k)) ) {
              sum += localV * SIMPDScale * (*this->gradgradMats0_)(i, j, k) * localDisp(k);
            }
          }
          Jv_overlap -> sumIntoGlobalValue(this->cellDofs_(this->myCellIds_[i], j), 0, sum);
        }
      }
    }
    Tpetra::Export<> exporter2(Jv_overlap->getMap(), Jv->getMap());
    //Jv->doExport(*Jv_overlap, exporter2, Tpetra::REPLACE);
    Jv->doExport(*Jv_overlap, exporter2, Tpetra::ADD);
  }

  void ApplyAdjointJacobian2ToVec (ROL::Ptr<Tpetra::MultiVector<> > &Jv,
                             const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                             const ROL::Ptr<const Tpetra::MultiVector<> > &v) {
    // Jv is myCellMap_
    // u should be myUniqueMap_
    // v should be myUniqueMap_
    //
    #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ConstraintDerivativeTime_example_PDEOPT_TOOLS_PDEFEM_GLOB);
    #endif
    ROL::Ptr<Tpetra::MultiVector<> > u_overlap = ROL::makePtr<Tpetra::MultiVector<>>(this->myOverlapMap_, 1, true);
    Tpetra::Import<> importer1(u->getMap(), u_overlap->getMap());
    u_overlap->doImport(*u, importer1, Tpetra::REPLACE);
    //only apply BC to u
    this->ApplyBCToVec (u_overlap);
    Teuchos::ArrayRCP<const Real> uData = u_overlap->get1dView();

    ROL::Ptr<Tpetra::MultiVector<> > v_overlap = ROL::makePtr<Tpetra::MultiVector<>>(this->myOverlapMap_, 1, true);
    Tpetra::Import<> importer2(v->getMap(), v_overlap->getMap());
    v_overlap->doImport(*v, importer2, Tpetra::REPLACE);
    //ApplyBCToVec (v_overlap);
    Teuchos::ArrayRCP<const Real> vData = v_overlap->get1dView();

    Intrepid::FieldContainer<Real> u_local(this->numLocalDofs_);
    Intrepid::FieldContainer<Real> v_local(this->numLocalDofs_);

    Real SIMPDScale(0), sum(0), zero(0);
    for(int i=0; i<this->numCells_; i++) {
      SIMPDScale = this->SIMPmaterial_[i]->getSIMPFirstDerivativeScaleFactor();
      u_local.initialize(zero);
      v_local.initialize(zero);

      for(int j=0; j<this->numLocalDofs_; j++) {
        u_local(j) = uData[u_overlap->getMap()->getLocalElement(this->cellDofs_(this->myCellIds_[i], j))];
        v_local(j) = vData[v_overlap->getMap()->getLocalElement(this->cellDofs_(this->myCellIds_[i], j))];
      }

      sum = zero;
      for(int j=0; j<this->numLocalDofs_; j++) {
        if(this->check_myGlobalDof_on_boundary(this->cellDofs_(this->myCellIds_[i],j))) {
          sum += u_local(j) * v_local(j);
        }
        else {
          for(int k=0; k<this->numLocalDofs_; k++) {
            if( !this->check_myGlobalDof_on_boundary(this->cellDofs_(this->myCellIds_[i],k)) ) {
              sum += SIMPDScale * (*this->gradgradMats0_)(i, j, k) * v_local(j) * u_local(k);
            }
          }
        }
      }
      //put value into myDensityDerivative_
      Jv -> replaceGlobalValue(this->myCellIds_[i], 0, sum);
    }
  }

  void ApplyAdjointHessian22ToVec(ROL::Ptr<Tpetra::MultiVector<> > &Hvw,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &v,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &w) {
    // Hvw is myCellMap_
    // u should be myUniqueMap_
    // v should be myCellMap_
    // w shluld be myUniqueMapi_
    //
    #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ConstraintDerivativeTime_example_PDEOPT_TOOLS_PDEFEM_GLOB);
    #endif
    ROL::Ptr<Tpetra::MultiVector<> > u_overlap = ROL::makePtr<Tpetra::MultiVector<>>(this->myOverlapMap_, 1, true);
    Tpetra::Import<> importer1(u->getMap(), u_overlap->getMap());
    u_overlap->doImport(*u, importer1, Tpetra::REPLACE);
    //only apply BC to U
    this->ApplyBCToVec (u_overlap);
    Teuchos::ArrayRCP<const Real> uData = u_overlap->get1dView();

    Tpetra::Export<> exporter(v->getMap(), this->myCellMap_);
    ROL::Ptr<Tpetra::MultiVector<> > v_local = ROL::makePtr<Tpetra::MultiVector<>>(this->myCellMap_, 1, true);
    v_local->doExport(*v, exporter, Tpetra::REPLACE);
    Teuchos::ArrayRCP<const Real> vData = v_local->get1dView();

    ROL::Ptr<Tpetra::MultiVector<> > w_overlap = ROL::makePtr<Tpetra::MultiVector<>>(this->myOverlapMap_, 1, true);
    Tpetra::Import<> importer2(w->getMap(), w_overlap->getMap());
    w_overlap->doImport(*w, importer2, Tpetra::REPLACE);
    //ApplyBCToVec (w_overlap);
    Teuchos::ArrayRCP<const Real> wData = w_overlap->get1dView();

    Intrepid::FieldContainer<Real> u_local(this->numLocalDofs_);
    Intrepid::FieldContainer<Real> w_local(this->numLocalDofs_);

    Real SIMPDDScale(0), localV(0), sum(0), zero(0);

    for(int i=0; i<this->numCells_; i++) {
      SIMPDDScale = this->SIMPmaterial_[i]->getSIMPSecondDerivativeScaleFactor();
      localV = vData[v_local->getMap()->getLocalElement(this->myCellIds_[i])];

      u_local.initialize(zero);
      w_local.initialize(zero);

      for(int j=0; j<this->numLocalDofs_; j++) {
        u_local(j) = uData[u_overlap->getMap()->getLocalElement(this->cellDofs_(this->myCellIds_[i], j))];
        w_local(j) = wData[w_overlap->getMap()->getLocalElement(this->cellDofs_(this->myCellIds_[i], j))];
      }

      sum = zero;
      for(int j=0; j<this->numLocalDofs_; j++) {
        if(this->check_myGlobalDof_on_boundary(this->cellDofs_(this->myCellIds_[i],j))) {
          sum += u_local(j) * w_local(j);
        }
        else {
          for(int k=0; k<this->numLocalDofs_; k++) {
            if( !this->check_myGlobalDof_on_boundary(this->cellDofs_(this->myCellIds_[i],k)) ) {
              sum += localV * SIMPDDScale * (*this->gradgradMats0_)(i, j, k) * w_local(j) * u_local(k);
            }
          }
        }
      }
      //put value into myDensityDerivative_
      Hvw -> replaceGlobalValue(this->myCellIds_[i], 0, sum);
    }
  }

}; // class ElasticitySIMPOperators

#endif
