// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_BELOSSOLVERTPETRA_DEF_HPP
#define _FROSCH_BELOSSOLVERTPETRA_DEF_HPP

#include <FROSch_BelosSolverTpetra_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template<class SC,class LO,class GO,class NO>
    int BelosSolverTpetra<SC,LO,GO,NO>::initialize()
    {
        FROSCH_TIMER_START_SOLVER(initializeTime,"BelosSolverTpetra::initialize");
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int BelosSolverTpetra<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_SOLVER(computeTime,"BelosSolverTpetra::compute");
        FROSCH_ASSERT(this->IsInitialized_,"FROSch::BelosSolverTpetra: !this->IsInitialized_");
        this->IsComputed_ = true;
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    void BelosSolverTpetra<SC,LO,GO,NO>::apply(const XMultiVector &x,
                                               XMultiVector &y,
                                               ETransp mode,
                                               SC alpha,
                                               SC beta) const
    {
        FROSCH_TIMER_START_SOLVER(applyTime,"BelosSolverTpetra::apply");
        FROSCH_ASSERT(this->IsComputed_,"FROSch::BelosSolverTpetra: !this->IsComputed_.");

        const TpetraMultiVector<SC,LO,GO,NO> * xTpetraMultiVectorX = dynamic_cast<const TpetraMultiVector<SC,LO,GO,NO> *>(&x);
        FROSCH_ASSERT(xTpetraMultiVectorX,"FROSch::BelosSolverTpetra: dynamic_cast failed.");
        TMultiVectorPtr tpetraMultiVectorX = xTpetraMultiVectorX->getTpetra_MultiVector();

        if (Y_.is_null()) Y_ = XMultiVectorFactory::Build(y.getMap(),y.getNumVectors());
        const TpetraMultiVector<SC,LO,GO,NO> * xTpetraMultiVectorY = dynamic_cast<const TpetraMultiVector<SC,LO,GO,NO> *>(Y_.get());
        FROSCH_ASSERT(xTpetraMultiVectorY,"FROSch::BelosSolverTpetra: dynamic_cast failed.");
        TMultiVectorPtr tpetraMultiVectorY = xTpetraMultiVectorY->getTpetra_MultiVector();

        BelosLinearProblem_->setProblem(tpetraMultiVectorY,tpetraMultiVectorX);
        BelosSolver_->solve();

        y.update(alpha,*Y_,beta);
    }

    template<class SC,class LO,class GO,class NO>
    int BelosSolverTpetra<SC,LO,GO,NO>::updateMatrix(ConstXMatrixPtr k,
                                                     bool reuseInitialize)
    {
        FROSCH_TIMER_START_SOLVER(updateMatrixTime,"BelosSolverTpetra::updateMatrix");
        this->K_ = k;
        FROSCH_ASSERT(!this->K_.is_null(),"FROSch::BelosSolverTpetra: K_ is null.");

        const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*this->K_);
        const TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<const TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());
        ConstTCrsMatrixPtr tpetraMat = xTpetraMat.getTpetra_CrsMatrix();
        TEUCHOS_TEST_FOR_EXCEPT(tpetraMat.is_null());

        BelosLinearProblem_->setOperator(tpetraMat);

        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    BelosSolverTpetra<SC,LO,GO,NO>::BelosSolverTpetra(ConstXMatrixPtr k,
                                                      ParameterListPtr parameterList,
                                                      string description) :
    Solver<SC,LO,GO,NO> (k,parameterList,description)
    {
        FROSCH_TIMER_START_SOLVER(BelosSolverTpetraTime,"BelosSolverTpetra::BelosSolverTpetra");
        FROSCH_ASSERT(!this->K_.is_null(),"FROSch::BelosSolverTpetra: K_ is null.");
        FROSCH_ASSERT(this->K_->getRowMap()->lib()==UseTpetra,"FROSch::BelosSolverTpetra: Not compatible with Epetra.")

        const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*this->K_);
        const TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<const TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());
        ConstTCrsMatrixPtr tpetraMat = xTpetraMat.getTpetra_CrsMatrix();
        TEUCHOS_TEST_FOR_EXCEPT(tpetraMat.is_null());

        TMultiVectorPtr xTmp;
        TMultiVectorPtr bTmp;

        BelosLinearProblem_.reset(new BelosLinearProblem(tpetraMat,xTmp,bTmp));

        ParameterListPtr belosParameterList = sublist(this->ParameterList_,"Belos");
        if (belosParameterList->isSublist(this->ParameterList_->get("Solver","GMRES"))) belosParameterList = sublist(belosParameterList,this->ParameterList_->get("Solver","GMRES"));

        BelosSolverFactory belosFactory;
        BelosSolver_ = belosFactory.create(this->ParameterList_->get("Solver","GMRES"),belosParameterList);
        BelosSolver_->setProblem(BelosLinearProblem_);
    }

}

#endif
