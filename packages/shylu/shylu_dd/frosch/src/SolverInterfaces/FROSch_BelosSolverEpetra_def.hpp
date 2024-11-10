// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_BELOSSOLVEREPETRA_DEF_HPP
#define _FROSCH_BELOSSOLVEREPETRA_DEF_HPP

#include <FROSch_BelosSolverEpetra_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template<class SC,class LO,class GO,class NO>
    int BelosSolverEpetra<SC,LO,GO,NO>::initialize()
    {
        FROSCH_TIMER_START_SOLVER(initializeTime,"BelosSolverEpetra::initialize");
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int BelosSolverEpetra<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_SOLVER(computeTime,"BelosSolverEpetra::compute");
        FROSCH_ASSERT(this->IsInitialized_,"FROSch::BelosSolverEpetra: !this->IsInitialized_");
        this->IsComputed_ = true;
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    void BelosSolverEpetra<SC,LO,GO,NO>::apply(const XMultiVector &x,
                                               XMultiVector &y,
                                               ETransp mode,
                                               SC alpha,
                                               SC beta) const
    {
        FROSCH_TIMER_START_SOLVER(applyTime,"BelosSolverEpetra::apply");
        FROSCH_ASSERT(this->IsComputed_,"FROSch::BelosSolverEpetra: !this->IsComputed_.");

        const EpetraMultiVectorT<GO,NO> * xEpetraMultiVectorX = dynamic_cast<const EpetraMultiVectorT<GO,NO> *>(&x);
        FROSCH_ASSERT(xEpetraMultiVectorX,"FROSch::BelosSolverEpetra: dynamic_cast failed.");
        RCP<EMultiVector> epetraMultiVectorX = xEpetraMultiVectorX->getEpetra_MultiVector();

        if (Y_.is_null()) Y_ = XMultiVectorFactory::Build(y.getMap(),y.getNumVectors());
        EpetraMultiVectorT<GO,NO> * xEpetraMultiVectorY = dynamic_cast<EpetraMultiVectorT<GO,NO> *>(Y_.get());
        FROSCH_ASSERT(xEpetraMultiVectorY,"FROSch::BelosSolverEpetra: dynamic_cast failed.");
        RCP<EMultiVector> epetraMultiVectorY = xEpetraMultiVectorY->getEpetra_MultiVector();

        BelosLinearProblem_->setProblem(epetraMultiVectorY,epetraMultiVectorX);
        BelosSolver_->solve();

        y.update(alpha,*Y_,beta);
    }

    template<class SC,class LO,class GO,class NO>
    int BelosSolverEpetra<SC,LO,GO,NO>::updateMatrix(ConstXMatrixPtr k,
                                                     bool reuseInitialize)
    {
        FROSCH_TIMER_START_SOLVER(updateMatrixTime,"BelosSolverEpetra::updateMatrix");
        this->K_ = k;
        FROSCH_ASSERT(!this->K_.is_null(),"FROSch::BelosSolverEpetra: K_ is null.");

        const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*this->K_);
        const EpetraCrsMatrixT<GO,NO>& xEpetraMat = dynamic_cast<const EpetraCrsMatrixT<GO,NO>&>(*crsOp.getCrsMatrix());
        ECrsMatrixPtr epetraMat = xEpetraMat.getEpetra_CrsMatrixNonConst();
        TEUCHOS_TEST_FOR_EXCEPT(epetraMat.is_null());

        BelosLinearProblem_->setOperator(epetraMat);

        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    BelosSolverEpetra<SC,LO,GO,NO>::BelosSolverEpetra(ConstXMatrixPtr k,
                                                      ParameterListPtr parameterList,
                                                      string description) :
    Solver<SC,LO,GO,NO> (k,parameterList,description)
    {
        FROSCH_TIMER_START_SOLVER(BelosSolverEpetraTime,"BelosSolverEpetra::BelosSolverEpetra");
        FROSCH_ASSERT(!this->K_.is_null(),"FROSch::BelosSolverEpetra: K_ is null.");
        FROSCH_ASSERT(this->K_->getRowMap()->lib()==UseEpetra,"FROSch::BelosSolverEpetra: Not compatible with Tpetra.")

        const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*this->K_);
        const EpetraCrsMatrixT<GO,NO>& xEpetraMat = dynamic_cast<const EpetraCrsMatrixT<GO,NO>&>(*crsOp.getCrsMatrix());
        ECrsMatrixPtr epetraMat = xEpetraMat.getEpetra_CrsMatrixNonConst();
        TEUCHOS_TEST_FOR_EXCEPT(epetraMat.is_null());

        EMultiVectorPtr xTmp;
        EMultiVectorPtr bTmp;

        BelosLinearProblem_.reset(new BelosLinearProblem(epetraMat,xTmp,bTmp));

        ParameterListPtr belosParameterList = sublist(this->ParameterList_,"Belos");
        if (belosParameterList->isSublist(this->ParameterList_->get("Solver","GMRES"))) belosParameterList = sublist(belosParameterList,this->ParameterList_->get("Solver","GMRES"));

        BelosSolverFactory belosFactory;
        BelosSolver_ = belosFactory.create(this->ParameterList_->get("Solver","GMRES"),belosParameterList);
        BelosSolver_->setProblem(BelosLinearProblem_);
    }

}

#endif
