// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_AMESOSSOLVEREPETRA_DEF_HPP
#define _FROSCH_AMESOSSOLVEREPETRA_DEF_HPP

#include <FROSch_AmesosSolverEpetra_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template<class SC,class LO,class GO,class NO>
    AmesosSolverEpetra<SC,LO,GO,NO>::~AmesosSolverEpetra()
    {
        AmesosSolver_.reset();
        EpetraLinearProblem_.reset();
    }

    template<class SC,class LO,class GO,class NO>
    int AmesosSolverEpetra<SC,LO,GO,NO>::initialize()
    {
        FROSCH_TIMER_START_SOLVER(initializeTime,"AmesosSolverEpetra::initialize");
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        AMESOS_CHK_ERR(AmesosSolver_->SymbolicFactorization());
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int AmesosSolverEpetra<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_SOLVER(computeTime,"AmesosSolverEpetra::compute");
        FROSCH_ASSERT(this->IsInitialized_,"FROSch::AmesosSolverEpetra: !this->IsInitialized_");
        this->IsComputed_ = true;
        AMESOS_CHK_ERR(AmesosSolver_->NumericFactorization());
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    void AmesosSolverEpetra<SC,LO,GO,NO>::apply(const XMultiVector &x,
                                                XMultiVector &y,
                                                ETransp mode,
                                                SC alpha,
                                                SC beta) const
    {
        FROSCH_TIMER_START_SOLVER(applyTime,"AmesosSolverEpetra::apply");
        FROSCH_ASSERT(this->IsComputed_,"FROSch::AmesosSolverEpetra: !this->IsComputed_.");

        const EpetraMultiVectorT<GO,NO> * xEpetraMultiVectorX = dynamic_cast<const EpetraMultiVectorT<GO,NO> *>(&x);
        FROSCH_ASSERT(xEpetraMultiVectorX,"FROSch::AmesosSolverEpetra: dynamic_cast failed.");
        RCP<EMultiVector> epetraMultiVectorX = xEpetraMultiVectorX->getEpetra_MultiVector();

        if (Y_.is_null()) Y_ = XMultiVectorFactory::Build(y.getMap(),x.getNumVectors());
        EpetraMultiVectorT<GO,NO> * xEpetraMultiVectorY = dynamic_cast<EpetraMultiVectorT<GO,NO> *>(Y_.get());
        FROSCH_ASSERT(xEpetraMultiVectorY,"FROSch::AmesosSolverEpetra: dynamic_cast failed.");
        RCP<EMultiVector> epetraMultiVectorY = xEpetraMultiVectorY->getEpetra_MultiVector();

        EpetraLinearProblem_->SetLHS(epetraMultiVectorY.get());
        EpetraLinearProblem_->SetRHS(epetraMultiVectorX.get());

        FROSCH_ASSERT(mode==NO_TRANS,"FROSch::AmesosSolverEpetra: mode!=NO_TRANS");
        EpetraLinearProblem_->GetMatrix()->SetUseTranspose(mode==TRANS);
        AmesosSolver_->Solve();

        y.update(alpha,*Y_,beta);
    }

    template<class SC,class LO,class GO,class NO>
    int AmesosSolverEpetra<SC,LO,GO,NO>::updateMatrix(ConstXMatrixPtr k,
                                                      bool reuseInitialize)
    {
        FROSCH_ASSERT(false,"FROSch::AmesosSolverEpetra: updateMatrix() is not implemented for the AmesosSolverEpetra yet.");
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    AmesosSolverEpetra<SC,LO,GO,NO>::AmesosSolverEpetra(ConstXMatrixPtr k,
                                                          ParameterListPtr parameterList,
                                                          string description) :
    Solver<SC,LO,GO,NO> (k,parameterList,description)
    {
        FROSCH_TIMER_START_SOLVER(AmesosSolverEpetraTime,"AmesosSolverEpetra::AmesosSolverEpetra");
        FROSCH_ASSERT(!this->K_.is_null(),"FROSch::AmesosSolverEpetra: K_ is null.");
        FROSCH_ASSERT(this->K_->getRowMap()->lib()==UseEpetra,"FROSch::AmesosSolverEpetra: Not compatible with Tpetra.")

        const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*this->K_);
        const EpetraCrsMatrixT<GO,NO>& xEpetraMat = dynamic_cast<const EpetraCrsMatrixT<GO,NO>&>(*crsOp.getCrsMatrix());
        ECrsMatrixPtr epetraMat = xEpetraMat.getEpetra_CrsMatrixNonConst();
        TEUCHOS_TEST_FOR_EXCEPT(epetraMat.is_null());

        EMultiVectorPtr xTmp;
        EMultiVectorPtr bTmp;

        EpetraLinearProblem_.reset(new ELinearProblem(epetraMat.get(),xTmp.get(),bTmp.get()));

        Amesos amesosFactory;
        AmesosSolver_.reset(amesosFactory.Create(this->ParameterList_->get("Solver","Klu"),*EpetraLinearProblem_));

        ParameterListPtr amesosParameterList = sublist(this->ParameterList_,"Amesos");
        if (amesosParameterList->isSublist(this->ParameterList_->get("Solver","Klu"))) amesosParameterList = sublist(amesosParameterList,this->ParameterList_->get("Solver","Klu"));
        AmesosSolver_->SetParameters(*amesosParameterList);
    }

}

#endif
