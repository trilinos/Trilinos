// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_AMESOS2SOLVERTPETRA_DEF_HPP
#define _FROSCH_AMESOS2SOLVERTPETRA_DEF_HPP

#include <FROSch_Amesos2SolverTpetra_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template<class SC,class LO,class GO,class NO>
    Amesos2SolverTpetra<SC,LO,GO,NO>::~Amesos2SolverTpetra()
    {
        Amesos2Solver_.reset();
    }

    template<class SC,class LO,class GO,class NO>
    int Amesos2SolverTpetra<SC,LO,GO,NO>::initialize()
    {
        FROSCH_TIMER_START_SOLVER(initializeTime,"Amesos2SolverTpetra::initialize");
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        Amesos2Solver_->symbolicFactorization();
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int Amesos2SolverTpetra<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_SOLVER(computeTime,"Amesos2SolverTpetra::compute");
        FROSCH_ASSERT(this->IsInitialized_,"FROSch::Amesos2SolverTpetra: !this->IsInitialized_");
        this->IsComputed_ = true;
        Amesos2Solver_->numericFactorization();
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    void Amesos2SolverTpetra<SC,LO,GO,NO>::apply(const XMultiVector &x,
                                                 XMultiVector &y,
                                                 ETransp mode,
                                                 SC alpha,
                                                 SC beta) const
    {
        FROSCH_TIMER_START_SOLVER(applyTime,"Amesos2SolverTpetra::apply");
        FROSCH_ASSERT(this->IsComputed_,"FROSch::Amesos2SolverTpetra: !this->IsComputed_.");

        const TpetraMultiVector<SC,LO,GO,NO> * xTpetraMultiVectorX = dynamic_cast<const TpetraMultiVector<SC,LO,GO,NO> *>(&x);
        FROSCH_ASSERT(xTpetraMultiVectorX,"FROSch::Amesos2SolverTpetra: dynamic_cast failed.");
        TMultiVectorPtr tpetraMultiVectorX = xTpetraMultiVectorX->getTpetra_MultiVector();

        if (Y_.is_null() || Y_->getNumVectors() != x.getNumVectors()) {
            Y_ = XMultiVectorFactory::Build(y.getMap(),x.getNumVectors());
        }
        const TpetraMultiVector<SC,LO,GO,NO> * xTpetraMultiVectorY = dynamic_cast<const TpetraMultiVector<SC,LO,GO,NO> *>(Y_.get());
        FROSCH_ASSERT(xTpetraMultiVectorY,"FROSch::Amesos2SolverTpetra: dynamic_cast failed.");
        TMultiVectorPtr tpetraMultiVectorY = xTpetraMultiVectorY->getTpetra_MultiVector();

        Amesos2Solver_->setX(tpetraMultiVectorY);
        Amesos2Solver_->setB(tpetraMultiVectorX);

        FROSCH_ASSERT(mode==NO_TRANS,"FROSch::Amesos2SolverTpetra: mode!=NO_TRANS");
        Amesos2Solver_->solve(); // What about solving with transposed matrices?

        y.update(alpha,*Y_,beta);
    }

    template<class SC,class LO,class GO,class NO>
    int Amesos2SolverTpetra<SC,LO,GO,NO>::updateMatrix(ConstXMatrixPtr k,
                                                       bool reuseInitialize)
    {
        FROSCH_TIMER_START_SOLVER(updateMatrixTime,"Amesos2SolverTpetra::updateMatrix");
        this->K_ = k;
        FROSCH_ASSERT(!this->K_.is_null(),"FROSch::Amesos2SolverTpetra: K_ is null.");

        const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*this->K_);
        const TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<const TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());
        ConstTCrsMatrixPtr tpetraMat = xTpetraMat.getTpetra_CrsMatrix();
        TEUCHOS_TEST_FOR_EXCEPT(tpetraMat.is_null());

        if (reuseInitialize) {
            Amesos2Solver_->setA(tpetraMat,Amesos2::SYMBFACT);
        } else {
            Amesos2Solver_->setA(tpetraMat,Amesos2::CLEAN);
        }
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    Amesos2SolverTpetra<SC,LO,GO,NO>::Amesos2SolverTpetra(ConstXMatrixPtr k,
                                                          ParameterListPtr parameterList,
                                                          string description) :
    Solver<SC,LO,GO,NO> (k,parameterList,description)
    {
        FROSCH_TIMER_START_SOLVER(Amesos2SolverTpetraTime,"Amesos2SolverTpetra::Amesos2SolverTpetra");
        FROSCH_ASSERT(!this->K_.is_null(),"FROSch::Amesos2SolverTpetra: K_ is null.");
        FROSCH_ASSERT(this->K_->getRowMap()->lib()==UseTpetra,"FROSch::Amesos2SolverTpetra: Not compatible with Epetra.")

        const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*this->K_);
        const TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<const TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());
        ConstTCrsMatrixPtr tpetraMat = xTpetraMat.getTpetra_CrsMatrix();
        TEUCHOS_TEST_FOR_EXCEPT(tpetraMat.is_null());

        TMultiVectorPtr xTmp;
        TMultiVectorPtr bTmp;

        ParameterListPtr amesos2ParameterList = sublist(this->ParameterList_,"Amesos2");
        if (amesos2ParameterList->isSublist(this->ParameterList_->get("Solver","Klu"))) amesos2ParameterList = sublist(amesos2ParameterList,this->ParameterList_->get("Solver","Klu"));
        amesos2ParameterList->setName("Amesos2");

        Amesos2Solver_ = Amesos2::create<TCrsMatrix,TMultiVector>(this->ParameterList_->get("Solver","Klu"),tpetraMat,xTmp,bTmp);
        Amesos2Solver_->setParameters(amesos2ParameterList);
    }

}

#endif
