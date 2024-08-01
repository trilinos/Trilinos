// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_AMESOS2SOLVEREPETRA_DEF_HPP
#define _FROSCH_AMESOS2SOLVEREPETRA_DEF_HPP

#include <FROSch_Amesos2SolverEpetra_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template<class SC,class LO,class GO,class NO>
    Amesos2SolverEpetra<SC,LO,GO,NO>::~Amesos2SolverEpetra()
    {
        Amesos2Solver_.reset();
    }

    template<class SC,class LO,class GO,class NO>
    int Amesos2SolverEpetra<SC,LO,GO,NO>::initialize()
    {
        FROSCH_TIMER_START_SOLVER(initializeTime,"Amesos2SolverEpetra::initialize");
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        Amesos2Solver_->symbolicFactorization();
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int Amesos2SolverEpetra<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_SOLVER(computeTime,"Amesos2SolverEpetra::compute");
        FROSCH_ASSERT(this->IsInitialized_,"FROSch::Amesos2SolverEpetra: !this->IsInitialized_");
        this->IsComputed_ = true;
        Amesos2Solver_->numericFactorization();
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    void Amesos2SolverEpetra<SC,LO,GO,NO>::apply(const XMultiVector &x,
                                                 XMultiVector &y,
                                                 ETransp mode,
                                                 SC alpha,
                                                 SC beta) const
    {
        FROSCH_TIMER_START_SOLVER(applyTime,"Amesos2SolverEpetra::apply");
        FROSCH_ASSERT(this->IsComputed_,"FROSch::Amesos2SolverEpetra: !this->IsComputed_.");

        const EpetraMultiVectorT<GO,NO> * xEpetraMultiVectorX = dynamic_cast<const EpetraMultiVectorT<GO,NO> *>(&x);
        FROSCH_ASSERT(xEpetraMultiVectorX,"FROSch::Amesos2SolverEpetra: dynamic_cast failed.");
        RCP<EMultiVector> epetraMultiVectorX = xEpetraMultiVectorX->getEpetra_MultiVector();

        if (Y_.is_null()) Y_ = XMultiVectorFactory::Build(y.getMap(),x.getNumVectors());
        EpetraMultiVectorT<GO,NO> * xEpetraMultiVectorY = dynamic_cast<EpetraMultiVectorT<GO,NO> *>(Y_.get());
        FROSCH_ASSERT(xEpetraMultiVectorY,"FROSch::Amesos2SolverEpetra: dynamic_cast failed.");
        RCP<EMultiVector> epetraMultiVectorY = xEpetraMultiVectorY->getEpetra_MultiVector();

        Amesos2Solver_->setX(epetraMultiVectorY);
        Amesos2Solver_->setB(epetraMultiVectorX);

        FROSCH_ASSERT(mode==NO_TRANS,"FROSch::Amesos2SolverEpetra: mode!=NO_TRANS");
        Amesos2Solver_->solve(); // What about solving with transposed matrices?

        y.update(alpha,*Y_,beta);
    }

    template<class SC,class LO,class GO,class NO>
    int Amesos2SolverEpetra<SC,LO,GO,NO>::updateMatrix(ConstXMatrixPtr k,
                                                      bool reuseInitialize)
    {
        FROSCH_TIMER_START_SOLVER(updateMatrixTime,"Amesos2SolverEpetra::updateMatrix");
        this->K_ = k;
        FROSCH_ASSERT(!this->K_.is_null(),"FROSch::Amesos2SolverEpetra: K_ is null.");

        const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*this->K_);
        const EpetraCrsMatrixT<GO,NO>& xEpetraMat = dynamic_cast<const EpetraCrsMatrixT<GO,NO>&>(*crsOp.getCrsMatrix());
        ECrsMatrixPtr epetraMat = xEpetraMat.getEpetra_CrsMatrixNonConst();
        TEUCHOS_TEST_FOR_EXCEPT(epetraMat.is_null());

        if (reuseInitialize) {
            Amesos2Solver_->setA(epetraMat,Amesos2::SYMBFACT);
        } else {
            Amesos2Solver_->setA(epetraMat,Amesos2::CLEAN);
        }
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    Amesos2SolverEpetra<SC,LO,GO,NO>::Amesos2SolverEpetra(ConstXMatrixPtr k,
                                                          ParameterListPtr parameterList,
                                                          string description) :
    Solver<SC,LO,GO,NO> (k,parameterList,description)
    {
        FROSCH_TIMER_START_SOLVER(Amesos2SolverEpetraTime,"Amesos2SolverEpetra::Amesos2SolverEpetra");
        FROSCH_ASSERT(!this->K_.is_null(),"FROSch::Amesos2SolverEpetra: K_ is null.");
        FROSCH_ASSERT(this->K_->getRowMap()->lib()==UseEpetra,"FROSch::Amesos2SolverEpetra: Not compatible with Tpetra.")

        const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*this->K_);
        const EpetraCrsMatrixT<GO,NO>& xEpetraMat = dynamic_cast<const EpetraCrsMatrixT<GO,NO>&>(*crsOp.getCrsMatrix());
        ECrsMatrixPtr epetraMat = xEpetraMat.getEpetra_CrsMatrixNonConst();
        TEUCHOS_TEST_FOR_EXCEPT(epetraMat.is_null());

        EMultiVectorPtr xTmp;
        EMultiVectorPtr bTmp;

        ParameterListPtr amesos2ParameterList = sublist(this->ParameterList_,"Amesos2");
        if (amesos2ParameterList->isSublist(this->ParameterList_->get("Solver","Klu"))) amesos2ParameterList = sublist(amesos2ParameterList,this->ParameterList_->get("Solver","Klu"));
        amesos2ParameterList->setName("Amesos2");

        Amesos2Solver_ = Amesos2::create<ECrsMatrix,EMultiVector>(this->ParameterList_->get("Solver","Klu"),epetraMat,xTmp,bTmp);
        Amesos2Solver_->setParameters(amesos2ParameterList);
    }

}

#endif
