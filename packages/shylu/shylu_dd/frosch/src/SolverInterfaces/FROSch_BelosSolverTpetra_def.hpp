//@HEADER
// ************************************************************************
//
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Alexander Heinlein (alexander.heinlein@uni-koeln.de)
//
// ************************************************************************
//@HEADER

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
