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
