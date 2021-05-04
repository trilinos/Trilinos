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

#ifndef _FROSCH_IFPACK2PRECONDITIONERTPETRA_DEF_HPP
#define _FROSCH_IFPACK2PRECONDITIONERTPETRA_DEF_HPP

#include <FROSch_Ifpack2PreconditionerTpetra_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template<class SC,class LO,class GO,class NO>
    int Ifpack2PreconditionerTpetra<SC,LO,GO,NO>::initialize()
    {
        FROSCH_TIMER_START_SOLVER(initializeTime,"Ifpack2PreconditionerTpetra::initialize");
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        Ifpack2Preconditioner_->initialize();
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int Ifpack2PreconditionerTpetra<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_SOLVER(computeTime,"Ifpack2PreconditionerTpetra::compute");
        FROSCH_ASSERT(this->IsInitialized_,"FROSch::Ifpack2PreconditionerTpetra: !this->IsInitialized_");
        this->IsComputed_ = true;
        Ifpack2Preconditioner_->compute();
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    void Ifpack2PreconditionerTpetra<SC,LO,GO,NO>::apply(const XMultiVector &x,
                                                         XMultiVector &y,
                                                         ETransp mode,
                                                         SC alpha,
                                                         SC beta) const
    {
        FROSCH_TIMER_START_SOLVER(applyTime,"Ifpack2PreconditionerTpetra::apply");
        FROSCH_ASSERT(this->IsComputed_,"FROSch::Ifpack2PreconditionerTpetra: !this->IsComputed_.");

        const TpetraMultiVector<SC,LO,GO,NO> * xTpetraMultiVectorX = dynamic_cast<const TpetraMultiVector<SC,LO,GO,NO> *>(&x);
        FROSCH_ASSERT(xTpetraMultiVectorX,"FROSch::Ifpack2PreconditionerTpetra: dynamic_cast failed.");
        TMultiVectorPtr tpetraMultiVectorX = xTpetraMultiVectorX->getTpetra_MultiVector();

        const TpetraMultiVector<SC,LO,GO,NO> * xTpetraMultiVectorY = dynamic_cast<const TpetraMultiVector<SC,LO,GO,NO> *>(&y);
        FROSCH_ASSERT(xTpetraMultiVectorY,"FROSch::Ifpack2PreconditionerTpetra: dynamic_cast failed.");
        TMultiVectorPtr tpetraMultiVectorY = xTpetraMultiVectorY->getTpetra_MultiVector();

        Ifpack2Preconditioner_->apply(*tpetraMultiVectorX,*tpetraMultiVectorY,mode,alpha,beta);
    }

    template<class SC,class LO,class GO,class NO>
    int Ifpack2PreconditionerTpetra<SC,LO,GO,NO>::updateMatrix(ConstXMatrixPtr k,
                                                               bool reuseInitialize)
    {
        FROSCH_ASSERT(false,"FROSch::Ifpack2PreconditionerTpetra: updateMatrix() is not implemented for the Ifpack2PreconditionerTpetra yet.");
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    Ifpack2PreconditionerTpetra<SC,LO,GO,NO>::Ifpack2PreconditionerTpetra(ConstXMatrixPtr k,
                                                                          ParameterListPtr parameterList,
                                                                          string description) :
    Solver<SC,LO,GO,NO> (k,parameterList,description)
    {
        FROSCH_TIMER_START_SOLVER(Ifpack2PreconditionerTpetraTime,"Ifpack2PreconditionerTpetra::Ifpack2PreconditionerTpetra");
        FROSCH_ASSERT(!this->K_.is_null(),"FROSch::Ifpack2PreconditionerTpetra: K_ is null.");
        FROSCH_ASSERT(this->K_->getRowMap()->lib()==UseTpetra,"FROSch::Ifpack2PreconditionerTpetra: Not compatible with Epetra.")

        const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*this->K_);
        const TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<const TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());
        ConstTCrsMatrixPtr tpetraMat = xTpetraMat.getTpetra_CrsMatrix();
        TEUCHOS_TEST_FOR_EXCEPT(tpetraMat.is_null());

        Ifpack2::Details::OneLevelFactory<TRowMatrix> ifpack2Factory;
        Ifpack2Preconditioner_ = ifpack2Factory.create(this->ParameterList_->get("Solver","RILUK"),tpetraMat);

        ParameterListPtr ifpack2ParameterList = sublist(this->ParameterList_,"Ifpack2");
        if (ifpack2ParameterList->isSublist(this->ParameterList_->get("Solver","RILUK"))) ifpack2ParameterList = sublist(ifpack2ParameterList,this->ParameterList_->get("Solver","RILUK"));
        ifpack2ParameterList->setName("Ifpack2");
        Ifpack2Preconditioner_->setParameters(*ifpack2ParameterList);
    }

}

#endif
