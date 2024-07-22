// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_IFPACK2PRECONDITIONERTPETRA_DEF_HPP
#define _FROSCH_IFPACK2PRECONDITIONERTPETRA_DEF_HPP

#include <FROSch_Ifpack2PreconditionerTpetra_decl.hpp>

#include "Ifpack2_Details_getCrsMatrix.hpp"
#include "Ifpack2_RILUK_decl.hpp"

#ifdef HAVE_SHYLU_DDFROSCH_ZOLTAN2
#include "Zoltan2_TpetraRowGraphAdapter.hpp"
#include "Zoltan2_OrderingProblem.hpp"
#include "Zoltan2_OrderingSolution.hpp"
#endif

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

#ifdef HAVE_SHYLU_DDFROSCH_ZOLTAN2
        if (this->useZoltan2 && this->useRILUK) {
            auto A = Ifpack2Preconditioner_->getMatrix();
            auto filteredA = rcp_dynamic_cast<const TRowMatrixFilterType> (A);

            // make copies of X & Y
            TMultiVector ReorderedX (*tpetraMultiVectorX, Teuchos::Copy);
            TMultiVector ReorderedY (*tpetraMultiVectorY, Teuchos::Copy);

            // permute X & Y
            filteredA->permuteOriginalToReordered (*tpetraMultiVectorX, ReorderedX);
            filteredA->permuteOriginalToReordered (*tpetraMultiVectorY, ReorderedY);

            // solve
            Ifpack2Preconditioner_->apply(ReorderedX,ReorderedY,mode,alpha,beta);

            // permute X back
            filteredA->permuteReorderedToOriginal (ReorderedY, *tpetraMultiVectorY);
        } else
#endif
        {
            Ifpack2Preconditioner_->apply(*tpetraMultiVectorX,*tpetraMultiVectorY,mode,alpha,beta);
        }
    }

    template<class SC,class LO,class GO,class NO>
    int Ifpack2PreconditionerTpetra<SC,LO,GO,NO>::updateMatrix(ConstXMatrixPtr k,
                                                               bool reuseInitialize)
    {
        if (this->useRILUK) {
            const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*this->K_);
            const TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<const TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());
            ConstTCrsMatrixPtr tpetraMat = xTpetraMat.getTpetra_CrsMatrix();

            auto RILUPreconditioner = rcp_dynamic_cast<Ifpack2::RILUK<TRowMatrix>>(Ifpack2Preconditioner_);
#ifdef HAVE_SHYLU_DDFROSCH_ZOLTAN2
            if (this->useZoltan2) {
                // if K is replaced, then we need to re-wrap into matrix filter
                auto rowMat = rcp_dynamic_cast<const TRowMatrix>(tpetraMat);
                auto filteredMat = rcp(new TRowMatrixFilterType(rowMat, this->perm, this->revperm));
                RILUPreconditioner->setMatrix(filteredMat);
            } else
#endif
            {
                RILUPreconditioner->setMatrix(tpetraMat);
            }
            return 0;
        }
        FROSCH_ASSERT(false,"FROSch::Ifpack2PreconditionerTpetra: updateMatrix() is not implemented for the Ifpack2PreconditionerTpetra yet.");
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

        auto solverName = this->ParameterList_->get("Solver","RILUK");
        this->useRILUK = (solverName == "RILUK");
        this->useZoltan2 = this->ParameterList_->get("RILUK: use reordering", false);

        Ifpack2::Details::OneLevelFactory<TRowMatrix> ifpack2Factory;
#ifdef HAVE_SHYLU_DDFROSCH_ZOLTAN2
        if (this->useZoltan2 && this->useRILUK) {
            // pre-ordering matrix, before constructing Ifpack2Preconditioner_ with the matrix
            {
                typedef Tpetra::RowGraph<LO, GO, NO> row_graph_type;
                typedef Zoltan2::TpetraRowGraphAdapter<row_graph_type> z2_adapter_type;
                typedef Zoltan2::OrderingProblem<z2_adapter_type> z2_ordering_problem_type;
                typedef Zoltan2::LocalOrderingSolution<LO> z2_ordering_solution_type;

                auto comm = tpetraMat->getRowMap()->getComm();
                auto constActiveGraph = Teuchos::rcp_const_cast<const row_graph_type>(tpetraMat->getGraph());
                z2_adapter_type Zoltan2Graph (constActiveGraph);

                Teuchos::ParameterList zoltan2_params = this->ParameterList_->sublist ("RILUK: reordering list");
                z2_ordering_problem_type MyOrderingProblem (&Zoltan2Graph, &zoltan2_params, comm);

                MyOrderingProblem.solve ();
                z2_ordering_solution_type sol (*MyOrderingProblem.getLocalOrderingSolution());
                this->perm = sol.getPermutationRCPConst (true);
                this->revperm = sol.getPermutationRCPConst ();
            }

            // wrap into matrix filter
            auto rowMat = rcp_dynamic_cast<const TRowMatrix>(tpetraMat);
            auto filteredMat = rcp(new TRowMatrixFilterType(rowMat, this->perm, this->revperm));
            Ifpack2Preconditioner_ = ifpack2Factory.create(solverName,filteredMat);
        } else
#endif
        {
            Ifpack2Preconditioner_ = ifpack2Factory.create(solverName,tpetraMat);
        }

        // set ifpack2 parameters
        ParameterListPtr ifpack2ParameterList = sublist(this->ParameterList_,"Ifpack2");
        if (ifpack2ParameterList->isSublist(solverName)) {
            ifpack2ParameterList = sublist(ifpack2ParameterList, solverName);
        }
        ifpack2ParameterList->setName("Ifpack2");
        Ifpack2Preconditioner_->setParameters(*ifpack2ParameterList);
    }

}

#endif
