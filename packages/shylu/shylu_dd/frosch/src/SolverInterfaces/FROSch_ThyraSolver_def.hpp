// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_THYRASOLVER_DEF_HPP
#define _FROSCH_THYRASOLVER_DEF_HPP

#include <FROSch_ThyraSolver_decl.hpp>


namespace FROSch {

    using namespace Stratimikos;
    using namespace Teuchos;
    using namespace Thyra;
    using namespace Xpetra;

    template<class SC,class LO,class GO,class NO>
    int ThyraSolver<SC,LO,GO,NO>::initialize()
    {
        FROSCH_TIMER_START_SOLVER(initializeTime,"ThyraSolver::initialize");
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int ThyraSolver<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_SOLVER(computeTime,"ThyraSolver::compute");
        FROSCH_ASSERT(this->IsInitialized_,"FROSch::ThyraSolver: !this->IsInitialized_");
        this->IsComputed_ = true;
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    void ThyraSolver<SC,LO,GO,NO>::apply(const XMultiVector &x,
                                         XMultiVector &y,
                                         ETransp mode,
                                         SC alpha,
                                         SC beta) const
    {
        FROSCH_TIMER_START_SOLVER(applyTime,"ThyraSolver::apply");
        FROSCH_ASSERT(this->IsComputed_,"FROSch::ThyraSolver: !this->IsComputed_.");

        ConstMultiVectorBasePtr xThyra = ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(rcpFromRef(x));
        if (YX_.is_null()) YX_ = XMultiVectorFactory::Build(y.getMap(),y.getNumVectors());
        if (YT_.is_null()) YT_ = rcp_const_cast<MultiVectorBase<SC> >(ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(YX_));

        EOpTransp tMode = NOTRANS;
        switch (mode) {
            case Teuchos::NO_TRANS:   tMode = Thyra::NOTRANS;       break;
            case Teuchos::TRANS:      tMode = Thyra::TRANS;         break;
            case Teuchos::CONJ_TRANS: tMode = Thyra::CONJTRANS;     break;
            default: FROSCH_ASSERT(false,"FROSch::ThyraPreconditionerTpetra: mode unknown.");
        }

        SolveStatus<SC> status = solve<SC>(*ThyraSolver_,tMode,*xThyra,YT_.ptr());

        // It seems that we have to convert the Thyra vector back to Xpetra. Is there a cheaper/more elegant way?
        // Same for ThyraPreconditioner
        XMultiVectorPtr yXpetra = ThyraUtils<SC,LO,GO,NO>::toXpetra(YT_,y.getMap()->getComm());
        y = *yXpetra;

        y.update(alpha,*YX_,beta);
    }

    template<class SC,class LO,class GO,class NO>
    int ThyraSolver<SC,LO,GO,NO>::updateMatrix(ConstXMatrixPtr k,
                                               bool reuseInitialize)
    {
        FROSCH_ASSERT(false,"FROSch::ThyraSolver: updateMatrix() is not implemented for the ThyraSolver yet.");
    }

    template<class SC,class LO,class GO,class NO>
    ThyraSolver<SC,LO,GO,NO>::ThyraSolver(ConstXMatrixPtr k,
                                          ParameterListPtr parameterList,
                                          string description) :
    Solver<SC,LO,GO,NO> (k,parameterList,description)
    {
        FROSCH_TIMER_START_SOLVER(ThyraSolverTime,"ThyraSolver::ThyraSolver");
        FROSCH_ASSERT(!this->K_.is_null(),"FROSch::ThyraSolver: K_ is null.");

        const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*this->K_);
        RCP<const Thyra::LinearOpBase<SC> > thyraOp = ThyraUtils<SC,LO,GO,NO>::toThyra(crsOp.getCrsMatrix());

        LinearSolverBuilder<SC> linearSolverBuilder;
        enableFROSch<SC,LO,GO,NO>(linearSolverBuilder);
        linearSolverBuilder.setParameterList(sublist(this->ParameterList_,"ThyraSolver"));

        LinearOpWithSolveFactoryBasePtr thyraFactory = linearSolverBuilder.createLinearSolveStrategy("");
        ThyraSolver_ = linearOpWithSolve(*thyraFactory,thyraOp);
    }

}

#endif
