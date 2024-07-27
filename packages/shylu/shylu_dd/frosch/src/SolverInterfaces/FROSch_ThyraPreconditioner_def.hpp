// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_THYRAPRECONDITIONER_DEF_HPP
#define _FROSCH_THYRAPRECONDITIONER_DEF_HPP

#include <FROSch_ThyraPreconditioner_decl.hpp>


namespace FROSch {

    using namespace Stratimikos;
    using namespace Teuchos;
    using namespace Thyra;
    using namespace Xpetra;

    template<class SC,class LO,class GO,class NO>
    int ThyraPreconditioner<SC,LO,GO,NO>::initialize()
    {
        FROSCH_TIMER_START_SOLVER(initializeTime,"ThyraPreconditioner::initialize");
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int ThyraPreconditioner<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_SOLVER(computeTime,"ThyraPreconditioner::compute");
        FROSCH_ASSERT(this->IsInitialized_,"FROSch::ThyraPreconditioner: !this->IsInitialized_");
        this->IsComputed_ = true;
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    void ThyraPreconditioner<SC,LO,GO,NO>::apply(const XMultiVector &x,
                                                 XMultiVector &y,
                                                 ETransp mode,
                                                 SC alpha,
                                                 SC beta) const
    {
        FROSCH_TIMER_START_SOLVER(applyTime,"ThyraPreconditioner::apply");
        FROSCH_ASSERT(this->IsComputed_,"FROSch::ThyraPreconditioner: !this->IsComputed_.");

        ConstMultiVectorBasePtr xThyra = ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(rcpFromRef(x));
        MultiVectorBasePtr yThyra = rcp_const_cast<MultiVectorBase<SC> >(ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(rcpFromRef(y)));

        EOpTransp tMode = NOTRANS;
        switch (mode) {
            case Teuchos::NO_TRANS:   tMode = Thyra::NOTRANS;       break;
            case Teuchos::TRANS:      tMode = Thyra::TRANS;         break;
            case Teuchos::CONJ_TRANS: tMode = Thyra::CONJTRANS;     break;
            default: FROSCH_ASSERT(false,"FROSch::ThyraPreconditioner: mode unknown.");
        }

        ThyraPreconditioner_->getUnspecifiedPrecOp()->apply(tMode,*xThyra,yThyra.ptr(),alpha,beta);

        // It seems that we have to convert the Thyra vector back to Xpetra. Is there a cheaper/more elegant way?
        // Same for ThyraSolver
        XMultiVectorPtr yXpetra = ThyraUtils<SC,LO,GO,NO>::toXpetra(yThyra,y.getMap()->getComm());
        y = *yXpetra;
    }

    template<class SC,class LO,class GO,class NO>
    int ThyraPreconditioner<SC,LO,GO,NO>::updateMatrix(ConstXMatrixPtr k,
                                               bool reuseInitialize)
    {
        FROSCH_ASSERT(false,"FROSch::ThyraPreconditioner: updateMatrix() is not implemented for the ThyraPreconditioner yet.");
    }

    template<class SC,class LO,class GO,class NO>
    ThyraPreconditioner<SC,LO,GO,NO>::ThyraPreconditioner(ConstXMatrixPtr k,
                                                          ParameterListPtr parameterList,
                                                          string description) :
    Solver<SC,LO,GO,NO> (k,parameterList,description)
    {
        FROSCH_TIMER_START_SOLVER(ThyraPreconditionerTime,"ThyraPreconditioner::ThyraPreconditioner");
        FROSCH_ASSERT(!this->K_.is_null(),"FROSch::ThyraPreconditioner: K_ is null.");

        const CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<const CrsMatrixWrap<SC,LO,GO,NO>&>(*this->K_);
        RCP<const Thyra::LinearOpBase<SC> > thyraOp = ThyraUtils<SC,LO,GO,NO>::toThyra(crsOp.getCrsMatrix());

        LinearSolverBuilder<SC> linearSolverBuilder;
        enableFROSch<SC,LO,GO,NO>(linearSolverBuilder);
        linearSolverBuilder.setParameterList(sublist(this->ParameterList_,"ThyraPreconditioner"));

        PreconditionerFactoryBasePtr thyraFactory = linearSolverBuilder.createPreconditioningStrategy("");
        ThyraPreconditioner_ = prec(*thyraFactory,thyraOp);
    }

}

#endif
