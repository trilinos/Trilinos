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
        return 0;
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

        DefaultLinearSolverBuilder linearSolverBuilder;
        enableFROSch<LO,GO,NO>(linearSolverBuilder);
        linearSolverBuilder.setParameterList(sublist(this->ParameterList_,"ThyraPreconditioner"));

        PreconditionerFactoryBasePtr thyraFactory = linearSolverBuilder.createPreconditioningStrategy("");
        ThyraPreconditioner_ = prec(*thyraFactory,thyraOp);
    }

}

#endif
