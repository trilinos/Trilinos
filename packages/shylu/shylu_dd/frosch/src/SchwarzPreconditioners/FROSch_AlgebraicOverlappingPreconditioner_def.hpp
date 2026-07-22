// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_ALGEBRAICOVERLAPPINGPRECONDITIONER_DEF_HPP
#define _FROSCH_ALGEBRAICOVERLAPPINGPRECONDITIONER_DEF_HPP

#include <FROSch_AlgebraicOverlappingPreconditioner_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>::AlgebraicOverlappingPreconditioner(ConstXMatrixPtr k,
                                                                                        ParameterListPtr parameterList) :
    SchwarzPreconditioner<SC,LO,GO,NO> (parameterList,k->getRangeMap()->getComm()),
    K_ (k),
    SumOperator_ (new SumOperator<SC,LO,GO,NO>(k->getRangeMap()->getComm()))
    {
        FROSCH_DETAILTIMER_START_LEVELID(algebraicOverlappingPreconditionerTime,"AlgebraicOverlappingPreconditioner::AlgebraicOverlappingPreconditioner");
        // Set the LevelID in the sublist
        parameterList->sublist("AlgebraicOverlappingOperator").set("Level ID",this->LevelID_);
        OverlappingOperator_.reset(new AlgebraicOverlappingOperator<SC,LO,GO,NO>(k,sublist(parameterList,"AlgebraicOverlappingOperator")));
        SumOperator_->addOperator(OverlappingOperator_);
    }

    template <class SC,class LO,class GO,class NO>
    int AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>::initialize(bool useDefaultParameters)
    {
        return initialize(1,null);
    }

    template <class SC,class LO,class GO,class NO>
    int AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>::initialize(int overlap,
                                                                    ConstXMapPtr repeatedMap)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"AlgebraicOverlappingPreconditioner::initialize");
        return OverlappingOperator_->initialize(overlap,repeatedMap);
    }

    template <class SC,class LO,class GO,class NO>
    int AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_LEVELID(computeTime,"AlgebraicOverlappingPreconditioner::compute");
        return OverlappingOperator_->compute();
    }

    template <class SC,class LO,class GO,class NO>
    void AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>::apply(const XMultiVector &x,
                                                                XMultiVector &y,
                                                                ETransp mode,
                                                                SC alpha,
                                                                SC beta) const
    {
        FROSCH_TIMER_START_LEVELID(applyTime,"AlgebraicOverlappingPreconditioner::apply");
        return SumOperator_->apply(x,y,true,mode,alpha,beta);
    }

    template <class SC,class LO,class GO,class NO>
    const typename AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>::ConstXMapPtr AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>::getDomainMap() const
    {
        return K_->getDomainMap();
    }

    template <class SC,class LO,class GO,class NO>
    const typename AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>::ConstXMapPtr AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>::getRangeMap() const
    {
        return K_->getRangeMap();
    }

    template <class SC,class LO,class GO,class NO>
    void AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>::describe(FancyOStream &out,
                                                                   const EVerbosityLevel verbLevel) const
    {
        SumOperator_->describe(out,verbLevel);
    }

    template <class SC,class LO,class GO,class NO>
    string AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>::description() const
    {
        return "Algebraic Overlapping Preconditioner";
    }

    template <class SC,class LO,class GO,class NO>
    int AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>::resetMatrix(ConstXMatrixPtr &k)
    {
        FROSCH_DETAILTIMER_START_LEVELID(resetMatrixTime,"TwoLevelPreconditioner::resetMatrix");
        this->K_ = k;
        OverlappingOperator_->resetMatrix(this->K_);
        return 0;
    }
}

#endif
