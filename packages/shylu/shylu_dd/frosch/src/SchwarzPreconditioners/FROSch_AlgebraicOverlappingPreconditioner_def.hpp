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

#ifndef _FROSCH_ALGEBRAICOVERLAPPINGPRECONDITIONER_DEF_HPP
#define _FROSCH_ALGEBRAICOVERLAPPINGPRECONDITIONER_DEF_HPP

#include <FROSch_AlgebraicOverlappingPreconditioner_decl.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>::AlgebraicOverlappingPreconditioner(ConstXMatrixPtr k,
                                                                                        ParameterListPtr parameterList) :
    SchwarzPreconditioner<SC,LO,GO,NO> (parameterList,k->getRangeMap()->getComm()),
    K_ (k),
    SumOperator_ (new SumOperator<SC,LO,GO,NO>(k->getRangeMap()->getComm())),
    FirstLevelOperator_ ()
    {
        FROSCH_TIMER_START_LEVELID(algebraicOverlappingPreconditionerTime,"AlgebraicOverlappingPreconditioner::AlgebraicOverlappingPreconditioner");
        // Set the LevelID in the sublist
        parameterList->sublist("GDSWCoarseOperator").set("Level ID",this->LevelID_);
        FirstLevelOperator_.reset(new AlgebraicOverlappingOperator<SC,LO,GO,NO>(k,sublist(parameterList,"AlgebraicOverlappingOperator")));
        SumOperator_->addOperator(FirstLevelOperator_);
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
        return FirstLevelOperator_->initialize(overlap,repeatedMap);
    }

    template <class SC,class LO,class GO,class NO>
    int AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_LEVELID(computeTime,"AlgebraicOverlappingPreconditioner::compute");
        return FirstLevelOperator_->compute();
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
    typename AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>::ConstXMapPtr AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>::getDomainMap() const
    {
        return K_->getDomainMap();
    }

    template <class SC,class LO,class GO,class NO>
    typename AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>::ConstXMapPtr AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>::getRangeMap() const
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
    std::string AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>::description() const
    {
        return "Algebraic Overlapping Preconditioner";
    }

}

#endif
