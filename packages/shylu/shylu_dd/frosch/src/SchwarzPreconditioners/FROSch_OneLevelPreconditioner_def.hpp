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

#ifndef _FROSCH_ONELEVELPRECONDITIONER_DEF_HPP
#define _FROSCH_ONELEVELPRECONDITIONER_DEF_HPP

#include <FROSch_OneLevelPreconditioner_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    OneLevelPreconditioner<SC,LO,GO,NO>::OneLevelPreconditioner(ConstXMatrixPtr k,
                                                                ParameterListPtr parameterList) :
    SchwarzPreconditioner<SC,LO,GO,NO> (parameterList,k->getRangeMap()->getComm()),
    K_ (k),
    SumOperator_ (new SumOperator<SC,LO,GO,NO>(k->getRangeMap()->getComm())),
    MultiplicativeOperator_ (new MultiplicativeOperator<SC,LO,GO,NO>(k,parameterList)),
    OverlappingOperator_ ()
    {
        FROSCH_DETAILTIMER_START_LEVELID(oneLevelPreconditionerTime,"OneLevelPreconditioner::OneLevelPreconditioner");
        if (!this->ParameterList_->get("OverlappingOperator Type","AlgebraicOverlappingOperator").compare("AlgebraicOverlappingOperator")) {
            // Set the LevelID in the sublist
            parameterList->sublist("AlgebraicOverlappingOperator").set("Level ID",this->LevelID_);
            OverlappingOperator_ = AlgebraicOverlappingOperatorPtr(new AlgebraicOverlappingOperator<SC,LO,GO,NO>(k,sublist(parameterList,"AlgebraicOverlappingOperator")));
        } else {
            FROSCH_ASSERT(false,"OverlappingOperator Type unkown.");
        }
        if (!this->ParameterList_->get("Level Combination","Additive").compare("Multiplicative")) {
            UseMultiplicative_ = true;
        }
        if (UseMultiplicative_) {
            MultiplicativeOperator_->addOperator(OverlappingOperator_);
        } else {
            SumOperator_->addOperator(OverlappingOperator_);
        }

    }

    template <class SC,class LO,class GO,class NO>
    int OneLevelPreconditioner<SC,LO,GO,NO>::initialize(bool useDefaultParameters)
    {
        if (useDefaultParameters) {
            return initialize(1,false);
        } else {
            return initialize(this->ParameterList_->get("Overlap",1),false);
        }
    }

    template <class SC,class LO,class GO,class NO>
    int OneLevelPreconditioner<SC,LO,GO,NO>::initialize(int overlap,
                                                        bool buildRepeatedMap)
    {
        ConstXMapPtr repeatedMap;
        if (buildRepeatedMap) {
            repeatedMap = BuildRepeatedMap(this->K_->getCrsGraph());
        }
        return initialize(overlap,repeatedMap);
    }

    template <class SC,class LO,class GO,class NO>
    int OneLevelPreconditioner<SC,LO,GO,NO>::initialize(int overlap,
                                                        ConstXMapPtr repeatedMap)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"OneLevelPreconditioner::initialize");
        int ret = 0;
        if (overlap<0) {
            overlap = this->ParameterList_->get("Overlap",1);
        }
        if (!this->ParameterList_->get("OverlappingOperator Type","AlgebraicOverlappingOperator").compare("AlgebraicOverlappingOperator")) {
            AlgebraicOverlappingOperatorPtr algebraicOverlappigOperator = rcp_static_cast<AlgebraicOverlappingOperator<SC,LO,GO,NO> >(OverlappingOperator_);
            ret = algebraicOverlappigOperator->initialize(overlap,repeatedMap);
        } else {
            FROSCH_ASSERT(false,"OverlappingOperator Type unkown.");
        }
        return ret;
    }

    template <class SC,class LO,class GO,class NO>
    int OneLevelPreconditioner<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_LEVELID(computeTime,"OneLevelPreconditioner::compute");
        return OverlappingOperator_->compute();
    }

    template <class SC,class LO,class GO,class NO>
    void OneLevelPreconditioner<SC,LO,GO,NO>::apply(const XMultiVector &x,
                                                    XMultiVector &y,
                                                    ETransp mode,
                                                    SC alpha,
                                                    SC beta) const
    {
        FROSCH_TIMER_START_LEVELID(applyTime,"OneLevelPreconditioner::apply");
        if (UseMultiplicative_) {
            return MultiplicativeOperator_->apply(x,y,true,mode,alpha,beta);
        }
        else{
            return SumOperator_->apply(x,y,true,mode,alpha,beta);
        }
    }

    template <class SC,class LO,class GO,class NO>
    typename OneLevelPreconditioner<SC,LO,GO,NO>::ConstXMapPtr OneLevelPreconditioner<SC,LO,GO,NO>::getDomainMap() const
    {
        return K_->getDomainMap();
    }

    template <class SC,class LO,class GO,class NO>
    typename OneLevelPreconditioner<SC,LO,GO,NO>::ConstXMapPtr OneLevelPreconditioner<SC,LO,GO,NO>::getRangeMap() const
    {
        return K_->getRangeMap();
    }

    template <class SC,class LO,class GO,class NO>
    void OneLevelPreconditioner<SC,LO,GO,NO>::describe(FancyOStream &out,
                                                       const EVerbosityLevel verbLevel) const
    {
        if (UseMultiplicative_) {
            MultiplicativeOperator_->describe(out,verbLevel);
        }
        else{
            SumOperator_->describe(out,verbLevel);
        }

    }

    template <class SC,class LO,class GO,class NO>
    string OneLevelPreconditioner<SC,LO,GO,NO>::description() const
    {
        return "One-Level Preconditioner";
    }

    template <class SC,class LO,class GO,class NO>
    int OneLevelPreconditioner<SC,LO,GO,NO>::resetMatrix(ConstXMatrixPtr &k)
    {
        FROSCH_TIMER_START_LEVELID(resetMatrixTime,"OneLevelPreconditioner::resetMatrix");
        K_ = k;
        OverlappingOperator_->resetMatrix(K_);
        return 0;
    }

}

#endif
