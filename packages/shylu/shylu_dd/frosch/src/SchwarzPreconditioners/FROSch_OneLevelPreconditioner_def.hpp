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
    
    template <class SC,class LO,class GO,class NO>
    OneLevelPreconditioner<SC,LO,GO,NO>::OneLevelPreconditioner(CrsMatrixPtr k,
                                                                            ParameterListPtr parameterList) :
    SchwarzPreconditioner<SC,LO,GO,NO> (parameterList,k->getRangeMap()->getComm()),
    K_ (k),
    SumOperator_ (new SumOperator<SC,LO,GO,NO>(k->getRangeMap()->getComm())),
    OverlappingOperator_ ()
    {
        if (!this->ParameterList_->get("OverlappingOperator Type","AlgebraicOverlappingOperator").compare("AlgebraicOverlappingOperator")) {
            OverlappingOperator_ = AlgebraicOverlappingOperatorPtr(new AlgebraicOverlappingOperator<SC,LO,GO,NO>(k,sublist(parameterList,"AlgebraicOverlappingOperator")));
        } else {
            FROSCH_ASSERT(0!=0,"OverlappingOperator Type unkown.");
        }
        SumOperator_->addOperator(OverlappingOperator_);
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
        MapPtr repeatedMap;
        if (buildRepeatedMap) {
            repeatedMap = BuildRepeatedMap(this->K_);
        }
        return initialize(overlap,repeatedMap);
    }
    
    template <class SC,class LO,class GO,class NO>
    int OneLevelPreconditioner<SC,LO,GO,NO>::initialize(int overlap,
                                                        MapPtr repeatedMap)
    {
        int ret = 0;
        if (overlap<0) {
            overlap = this->ParameterList_->get("Overlap",1);
        }
        if (!this->ParameterList_->get("OverlappingOperator Type","AlgebraicOverlappingOperator").compare("AlgebraicOverlappingOperator")) {
            AlgebraicOverlappingOperatorPtr algebraicOverlappigOperator = Teuchos::rcp_static_cast<AlgebraicOverlappingOperator<SC,LO,GO,NO> >(OverlappingOperator_);
            ret = algebraicOverlappigOperator->initialize(overlap,repeatedMap);
        } else {
            FROSCH_ASSERT(0!=0,"OverlappingOperator Type unkown.");
        }
        return ret;
    }
    
    template <class SC,class LO,class GO,class NO>
    int OneLevelPreconditioner<SC,LO,GO,NO>::compute()
    {
        return OverlappingOperator_->compute();
    }
    
    template <class SC,class LO,class GO,class NO>
    void OneLevelPreconditioner<SC,LO,GO,NO>::apply(const MultiVector &x,
                                                    MultiVector &y,
                                                    Teuchos::ETransp mode,
                                                    SC alpha,
                                                    SC beta) const
    {
        return SumOperator_->apply(x,y,true,mode,alpha,beta);
    }
    
    template <class SC,class LO,class GO,class NO>
    typename OneLevelPreconditioner<SC,LO,GO,NO>::ConstMapPtr OneLevelPreconditioner<SC,LO,GO,NO>::getDomainMap() const
    {
        return K_->getDomainMap();
    }
    
    template <class SC,class LO,class GO,class NO>
    typename OneLevelPreconditioner<SC,LO,GO,NO>::ConstMapPtr OneLevelPreconditioner<SC,LO,GO,NO>::getRangeMap() const
    {
        return K_->getRangeMap();
    }
    
    template <class SC,class LO,class GO,class NO>
    void OneLevelPreconditioner<SC,LO,GO,NO>::describe(Teuchos::FancyOStream &out,
                                                        const Teuchos::EVerbosityLevel verbLevel) const
    {
        SumOperator_->describe(out,verbLevel);
    }
    
    template <class SC,class LO,class GO,class NO>
    std::string OneLevelPreconditioner<SC,LO,GO,NO>::description() const
    {
        return "One-Level Preconditioner";
    }
    
}

#endif
