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
// Questions? Contact Christian Hochmuth (c.hochmuth@uni-koeln.de)
//
// ************************************************************************
//@HEADER

#ifndef _FROSCH_MULTIPLICATIVEOPERATOR_DEF_HPP
#define _FROSCH_MULTIPLICATIVEOPERATOR_DEF_HPP

#include <FROSch_MultiplicativeOperator_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    MultiplicativeOperator<SC,LO,GO,NO>::MultiplicativeOperator(ConstXMatrixPtr k,
                                                                ParameterListPtr parameterList) :
    SchwarzOperator<SC,LO,GO,NO> (k, parameterList)
    {
        FROSCH_DETAILTIMER_START_LEVELID(multiplicativeOperatorTime,"MultiplicativeOperator::MultiplicativeOperator");
    }

    template <class SC,class LO,class GO,class NO>
    MultiplicativeOperator<SC,LO,GO,NO>::MultiplicativeOperator(ConstXMatrixPtr k,
                                                                SchwarzOperatorPtrVecPtr operators,
                                                                ParameterListPtr parameterList) :
    SchwarzOperator<SC,LO,GO,NO> (k, parameterList)
    {
        FROSCH_DETAILTIMER_START_LEVELID(multiplicativeOperatorTime,"MultiplicativeOperator::MultiplicativeOperator");
        OperatorVector_.push_back(operators.at(0));
        for (unsigned i=1; i<operators.size(); i++) {
            FROSCH_ASSERT(operators[i]->OperatorDomainMap().SameAs(OperatorVector_[i]->OperatorDomainMap()),"The DomainMaps of the operators are not identical.");
            FROSCH_ASSERT(operators[i]->OperatorRangeMap().SameAs(OperatorVector_[i]->OperatorRangeMap()),"The RangeMaps of the operators are not identical.");

            OperatorVector_.push_back(operators[i]);
            EnableOperators_.push_back(true);
        }
    }

    template <class SC,class LO,class GO,class NO>
    MultiplicativeOperator<SC,LO,GO,NO>::~MultiplicativeOperator()
    {

    }

    template <class SC,class LO,class GO,class NO>
    void MultiplicativeOperator<SC,LO,GO,NO>::preApplyCoarse(XMultiVector &x,
                                                             XMultiVector &y)
    {
        FROSCH_DETAILTIMER_START_LEVELID(preApplyCoarseTime,"MultiplicativeOperator::preApplyCoarse");
        FROSCH_ASSERT(this->OperatorVector_.size()==2,"Should be a Two-Level Operator.");
        this->OperatorVector_[1]->apply(x,y,true);
    }

    // Y = alpha * A^mode * X + beta * Y
    template <class SC,class LO,class GO,class NO>
    void MultiplicativeOperator<SC,LO,GO,NO>::apply(const XMultiVector &x,
                                                    XMultiVector &y,
                                                    bool usePreconditionerOnly,
                                                    ETransp mode,
                                                    SC alpha,
                                                    SC beta) const
    {
        FROSCH_TIMER_START_LEVELID(applyTime,"MultiplicativeOperator::apply");
        FROSCH_ASSERT(usePreconditionerOnly,"MultiplicativeOperator can only be used as a preconditioner.");
        FROSCH_ASSERT(this->OperatorVector_.size()==2,"Should be a Two-Level Operator.");


        if (XTmp_.is_null()) XTmp_ = MultiVectorFactory<SC,LO,GO,NO>::Build(x.getMap(),x.getNumVectors());
        *XTmp_ = x; // Need this for the case when x aliases y

        if (YTmp_.is_null()) XMultiVectorPtr YTmp_ = MultiVectorFactory<SC,LO,GO,NO>::Build(y.getMap(),y.getNumVectors());
        *YTmp_ = y; // for the second apply

        this->OperatorVector_[0]->apply(*XTmp_,*YTmp_,true);

        this->K_->apply(*YTmp_,*XTmp_);

        this->OperatorVector_[1]->apply(*XTmp_,*XTmp_,true);

        YTmp_->update(ScalarTraits<SC>::one(),*XTmp_,-ScalarTraits<SC>::one());
        y.update(alpha,*YTmp_,beta);
    }

    template <class SC,class LO,class GO,class NO>
    int MultiplicativeOperator<SC,LO,GO,NO>::initialize()
    {
        if (this->Verbose_) {
            FROSCH_ASSERT(false,"ERROR: Each of the Operators has to be initialized manually.");
        }
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int MultiplicativeOperator<SC,LO,GO,NO>::initialize(ConstXMapPtr repeatedMap)
    {
        if (this->Verbose_) {
            FROSCH_ASSERT(false,"ERROR: Each of the Operators has to be initialized manually.");
        }
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int MultiplicativeOperator<SC,LO,GO,NO>::compute()
    {
        if (this->Verbose_) {
            FROSCH_ASSERT(false,"ERROR: Each of the Operators has to be computed manually.");
        }
        return 0;
    }


    template <class SC,class LO,class GO,class NO>
    typename MultiplicativeOperator<SC,LO,GO,NO>::ConstXMapPtr MultiplicativeOperator<SC,LO,GO,NO>::getDomainMap() const
    {
        return OperatorVector_[0]->getDomainMap();
    }

    template <class SC,class LO,class GO,class NO>
    typename MultiplicativeOperator<SC,LO,GO,NO>::ConstXMapPtr MultiplicativeOperator<SC,LO,GO,NO>::getRangeMap() const
    {
        return OperatorVector_[0]->getRangeMap();
    }

    template <class SC,class LO,class GO,class NO>
    void MultiplicativeOperator<SC,LO,GO,NO>::describe(FancyOStream &out,
                                                       const EVerbosityLevel verbLevel) const
    {
        FROSCH_ASSERT(false,"describe() has to be implemented properly...");
    }

    template <class SC,class LO,class GO,class NO>
    string MultiplicativeOperator<SC,LO,GO,NO>::description() const
    {
        string labelString = "Level operator: ";

        for (UN i=0; i<OperatorVector_.size(); i++) {
            labelString += OperatorVector_.at(i)->description();
            if (i<OperatorVector_.size()-1) {
                labelString += ",";
            }
        }
        return labelString;
    }

    template <class SC,class LO,class GO,class NO>
    int MultiplicativeOperator<SC,LO,GO,NO>::addOperator(SchwarzOperatorPtr op)
    {
        FROSCH_DETAILTIMER_START_LEVELID(addOperatorTime,"MultiplicativeOperator::addOperator");
        int ret = 0;
        if (OperatorVector_.size()>0) {
            if (!op->getDomainMap()->isSameAs(*OperatorVector_[0]->getDomainMap())) {
                if (this->Verbose_) cerr << "MultiplicativeOperator<SC,LO,GO,NO>::addOperator(SchwarzOperatorPtr op)\t\t!op->getDomainMap().isSameAs(OperatorVector_[0]->getDomainMap())\n";
                ret -= 1;
            }
            if (!op->getRangeMap()->isSameAs(*OperatorVector_[0]->getRangeMap())) {
                if (this->Verbose_) cerr << "MultiplicativeOperator<SC,LO,GO,NO>::addOperator(SchwarzOperatorPtr op)\t\t!op->getRangeMap().isSameAs(OperatorVector_[0]->getRangeMap())\n";
                ret -= 10;
            }
            //FROSCH_ASSERT(op->OperatorDomainMap().SameAs(OperatorVector_.at(0)->OperatorDomainMap()),"The DomainMaps of the operators are not identical.");
            //FROSCH_ASSERT(op->OperatorRangeMap().SameAs(OperatorVector_.at(0)->OperatorRangeMap()),"The RangeMaps of the operators are not identical.");
        }
        OperatorVector_.push_back(op);
        EnableOperators_.push_back(true);
        return ret;
    }

    template <class SC,class LO,class GO,class NO>
    int MultiplicativeOperator<SC,LO,GO,NO>::addOperators(SchwarzOperatorPtrVecPtr operators)
    {
        FROSCH_DETAILTIMER_START_LEVELID(addOperatorsTime,"MultiplicativeOperator::addOperators");
        int ret = 0;
        for (UN i=1; i<operators.size(); i++) {
            if (0>addOperator(operators[i])) ret -= pow(10,i);
        }
        return ret;
    }

    template <class SC,class LO,class GO,class NO>
    int MultiplicativeOperator<SC,LO,GO,NO>::resetOperator(UN iD,
                                                           SchwarzOperatorPtr op)
    {
        FROSCH_DETAILTIMER_START_LEVELID(resetOperatorTime,"MultiplicativeOperator::resetOperator");
        FROSCH_ASSERT(iD<OperatorVector_.size(),"iD exceeds the length of the OperatorVector_");
        int ret = 0;
        if (!op->getDomainMap().isSameAs(OperatorVector_[0]->getDomainMap())) {
            if (this->Verbose_) cerr << "MultiplicativeOperator<SC,LO,GO,NO>::addOperator(SchwarzOperatorPtr op)\t\t!op->getDomainMap().isSameAs(OperatorVector_[0]->getDomainMap())\n";
            ret -= 1;
        }
        if (!op->getRangeMap().isSameAs(OperatorVector_[0]->getRangeMap())) {
            if (this->Verbose_) cerr << "MultiplicativeOperator<SC,LO,GO,NO>::addOperator(SchwarzOperatorPtr op)\t\t!op->getRangeMap().isSameAs(OperatorVector_[0]->getRangeMap())\n";
            ret -= 10;
        }
        OperatorVector_[iD] = op;
        return ret;
    }

    template <class SC,class LO,class GO,class NO>
    int MultiplicativeOperator<SC,LO,GO,NO>::enableOperator(UN iD,
                                                            bool enable)
    {
        FROSCH_DETAILTIMER_START_LEVELID(enableOperatorTime,"MultiplicativeOperator::enableOperator");
        EnableOperators_[iD] = enable;
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    typename MultiplicativeOperator<SC,LO,GO,NO>::UN MultiplicativeOperator<SC,LO,GO,NO>::getNumOperators()
    {
        return OperatorVector_.size();
    }

}

#endif
