// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_SCHWARZOPERATOR_DEF_HPP
#define _FROSCH_SCHWARZOPERATOR_DEF_HPP

#include <FROSch_SchwarzOperator_decl.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template<class SC,class LO,class GO,class NO>
    SchwarzOperator<SC,LO,GO,NO>::SchwarzOperator(CommPtr comm) :
    MpiComm_ (comm),
    Verbose_ (comm->getRank()==0)
    {

    }

    template<class SC,class LO,class GO,class NO>
    SchwarzOperator<SC,LO,GO,NO>::SchwarzOperator(ConstXMatrixPtr k,
                                                  ParameterListPtr parameterList) :
    MpiComm_ (k->getRangeMap()->getComm()),
    K_ (k),
    ParameterList_ (parameterList),
    Verbose_ (MpiComm_->getRank()==0),
    LevelID_ (ParameterList_->get("Level ID",UN(1)))
    {
        FROSCH_ASSERT(getDomainMap()->isSameAs(*getRangeMap()),"SchwarzOperator assumes DomainMap==RangeMap");
    }

    template<class SC,class LO,class GO,class NO>
    SchwarzOperator<SC,LO,GO,NO>::~SchwarzOperator()
    {

    }

    template<class SC,class LO,class GO,class NO>
    void SchwarzOperator<SC,LO,GO,NO>::apply(const XMultiVector &x,
                                             XMultiVector &y,
                                             ETransp mode,
                                             SC alpha,
                                             SC beta) const
    {
        return apply(x,y,false,mode,alpha,beta);
    }

    template<class SC,class LO,class GO,class NO>
    const typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtr SchwarzOperator<SC,LO,GO,NO>::getDomainMap() const
    {
        return K_->getDomainMap();
    }

    template<class SC,class LO,class GO,class NO>
    const typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtr SchwarzOperator<SC,LO,GO,NO>::getRangeMap() const
    {
        return K_->getRangeMap();
    }

    template<class SC,class LO,class GO,class NO>
    bool SchwarzOperator<SC,LO,GO,NO>::isInitialized() const
    {
        return IsInitialized_;
    }

    template<class SC,class LO,class GO,class NO>
    bool SchwarzOperator<SC,LO,GO,NO>::isComputed() const
    {
        return IsComputed_;
    }

    template<class SC,class LO,class GO,class NO>
    int SchwarzOperator<SC,LO,GO,NO>::resetMatrix(ConstXMatrixPtr &k) {
    // Maybe set IsComputed_ = false ? -> Go through code to be saver/cleaner
    // This function must be actively called by the user, and is only for recycling purposes.
    // The preconditioner is still computed and this point or the preconditioner was never computed and can now be computed with the matrix k now.
        K_ = k;
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    void SchwarzOperator<SC,LO,GO,NO>::residual(const XMultiVector & X,
                                                const XMultiVector & B,
                                                XMultiVector& R) const
    {
        SC one = ScalarTraits<SC>::one(), negone = -one;
        apply(X,R);
        R.update(one,B,negone);
    }
}

#endif
