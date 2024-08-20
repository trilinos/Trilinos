// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_SOLVER_DEF_HPP
#define _FROSCH_SOLVER_DEF_HPP

#include <FROSch_Solver_decl.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template<class SC,class LO,class GO,class NO>
    const typename Solver<SC,LO,GO,NO>::ConstXMapPtr Solver<SC,LO,GO,NO>::getDomainMap() const
    {
        return K_->getDomainMap();
    }

    template<class SC,class LO,class GO,class NO>
    const typename Solver<SC,LO,GO,NO>::ConstXMapPtr Solver<SC,LO,GO,NO>::getRangeMap() const
    {
        return K_->getRangeMap();
    }

    template<class SC,class LO,class GO,class NO>
    bool Solver<SC,LO,GO,NO>::isInitialized() const
    {
        return IsInitialized_;
    }

    template<class SC,class LO,class GO,class NO>
    bool Solver<SC,LO,GO,NO>::isComputed() const
    {
        return IsComputed_;
    }

    template<class SC,class LO,class GO,class NO>
    void Solver<SC,LO,GO,NO>::residual(const XMultiVector & x,
                                       const XMultiVector & b,
                                       XMultiVector& r) const
    {
        apply(x,r);
        r.update(ScalarTraits<SC>::one(),b,-ScalarTraits<SC>::one());
    }

    template<class SC,class LO,class GO,class NO>
    Solver<SC,LO,GO,NO>::Solver(ConstXMatrixPtr k,
                                ParameterListPtr parameterList,
                                string description) :
    K_ (k),
    ParameterList_ (parameterList),
    Description_ (description)
    {}

}

#endif
