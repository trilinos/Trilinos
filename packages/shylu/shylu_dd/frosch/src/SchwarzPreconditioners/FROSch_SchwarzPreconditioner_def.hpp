// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_SCHWARZPRECONDITIONER_DEF_HPP
#define _FROSCH_SCHWARZPRECONDITIONER_DEF_HPP

#include <FROSch_SchwarzPreconditioner_decl.hpp>

#include <FROSch_SumOperator_def.hpp>
#include <FROSch_MultiplicativeOperator_def.hpp>
#include <FROSch_AlgebraicOverlappingOperator_def.hpp>
#include <FROSch_GDSWCoarseOperator_def.hpp>
#include <FROSch_RGDSWCoarseOperator_def.hpp>
#include <FROSch_IPOUHarmonicCoarseOperator_def.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    SchwarzPreconditioner<SC,LO,GO,NO>::SchwarzPreconditioner(ParameterListPtr parameterList,
                                                              CommPtr comm) :
    MpiComm_ (comm),
    ParameterList_ (parameterList),
    Verbose_ (comm->getRank()==0),
    LevelID_ (ParameterList_->get("Level ID",UN(1)))
    {

    }

    template <class SC,class LO,class GO,class NO>
    SchwarzPreconditioner<SC,LO,GO,NO>::~SchwarzPreconditioner()
    {

    }

    template <class SC,class LO,class GO,class NO>
    bool SchwarzPreconditioner<SC,LO,GO,NO>::isInitialized() const
    {
        return IsInitialized_; // TODO: Das hat noch keine Bedeutung
    }

    template <class SC,class LO,class GO,class NO>
    bool SchwarzPreconditioner<SC,LO,GO,NO>::isComputed() const
    {
        return IsComputed_; // TODO: Das hat noch keine Bedeutung
    }

    template <class SC,class LO,class GO,class NO>
    void SchwarzPreconditioner<SC,LO,GO,NO>::residual(const XMultiVector & X,
                                                      const XMultiVector & B,
                                                      XMultiVector& R) const
    {
        SC one = Teuchos::ScalarTraits<SC>::one(), negone = -one;
        apply(X,R);
        R.update(one,B,negone);
    }

}

#endif
