// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_RGDSWPRECONDITIONER_DECL_HPP
#define _FROSCH_RGDSWPRECONDITIONER_DECL_HPP

#include <FROSch_OneLevelPreconditioner_def.hpp>
#include <FROSch_AlgebraicOverlappingPreconditioner_def.hpp>
#include <FROSch_RGDSWCoarseOperator_def.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class RGDSWPreconditioner : public AlgebraicOverlappingPreconditioner<SC,LO,GO,NO> {

    protected:

        using XMapPtr                   = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr              = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMapPtr;
        using XMapPtrVecPtr             = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMapPtrVecPtr;
        using ConstXMapPtrVecPtr        = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMapPtrVecPtr;

        using XMatrixPtr                = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr           = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMatrixPtr;

        using XMultiVectorPtr           = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMultiVectorPtr;
        using ConstXMultiVectorPtr      = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMultiVectorPtr;

        using ParameterListPtr          = typename SchwarzPreconditioner<SC,LO,GO,NO>::ParameterListPtr;

        using RGDSWCoarseOperatorPtr    = typename SchwarzPreconditioner<SC,LO,GO,NO>::RGDSWCoarseOperatorPtr;

        using UN                        = typename SchwarzPreconditioner<SC,LO,GO,NO>::UN;

        using GOVecPtr                  = typename SchwarzPreconditioner<SC,LO,GO,NO>::GOVecPtr;

    public:
        using AlgebraicOverlappingPreconditioner<SC,LO,GO,NO>::initialize;

        RGDSWPreconditioner(ConstXMatrixPtr k,
                            ParameterListPtr parameterList);

        int initialize(bool useDefaultParameters = true);

        int initialize(ConstXMapPtr repeatedMap,
                       bool useDefaultParameters = true);

        int initialize(GOVecPtr &dirichletBoundaryDofs,
                       bool useDefaultParameters = true);

        int initialize(ConstXMapPtr repeatedMap,
                       GOVecPtr &dirichletBoundaryDofs,
                       bool useDefaultParameters = true);

        int initialize(UN dimension,
                       int overlap);

        int initialize(UN dimension,
                       int overlap,
                       ConstXMapPtr repeatedMap);

        int initialize(UN dimension,
                       int overlap,
                       ConstXMapPtr repeatedMap,
                       GOVecPtr &dirichletBoundaryDofs);

        int initialize(UN dimension,
                       UN dofsPerNode,
                       DofOrdering dofOrdering,
                       int overlap,
                       ConstXMapPtr repeatedMap);

        int initialize(UN dimension,
                       UN dofsPerNode,
                       DofOrdering dofOrdering,
                       int overlap,
                       ConstXMapPtr repeatedMap,
                       GOVecPtr &dirichletBoundaryDofs);

        int initialize(UN dimension,
                       UN dofsPerNode,
                       DofOrdering dofOrdering,
                       int overlap,
                       ConstXMapPtr repeatedMap,
                       ConstXMultiVectorPtr &nodeList);

        int initialize(UN dimension,
                       UN dofsPerNode,
                       DofOrdering dofOrdering,
                       int overlap,
                       ConstXMapPtr repeatedMap,
                       GOVecPtr &dirichletBoundaryDofs,
                       ConstXMultiVectorPtr &nodeList);

        int compute();

        void describe(FancyOStream &out,
                      const EVerbosityLevel verbLevel=Describable::verbLevel_default) const;

        string description() const; // @suppress("Type cannot be resolved")

        virtual int resetMatrix(ConstXMatrixPtr &k);

    protected:

        RGDSWCoarseOperatorPtr CoarseOperator_;
    };

}

#endif
