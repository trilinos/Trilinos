// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_TWOLEVELPRECONDITIONER_DECL_HPP
#define _FROSCH_TWOLEVELPRECONDITIONER_DECL_HPP

#include <FROSch_OneLevelPreconditioner_def.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class TwoLevelPreconditioner : public OneLevelPreconditioner<SC,LO,GO,NO> {

    protected:

        using XMapPtr                           = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr                      = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMapPtr;
        using XMapPtrVecPtr                     = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMapPtrVecPtr;
        using ConstXMapPtrVecPtr                = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMapPtrVecPtr;

        using XMatrixPtr                        = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr                   = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMatrixPtr;

        using XMultiVectorPtr                   = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMultiVectorPtr;
        using ConstXMultiVectorPtr              = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMultiVectorPtr;

        using ParameterListPtr                  = typename SchwarzPreconditioner<SC,LO,GO,NO>::ParameterListPtr;

        using AlgebraicOverlappingOperatorPtr   = typename SchwarzPreconditioner<SC,LO,GO,NO>::AlgebraicOverlappingOperatorPtr;
        using CoarseOperatorPtr                 = typename SchwarzPreconditioner<SC,LO,GO,NO>::CoarseOperatorPtr;
        using GDSWCoarseOperatorPtr             = typename SchwarzPreconditioner<SC,LO,GO,NO>::GDSWCoarseOperatorPtr;
        using RGDSWCoarseOperatorPtr            = typename SchwarzPreconditioner<SC,LO,GO,NO>::RGDSWCoarseOperatorPtr;
        using IPOUHarmonicCoarseOperatorPtr     = typename SchwarzPreconditioner<SC,LO,GO,NO>::IPOUHarmonicCoarseOperatorPtr;

        using UN                                = typename SchwarzPreconditioner<SC,LO,GO,NO>::UN;

        using GOVecPtr                          = typename SchwarzPreconditioner<SC,LO,GO,NO>::GOVecPtr;

    public:
        using OneLevelPreconditioner<SC,LO,GO,NO>::initialize;

        TwoLevelPreconditioner(ConstXMatrixPtr k,
                               ParameterListPtr parameterList);

        int initialize(bool useDefaultParameters = true);

        int initialize(UN dimension,
                       int overlap,
                       UN dofsPerNode,
                       DofOrdering dofOrdering);

        int initialize(UN dimension,
                       int overlap,
                       ConstXMapPtr repeatedMap,
                       UN dofsPerNode,
                       DofOrdering dofOrdering,
                       ConstXMultiVectorPtr nodeList = null);

        int initialize(UN dimension,
                       UN dofsPerNode,
                       int overlap = -1,
                       ConstXMultiVectorPtr nullSpaceBasis = null,
                       ConstXMultiVectorPtr nodeList = null,
                       DofOrdering dofOrdering = NodeWise,
                       ConstXMapPtr repeatedMap = null,
                       ConstXMapPtrVecPtr dofsMaps = null,
                       GOVecPtr dirichletBoundaryDofs = null);

        int compute();

        void describe(FancyOStream &out,
                      const EVerbosityLevel verbLevel=Describable::verbLevel_default) const;

        string description() const;

        int resetMatrix(ConstXMatrixPtr &k);

    protected:

        CoarseOperatorPtr CoarseOperator_;

    };

}

#endif
