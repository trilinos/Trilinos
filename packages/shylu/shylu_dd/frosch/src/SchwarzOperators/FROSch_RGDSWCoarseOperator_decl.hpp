// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_RGDSWCOARSEOPERATOR_DECL_HPP
#define _FROSCH_RGDSWCOARSEOPERATOR_DECL_HPP

#include <FROSch_GDSWCoarseOperator_def.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class RGDSWCoarseOperator : public GDSWCoarseOperator<SC,LO,GO,NO> {

    protected:

        using CommPtr                 = typename SchwarzOperator<SC,LO,GO,NO>::CommPtr;

        using XMapPtr                 = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr            = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtr;
        using XMapPtrVecPtr           = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtrVecPtr;
        using ConstXMapPtrVecPtr      = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtrVecPtr;

        using XMatrixPtr              = typename SchwarzOperator<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr         = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMatrixPtr;

        using XMultiVectorPtr         = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVectorPtr;
        using ConstXMultiVectorPtr    = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMultiVectorPtr;
        using XMultiVectorPtrVecPtr   = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVectorPtrVecPtr;

        using ParameterListPtr        = typename SchwarzOperator<SC,LO,GO,NO>::ParameterListPtr;

        using DDInterfacePtr          = typename SchwarzOperator<SC,LO,GO,NO>::DDInterfacePtr;

        using EntitySetPtr            = typename SchwarzOperator<SC,LO,GO,NO>::EntitySetPtr;
        using EntitySetPtrVecPtr      = typename SchwarzOperator<SC,LO,GO,NO>::EntitySetPtrVecPtr;

        using InterfaceEntityPtr      = typename SchwarzOperator<SC,LO,GO,NO>::InterfaceEntityPtr;

        using UN                      = typename SchwarzOperator<SC,LO,GO,NO>::UN;

        using LOVec                   = typename SchwarzOperator<SC,LO,GO,NO>::LOVec;
        using LOVecPtr                = typename SchwarzOperator<SC,LO,GO,NO>::LOVecPtr;
        using LOVecPtr2D              = typename SchwarzOperator<SC,LO,GO,NO>::LOVecPtr2D;

        using GOVec                   = typename SchwarzOperator<SC,LO,GO,NO>::GOVec;
        using GOVecPtr                = typename SchwarzOperator<SC,LO,GO,NO>::GOVecPtr;
        using GOVecPtr2D              = typename SchwarzOperator<SC,LO,GO,NO>::GOVecPtr2D;

        using SCVecPtr                = typename SchwarzOperator<SC,LO,GO,NO>::SCVecPtr;

    public:

        RGDSWCoarseOperator(ConstXMatrixPtr k,
                            ParameterListPtr parameterList);

        virtual int resetCoarseSpaceBlock(UN blockId,
                                          UN dimension,
                                          UN dofsPerNode,
                                          ConstXMapPtr nodesMap,
                                          ConstXMapPtrVecPtr dofsMaps,
                                          GOVecPtr dirichletBoundaryDofs,
                                          ConstXMultiVectorPtr nodeList);

        virtual XMapPtr BuildRepeatedMapCoarseLevel(ConstXMapPtr &nodesMap,
                                                    UN dofsPerNode,
                                                    ConstXMapPtrVecPtr dofsMaps,
                                                   UN partitionType);


    protected:

        virtual XMultiVectorPtrVecPtr computeTranslations(UN blockId,
                                                          EntitySetPtr Roots,
                                                          EntitySetPtrVecPtr entitySetVector,
                                                          DistanceFunction distanceFunction = ConstantDistanceFunction);

        virtual XMultiVectorPtrVecPtr computeRotations(UN blockId,
                                                       UN dimension,
                                                       ConstXMultiVectorPtr nodeList,
                                                       EntitySetPtr Roots,
                                                       EntitySetPtrVecPtr entitySetVector,
                                                       DistanceFunction distanceFunction = ConstantDistanceFunction);
    };

}

#endif
