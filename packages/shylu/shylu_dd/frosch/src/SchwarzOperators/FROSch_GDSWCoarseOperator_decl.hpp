// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_GDSWCOARSEOPERATOR_DECL_HPP
#define _FROSCH_GDSWCOARSEOPERATOR_DECL_HPP

#include <FROSch_HarmonicCoarseOperator_def.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class GDSWCoarseOperator : public HarmonicCoarseOperator<SC,LO,GO,NO> {

    protected:

        using CommPtr                       = typename SchwarzOperator<SC,LO,GO,NO>::CommPtr;

        using XMapPtr                       = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr                  = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtr;
        using XMapPtrVecPtr                 = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtrVecPtr;
        using ConstXMapPtrVecPtr            = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtrVecPtr;
        using XMapPtrVecPtr2D               = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtrVecPtr2D;
        using ConstXMapPtrVecPtr2D          = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtrVecPtr2D;

        using XMatrixPtr                    = typename SchwarzOperator<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr               = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMatrixPtr;

        using XMultiVectorPtr               = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVectorPtr;
        using ConstXMultiVectorPtr          = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMultiVectorPtr;
        using XMultiVectorPtrVecPtr         = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVectorPtrVecPtr;
        using ConstXMultiVectorPtrVecPtr    = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMultiVectorPtrVecPtr;

        using ParameterListPtr              = typename SchwarzOperator<SC,LO,GO,NO>::ParameterListPtr;

        using DDInterfacePtr                = typename SchwarzOperator<SC,LO,GO,NO>::DDInterfacePtr;

        using EntitySetPtr                  = typename SchwarzOperator<SC,LO,GO,NO>::EntitySetPtr;

        using SolverPtr                     = typename SchwarzOperator<SC,LO,GO,NO>::SolverPtr;
        using SolverFactoryPtr              = typename SchwarzOperator<SC,LO,GO,NO>::SolverFactoryPtr;

        using UN                            = typename SchwarzOperator<SC,LO,GO,NO>::UN;
        using UNVecPtr                      = typename SchwarzOperator<SC,LO,GO,NO>::UNVecPtr;

        using LOVec                         = typename SchwarzOperator<SC,LO,GO,NO>::LOVec;
        using LOVecPtr                      = typename SchwarzOperator<SC,LO,GO,NO>::LOVecPtr;
        using LOVecPtr2D                    = typename SchwarzOperator<SC,LO,GO,NO>::LOVecPtr2D;

        using GOVec                         = typename SchwarzOperator<SC,LO,GO,NO>::GOVec;
        using GOVecPtr                      = typename SchwarzOperator<SC,LO,GO,NO>::GOVecPtr;
        using GOVecView                     = typename SchwarzOperator<SC,LO,GO,NO>::GOVecView;
        using GOVecPtr2D                    = typename SchwarzOperator<SC,LO,GO,NO>::GOVecPtr2D;

        using SCVec                         = typename SchwarzOperator<SC,LO,GO,NO>::SCVec;
        using SCVecPtr                      = typename SchwarzOperator<SC,LO,GO,NO>::SCVecPtr;

        using BoolVecPtr                    = typename SchwarzOperator<SC,LO,GO,NO>::BoolVecPtr;

    public:

        GDSWCoarseOperator(ConstXMatrixPtr k,
                           ParameterListPtr parameterList);

        virtual int initialize()
        {
            FROSCH_ASSERT(false,"GDSWCoarseOperator cannot be built without a repeated Map");
        };

        int initialize(UN dimension,
                       ConstXMapPtr repeatedMap);

        int initialize(UN dimension,
                       ConstXMapPtr repeatedMap,
                       GOVecPtr dirichletBoundaryDofs);

        int initialize(UN dimension,
                       UN dofsPerNode,
                       ConstXMapPtr repeatedNodesMap,
                       ConstXMapPtrVecPtr RepeatedDofMaps);

        int initialize(UN dimension,
                       UN dofsPerNode,
                       ConstXMapPtr repeatedNodesMap,
                       ConstXMapPtrVecPtr RepeatedDofMaps,
                       GOVecPtr dirichletBoundaryDofs);

        int initialize(UN dimension,
                       UN dofsPerNode,
                       ConstXMapPtr repeatedNodesMap,
                       ConstXMapPtrVecPtr RepeatedDofMaps,
                       ConstXMultiVectorPtr nodeList);

        int initialize(UN dimension,
                       UN dofsPerNode,
                       ConstXMapPtr repeatedNodesMap,
                       ConstXMapPtrVecPtr RepeatedDofMaps,
                       GOVecPtr dirichletBoundaryDofs,
                       ConstXMultiVectorPtr nodeList);

        int initialize(UN dimension,
                       UNVecPtr dofsPerNodeVec,
                       ConstXMapPtrVecPtr repeatedNodesMapVec,
                       ConstXMapPtrVecPtr2D repeatedDofMapsVec,
                       GOVecPtr2D dirichletBoundaryDofsVec,
                       ConstXMultiVectorPtrVecPtr nodeListVec);

        void describe(FancyOStream &out,
                      const EVerbosityLevel verbLevel=Describable::verbLevel_default) const;

        string description() const;

        // AH: Could this be moved to protected?
        virtual XMapPtr BuildRepeatedMapCoarseLevel(ConstXMapPtr &nodesMap,
                                                    UN dofsPerNode,
                                                    ConstXMapPtrVecPtr dofsMaps,
                                                    UN partitionType );

    protected:

        int buildCoarseSpace(UN dimension,
                             ConstXMapPtr nodesMap);

        int buildCoarseSpace(UN dimension,
                             ConstXMapPtr nodesMap,
                             GOVecPtr dirichletBoundaryDofs); // Das kann man auch mit in den Fall davor reinnehmen ?!

        int buildCoarseSpace(UN dimension,
                             UN dofsPerNode,
                             ConstXMapPtr nodesMap,
                             ConstXMapPtrVecPtr dofsMaps);

        int buildCoarseSpace(UN dimension,
                             UN dofsPerNode,
                             ConstXMapPtr nodesMap,
                             ConstXMapPtrVecPtr dofsMaps,
                             GOVecPtr dirichletBoundaryDofs);

        int buildCoarseSpace(UN dimension,
                             UN dofsPerNode,
                             ConstXMapPtr nodesMap,
                             ConstXMapPtrVecPtr dofsMaps,
                             ConstXMultiVectorPtr nodeList);

        int buildCoarseSpace(UN dimension,
                             UN dofsPerNode,
                             ConstXMapPtr nodesMap,
                             ConstXMapPtrVecPtr dofsMaps,
                             GOVecPtr dirichletBoundaryDofs,
                             ConstXMultiVectorPtr nodeList);

        int buildCoarseSpace(UN dimension,
                             UNVecPtr dofsPerNodeVec,
                             ConstXMapPtrVecPtr repeatedNodesMapVec,
                             ConstXMapPtrVecPtr2D repeatedDofMapsVec,
                             GOVecPtr2D dirichletBoundaryDofsVec,
                             ConstXMultiVectorPtrVecPtr nodeListVec);

        virtual int resetCoarseSpaceBlock(UN blockId,
                                          UN dimension,
                                          UN dofsPerNode,
                                          ConstXMapPtr nodesMap,
                                          ConstXMapPtrVecPtr dofsMaps,
                                          GOVecPtr dirichletBoundaryDofs,
                                          ConstXMultiVectorPtr nodeList);

        DDInterfacePtr DDInterface_;
    };

}

#endif
