// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_IPOUHARMONICCOARSEOPERATOR_DECL_HPP
#define _FROSCH_IPOUHARMONICCOARSEOPERATOR_DECL_HPP

#include <FROSch_ConstantPartitionOfUnity_def.hpp>
#include <FROSch_GDSWInterfacePartitionOfUnity_def.hpp>
#include <FROSch_GDSWStarInterfacePartitionOfUnity_def.hpp>
#include <FROSch_RGDSWInterfacePartitionOfUnity_def.hpp>

#include <FROSch_HarmonicCoarseOperator_def.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class  IPOUHarmonicCoarseOperator : public HarmonicCoarseOperator<SC,LO,GO,NO> {

    protected:

        using CommPtr                           = typename SchwarzOperator<SC,LO,GO,NO>::CommPtr;

        using XMapPtr                           = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr                      = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtr;
        using XMapPtrVecPtr                     = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtrVecPtr;
        using ConstXMapPtrVecPtr                = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtrVecPtr;
        using XMapPtrVecPtr2D                   = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtrVecPtr2D;
        using ConstXMapPtrVecPtr2D              = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtrVecPtr2D;

        using XMatrixPtr                        = typename SchwarzOperator<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr                   = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMatrixPtr;

        using XMultiVectorPtr                   = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVectorPtr;
        using ConstXMultiVectorPtr              = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMultiVectorPtr;
        using XMultiVectorPtrVecPtr             = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVectorPtrVecPtr;
        using ConstXMultiVectorPtrVecPtr        = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMultiVectorPtrVecPtr;

        using XCrsGraph                         = typename SchwarzOperator<SC,LO,GO,NO>::XCrsGraph;
        using GraphPtr                          = typename SchwarzOperator<SC,LO,GO,NO>::GraphPtr;
        using ConstXCrsGraphPtr                 = typename SchwarzOperator<SC,LO,GO,NO>::ConstXCrsGraphPtr;
        using ParameterListPtr                  = typename SchwarzOperator<SC,LO,GO,NO>::ParameterListPtr;

        using DDInterfacePtr                    = typename SchwarzOperator<SC,LO,GO,NO>::DDInterfacePtr;

        using EntitySetPtr                      = typename SchwarzOperator<SC,LO,GO,NO>::EntitySetPtr;

        using InterfaceEntityPtr                = typename SchwarzOperator<SC,LO,GO,NO>::InterfaceEntityPtr;

        using CoarseSpacePtr                    = typename SchwarzOperator<SC,LO,GO,NO>::CoarseSpacePtr;

        using PartitionOfUnityPtr               = typename SchwarzOperator<SC,LO,GO,NO>::PartitionOfUnityPtr;
        using InterfacePartitionOfUnityPtr      = typename SchwarzOperator<SC,LO,GO,NO>::InterfacePartitionOfUnityPtr;

        using LocalPartitionOfUnityBasisPtr     = typename SchwarzOperator<SC,LO,GO,NO>::LocalPartitionOfUnityBasisPtr;

        using SolverPtr                         = typename SchwarzOperator<SC,LO,GO,NO>::SolverPtr;
        using SolverFactoryPtr                  = typename SchwarzOperator<SC,LO,GO,NO>::SolverFactoryPtr;

        using UN                                = typename SchwarzOperator<SC,LO,GO,NO>::UN;
        using UNVecPtr                          = typename SchwarzOperator<SC,LO,GO,NO>::UNVecPtr;

        using LOVec                             = typename SchwarzOperator<SC,LO,GO,NO>::LOVec;
        using LOVecPtr                          = typename SchwarzOperator<SC,LO,GO,NO>::LOVecPtr;
        using LOVecPtr2D                        = typename SchwarzOperator<SC,LO,GO,NO>::LOVecPtr2D;

        using GOVec                             = typename SchwarzOperator<SC,LO,GO,NO>::GOVec;
        using GOVecPtr                          = typename SchwarzOperator<SC,LO,GO,NO>::GOVecPtr;
        using GOVecView                         = typename SchwarzOperator<SC,LO,GO,NO>::GOVecView;
        using GOVecPtr2D                        = typename SchwarzOperator<SC,LO,GO,NO>::GOVecPtr2D;

        using SCVec                             = typename SchwarzOperator<SC,LO,GO,NO>::SCVec;
        using SCVecPtr                          = typename SchwarzOperator<SC,LO,GO,NO>::SCVecPtr;
        using ConstSCVecPtr                     = typename SchwarzOperator<SC,LO,GO,NO>::ConstSCVecPtr;

        using BoolVecPtr                        = typename SchwarzOperator<SC,LO,GO,NO>::BoolVecPtr;

    public:

         IPOUHarmonicCoarseOperator(ConstXMatrixPtr k,
                                    ParameterListPtr parameterList);

        virtual int initialize()
        {
            FROSCH_ASSERT(false," IPOUHarmonicCoarseOperator cannot be built without a repeated Map");
        };

        int initialize(UN dimension,
                       UN dofsPerNode,
                       ConstXMapPtr nodesMap,
                       ConstXMapPtrVecPtr dofsMaps,
                       ConstXMultiVectorPtr nullSpaceBasis,
                       ConstXMultiVectorPtr nodeList,
                       GOVecPtr dirichletBoundaryDofs);

        int initialize(UN dimension,
                       UNVecPtr dofsPerNodeVec,
                       ConstXMapPtrVecPtr repeatedNodesMapVec,
                       ConstXMapPtrVecPtr2D repeatedDofMapsVec,
                       ConstXMultiVectorPtrVecPtr nullSpaceBasisVec,
                       ConstXMultiVectorPtrVecPtr nodeListVec,
                       GOVecPtr2D dirichletBoundaryDofsVec);

        void describe(FancyOStream &out,
                      const EVerbosityLevel verbLevel=Describable::verbLevel_default) const;

        string description() const;

        virtual XMapPtr BuildRepeatedMapCoarseLevel(ConstXMapPtr &nodesMap,
                                                    UN dofsPerNode,
                                                    ConstXMapPtrVecPtr dofsMaps,
                                                    UN partitionType);

    protected:

        int buildCoarseSpace(UN dimension,
                             UN dofsPerNode,
                             ConstXMapPtr nodesMap,
                             ConstXMapPtrVecPtr dofsMaps,
                             ConstXMultiVectorPtr nullSpaceBasis,
                             GOVecPtr dirichletBoundaryDofs,
                             ConstXMultiVectorPtr nodeList);


        int buildCoarseSpace(UN dimension,
                             UNVecPtr dofsPerNodeVec,
                             ConstXMapPtrVecPtr repeatedNodesMapVec,
                             ConstXMapPtrVecPtr2D repeatedDofMapsVec,
                             ConstXMultiVectorPtrVecPtr nullSpaceBasisVec,
                             GOVecPtr2D dirichletBoundaryDofsVec,
                             ConstXMultiVectorPtrVecPtr nodeListVec);

        virtual int resetCoarseSpaceBlock(UN blockId,
                                          UN dimension,
                                          UN dofsPerNode,
                                          ConstXMapPtr nodesMap,
                                          ConstXMapPtrVecPtr dofsMaps,
                                          ConstXMultiVectorPtr nullSpaceBasis,
                                          GOVecPtr dirichletBoundaryDofs,
                                          ConstXMultiVectorPtr nodeList);


        /*
         Todo: This should be vectors!
         vvvvvvvvvv
         */
        PartitionOfUnityPtr PartitionOfUnity_;

        LocalPartitionOfUnityBasisPtr LocalPartitionOfUnityBasis_;

        GOVec NumEnt_;
        /*
         ^^^^^^^^^^
         */
    };

}

#endif
