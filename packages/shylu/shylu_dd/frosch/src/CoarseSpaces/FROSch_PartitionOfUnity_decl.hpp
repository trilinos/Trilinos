// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_PARTITIONOFUNITY_DECL_HPP
#define _FROSCH_PARTITIONOFUNITY_DECL_HPP

#include <FROSch_DDInterface_def.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class PartitionOfUnity {

    protected:

        using CommPtr                       = RCP<const Comm<int> >;

        using XMap                          = Map<LO,GO,NO>;
        using XMapPtr                       = RCP<XMap>;
        using ConstXMapPtr                  = RCP<const XMap>;
        using XMapPtrVec                    = Array<XMapPtr>;
        using XMapPtrVecPtr                 = ArrayRCP<XMapPtr>;
        using ConstXMapPtrVecPtr            = ArrayRCP<ConstXMapPtr>;

        using XMatrix                       = Matrix<SC,LO,GO,NO>;
        using XMatrixPtr                    = RCP<XMatrix>;
        using ConstXMatrixPtr               = RCP<const XMatrix>;

        using XMultiVector                  = MultiVector<SC,LO,GO,NO>;
        using ConstXMultiVectorPtr          = RCP<const XMultiVector>;
        using XMultiVectorPtr               = RCP<XMultiVector>;
        using XMultiVectorPtrVecPtr         = ArrayRCP<XMultiVectorPtr>;
        using ConstXMultiVectorPtrVecPtr    = ArrayRCP<ConstXMultiVectorPtr>;

        using ParameterListPtr              = RCP<ParameterList>;

        using DDInterfacePtr                = RCP<DDInterface<SC,LO,GO,NO> >;
        using ConstDDInterfacePtr           = RCP<const DDInterface<SC,LO,GO,NO> >;

        using EntitySetPtr                  = RCP<EntitySet<SC,LO,GO,NO> >;
        using EntitySetPtrVecPtr            = ArrayRCP<EntitySetPtr>;

        using InterfaceEntityPtr            = RCP<InterfaceEntity<SC,LO,GO,NO> >;

        using UN                            = unsigned;
        using ConstUN                       = const UN;

        using LOVec                         = Array<LO>;
        using LOVecPtr                      = ArrayRCP<LO>;
        using LOVecPtr2D                    = ArrayRCP<LOVecPtr>;

        using GOVec                         = Array<GO>;
        using GOVecView                     = ArrayView<GO>;

        using SCVec                         = Array<SC>;
        using SCVecPtr                      = ArrayRCP<SC>;

    public:

        PartitionOfUnity(CommPtr mpiComm,
                         CommPtr serialComm,
                         UN dofsPerNode,
                         ConstXMapPtr nodesMap,
                         ConstXMapPtrVecPtr dofsMaps,
                         ParameterListPtr parameterList,
                         Verbosity verbosity = All,
                         UN levelID = 1);

        virtual ~PartitionOfUnity();

        virtual int removeDirichletNodes(GOVecView dirichletBoundaryDofs,
                                         ConstXMultiVectorPtr nodeList = null) = 0;

        virtual int computePartitionOfUnity(ConstXMultiVectorPtr nodeList = null) = 0;

        int assembledPartitionOfUnityMaps();

        ConstXMultiVectorPtrVecPtr getLocalPartitionOfUnity() const;

        ConstXMapPtrVecPtr getPartitionOfUnityMaps() const;

        ConstXMapPtr getAssembledPartitionOfUnityMap() const;

    protected:

        CommPtr MpiComm_;
        CommPtr SerialComm_;

        ParameterListPtr ParameterList_;

        ConstXMultiVectorPtrVecPtr LocalPartitionOfUnity_;

        ConstXMapPtrVecPtr PartitionOfUnityMaps_;

        ConstXMapPtr AssmbledPartitionOfUnityMap_;

        bool Verbose_ = false;

        Verbosity Verbosity_ = All;

        const UN LevelID_;
    };

}

#endif
