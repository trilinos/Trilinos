// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_CONSTANTPARTITIONOFUNITY_DECL_HPP
#define _FROSCH_CONSTANTPARTITIONOFUNITY_DECL_HPP

#include <FROSch_PartitionOfUnity_def.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class ConstantPartitionOfUnity : public PartitionOfUnity<SC,LO,GO,NO> {

    protected:

        using CommPtr                       = typename PartitionOfUnity<SC,LO,GO,NO>::CommPtr;

        using XMapPtr                       = typename PartitionOfUnity<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr                  = typename PartitionOfUnity<SC,LO,GO,NO>::ConstXMapPtr;
        using XMapPtrVecPtr                 = typename PartitionOfUnity<SC,LO,GO,NO>::XMapPtrVecPtr;
        using ConstXMapPtrVecPtr            = typename PartitionOfUnity<SC,LO,GO,NO>::ConstXMapPtrVecPtr;

        using XMultiVectorPtr               = typename PartitionOfUnity<SC,LO,GO,NO>::XMultiVectorPtr;
        using ConstXMultiVectorPtr          = typename PartitionOfUnity<SC,LO,GO,NO>::ConstXMultiVectorPtr;
        using XMultiVectorPtrVecPtr         = typename PartitionOfUnity<SC,LO,GO,NO>::XMultiVectorPtrVecPtr;
        using ConstXMultiVectorPtrVecPtr    = typename PartitionOfUnity<SC,LO,GO,NO>::ConstXMultiVectorPtrVecPtr;

        using ParameterListPtr              = typename PartitionOfUnity<SC,LO,GO,NO>::ParameterListPtr;

        using DDInterfacePtr                = typename PartitionOfUnity<SC,LO,GO,NO>::DDInterfacePtr;

        using EntitySetPtr                  = typename PartitionOfUnity<SC,LO,GO,NO>::EntitySetPtr;

        using UN                            = typename PartitionOfUnity<SC,LO,GO,NO>::UN;

        using LOVec                         = typename PartitionOfUnity<SC,LO,GO,NO>::LOVec;

        using GOVec                         = typename PartitionOfUnity<SC,LO,GO,NO>::GOVec;
        using GOVecView                     = typename PartitionOfUnity<SC,LO,GO,NO>::GOVecView;

        using SCVec                         = typename PartitionOfUnity<SC,LO,GO,NO>::SCVec;

    public:

        ConstantPartitionOfUnity(CommPtr mpiComm,
                                 CommPtr serialComm,
                                 UN dimension,
                                 UN dofsPerNode,
                                 ConstXMapPtr nodesMap,
                                 ConstXMapPtrVecPtr dofsMaps,
                                 ParameterListPtr parameterList,
                                 Verbosity verbosity = All,
                                 UN levelID = 1,
                                 DDInterfacePtr ddInterface = null);

        virtual ~ConstantPartitionOfUnity();

        virtual int removeDirichletNodes(GOVecView dirichletBoundaryDofs,
                                         ConstXMultiVectorPtr nodeList = null);

        virtual int computePartitionOfUnity(ConstXMultiVectorPtr nodeList = null);

    protected:

        DDInterfacePtr DDInterface_;

        bool UseVolumes_ = false;

        EntitySetPtr Volumes_;
    };

}

#endif
