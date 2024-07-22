// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_RGDSWINTERFACEPARTITIONOFUNITY_DECL_HPP
#define _FROSCH_RGDSWINTERFACEPARTITIONOFUNITY_DECL_HPP

#include <FROSch_GDSWInterfacePartitionOfUnity_def.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class RGDSWInterfacePartitionOfUnity : public GDSWInterfacePartitionOfUnity<SC,LO,GO,NO> {

    protected:

        using CommPtr                       = typename PartitionOfUnity<SC,LO,GO,NO>::CommPtr;

        using XMap                          = typename PartitionOfUnity<SC,LO,GO,NO>::XMap;
        using XMapPtr                       = typename PartitionOfUnity<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr                  = typename PartitionOfUnity<SC,LO,GO,NO>::ConstXMapPtr;
        using XMapPtrVecPtr                 = typename PartitionOfUnity<SC,LO,GO,NO>::XMapPtrVecPtr;
        using ConstXMapPtrVecPtr            = typename PartitionOfUnity<SC,LO,GO,NO>::ConstXMapPtrVecPtr;

        using XMatrix                       = typename PartitionOfUnity<SC,LO,GO,NO>::XMatrix;
        using XMatrixPtr                    = typename PartitionOfUnity<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr               = typename PartitionOfUnity<SC,LO,GO,NO>::ConstXMatrixPtr;

        using XMultiVector                  = typename PartitionOfUnity<SC,LO,GO,NO>::XMultiVector;
        using ConstXMultiVectorPtr          = typename PartitionOfUnity<SC,LO,GO,NO>::ConstXMultiVectorPtr;
        using XMultiVectorPtr               = typename PartitionOfUnity<SC,LO,GO,NO>::XMultiVectorPtr;
        using XMultiVectorPtrVecPtr         = typename PartitionOfUnity<SC,LO,GO,NO>::XMultiVectorPtrVecPtr;
        using ConstXMultiVectorPtrVecPtr    = typename PartitionOfUnity<SC,LO,GO,NO>::ConstXMultiVectorPtrVecPtr;

        using ParameterListPtr              = typename PartitionOfUnity<SC,LO,GO,NO>::ParameterListPtr;

        using DDInterfacePtr                = typename PartitionOfUnity<SC,LO,GO,NO>::DDInterfacePtr;

        using EntitySetPtr                  = typename PartitionOfUnity<SC,LO,GO,NO>::EntitySetPtr;
        using EntitySetPtrVecPtr            = typename PartitionOfUnity<SC,LO,GO,NO>::EntitySetPtrVecPtr;

        using InterfaceEntityPtr            = typename PartitionOfUnity<SC,LO,GO,NO>::InterfaceEntityPtr;

        using UN                            = typename PartitionOfUnity<SC,LO,GO,NO>::UN;

        using GOVec                         = typename PartitionOfUnity<SC,LO,GO,NO>::GOVec;
        using GOVecView                     = typename PartitionOfUnity<SC,LO,GO,NO>::GOVecView;

    public:

        RGDSWInterfacePartitionOfUnity(CommPtr mpiComm,
                                       CommPtr serialComm,
                                       UN dimension,
                                       UN dofsPerNode,
                                       ConstXMapPtr nodesMap,
                                       ConstXMapPtrVecPtr dofsMaps,
                                       ParameterListPtr parameterList,
                                       Verbosity verbosity = All,
                                       UN levelID = 1);

        virtual int computePartitionOfUnity(ConstXMultiVectorPtr nodeList);

    protected:

        bool UseRoots_ = false;

        EntitySetPtr Roots_;

        EntitySetPtrVecPtr EntitySetVector_ = EntitySetPtrVecPtr(0);

        DistanceFunction DistanceFunction_ = ConstantDistanceFunction;
    };

}

#endif
