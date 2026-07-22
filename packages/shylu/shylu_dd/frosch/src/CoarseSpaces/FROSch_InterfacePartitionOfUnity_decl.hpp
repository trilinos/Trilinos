// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_INTERFACEPARTITIONOFUNITY_DECL_HPP
#define _FROSCH_INTERFACEPARTITIONOFUNITY_DECL_HPP

#include <FROSch_PartitionOfUnity_def.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class InterfacePartitionOfUnity : public PartitionOfUnity<SC,LO,GO,NO> {

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
        using ConstDDInterfacePtr           = typename PartitionOfUnity<SC,LO,GO,NO>::ConstDDInterfacePtr;

        using EntitySetPtr                  = typename PartitionOfUnity<SC,LO,GO,NO>::EntitySetPtr;
        using EntitySetPtrVecPtr            = typename PartitionOfUnity<SC,LO,GO,NO>::EntitySetPtrVecPtr;

        using InterfaceEntityPtr            = typename PartitionOfUnity<SC,LO,GO,NO>::InterfaceEntityPtr;

        using UN                            = typename PartitionOfUnity<SC,LO,GO,NO>::UN;
        using ConstUN                       = typename PartitionOfUnity<SC,LO,GO,NO>::ConstUN;

        using GOVec                         = typename PartitionOfUnity<SC,LO,GO,NO>::GOVec;
        using GOVecView                     = typename PartitionOfUnity<SC,LO,GO,NO>::GOVecView;

        using SCVecPtr                      = typename PartitionOfUnity<SC,LO,GO,NO>::SCVecPtr;

    public:

        InterfacePartitionOfUnity(CommPtr mpiComm,
                                  CommPtr serialComm,
                                  UN dimension,
                                  UN dofsPerNode,
                                  ConstXMapPtr nodesMap,
                                  ConstXMapPtrVecPtr dofsMaps,
                                  ParameterListPtr parameterList,
                                  Verbosity verbosity = All,
                                  UN levelID = 1);

        virtual ~InterfacePartitionOfUnity();        

        virtual int sortInterface(ConstXMatrixPtr matrix = null,
                                  ConstXMultiVectorPtr nodeList = null) = 0;

        ConstDDInterfacePtr getDDInterface() const;
        
        DDInterfacePtr getDDInterfaceNonConst() const;
        
        virtual int computePartitionOfUnity(ConstXMultiVectorPtr nodeList = null) = 0;

    protected:

        DDInterfacePtr DDInterface_;
    };

}

#endif
