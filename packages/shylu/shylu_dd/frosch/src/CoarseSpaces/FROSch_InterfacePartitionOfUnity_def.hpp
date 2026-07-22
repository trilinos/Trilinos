// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_INTERFACEPARTITIONOFUNITY_DEF_HPP
#define _FROSCH_INTERFACEPARTITIONOFUNITY_DEF_HPP

#include <FROSch_InterfacePartitionOfUnity_decl.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    InterfacePartitionOfUnity<SC,LO,GO,NO>::InterfacePartitionOfUnity(CommPtr mpiComm,
                                                                      CommPtr serialComm,
                                                                      UN dimension,
                                                                      UN dofsPerNode,
                                                                      ConstXMapPtr nodesMap,
                                                                      ConstXMapPtrVecPtr dofsMaps,
                                                                      ParameterListPtr parameterList,
                                                                      Verbosity verbosity,
                                                                      UN levelID) :
    PartitionOfUnity<SC,LO,GO,NO> (mpiComm,serialComm,dofsPerNode,nodesMap,dofsMaps,parameterList,verbosity,levelID)
    {
        FROSCH_DETAILTIMER_START_LEVELID(interfacePartitionOfUnityTime,"InterfacePartitionOfUnity::InterfacePartitionOfUnity");
        CommunicationStrategy communicationStrategy = CreateOneToOneMap;
        if (!this->ParameterList_->get("Interface Communication Strategy","CreateOneToOneMap").compare("CrsMatrix")) {
            communicationStrategy = CommCrsMatrix;
        } else if (!this->ParameterList_->get("Interface Communication Strategy","CreateOneToOneMap").compare("CrsGraph")) {
            communicationStrategy = CommCrsGraph;
        } else if (!this->ParameterList_->get("Interface Communication Strategy","CreateOneToOneMap").compare("CreateOneToOneMap")) {
            communicationStrategy = CreateOneToOneMap;
        } else {
            FROSCH_ASSERT(false,"FROSch::InterfacePartitionOfUnity: Specify a valid communication strategy for the identification of the interface components.");
        }

        DDInterface_.reset(new DDInterface<SC,LO,GO,NO>(dimension,dofsPerNode,nodesMap.getConst(),this->Verbosity_,this->LevelID_,communicationStrategy));
        DDInterface_->resetGlobalDofs(dofsMaps);
    }

    template <class SC,class LO,class GO,class NO>
    InterfacePartitionOfUnity<SC,LO,GO,NO>::~InterfacePartitionOfUnity()
    {

    }

    template <class SC,class LO,class GO,class NO>
    typename InterfacePartitionOfUnity<SC,LO,GO,NO>::ConstDDInterfacePtr InterfacePartitionOfUnity<SC,LO,GO,NO>::getDDInterface() const
    {
        return DDInterface_.getConst();
    }

    template <class SC,class LO,class GO,class NO>
    typename InterfacePartitionOfUnity<SC,LO,GO,NO>::DDInterfacePtr InterfacePartitionOfUnity<SC,LO,GO,NO>::getDDInterfaceNonConst() const
    {
        return DDInterface_;
    }
}

#endif
