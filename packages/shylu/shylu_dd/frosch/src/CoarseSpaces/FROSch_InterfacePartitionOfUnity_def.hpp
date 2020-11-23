//@HEADER
// ************************************************************************
//
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Alexander Heinlein (alexander.heinlein@uni-koeln.de)
//
// ************************************************************************
//@HEADER

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
            FROSCH_ASSERT(false,"FROSch::InterfacePartitionOfUnity : ERROR: Specify a valid communication strategy for the identification of the interface components.");
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
