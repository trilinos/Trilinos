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

#ifndef _FROSCH_CONSTANTPARTITIONOFUNITY_DEF_HPP
#define _FROSCH_CONSTANTPARTITIONOFUNITY_DEF_HPP

#include <FROSch_ConstantPartitionOfUnity_decl.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    ConstantPartitionOfUnity<SC,LO,GO,NO>::ConstantPartitionOfUnity(CommPtr mpiComm,
                                                                    CommPtr serialComm,
                                                                    UN dimension,
                                                                    UN dofsPerNode,
                                                                    ConstXMapPtr nodesMap,
                                                                    ConstXMapPtrVecPtr dofsMaps,
                                                                    ParameterListPtr parameterList,
                                                                    Verbosity verbosity,
                                                                    UN levelID,
                                                                    DDInterfacePtr ddInterface) :
    PartitionOfUnity<SC,LO,GO,NO> (mpiComm,serialComm,dofsPerNode,nodesMap,dofsMaps,parameterList,verbosity,levelID),
    DDInterface_ (ddInterface)
    {
        FROSCH_TIMER_START_LEVELID(constantPartitionOfUnityTime,"ConstantPartitionOfUnity::ConstantPartitionOfUnity");
        
        if (!this->ParameterList_->get("Type","Full").compare("Full")) {
            UseVolumes_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("Volumes")) {
            UseVolumes_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("Custom")) {
            UseVolumes_ = this->ParameterList_->sublist("Custom").get("Volumes",false);
        } else {
            FROSCH_ASSERT(false,"FROSch::ConstantPartitionOfUnity : ERROR: Specify a valid Type.");
        }
        
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
        
        if (DDInterface_.is_null()) DDInterface_.reset(new DDInterface<SC,LO,GO,NO>(dimension,dofsPerNode,nodesMap.getConst(),this->Verbosity_,this->LevelID_,communicationStrategy));
        FROSCH_ASSERT(DDInterface_->getInterface()->getEntity(0)->getNumNodes()==0,"FROSch::ConstantPartitionOfUnity : ERROR: Is only reasonable if there is no interface.");
        DDInterface_->resetGlobalDofs(dofsMaps);
        Volumes_ = DDInterface_->getInterior()->deepCopy();
        
        this->LocalPartitionOfUnity_ = XMultiVectorPtrVecPtr(1);
        this->PartitionOfUnityMaps_ = XMapPtrVecPtr(1);
    }

    template <class SC,class LO,class GO,class NO>
    ConstantPartitionOfUnity<SC,LO,GO,NO>::~ConstantPartitionOfUnity()
    {

    }

    template <class SC,class LO,class GO,class NO>
    int ConstantPartitionOfUnity<SC,LO,GO,NO>::removeDirichletNodes(GOVecView dirichletBoundaryDofs,
                                                                    ConstXMultiVectorPtr nodeList)
    {
        FROSCH_TIMER_START_LEVELID(removeDirichletNodesTime,"ConstantPartitionOfUnity::removeDirichletNodes");
        if (!dirichletBoundaryDofs.is_null()) {
            GOVec tmpDirichletBoundaryDofs(dirichletBoundaryDofs());
            sortunique(tmpDirichletBoundaryDofs);
            Volumes_->removeNodesWithDofs(tmpDirichletBoundaryDofs());
            Volumes_->removeEmptyEntities();
            Volumes_->setUniqueIDToFirstGlobalNodeID();
        }
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int ConstantPartitionOfUnity<SC,LO,GO,NO>::computePartitionOfUnity(ConstXMultiVectorPtr nodeList)
    {
        FROSCH_TIMER_START_LEVELID(computePartitionOfUnityTime,"ConstantPartitionOfUnity::computePartitionOfUnity");
        // Interface
        UN dofsPerNode = DDInterface_->getInterior()->getEntity(0)->getDofsPerNode();
        UN numInteriorDofs = dofsPerNode*DDInterface_->getInterior()->getEntity(0)->getNumNodes();
        
        if (UseVolumes_) Volumes_->buildEntityMap(DDInterface_->getNodesMap());
        
        if (this->Verbosity_==All) {
            // Count entities
            GOVec global(1);
            LOVec local(1);
            LOVec sum(1);
            SCVec avg(1);
            LOVec min(1);
            LOVec max(1);
            if (UseVolumes_) {
                global[0] = Volumes_->getEntityMap()->getMaxAllGlobalIndex();
                if (DDInterface_->getNodesMap()->lib()==UseEpetra || Volumes_->getEntityMap()->getGlobalNumElements()>0) {
                    global[0] += 1;
                }
                if (global[0]<0) global[0] = 0;
                local[0] = (LO) std::max((LO) Volumes_->getEntityMap()->getNodeNumElements(),(LO) 0);
                reduceAll(*this->MpiComm_,REDUCE_SUM,local[0],ptr(&sum[0]));
                avg[0] = std::max(sum[0]/double(this->MpiComm_->getSize()),0.0);
                reduceAll(*this->MpiComm_,REDUCE_MIN,local[0],ptr(&min[0]));
                reduceAll(*this->MpiComm_,REDUCE_MAX,local[0],ptr(&max[0]));
            } else {
                global[0] = -1;
                local[0] = -1;
                avg[0] = -1;
                min[0] = -1;
                max[0] = -1;
            }
            
            if (global[0]<0) {
                global[0] = -1;
            }
            
            if (this->Verbose_) {
                std::cout << "\n\
    ------------------------------------------------------------------------------\n\
     Volumes statistics\n\
    ------------------------------------------------------------------------------\n\
      Volumes:        total / avg / min / max     ---  " << global[0] << " / " << avg[0] << " / " << min[0] << " / " << max[0] << "\n\
    ------------------------------------------------------------------------------\n";
            }
        }

        // Maps
        if (UseVolumes_) {
            this->PartitionOfUnityMaps_[0] = Volumes_->getEntityMap();
        }

        if (this->Verbose_) {
            std::cout << std::boolalpha << "\n\
    ------------------------------------------------------------------------------\n\
     Constant Partition Of Unity \n\
    ------------------------------------------------------------------------------\n\
      Volumes                                     --- " << UseVolumes_ << "\n\
    ------------------------------------------------------------------------------\n" << std::noboolalpha;
        }

        // Build Partition Of Unity Vectors
        XMapPtr serialMap = MapFactory<LO,GO,NO>::Build(DDInterface_->getNodesMap()->lib(),numInteriorDofs,0,this->SerialComm_);

        if (UseVolumes_ && Volumes_->getNumEntities()>0) {
            XMultiVectorPtr tmpVector = MultiVectorFactory<SC,LO,GO,NO>::Build(serialMap,Volumes_->getNumEntities());

            for (UN i=0; i<Volumes_->getNumEntities(); i++) {
                for (UN j=0; j<Volumes_->getEntity(i)->getNumNodes(); j++) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        tmpVector->replaceLocalValue(Volumes_->getEntity(i)->getGammaDofID(j,k),i,ScalarTraits<SC>::one());
                    }
                }
            }

            this->LocalPartitionOfUnity_[0] = tmpVector;
        }
        return 0;
    }
}

#endif
