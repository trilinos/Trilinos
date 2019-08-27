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

#ifndef _FROSCH_RGDSWINTERFACEPARTITIONOFUNITY_DEF_HPP
#define _FROSCH_RGDSWINTERFACEPARTITIONOFUNITY_DEF_HPP

#include <FROSch_RGDSWInterfacePartitionOfUnity_decl.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    RGDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::RGDSWInterfacePartitionOfUnity(CommPtr mpiComm,
                                                                                CommPtr serialComm,
                                                                                UN dimension,
                                                                                UN dofsPerNode,
                                                                                ConstXMapPtr nodesMap,
                                                                                ConstXMapPtrVecPtr dofsMaps,
                                                                                ParameterListPtr parameterList,
                                                                                Verbosity verbosity,
                                                                                UN levelID) :
    GDSWInterfacePartitionOfUnity<SC,LO,GO,NO> (mpiComm,serialComm,dimension,dofsPerNode,nodesMap,dofsMaps,parameterList,verbosity,levelID),
    UseCoarseNodes_ (false),
    CoarseNodes_ (),
    EntitySetVector_ (),
    DistanceFunction_ (ConstantDistanceFunction)
    {
        FROSCH_TIMER_START_LEVELID(rGDSWInterfacePartitionOfUnityTime,"RGDSWInterfacePartitionOfUnity::RGDSWInterfacePartitionOfUnity");
        this->UseVertices_ = false;
        this->UseShortEdges_ = false;
        this->UseStraightEdges_ = false;
        this->UseEdges_ = false;
        this->UseFaces_ = false;

        if (!this->ParameterList_->get("Type","Full").compare("Full")) {
            UseCoarseNodes_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("CoarseNodes")) {
            UseCoarseNodes_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("Custom")) {
            UseCoarseNodes_ = this->ParameterList_->sublist("Custom").get("CoarseNodes",false);
        } else {
            FROSCH_ASSERT(false,"FROSch::RGDSWInterfacePartitionOfUnity : ERROR: Specify a valid Type.");
        }

        if (!this->ParameterList_->get("Distance Function","Constant").compare("Constant")) {
            DistanceFunction_ = ConstantDistanceFunction;
        } else if (!this->ParameterList_->get("Distance Function","Constant").compare("Inverse Euclidean")) {
            DistanceFunction_ = InverseEuclideanDistanceFunction;
        } else {
            FROSCH_ASSERT(false,"FROSch::RGDSWInterfacePartitionOfUnity : ERROR: Specify a valid Distance Function.");
        }
        this->LocalPartitionOfUnity_ = XMultiVectorPtrVecPtr(1);
        this->PartitionOfUnityMaps_ = XMapPtrVecPtr(1);
    }

    template <class SC,class LO,class GO,class NO>
    int RGDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::computePartitionOfUnity(ConstXMultiVectorPtr nodeList)
    {
        FROSCH_TIMER_START_LEVELID(computePartitionOfUnityTime,"RGDSWInterfacePartitionOfUnity::computePartitionOfUnity");
        // Interface
        UN dofsPerNode = this->DDInterface_->getInterface()->getEntity(0)->getDofsPerNode();
        UN numInterfaceDofs = dofsPerNode*this->DDInterface_->getInterface()->getEntity(0)->getNumNodes();

        this->DDInterface_->buildEntityHierarchy();

        this->DDInterface_->computeDistancesToCoarseNodes(this->DDInterface_->getDimension(),nodeList,DistanceFunction_);

        this->DDInterface_->buildEntityMaps(false,
                                            false,
                                            false,
                                            false,
                                            false,
                                            UseCoarseNodes_);

        EntitySetVector_ = this->DDInterface_->getEntitySetVector();

        // Map
        if (UseCoarseNodes_) {
            CoarseNodes_ = this->DDInterface_->getCoarseNodes();
            this->PartitionOfUnityMaps_[0] = CoarseNodes_->getEntityMap();
        }

        if (this->MpiComm_->getRank() == 0) {
            std::cout << std::boolalpha << "\n\
    ------------------------------------------------------------------------------\n\
     RGDSW Interface Partition Of Unity (RGDSW IPOU)\n\
    ------------------------------------------------------------------------------\n\
      Coarse nodes                               --- " << UseCoarseNodes_ << "\n\
    ------------------------------------------------------------------------------\n" << std::noboolalpha;
        }

        // Build Partition Of Unity Vectors
        XMapPtr serialInterfaceMap = MapFactory<LO,GO,NO>::Build(this->DDInterface_->getNodesMap()->lib(),numInterfaceDofs,0,this->SerialComm_);

        if (UseCoarseNodes_ && CoarseNodes_->getNumEntities()>0) {
            XMultiVectorPtr tmpVector = MultiVectorFactory<SC,LO,GO,NO>::Build(serialInterfaceMap,CoarseNodes_->getNumEntities());

            // Loop over EntitySetVector_
            for (UN i=0; i<EntitySetVector_.size(); i++) {
                // Loop over entities
                for (UN j=0; j<EntitySetVector_[i]->getNumEntities(); j++) {
                    InterfaceEntityPtr tmpEntity = EntitySetVector_[i]->getEntity(j);
                    LO coarseNodeID = tmpEntity->getCoarseNodeID();
                    UN numCoarseNodes = tmpEntity->getCoarseNodes()->getNumEntities();
                    if (coarseNodeID==-1) {
                        FROSCH_ASSERT(numCoarseNodes!=0,"coarseNodeID==-1 but numCoarseNodes==0!");
                        for (UN m=0; m<numCoarseNodes; m++) {
                            InterfaceEntityPtr tmpCoarseNode = tmpEntity->getCoarseNodes()->getEntity(m);
                            LO index = tmpCoarseNode->getCoarseNodeID();
                            // Offspring: loop over nodes
                            for (UN l=0; l<tmpEntity->getNumNodes(); l++) {
                                SC value = tmpEntity->getDistanceToCoarseNode(l,m)/tmpEntity->getDistanceToCoarseNode(l,numCoarseNodes);
                                for (UN k=0; k<dofsPerNode; k++) {
                                    tmpVector->replaceLocalValue(tmpEntity->getGammaDofID(l,k),index,value*ScalarTraits<SC>::one());
                                }
                            }
                        }
                    } else {
                        // Coarse node: loop over nodes
                        for (UN l=0; l<EntitySetVector_[i]->getEntity(j)->getNumNodes(); l++) {
                            for (UN k=0; k<dofsPerNode; k++) {
                                tmpVector->replaceLocalValue(tmpEntity->getGammaDofID(l,k),coarseNodeID,ScalarTraits<SC>::one());
                            }
                        }
                    }
                }
            }
            this->LocalPartitionOfUnity_[0] = tmpVector;
        }

        return 0;
    }

}

#endif
