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

    using namespace std;
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
        FROSCH_DETAILTIMER_START_LEVELID(constantPartitionOfUnityTime,"ConstantPartitionOfUnity::ConstantPartitionOfUnity");

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

        this->LocalPartitionOfUnity_ = ConstXMultiVectorPtrVecPtr(1);
        this->PartitionOfUnityMaps_ = ConstXMapPtrVecPtr(1);
    }

    template <class SC,class LO,class GO,class NO>
    ConstantPartitionOfUnity<SC,LO,GO,NO>::~ConstantPartitionOfUnity()
    {

    }

    template <class SC,class LO,class GO,class NO>
    int ConstantPartitionOfUnity<SC,LO,GO,NO>::removeDirichletNodes(GOVecView dirichletBoundaryDofs,
                                                                    ConstXMultiVectorPtr nodeList)
    {
        FROSCH_DETAILTIMER_START_LEVELID(removeDirichletNodesTime,"ConstantPartitionOfUnity::removeDirichletNodes");
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
        FROSCH_DETAILTIMER_START_LEVELID(computePartitionOfUnityTime,"ConstantPartitionOfUnity::computePartitionOfUnity");
        // Interface
        UN dofsPerNode = DDInterface_->getInterior()->getEntity(0)->getDofsPerNode();
        UN numInteriorDofs = dofsPerNode*DDInterface_->getInterior()->getEntity(0)->getNumNodes();

        if (UseVolumes_) Volumes_->buildEntityMap(DDInterface_->getNodesMap());

        if (this->Verbosity_==All) {
            FROSCH_DETAILTIMER_START_LEVELID(printStatisticsTime,"print statistics");
            // Count entities
            GOVec globalVec(1);
            LOVec localVec(1);
            LOVec sumVec(1);
            SCVec avgVec(1);
            LOVec minVec(1);
            LOVec maxVec(1);
            if (UseVolumes_) {
                globalVec[0] = Volumes_->getEntityMap()->getMaxAllGlobalIndex();
                if (DDInterface_->getNodesMap()->lib()==UseEpetra || Volumes_->getEntityMap()->getGlobalNumElements()>0) {
                    globalVec[0] += 1;
                }
                if (globalVec[0]<0) globalVec[0] = 0;
                localVec[0] = (LO) max((LO) Volumes_->getEntityMap()->getNodeNumElements(),(LO) 0);
                reduceAll(*this->MpiComm_,REDUCE_SUM,localVec[0],ptr(&sumVec[0]));
                avgVec[0] = max(sumVec[0]/double(this->MpiComm_->getSize()),0.0);
                reduceAll(*this->MpiComm_,REDUCE_MIN,localVec[0],ptr(&minVec[0]));
                reduceAll(*this->MpiComm_,REDUCE_MAX,localVec[0],ptr(&maxVec[0]));
            } else {
                globalVec[0] = -1;
                localVec[0] = -1;
                avgVec[0] = -1;
                minVec[0] = -1;
                maxVec[0] = -1;
            }

            if (globalVec[0]<0) {
                globalVec[0] = -1;
            }

            if (this->Verbose_) {
                cout
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << setw(89) << "-----------------------------------------------------------------------------------------"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| "
                << left << setw(74) << "> Volumes statistics " << right << setw(8) << "(Level " << setw(2) << this->LevelID_ << ")" << right
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << setw(89) << "========================================================================================="
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| " << left << setw(20) << " " << right
                << " | " << setw(10) << "total"
                << " | " << setw(10) << "avg"
                << " | " << setw(10) << "min"
                << " | " << setw(10) << "max"
                << " | " << setw(10) << "global sum"
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << setw(89) << "-----------------------------------------------------------------------------------------"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| " << left << setw(20) << "Volumes" << right
                << " | "; globalVec[0]<0 ? cout << setw(10) << " " : cout << setw(10) << globalVec[0]; cout
                << " | "; avgVec[0]<0 ? cout << setw(10) << " " : cout << setw(10) << setprecision(5) << avgVec[0]; cout
                << " | "; minVec[0]<0 ? cout << setw(10) << " " : cout << setw(10) << minVec[0]; cout
                << " | "; maxVec[0]<0 ? cout << setw(10) << " " : cout << setw(10) << maxVec[0]; cout
                << " | "; sumVec[0]<0 ? cout << setw(10) << " " : cout << setw(10) << sumVec[0]; cout
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << setw(89) << "-----------------------------------------------------------------------------------------"
                << endl;
            }
        }

        // Maps
        if (UseVolumes_) {
            this->PartitionOfUnityMaps_[0] = Volumes_->getEntityMap();
        }

        if (this->Verbose_) {
            cout
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << setw(89) << "-----------------------------------------------------------------------------------------"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << "| "
            << left << setw(74) << "> Constant Partition Of Unity " << right << setw(8) << "(Level " << setw(2) << this->LevelID_ << ")"
            << " |"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << setw(89) << "========================================================================================="
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << "| " << left << setw(41) << "Volumes" << right
            << " | " << setw(41) << boolalpha << UseVolumes_ << noboolalpha
            << " |"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << setw(89) << "-----------------------------------------------------------------------------------------"
            << endl;
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
