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

#ifndef _FROSCH_RGDSWCOARSEOPERATOR_DEF_HPP
#define _FROSCH_RGDSWCOARSEOPERATOR_DEF_HPP

#include <FROSch_RGDSWCoarseOperator_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    RGDSWCoarseOperator<SC,LO,GO,NO>::RGDSWCoarseOperator(ConstXMatrixPtr k,
                                                          ParameterListPtr parameterList) :
    GDSWCoarseOperator<SC,LO,GO,NO> (k,parameterList)
    {
        FROSCH_DETAILTIMER_START_LEVELID(rGDSWCoarseOperatorTime,"RGDSWCoarseOperator::RGDSWCoarseOperator");
    }

    template <class SC,class LO,class GO,class NO>
    int RGDSWCoarseOperator<SC,LO,GO,NO>::resetCoarseSpaceBlock(UN blockId,
                                                                UN dimension,
                                                                UN dofsPerNode,
                                                                ConstXMapPtr nodesMap,
                                                                ConstXMapPtrVecPtr dofsMaps,
                                                                GOVecPtr dirichletBoundaryDofs,
                                                                ConstXMultiVectorPtr nodeList)
    {
        FROSCH_DETAILTIMER_START_LEVELID(resetCoarseSpaceBlockTime,"RGDSWCoarseOperator::resetCoarseSpaceBlock");
        FROSCH_ASSERT(dofsMaps.size()==dofsPerNode,"dofsMaps.size()!=dofsPerNode");
        FROSCH_ASSERT(blockId<this->NumberOfBlocks_,"Block does not exist yet and can therefore not be reset.");

        if (!this->DistributionList_->get("Type","linear").compare("ZoltanDual")) {
            FROSCH_ASSERT(false,"RGDSWCoarseOperator:: Distribution Type ZoltanDual only works for IPOUHarmonicCoarseOperator");
        }

        // Process the parameter list
        stringstream blockIdStringstream;
        blockIdStringstream << blockId+1;
        string blockIdString = blockIdStringstream.str();
        RCP<ParameterList> coarseSpaceList = sublist(sublist(this->ParameterList_,"Blocks"),blockIdString.c_str());

        CommunicationStrategy communicationStrategy = CreateOneToOneMap;
        if (!coarseSpaceList->get("Interface Communication Strategy","CreateOneToOneMap").compare("CrsMatrix")) {
            communicationStrategy = CommCrsMatrix;
        } else if (!coarseSpaceList->get("Interface Communication Strategy","CreateOneToOneMap").compare("CrsGraph")) {
            communicationStrategy = CommCrsGraph;
        } else if (!coarseSpaceList->get("Interface Communication Strategy","CreateOneToOneMap").compare("CreateOneToOneMap")) {
            communicationStrategy = CreateOneToOneMap;
        } else {
            FROSCH_ASSERT(false,"FROSch::RGDSWCoarseOperator : ERROR: Specify a valid communication strategy for the identification of the interface components.");
        }

        Verbosity verbosity = All;
        if (!coarseSpaceList->get("Verbosity","All").compare("None")) {
            verbosity = None;
        } else if (!coarseSpaceList->get("Verbosity","All").compare("All")) {
            verbosity = All;
        } else {
            FROSCH_ASSERT(false,"FROSch::RGDSWCoarseOperator : ERROR: Specify a valid verbosity level.");
        }

        bool useForCoarseSpace = coarseSpaceList->get("Use For Coarse Space",false);
        string option = coarseSpaceList->get("Option","1");
        DistanceFunction distanceFunction = ConstantDistanceFunction;
        if (!option.compare("1")) {

        } else if (!option.compare("2.2")) {
            distanceFunction = InverseEuclideanDistanceFunction;
        } else {
            FROSCH_ASSERT(false,"Option is unknown!");
        }

        bool useRotations = coarseSpaceList->get("Rotations",true);
        if (useRotations && nodeList.is_null()) {
            //FROSCH_ASSERT(option==1,"Only option 1 can be constructed without a valid node list.");
            useRotations = false;
            FROSCH_WARNING("FROSch::RGDSWCoarseOperator",this->Verbose_,"Rotations cannot be used since nodeList.is_null().");
        }

        this->DofsMaps_[blockId] = dofsMaps;
        this->DofsPerNode_[blockId] = dofsPerNode;

        Array<GO> tmpDirichletBoundaryDofs(dirichletBoundaryDofs()); // Here, we do a copy. Maybe, this is not necessary
        sortunique(tmpDirichletBoundaryDofs);

        this->DDInterface_.reset(new DDInterface<SC,LO,GO,NO>(dimension,this->DofsPerNode_[blockId],nodesMap.getConst(),verbosity,this->LevelID_,communicationStrategy));
        this->DDInterface_->resetGlobalDofs(dofsMaps);
        this->DDInterface_->removeDirichletNodes(tmpDirichletBoundaryDofs);

        EntitySetPtr interface = this->DDInterface_->getInterface();
        EntitySetPtr interior = this->DDInterface_->getInterior();

        if (useForCoarseSpace) {

            if (this->Verbose_) {
                cout
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << setw(89) << "-----------------------------------------------------------------------------------------"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| "
                << left << setw(74) << "RGDSWCoarseOperator " << right << setw(8) << "(Level " << setw(2) << this->LevelID_ << ")"
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << setw(89) << "========================================================================================="
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| " << left << setw(41) << "Block" << right
                << " | " << setw(41) << blockId
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| " << left << setw(41) << "Numer of degrees of freedom per node" << right
                << " | " << setw(41) << dimension
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| " << left << setw(41) << "Numer of degrees of freedom per node" << right
                << " | " << setw(41) << dofsPerNode
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| " << left << setw(41) << "Option" << right
                << " | " << setw(41) << option
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << setw(89) << "-----------------------------------------------------------------------------------------"
                << endl;
            }

            // Check for interface
            if (interface->getEntity(0)->getNumNodes()==0) {
                FROSCH_NOTIFICATION("FROSch::RGDSWCoarseOperator",this->Verbose_,"No interface found => Volume functions will be used instead.");
                this->computeVolumeFunctions(blockId,dimension,nodesMap,nodeList,interior);
            } else {
                this->GammaDofs_[blockId] = LOVecPtr(this->DofsPerNode_[blockId]*interface->getEntity(0)->getNumNodes());
                this->IDofs_[blockId] = LOVecPtr(this->DofsPerNode_[blockId]*interior->getEntity(0)->getNumNodes());
                for (UN k=0; k<this->DofsPerNode_[blockId]; k++) {
                    for (UN i=0; i<interface->getEntity(0)->getNumNodes(); i++) {
                        this->GammaDofs_[blockId][this->DofsPerNode_[blockId]*i+k] = interface->getEntity(0)->getLocalDofID(i,k);
                    }
                    for (UN i=0; i<interior->getEntity(0)->getNumNodes(); i++) {
                        this->IDofs_[blockId][this->DofsPerNode_[blockId]*i+k] = interior->getEntity(0)->getLocalDofID(i,k);
                    }
                }

                this->InterfaceCoarseSpaces_[blockId].reset(new CoarseSpace<SC,LO,GO,NO>(this->MpiComm_,this->SerialComm_));

                if (this->ParameterList_->get("Test Unconnected Interface",true)) {
                    this->DDInterface_->divideUnconnectedEntities(this->K_);
                }

                this->DDInterface_->buildEntityHierarchy();

                this->DDInterface_->computeDistancesToRoots(dimension,nodeList,distanceFunction);

                //////////////////////////
                // Root Basis Functions //
                //////////////////////////
                EntitySetPtrVecPtr entitySetVector = this->DDInterface_->getEntitySetVector();
                this->DDInterface_->buildEntityMaps(false,
                                                    false,
                                                    false,
                                                    false,
                                                    false,
                                                    true);

                XMultiVectorPtrVecPtr translations = this->computeTranslations(blockId,this->DDInterface_->getRoots(),entitySetVector,distanceFunction);
                for (UN i=0; i<translations.size(); i++) {
                    this->InterfaceCoarseSpaces_[blockId]->addSubspace(this->DDInterface_->getRoots()->getEntityMap(),null,translations[i]);
                }

                if (useRotations) {
                    XMultiVectorPtrVecPtr rotations = this->computeRotations(blockId,dimension,nodeList,this->DDInterface_->getRoots(),entitySetVector,distanceFunction);
                    for (UN i=0; i<rotations.size(); i++) {                        this->InterfaceCoarseSpaces_[blockId]->addSubspace(this->DDInterface_->getRoots()->getEntityMap(),null,rotations[i]);
                    }
                }

                this->InterfaceCoarseSpaces_[blockId]->assembleCoarseSpace();

                if (this->Verbose_) {
                    cout
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << setw(89) << "-----------------------------------------------------------------------------------------"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << "| "
                    << left << setw(74) << "> RGDSW coarse space " << right << setw(8) << "(Level " << setw(2) << this->LevelID_ << ")"
                    << " |"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << setw(89) << "========================================================================================="
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << "| " << left << setw(20) << "Coarse nodes " << " | " << setw(19) << " Translations" << right
                    << " | " << setw(41) << boolalpha << useForCoarseSpace << noboolalpha
                    << " |"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << "| " << left << setw(20) << "Coarse nodes " << " | " << setw(19) << " Rotations" << right
                    << " | " << setw(41) << boolalpha << useRotations << noboolalpha
                    << " |"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << setw(89) << "-----------------------------------------------------------------------------------------"
                    << endl;
                }
            }
        }
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    typename RGDSWCoarseOperator<SC,LO,GO,NO>::XMapPtr RGDSWCoarseOperator<SC,LO,GO,NO>::BuildRepeatedMapCoarseLevel(ConstXMapPtr &nodesMap,
                                                UN dofsPerNode,
                                                ConstXMapPtrVecPtr dofsMaps,
                                               UN partitionType)
    {
        FROSCH_ASSERT(false,"For RGDSWCoarseOperator the ZoltanDual Option is not implemented!");
    }

    template <class SC,class LO,class GO,class NO>
    typename RGDSWCoarseOperator<SC,LO,GO,NO>::XMultiVectorPtrVecPtr RGDSWCoarseOperator<SC,LO,GO,NO>::computeTranslations(UN blockId,
                                                                                                                           EntitySetPtr Roots,
                                                                                                                           EntitySetPtrVecPtr entitySetVector,
                                                                                                                           DistanceFunction distanceFunction)
    {
        FROSCH_DETAILTIMER_START_LEVELID(computeTranslationsTime,"RGDSWCoarseOperator::computeTranslations");
        XMultiVectorPtrVecPtr translations(this->DofsPerNode_[blockId]);
        XMapPtr serialGammaMap = MapFactory<LO,GO,NO>::Build(this->K_->getRangeMap()->lib(),this->GammaDofs_[blockId].size(),0,this->SerialComm_);
        for (UN i=0; i<this->DofsPerNode_[blockId]; i++) {
            if (Roots->getNumEntities()>0) {
                translations[i] = MultiVectorFactory<SC,LO,GO,NO>::Build(serialGammaMap,Roots->getNumEntities());
            } else {
                translations[i] = null;
            }
        }

        // Loop over Dofs
        for (UN k=0; k<this->DofsPerNode_[blockId]; k++) {
            // Loop over entitySetVector
            for (UN i=0; i<entitySetVector.size(); i++) {
                // Loop over entities
                for (UN j=0; j<entitySetVector[i]->getNumEntities(); j++) {
                    InterfaceEntityPtr tmpEntity = entitySetVector[i]->getEntity(j);
                    LO rootID = tmpEntity->getRootID();
                    UN numRoots = tmpEntity->getRoots()->getNumEntities();
                    if (rootID==-1) {
                        //if (numRoots==0) cout << rootID << " " << numRoots << " " << tmpEntity->getAncestors()->getNumEntities() << endl;
                        FROSCH_ASSERT(numRoots!=0,"rootID==-1 but numRoots==0!");
                        for (UN m=0; m<numRoots; m++) {
                            InterfaceEntityPtr tmpRoot = tmpEntity->getRoots()->getEntity(m);
                            LO index = tmpRoot->getRootID();
                            // Offspring: loop over nodes
                            for (UN l=0; l<tmpEntity->getNumNodes(); l++) {
                                SC value = tmpEntity->getDistanceToRoot(l,m)/tmpEntity->getDistanceToRoot(l,numRoots);
                                translations[k]->replaceLocalValue(tmpEntity->getGammaDofID(l,k),index,value);
                            }
                        }
                    } else {
                        // Coarse node: loop over nodes
                        for (UN l=0; l<entitySetVector[i]->getEntity(j)->getNumNodes(); l++) {
                            translations[k]->replaceLocalValue(tmpEntity->getGammaDofID(l,k),rootID,ScalarTraits<SC>::one());
                        }
                    }
                }
            }
        }
        return translations;
    }

    template <class SC,class LO,class GO,class NO>
    typename RGDSWCoarseOperator<SC,LO,GO,NO>::XMultiVectorPtrVecPtr RGDSWCoarseOperator<SC,LO,GO,NO>::computeRotations(UN blockId,
                                                                                                                        UN dimension,
                                                                                                                        ConstXMultiVectorPtr nodeList,
                                                                                                                        EntitySetPtr Roots,
                                                                                                                        EntitySetPtrVecPtr entitySetVector,
                                                                                                                        DistanceFunction distanceFunction)
    {
        FROSCH_DETAILTIMER_START_LEVELID(computeRotationsTime,"RGDSWCoarseOperator::computeRotations");
        FROSCH_ASSERT(nodeList->getNumVectors()==dimension,"dimension of the nodeList is wrong.");
        FROSCH_ASSERT(dimension==this->DofsPerNode_[blockId],"dimension!=this->DofsPerNode_[blockId]");

        UN rotationsPerEntity = 0;
        switch (dimension) {
            case 1:
                return null;
                break;
            case 2:
                rotationsPerEntity = 1;
                break;
            case 3:
                rotationsPerEntity = 3;
                break;
            default:
                FROSCH_ASSERT(false,"The dimension is neither 2 nor 3!");
                break;
        }

        XMultiVectorPtrVecPtr rotations(rotationsPerEntity);
        XMapPtr serialGammaMap = MapFactory<LO,GO,NO>::Build(this->K_->getRangeMap()->lib(),this->GammaDofs_[blockId].size(),0,this->SerialComm_);
        for (UN i=0; i<rotationsPerEntity; i++) {
            if (Roots->getNumEntities()>0) {
                rotations[i] = MultiVectorFactory<SC,LO,GO,NO>::Build(serialGammaMap,Roots->getNumEntities());
            } else {
                rotations[i] = null;
            }
        }

        SC x,y,z,rx,ry,rz;
        // Loop over entitySetVector
        for (UN i=0; i<entitySetVector.size(); i++) {
            // Loop over entities
            for (UN j=0; j<entitySetVector[i]->getNumEntities(); j++) {
                InterfaceEntityPtr tmpEntity = entitySetVector[i]->getEntity(j);
                LO rootID = tmpEntity->getRootID();
                UN numRoots = tmpEntity->getRoots()->getNumEntities();
                if (rootID==-1) {
                    FROSCH_ASSERT(numRoots!=0,"rootID==-1 but numRoots==0!");
                    for (UN m=0; m<numRoots; m++) {
                        InterfaceEntityPtr tmpRoot = tmpEntity->getRoots()->getEntity(m);
                        LO index = tmpRoot->getRootID();
                        // Offspring: loop over nodes
                        for (UN l=0; l<tmpEntity->getNumNodes(); l++) {
                            SC value = tmpEntity->getDistanceToRoot(l,m)/tmpEntity->getDistanceToRoot(l,numRoots);

                            // Rotations
                            x = nodeList->getData(0)[tmpEntity->getLocalNodeID(l)];
                            y = nodeList->getData(1)[tmpEntity->getLocalNodeID(l)];

                            // Rotation 1
                            rx = y;
                            ry = -x;
                            rz = 0;
                            rotations[0]->replaceLocalValue(tmpEntity->getGammaDofID(l,0),index,value*rx);
                            rotations[0]->replaceLocalValue(tmpEntity->getGammaDofID(l,1),index,value*ry);
                            if (dimension == 3) {
                                z = nodeList->getData(2)[tmpEntity->getLocalNodeID(l)];

                                rotations[0]->replaceLocalValue(tmpEntity->getGammaDofID(l,2),index,value*rz);

                                // Rotation 2
                                rx = -z;
                                ry = 0;
                                rz = x;
                                rotations[1]->replaceLocalValue(tmpEntity->getGammaDofID(l,0),index,value*rx);
                                rotations[1]->replaceLocalValue(tmpEntity->getGammaDofID(l,1),index,value*ry);
                                rotations[1]->replaceLocalValue(tmpEntity->getGammaDofID(l,2),index,value*rz);

                                // Rotation 3
                                rx = 0;
                                ry = z;
                                rz = -y;
                                rotations[2]->replaceLocalValue(tmpEntity->getGammaDofID(l,0),index,value*rx);
                                rotations[2]->replaceLocalValue(tmpEntity->getGammaDofID(l,1),index,value*ry);
                                rotations[2]->replaceLocalValue(tmpEntity->getGammaDofID(l,2),index,value*rz);
                            }
                        }
                    }
                } else {
                    // Coarse node: loop over nodes
                    for (UN l=0; l<entitySetVector[i]->getEntity(j)->getNumNodes(); l++) {
                        // Rotations
                        x = nodeList->getData(0)[tmpEntity->getLocalNodeID(l)];
                        y = nodeList->getData(1)[tmpEntity->getLocalNodeID(l)];

                        // Rotation 1
                        rx = y;
                        ry = -x;
                        rz = 0;
                        rotations[0]->replaceLocalValue(tmpEntity->getGammaDofID(l,0),rootID,rx);
                        rotations[0]->replaceLocalValue(tmpEntity->getGammaDofID(l,1),rootID,ry);
                        if (dimension == 3) {
                            z = nodeList->getData(2)[tmpEntity->getLocalNodeID(l)];

                            rotations[0]->replaceLocalValue(tmpEntity->getGammaDofID(l,2),rootID,rz);

                            // Rotation 2
                            rx = -z;
                            ry = 0;
                            rz = x;
                            rotations[1]->replaceLocalValue(tmpEntity->getGammaDofID(l,0),rootID,rx);
                            rotations[1]->replaceLocalValue(tmpEntity->getGammaDofID(l,1),rootID,ry);
                            rotations[1]->replaceLocalValue(tmpEntity->getGammaDofID(l,2),rootID,rz);

                            // Rotation 3
                            rx = 0;
                            ry = z;
                            rz = -y;
                            rotations[2]->replaceLocalValue(tmpEntity->getGammaDofID(l,0),rootID,rx);
                            rotations[2]->replaceLocalValue(tmpEntity->getGammaDofID(l,1),rootID,ry);
                            rotations[2]->replaceLocalValue(tmpEntity->getGammaDofID(l,2),rootID,rz);
                        }
                    }
                }
            }
        }
        return rotations;
    }
}

#endif
