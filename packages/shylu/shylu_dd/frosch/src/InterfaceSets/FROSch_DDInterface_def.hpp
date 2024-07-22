// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_DDINTERFACE_DEF_HPP
#define _FROSCH_DDINTERFACE_DEF_HPP

#include <FROSch_DDInterface_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    DDInterface<SC,LO,GO,NO>::DDInterface(UN dimension,
                                          UN dofsPerNode,
                                          ConstXMapPtr localToGlobalMap,
                                          Verbosity verbosity,
                                          UN levelID,
                                          CommunicationStrategy commStrategy) :
    MpiComm_ (localToGlobalMap->getComm()),
    Dimension_ (dimension),
    DofsPerNode_ (dofsPerNode),
    NumMyNodes_ (localToGlobalMap->getLocalNumElements()),
    NodesMap_ (localToGlobalMap),
    CommStrategy_ (commStrategy),
    Verbose_ (MpiComm_->getRank()==0),
    Verbosity_ (verbosity),
    LevelID_ (levelID)
    {
        FROSCH_DETAILTIMER_START_LEVELID(dDInterfaceTime,"DDInterface::DDInterface");
        FROSCH_ASSERT(((Dimension_==2)||(Dimension_==3)),"FROSch::DDInterface: Only dimension 2 and 3 are available");

        //if (Verbose_ && Verbosity_==All) cout << "FROSch::DDInterface" << endl;

        IntVecVecPtr componentsSubdomains;
        IntVecVec componentsSubdomainsUnique;

        communicateLocalComponents(componentsSubdomains,componentsSubdomainsUnique);

        identifyLocalComponents(componentsSubdomains,componentsSubdomainsUnique);
    }

    template <class SC,class LO,class GO,class NO>
    DDInterface<SC,LO,GO,NO>::~DDInterface()
    {

    } // Do we need sth here?

    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::resetGlobalDofs(ConstXMapPtrVecPtr dofsMaps)
    {
        FROSCH_DETAILTIMER_START_LEVELID(resetGlobalDofsTime,"DDInterface::resetGlobalDofs");
        //if (Verbose_ && Verbosity_==All) cout << "FROSch::DDInterface : Resetting Global IDs" << endl;

        // EntityVector
        for (UN l=0; l<EntitySetVector_.size(); l++) {
            for (UN i=0; i<EntitySetVector_[l]->getNumEntities(); i++) {
                for (UN j=0; j<EntitySetVector_[l]->getEntity(i)->getNumNodes(); j++) {
                    LO localID = EntitySetVector_[l]->getEntity(i)->getLocalNodeID(j);
                    UNVecPtr dofIDs(DofsPerNode_);
                    GOVecPtr dofsGlobal(DofsPerNode_);
                    for (UN k=0; k<DofsPerNode_; k++) {
                        dofIDs[k] = k;
                        dofsGlobal[k] = dofsMaps[k]->getGlobalElement(localID);
                    }
                    EntitySetVector_[l]->getEntity(i)->resetGlobalDofs(j,DofsPerNode_,&(dofIDs[0]),&(dofsGlobal[0]));
                }
            }
        }

        // Interface
        for (UN i=0; i<Interface_->getNumEntities(); i++) {
            for (UN j=0; j<Interface_->getEntity(i)->getNumNodes(); j++) {
                LO localID = Interface_->getEntity(i)->getLocalNodeID(j);
                UNVecPtr dofIDs(DofsPerNode_);
                GOVecPtr dofsGlobal(DofsPerNode_);
                for (UN k=0; k<DofsPerNode_; k++) {
                    dofIDs[k] = k;
                    dofsGlobal[k] = dofsMaps[k]->getGlobalElement(localID);
                }
                Interface_->getEntity(i)->resetGlobalDofs(j,DofsPerNode_,&(dofIDs[0]),&(dofsGlobal[0]));
            }
        }

        // Interior
        for (UN i=0; i<Interior_->getNumEntities(); i++) {
            for (UN j=0; j<Interior_->getEntity(i)->getNumNodes(); j++) {
                LO localID = Interior_->getEntity(i)->getLocalNodeID(j);
                UNVecPtr dofIDs(DofsPerNode_);
                GOVecPtr dofsGlobal(DofsPerNode_);
                for (UN k=0; k<DofsPerNode_; k++) {
                    dofIDs[k] = k;
                    dofsGlobal[k] = dofsMaps[k]->getGlobalElement(localID);
                }
                Interior_->getEntity(i)->resetGlobalDofs(j,DofsPerNode_,&(dofIDs[0]),&(dofsGlobal[0]));
            }
        }

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::removeDirichletNodes(GOVecView dirichletBoundaryDofs)
    {
        FROSCH_DETAILTIMER_START_LEVELID(removeDirichletNodesTime,"DDInterface::removeDirichletNodes");
        //if (Verbose_ && Verbosity_==All) cout << "FROSch::DDInterface : Removing Dirichlet Nodes from the domain decomposition interface" << endl;

        // EntityVector
        for (UN l=0; l<EntitySetVector_.size(); l++) {
            EntitySetVector_[l]->removeNodesWithDofs(dirichletBoundaryDofs);
        }
        removeEmptyEntities();
        for (UN l=0; l<EntitySetVector_.size(); l++) {
            EntitySetVector_[l]->setUniqueIDToFirstGlobalNodeID();
        }
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::divideUnconnectedEntities(ConstXMatrixPtr matrix)
    {
        FROSCH_DETAILTIMER_START_LEVELID(divideUnconnectedEntitiesTime,"DDInterface::divideUnconnectedEntities");
        //if (Verbose_ && Verbosity_==All) cout << "FROSch::DDInterface : Decomposing unconnected interface components" << endl;

        GOVecPtr indicesGammaDofs(DofsPerNode_*Interface_->getEntity(0)->getNumNodes());
        for (UN k=0; k<DofsPerNode_; k++) {
            for (UN i=0; i<Interface_->getEntity(0)->getNumNodes(); i++) {
                indicesGammaDofs[Interface_->getEntity(0)->getGammaDofID(i,k)] = Interface_->getEntity(0)->getGlobalDofID(i,k);
            }
        }

        const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();
        XMapPtr map = MapFactory<LO,GO,NO>::Build(matrix->getRowMap()->lib(),INVALID,indicesGammaDofs(),0,MpiComm_);
        matrix = FROSch::ExtractLocalSubdomainMatrix(matrix.getConst(),map.getConst(),ScalarTraits<SC>::one());

        // Operate on hierarchy
        for (UN i=0; i<EntitySetVector_.size(); i++) {
            EntitySetVector_[i]->divideUnconnectedEntities(matrix,MpiComm_->getRank());
        }

        /*
        LO numSeparateEdges = Edges_->divideUnconnectedEntities(matrix,MpiComm_->getRank());
        LO numSeparateFaces = Faces_->divideUnconnectedEntities(matrix,MpiComm_->getRank());

        if (Verbose_ && Verbosity_==All) {
            cout << "\n\
            --------------------------------------------\n\
            # separate edges:     --- " << numSeparateEdges << "\n\
            # separate faces:     --- " << numSeparateFaces << "\n\
            --------------------------------------------\n";
        }
        */

        removeEmptyEntities();

        // We need to set the unique ID; otherwise, we cannot sort entities
        for (UN i=0; i<EntitySetVector_.size(); i++) {
            EntitySetVector_[i]->setUniqueIDToFirstGlobalNodeID();
        }
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::flagEntities(ConstXMultiVectorPtr nodeList)
    {
        FROSCH_DETAILTIMER_START_LEVELID(flagEntitiesTime,"DDInterface::flagEntities");
        for (UN l=0; l<EntitySetVector_.size(); l++) {
            EntitySetVector_[l]->flagNodes();
            EntitySetVector_[l]->flagShortEntities();
        }
        if (!nodeList.is_null()) {
            for (UN l=0; l<EntitySetVector_.size(); l++) {
                EntitySetVector_[l]->flagStraightEntities(Dimension_,nodeList);
            }
        }
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::removeEmptyEntities()
    {
        FROSCH_DETAILTIMER_START_LEVELID(removeEmptyEntitiesTime,"DDInterface::removeEmptyEntities");
        //if (Verbose_ && Verbosity_==All) cout << "FROSch::DDInterface : Removing empty interface components" << endl;

        for (UN l=0; l<EntitySetVector_.size(); l++) {
            EntitySetVector_[l]->removeEmptyEntities();
        }
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::sortVerticesEdgesFaces(ConstXMultiVectorPtr nodeList)
    {
        FROSCH_DETAILTIMER_START_LEVELID(sortVerticesEdgesFacesTime,"DDInterface::sortVerticesEdgesFaces");
        //if (Verbose_ && Verbosity_==All) cout << "FROSch::DDInterface : Sorting interface components" << endl;

        // Clear EntitySets if non-empty
        if (Vertices_->getNumEntities()>0) Vertices_.reset(new EntitySet<SC,LO,GO,NO>(VertexType));
        if (ShortEdges_->getNumEntities()>0) ShortEdges_.reset(new EntitySet<SC,LO,GO,NO>(EdgeType));
        if (StraightEdges_->getNumEntities()>0) StraightEdges_.reset(new EntitySet<SC,LO,GO,NO>(EdgeType));
        if (Edges_->getNumEntities()>0) Edges_.reset(new EntitySet<SC,LO,GO,NO>(EdgeType));
        if (Faces_->getNumEntities()>0) Faces_.reset(new EntitySet<SC,LO,GO,NO>(FaceType));

        flagEntities(nodeList);

        // Make sure that we do not sort any empty entities
        removeEmptyEntities();

        for (UN l=0; l<EntitySetVector_.size(); l++) {
            switch (l) {
                case 0:
                    FROSCH_ASSERT(EntitySetVector_[l]->getNumEntities()==0,"FROSch::DDInterface: This case is impossible.");
                    break;
                case 1:
                    FROSCH_ASSERT(EntitySetVector_[l]->getNumEntities()==0,"FROSch::DDInterface: In this case, the entity is interior to the subdomain.");
                    break;
                case 2:
                    for (UN i=0; i<EntitySetVector_[l]->getNumEntities(); i++) {
                        switch (EntitySetVector_[l]->getEntity(i)->getEntityFlag()) {
                            case DefaultFlag: // By default, an entity which belongs to 2 subdomains is a face
                                EntitySetVector_[l]->getEntity(i)->resetEntityType(FaceType);
                                Faces_->addEntity(EntitySetVector_[l]->getEntity(i));
                                break;
                            case StraightFlag: // If an entity is straight, it is always a straight edge
                                EntitySetVector_[l]->getEntity(i)->resetEntityType(EdgeType);
                                StraightEdges_->addEntity(EntitySetVector_[l]->getEntity(i));
                                break;
                            case ShortFlag: // If an entity is a short, it is always a short edge
                                EntitySetVector_[l]->getEntity(i)->resetEntityType(EdgeType);
                                ShortEdges_->addEntity(EntitySetVector_[l]->getEntity(i));
                                break;
                            case NodeFlag: // If an entity is a node, it is always a vertex
                                EntitySetVector_[l]->getEntity(i)->resetEntityType(VertexType);
                                Vertices_->addEntity(EntitySetVector_[l]->getEntity(i));
                                break;
                            default:
                                break;
                        }
                    }
                    break;
                default:
                    for (UN i=0; i<EntitySetVector_[l]->getNumEntities(); i++) {
                        switch (EntitySetVector_[l]->getEntity(i)->getEntityFlag()) {
                            case DefaultFlag: // By default, an entity which belongs to more than 2 subdomains is an edge
                                EntitySetVector_[l]->getEntity(i)->resetEntityType(EdgeType);
                                Edges_->addEntity(EntitySetVector_[l]->getEntity(i));
                                break;
                            case StraightFlag:  // If an entity is straight, it is always a straight edge
                                EntitySetVector_[l]->getEntity(i)->resetEntityType(EdgeType);
                                StraightEdges_->addEntity(EntitySetVector_[l]->getEntity(i));
                                break;
                            case ShortFlag: // If an entity is a short, it is always a short edge
                                EntitySetVector_[l]->getEntity(i)->resetEntityType(EdgeType);
                                ShortEdges_->addEntity(EntitySetVector_[l]->getEntity(i));
                                break;
                            case NodeFlag: // If an entity is a node, it is always a vertex
                                EntitySetVector_[l]->getEntity(i)->resetEntityType(VertexType);
                                Vertices_->addEntity(EntitySetVector_[l]->getEntity(i));
                                break;
                            default:
                                break;
                        }
                    }
                    break;
            }

        }
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::buildEntityMaps(bool buildVerticesMap,
                                                  bool buildShortEdgesMap,
                                                  bool buildStraightEdgesMap,
                                                  bool buildEdgesMap,
                                                  bool buildFacesMap,
                                                  bool buildRootsMap,
                                                  bool buildLeafsMap)
    {
        FROSCH_DETAILTIMER_START_LEVELID(buildEntityMapsTime,"DDInterface::buildEntityMaps");
        //if (Verbose_ && Verbosity_==All) cout << "FROSch::DDInterface : Building global interface component maps" << endl;

        if (buildVerticesMap) Vertices_->buildEntityMap(NodesMap_);
        if (buildShortEdgesMap) ShortEdges_->buildEntityMap(NodesMap_);
        if (buildStraightEdgesMap) StraightEdges_->buildEntityMap(NodesMap_);
        if (buildEdgesMap) Edges_->buildEntityMap(NodesMap_);
        if (buildFacesMap) Faces_->buildEntityMap(NodesMap_);
        if (buildRootsMap) Roots_->buildEntityMap(NodesMap_);
        if (buildLeafsMap) Leafs_->buildEntityMap(NodesMap_);

        if (Verbosity_==All) {
            FROSCH_DETAILTIMER_START_LEVELID(printStatisticsTime,"print statistics");
            // Count entities
            GOVec globalVec(7);
            LOVec localVec(7);
            LOVec sumVec(7);
            SCVec avgVec(7);
            LOVec minVec(7);
            LOVec maxVec(7);
            if (buildVerticesMap) {
                globalVec[0] = Vertices_->getEntityMap()->getMaxAllGlobalIndex();
                if (NodesMap_->lib()==UseEpetra || Vertices_->getEntityMap()->getGlobalNumElements()>0) {
                    globalVec[0] += 1;
                }
                if (globalVec[0]<0) globalVec[0] = 0;
                localVec[0] = (LO) max((LO) Vertices_->getEntityMap()->getLocalNumElements(),(LO) 0);
                reduceAll(*this->MpiComm_,REDUCE_SUM,localVec[0],ptr(&sumVec[0]));
                avgVec[0] = max(sumVec[0]/double(MpiComm_->getSize()),0.0);
                reduceAll(*MpiComm_,REDUCE_MIN,localVec[0],ptr(&minVec[0]));
                reduceAll(*MpiComm_,REDUCE_MAX,localVec[0],ptr(&maxVec[0]));
            } else {
                globalVec[0] = -1;
                localVec[0] = -1;
                avgVec[0] = -1;
                minVec[0] = -1;
                maxVec[0] = -1;
                sumVec[0] = -1;
            }
            if (buildShortEdgesMap) {
                globalVec[1] = ShortEdges_->getEntityMap()->getMaxAllGlobalIndex();
                if (NodesMap_->lib()==UseEpetra || ShortEdges_->getEntityMap()->getGlobalNumElements()>0) {
                    globalVec[1] += 1;
                }
                if (globalVec[1]<0) globalVec[1] = 0;
                localVec[1] = (LO) max((LO) ShortEdges_->getEntityMap()->getLocalNumElements(),(LO) 0);
                reduceAll(*this->MpiComm_,REDUCE_SUM,localVec[1],ptr(&sumVec[1]));
                avgVec[1] = max(sumVec[1]/double(MpiComm_->getSize()),0.0);
                reduceAll(*MpiComm_,REDUCE_MIN,localVec[1],ptr(&minVec[1]));
                reduceAll(*MpiComm_,REDUCE_MAX,localVec[1],ptr(&maxVec[1]));
            } else {
                globalVec[1] = -1;
                localVec[1] = -1;
                avgVec[1] = -1;
                minVec[1] = -1;
                maxVec[1] = -1;
                sumVec[1] = -1;
            }
            if (buildStraightEdgesMap) {
                globalVec[2] = StraightEdges_->getEntityMap()->getMaxAllGlobalIndex();
                if (NodesMap_->lib()==UseEpetra || StraightEdges_->getEntityMap()->getGlobalNumElements()>0) {
                    globalVec[2] += 1;
                }
                if (globalVec[2]<0) globalVec[2] = 0;
                localVec[2] = (LO) max((LO) StraightEdges_->getEntityMap()->getLocalNumElements(),(LO) 0);
                reduceAll(*this->MpiComm_,REDUCE_SUM,localVec[2],ptr(&sumVec[2]));
                avgVec[2] = max(sumVec[2]/double(MpiComm_->getSize()),0.0);
                reduceAll(*MpiComm_,REDUCE_MIN,localVec[2],ptr(&minVec[2]));
                reduceAll(*MpiComm_,REDUCE_MAX,localVec[2],ptr(&maxVec[2]));
            } else {
                globalVec[2] = -1;
                localVec[2] = -1;
                avgVec[2] = -1;
                minVec[2] = -1;
                maxVec[2] = -1;
                sumVec[2] = -1;
            }
            if (buildEdgesMap) {
                globalVec[3] = (LO) Edges_->getEntityMap()->getMaxAllGlobalIndex();
                if (NodesMap_->lib()==UseEpetra || Edges_->getEntityMap()->getGlobalNumElements()>0) {
                    globalVec[3] += 1;
                }
                if (globalVec[3]<0) globalVec[3] = 0;
                localVec[3] = max((LO) Edges_->getEntityMap()->getLocalNumElements(),(LO) 0);
                reduceAll(*this->MpiComm_,REDUCE_SUM,localVec[3],ptr(&sumVec[3]));
                avgVec[3] = max(sumVec[3]/double(MpiComm_->getSize()),0.0);
                reduceAll(*MpiComm_,REDUCE_MIN,localVec[3],ptr(&minVec[3]));
                reduceAll(*MpiComm_,REDUCE_MAX,localVec[3],ptr(&maxVec[3]));
            } else {
                globalVec[3] = -1;
                localVec[3] = -1;
                avgVec[3] = -1;
                minVec[3] = -1;
                maxVec[3] = -1;
                sumVec[3] = -1;
            }
            if (buildFacesMap) {
                globalVec[4] = (LO) Faces_->getEntityMap()->getMaxAllGlobalIndex();
                if (NodesMap_->lib()==UseEpetra || Faces_->getEntityMap()->getGlobalNumElements()>0) {
                    globalVec[4] += 1;
                }
                if (globalVec[4]<0) globalVec[4] = 0;
                localVec[4] = max((LO) Faces_->getEntityMap()->getLocalNumElements(),(LO) 0);
                reduceAll(*this->MpiComm_,REDUCE_SUM,localVec[4],ptr(&sumVec[4]));
                avgVec[4] = max(sumVec[4]/double(MpiComm_->getSize()),0.0);
                reduceAll(*MpiComm_,REDUCE_MIN,localVec[4],ptr(&minVec[4]));
                reduceAll(*MpiComm_,REDUCE_MAX,localVec[4],ptr(&maxVec[4]));
            } else {
                globalVec[4] = -1;
                localVec[4] = -1;
                avgVec[4] = -1;
                minVec[4] = -1;
                maxVec[4] = -1;
                sumVec[4] = -1;
            }
            if (buildRootsMap) {
                globalVec[5] = Roots_->getEntityMap()->getMaxAllGlobalIndex();
                if (NodesMap_->lib()==UseEpetra || Roots_->getEntityMap()->getGlobalNumElements()>0) {
                    globalVec[5] += 1;
                }
                if (globalVec[5]<0) globalVec[5] = 0;
                localVec[5] = (LO) max((LO) Roots_->getEntityMap()->getLocalNumElements(),(LO) 0);
                reduceAll(*this->MpiComm_,REDUCE_SUM,localVec[5],ptr(&sumVec[5]));
                avgVec[5] = max(sumVec[5]/double(MpiComm_->getSize()),0.0);
                reduceAll(*MpiComm_,REDUCE_MIN,localVec[5],ptr(&minVec[5]));
                reduceAll(*MpiComm_,REDUCE_MAX,localVec[5],ptr(&maxVec[5]));
            } else {
                globalVec[5] = -1;
                localVec[5] = -1;
                avgVec[5] = -1;
                minVec[5] = -1;
                maxVec[5] = -1;
                sumVec[5] = -1;
            }
            if (buildLeafsMap) {
                globalVec[6] = Leafs_->getEntityMap()->getMaxAllGlobalIndex();
                if (NodesMap_->lib()==UseEpetra || Leafs_->getEntityMap()->getGlobalNumElements()>0) {
                    globalVec[6] += 1;
                }
                if (globalVec[6]<0) globalVec[6] = 0;
                localVec[6] = (LO) max((LO) Leafs_->getEntityMap()->getLocalNumElements(),(LO) 0);
                reduceAll(*this->MpiComm_,REDUCE_SUM,localVec[6],ptr(&sumVec[6]));
                avgVec[6] = max(sumVec[6]/double(MpiComm_->getSize()),0.0);
                reduceAll(*MpiComm_,REDUCE_MIN,localVec[6],ptr(&minVec[6]));
                reduceAll(*MpiComm_,REDUCE_MAX,localVec[6],ptr(&maxVec[6]));
            } else {
                globalVec[6] = -1;
                localVec[6] = -1;
                avgVec[6] = -1;
                minVec[6] = -1;
                maxVec[6] = -1;
                sumVec[6] = -1;
            }

            for (UN i=0; i<globalVec.size(); i++) {
                if (globalVec[i]<0) {
                    globalVec[i] = -1;
                }
            }

            NumEntity_ = globalVec;

            if (Verbose_) {
                cout
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << setw(89) << "-----------------------------------------------------------------------------------------"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| "
                << left << setw(74) << "> Interface Statistics " << right << setw(8) << "(Level " << setw(2) << LevelID_ << ")"
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << setw(89) << "========================================================================================="
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| " << left << setw(41) << "Interface communication strategy" << right
                << " | " << setw(41) << CommStrategy_
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << setw(89) << "-----------------------------------------------------------------------------------------"
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
                << "| " << left << setw(20) << "Vertices" << right
                << " | "; globalVec[0]<0 ? cout << setw(10) << " " : cout << setw(10) << globalVec[0]; cout
                << " | "; avgVec[0]<0 ? cout << setw(10) << " " : cout << setw(10) << setprecision(5) << avgVec[0]; cout
                << " | "; minVec[0]<0 ? cout << setw(10) << " " : cout << setw(10) << minVec[0]; cout
                << " | "; maxVec[0]<0 ? cout << setw(10) << " " : cout << setw(10) << maxVec[0]; cout
                << " | "; sumVec[0]<0 ? cout << setw(10) << " " : cout << setw(10) << sumVec[0]; cout
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| " << left << setw(20) << "Short edges" << right
                << " | "; globalVec[1]<0 ? cout << setw(10) << " " : cout << setw(10) << globalVec[1]; cout
                << " | "; avgVec[1]<0 ? cout << setw(10) << " " : cout << setw(10) << setprecision(5) << avgVec[1]; cout
                << " | "; minVec[1]<0 ? cout << setw(10) << " " : cout << setw(10) << minVec[1]; cout
                << " | "; maxVec[1]<0 ? cout << setw(10) << " " : cout << setw(10) << maxVec[1]; cout
                << " | "; sumVec[1]<0 ? cout << setw(10) << " " : cout << setw(10) << sumVec[1]; cout
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| " << left << setw(20) << "Straight edges" << right
                << " | "; globalVec[2]<0 ? cout << setw(10) << " " : cout << setw(10) << globalVec[2]; cout
                << " | "; avgVec[2]<0 ? cout << setw(10) << " " : cout << setw(10) << setprecision(5) << avgVec[2]; cout
                << " | "; minVec[2]<0 ? cout << setw(10) << " " : cout << setw(10) << minVec[2]; cout
                << " | "; maxVec[2]<0 ? cout << setw(10) << " " : cout << setw(10) << maxVec[2]; cout
                << " | "; sumVec[2]<0 ? cout << setw(10) << " " : cout << setw(10) << sumVec[2]; cout
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| " << left << setw(20) << "Edges" << right
                << " | "; globalVec[3]<0 ? cout << setw(10) << " " : cout << setw(10) << globalVec[3]; cout
                << " | "; avgVec[3]<0 ? cout << setw(10) << " " : cout << setw(10) << setprecision(5) << avgVec[3]; cout
                << " | "; minVec[3]<0 ? cout << setw(10) << " " : cout << setw(10) << minVec[3]; cout
                << " | "; maxVec[3]<0 ? cout << setw(10) << " " : cout << setw(10) << maxVec[3]; cout
                << " | "; sumVec[3]<0 ? cout << setw(10) << " " : cout << setw(10) << sumVec[3]; cout
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| " << left << setw(20) << "Faces" << right
                << " | "; globalVec[4]<0 ? cout << setw(10) << " " : cout << setw(10) << globalVec[4]; cout
                << " | "; avgVec[4]<0 ? cout << setw(10) << " " : cout << setw(10) << setprecision(5) << avgVec[4]; cout
                << " | "; minVec[4]<0 ? cout << setw(10) << " " : cout << setw(10) << minVec[4]; cout
                << " | "; maxVec[4]<0 ? cout << setw(10) << " " : cout << setw(10) << maxVec[4]; cout
                << " | "; sumVec[4]<0 ? cout << setw(10) << " " : cout << setw(10) << sumVec[4]; cout
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| " << left << setw(20) << "Roots" << right
                << " | "; globalVec[5]<0 ? cout << setw(10) << " " : cout << setw(10) << globalVec[5]; cout
                << " | "; avgVec[5]<0 ? cout << setw(10) << " " : cout << setw(10) << setprecision(5) << avgVec[5]; cout
                << " | "; minVec[5]<0 ? cout << setw(10) << " " : cout << setw(10) << minVec[5]; cout
                << " | "; maxVec[5]<0 ? cout << setw(10) << " " : cout << setw(10) << maxVec[5]; cout
                << " | "; sumVec[5]<0 ? cout << setw(10) << " " : cout << setw(10) << sumVec[5]; cout
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| " << left << setw(20) << "Leafs" << right
                << " | "; globalVec[6]<0 ? cout << setw(10) << " " : cout << setw(10) << globalVec[6]; cout
                << " | "; avgVec[6]<0 ? cout << setw(10) << " " : cout << setw(10) << setprecision(5) << avgVec[6]; cout
                << " | "; minVec[6]<0 ? cout << setw(10) << " " : cout << setw(10) << minVec[6]; cout
                << " | "; maxVec[6]<0 ? cout << setw(10) << " " : cout << setw(10) << maxVec[6]; cout
                << " | "; sumVec[6]<0 ? cout << setw(10) << " " : cout << setw(10) << sumVec[6]; cout
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << setw(89) << "-----------------------------------------------------------------------------------------"
                << endl;
            }
        }

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::buildEntityHierarchy()
    {
        FROSCH_DETAILTIMER_START_LEVELID(buildEntityHierarchyTime,"DDInterface::buildEntityHierarchy");
        //if (Verbose_ && Verbosity_==All) cout << "FROSch::DDInterface : Building hierarchy of interface components" << endl;

        // Build hierarchy
        for (UN i=0; i<EntitySetVector_.size(); i++) {
            for (UN j=i+1; j<EntitySetVector_.size(); j++) {
                EntitySetVector_[i]->findAncestorsInSet(EntitySetVector_[j]);
            }
        }

        // Find roots
        for (UN i=0; i<EntitySetVector_.size(); i++) {
            EntitySetPtr tmpRoots = EntitySetVector_[i]->findRoots();
            Roots_->addEntitySet(tmpRoots);
        }
        Roots_->sortUnique();
        Roots_->setRootID();

        // Find Leafs
        for (UN i=0; i<EntitySetVector_.size(); i++) {
            EntitySetPtr tmpLeafs = EntitySetVector_[i]->findLeafs();
            Leafs_->addEntitySet(tmpLeafs);
        }
        Leafs_->sortUnique();
        Leafs_->setLeafID();
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::computeDistancesToRoots(UN dimension,
                                                          ConstXMultiVectorPtr &nodeList,
                                                          DistanceFunction distanceFunction)
    {
        FROSCH_DETAILTIMER_START_LEVELID(computeDistancesToRootsTime,"DDInterface::computeDistancesToRoots");
        //if (Verbose_ && Verbosity_==All) cout << "FROSch::DDInterface : Computing distances to the coarse nodes" << endl;

        for (UN i=0; i<EntitySetVector_.size(); i++) {
            EntitySetVector_[i]->computeDistancesToRoots(dimension,nodeList,distanceFunction);
        }
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::identifyConnectivityEntities(UNVecPtr multiplicities,
                                                               EntityFlagVecPtr flags)
    {
        FROSCH_DETAILTIMER_START_LEVELID(identifyConnectivityEntitiesTime,"DDInterface::identifyConnectivityEntities");
        //if (Verbose_ && Verbosity_==All) cout << "FROSch::DDInterface : Preparing subdomain graph" << endl;

        if (multiplicities.is_null()) {
            multiplicities = UNVecPtr(1,2);
        }
        if (flags.is_null()) {
            flags = EntityFlagVecPtr(4);
            flags[0] = DefaultFlag;
            flags[1] = StraightFlag;
            flags[2] = ShortFlag;
            flags[3] = NodeFlag;
        }

        for (UN j=0; j<multiplicities.size(); j++) {
            for (UN i=0; i<EntitySetVector_[multiplicities[j]]->getNumEntities(); i++) {
                if (binary_search(flags.begin(),flags.end(),EntitySetVector_[multiplicities[j]]->getEntity(i)->getEntityFlag())) {
                    ConnectivityEntities_->addEntity(EntitySetVector_[multiplicities[j]]->getEntity(i));
                }
            }
        }
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    typename DDInterface<SC,LO,GO,NO>::UN DDInterface<SC,LO,GO,NO>::getDimension() const
    {
        return Dimension_;
    }

    template <class SC,class LO,class GO,class NO>
    typename DDInterface<SC,LO,GO,NO>::UN DDInterface<SC,LO,GO,NO>::getDofsPerNode() const
    {
        return DofsPerNode_;
    }

    template <class SC,class LO,class GO,class NO>
    LO DDInterface<SC,LO,GO,NO>::getNumMyNodes() const
    {
        return NumMyNodes_;
    }

    template <class SC,class LO,class GO,class NO>
    typename DDInterface<SC,LO,GO,NO>::EntitySetConstPtr & DDInterface<SC,LO,GO,NO>::getVertices() const
    {
        return Vertices_;
    }

    template <class SC,class LO,class GO,class NO>
    typename DDInterface<SC,LO,GO,NO>::EntitySetConstPtr & DDInterface<SC,LO,GO,NO>::getShortEdges() const
    {
        return ShortEdges_;
    }

    template <class SC,class LO,class GO,class NO>
    typename DDInterface<SC,LO,GO,NO>::EntitySetConstPtr & DDInterface<SC,LO,GO,NO>::getStraightEdges() const
    {
        return StraightEdges_;
    }

    template <class SC,class LO,class GO,class NO>
    typename DDInterface<SC,LO,GO,NO>::EntitySetConstPtr & DDInterface<SC,LO,GO,NO>::getEdges() const
    {
        return Edges_;
    }

    template <class SC,class LO,class GO,class NO>
    typename DDInterface<SC,LO,GO,NO>::EntitySetConstPtr & DDInterface<SC,LO,GO,NO>::getFaces() const
    {
        return Faces_;
    }

    template <class SC,class LO,class GO,class NO>
    typename DDInterface<SC,LO,GO,NO>::EntitySetConstPtr & DDInterface<SC,LO,GO,NO>::getInterface() const
    {
        return Interface_;
    }

    template <class SC,class LO,class GO,class NO>
    typename DDInterface<SC,LO,GO,NO>::EntitySetConstPtr & DDInterface<SC,LO,GO,NO>::getInterior() const
    {
        return Interior_;
    }

    template <class SC,class LO,class GO,class NO>
    typename DDInterface<SC,LO,GO,NO>::EntitySetConstPtr & DDInterface<SC,LO,GO,NO>::getRoots() const
    {
        return Roots_;
    }

    template <class SC,class LO,class GO,class NO>
    typename DDInterface<SC,LO,GO,NO>::EntitySetConstPtr & DDInterface<SC,LO,GO,NO>::getLeafs() const
    {
        return Leafs_;
    }

    template <class SC,class LO,class GO,class NO>
    typename DDInterface<SC,LO,GO,NO>::EntitySetPtrConstVecPtr & DDInterface<SC,LO,GO,NO>::getEntitySetVector() const
    {
        return EntitySetVector_;
    }

    template <class SC,class LO,class GO,class NO>
    typename DDInterface<SC,LO,GO,NO>::GOVec DDInterface<SC,LO,GO,NO>::getNumEnt() const
    {
        return NumEntity_;
    }

    template <class SC,class LO,class GO,class NO>
    typename DDInterface<SC,LO,GO,NO>::EntitySetConstPtr & DDInterface<SC,LO,GO,NO>::getConnectivityEntities() const
    {
        return ConnectivityEntities_;
    }

    template <class SC,class LO,class GO,class NO>
    typename DDInterface<SC,LO,GO,NO>::ConstXMapPtr DDInterface<SC,LO,GO,NO>::getNodesMap() const
    {
        return NodesMap_.getConst();
    }

    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::communicateLocalComponents(IntVecVecPtr &componentsSubdomains,
                                                             IntVecVec &componentsSubdomainsUnique)
    {
        FROSCH_DETAILTIMER_START_LEVELID(communicateLocalComponentsTime,"DDInterface::communicateLocalComponents");
        //if (Verbose_ && Verbosity_==All) cout << "FROSch::DDInterface : Communicating nodes" << endl;

        if (NodesMap_->lib() == UseEpetra && CommStrategy_ == CreateOneToOneMap) {
            FROSCH_WARNING("FROSch::DDInterface",Verbose_,"CreateOneToOneMap communication strategy does not work for Epetra => Switching to CommCrsGraph.");
            CommStrategy_ = CommCrsGraph;
        }

        // Different communication strategies
        const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();
        switch (CommStrategy_) {
            case CommCrsMatrix:
                {
                    UniqueNodesMap_ = BuildUniqueMap<LO,GO,NO>(NodesMap_);
                    RCP<Matrix<SC,LO,GO,NO> > commMat = MatrixFactory<SC,LO,GO,NO>::Build(NodesMap_,10);
                    RCP<Matrix<SC,LO,GO,NO> > commMatTmp = MatrixFactory<SC,LO,GO,NO>::Build(UniqueNodesMap_,10);
                    XExportPtr commExporter = ExportFactory<LO,GO,NO>::Build(NodesMap_,UniqueNodesMap_);

                    Array<SC> one(1,ScalarTraits<SC>::one());
                    Array<GO> myPID(1,MpiComm_->getRank());
                    for (int i=0; i<NumMyNodes_; i++) {
                        commMat->insertGlobalValues(NodesMap_->getGlobalElement(i),myPID(),one());
                    }
                    XMapPtr domainMap = MapFactory<LO,GO,NO>::Build(NodesMap_->lib(),INVALID,myPID(),0,NodesMap_->getComm());

                    commMat->fillComplete(domainMap,NodesMap_);
                    commMatTmp->doExport(*commMat,*commExporter,INSERT);
                    commMatTmp->fillComplete(domainMap,UniqueNodesMap_);
                    commMat = MatrixFactory<SC,LO,GO,NO>::Build(NodesMap_,LO(0));
                    commMat->doImport(*commMatTmp,*commExporter,INSERT);

                    componentsSubdomains = IntVecVecPtr(NumMyNodes_);

                    ArrayView<const GO> indices;
                    ArrayView<const SC> values;
                    for (LO i=0; i<NumMyNodes_; i++) {
                        commMat->getGlobalRowView(NodesMap_->getGlobalElement(i),indices,values);
                        componentsSubdomains[i].resize(indices.size());
                        for (LO j=0; j<indices.size(); j++) {
                            componentsSubdomains[i][j] = as<int>(indices[j]);
                        }
                    }
                }
                break;

            case CommCrsGraph:
                {
                    UniqueNodesMap_ = BuildUniqueMap<LO,GO,NO>(NodesMap_);

                    XCrsGraphPtr commGraph = CrsGraphFactory<LO,GO,NO>::Build(NodesMap_,10); // AH 08/07/2019: Can we put 1 instead of 10 here?
                    XCrsGraphPtr commGraphTmp = CrsGraphFactory<LO,GO,NO>::Build(UniqueNodesMap_,10); // We assume that any node is part of no more than 10 subdomains
                    XExportPtr commExporter = ExportFactory<LO,GO,NO>::Build(NodesMap_,UniqueNodesMap_);

                    Array<GO> myPID(1,MpiComm_->getRank());
                    for (int i=0; i<NumMyNodes_; i++) {
                        commGraph->insertGlobalIndices(NodesMap_->getGlobalElement(i),myPID());
                    }
                    XMapPtr domainMap = MapFactory<LO,GO,NO>::Build(NodesMap_->lib(),INVALID,myPID(),0,NodesMap_->getComm());

                    commGraph->fillComplete(domainMap,NodesMap_); // AH 08/07/2019: Can we remove some fillComplete?
                    commGraphTmp->doExport(*commGraph,*commExporter,INSERT);
                    commGraphTmp->fillComplete(domainMap,UniqueNodesMap_);
                    commGraph = CrsGraphFactory<LO,GO,NO>::Build(NodesMap_);
                    commGraph->doImport(*commGraphTmp,*commExporter,INSERT);

                    componentsSubdomains = IntVecVecPtr(NumMyNodes_);

                    ArrayView<const GO> indices;
                    for (LO i=0; i<NumMyNodes_; i++) {
                        commGraph->getGlobalRowView(NodesMap_->getGlobalElement(i),indices);
                        componentsSubdomains[i].resize(indices.size());
                        for (LO j=0; j<indices.size(); j++) {
                            componentsSubdomains[i][j] = as<int>(indices[j]);
                        }
                    }
                }
                break;

            case CreateOneToOneMap:
                {
                    RCP<LowerPIDTieBreak<LO,GO,NO> > lowerPIDTieBreak(new LowerPIDTieBreak<LO,GO,NO>(MpiComm_,NodesMap_,Dimension_,LevelID_));
                    UniqueNodesMap_ = BuildUniqueMap<LO,GO,NO>(NodesMap_,true,lowerPIDTieBreak);
                    lowerPIDTieBreak->sendDataToOriginalMap();
                    componentsSubdomains = lowerPIDTieBreak->getComponents();
                }
                break;

            default:
                FROSCH_ASSERT(false,"FROSch::DDInterface: Specify a valid communication strategy.");
        }

        componentsSubdomainsUnique = IntVecVec(NumMyNodes_);
        for (LO i=0; i<NumMyNodes_; i++) {
            sortunique(componentsSubdomains[i]);
            if (componentsSubdomains[i].size() == 0) componentsSubdomains[i].push_back(MpiComm_->getRank()); // For Tpetra this is empty if the repeatedMap is already unique. In this case, we have to add the local rank. Otherwise, we obtain nodes with multiplicity 0.
            componentsSubdomainsUnique[i] = componentsSubdomains[i];
//            if (MpiComm_->getRank() == 0) cout << MpiComm_->getRank() << ": " << i << " " << componentsSubdomains[i] << endl;
        }
        sortunique(componentsSubdomainsUnique);

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::identifyLocalComponents(IntVecVecPtr &componentsSubdomains,
                                                          IntVecVec &componentsSubdomainsUnique)
    {
        FROSCH_DETAILTIMER_START_LEVELID(identifyLocalComponentsTime,"DDInterface::identifyLocalComponents");
        //if (Verbose_ && Verbosity_==All) cout << "FROSch::DDInterface : Classifying interface components based on equivalence classes" << endl;

        // Hier herausfinden, ob Ecke, Kante oder Fläche
        UNVecPtr componentsMultiplicity(componentsSubdomainsUnique.size());
        IntVecVecPtr components(componentsSubdomainsUnique.size());
        IntVecVecPtr componentsGamma(componentsSubdomainsUnique.size());
        UN maxMultiplicity = 0;
        for (UN i=0; i<componentsSubdomainsUnique.size(); i++) {
            componentsMultiplicity[i] = componentsSubdomainsUnique[i].size();
            maxMultiplicity = max(maxMultiplicity,componentsMultiplicity[i]);
        }
        EntitySetVector_ = EntitySetPtrVecPtr(maxMultiplicity+1);
        for (UN i=0; i<maxMultiplicity+1; i++) {
            EntitySetVector_[i].reset(new EntitySet<SC,LO,GO,NO>(DefaultType));
        }

        typename IntVecVec::iterator classIterator;
        LOVecPtr localComponentIndices(NumMyNodes_);
        for (int i=0; i<NumMyNodes_; i++) {
            classIterator = lower_bound(componentsSubdomainsUnique.begin(),componentsSubdomainsUnique.end(),componentsSubdomains[i]);
            localComponentIndices[i] = classIterator - componentsSubdomainsUnique.begin();
        }

        LO tmp1 = 0; // The interface and interior have multiplicity 0 in our construction
        int *tmp2 = NULL;
        RCP<InterfaceEntity<SC,LO,GO,NO> > interior(new InterfaceEntity<SC,LO,GO,NO>(InteriorType,DofsPerNode_,tmp1,tmp2));
        RCP<InterfaceEntity<SC,LO,GO,NO> > interface(new InterfaceEntity<SC,LO,GO,NO>(InterfaceType,DofsPerNode_,tmp1,tmp2));
        for (LO i=0; i<NumMyNodes_; i++) {
            if (componentsMultiplicity[localComponentIndices[i]] == 1) {
                LO nodeIDI = interior->getNumNodes();
                LO nodeIDLocal = i;
                GO nodeIDGlobal = NodesMap_->getGlobalElement(nodeIDLocal);
                LOVecPtr dofsI(DofsPerNode_);
                LOVecPtr dofsLocal(DofsPerNode_);
                GOVecPtr dofsGlobal(DofsPerNode_);
                for (UN k=0; k<DofsPerNode_; k++) {
                    dofsI[k] = DofsPerNode_*nodeIDI+k;
                    dofsLocal[k] = DofsPerNode_*nodeIDLocal+k;
                    dofsGlobal[k] = DofsPerNode_*nodeIDGlobal+k;
                }
                interior->addNode(nodeIDI,nodeIDLocal,nodeIDGlobal,DofsPerNode_,dofsI,dofsLocal,dofsGlobal);
            } else {
                FROSCH_ASSERT(componentsMultiplicity[localComponentIndices[i]]>1,"FROSch::DDInterface: There cannot be any nodes with multiplicity 0.");
                LO nodeIDGamma = interface->getNumNodes();
                LO nodeIDLocal = i;
                GO nodeIDGlobal = NodesMap_->getGlobalElement(nodeIDLocal);
                LOVecPtr dofsGamma(DofsPerNode_);
                LOVecPtr dofsLocal(DofsPerNode_);
                GOVecPtr dofsGlobal(DofsPerNode_);
                for (UN k=0; k<DofsPerNode_; k++) {
                    dofsGamma[k] = DofsPerNode_*nodeIDGamma+k;
                    dofsLocal[k] = DofsPerNode_*nodeIDLocal+k;
                    dofsGlobal[k] = DofsPerNode_*nodeIDGlobal+k;
                }
                interface->addNode(nodeIDGamma,nodeIDLocal,nodeIDGlobal,DofsPerNode_,dofsGamma,dofsLocal,dofsGlobal);

                components[localComponentIndices[i]].push_back(i);
                componentsGamma[localComponentIndices[i]].push_back(interface->getNumNodes()-1);
            }
        }
        Interior_->addEntity(interior);
        Interface_->addEntity(interface);

        for (UN i=0; i<componentsSubdomainsUnique.size(); i++) {
            FROSCH_ASSERT(componentsMultiplicity[i]>0,"FROSch::DDInterface: There cannot be any component with multiplicity 0.");
            RCP<InterfaceEntity<SC,LO,GO,NO> > tmpEntity(new InterfaceEntity<SC,LO,GO,NO>(VertexType,DofsPerNode_,componentsMultiplicity[i],&(componentsSubdomainsUnique[i][0])));
            LO nodeIDGamma;
            LO nodeIDLocal;
            GO nodeIDGlobal;
            LOVecPtr dofsGamma(DofsPerNode_);
            LOVecPtr dofsLocal(DofsPerNode_);
            GOVecPtr dofsGlobal(DofsPerNode_);

            sortunique(components[i]);

            //
            for (UN j=0; j<components[i].size(); j++) {
                nodeIDGamma = componentsGamma[i][j];
                nodeIDLocal = components[i][j];
                nodeIDGlobal = NodesMap_->getGlobalElement(nodeIDLocal);
                for (UN k=0; k<DofsPerNode_; k++) {
                    dofsGamma[k] = DofsPerNode_*nodeIDGamma+k;
                    dofsLocal[k] = DofsPerNode_*nodeIDLocal+k;
                    dofsGlobal[k] = DofsPerNode_*nodeIDGlobal+k;
                }

                tmpEntity->addNode(nodeIDGamma,nodeIDLocal,nodeIDGlobal,DofsPerNode_,dofsGamma,dofsLocal,dofsGlobal);
            }
            tmpEntity->resetEntityType(DefaultType);
            EntitySetVector_[componentsMultiplicity[i]]->addEntity(tmpEntity);
        }

        // Remove the empty entity stemming from the interior nodes
        removeEmptyEntities();

        // We need to set the unique ID; otherwise, we cannot sort entities
        for (UN i=0; i<EntitySetVector_.size(); i++) {
            EntitySetVector_[i]->setUniqueIDToFirstGlobalNodeID();
        }
        return 0;
    }

}

#endif
