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

#ifndef _FROSCH_DDINTERFACE_DEF_HPP
#define _FROSCH_DDINTERFACE_DEF_HPP

#include <FROSch_DDInterface_decl.hpp>


namespace FROSch {

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
    NumMyNodes_ (localToGlobalMap->getNodeNumElements()),
    Vertices_ (new EntitySet<SC,LO,GO,NO>(VertexType)),
    ShortEdges_ (new EntitySet<SC,LO,GO,NO>(EdgeType)),
    StraightEdges_ (new EntitySet<SC,LO,GO,NO>(EdgeType)),
    Edges_ (new EntitySet<SC,LO,GO,NO>(EdgeType)),
    Faces_ (new EntitySet<SC,LO,GO,NO>(FaceType)),
    Interface_ (new EntitySet<SC,LO,GO,NO>(InterfaceType)),
    Interior_ (new EntitySet<SC,LO,GO,NO>(InteriorType)),
    Roots_ (new EntitySet<SC,LO,GO,NO>(DefaultType)),
    Leafs_ (new EntitySet<SC,LO,GO,NO>(DefaultType)),
    ConnectivityEntities_ (new EntitySet<SC,LO,GO,NO>(DefaultType)),
    EntitySetVector_ (),
    NodesMap_ (localToGlobalMap),
    UniqueNodesMap_ (),
    Verbose_ (MpiComm_->getRank()==0),
    Verbosity_ (verbosity),
    LevelID_ (levelID)
    {
        FROSCH_TIMER_START_LEVELID(dDInterfaceTime,"DDInterface::DDInterface");
        FROSCH_ASSERT(((Dimension_==2)||(Dimension_==3)),"FROSch::DDInterface : ERROR: Only dimension 2 and 3 are available");

        //if (Verbose_ && Verbosity_==All) std::cout << "FROSch::DDInterface" << std::endl;

        IntVecVecPtr componentsSubdomains;
        IntVecVec componentsSubdomainsUnique;

        communicateLocalComponents(componentsSubdomains,componentsSubdomainsUnique,commStrategy);

        identifyLocalComponents(componentsSubdomains,componentsSubdomainsUnique);
    }

    template <class SC,class LO,class GO,class NO>
    DDInterface<SC,LO,GO,NO>::~DDInterface()
    {

    } // Do we need sth here?

    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::resetGlobalDofs(ConstXMapPtrVecPtr dofsMaps)
    {
        FROSCH_TIMER_START_LEVELID(resetGlobalDofsTime,"DDInterface::resetGlobalDofs");
        //if (Verbose_ && Verbosity_==All) std::cout << "FROSch::DDInterface : Resetting Global IDs" << std::endl;

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
        FROSCH_TIMER_START_LEVELID(removeDirichletNodesTime,"DDInterface::removeDirichletNodes");
        //if (Verbose_ && Verbosity_==All) std::cout << "FROSch::DDInterface : Removing Dirichlet Nodes from the domain decomposition interface" << std::endl;

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
        FROSCH_TIMER_START_LEVELID(divideUnconnectedEntitiesTime,"DDInterface::divideUnconnectedEntities");
        //if (Verbose_ && Verbosity_==All) std::cout << "FROSch::DDInterface : Decomposing unconnected interface components" << std::endl;

        GOVecPtr indicesGammaDofs(DofsPerNode_*Interface_->getEntity(0)->getNumNodes());
        for (UN k=0; k<DofsPerNode_; k++) {
            for (UN i=0; i<Interface_->getEntity(0)->getNumNodes(); i++) {
                indicesGammaDofs[Interface_->getEntity(0)->getGammaDofID(i,k)] = Interface_->getEntity(0)->getGlobalDofID(i,k);
            }
        }
        XMapPtr map = MapFactory<LO,GO,NO>::Build(matrix->getRowMap()->lib(),-1,indicesGammaDofs(),0,MpiComm_);
        matrix = FROSch::ExtractLocalSubdomainMatrix(matrix.getConst(),map.getConst(),ScalarTraits<SC>::one());

        // Operate on hierarchy
        for (UN i=0; i<EntitySetVector_.size(); i++) {
            EntitySetVector_[i]->divideUnconnectedEntities(matrix,MpiComm_->getRank());
        }

        /*
        LO numSeparateEdges = Edges_->divideUnconnectedEntities(matrix,MpiComm_->getRank());
        LO numSeparateFaces = Faces_->divideUnconnectedEntities(matrix,MpiComm_->getRank());

        if (Verbose_ && Verbosity_==All) {
            std::cout << "\n\
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
        FROSCH_TIMER_START_LEVELID(flagEntitiesTime,"DDInterface::flagEntities");
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
        FROSCH_TIMER_START_LEVELID(removeEmptyEntitiesTime,"DDInterface::removeEmptyEntities");
        //if (Verbose_ && Verbosity_==All) std::cout << "FROSch::DDInterface : Removing empty interface components" << std::endl;

        for (UN l=0; l<EntitySetVector_.size(); l++) {
            EntitySetVector_[l]->removeEmptyEntities();
        }
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::sortVerticesEdgesFaces(ConstXMultiVectorPtr nodeList)
    {
        FROSCH_TIMER_START_LEVELID(sortVerticesEdgesFacesTime,"DDInterface::sortVerticesEdgesFaces");
        //if (Verbose_ && Verbosity_==All) std::cout << "FROSch::DDInterface : Sorting interface components" << std::endl;

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
                    FROSCH_ASSERT(EntitySetVector_[l]->getNumEntities()==0,"FROSch::DDInterface : ERROR: This case is impossible.");
                    break;
                case 1:
                    FROSCH_ASSERT(EntitySetVector_[l]->getNumEntities()==0,"FROSch::DDInterface : ERROR: In this case, the entity is interior to the subdomain.");
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
        FROSCH_TIMER_START_LEVELID(buildEntityMapsTime,"DDInterface::buildEntityMaps");
        //if (Verbose_ && Verbosity_==All) std::cout << "FROSch::DDInterface : Building global interface component maps" << std::endl;

        if (buildVerticesMap) Vertices_->buildEntityMap(NodesMap_);
        if (buildShortEdgesMap) ShortEdges_->buildEntityMap(NodesMap_);
        if (buildStraightEdgesMap) StraightEdges_->buildEntityMap(NodesMap_);
        if (buildEdgesMap) Edges_->buildEntityMap(NodesMap_);
        if (buildFacesMap) Faces_->buildEntityMap(NodesMap_);
        if (buildRootsMap) Roots_->buildEntityMap(NodesMap_);
        if (buildLeafsMap) Leafs_->buildEntityMap(NodesMap_);

        if (Verbosity_==All) {
            // Count entities
            GOVec global(7);
            LOVec local(7);
            LOVec sum(7);
            SCVec avg(7);
            LOVec min(7);
            LOVec max(7);
            if (buildVerticesMap) {
                global[0] = Vertices_->getEntityMap()->getMaxAllGlobalIndex();
                if (NodesMap_->lib()==UseEpetra || Vertices_->getEntityMap()->getGlobalNumElements()>0) {
                    global[0] += 1;
                }
                if (global[0]<0) global[0] = 0;
                local[0] = (LO) std::max((LO) Vertices_->getEntityMap()->getNodeNumElements(),(LO) 0);
                reduceAll(*this->MpiComm_,REDUCE_SUM,local[0],ptr(&sum[0]));
                avg[0] = std::max(sum[0]/double(MpiComm_->getSize()),0.0);
                reduceAll(*MpiComm_,REDUCE_MIN,local[0],ptr(&min[0]));
                reduceAll(*MpiComm_,REDUCE_MAX,local[0],ptr(&max[0]));
            } else {
                global[0] = -1;
                local[0] = -1;
                avg[0] = -1;
                min[0] = -1;
                max[0] = -1;
            }
            if (buildShortEdgesMap) {
                global[1] = ShortEdges_->getEntityMap()->getMaxAllGlobalIndex();
                if (NodesMap_->lib()==UseEpetra || ShortEdges_->getEntityMap()->getGlobalNumElements()>0) {
                    global[1] += 1;
                }
                if (global[1]<0) global[1] = 0;
                local[1] = (LO) std::max((LO) ShortEdges_->getEntityMap()->getNodeNumElements(),(LO) 0);
                reduceAll(*this->MpiComm_,REDUCE_SUM,local[1],ptr(&sum[1]));
                avg[1] = std::max(sum[1]/double(MpiComm_->getSize()),0.0);
                reduceAll(*MpiComm_,REDUCE_MIN,local[1],ptr(&min[1]));
                reduceAll(*MpiComm_,REDUCE_MAX,local[1],ptr(&max[1]));
            } else {
                global[1] = -1;
                local[1] = -1;
                avg[1] = -1;
                min[1] = -1;
                max[1] = -1;
            }
            if (buildStraightEdgesMap) {
                global[2] = StraightEdges_->getEntityMap()->getMaxAllGlobalIndex();
                if (NodesMap_->lib()==UseEpetra || StraightEdges_->getEntityMap()->getGlobalNumElements()>0) {
                    global[2] += 1;
                }
                if (global[2]<0) global[2] = 0;
                local[2] = (LO) std::max((LO) StraightEdges_->getEntityMap()->getNodeNumElements(),(LO) 0);
                reduceAll(*this->MpiComm_,REDUCE_SUM,local[2],ptr(&sum[2]));
                avg[2] = std::max(sum[2]/double(MpiComm_->getSize()),0.0);
                reduceAll(*MpiComm_,REDUCE_MIN,local[2],ptr(&min[2]));
                reduceAll(*MpiComm_,REDUCE_MAX,local[2],ptr(&max[2]));
            } else {
                global[2] = -1;
                local[2] = -1;
                avg[2] = -1;
                min[2] = -1;
                max[2] = -1;
            }
            if (buildEdgesMap) {
                global[3] = (LO) Edges_->getEntityMap()->getMaxAllGlobalIndex();
                if (NodesMap_->lib()==UseEpetra || Edges_->getEntityMap()->getGlobalNumElements()>0) {
                    global[3] += 1;
                }
                if (global[3]<0) global[3] = 0;
                local[3] = std::max((LO) Edges_->getEntityMap()->getNodeNumElements(),(LO) 0);
                reduceAll(*this->MpiComm_,REDUCE_SUM,local[3],ptr(&sum[3]));
                avg[3] = std::max(sum[3]/double(MpiComm_->getSize()),0.0);
                reduceAll(*MpiComm_,REDUCE_MIN,local[3],ptr(&min[3]));
                reduceAll(*MpiComm_,REDUCE_MAX,local[3],ptr(&max[3]));
            } else {
                global[3] = -1;
                local[3] = -1;
                avg[3] = -1;
                min[3] = -1;
                max[3] = -1;
            }
            if (buildFacesMap) {
                global[4] = (LO) Faces_->getEntityMap()->getMaxAllGlobalIndex();
                if (NodesMap_->lib()==UseEpetra || Faces_->getEntityMap()->getGlobalNumElements()>0) {
                    global[4] += 1;
                }
                if (global[4]<0) global[4] = 0;
                local[4] = std::max((LO) Faces_->getEntityMap()->getNodeNumElements(),(LO) 0);
                reduceAll(*this->MpiComm_,REDUCE_SUM,local[4],ptr(&sum[4]));
                avg[4] = std::max(sum[4]/double(MpiComm_->getSize()),0.0);
                reduceAll(*MpiComm_,REDUCE_MIN,local[4],ptr(&min[4]));
                reduceAll(*MpiComm_,REDUCE_MAX,local[4],ptr(&max[4]));
            } else {
                global[4] = -1;
                local[4] = -1;
                avg[4] = -1;
                min[4] = -1;
                max[4] = -1;
            }
            if (buildRootsMap) {
                global[5] = Roots_->getEntityMap()->getMaxAllGlobalIndex();
                if (NodesMap_->lib()==UseEpetra || Roots_->getEntityMap()->getGlobalNumElements()>0) {
                    global[5] += 1;
                }
                if (global[5]<0) global[5] = 0;
                local[5] = (LO) std::max((LO) Roots_->getEntityMap()->getNodeNumElements(),(LO) 0);
                reduceAll(*this->MpiComm_,REDUCE_SUM,local[5],ptr(&sum[5]));
                avg[5] = std::max(sum[5]/double(MpiComm_->getSize()),0.0);
                reduceAll(*MpiComm_,REDUCE_MIN,local[5],ptr(&min[5]));
                reduceAll(*MpiComm_,REDUCE_MAX,local[5],ptr(&max[5]));
            } else {
                global[5] = -1;
                local[5] = -1;
                avg[5] = -1;
                min[5] = -1;
                max[5] = -1;
            }
            if (buildLeafsMap) {
                global[6] = Leafs_->getEntityMap()->getMaxAllGlobalIndex();
                if (NodesMap_->lib()==UseEpetra || Leafs_->getEntityMap()->getGlobalNumElements()>0) {
                    global[6] += 1;
                }
                if (global[6]<0) global[6] = 0;
                local[6] = (LO) std::max((LO) Leafs_->getEntityMap()->getNodeNumElements(),(LO) 0);
                reduceAll(*this->MpiComm_,REDUCE_SUM,local[6],ptr(&sum[6]));
                avg[6] = std::max(sum[6]/double(MpiComm_->getSize()),0.0);
                reduceAll(*MpiComm_,REDUCE_MIN,local[6],ptr(&min[6]));
                reduceAll(*MpiComm_,REDUCE_MAX,local[6],ptr(&max[6]));
            } else {
                global[6] = -1;
                local[6] = -1;
                avg[6] = -1;
                min[6] = -1;
                max[6] = -1;
            }

            for (UN i=0; i<global.size(); i++) {
                if (global[i]<0) {
                    global[i] = -1;
                }
            }

            if (Verbose_) {
                std::cout << "\n\
    ------------------------------------------------------------------------------\n\
     Interface statistics\n\
    ------------------------------------------------------------------------------\n\
      Vertices:       total / avg / min / max     ---  " << global[0] << " / " << avg[0] << " / " << min[0] << " / " << max[0] << "\n\
      ShortEdges:     total / avg / min / max     ---  " << global[1] << " / " << avg[1] << " / " << min[1] << " / " << max[1] << "\n\
      StraightEdges:  total / avg / min / max     ---  " << global[2] << " / " << avg[2] << " / " << min[2] << " / " << max[2] << "\n\
      Edges:          total / avg / min / max     ---  " << global[3] << " / " << avg[3] << " / " << min[3] << " / " << max[3] << "\n\
      Faces:          total / avg / min / max     ---  " << global[4] << " / " << avg[4] << " / " << min[4] << " / " << max[4] << "\n\
      Roots:          total / avg / min / max     ---  " << global[5] << " / " << avg[5] << " / " << min[5] << " / " << max[5] << "\n\
      Leafs:          total / avg / min / max     ---  " << global[6] << " / " << avg[6] << " / " << min[6] << " / " << max[6] << "\n\
    ------------------------------------------------------------------------------\n";
            }
        }

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::buildEntityHierarchy()
    {
        FROSCH_TIMER_START_LEVELID(buildEntityHierarchyTime,"DDInterface::buildEntityHierarchy");
        //if (Verbose_ && Verbosity_==All) std::cout << "FROSch::DDInterface : Building hierarchy of interface components" << std::endl;

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
        FROSCH_TIMER_START_LEVELID(computeDistancesToRootsTime,"DDInterface::computeDistancesToRoots");
        //if (Verbose_ && Verbosity_==All) std::cout << "FROSch::DDInterface : Computing distances to the coarse nodes" << std::endl;

        for (UN i=0; i<EntitySetVector_.size(); i++) {
            EntitySetVector_[i]->computeDistancesToRoots(dimension,nodeList,distanceFunction);
        }
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::identifyConnectivityEntities(UNVecPtr multiplicities,
                                                               EntityFlagVecPtr flags)
    {
        FROSCH_TIMER_START_LEVELID(identifyConnectivityEntitiesTime,"DDInterface::identifyConnectivityEntities");
        //if (Verbose_ && Verbosity_==All) std::cout << "FROSch::DDInterface : Preparing subdomain graph" << std::endl;

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
                if (std::binary_search(flags.begin(),flags.end(),EntitySetVector_[multiplicities[j]]->getEntity(i)->getEntityFlag())) {
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
                                                             IntVecVec &componentsSubdomainsUnique,
                                                             CommunicationStrategy commStrategy)
    {
        FROSCH_TIMER_START_LEVELID(communicateLocalComponentsTime,"DDInterface::communicateLocalComponents");
        //if (Verbose_ && Verbosity_==All) std::cout << "FROSch::DDInterface : Communicating nodes" << std::endl;

        if (NodesMap_->lib() == UseEpetra && commStrategy == CreateOneToOneMap) {
            if (Verbose_) std::cout << "FROSch::DDInterface : WARNING: CreateOneToOneMap communication strategy does not work for Epetra => Switching to CommCrsGraph" << std::endl;
            commStrategy = CommCrsGraph;
        }

        // Different communication strategies
        switch (commStrategy) {
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
                    XMapPtr rangeMap = MapFactory<LO,GO,NO>::Build(NodesMap_->lib(),-1,myPID(),0,NodesMap_->getComm());

                    commMat->fillComplete(NodesMap_,rangeMap);
                    commMatTmp->doExport(*commMat,*commExporter,INSERT);
                    commMatTmp->fillComplete(UniqueNodesMap_,rangeMap);
                    commMat = MatrixFactory<SC,LO,GO,NO>::Build(NodesMap_,10);
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
                    XMapPtr rangeMap = MapFactory<LO,GO,NO>::Build(NodesMap_->lib(),-1,myPID(),0,NodesMap_->getComm());

                    commGraph->fillComplete(NodesMap_,rangeMap); // AH 08/07/2019: Can we remove some fillComplete?
                    commGraphTmp->doExport(*commGraph,*commExporter,INSERT);
                    commGraphTmp->fillComplete(UniqueNodesMap_,rangeMap);
                    commGraph = CrsGraphFactory<LO,GO,NO>::Build(NodesMap_,10);
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
                FROSCH_ASSERT(false,"FROSch::DDInterface : ERROR: Specify a valid communication strategy.");
                break;
        }

        componentsSubdomainsUnique = IntVecVec(NumMyNodes_);
        for (LO i=0; i<NumMyNodes_; i++) {
            sortunique(componentsSubdomains[i]);
            if (componentsSubdomains[i].size() == 0) componentsSubdomains[i].push_back(MpiComm_->getRank()); // For Tpetra this is empty if the repeatedMap is already unique. In this case, we have to add the local rank. Otherwise, we obtain nodes with multiplicity 0.
            componentsSubdomainsUnique[i] = componentsSubdomains[i];
//            if (MpiComm_->getRank() == 0) std::cout << MpiComm_->getRank() << ": " << i << " " << componentsSubdomains[i] << std::endl;
        }
        sortunique(componentsSubdomainsUnique);

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::identifyLocalComponents(IntVecVecPtr &componentsSubdomains,
                                                          IntVecVec &componentsSubdomainsUnique)
    {
        FROSCH_TIMER_START_LEVELID(identifyLocalComponentsTime,"DDInterface::identifyLocalComponents");
        //if (Verbose_ && Verbosity_==All) std::cout << "FROSch::DDInterface : Classifying interface components based on equivalence classes" << std::endl;

        // Hier herausfinden, ob Ecke, Kante oder Fläche
        UNVecPtr componentsMultiplicity(componentsSubdomainsUnique.size());
        IntVecVecPtr components(componentsSubdomainsUnique.size());
        IntVecVecPtr componentsGamma(componentsSubdomainsUnique.size());
        UN maxMultiplicity = 0;
        for (UN i=0; i<componentsSubdomainsUnique.size(); i++) {
            componentsMultiplicity[i] = componentsSubdomainsUnique[i].size();
            maxMultiplicity = std::max(maxMultiplicity,componentsMultiplicity[i]);
        }
        EntitySetVector_ = EntitySetPtrVecPtr(maxMultiplicity+1);
        for (UN i=0; i<maxMultiplicity+1; i++) {
            EntitySetVector_[i].reset(new EntitySet<SC,LO,GO,NO>(DefaultType));
        }

        typename IntVecVec::iterator classIterator;
        LOVecPtr localComponentIndices(NumMyNodes_);
        for (int i=0; i<NumMyNodes_; i++) {
            classIterator = std::lower_bound(componentsSubdomainsUnique.begin(),componentsSubdomainsUnique.end(),componentsSubdomains[i]);
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
                FROSCH_ASSERT(componentsMultiplicity[localComponentIndices[i]]>1,"FROSch::DDInterface : ERROR: There cannot be any nodes with multiplicity 0.");
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
            FROSCH_ASSERT(componentsMultiplicity[i]>0,"FROSch::DDInterface : ERROR: There cannot be any component with multiplicity 0.");
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
