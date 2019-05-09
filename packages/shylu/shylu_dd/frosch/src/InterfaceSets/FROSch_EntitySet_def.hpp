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

#ifndef _FROSCH_ENTITYSET_DEF_HPP
#define _FROSCH_ENTITYSET_DEF_HPP

#include <FROSch_EntitySet_decl.hpp>

namespace FROSch {
    
    template<class SC,class LO,class GO,class NO>
    EntitySet<SC,LO,GO,NO>::EntitySet(EntityType type):
    Type_ (type),
    EntityVector_ (0),
    EntityMapIsUpToDate_ (false),
    EntityMap_ ()
    {
        
    }
    
    template<class SC,class LO,class GO,class NO>
    EntitySet<SC,LO,GO,NO>::EntitySet(const EntitySet<SC,LO,GO,NO> &entitySet) :
    Type_ (entitySet.getEntityType()),
    EntityVector_ (entitySet.getEntityVector()),
    EntityMapIsUpToDate_ (false),
    EntityMap_ ()
    {
        
    }
    
    template<class SC,class LO,class GO,class NO>
    EntitySet<SC,LO,GO,NO>::~EntitySet()
    {
        
    } // Do we need sth here?
    
    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::addEntity(InterfaceEntityPtr entity)
    {
        FROSCH_ASSERT(Type_==DefaultType||entity->getEntityType()==Type_,"Entity to add is of wrong type.");
        EntityVector_.push_back(entity);
        EntityMapIsUpToDate_ = false;
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::addEntitySet(EntitySetPtr entitySet)
    {
        for (UN i=0; i<entitySet->getNumEntities(); i++) {
            addEntity(entitySet->getEntity(i));
        }
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::buildEntityMap(ConstMapPtr localToGlobalNodesMap)
    {
        if (!EntityMapIsUpToDate_) {
            LO localNumberEntities = getNumEntities();
            LO globalNumberEntities = 0; // AH 10/13/2017: Can we stick with LO here
            LO maxLocalNumberEntities = 0;
            reduceAll(*localToGlobalNodesMap->getComm(),Teuchos::REDUCE_SUM,localNumberEntities,Teuchos::ptr(&globalNumberEntities));
            reduceAll(*localToGlobalNodesMap->getComm(),Teuchos::REDUCE_MAX,localNumberEntities,Teuchos::ptr(&maxLocalNumberEntities));
            
            GOVec localToGlobalVector(0);
            if (globalNumberEntities>0) {
                // Set the Unique iD
                setUniqueIDToFirstGlobalNodeID();
                
                GOVec entities(maxLocalNumberEntities);
                for (UN i=0; i<getNumEntities(); i++) {
                    entities[i] = EntityVector_[i]->getUniqueID()+1;
                    EntityVector_[i]->setLocalID(i);
                }
                MapPtr entityMapping = Xpetra::MapFactory<LO,GO,NO>::Build(localToGlobalNodesMap->lib(),-1,entities(),0,localToGlobalNodesMap->getComm());
                
                GOVec allEntities(maxLocalNumberEntities*localToGlobalNodesMap->getComm()->getSize(),0);
                //localToGlobalNodesMap->getComm().GatherAll(&(entities->at(0)),&(allEntities->at(0)),maxLocalNumberEntities);
                gatherAll(*localToGlobalNodesMap->getComm(),maxLocalNumberEntities,entities.getRawPtr(),maxLocalNumberEntities*localToGlobalNodesMap->getComm()->getSize(),allEntities.getRawPtr());
                
                allEntities.push_back(0); // Um sicherzugehen, dass der erste Eintrag nach sort_unique eine 0 ist.
                
                sortunique(allEntities);
                
                localToGlobalVector.resize(localNumberEntities);
                int LocalID;
                for (UN i=1; i<allEntities.size(); i++) { // Wir fangen bei 1 an, weil wir am Anfang 1 auf die ID addiert haben
                    LocalID = entityMapping->getLocalElement(allEntities[i]);
                    if ( LocalID != -1) {
                        localToGlobalVector[LocalID] = i-1;
                    }
                }
                
            }
            EntityMap_ = Xpetra::MapFactory<LO,GO,NO>::Build(localToGlobalNodesMap->lib(),-1,localToGlobalVector(),0,localToGlobalNodesMap->getComm());
            EntityMapIsUpToDate_ = true;
        }
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::findAncestorsInSet(EntitySetPtr entitySet)
    {
        for (UN i=0; i<getNumEntities(); i++) {            
            getEntity(i)->findAncestorsInSet(entitySet);
        }
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::clearAncestors()
    {
        for (UN i=0; i<getNumEntities(); i++) {
            EntityVector_[i]->clearAncestors();
        }
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::clearOffspring()
    {
        for (UN i=0; i<getNumEntities(); i++) {
            EntityVector_[i]->clearOffspring();
        }
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    typename EntitySet<SC,LO,GO,NO>::EntitySetPtr EntitySet<SC,LO,GO,NO>::findCoarseNodes()
    {
        EntitySetPtr coarseNodes(new EntitySet<SC,LO,GO,NO>(DefaultType));
        for (UN i=0; i<getNumEntities(); i++) {
            EntitySetPtr tmpCoarseNodes = EntityVector_[i]->findCoarseNodes();
            if (tmpCoarseNodes.is_null()) {
                FROSCH_ASSERT(EntityVector_[i]->getAncestors()->getNumEntities()==0,"EntityVector_[i]->getAncestors()->getNumEntities()!=0");
                coarseNodes->addEntity(EntityVector_[i]);
            } else {
                FROSCH_ASSERT(EntityVector_[i]->getAncestors()->getNumEntities()!=0,"EntityVector_[i]->getAncestors()->getNumEntities()==0");
                FROSCH_ASSERT(tmpCoarseNodes->getNumEntities()>0,"tmpCoarseNodes->getNumEntities()<=0");
                coarseNodes->addEntitySet(tmpCoarseNodes);
            }
        }
        coarseNodes->sortUnique();
        return coarseNodes;
    }
    
    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::clearCoarseNodes()
    {
        for (UN i=0; i<getNumEntities(); i++) {
            EntityVector_[i]->clearCoarseNodes();
        }
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::computeDistancesToCoarseNodes(UN dimension,
                                                              MultiVectorPtr &nodeList,
                                                              DistanceFunction distanceFunction)
    {
        for (UN i=0; i<getNumEntities(); i++) {
            EntityVector_[i]->computeDistancesToCoarseNodes(dimension,nodeList,distanceFunction);
        }
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::divideUnconnectedEntities(CrsMatrixPtr matrix,
                                                          int pID)
    {
        UN before = getNumEntities();
        UN i=0;
        while (i<getNumEntities()) {
            InterfaceEntityPtr tmpEntity = getEntity(i)->divideEntity(matrix,pID);
            if (tmpEntity->getNumNodes()>0) {
                addEntity(tmpEntity);
            }
            i++;
        }
        if (getNumEntities()-before>0) EntityMapIsUpToDate_ = false;
        return getNumEntities()-before;
    }
    
    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::flagNodes()
    {
        for (UN i=0; i<getNumEntities(); i++) {
            if (EntityVector_[i]->getNumNodes()==1) {
                EntityVector_[i]->resetEntityFlag(NodeFlag);
            }
        }
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::flagShortEntities()
    {
        for (UN i=0; i<getNumEntities(); i++) {
            if (EntityVector_[i]->getNumNodes()==2) {
                EntityVector_[i]->resetEntityFlag(ShortFlag);
            }
        }
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::flagStraightEntities(UN dimension,
                                                     MultiVectorPtr &nodeList)
    {
        FROSCH_ASSERT(dimension==nodeList->getNumVectors(),"Inconsistent Dimension.");
        
        bool straight;
        LO length,j;
        SCVec pt1(dimension);
        SCVec dir1(dimension);
        SCVec dir2(dimension);
        
        for (UN i=0; i<getNumEntities(); i++) {
            straight = true;
            length = EntityVector_[i]->getNumNodes();
            
            j=2;
            
            if (length>2) {
                // Anfangssteigung berechnen
                for (UN k=0; k<dimension; k++) {
                    pt1[k] = nodeList->getData(k)[EntityVector_[i]->getLocalNodeID(0)];
                }
                
                for (UN k=0; k<dimension; k++) {
                    dir1[k] = nodeList->getData(k)[EntityVector_[i]->getLocalNodeID(1)]-pt1[k];
                }
                
                while (j<length) {
                    // Steigung zum zweiten Punkt berechnen
                    for (UN k=0; k<dimension; k++) {
                        dir2[k] = nodeList->getData(k)[EntityVector_[i]->getLocalNodeID(j)]-pt1[k];
                    }
                    
                    if (!ismultiple<SC,LO>(dir1(),dir2())) {
                        straight = false;
                        break;
                    }
                    j++;
                }
                if (straight) {
                    EntityVector_[i]->resetEntityFlag(StraightFlag);
                }
            }
        }
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    typename EntitySet<SC,LO,GO,NO>::EntitySetPtr EntitySet<SC,LO,GO,NO>::sortOutEntities(EntityFlag flag)
    {
        UN before = getNumEntities();
        EntitySetPtr removedEntities(new EntitySet<SC,LO,GO,NO>(DefaultType));
        for (UN i=0; i<getNumEntities(); i++) {
            if (EntityVector_[i]->getEntityFlag()==flag) {
                removedEntities->addEntities(EntityVector_[i]);
                EntityVector_.erase(EntityVector_.begin()+i);
                i--;
            }
        }
        if (getNumEntities()-before>0) EntityMapIsUpToDate_ = false;
        return arcp(removedEntities);
    }
    
    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::removeEntity(UN iD)
    {
        FROSCH_ASSERT(iD<getNumEntities(),"iD is larger than the number of entities.");
        EntityVector_.erase(EntityVector_.begin()+iD);
        EntityMapIsUpToDate_ = false;
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::removeEmptyEntities()
    {
        UN before = getNumEntities();
        for (UN i=0; i<getNumEntities(); i++) {
            if (EntityVector_[i]->getNumNodes()==0) {
                EntityVector_.erase(EntityVector_.begin()+i);
                i--;
            }
        }
        if (getNumEntities()-before>0) EntityMapIsUpToDate_ = false;
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::sortUnique()
    {
        for (UN i=0; i<getNumEntities(); i++) {
            EntityVector_[i]->sortByGlobalID();
        }
        
        std::sort(EntityVector_.begin(),EntityVector_.end(),compareInterfaceEntities<SC,LO,GO,NO>);
        EntityVector_.erase(std::unique(EntityVector_.begin(),EntityVector_.end(),equalInterfaceEntities<SC,LO,GO,NO>),EntityVector_.end());
        EntityMapIsUpToDate_ = false;
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    bool EntitySet<SC,LO,GO,NO>::checkForVertices()
    {
        for (UN i=0; i<getNumEntities(); i++) {
            if (EntityVector_[i]->getNumNodes()==1) {
                i--;
                return true;
            }
        }
        return false;
    }
    
    template<class SC,class LO,class GO,class NO>
    bool EntitySet<SC,LO,GO,NO>::checkForShortEdges()
    {
        for (UN i=0; i<getNumEntities(); i++) {
            if (EntityVector_[i]->getNumNodes()==2) {
                i--;
                return true;
            }
        }
        return false;
    }
    
    template<class SC,class LO,class GO,class NO>
    bool EntitySet<SC,LO,GO,NO>::checkForStraightEdges(UN dimension,
                                                       MultiVectorPtr &nodeList)
    {
        FROSCH_ASSERT(dimension==nodeList->getNumVectors(),"Inconsistent Dimension.");
        
        bool straight;
        LO length,j;
        SCVec pt1(dimension);
        SCVec dir1(dimension);
        SCVec dir2(dimension);
        
        for (UN i=0; i<getNumEntities(); i++) {
            straight = true;
            length = EntityVector_[i]->getNumNodes();
            j=2;
            if (length>2) {
                // Anfangssteigung berechnen
                for (UN k=0; k<dimension; k++) {
                    pt1[k] = nodeList->getData(k)[EntityVector_[i]->getLocalNodeID(0)];
                }
                for (UN k=0; k<dimension; k++) {
                    dir1[k] = nodeList->getData(k)[EntityVector_[i]->getLocalNodeID(1)]-pt1[k];
                }
                
                while (j<length) {
                    for (UN k=0; k<dimension; k++) {
                        dir2[k] = nodeList->getData(k)[EntityVector_[i]->getLocalNodeID(j)]-pt1[k];
                    }
                    if (!ismultiple<SC,LO>(dir1(),dir2())) {
                        straight = false;
                        break;
                    }
                    j++;
                }
                if (straight) {
                    i--;
                    return true;
                }
            }
        }
        return false;
    }
    
    template<class SC,class LO,class GO,class NO>
    bool EntitySet<SC,LO,GO,NO>::checkForEmptyEntities()
    {
        for (UN i=0; i<getNumEntities(); i++) {
            if (EntityVector_[i]->getNumNodes()==0) {
                i--;
                return true;
            }
        }
        return false;
    }
    
    /////////////////
    // Set Methods //
    /////////////////
    
    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::setUniqueIDToFirstGlobalNodeID()
    {
        for (UN i=0; i<getNumEntities(); i++) {
            EntityVector_[i]->sortByGlobalID();
            EntityVector_[i]->setUniqueIDToFirstGlobalID();
        }
        EntityMapIsUpToDate_ = false;
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::setCoarseNodeID()
    {
        for (UN i=0; i<getNumEntities(); i++) {
            EntityVector_[i]->setCoarseNodeID(i);
        }
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::resetEntityType(EntityType type)
    {
        Type_ = type;
        for (UN i=0; i<getNumEntities(); i++) {
            EntityVector_[i]->resetEntityType(type);
        }
        return 0;
    }
    
    /////////////////
    // Get Methods //
    /////////////////
    
    template<class SC,class LO,class GO,class NO>
    EntityType EntitySet<SC,LO,GO,NO>::getEntityType() const
    {
        return Type_;
    }
    
    template<class SC,class LO,class GO,class NO>
    typename EntitySet<SC,LO,GO,NO>::UN EntitySet<SC,LO,GO,NO>::getNumEntities() const
    {
        return EntityVector_.size();
    }
    
    template<class SC,class LO,class GO,class NO>
    const typename EntitySet<SC,LO,GO,NO>::InterfaceEntityPtrVec & EntitySet<SC,LO,GO,NO>::getEntityVector() const
    {
        return EntityVector_;
    }
    
    template<class SC,class LO,class GO,class NO>
    const typename EntitySet<SC,LO,GO,NO>::InterfaceEntityPtr EntitySet<SC,LO,GO,NO>::getEntity(UN iD) const
    {
        FROSCH_ASSERT(iD<getNumEntities(),"iD>=getNumEntities().");
        return EntityVector_[iD];
    }
    
    template<class SC,class LO,class GO,class NO>
    const typename EntitySet<SC,LO,GO,NO>::MapPtr EntitySet<SC,LO,GO,NO>::getEntityMap() const
    {
        FROSCH_ASSERT(EntityMapIsUpToDate_,"EntitySet: the entity map has not been built or is not up to date.");
        return EntityMap_;
    }
    
    template<class SC,class LO,class GO,class NO>
    const typename EntitySet<SC,LO,GO,NO>::SCVecPtr EntitySet<SC,LO,GO,NO>::getDirection(UN dimension,
                                                                                         MultiVectorPtr &nodeList,
                                                                                         UN iD) const
    {
        FROSCH_ASSERT(iD<getNumEntities(),"iD>=getNumEntities().");
        
        if (EntityVector_[iD]->getEntityFlag()==StraightFlag) {
            
            LO length = EntityVector_[iD]->getNumNodes();
            LO j=2;
            
            FROSCH_ASSERT(length>2,"Edge is not a straight edge!");
            
            bool straight=true;
            
            SCVec pt1(dimension);
            SCVec dir2(dimension);
            SCVecPtr dir1(dimension);
            
            // Anfangssteigung berechnen
            for (UN k=0; k<dimension; k++) {
                pt1[k] = nodeList->getData(k)[EntityVector_[iD]->getLocalNodeID(0)];
            }
            for (UN k=0; k<dimension; k++) {
                dir1[k] = nodeList->getData(k)[EntityVector_[iD]->getLocalNodeID(1)]-pt1[k];
            }
            
            while (j<length) {
                for (UN k=0; k<dimension; k++) {
                    dir2[k] = nodeList->getData(k)[EntityVector_[iD]->getLocalNodeID(j)]-pt1[k];
                }
                
                if (!ismultiple<SC,LO>(dir1(),dir2())) {
                    straight = false;
                    break;
                }
                j++;
            }
            
            FROSCH_ASSERT(straight,"Edge is not straight!");
            
            return dir1;
            
        } else if (EntityVector_[iD]->getEntityFlag()==ShortFlag) {
            
            int length = EntityVector_[iD]->getNumNodes();
            
            FROSCH_ASSERT(length==2,"Edge is not a short edge!");
            
            SCVec pt1(dimension);
            SCVecPtr dir1(dimension);
            
            // Anfangssteigung berechnen
            for (UN k=0; k<dimension; k++) {
                pt1[k] = nodeList->getData(k)[EntityVector_[iD]->getLocalNodeID(0)];
            }
            for (UN k=0; k<dimension; k++) {
                dir1[k] = nodeList->getData(k)[EntityVector_[iD]->getLocalNodeID(1)]-pt1[k];
            }
            
            return dir1;
            
        } else {
            FROSCH_ASSERT(false,"There is a problem while computing the direction of an edge!");
        }
    }
}

#endif
