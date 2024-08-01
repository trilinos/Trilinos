// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_ENTITYSET_DEF_HPP
#define _FROSCH_ENTITYSET_DEF_HPP

#include <FROSch_EntitySet_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template<class SC,class LO,class GO,class NO>
    EntitySet<SC,LO,GO,NO>::EntitySet(EntityType type) :
    Type_ (type)
    {

    }

    template<class SC,class LO,class GO,class NO>
    EntitySet<SC,LO,GO,NO>::EntitySet(const EntitySet<SC,LO,GO,NO> &entitySet) :
    Type_ (entitySet.getEntityType()),
    EntityVector_ (entitySet.getEntityVector())
    {

    }

    template<class SC,class LO,class GO,class NO>
    EntitySet<SC,LO,GO,NO>::~EntitySet()
    {

    } // Do we need sth here?

    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::addEntity(InterfaceEntityPtr entity)
    {
        FROSCH_ASSERT(Type_==DefaultType||entity->getEntityType()==Type_,"FROSch::EntitySet: Entity to add is of wrong type.");
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
    typename EntitySet<SC,LO,GO,NO>::EntitySetPtr EntitySet<SC,LO,GO,NO>::deepCopy()
    {
        EntitySetPtr copy(new EntitySet<SC,LO,GO,NO>(Type_));
        for (UN i=0; i<getNumEntities(); i++) {
            InterfaceEntityPtr entity(new InterfaceEntity<SC,LO,GO,NO>(getEntity(i)->getEntityType(),
                                                                       getEntity(i)->getDofsPerNode(),
                                                                       getEntity(i)->getMultiplicity(),
                                                                       getEntity(i)->getSubdomainsVector().getRawPtr()));
            for (UN j=0; j<getEntity(i)->getNumNodes(); j++) {
                entity->addNode(getEntity(i)->getNode(j).NodeIDGamma_,
                                getEntity(i)->getNode(j).NodeIDLocal_,
                                getEntity(i)->getNode(j).NodeIDGlobal_,
                                getEntity(i)->getNode(j).DofsGamma_.size(),
                                getEntity(i)->getNode(j).DofsGamma_,
                                getEntity(i)->getNode(j).DofsLocal_,
                                getEntity(i)->getNode(j).DofsGlobal_);
            }
            copy->addEntity(entity);
        }
        return copy;
    }

    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::buildEntityMap(ConstXMapPtr localToGlobalNodesMap)
    {
        if (!EntityMapIsUpToDate_) {
            LO localNumberEntities = getNumEntities();
            LO globalNumberEntities = 0; // AH 10/13/2017: Can we stick with LO here
            LO maxLocalNumberEntities = 0;
            reduceAll(*localToGlobalNodesMap->getComm(),REDUCE_SUM,localNumberEntities,ptr(&globalNumberEntities));
            reduceAll(*localToGlobalNodesMap->getComm(),REDUCE_MAX,localNumberEntities,ptr(&maxLocalNumberEntities));

            GOVec localToGlobalVector(0);
            const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();
            if (globalNumberEntities>0) {
                // Set the Unique iD
                setUniqueIDToFirstGlobalNodeID();

                GOVec entities(maxLocalNumberEntities);
                for (UN i=0; i<getNumEntities(); i++) {
                    entities[i] = getEntity(i)->getUniqueID()+1;
                    getEntity(i)->setLocalID(i);
                }
                XMapPtr entityMapping = MapFactory<LO,GO,NO>::Build(localToGlobalNodesMap->lib(),INVALID,entities(),0,localToGlobalNodesMap->getComm());

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
            EntityMap_ = MapFactory<LO,GO,NO>::Build(localToGlobalNodesMap->lib(),INVALID,localToGlobalVector(),0,localToGlobalNodesMap->getComm());
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
    typename EntitySet<SC,LO,GO,NO>::EntitySetPtr EntitySet<SC,LO,GO,NO>::findRoots()
    {
        EntitySetPtr Roots(new EntitySet<SC,LO,GO,NO>(DefaultType));
        for (UN i=0; i<getNumEntities(); i++) {
            EntitySetPtr tmpRoots = getEntity(i)->findRoots();
            if (tmpRoots.is_null()) {
                FROSCH_ASSERT(getEntity(i)->getAncestors()->getNumEntities()==0,"FROSch::EntitySet: getEntity(i)->getAncestors()->getNumEntities()!=0");
                Roots->addEntity(getEntity(i));
            } else {
                FROSCH_ASSERT(getEntity(i)->getAncestors()->getNumEntities()!=0,"FROSch::EntitySet: getEntity(i)->getAncestors()->getNumEntities()==0");
                FROSCH_ASSERT(tmpRoots->getNumEntities()>0,"FROSch::EntitySet: tmpRoots->getNumEntities()<=0");
                Roots->addEntitySet(tmpRoots);
            }
        }
        Roots->sortUnique();
        return Roots;
    }

    template<class SC,class LO,class GO,class NO>
    typename EntitySet<SC,LO,GO,NO>::EntitySetPtr EntitySet<SC,LO,GO,NO>::findLeafs()
    {
        EntitySetPtr Leafs(new EntitySet<SC,LO,GO,NO>(DefaultType));
        for (UN i=0; i<getNumEntities(); i++) {
            if (getEntity(i)->getOffspring()->getNumEntities()==0) {
                Leafs->addEntity(getEntity(i));
            }
        }
        Leafs->sortUnique();
        return Leafs;
    }

    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::clearAncestors()
    {
        for (UN i=0; i<getNumEntities(); i++) {
            getEntity(i)->clearAncestors();
        }
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::clearOffspring()
    {
        for (UN i=0; i<getNumEntities(); i++) {
            getEntity(i)->clearOffspring();
        }
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::clearRoots()
    {
        for (UN i=0; i<getNumEntities(); i++) {
            getEntity(i)->clearRoots();
        }
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::clearLeafs()
    {
        for (UN i=0; i<getNumEntities(); i++) {
            getEntity(i)->clearLeafs();
        }
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::computeDistancesToRoots(UN dimension,
                                                        ConstXMultiVectorPtr &nodeList,
                                                        DistanceFunction distanceFunction)
    {
        for (UN i=0; i<getNumEntities(); i++) {
            getEntity(i)->computeDistancesToRoots(dimension,nodeList,distanceFunction);
        }
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::divideUnconnectedEntities(ConstXMatrixPtr matrix,
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
            if (getEntity(i)->getNumNodes()==1) {
                getEntity(i)->resetEntityFlag(NodeFlag);
            }
        }
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::flagShortEntities()
    {
        for (UN i=0; i<getNumEntities(); i++) {
            if (getEntity(i)->getNumNodes()==2) {
                getEntity(i)->resetEntityFlag(ShortFlag);
            }
        }
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::flagStraightEntities(UN dimension,
                                                     ConstXMultiVectorPtr &nodeList)
    {
        FROSCH_ASSERT(dimension==nodeList->getNumVectors(),"FROSch::EntitySet: Inconsistent Dimension.");

        bool straight;
        LO length,j;
        SCVec pt1(dimension);
        SCVec dir1(dimension);
        SCVec dir2(dimension);

        for (UN i=0; i<getNumEntities(); i++) {
            straight = true;
            length = getEntity(i)->getNumNodes();

            j=2;

            if (length>2) {
                // Anfangssteigung berechnen
                for (UN k=0; k<dimension; k++) {
                    pt1[k] = nodeList->getData(k)[getEntity(i)->getLocalNodeID(0)];
                }

                for (UN k=0; k<dimension; k++) {
                    dir1[k] = nodeList->getData(k)[getEntity(i)->getLocalNodeID(1)]-pt1[k];
                }

                while (j<length) {
                    // Steigung zum zweiten Punkt berechnen
                    for (UN k=0; k<dimension; k++) {
                        dir2[k] = nodeList->getData(k)[getEntity(i)->getLocalNodeID(j)]-pt1[k];
                    }

                    if (!ismultiple<SC,LO>(dir1(),dir2())) {
                        straight = false;
                        break;
                    }
                    j++;
                }
                if (straight) {
                    getEntity(i)->resetEntityFlag(StraightFlag);
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
            if (getEntity(i)->getEntityFlag()==flag) {
                removedEntities->addEntities(getEntity(i));
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
        FROSCH_ASSERT(iD<getNumEntities(),"FROSch::EntitySet: Cannot access Entity because iD>=getNumEntities().");
        EntityVector_.erase(EntityVector_.begin()+iD);
        EntityMapIsUpToDate_ = false;
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::removeNodesWithDofs(GOVecView dirichletBoundaryDofs)
    {
        UN dofsPerNode = 0;
        if (getNumEntities()>0) dofsPerNode = EntityVector_[0]->getDofsPerNode();
        for (UN i=0; i<getNumEntities(); i++) {
            UN length = getEntity(i)->getNumNodes();
            for (UN j=0; j<length; j++) {
                UN itmp = length-1-j;
                UN k = 0;
                while (k<dofsPerNode) {
                    GO dofGlobal = getEntity(i)->getGlobalDofID(itmp,k);
                    if (binary_search(dirichletBoundaryDofs.begin(),dirichletBoundaryDofs.end(),dofGlobal)) {
                        getEntity(i)->removeNode(itmp);
                        break;
                    }
                    k++;
                }
            }
        }
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::removeEmptyEntities()
    {
        UN before = getNumEntities();
        for (UN i=0; i<getNumEntities(); i++) {
            if (getEntity(i)->getNumNodes()==0) {
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
            getEntity(i)->sortByGlobalID();
        }

        std::sort(EntityVector_.begin(),EntityVector_.end(),compareInterfaceEntities<SC,LO,GO,NO>);
        EntityVector_.erase(unique(EntityVector_.begin(),EntityVector_.end(),equalInterfaceEntities<SC,LO,GO,NO>),EntityVector_.end());
        EntityMapIsUpToDate_ = false;
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    bool EntitySet<SC,LO,GO,NO>::checkForVertices()
    {
        for (UN i=0; i<getNumEntities(); i++) {
            if (getEntity(i)->getNumNodes()==1) {
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
            if (getEntity(i)->getNumNodes()==2) {
                i--;
                return true;
            }
        }
        return false;
    }

    template<class SC,class LO,class GO,class NO>
    bool EntitySet<SC,LO,GO,NO>::checkForStraightEdges(UN dimension,
                                                       ConstXMultiVectorPtr &nodeList)
    {
        FROSCH_ASSERT(dimension==nodeList->getNumVectors(),"FROSch::EntitySet: Inconsistent Dimension.");

        bool straight;
        LO length,j;
        SCVec pt1(dimension);
        SCVec dir1(dimension);
        SCVec dir2(dimension);

        for (UN i=0; i<getNumEntities(); i++) {
            straight = true;
            length = getEntity(i)->getNumNodes();
            j=2;
            if (length>2) {
                // Anfangssteigung berechnen
                for (UN k=0; k<dimension; k++) {
                    pt1[k] = nodeList->getData(k)[getEntity(i)->getLocalNodeID(0)];
                }
                for (UN k=0; k<dimension; k++) {
                    dir1[k] = nodeList->getData(k)[getEntity(i)->getLocalNodeID(1)]-pt1[k];
                }

                while (j<length) {
                    for (UN k=0; k<dimension; k++) {
                        dir2[k] = nodeList->getData(k)[getEntity(i)->getLocalNodeID(j)]-pt1[k];
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
            if (getEntity(i)->getNumNodes()==0) {
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
            getEntity(i)->sortByGlobalID();
            getEntity(i)->setUniqueIDToFirstGlobalID();
        }
        EntityMapIsUpToDate_ = false;
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::setRootID()
    {
        for (UN i=0; i<getNumEntities(); i++) {
            getEntity(i)->setRootID(i);
        }
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::setLeafID()
    {
        for (UN i=0; i<getNumEntities(); i++) {
            getEntity(i)->setLeafID(i);
        }
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int EntitySet<SC,LO,GO,NO>::resetEntityType(EntityType type)
    {
        Type_ = type;
        for (UN i=0; i<getNumEntities(); i++) {
            getEntity(i)->resetEntityType(type);
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
        FROSCH_ASSERT(iD<getNumEntities(),"FROSch::EntitySet: Cannot access Entity because iD>=getNumEntities().");
        return EntityVector_[iD];
    }

    template<class SC,class LO,class GO,class NO>
    const typename EntitySet<SC,LO,GO,NO>::XMapPtr EntitySet<SC,LO,GO,NO>::getEntityMap() const
    {
        FROSCH_ASSERT(EntityMapIsUpToDate_,"FROSch::EntitySet:  the entity map has not been built or is not up to date.");
        return EntityMap_;
    }

    template<class SC,class LO,class GO,class NO>
    const typename EntitySet<SC,LO,GO,NO>::SCVecPtr EntitySet<SC,LO,GO,NO>::getDirection(UN dimension,
                                                                                         ConstXMultiVectorPtr &nodeList,
                                                                                         UN iD) const
    {
        FROSCH_ASSERT(iD<getNumEntities(),"FROSch::EntitySet: Cannot access Entity because iD>=getNumEntities().");

        if (getEntity(iD)->getEntityFlag()==StraightFlag) {

            LO length = getEntity(iD)->getNumNodes();
            LO j=2;

            FROSCH_ASSERT(length>2,"Edge is not a straight edge!");

            bool straight=true;

            SCVec pt1(dimension);
            SCVec dir2(dimension);
            SCVecPtr dir1(dimension);

            // Anfangssteigung berechnen
            for (UN k=0; k<dimension; k++) {
                pt1[k] = nodeList->getData(k)[getEntity(iD)->getLocalNodeID(0)];
            }
            for (UN k=0; k<dimension; k++) {
                dir1[k] = nodeList->getData(k)[getEntity(iD)->getLocalNodeID(1)]-pt1[k];
            }

            while (j<length) {
                for (UN k=0; k<dimension; k++) {
                    dir2[k] = nodeList->getData(k)[getEntity(iD)->getLocalNodeID(j)]-pt1[k];
                }

                if (!ismultiple<SC,LO>(dir1(),dir2())) {
                    straight = false;
                    break;
                }
                j++;
            }

            FROSCH_ASSERT(straight,"FROSch::EntitySet: Edge is not straight!");

            return dir1;

        } else if (getEntity(iD)->getEntityFlag()==ShortFlag) {

            int length = getEntity(iD)->getNumNodes();

            FROSCH_ASSERT(length==2,"FROSch::EntitySet: Edge is not a short edge!");

            SCVec pt1(dimension);
            SCVecPtr dir1(dimension);

            // Anfangssteigung berechnen
            for (UN k=0; k<dimension; k++) {
                pt1[k] = nodeList->getData(k)[getEntity(iD)->getLocalNodeID(0)];
            }
            for (UN k=0; k<dimension; k++) {
                dir1[k] = nodeList->getData(k)[getEntity(iD)->getLocalNodeID(1)]-pt1[k];
            }

            return dir1;

        } else {
            FROSCH_ASSERT(false,"FROSch::EntitySet: There is a problem while computing the direction of an edge!");
        }
    }
}

#endif
