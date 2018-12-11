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

#ifndef _FROSCH_GDSWCOARSEOPERATOR_DEF_HPP
#define _FROSCH_GDSWCOARSEOPERATOR_DEF_HPP

#include <FROSch_GDSWCoarseOperator_decl.hpp>

namespace FROSch {
    
    template <class SC,class LO,class GO,class NO>
    GDSWCoarseOperator<SC,LO,GO,NO>::GDSWCoarseOperator(CrsMatrixPtr k,
                                                        ParameterListPtr parameterList) :
    HarmonicCoarseOperator<SC,LO,GO,NO> (k,parameterList),
    DDInterface_ ()
    {
        
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                    MapPtr repeatedMap)
    {
        buildCoarseSpace(dimension,repeatedMap);
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                    MapPtr repeatedMap,
                                                    GOVecPtr dirichletBoundaryDofs)
    {
        buildCoarseSpace(dimension,repeatedMap,dirichletBoundaryDofs);
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UN dofsPerNode,
                                                    MapPtr repeatedNodesMap,
                                                    MapPtrVecPtr repeatedDofMaps)
    {
        buildCoarseSpace(dimension,dofsPerNode,repeatedNodesMap,repeatedDofMaps);
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UN dofsPerNode,
                                                    MapPtr repeatedNodesMap,
                                                    MapPtrVecPtr repeatedDofMaps,
                                                    GOVecPtr dirichletBoundaryDofs)
    {
        buildCoarseSpace(dimension,dofsPerNode,repeatedNodesMap,repeatedDofMaps,dirichletBoundaryDofs);
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UN dofsPerNode,
                                                    MapPtr repeatedNodesMap,
                                                    MapPtrVecPtr repeatedDofMaps,
                                                    MultiVectorPtr nodeList)
    {
        buildCoarseSpace(dimension,dofsPerNode,repeatedNodesMap,repeatedDofMaps,nodeList);
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UN dofsPerNode,
                                                    MapPtr repeatedNodesMap,
                                                    MapPtrVecPtr repeatedDofMaps,
                                                    GOVecPtr dirichletBoundaryDofs,
                                                    MultiVectorPtr nodeList)
    {
        buildCoarseSpace(dimension,dofsPerNode,repeatedNodesMap,repeatedDofMaps,dirichletBoundaryDofs,nodeList);
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UNVecPtr dofsPerNodeVec,
                                                    MapPtrVecPtr repeatedNodesMapVec,
                                                    MapPtrVecPtr2D repeatedDofMapsVec,
                                                    GOVecPtr2D dirichletBoundaryDofsVec,
                                                    MultiVectorPtrVecPtr nodeListVec)
    {
        buildCoarseSpace(dimension,dofsPerNodeVec,repeatedNodesMapVec,repeatedDofMapsVec,dirichletBoundaryDofsVec,nodeListVec);
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }

    
    template <class SC,class LO,class GO,class NO>
    void GDSWCoarseOperator<SC,LO,GO,NO>::describe(Teuchos::FancyOStream &out,
                                                   const Teuchos::EVerbosityLevel verbLevel) const
    {
        FROSCH_ASSERT(false,"describe() has be implemented properly...");
    }
    
    template <class SC,class LO,class GO,class NO>
    std::string GDSWCoarseOperator<SC,LO,GO,NO>::description() const
    {
        return "GDSW Coarse Operator";
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                          MapPtr nodesMap)
    {
        MapPtrVecPtr dofsMaps(1);
        dofsMaps[0] = nodesMap;
        buildCoarseSpace(dimension,1,nodesMap,dofsMaps);
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                          MapPtr nodesMap,
                                                          GOVecPtr dirichletBoundaryDofs)
    {
        MapPtrVecPtr dofsMaps(1);
        dofsMaps[0] = nodesMap;
        buildCoarseSpace(dimension,1,nodesMap,dofsMaps,dirichletBoundaryDofs);
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                          UN dofsPerNode,
                                                          MapPtr nodesMap,
                                                          MapPtrVecPtr dofsMaps)
    {
        GOVecPtr dirichletBoundaryDofs = FindOneEntryOnlyRowsGlobal(this->K_,nodesMap);
        buildCoarseSpace(dimension,dofsPerNode,nodesMap,dofsMaps,dirichletBoundaryDofs);
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                          UN dofsPerNode,
                                                          MapPtr nodesMap,
                                                          MapPtrVecPtr dofsMaps,
                                                          GOVecPtr dirichletBoundaryDofs)
    {        
        MultiVectorPtr nodeList;
        buildCoarseSpace(dimension,dofsPerNode,nodesMap,dofsMaps,dirichletBoundaryDofs,nodeList);
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                          UN dofsPerNode,
                                                          MapPtr nodesMap,
                                                          MapPtrVecPtr dofsMaps,
                                                          MultiVectorPtr nodeList)
    {
        GOVecPtr dirichletBoundaryDofs = FindOneEntryOnlyRowsGlobal(this->K_,nodesMap);
        buildCoarseSpace(dimension,dofsPerNode,nodesMap,dofsMaps,dirichletBoundaryDofs,nodeList);
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                          UN dofsPerNode,
                                                          MapPtr nodesMap,
                                                          MapPtrVecPtr dofsMaps,
                                                          GOVecPtr dirichletBoundaryDofs,
                                                          MultiVectorPtr nodeList)
    {
        FROSCH_ASSERT(dofsMaps.size()==dofsPerNode,"dofsMaps.size()!=dofsPerNode");
        
        // Das könnte man noch ändern
        // TODO: DAS SOLLTE ALLES IN EINE FUNKTION IN HARMONICCOARSEOPERATOR
        this->GammaDofs_.resize(this->GammaDofs_.size()+1);
        this->IDofs_.resize(this->IDofs_.size()+1);
        this->InterfaceCoarseSpaces_.resize(this->InterfaceCoarseSpaces_.size()+1);
        this->DofsMaps_.resize(this->DofsMaps_.size()+1);
        this->DofsPerNode_.resize(this->DofsPerNode_.size()+1);
        this->BlockCoarseDimension_.resize(this->BlockCoarseDimension_.size()+1);
        this->NumberOfBlocks_++;
        
        resetCoarseSpaceBlock(this->NumberOfBlocks_-1,dimension,dofsPerNode,nodesMap,dofsMaps,dirichletBoundaryDofs,nodeList);
        
        return 0;
    }
    
    
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                          UNVecPtr dofsPerNodeVec,
                                                          MapPtrVecPtr repeatedNodesMapVec,
                                                          MapPtrVecPtr2D repeatedDofMapsVec,
                                                          GOVecPtr2D dirichletBoundaryDofsVec,
                                                          MultiVectorPtrVecPtr nodeListVec)
    {
        // Das könnte man noch ändern
        // TODO: DAS SOLLTE ALLES IN EINE FUNKTION IN HARMONICCOARSEOPERATOR
        for (UN i=0; i<repeatedNodesMapVec.size(); i++) {
            this->GammaDofs_.resize(this->GammaDofs_.size()+1);
            this->IDofs_.resize(this->IDofs_.size()+1);
            this->InterfaceCoarseSpaces_.resize(this->InterfaceCoarseSpaces_.size()+1);
            this->DofsMaps_.resize(this->DofsMaps_.size()+1);
            this->DofsPerNode_.resize(this->DofsPerNode_.size()+1);
            this->BlockCoarseDimension_.resize(this->BlockCoarseDimension_.size()+1);
            this->NumberOfBlocks_++;
            resetCoarseSpaceBlock(this->NumberOfBlocks_-1,dimension,dofsPerNodeVec[i],repeatedNodesMapVec[i],repeatedDofMapsVec[i],dirichletBoundaryDofsVec[i],nodeListVec[i]);
        }
        return 0;
    }
    
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::resetCoarseSpaceBlock(UN blockId,
                                                               UN dimension,
                                                               UN dofsPerNode,
                                                               MapPtr nodesMap,
                                                               MapPtrVecPtr dofsMaps,
                                                               GOVecPtr dirichletBoundaryDofs,
                                                               MultiVectorPtr nodeList)
    {
        
        FROSCH_ASSERT(dofsMaps.size()==dofsPerNode,"dofsMaps.size()!=dofsPerNode");
        FROSCH_ASSERT(blockId<this->NumberOfBlocks_,"Block does not exist yet and can therefore not be reset.");
        
        // Process the parameter list
        std::stringstream blockIdStringstream;
        blockIdStringstream << blockId+1;
        std::string blockIdString = blockIdStringstream.str();
        Teuchos::RCP<Teuchos::ParameterList> coarseSpaceList = sublist(sublist(this->ParameterList_,"Blocks"),blockIdString.c_str());
        
        bool useForCoarseSpace = coarseSpaceList->get("Use For Coarse Space",true);
        
        bool useVertexTranslations = coarseSpaceList->sublist("Custom").get("Vertices: translations",true);
        
        bool useShortEdgeTranslations = coarseSpaceList->sublist("Custom").get("ShortEdges: translations",true);
        bool useShortEdgeRotations = coarseSpaceList->sublist("Custom").get("ShortEdges: rotations",true);
        
        bool useStraightEdgeTranslations = coarseSpaceList->sublist("Custom").get("StraightEdges: translations",true);
        bool useStraightEdgeRotations = coarseSpaceList->sublist("Custom").get("StraightEdges: rotations",true);
        
        bool useEdgeTranslations = coarseSpaceList->sublist("Custom").get("Edges: translations",true);
        bool useEdgeRotations = coarseSpaceList->sublist("Custom").get("Edges: rotations",true);
        
        bool useFaceTranslations = coarseSpaceList->sublist("Custom").get("Faces: translations",true);
        bool useFaceRotations = coarseSpaceList->sublist("Custom").get("Faces: rotations",true);
        
        bool useRotations = coarseSpaceList->get("Rotations",true);
        if (useRotations && nodeList.is_null()) {
            useRotations = false;
            if (this->Verbose_) std::cout << "\nWarning: Rotations cannot be used!\n";
        }
        if (!useRotations) {
            useShortEdgeRotations = false;
            useStraightEdgeRotations = false;
            useEdgeRotations = false;
            useFaceRotations = false;
        }
        
        this->DofsMaps_[blockId] = dofsMaps;
        this->DofsPerNode_[blockId] = dofsPerNode;

        Teuchos::Array<GO> tmpDirichletBoundaryDofs(dirichletBoundaryDofs()); // Here, we do a copy. Maybe, this is not necessary
        sortunique(tmpDirichletBoundaryDofs);

        DDInterface_.reset(new DDInterface<SC,LO,GO,NO>(dimension,dofsPerNode,nodesMap));
        DDInterface_->resetGlobalDofs(dofsMaps);
        DDInterface_->removeDirichletNodes(tmpDirichletBoundaryDofs());
        if (this->ParameterList_->get("Test Unconnected Interface",true)) {
            DDInterface_->divideUnconnectedEntities(this->K_);
        }

        DDInterface_->sortEntities(nodeList);
        
        EntitySetPtr vertices,shortEdges,straightEdges,edges,faces,interface,interior;
        
        interface = DDInterface_->getInterface();
        interior = DDInterface_->getInterior();
        
        this->GammaDofs_[blockId] = LOVecPtr(dofsPerNode*interface->getEntity(0)->getNumNodes());
        this->IDofs_[blockId] = LOVecPtr(dofsPerNode*interior->getEntity(0)->getNumNodes());
        for (UN k=0; k<dofsPerNode; k++) {
            for (UN i=0; i<interface->getEntity(0)->getNumNodes(); i++) {
                this->GammaDofs_[blockId][dofsPerNode*i+k] = interface->getEntity(0)->getLocalDofID(i,k);
            }
            for (UN i=0; i<interior->getEntity(0)->getNumNodes(); i++) {
                this->IDofs_[blockId][dofsPerNode*i+k] = interior->getEntity(0)->getLocalDofID(i,k);
            }
        }
        
        this->InterfaceCoarseSpaces_[blockId].reset(new CoarseSpace<SC,LO,GO,NO>());
        
        if (useForCoarseSpace && (useVertexTranslations||useShortEdgeTranslations||useShortEdgeRotations||useStraightEdgeTranslations||useStraightEdgeRotations||useEdgeTranslations||useEdgeRotations||useFaceTranslations||useFaceRotations)) {
            
            ////////////////////////////////
            // Build Processor Map Coarse //
            ////////////////////////////////
            // Vertices
            if (useVertexTranslations) {
                vertices = DDInterface_->getVertices();
                vertices->buildEntityMap(nodesMap);
                
                MultiVectorPtrVecPtr translations = computeTranslations(blockId,this->DofsPerNode_[blockId],vertices);
                for (UN i=0; i<translations.size(); i++) {
                    this->InterfaceCoarseSpaces_[blockId]->addSubspace(vertices->getEntityMap(),translations[i]);
                }
            }
            // ShortEdges
            if (useShortEdgeTranslations || useShortEdgeRotations) {
                shortEdges = DDInterface_->getShortEdges();
                shortEdges->buildEntityMap(nodesMap);
                
                if (useShortEdgeTranslations) {
                    MultiVectorPtrVecPtr translations = computeTranslations(blockId,this->DofsPerNode_[blockId],shortEdges);
                    for (UN i=0; i<translations.size(); i++) {
                        this->InterfaceCoarseSpaces_[blockId]->addSubspace(shortEdges->getEntityMap(),translations[i]);
                    }
                }
                if (useShortEdgeRotations) {
                    MultiVectorPtrVecPtr rotations = computeRotations(blockId,dimension,this->DofsPerNode_[blockId],nodeList,shortEdges);
                    for (UN i=0; i<rotations.size(); i++) {
                        this->InterfaceCoarseSpaces_[blockId]->addSubspace(shortEdges->getEntityMap(),rotations[i]);
                    }
                }
            }
            // StraightEdges
            if (useStraightEdgeTranslations || useStraightEdgeRotations) {
                straightEdges = DDInterface_->getStraightEdges();
                straightEdges->buildEntityMap(nodesMap);
                
                if (useShortEdgeTranslations) {
                    MultiVectorPtrVecPtr translations = computeTranslations(blockId,this->DofsPerNode_[blockId],straightEdges);
                    for (UN i=0; i<translations.size(); i++) {
                        this->InterfaceCoarseSpaces_[blockId]->addSubspace(straightEdges->getEntityMap(),translations[i]);
                    }
                }
                if (useShortEdgeRotations) {
                    MultiVectorPtrVecPtr rotations = computeRotations(blockId,dimension,this->DofsPerNode_[blockId],nodeList,straightEdges);
                    for (UN i=0; i<rotations.size(); i++) {
                        this->InterfaceCoarseSpaces_[blockId]->addSubspace(straightEdges->getEntityMap(),rotations[i]);
                    }
                }
            }
            // Edges
            if (useEdgeTranslations || useEdgeRotations) {
                edges = DDInterface_->getEdges();
                edges->buildEntityMap(nodesMap);
                
                if (useShortEdgeTranslations) {
                    MultiVectorPtrVecPtr translations = computeTranslations(blockId,this->DofsPerNode_[blockId],edges);
                    for (UN i=0; i<translations.size(); i++) {
                        this->InterfaceCoarseSpaces_[blockId]->addSubspace(edges->getEntityMap(),translations[i]);
                    }
                }
                if (useShortEdgeRotations) {
                    MultiVectorPtrVecPtr rotations = computeRotations(blockId,dimension,this->DofsPerNode_[blockId],nodeList,edges);
                    for (UN i=0; i<rotations.size(); i++) {
                        this->InterfaceCoarseSpaces_[blockId]->addSubspace(edges->getEntityMap(),rotations[i]);
                    }
                }
            }
            // Faces
            if (useFaceTranslations || useFaceRotations) {
                faces = DDInterface_->getFaces();
                faces->buildEntityMap(nodesMap);
                
                if (useShortEdgeTranslations) {
                    MultiVectorPtrVecPtr translations = computeTranslations(blockId,this->DofsPerNode_[blockId],faces);
                    for (UN i=0; i<translations.size(); i++) {
                        this->InterfaceCoarseSpaces_[blockId]->addSubspace(faces->getEntityMap(),translations[i]);
                    }
                }
                if (useShortEdgeRotations) {
                    MultiVectorPtrVecPtr rotations = computeRotations(blockId,dimension,this->DofsPerNode_[blockId],nodeList,faces);
                    for (UN i=0; i<rotations.size(); i++) {
                        this->InterfaceCoarseSpaces_[blockId]->addSubspace(faces->getEntityMap(),rotations[i]);
                    }
                }
            }
            
            this->InterfaceCoarseSpaces_[blockId]->assembleCoarseSpace();
            
            // Count entities
            LOVec numEntitiesGlobal(5);
            if (useVertexTranslations) {
                numEntitiesGlobal[0] = vertices->getEntityMap()->getMaxAllGlobalIndex();
                if (vertices->getEntityMap()->lib()==Xpetra::UseEpetra || vertices->getEntityMap()->getGlobalNumElements()>0) {
                    numEntitiesGlobal[0] += 1;
                }
            } else {
                numEntitiesGlobal[0] = -1;
            }
            if (useShortEdgeTranslations || useShortEdgeRotations) {
                numEntitiesGlobal[1] = shortEdges->getEntityMap()->getMaxAllGlobalIndex();
                if (shortEdges->getEntityMap()->lib()==Xpetra::UseEpetra || shortEdges->getEntityMap()->getGlobalNumElements()>0) {
                    numEntitiesGlobal[1] += 1;
                }
            } else {
                numEntitiesGlobal[1] = -1;
            }
            if (useStraightEdgeTranslations || useStraightEdgeRotations) {
                numEntitiesGlobal[2] = straightEdges->getEntityMap()->getMaxAllGlobalIndex();
                if (straightEdges->getEntityMap()->lib()==Xpetra::UseEpetra || straightEdges->getEntityMap()->getGlobalNumElements()>0) {
                    numEntitiesGlobal[2] += 1;
                }
            } else {
                numEntitiesGlobal[2] = -1;
            }
            if (useEdgeTranslations || useEdgeRotations) {
                numEntitiesGlobal[3] = edges->getEntityMap()->getMaxAllGlobalIndex();
                if (edges->getEntityMap()->lib()==Xpetra::UseEpetra || edges->getEntityMap()->getGlobalNumElements()>0) {
                    numEntitiesGlobal[3] += 1;
                }
            } else {
                numEntitiesGlobal[3] = -1;
            }
            if (useFaceTranslations || useFaceRotations) {
                numEntitiesGlobal[4] = faces->getEntityMap()->getMaxAllGlobalIndex();
                if (faces->getEntityMap()->lib()==Xpetra::UseEpetra || faces->getEntityMap()->getGlobalNumElements()>0) {
                    numEntitiesGlobal[4] += 1;
                }
            } else {
                numEntitiesGlobal[4] = -1;
            }
            
            for (UN i=0; i<numEntitiesGlobal.size(); i++) {
                if (numEntitiesGlobal[i]<0) {
                    numEntitiesGlobal[i] = 0;
                }
            }

            if (this->Verbose_) {
                
                std::cout << "\n\
                --------------------------------------------\n\
                # vertices:       --- " << numEntitiesGlobal[0] << "\n\
                # shortEdges:     --- " << numEntitiesGlobal[1] << "\n\
                # straightEdges:  --- " << numEntitiesGlobal[2] << "\n\
                # edges:          --- " << numEntitiesGlobal[3] << "\n\
                # faces:          --- " << numEntitiesGlobal[4] << "\n\
                --------------------------------------------\n\
                Coarse space:\n\
                --------------------------------------------\n\
                vertices: translations      --- " << useVertexTranslations << "\n\
                shortEdges: translations    --- " << useShortEdgeTranslations << "\n\
                shortEdges: rotations       --- " << useShortEdgeRotations << "\n\
                straightEdges: translations --- " << useStraightEdgeTranslations << "\n\
                straightEdges: rotations    --- " << useStraightEdgeRotations << "\n\
                edges: translations         --- " << useEdgeTranslations << "\n\
                edges: rotations            --- " << useEdgeRotations << "\n\
                faces: translations         --- " << useFaceTranslations << "\n\
                faces: rotations            --- " << useFaceRotations << "\n\
                --------------------------------------------\n";
            }

            this->BlockCoarseDimension_[blockId] = 0;
            for (UN i=0; i<numEntitiesGlobal.size(); i++) {
                this->BlockCoarseDimension_[blockId] += numEntitiesGlobal[i];
            }
        }
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    typename GDSWCoarseOperator<SC,LO,GO,NO>::MultiVectorPtrVecPtr GDSWCoarseOperator<SC,LO,GO,NO>::computeTranslations(UN blockId,
                                                                                                                        UN dofsPerNode,
                                                                                                                        EntitySetPtr entitySet)
    {
        MultiVectorPtrVecPtr translations(dofsPerNode);
        MapPtr serialGammaMap = Xpetra::MapFactory<LO,GO,NO>::Build(this->K_->getRangeMap()->lib(),this->GammaDofs_[blockId].size(),0,this->SerialComm_);
        for (UN i=0; i<dofsPerNode; i++) {
            if (entitySet->getNumEntities()>0) {
                translations[i] = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(serialGammaMap,entitySet->getNumEntities());
            } else {
                translations[i] = Teuchos::null;
            }
        }
        
        for (UN k=0; k<dofsPerNode; k++) {
            UN ii=0;
            for (UN i=0; i<entitySet->getNumEntities(); i++) {
                for (UN j=0; j<entitySet->getEntity(i)->getNumNodes(); j++) {
                    translations[k]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,k),ii,1.0);
                }
                ii++;
            }
        }
        return translations;
    }
    
    template <class SC,class LO,class GO,class NO>
    typename GDSWCoarseOperator<SC,LO,GO,NO>::MultiVectorPtrVecPtr GDSWCoarseOperator<SC,LO,GO,NO>::computeRotations(UN blockId,
                                                                                                                     UN dimension,
                                                                                                                     UN dofsPerNode,
                                                                                                                     MultiVectorPtr nodeList,
                                                                                                                     EntitySetPtr entitySet)
    {
        FROSCH_ASSERT(nodeList->getNumVectors()==dimension,"dimension of the nodeList is wrong.");
        FROSCH_ASSERT(dimension==dofsPerNode,"dimension!=dofsPerNode");
        
        UN rotationsPerEntity = 0;
        switch (dimension) {
            case 1:
                return Teuchos::null;
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
        
        MultiVectorPtrVecPtr rotations(rotationsPerEntity);
        MapPtr serialGammaMap = Xpetra::MapFactory<LO,GO,NO>::Build(this->K_->getRangeMap()->lib(),this->GammaDofs_[blockId].size(),0,this->SerialComm_);
        for (UN i=0; i<rotationsPerEntity; i++) {
            if (entitySet->getNumEntities()>0) {
                rotations[i] = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(serialGammaMap,entitySet->getNumEntities());
            } else {
                rotations[i] = Teuchos::null;
            }
        }
        
        UN ii=0;
        SC x,y,z,rx,ry,rz;
        SCVecPtr dir;
        for (UN i=0; i<entitySet->getNumEntities(); i++) {
            for (UN j=0; j<entitySet->getEntity(i)->getNumNodes(); j++) {
                x = nodeList->getData(0)[entitySet->getEntity(i)->getLocalNodeID(j)];
                y = nodeList->getData(1)[entitySet->getEntity(i)->getLocalNodeID(j)];
                
                // Rotation 1
                rx = y;
                ry = -x;
                rz = 0;
                rotations[0]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,0),ii,rx);
                rotations[0]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,1),ii,ry);
                if (dimension == 3) {
                    z = nodeList->getData(2)[entitySet->getEntity(i)->getLocalNodeID(j)];
                    
                    rotations[0]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,2),ii,rz);
                
                    // Rotation 2
                    rx = -z;
                    ry = 0;
                    rz = x;
                    rotations[1]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,0),ii,rx);
                    rotations[1]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,1),ii,ry);
                    rotations[1]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,2),ii,rz);
                
                    // Rotation 3
                    rx = 0;
                    ry = z;
                    rz = -y;
                    rotations[2]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,0),ii,rx);
                    rotations[2]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,1),ii,ry);
                    rotations[2]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,2),ii,rz);
                }
            }
            ii++;
        }
        return rotations;
    }
    
}

#endif
