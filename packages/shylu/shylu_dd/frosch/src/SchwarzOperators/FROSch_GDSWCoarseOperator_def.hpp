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
    void GDSWCoarseOperator<SC,LO,GO,NO>::describe(Teuchos::FancyOStream &out,
                                                   const Teuchos::EVerbosityLevel verbLevel) const
    {
        FROSCH_ASSERT(0!=0,"describe() has be implemented properly...");
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
        this->BlockCoarseMaps_.resize(this->BlockCoarseMaps_.size()+1);
        this->MVPhiGamma_.resize(this->MVPhiGamma_.size()+1);
        this->DofsMaps_.resize(this->DofsMaps_.size()+1);
        this->DofsPerNode_.resize(this->DofsPerNode_.size()+1);
        
        this->NumberOfBlocks_++;
        
        resetCoarseSpaceBlock(this->NumberOfBlocks_-1,dimension,dofsPerNode,nodesMap,dofsMaps,dirichletBoundaryDofs,nodeList);
        
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
        
        Teuchos::ArrayRCP<bool> coarseSpaceFunctions(9);
        
        coarseSpaceFunctions[0] = coarseSpaceList->sublist("Custom").get("Vertices: translations",true);
        
        coarseSpaceFunctions[1] = coarseSpaceList->sublist("Custom").get("ShortEdges: translations",true);
        coarseSpaceFunctions[2] = coarseSpaceList->sublist("Custom").get("ShortEdges: rotations",true);
        
        coarseSpaceFunctions[3] = coarseSpaceList->sublist("Custom").get("StraightEdges: translations",true);
        coarseSpaceFunctions[4] = coarseSpaceList->sublist("Custom").get("StraightEdges: rotations",true);
        
        coarseSpaceFunctions[5] = coarseSpaceList->sublist("Custom").get("Edges: translations",true);
        coarseSpaceFunctions[6] = coarseSpaceList->sublist("Custom").get("Edges: rotations",true);
        
        coarseSpaceFunctions[7] = coarseSpaceList->sublist("Custom").get("Faces: translations",true);
        coarseSpaceFunctions[8] = coarseSpaceList->sublist("Custom").get("Faces: rotations",true);
        
        
        bool useRotations = coarseSpaceList->get("Rotations",true);
        if (useRotations && nodeList.is_null()) {
            useRotations = false;
            if (this->Verbose_) std::cout << "\nWarning: Rotations cannot be used!\n";
        }
        if (!useRotations) {
            coarseSpaceFunctions[2] = false;
            coarseSpaceFunctions[4] = false;
            coarseSpaceFunctions[6] = false;
            coarseSpaceFunctions[8] = false;
        }
        
        this->DofsMaps_[blockId] = dofsMaps;
        this->DofsPerNode_[blockId] = dofsPerNode;
        
        Teuchos::Array<GO> tmpDirichletBoundaryDofs(dirichletBoundaryDofs()); // Here, we do a copy. Maybe, this is not necessary
        sortunique(tmpDirichletBoundaryDofs);
        
        DDInterface_.reset(new DDInterface<SC,LO,GO,NO>(dimension,dofsPerNode,nodesMap));
        DDInterface_->resetGlobalDofs(dofsMaps);
        DDInterface_->removeDirichletNodes(tmpDirichletBoundaryDofs());
        DDInterface_->divideUnconnectedEntities(this->K_);
        
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
        
        if (useForCoarseSpace && (coarseSpaceFunctions[0]||coarseSpaceFunctions[1]||coarseSpaceFunctions[2]||coarseSpaceFunctions[3]||coarseSpaceFunctions[4]||coarseSpaceFunctions[5]||coarseSpaceFunctions[6]||coarseSpaceFunctions[7]||coarseSpaceFunctions[8])) {
            
            ////////////////////////////////
            // Build Processor Map Coarse //
            ////////////////////////////////
            MapPtrVecPtr mapVector( dofsPerNode*(coarseSpaceFunctions[0]+coarseSpaceFunctions[1]+coarseSpaceFunctions[3]+coarseSpaceFunctions[5]+coarseSpaceFunctions[7])+(dofsPerNode-1)*(coarseSpaceFunctions[2]+coarseSpaceFunctions[4]+coarseSpaceFunctions[6]+coarseSpaceFunctions[8])+((dimension==3) && (dofsPerNode==3))*coarseSpaceFunctions[6]+((dimension==3)&&(dofsPerNode==3))*coarseSpaceFunctions[8] ); // Beachte: In 2D gibt es sowieso keine faces
            
            if (coarseSpaceFunctions[0]) {
                vertices = DDInterface_->getVertices();
                vertices->buildEntityMap(nodesMap);
            }
            if (coarseSpaceFunctions[1] || coarseSpaceFunctions[2]) {
                shortEdges = DDInterface_->getShortEdges();
                shortEdges->buildEntityMap(nodesMap);
            }
            if (coarseSpaceFunctions[3] || coarseSpaceFunctions[4]) {
                straightEdges = DDInterface_->getStraightEdges();
                straightEdges->buildEntityMap(nodesMap);
            }
            if (coarseSpaceFunctions[5] || coarseSpaceFunctions[6]) {
                edges = DDInterface_->getEdges();
                edges->buildEntityMap(nodesMap);
            }
            if (coarseSpaceFunctions[7] || coarseSpaceFunctions[8]) {
                faces = DDInterface_->getFaces();
                faces->buildEntityMap(nodesMap);
            }
            
            // Vertices
            int ii=0;
            if (coarseSpaceFunctions[0]) {
                for (UN i=0; i<dofsPerNode; i++) {
                    mapVector[ii] = vertices->getEntityMap();
                    ii++;
                }
            }
            // ShortEdges
            if (coarseSpaceFunctions[1]) {
                for (UN i=0; i<dofsPerNode; i++) {
                    mapVector[ii] = shortEdges->getEntityMap();
                    ii++;
                }
            }
            if (coarseSpaceFunctions[2]) {
                for (UN i=0; i<dofsPerNode-1; i++) {
                    mapVector[ii] = shortEdges->getEntityMap();
                    ii++;
                }
            }
            // StraightEdges
            if (coarseSpaceFunctions[3]) {
                for (UN i=0; i<dofsPerNode; i++) {
                    
                    mapVector[ii] = straightEdges->getEntityMap();
                    ii++;
                }
            }
            if (coarseSpaceFunctions[4]) {
                for (UN i=0; i<dofsPerNode-1; i++) {
                    mapVector[ii] = straightEdges->getEntityMap();
                    ii++;
                    
                }
            }
            // Edges
            if (coarseSpaceFunctions[5]) {
                for (UN i=0; i<dofsPerNode; i++) {
                    mapVector[ii] = edges->getEntityMap();
                    ii++;
                }
            }
            if (coarseSpaceFunctions[6]) {
                for (UN i=0; i<dofsPerNode-1+((dimension==3)&&(dofsPerNode==3)); i++) {
                    mapVector[ii] = edges->getEntityMap();
                    ii++;
                }
            }
            // Faces
            if (coarseSpaceFunctions[7]) {
                for (UN i=0; i<dofsPerNode; i++) {
                    mapVector[ii] = faces->getEntityMap();
                    ii++;
                }
            }
            if (coarseSpaceFunctions[8]) {
                for (UN i=0; i<dofsPerNode-1+((dimension==3)&&(dofsPerNode==3)); i++) { // Beachte: In 2D gibt es sowieso keine faces
                    mapVector[ii] = faces->getEntityMap();
                    ii++;
                }
            }
            
            LOVec numEntitiesGlobal(5);
            if (coarseSpaceFunctions[0]) {
                numEntitiesGlobal[0] = vertices->getEntityMap()->getMaxAllGlobalIndex();
                if (vertices->getEntityMap()->lib()==Xpetra::UseEpetra || vertices->getEntityMap()->getGlobalNumElements()>0) {
                    numEntitiesGlobal[0] += 1;
                }
            } else {
                numEntitiesGlobal[0] = -1;
            }
            if (coarseSpaceFunctions[1] || coarseSpaceFunctions[2]) {
                numEntitiesGlobal[1] = shortEdges->getEntityMap()->getMaxAllGlobalIndex();
                if (shortEdges->getEntityMap()->lib()==Xpetra::UseEpetra || shortEdges->getEntityMap()->getGlobalNumElements()>0) {
                    numEntitiesGlobal[1] += 1;
                }
            } else {
                numEntitiesGlobal[1] = -1;
            }
            if (coarseSpaceFunctions[3] || coarseSpaceFunctions[4]) {
                numEntitiesGlobal[2] = straightEdges->getEntityMap()->getMaxAllGlobalIndex();
                if (straightEdges->getEntityMap()->lib()==Xpetra::UseEpetra || straightEdges->getEntityMap()->getGlobalNumElements()>0) {
                    numEntitiesGlobal[2] += 1;
                }
            } else {
                numEntitiesGlobal[2] = -1;
            }
            if (coarseSpaceFunctions[5] || coarseSpaceFunctions[6]) {
                numEntitiesGlobal[3] = edges->getEntityMap()->getMaxAllGlobalIndex();
                if (edges->getEntityMap()->lib()==Xpetra::UseEpetra || edges->getEntityMap()->getGlobalNumElements()>0) {
                    numEntitiesGlobal[3] += 1;
                }
            } else {
                numEntitiesGlobal[3] = -1;
            }
            if (coarseSpaceFunctions[7] || coarseSpaceFunctions[8]) {
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
                vertices: translations      --- " << coarseSpaceFunctions[0] << "\n\
                shortEdges: translations    --- " << coarseSpaceFunctions[1] << "\n\
                shortEdges: rotations       --- " << coarseSpaceFunctions[2] << "\n\
                straightEdges: translations --- " << coarseSpaceFunctions[3] << "\n\
                straightEdges: rotations    --- " << coarseSpaceFunctions[4] << "\n\
                edges: translations         --- " << coarseSpaceFunctions[5] << "\n\
                edges: rotations            --- " << coarseSpaceFunctions[6] << "\n\
                faces: translations         --- " << coarseSpaceFunctions[7] << "\n\
                faces: rotations            --- " << coarseSpaceFunctions[8] << "\n\
                --------------------------------------------\n";
            }

            LOVecPtr2D partMappings;
            this->BlockCoarseMaps_[blockId] = AssembleMaps(mapVector(),partMappings);

            ////////////////////
            // Build PhiGamma //
            ////////////////////
            phiGammaGDSW(blockId,useRotations,dimension,dofsPerNode,nodeList,partMappings,vertices,shortEdges,straightEdges,edges,faces,coarseSpaceFunctions);
        }
        
        return 0;
    }
    
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::phiGammaGDSW(UN blockId,
                                                      bool buildRotations,
                                                      UN dimension,
                                                      UN dofsPerNode,
                                                      MultiVectorPtr nodeList,
                                                      LOVecPtr2D partMappings, // TODO Kann man hier den 2D Vektor auch ersetzen
                                                      EntitySetPtr vertices,
                                                      EntitySetPtr shortEdges,
                                                      EntitySetPtr straightEdges,
                                                      EntitySetPtr edges,
                                                      EntitySetPtr faces,
                                                      BoolVecPtr coarseSpaceFunctions)
    {
        if (buildRotations) {
            FROSCH_ASSERT(nodeList->getNumVectors()==dimension,"dimension of the nodeList is wrong.");
        }

        //Epetra_SerialComm serialComm;
        MapPtr serialGammaMap = Xpetra::MapFactory<LO,GO,NO>::Build(this->BlockCoarseMaps_[blockId]->lib(),this->GammaDofs_[blockId].size(),0,this->SerialComm_);
        this->MVPhiGamma_[blockId] = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(serialGammaMap,this->BlockCoarseMaps_[blockId]->getNodeNumElements());

        LO ii=0;
        SC x,y,z,rx,ry,rz;
        SCVecPtr dir;
        LO tmp=0;
        // vertices
        if (coarseSpaceFunctions[0]) {
            for (UN k=0; k<dofsPerNode; k++) {
                for (UN i=0; i<vertices->getNumEntities(); i++) {
                    this->MVPhiGamma_[blockId]->replaceLocalValue(vertices->getEntity(i)->getGammaDofID(0,k),partMappings[ii][i],1.0);
                }
                ii++;
            }
        }
        
        // Short edges
        if (coarseSpaceFunctions[1]) { // Translations
            for (UN k=0; k<dofsPerNode; k++) {
                for (UN i=0; i<shortEdges->getNumEntities(); i++) {
                    for (UN j=0; j<shortEdges->getEntity(i)->getNumNodes(); j++) {
                        this->MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,k),partMappings[ii][i],1.0);
                    }
                }
                ii++;
            }
        }

        if (coarseSpaceFunctions[2]) { // Rotations
            
            FROSCH_ASSERT(dofsPerNode>1,"Dofs<2 => Rotations cannot be built.");
            
            if (dimension == 2) {
                for (UN i=0; i<shortEdges->getNumEntities(); i++) {
                    // Rotation 1
                    for (UN j=0; j<shortEdges->getEntity(i)->getNumNodes(); j++) {
                        x = nodeList->getData(0)[shortEdges->getEntity(i)->getLocalNodeID(j)];
                        y = nodeList->getData(1)[shortEdges->getEntity(i)->getLocalNodeID(j)];
                        rx = -y;
                        ry = x;
                        this->MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,0),partMappings[ii][i],rx);
                        this->MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,1),partMappings[ii][i],ry);
                    }
                }
                ii++;
            } else if (dimension == 3) {
                for (UN i=0; i<shortEdges->getNumEntities(); i++) {
                    for (UN j=0; j<shortEdges->getEntity(i)->getNumNodes(); j++) {
                        // Get the direction of the short edge
                        dir = shortEdges->getDirection(dimension,nodeList,i);
                        tmp = 0;
                        
                        FROSCH_ASSERT(sqrt(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2])>1.0e-12,"The direction vector is 0. ERROR!");
                        
                        x = nodeList->getData(0)[shortEdges->getEntity(i)->getLocalNodeID(j)];
                        y = nodeList->getData(1)[shortEdges->getEntity(i)->getLocalNodeID(j)];
                        z = nodeList->getData(2)[shortEdges->getEntity(i)->getLocalNodeID(j)];
                        
                        // Rotation 1
                        if ((fabs(dir[0])>1.0e-12) || (fabs(dir[1])>1.0e-12)) {
                            rx = y;
                            ry = -x;
                            rz = 0;
                            this->MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,0),partMappings[ii+tmp][i],rx);
                            this->MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,1),partMappings[ii+tmp][i],ry);
                            this->MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,2),partMappings[ii+tmp][i],rz);
                            tmp++;
                        }
                        
                        // Rotation 2
                        if ((fabs(dir[0])>1.0e-12) || (fabs(dir[2])>1.0e-12)) {
                            rx = -z;
                            ry = 0;
                            rz = x;
                            this->MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,0),partMappings[ii+tmp][i],rx);
                            this->MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,1),partMappings[ii+tmp][i],ry);
                            this->MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,2),partMappings[ii+tmp][i],rz);
                            tmp++;
                        }
                        
                        // Rotation 3
                        if (((fabs(dir[1])>1.0e-12) || (fabs(dir[2])>1.0e-12)) && (tmp<2)) {
                            rx = 0;
                            ry = z;
                            rz = -y;
                            this->MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,0),partMappings[ii+tmp][i],rx);
                            this->MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,1),partMappings[ii+tmp][i],ry);
                            this->MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,2),partMappings[ii+tmp][i],rz);
                        }
                    }
                }
                ii+=2;
            } else {
                FROSCH_ASSERT(0!=0,"The dimension is neither 2 nor 3!");
            }
        }
        
        
        // Straight edges
        if (coarseSpaceFunctions[3]) { // Translations
            for (UN k=0; k<dofsPerNode; k++) {
                for (UN i=0; i<straightEdges->getNumEntities(); i++) {
                    for (UN j=0; j<straightEdges->getEntity(i)->getNumNodes(); j++) {
                        this->MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,k),partMappings[ii][i],1.0);
                    }
                }
                ii++;
            }
        }
        
        if (coarseSpaceFunctions[4]) { // Rotations
            
            FROSCH_ASSERT(dofsPerNode>1,"Dofs<2 => Rotations cannot be built.");
            
            if (dimension == 2) {
                for (UN i=0; i<straightEdges->getNumEntities(); i++) {
                    // Rotation 1
                    for (UN j=0; j<straightEdges->getEntity(i)->getNumNodes(); j++) {
                        x = nodeList->getData(0)[straightEdges->getEntity(i)->getLocalNodeID(j)];
                        y = nodeList->getData(1)[straightEdges->getEntity(i)->getLocalNodeID(j)];
                        rx = -y;
                        ry = x;
                        this->MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,0),partMappings[ii][i],rx);
                        this->MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,1),partMappings[ii][i],ry);
                    }
                }
                ii++;
            } else if (dimension == 3) {
                for (UN i=0; i<straightEdges->getNumEntities(); i++) {
                    for (UN j=0; j<straightEdges->getEntity(i)->getNumNodes(); j++) {
                        // Get the direction of the straight edge
                        dir = straightEdges->getDirection(dimension,nodeList,i);
                        tmp=0;
                        
                        FROSCH_ASSERT(sqrt(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2])>1.0e-12,"The direction vector is 0. ERROR!");
                        
                        x = nodeList->getData(0)[straightEdges->getEntity(i)->getLocalNodeID(j)];
                        y = nodeList->getData(1)[straightEdges->getEntity(i)->getLocalNodeID(j)];
                        z = nodeList->getData(2)[straightEdges->getEntity(i)->getLocalNodeID(j)];
                        
                        // Rotation 1
                        if ((fabs(dir[0])>1.0e-12) || (fabs(dir[1])>1.0e-12)) {
                            rx = y;
                            ry = -x;
                            rz = 0;
                            this->MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,0),partMappings[ii+tmp][i],rx);
                            this->MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,1),partMappings[ii+tmp][i],ry);
                            this->MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,2),partMappings[ii+tmp][i],rz);
                            tmp++;
                        }
                        
                        // Rotation 2
                        if ((fabs(dir[0])>1.0e-12) || (fabs(dir[2])>1.0e-12)) {
                            rx = -z;
                            ry = 0;
                            rz = x;
                            this->MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,0),partMappings[ii+tmp][i],rx);
                            this->MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,1),partMappings[ii+tmp][i],ry);
                            this->MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,2),partMappings[ii+tmp][i],rz);
                            tmp++;
                        }
                        
                        // Rotation 3
                        if (((fabs(dir[1])>1.0e-12) || (fabs(dir[2])>1.0e-12)) && (tmp<2)) {
                            rx = 0;
                            ry = z;
                            rz = -y;
                            this->MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,0),partMappings[ii+tmp][i],rx);
                            this->MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,1),partMappings[ii+tmp][i],ry);
                            this->MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,2),partMappings[ii+tmp][i],rz);
                        }
                    }
                }
                ii+=2;
            } else {
                FROSCH_ASSERT(0!=0,"The dimension is neither 2 nor 3!");
            }
        }
        
        
        // edges
        if (coarseSpaceFunctions[5]) { // Translations
            for (UN k=0; k<dofsPerNode; k++) {
                for (UN i=0; i<edges->getNumEntities(); i++) {
                    for (UN j=0; j<edges->getEntity(i)->getNumNodes(); j++) {
                        this->MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,k),partMappings[ii][i],1.0);
                    }
                }
                ii++;
            }
        }
        
        if (coarseSpaceFunctions[6]) { // Rotations
            
            FROSCH_ASSERT(dofsPerNode>1,"Dofs<2 => Rotations cannot be built.");
            
            if (dimension == 2) {
                for (UN i=0; i<edges->getNumEntities(); i++) {
                    // Rotation 1
                    for (UN j=0; j<edges->getEntity(i)->getNumNodes(); j++) {
                        x = nodeList->getData(0)[edges->getEntity(i)->getLocalNodeID(j)];
                        y = nodeList->getData(1)[edges->getEntity(i)->getLocalNodeID(j)];
                        rx = -y;
                        ry = x;
                        this->MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,0),partMappings[ii][i],rx);
                        this->MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,1),partMappings[ii][i],ry);
                    }
                }
                ii++;
            } else if (dimension == 3) {
                for (UN i=0; i<edges->getNumEntities(); i++) {
                    for (UN j=0; j<edges->getEntity(i)->getNumNodes(); j++) {
                        x = nodeList->getData(0)[edges->getEntity(i)->getLocalNodeID(j)];
                        y = nodeList->getData(1)[edges->getEntity(i)->getLocalNodeID(j)];
                        z = nodeList->getData(2)[edges->getEntity(i)->getLocalNodeID(j)];
                        
                        // Rotation 1
                        rx = y;
                        ry = -x;
                        rz = 0;
                        this->MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,0),partMappings[ii][i],rx);
                        this->MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,1),partMappings[ii][i],ry);
                        this->MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,2),partMappings[ii][i],rz);
                        
                        // Rotation 2
                        rx = -z;
                        ry = 0;
                        rz = x;
                        this->MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,0),partMappings[ii+1][i],rx);
                        this->MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,1),partMappings[ii+1][i],ry);
                        this->MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,2),partMappings[ii+1][i],rz);
                        
                        // Rotation 3
                        rx = 0;
                        ry = z;
                        rz = -y;
                        this->MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,0),partMappings[ii+2][i],rx);
                        this->MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,1),partMappings[ii+2][i],ry);
                        this->MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,2),partMappings[ii+2][i],rz);
                    }
                }
                ii+=3;
            } else {
                FROSCH_ASSERT(0!=0,"The dimension is neither 2 nor 3!");
            }
        }
        
        
        // faces
        if (coarseSpaceFunctions[7]) { // Translations
            for (UN k=0; k<dofsPerNode; k++) {
                for (UN i=0; i<faces->getNumEntities(); i++) {
                    for (UN j=0; j<faces->getEntity(i)->getNumNodes(); j++) {
                        this->MVPhiGamma_[blockId]->replaceLocalValue(faces->getEntity(i)->getGammaDofID(j,k),partMappings[ii][i],1.0);
                    }
                }
                ii++;
            }
        }
        
        if (coarseSpaceFunctions[8]) { // Rotations
            
            FROSCH_ASSERT(dofsPerNode>1,"Dofs<2 => Rotations cannot be built.");
            
            for (UN i=0; i<faces->getNumEntities(); i++) {
                for (UN j=0; j<faces->getEntity(i)->getNumNodes(); j++) {
                    x = nodeList->getData(0)[faces->getEntity(i)->getLocalNodeID(j)];
                    y = nodeList->getData(1)[faces->getEntity(i)->getLocalNodeID(j)];
                    z = nodeList->getData(2)[faces->getEntity(i)->getLocalNodeID(j)];
                    
                    // Rotation 1
                    rx = y;
                    ry = -x;
                    rz = 0;
                    this->MVPhiGamma_[blockId]->replaceLocalValue(faces->getEntity(i)->getGammaDofID(j,0),partMappings[ii][i],rx);
                    this->MVPhiGamma_[blockId]->replaceLocalValue(faces->getEntity(i)->getGammaDofID(j,1),partMappings[ii][i],ry);
                    this->MVPhiGamma_[blockId]->replaceLocalValue(faces->getEntity(i)->getGammaDofID(j,2),partMappings[ii][i],rz);
                    
                    // Rotation 2
                    rx = -z;
                    ry = 0;
                    rz = x;
                    this->MVPhiGamma_[blockId]->replaceLocalValue(faces->getEntity(i)->getGammaDofID(j,0),partMappings[ii+1][i],rx);
                    this->MVPhiGamma_[blockId]->replaceLocalValue(faces->getEntity(i)->getGammaDofID(j,1),partMappings[ii+1][i],ry);
                    this->MVPhiGamma_[blockId]->replaceLocalValue(faces->getEntity(i)->getGammaDofID(j,2),partMappings[ii+1][i],rz);
                    
                    // Rotation 3
                    rx = 0;
                    ry = z;
                    rz = -y;
                    this->MVPhiGamma_[blockId]->replaceLocalValue(faces->getEntity(i)->getGammaDofID(j,0),partMappings[ii+2][i],rx);
                    this->MVPhiGamma_[blockId]->replaceLocalValue(faces->getEntity(i)->getGammaDofID(j,1),partMappings[ii+2][i],ry);
                    this->MVPhiGamma_[blockId]->replaceLocalValue(faces->getEntity(i)->getGammaDofID(j,2),partMappings[ii+2][i],rz);
                }
            }
            ii++;
        }
        //cout << *this->MVPhiGamma_[blockId];
        return 0;
    }
    
}

#endif
