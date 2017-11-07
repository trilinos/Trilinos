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
    CoarseOperator<SC,LO,GO,NO> (k,parameterList),
    DDInterface_ (),
    ExtensionSolver_ (),
    MVPhiGamma_ (0),
    BlockCoarseMaps_ (0),
    Dimensions_ (0),
    DofsPerNode_ (0),
    IndicesGamma_ (0),
    IndicesI_ (0),
    IndicesGammaDofs_ (0),
    IndicesIDofs_ (0),
    DofsMaps_ (0),
    NumberOfBlocks_ (0)
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
                                                    GOVecPtr myGlobalDirichletBoundaryDofs)
    {
        buildCoarseSpace(dimension,repeatedMap,myGlobalDirichletBoundaryDofs);
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UN dofsPerNode,
                                                    MapPtr repeatedNodesMap,
                                                    MapPtrVecPtr &repeatedDofMaps)
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
                                                    MapPtrVecPtr &repeatedDofMaps,
                                                    GOVecPtr myGlobalDirichletBoundaryDofs)
    {
        buildCoarseSpace(dimension,dofsPerNode,repeatedNodesMap,repeatedDofMaps,myGlobalDirichletBoundaryDofs);
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UN dofsPerNode,
                                                    MapPtr repeatedNodesMap,
                                                    MapPtrVecPtr &repeatedDofMaps,
                                                    SCVecPtr2D &localNodeList)
    {
        buildCoarseSpace(dimension,dofsPerNode,repeatedNodesMap,repeatedDofMaps,localNodeList);
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UN dofsPerNode,
                                                    MapPtr repeatedNodesMap,
                                                    MapPtrVecPtr &repeatedDofMaps,
                                                    GOVecPtr myGlobalDirichletBoundaryDofs,
                                                    SCVecPtr2D &localNodeList)
    {
        buildCoarseSpace(dimension,dofsPerNode,repeatedNodesMap,repeatedDofMaps,myGlobalDirichletBoundaryDofs,localNodeList);
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::compute()
    {
        // This is not optimal yet... Some work could be moved to Initialize
        if (this->Verbose_) {
            cerr << "WARNING: Some of the operations could be moved from initialize() to Compute().\n";
        }
        this->computeBasis();
        this->setUpCoarseOperator(); 
        this->IsComputed_ = true;
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
                                                          MapPtr &nodesMap)
    {
        MapPtrVecPtr dofsMaps(1);
        dofsMaps[0] = nodesMap;
        buildCoarseSpace(dimension,1,nodesMap,dofsMaps);
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                          MapPtr &nodesMap,
                                                          GOVecPtr myGlobalDirichletBoundaryDofs)
    {
        MapPtrVecPtr dofsMaps(1);
        dofsMaps[0] = nodesMap;
        buildCoarseSpace(dimension,1,nodesMap,dofsMaps,myGlobalDirichletBoundaryDofs);
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                          UN dofsPerNode,
                                                          MapPtr &nodesMap,
                                                          MapPtrVecPtr &dofsMaps)
    {
        GOVecPtr myGlobalDirichletBoundaryDofs = FindOneEntryOnlyRowsGlobal(this->K_,nodesMap);
        buildCoarseSpace(dimension,dofsPerNode,nodesMap,dofsMaps,myGlobalDirichletBoundaryDofs);
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                          UN dofsPerNode,
                                                          MapPtr &nodesMap,
                                                          MapPtrVecPtr &dofsMaps,
                                                          GOVecPtr myGlobalDirichletBoundaryDofs)
    {
        
        SCVecPtr2D localNodeList(0);
        buildCoarseSpace(dimension,dofsPerNode,nodesMap,dofsMaps,myGlobalDirichletBoundaryDofs,localNodeList);
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                          UN dofsPerNode,
                                                          MapPtr &nodesMap,
                                                          MapPtrVecPtr &dofsMaps,
                                                          SCVecPtr2D &localNodeList)
    {
        GOVecPtr myGlobalDirichletBoundaryDofs = FindOneEntryOnlyRowsGlobal(this->K_,nodesMap);
        buildCoarseSpace(dimension,dofsPerNode,nodesMap,dofsMaps,myGlobalDirichletBoundaryDofs,localNodeList);
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                          UN dofsPerNode,
                                                          MapPtr &nodesMap,
                                                          MapPtrVecPtr &dofsMaps,
                                                          GOVecPtr myGlobalDirichletBoundaryDofs,
                                                          SCVecPtr2D &localNodeList)
    {
        FROSCH_ASSERT(dofsMaps.size()==dofsPerNode,"dofsMaps.size()!=dofsPerNode");
        
        // Das könnte man noch ändern
        IndicesGamma_.resize(IndicesGamma_.size()+1);
        IndicesI_.resize(IndicesI_.size()+1);
        IndicesGammaDofs_.resize(IndicesGammaDofs_.size()+1);
        IndicesIDofs_.resize(IndicesIDofs_.size()+1);
        BlockCoarseMaps_.resize(BlockCoarseMaps_.size()+1);
        MVPhiGamma_.resize(MVPhiGamma_.size()+1);
        DofsMaps_.resize(DofsMaps_.size()+1);
        DofsPerNode_.resize(DofsPerNode_.size()+1);
        
        NumberOfBlocks_++;
        
        resetCoarseSpaceBlock(NumberOfBlocks_-1,dimension,dofsPerNode,nodesMap,dofsMaps,myGlobalDirichletBoundaryDofs,localNodeList);
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::resetCoarseSpaceBlock(UN blockId,
                                                               UN dimension,
                                                               UN dofsPerNode,
                                                               MapPtr &nodesMap,
                                                               MapPtrVecPtr &dofsMaps,
                                                               GOVecPtr &myGlobalDirichletBoundaryDofs,
                                                               SCVecPtr2D &localNodeList)
    {
        FROSCH_ASSERT(dofsMaps.size()==dofsPerNode,"dofsMaps.size()!=dofsPerNode");
        FROSCH_ASSERT(blockId<NumberOfBlocks_,"Block does not exist yet and can therefore not be reset.");
        
        // Process the parameter list
        std::stringstream blockIdStringstream;
        blockIdStringstream << blockId+1;
        std::string blockIdString = blockIdStringstream.str();
        Teuchos::RCP<Teuchos::ParameterList> coarseSpaceList = sublist(sublist(this->ParameterList_,"Blocks"),blockIdString.c_str());
        
        bool useForCoarseSpace = coarseSpaceList->get("Use For Coarse Space",true);
        
        Teuchos::ArrayRCP<bool> coarseSpaceFunctions(9);
        
        coarseSpaceFunctions[0] = coarseSpaceList->sublist("Custom").get("vertices: translations",true);
        
        coarseSpaceFunctions[1] = coarseSpaceList->sublist("Custom").get("shortEdges: translations",true);
        coarseSpaceFunctions[2] = coarseSpaceList->sublist("Custom").get("shortEdges: rotations",true);
        
        coarseSpaceFunctions[3] = coarseSpaceList->sublist("Custom").get("straightEdges: translations",true);
        coarseSpaceFunctions[4] = coarseSpaceList->sublist("Custom").get("straightEdges: rotations",true);
        
        coarseSpaceFunctions[5] = coarseSpaceList->sublist("Custom").get("edges: translations",true);
        coarseSpaceFunctions[6] = coarseSpaceList->sublist("Custom").get("edges: rotations",true);
        
        coarseSpaceFunctions[7] = coarseSpaceList->sublist("Custom").get("faces: translations",true);
        coarseSpaceFunctions[8] = coarseSpaceList->sublist("Custom").get("faces: rotations",true);
        
        
        bool useRotations = coarseSpaceList->get("Rotations",true);
        if (useRotations && (localNodeList.size()==0) && this->Verbose_) {
            useRotations = false;
            if (this->Verbose_) std::cout << "\nWarning: Rotations cannot be used!\n";
        }
        if (!useRotations) {
            coarseSpaceFunctions[2] = false;
            coarseSpaceFunctions[4] = false;
            coarseSpaceFunctions[6] = false;
            coarseSpaceFunctions[8] = false;
        }
        
        DofsMaps_[blockId] = dofsMaps;
        DofsPerNode_[blockId] = dofsPerNode;
        
        Teuchos::Array<GO> globalDirichletBoundaryDofs(myGlobalDirichletBoundaryDofs()); // Here, we do a copy. Maybe, this is not necessary
        sortunique(globalDirichletBoundaryDofs);
        
        DDInterface_.reset(new DDInterface<SC,LO,GO,NO>(dimension,dofsPerNode,nodesMap));
        DDInterface_->resetGlobalDofs(dofsMaps);
        DDInterface_->removeDirichletNodes(globalDirichletBoundaryDofs());
        DDInterface_->divideUnconnectedEntities(this->K_);
        if (localNodeList.size()>0) {
            DDInterface_->sortEntities(localNodeList);
        } else {
            DDInterface_->sortEntities();
        }
        
        EntitySetPtr vertices,shortEdges,straightEdges,edges,faces,interface,interior;
        
        interface = DDInterface_->getInterface();
        interior = DDInterface_->getInterior();
        
        IndicesGammaDofs_[blockId] = LOVecPtr(dofsPerNode*interface->getEntity(0)->getNumNodes());
        IndicesIDofs_[blockId] = LOVecPtr(dofsPerNode*interior->getEntity(0)->getNumNodes());
        for (UN k=0; k<dofsPerNode; k++) {
            for (UN i=0; i<interface->getEntity(0)->getNumNodes(); i++) {
                IndicesGammaDofs_[blockId][dofsPerNode*i+k] = interface->getEntity(0)->getLocalDofID(i,k);
            }
            for (UN i=0; i<interior->getEntity(0)->getNumNodes(); i++) {
                IndicesIDofs_[blockId][dofsPerNode*i+k] = interior->getEntity(0)->getLocalDofID(i,k);
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
            BlockCoarseMaps_[blockId] = AssembleMaps(mapVector,partMappings);

            ////////////////////
            // Build PhiGamma //
            ////////////////////
            phiGammaGDSW(blockId,useRotations,dimension,dofsPerNode,localNodeList,partMappings,vertices,shortEdges,straightEdges,edges,faces,coarseSpaceFunctions);
        }
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::addZeroCoarseSpaceBlock(MapPtr &dofsMap)
    {
        // Das könnte man noch ändern
        IndicesGamma_->resize(IndicesGamma_.size()+1);
        IndicesI_->resize(IndicesI_.size()+1);
        IndicesGammaDofs_->resize(IndicesGammaDofs_.size()+1);
        IndicesIDofs_->resize(IndicesIDofs_.size()+1);
        BlockCoarseMaps_->resize(BlockCoarseMaps_.size()+1);
        MVPhiGamma_->resize(MVPhiGamma_.size()+1);
        DofsMaps_->resize(DofsMaps_.size()+1);
        DofsPerNode_->resize(DofsPerNode_.size()+1);
        
        NumberOfBlocks_++;
        
        /////
        int blockId = NumberOfBlocks_-1;
        
        // Process the parameter list
        std::stringstream blockIdStringstream;
        blockIdStringstream << blockId+1;
        std::string blockIdString = blockIdStringstream.str();
        Teuchos::RCP<Teuchos::ParameterList> coarseSpaceList = sublist(sublist(this->ParameterList_,"Blocks"),blockIdString.c_str());
        
        bool useForCoarseSpace = coarseSpaceList->get("Use For Coarse Space",true);
        
        IndicesGamma_[blockId] = LOVecPtr(0);
        IndicesGammaDofs_[blockId] = LOVecPtr(0);
        
        if (useForCoarseSpace) {
            //Epetra_SerialComm serialComm;
            MapPtr serialGammaMap = Xpetra::MapFactory<LO,GO,NO>::Build(dofsMap->lib(),dofsMap->getNodeNumElements(),0,this->SerialComm_);
            MVPhiGamma_[blockId] = Xpetra::MultiVectorFactory<LO,GO,NO>::Build(serialGammaMap,dofsMap->getNodeNumElements());
        }
        
        for (int i=0; i<dofsMap->getNodeNumElements(); i++) {
            IndicesGamma_[blockId]->push_back(i);
            IndicesGammaDofs_[blockId]->push_back(i);
            
            if (useForCoarseSpace) {
                MVPhiGamma_[blockId]->replaceLocalValue(i,i,1.0);
            }
        }
        
        IndicesI_[blockId] = LOVecPtr(0);
        IndicesIDofs_[blockId] = LOVecPtr(0);
        
        if (useForCoarseSpace) {
            BlockCoarseMaps_[blockId] = Xpetra::MapFactory<LO,GO,NO>::Build(dofsMap->lib(),-1,IndicesGamma_[blockId](),0,this->MpiComm_);
        }
        
        DofsMaps_[blockId] = MapPtrVecPtr(0);
        DofsMaps_[blockId].push_back(dofsMap);
        
        DofsPerNode_[blockId] = 1;
        
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::computeBasis()
    {
        // Build repeatedMap for the local saddle point problem // Todo:: Eigentlich gehört das in Initialize
        MapPtr repeatedMap = assembleRepeatedMap();
        
        // Build local saddle point problem
        CrsMatrixPtr repeatedMatrix = FROSch::ExtractLocalSubdomainMatrix(this->K_,repeatedMap);
        
        // Extract submatrices
        GOVec indicesGammaDofsAll(0);
        GOVec indicesIDofsAll(0);
        LO tmp = 0;
        
        for (UN i=0; i<NumberOfBlocks_; i++) {
            for (UN j=0; j<IndicesGammaDofs_[i].size(); j++) {
                indicesGammaDofsAll.push_back(tmp+IndicesGammaDofs_[i][j]);
            }
            for (UN j=0; j<IndicesIDofs_[i].size(); j++) {
                indicesIDofsAll.push_back(tmp+IndicesIDofs_[i][j]);
            }
            tmp += IndicesGammaDofs_[i].size()+IndicesIDofs_[i].size(); // Ist das mit tmp korrekt?
        }
        
        CrsMatrixPtr kII;
        CrsMatrixPtr kIGamma;
        CrsMatrixPtr kGammaI;
        CrsMatrixPtr kGammaGamma;
        
        FROSch::BuildSubmatrices(repeatedMatrix,indicesIDofsAll(),kII,kIGamma,kGammaI,kGammaGamma);
        //if (this->Verbose_) std::cout << kII->NumMyRows() << " " << kGammaI->NumMyRows();
        //if (this->Verbose_) std::cout << *kII << *kIGamma << *kGammaI << *kGammaGamma;
        //FROSCH_ASSERT(0!=0,"STOP!");
        
        // Assemble coarse map
        MapPtr coarseMap = assembleCoarseMap(); // Todo:: Eigentlich gehört das in Initialize
        
        // Build the saddle point harmonic extensions
        computeAndFillPhi(repeatedMatrix,repeatedMap,coarseMap,indicesGammaDofsAll(),indicesIDofsAll(),kII,kIGamma);
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    typename GDSWCoarseOperator<SC,LO,GO,NO>::MapPtr GDSWCoarseOperator<SC,LO,GO,NO>::assembleRepeatedMap()
    {
        FROSCH_ASSERT(DofsMaps_.size()==NumberOfBlocks_,"DofsMaps_.size()!=NumberOfBlocks_");
        FROSCH_ASSERT(DofsPerNode_.size()==NumberOfBlocks_,"DofsPerNode_.size()!=NumberOfBlocks_");
        
        GOVec mapVector(0);
        for (UN i=0; i<NumberOfBlocks_; i++) {
            FROSCH_ASSERT(DofsMaps_[i].size()==DofsPerNode_[i],"DofsMaps_[i].size()!=DofsPerNode_[i]");
            UN numMyElements = DofsMaps_[i][0]->getNodeNumElements();
            for (UN j=1; j<DofsPerNode_[i]; j++) {
                FROSCH_ASSERT(DofsMaps_[i][j]->getNodeNumElements()==(int) numMyElements,"DofsMaps_[i][j]->getNodeNumElements()==numMyElements");
            }
            for (UN j=0; j<numMyElements; j++) {
                for (UN k=0; k<DofsPerNode_[i]; k++) {
                    mapVector.push_back(DofsMaps_[i][k]->getGlobalElement(j));
                }
            }
        }
        return Xpetra::MapFactory<LO,GO,NO>::Build(DofsMaps_[0][0]->lib(),-1,mapVector(),0,this->MpiComm_);
    }
    
    template <class SC,class LO,class GO,class NO>
    typename GDSWCoarseOperator<SC,LO,GO,NO>::MapPtr GDSWCoarseOperator<SC,LO,GO,NO>::assembleCoarseMap()
    {
        GOVec mapVector(0);
        GO tmp = 0;
        for (UN i=0; i<NumberOfBlocks_; i++) {
            if (!BlockCoarseMaps_[i].is_null()) {
                for (int j=0; j<BlockCoarseMaps_[i]->getNodeNumElements(); j++) {
                    mapVector.push_back(BlockCoarseMaps_[i]->getGlobalElement(j)+tmp);
                }
                tmp += BlockCoarseMaps_[i]->getMaxAllGlobalIndex()+1;
            }
        }
        return Xpetra::MapFactory<LO,GO,NO>::Build(DofsMaps_[0][0]->lib(),-1,mapVector(),0,this->MpiComm_);
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::phiGammaGDSW(UN blockId,
                                                      bool buildRotations,
                                                      UN dimension,
                                                      UN dofsPerNode,
                                                      SCVecPtr2D &localNodeList,
                                                      LOVecPtr2D &partMappings,
                                                      EntitySetPtr &vertices,
                                                      EntitySetPtr &shortEdges,
                                                      EntitySetPtr &straightEdges,
                                                      EntitySetPtr &edges,
                                                      EntitySetPtr &faces,
                                                      BoolVecPtr &coarseSpaceFunctions)
    {
        if (buildRotations) {
            FROSCH_ASSERT(localNodeList[0].size()==dimension,"dimension of the localNodeList is wrong.");
        }
        
        //Epetra_SerialComm serialComm;
        MapPtr serialGammaMap = Xpetra::MapFactory<LO,GO,NO>::Build(BlockCoarseMaps_[blockId]->lib(),IndicesGammaDofs_[blockId].size(),0,this->SerialComm_);
        MVPhiGamma_[blockId] = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(serialGammaMap,BlockCoarseMaps_[blockId]->getNodeNumElements());
        
        LO ii=0;
        SC x,y,z,rx,ry,rz;
        SCVecPtr dir;
        LO tmp=0;
        
        // vertices
        if (coarseSpaceFunctions[0]) {
            for (UN k=0; k<dofsPerNode; k++) {
                for (UN i=0; i<vertices->getNumEntities(); i++) {
                    MVPhiGamma_[blockId]->replaceLocalValue(vertices->getEntity(i)->getGammaDofID(0,k),partMappings[ii][i],1.0);
                }
                ii++;
            }
        }
        
        
        // Short edges
        if (coarseSpaceFunctions[1]) { // Translations
            for (UN k=0; k<dofsPerNode; k++) {
                for (UN i=0; i<shortEdges->getNumEntities(); i++) {
                    for (UN j=0; j<shortEdges->getEntity(i)->getNumNodes(); j++) {
                        MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,k),partMappings[ii][i],1.0);
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
                        x = localNodeList[shortEdges->getEntity(i)->getLocalNodeID(j)][0];
                        y = localNodeList[shortEdges->getEntity(i)->getLocalNodeID(j)][1];
                        rx = -y;
                        ry = x;
                        MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,0),partMappings[ii][i],rx);
                        MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,1),partMappings[ii][i],ry);
                    }
                }
                ii++;
            } else if (dimension == 3) {
                for (UN i=0; i<shortEdges->getNumEntities(); i++) {
                    for (UN j=0; j<shortEdges->getEntity(i)->getNumNodes(); j++) {
                        // Get the direction of the short edge
                        dir = shortEdges->getDirection(dimension,localNodeList,i);
                        tmp = 0;
                        
                        FROSCH_ASSERT(sqrt(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2])>1.0e-12,"The direction vector is 0. ERROR!");
                        
                        x = localNodeList[shortEdges->getEntity(i)->getLocalNodeID(j)][0];
                        y = localNodeList[shortEdges->getEntity(i)->getLocalNodeID(j)][1];
                        z = localNodeList[shortEdges->getEntity(i)->getLocalNodeID(j)][2];
                        
                        // Rotation 1
                        if ((fabs(dir[0])>1.0e-12) || (fabs(dir[1])>1.0e-12)) {
                            rx = y;
                            ry = -x;
                            rz = 0;
                            MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,0),partMappings[ii+tmp][i],rx);
                            MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,1),partMappings[ii+tmp][i],ry);
                            MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,2),partMappings[ii+tmp][i],rz);
                            tmp++;
                        }
                        
                        // Rotation 2
                        if ((fabs(dir[0])>1.0e-12) || (fabs(dir[2])>1.0e-12)) {
                            rx = -z;
                            ry = 0;
                            rz = x;
                            MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,0),partMappings[ii+tmp][i],rx);
                            MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,1),partMappings[ii+tmp][i],ry);
                            MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,2),partMappings[ii+tmp][i],rz);
                            tmp++;
                        }
                        
                        // Rotation 3
                        if (((fabs(dir[1])>1.0e-12) || (fabs(dir[2])>1.0e-12)) && (tmp<2)) {
                            rx = 0;
                            ry = z;
                            rz = -y;
                            MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,0),partMappings[ii+tmp][i],rx);
                            MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,1),partMappings[ii+tmp][i],ry);
                            MVPhiGamma_[blockId]->replaceLocalValue(shortEdges->getEntity(i)->getGammaDofID(j,2),partMappings[ii+tmp][i],rz);
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
                        MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,k),partMappings[ii][i],1.0);
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
                        x = localNodeList[straightEdges->getEntity(i)->getLocalNodeID(j)][0];
                        y = localNodeList[straightEdges->getEntity(i)->getLocalNodeID(j)][1];
                        rx = -y;
                        ry = x;
                        MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,0),partMappings[ii][i],rx);
                        MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,1),partMappings[ii][i],ry);
                    }
                }
                ii++;
            } else if (dimension == 3) {
                for (UN i=0; i<straightEdges->getNumEntities(); i++) {
                    for (UN j=0; j<straightEdges->getEntity(i)->getNumNodes(); j++) {
                        // Get the direction of the straight edge
                        dir = straightEdges->getDirection(dimension,localNodeList,i);
                        tmp=0;
                        
                        FROSCH_ASSERT(sqrt(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2])>1.0e-12,"The direction vector is 0. ERROR!");
                        
                        x = localNodeList[straightEdges->getEntity(i)->getLocalNodeID(j)][0];
                        y = localNodeList[straightEdges->getEntity(i)->getLocalNodeID(j)][1];
                        z = localNodeList[straightEdges->getEntity(i)->getLocalNodeID(j)][2];
                        
                        // Rotation 1
                        if ((fabs(dir[0])>1.0e-12) || (fabs(dir[1])>1.0e-12)) {
                            rx = y;
                            ry = -x;
                            rz = 0;
                            MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,0),partMappings[ii+tmp][i],rx);
                            MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,1),partMappings[ii+tmp][i],ry);
                            MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,2),partMappings[ii+tmp][i],rz);
                            tmp++;
                        }
                        
                        // Rotation 2
                        if ((fabs(dir[0])>1.0e-12) || (fabs(dir[2])>1.0e-12)) {
                            rx = -z;
                            ry = 0;
                            rz = x;
                            MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,0),partMappings[ii+tmp][i],rx);
                            MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,1),partMappings[ii+tmp][i],ry);
                            MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,2),partMappings[ii+tmp][i],rz);
                            tmp++;
                        }
                        
                        // Rotation 3
                        if (((fabs(dir[1])>1.0e-12) || (fabs(dir[2])>1.0e-12)) && (tmp<2)) {
                            rx = 0;
                            ry = z;
                            rz = -y;
                            MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,0),partMappings[ii+tmp][i],rx);
                            MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,1),partMappings[ii+tmp][i],ry);
                            MVPhiGamma_[blockId]->replaceLocalValue(straightEdges->getEntity(i)->getGammaDofID(j,2),partMappings[ii+tmp][i],rz);
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
                        MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,k),partMappings[ii][i],1.0);
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
                        x = localNodeList[edges->getEntity(i)->getLocalNodeID(j)][0];
                        y = localNodeList[edges->getEntity(i)->getLocalNodeID(j)][1];
                        rx = -y;
                        ry = x;
                        MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,0),partMappings[ii][i],rx);
                        MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,1),partMappings[ii][i],ry);
                    }
                }
                ii++;
            } else if (dimension == 3) {
                for (UN i=0; i<edges->getNumEntities(); i++) {
                    for (UN j=0; j<edges->getEntity(i)->getNumNodes(); j++) {
                        x = localNodeList[edges->getEntity(i)->getLocalNodeID(j)][0];
                        y = localNodeList[edges->getEntity(i)->getLocalNodeID(j)][1];
                        z = localNodeList[edges->getEntity(i)->getLocalNodeID(j)][2];
                        
                        // Rotation 1
                        rx = y;
                        ry = -x;
                        rz = 0;
                        MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,0),partMappings[ii][i],rx);
                        MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,1),partMappings[ii][i],ry);
                        MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,2),partMappings[ii][i],rz);
                        
                        // Rotation 2
                        rx = -z;
                        ry = 0;
                        rz = x;
                        MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,0),partMappings[ii+1][i],rx);
                        MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,1),partMappings[ii+1][i],ry);
                        MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,2),partMappings[ii+1][i],rz);
                        
                        // Rotation 3
                        rx = 0;
                        ry = z;
                        rz = -y;
                        MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,0),partMappings[ii+2][i],rx);
                        MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,1),partMappings[ii+2][i],ry);
                        MVPhiGamma_[blockId]->replaceLocalValue(edges->getEntity(i)->getGammaDofID(j,2),partMappings[ii+2][i],rz);
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
                        MVPhiGamma_[blockId]->replaceLocalValue(faces->getEntity(i)->getGammaDofID(j,k),partMappings[ii][i],1.0);
                    }
                }
                ii++;
            }
        }
        
        if (coarseSpaceFunctions[8]) { // Rotations
            
            FROSCH_ASSERT(dofsPerNode>1,"Dofs<2 => Rotations cannot be built.");
            
            for (UN i=0; i<faces->getNumEntities(); i++) {
                for (UN j=0; j<faces->getEntity(i)->getNumNodes(); j++) {
                    x = localNodeList[faces->getEntity(i)->getLocalNodeID(j)][0];
                    y = localNodeList[faces->getEntity(i)->getLocalNodeID(j)][1];
                    z = localNodeList[faces->getEntity(i)->getLocalNodeID(j)][2];
                    
                    // Rotation 1
                    rx = y;
                    ry = -x;
                    rz = 0;
                    MVPhiGamma_[blockId]->replaceLocalValue(faces->getEntity(i)->getGammaDofID(j,0),partMappings[ii][i],rx);
                    MVPhiGamma_[blockId]->replaceLocalValue(faces->getEntity(i)->getGammaDofID(j,1),partMappings[ii][i],ry);
                    MVPhiGamma_[blockId]->replaceLocalValue(faces->getEntity(i)->getGammaDofID(j,2),partMappings[ii][i],rz);
                    
                    // Rotation 2
                    rx = -z;
                    ry = 0;
                    rz = x;
                    MVPhiGamma_[blockId]->replaceLocalValue(faces->getEntity(i)->getGammaDofID(j,0),partMappings[ii+1][i],rx);
                    MVPhiGamma_[blockId]->replaceLocalValue(faces->getEntity(i)->getGammaDofID(j,1),partMappings[ii+1][i],ry);
                    MVPhiGamma_[blockId]->replaceLocalValue(faces->getEntity(i)->getGammaDofID(j,2),partMappings[ii+1][i],rz);
                    
                    // Rotation 3
                    rx = 0;
                    ry = z;
                    rz = -y;
                    MVPhiGamma_[blockId]->replaceLocalValue(faces->getEntity(i)->getGammaDofID(j,0),partMappings[ii+2][i],rx);
                    MVPhiGamma_[blockId]->replaceLocalValue(faces->getEntity(i)->getGammaDofID(j,1),partMappings[ii+2][i],ry);
                    MVPhiGamma_[blockId]->replaceLocalValue(faces->getEntity(i)->getGammaDofID(j,2),partMappings[ii+2][i],rz);
                }
            }
            ii++;
        }
        //cout << *MVPhiGamma_[blockId];
        return 0;
    }
    
    
    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::computeAndFillPhi(CrsMatrixPtr &repeatedMatrix,
                                                           MapPtr &repeatedMap,
                                                           MapPtr &coarseMap,
                                                           GOVecView indicesGammaDofsAll,
                                                           GOVecView indicesIDofsAll,
                                                           CrsMatrixPtr kII,
                                                           CrsMatrixPtr kIGamma)
    {
        this->Phi_ = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(this->K_->getRangeMap(),coarseMap,coarseMap->getNodeNumElements()); // Nonzeroes abhängig von dim/dofs!!!
        MultiVectorPtr mVtmp = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(kII->getRowMap(),coarseMap->getNodeNumElements());
        MultiVectorPtr mVPhiI = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(kII->getRowMap(),coarseMap->getNodeNumElements());

        //Build mVPhiGamma
        MultiVectorPtr mVPhiGamma = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(kIGamma->getDomainMap(),coarseMap->getNodeNumElements());
        LO jj=0;
        LO kk=0;
        for (UN i=0; i<NumberOfBlocks_; i++) {
            LO j = 0;
            LO k = 0;
            //if (this->Verbose_) std::cout << MVPhiGamma_[i]->MyLength() << std::endl;
            if (!MVPhiGamma_[i].is_null()) {
                for (j=0; j<MVPhiGamma_[i]->getNumVectors(); j++) {
                    for (k=0; k<MVPhiGamma_[i]->getLocalLength(); k++) {
                        //if (this->Verbose_) std::cout << j << " " << k << " " <<  (*(*MVPhiGamma_[i])(j))[k] << std::endl;
                        mVPhiGamma->replaceLocalValue(k+kk,j+jj,MVPhiGamma_[i]->getData(j)[k]);
                    }
                }
            } else { // Das ist für den Fall, dass keine Basisfunktionen für einen Block gebaut werden sollen
                //mVPhiGamma->replaceLocalValue(k+kk,j+jj,1.0);
                k=IndicesGammaDofs_[i].size();
            }
            jj += j;
            kk += k;
        }
        
        LO iD;
        SC valueTmp;
        LOVec indices;
        SCVec values;
        // IST DAS LANGSAM??????? TESTEN!!!!!
        for (LO i=0; i<mVPhiGamma->getLocalLength(); i++) {
            indices.resize(0);
            values.resize(0);
            for (LO j=0; j<mVPhiGamma->getNumVectors(); j++) {
                valueTmp=mVPhiGamma->getData(j)[i];
                if (fabs(valueTmp) > 1.0e-8) {
                    indices.push_back(j);
                    values.push_back(valueTmp);
                }
            }
            
            // CHECK, IF OK?!?!?
            iD = this->K_->getRowMap()->getLocalElement(repeatedMap->getGlobalElement(indicesGammaDofsAll[i]));
            //cout << ProcessorMapRepeatedDofs->GID(IndicesGamma[0][i]) << " " << ID << std::endl;
            
            //ID = IndicesGamma[0][i];
            if (iD!=-1) {
                this->Phi_->insertLocalValues(iD,indices(),values());
            }
            
        }
        
        // Hier Multiplikation kIGamma*PhiGamma
        kIGamma->apply(*mVPhiGamma,*mVtmp);
        
        mVtmp->scale(-1.0);
        
        // Jetzt der solver für kII
        ExtensionSolver_.reset(new SubdomainSolver<SC,LO,GO,NO>(kII,sublist(this->ParameterList_,"ExtensionSolver")));
        // DAS MÜSSEN WIR NOCH ÄNDERN -> initialize, compute, apply...
        ExtensionSolver_->initialize();
        ExtensionSolver_->compute();
        ExtensionSolver_->apply(*mVtmp,*mVPhiI);
        
        // Now we have to insert the values to the global matrix Phi
        // IST DAS LANGSAM??????? TESTEN!!!!!
        for (LO i=0; i<mVPhiI->getLocalLength(); i++) {
            indices.resize(0);
            values.resize(0);
            for (LO j=0; j<mVPhiI->getNumVectors(); j++) {
                valueTmp=mVPhiI->getData(j)[i]; //if (this->Verbose_) std::cout << i << " " << this->K_->getRowMap()->getLocalElement(repeatedMap->getGlobalElement(indicesIDofsAll[i])) << " " << j << " " << valueTmp << std::endl;
                if (fabs(valueTmp) > 1.0e-8) {
                    indices.push_back(j);
                    values.push_back(valueTmp);
                }
            }
            
            // CHECK, IF OK?!?!?
            iD = this->K_->getRowMap()->getLocalElement(repeatedMap->getGlobalElement(indicesIDofsAll[i]));
            //ID = IndicesI[0][i];
            if (iD!=-1) {
                this->Phi_->insertLocalValues(iD,indices(),values());
            }
        }
        this->Phi_->fillComplete(coarseMap,this->Phi_->getRowMap());
        
        return 0;
    }
    
}

#endif
