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
    
    template <class SC,class LO,class GO,class NO>
    RGDSWCoarseOperator<SC,LO,GO,NO>::RGDSWCoarseOperator(CrsMatrixPtr k,
                                                          ParameterListPtr parameterList) :
    GDSWCoarseOperator<SC,LO,GO,NO> (k,parameterList)
    {
        
    }
    
    template <class SC,class LO,class GO,class NO>
    int RGDSWCoarseOperator<SC,LO,GO,NO>::resetCoarseSpaceBlock(UN blockId,
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
        
        bool useForCoarseSpace = coarseSpaceList->get("Use For Coarse Space",false);
        int option = coarseSpaceList->get("Option",1);
        bool useRotations = coarseSpaceList->get("Rotations",true);
        
        if (useRotations && nodeList.is_null()) {
            FROSCH_ASSERT(option==1,"Only option 1 can be constructed without a valid node list.");
            useRotations = false;
            if (this->Verbose_) std::cout << "\nWarning: Rotations cannot be used!\n";
        }
        
        this->DofsMaps_[blockId] = dofsMaps;
        this->DofsPerNode_[blockId] = dofsPerNode;
        
        Teuchos::Array<GO> tmpDirichletBoundaryDofs(dirichletBoundaryDofs()); // Here, we do a copy. Maybe, this is not necessary
        sortunique(tmpDirichletBoundaryDofs);
        
        this->DDInterface_.reset(new DDInterface<SC,LO,GO,NO>(dimension,dofsPerNode,nodesMap));
        this->DDInterface_->resetGlobalDofs(dofsMaps);
        this->DDInterface_->removeDirichletNodes(tmpDirichletBoundaryDofs);
        this->DDInterface_->divideUnconnectedEntities(this->K_);
        
        EntitySetPtr vertices,edges,faces,interface,interior,AncestorVertices,AncestorEdges,AncestorFaces;
        MapPtr AncestorVerticesMap,AncestorEdgesMap,AncestorFacesMap;
        
        interface = this->DDInterface_->getInterface();
        interior = this->DDInterface_->getInterior();
        
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
        
        if (useForCoarseSpace) {
            this->DDInterface_->findAncestors();
            
            ////////////////////////////////
            // Build Processor Map Coarse //
            ////////////////////////////////
            MapPtrVecPtr mapVector(dofsPerNode*3+useRotations*3*(dofsPerNode-1+((dimension==3)&&(dofsPerNode==3))));
            
            vertices = this->DDInterface_->getVertices();
            vertices->buildEntityMap(nodesMap);
            
            edges = this->DDInterface_->getEdges();
            edges->buildEntityMap(nodesMap);
            
            faces = this->DDInterface_->getFaces();
            faces->buildEntityMap(nodesMap);
            
            // HIER MUSS NOCH WAS GEÄNDERT WERDEN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            AncestorVertices = this->DDInterface_->getAncestorVertices();
            AncestorVertices->buildEntityMap(nodesMap);
            
            AncestorEdges = this->DDInterface_->getAncestorEdges();
            AncestorEdges->buildEntityMap(nodesMap);
            
            AncestorFaces = this->DDInterface_->getAncestorFaces();
            AncestorFaces->buildEntityMap(nodesMap);
            
            UN ii=0;
            for (UN i=0; i<dofsPerNode; i++) {
                mapVector[ii] = AncestorVertices->getEntityMap();
                ii++;
            }
            if (useRotations) {
                for (UN i=0; i<dofsPerNode-1+((dimension==3)&&(dofsPerNode==3)); i++) {
                    mapVector[ii] = AncestorVertices->getEntityMap();
                    ii++;
                }
            }
            for (UN i=0; i<dofsPerNode; i++) {
                mapVector[ii] = AncestorEdges->getEntityMap();
                ii++;
            }
            if (useRotations) {
                for (UN i=0; i<dofsPerNode-1+((dimension==3)&&(dofsPerNode==3)); i++) {
                    mapVector[ii] = AncestorEdges->getEntityMap();
                    ii++;
                }
            }
            for (UN i=0; i<dofsPerNode; i++) {
                mapVector[ii] = AncestorFaces->getEntityMap();
                ii++;
            }
            if (useRotations) {
                for (UN i=0; i<dofsPerNode-1+((dimension==3)&&(dofsPerNode==3)); i++) {
                    mapVector[ii] = AncestorFaces->getEntityMap();
                    ii++;
                }
            }
            
            LOVec numEntitiesGlobal(3);
            numEntitiesGlobal[0] = AncestorVertices->getEntityMap()->getMaxAllGlobalIndex();
            if (AncestorVertices->getEntityMap()->lib()==Xpetra::UseEpetra || AncestorVertices->getEntityMap()->getGlobalNumElements()>0) {
                numEntitiesGlobal[0] += 1;
            }
            numEntitiesGlobal[1] = AncestorEdges->getEntityMap()->getMaxAllGlobalIndex();
            if (AncestorEdges->getEntityMap()->lib()==Xpetra::UseEpetra || AncestorEdges->getEntityMap()->getGlobalNumElements()>0) {
                numEntitiesGlobal[1] += 1;
            }
            numEntitiesGlobal[2] = AncestorFaces->getEntityMap()->getMaxAllGlobalIndex();
            if (AncestorFaces->getEntityMap()->lib()==Xpetra::UseEpetra || AncestorFaces->getEntityMap()->getGlobalNumElements()>0) {
                numEntitiesGlobal[2] += 1;
            }
            
            for (UN i=0; i<numEntitiesGlobal.size(); i++) {
                if (numEntitiesGlobal[i]<0) {
                    numEntitiesGlobal[i] = 0;
                }
            }
            
            if (this->MpiComm_->getRank() == 0) {
                std::cout << "\n\
                --------------------------------------------\n\
                # vertices:       --- " << numEntitiesGlobal[0] << "\n\
                # edges:          --- " << numEntitiesGlobal[1] << "\n\
                # faces:          --- " << numEntitiesGlobal[2] << "\n\
                --------------------------------------------\n\
                Coarse space:\n\
                --------------------------------------------\n\
                vertices: translations      --- " << 1 << "\n\
                vertices: rotations         --- " << 1 << "\n\
                --------------------------------------------\n";
            }
            
            LOVecPtr2D partMappings;
            this->BlockCoarseMaps_[blockId] = AssembleMaps(mapVector(),partMappings);
            
            ////////////////////
            // Build PhiGamma //
            ////////////////////
            phiGammaReducedGDSW(blockId,option,useRotations,dimension,dofsPerNode,nodeList,partMappings,vertices,edges,faces);            
        }
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int RGDSWCoarseOperator<SC,LO,GO,NO>::phiGammaReducedGDSW(UN blockId,
                                                              int option,
                                                              bool buildRotations,
                                                              UN dimension,
                                                              UN dofsPerNode,
                                                              MultiVectorPtr nodeList,
                                                              LOVecPtr2D partMappings,
                                                              EntitySetPtr vertices,
                                                              EntitySetPtr edges,
                                                              EntitySetPtr faces)
    {
        if (buildRotations || (option == 3) ) {
            FROSCH_ASSERT(nodeList->getNumVectors()==dimension,"dimension of the nodeList is wrong.");
        }
        
        //Epetra_SerialComm serialComm;
        MapPtr serialGammaMap = Xpetra::MapFactory<LO,GO,NO>::Build(this->BlockCoarseMaps_[blockId]->lib(),this->GammaDofs_[blockId].size(),0,this->SerialComm_);
        this->MVPhiGamma_[blockId] = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(serialGammaMap,this->BlockCoarseMaps_[blockId]->getNodeNumElements());
        
        //int tmp=0;
        
        /*
         Die Schleife ist noch nicht korrekt. Vermutlich muss man zweimal durchgehen: Einmal um herauszufinden welche Werte auf die Kanten und Flächen gesetzt werden müssen
         und beim zweiten Mal, um die Werte zu setzen.
         Außerdem sind im Moment noch zu viel Nullen im Vektor und die Länge des Vektors sollte man überprüfen... Wobei die Länge ist wahrscheinlich korrekt -> 3D
         
         PartMappings[itmp]->at(i) dann steht itmp für vertices (translation1=0, translation2=1, translation3=2, rotation1=3, rotation2=4, rotation3=5), edges (translation1=6, ...), faces (...)
         und i für die Nummer des Vertex, der Edge, etc.
         */
        //vector<double> faceVert(faces->getNumEntities());
        //vector<double> edgeValues(edges->getNumEntities());
        
        LO itmp=0;
        SC x,y,z,rx,ry,rz;
        SC edgeValue;
        SC faceValue;
        
        switch (option) {
            case 1:
            {
                LOVec vertexAncestorsFace(0);
                
                // Vertices translations
                if (dimension==2) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        for (UN i=0; i<edges->getNumEntities(); i++) {
                            InterfaceEntityPtr edge = edges->getEntity(i);
                            EntitySetPtr AncestorVertices = edge->getAncestors();
                            edgeValue = 1.0/SC(AncestorVertices->getNumEntities());
                            for (UN ii=0; ii<AncestorVertices->getNumEntities(); ii++) {
                                InterfaceEntityPtr AncestorVertex = AncestorVertices->getEntity(ii);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,k),partMappings[itmp][AncestorVertex->getAncestorID()],1.0);
                                for (UN j=0; j<edge->getNumNodes(); j++) {
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,k),partMappings[itmp][AncestorVertex->getAncestorID()],edgeValue);
                                }
                            }
                        }
                        itmp++;
                    }
                } else if (dimension==3) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        for (UN i=0; i<faces->getNumEntities(); i++) {
                            vertexAncestorsFace.resize(0);
                            InterfaceEntityPtr face = faces->getEntity(i);
                            EntitySetPtr AncestorEdges = face->getAncestors();
                            for (UN ii=0; ii<AncestorEdges->getNumEntities(); ii++) {
                                InterfaceEntityPtr AncestorEdge = AncestorEdges->getEntity(ii);
                                EntitySetPtr AncestorVertices = AncestorEdge->getAncestors();
                                edgeValue = 1.0/SC(AncestorVertices->getNumEntities());
                                for (UN iii=0; iii<AncestorVertices->getNumEntities(); iii++) {
                                    InterfaceEntityPtr AncestorVertex = AncestorVertices->getEntity(iii);
                                    vertexAncestorsFace.push_back(AncestorVertex->getAncestorID());
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,k),partMappings[itmp][AncestorVertex->getAncestorID()],1.0);
                                    for (UN j=0; j<AncestorEdge->getNumNodes(); j++) {
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,k),partMappings[itmp][AncestorVertex->getAncestorID()],edgeValue);
                                    }
                                }
                            }
                            sortunique(vertexAncestorsFace);
                            faceValue = 1.0/SC(vertexAncestorsFace.size());
                            for (UN ii=0; ii<vertexAncestorsFace.size(); ii++) {
                                for (UN j=0; j<face->getNumNodes(); j++) {
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,k),partMappings[itmp][vertexAncestorsFace[ii]],faceValue);
                                }
                            }
                        }
                        itmp++;
                    }
                }
                
                // Vertices rotations
                if (buildRotations) {
                    if (dimension==2) {
                        for (UN i=0; i<edges->getNumEntities(); i++) {
                            InterfaceEntityPtr edge = edges->getEntity(i);
                            EntitySetPtr AncestorVertices = edge->getAncestors();
                            edgeValue = 1.0/SC(AncestorVertices->getNumEntities());
                            for (UN ii=0; ii<AncestorVertices->getNumEntities(); ii++) {
                                InterfaceEntityPtr AncestorVertex = AncestorVertices->getEntity(ii);
                                
                                x = nodeList->getData(0)[AncestorVertex->getLocalNodeID(0)];
                                y = nodeList->getData(1)[AncestorVertex->getLocalNodeID(0)];
                                rx = -y;
                                ry = x;
                                this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,0),partMappings[itmp][AncestorVertex->getAncestorID()],rx);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,1),partMappings[itmp][AncestorVertex->getAncestorID()],ry);
                                
                                for (UN j=0; j<edge->getNumNodes(); j++) {
                                    x = nodeList->getData(0)[edge->getLocalNodeID(j)];
                                    y = nodeList->getData(1)[edge->getLocalNodeID(j)];
                                    rx = -y;
                                    ry = x;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,0),partMappings[itmp][AncestorVertex->getAncestorID()],edgeValue*rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,1),partMappings[itmp][AncestorVertex->getAncestorID()],edgeValue*ry);
                                    
                                }
                            }
                        }
                        itmp++;
                    } else if (dimension==3) {
                        for (UN i=0; i<faces->getNumEntities(); i++) {
                            vertexAncestorsFace.resize(0);
                            InterfaceEntityPtr face = faces->getEntity(i);
                            EntitySetPtr AncestorEdges = face->getAncestors();
                            for (UN ii=0; ii<AncestorEdges->getNumEntities(); ii++) {
                                InterfaceEntityPtr AncestorEdge = AncestorEdges->getEntity(ii);
                                EntitySetPtr AncestorVertices = AncestorEdge->getAncestors();
                                edgeValue = 1.0/SC(AncestorVertices->getNumEntities());
                                for (UN iii=0; iii<AncestorVertices->getNumEntities(); iii++) {
                                    InterfaceEntityPtr AncestorVertex = AncestorVertices->getEntity(iii);
                                    vertexAncestorsFace.push_back(AncestorVertex->getAncestorID());
                                    
                                    x = nodeList->getData(0)[AncestorVertex->getLocalNodeID(0)];
                                    y = nodeList->getData(1)[AncestorVertex->getLocalNodeID(0)];
                                    z = nodeList->getData(2)[AncestorVertex->getLocalNodeID(0)];
                                    
                                    // Rotation 1
                                    rx = y;
                                    ry = -x;
                                    rz = 0;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,0),partMappings[itmp][AncestorVertex->getAncestorID()],rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,1),partMappings[itmp][AncestorVertex->getAncestorID()],ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,2),partMappings[itmp][AncestorVertex->getAncestorID()],rz);
                                    
                                    // Rotation 2
                                    rx = -z;
                                    ry = 0;
                                    rz = x;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,0),partMappings[itmp+1][AncestorVertex->getAncestorID()],rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,1),partMappings[itmp+1][AncestorVertex->getAncestorID()],ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,2),partMappings[itmp+1][AncestorVertex->getAncestorID()],rz);
                                    
                                    // Rotation 3
                                    rx = 0;
                                    ry = z;
                                    rz = -y;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,0),partMappings[itmp+2][AncestorVertex->getAncestorID()],rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,1),partMappings[itmp+2][AncestorVertex->getAncestorID()],ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,2),partMappings[itmp+2][AncestorVertex->getAncestorID()],rz);
                                    
                                    for (UN j=0; j<AncestorEdge->getNumNodes(); j++) {
                                        x = nodeList->getData(0)[AncestorEdge->getLocalNodeID(j)];
                                        y = nodeList->getData(1)[AncestorEdge->getLocalNodeID(j)];
                                        z = nodeList->getData(2)[AncestorEdge->getLocalNodeID(j)];
                                        
                                        // Rotation 1
                                        rx = y;
                                        ry = -x;
                                        rz = 0;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,0),partMappings[itmp][AncestorVertex->getAncestorID()],edgeValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,1),partMappings[itmp][AncestorVertex->getAncestorID()],edgeValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,2),partMappings[itmp][AncestorVertex->getAncestorID()],edgeValue*rz);
                                        
                                        // Rotation 2
                                        rx = -z;
                                        ry = 0;
                                        rz = x;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,0),partMappings[itmp+1][AncestorVertex->getAncestorID()],edgeValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,1),partMappings[itmp+1][AncestorVertex->getAncestorID()],edgeValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,2),partMappings[itmp+1][AncestorVertex->getAncestorID()],edgeValue*rz);
                                        
                                        // Rotation 3
                                        rx = 0;
                                        ry = z;
                                        rz = -y;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,0),partMappings[itmp+2][AncestorVertex->getAncestorID()],edgeValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,1),partMappings[itmp+2][AncestorVertex->getAncestorID()],edgeValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,2),partMappings[itmp+2][AncestorVertex->getAncestorID()],edgeValue*rz);
                                    }
                                    
                                }
                            }
                            sortunique(vertexAncestorsFace);
                            faceValue = 1.0/SC(vertexAncestorsFace.size());
                            for (UN ii=0; ii<vertexAncestorsFace.size(); ii++) {
                                for (UN j=0; j<face->getNumNodes(); j++) {
                                    x = nodeList->getData(0)[face->getLocalNodeID(j)];
                                    y = nodeList->getData(1)[face->getLocalNodeID(j)];
                                    z = nodeList->getData(2)[face->getLocalNodeID(j)];
                                    
                                    // Rotation 1
                                    rx = y;
                                    ry = -x;
                                    rz = 0;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp][vertexAncestorsFace[ii]],faceValue*rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp][vertexAncestorsFace[ii]],faceValue*ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp][vertexAncestorsFace[ii]],faceValue*rz);
                                    
                                    // Rotation 2
                                    rx = -z;
                                    ry = 0;
                                    rz = x;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp+1][vertexAncestorsFace[ii]],faceValue*rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp+1][vertexAncestorsFace[ii]],faceValue*ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp+1][vertexAncestorsFace[ii]],faceValue*rz);
                                    
                                    // Rotation 3
                                    rx = 0;
                                    ry = z;
                                    rz = -y;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp+2][vertexAncestorsFace[ii]],faceValue*rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp+2][vertexAncestorsFace[ii]],faceValue*ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp+2][vertexAncestorsFace[ii]],faceValue*rz);
                                }
                            }
                        }
                        itmp+=3;
                    }
                }
                
                // Edges translations
                if (dimension==2) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        for (UN i=0; i<edges->getNumEntities(); i++) {
                            InterfaceEntityPtr edge = edges->getEntity(i);
                            EntitySetPtr AncestorVertices = edge->getAncestors();
                            if (AncestorVertices->getNumEntities()==0) {
                                for (UN j=0; j<edge->getNumNodes(); j++) {
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,k),partMappings[itmp][edge->getAncestorID()],1.0);
                                }
                            }
                        }
                        itmp++;
                    }
                } else if (dimension==3) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        for (UN i=0; i<faces->getNumEntities(); i++) {
                            InterfaceEntityPtr face = faces->getEntity(i);
                            EntitySetPtr AncestorEdges = face->getAncestors();
                            faceValue = 1.0/SC(AncestorEdges->getNumEntities());
                            for (UN ii=0; ii<AncestorEdges->getNumEntities(); ii++) {
                                InterfaceEntityPtr AncestorEdge = AncestorEdges->getEntity(ii);
                                EntitySetPtr AncestorVertices = AncestorEdge->getAncestors();
                                if (AncestorVertices->getNumEntities()==0) {
                                    for (UN j=0; j<AncestorEdge->getNumNodes(); j++) {
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,k),partMappings[itmp][AncestorEdge->getAncestorID()],1.0);
                                    }
                                    for (UN j=0; j<face->getNumNodes(); j++) {
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,k),partMappings[itmp][AncestorEdge->getAncestorID()],faceValue);
                                    }
                                }
                            }
                        }
                        itmp++;
                    }
                }
                
                // Edges rotations
                if (buildRotations) {
                    if (dimension==2) {
                        for (UN i=0; i<edges->getNumEntities(); i++) {
                            InterfaceEntityPtr edge = edges->getEntity(i);
                            EntitySetPtr AncestorVertices = edge->getAncestors();
                            if (AncestorVertices->getNumEntities()==0) {
                                for (UN j=0; j<edge->getNumNodes(); j++) {
                                    x = nodeList->getData(0)[edge->getLocalNodeID(j)];
                                    y = nodeList->getData(1)[edge->getLocalNodeID(j)];
                                    rx = -y;
                                    ry = x;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,0),partMappings[itmp][edge->getAncestorID()],rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,1),partMappings[itmp][edge->getAncestorID()],ry);
                                }
                            }
                        }
                        itmp++;
                    } else if (dimension==3) {
                        for (UN i=0; i<faces->getNumEntities(); i++) {
                            InterfaceEntityPtr face = faces->getEntity(i);
                            EntitySetPtr AncestorEdges = face->getAncestors();
                            faceValue = 1.0/SC(AncestorEdges->getNumEntities());
                            for (UN ii=0; ii<AncestorEdges->getNumEntities(); ii++) {
                                InterfaceEntityPtr AncestorEdge = AncestorEdges->getEntity(ii);
                                EntitySetPtr AncestorVertices = AncestorEdge->getAncestors();
                                if (AncestorVertices->getNumEntities()==0) {
                                    for (UN j=0; j<AncestorEdge->getNumNodes(); j++) {
                                        x = nodeList->getData(0)[AncestorEdge->getLocalNodeID(j)];
                                        y = nodeList->getData(1)[AncestorEdge->getLocalNodeID(j)];
                                        z = nodeList->getData(2)[AncestorEdge->getLocalNodeID(j)];
                                        
                                        // Rotation 1
                                        rx = y;
                                        ry = -x;
                                        rz = 0;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,0),partMappings[itmp][AncestorEdge->getAncestorID()],rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,1),partMappings[itmp][AncestorEdge->getAncestorID()],ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,2),partMappings[itmp][AncestorEdge->getAncestorID()],rz);
                                        
                                        // Rotation 2
                                        rx = -z;
                                        ry = 0;
                                        rz = x;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,0),partMappings[itmp+1][AncestorEdge->getAncestorID()],rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,1),partMappings[itmp+1][AncestorEdge->getAncestorID()],ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,2),partMappings[itmp+1][AncestorEdge->getAncestorID()],rz);
                                        
                                        // Rotation 3
                                        rx = 0;
                                        ry = z;
                                        rz = -y;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,0),partMappings[itmp+2][AncestorEdge->getAncestorID()],rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,1),partMappings[itmp+2][AncestorEdge->getAncestorID()],ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,2),partMappings[itmp+2][AncestorEdge->getAncestorID()],rz);
                                    }
                                    for (UN j=0; j<face->getNumNodes(); j++) {
                                        x = nodeList->getData(0)[face->getLocalNodeID(j)];
                                        y = nodeList->getData(1)[face->getLocalNodeID(j)];
                                        z = nodeList->getData(2)[face->getLocalNodeID(j)];
                                        
                                        // Rotation 1
                                        rx = y;
                                        ry = -x;
                                        rz = 0;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp][AncestorEdge->getAncestorID()],faceValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp][AncestorEdge->getAncestorID()],faceValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp][AncestorEdge->getAncestorID()],faceValue*rz);
                                        
                                        // Rotation 2
                                        rx = -z;
                                        ry = 0;
                                        rz = x;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp+1][AncestorEdge->getAncestorID()],faceValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp+1][AncestorEdge->getAncestorID()],faceValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp+1][AncestorEdge->getAncestorID()],faceValue*rz);
                                        
                                        // Rotation 3
                                        rx = 0;
                                        ry = z;
                                        rz = -y;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp+2][AncestorEdge->getAncestorID()],faceValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp+2][AncestorEdge->getAncestorID()],faceValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp+2][AncestorEdge->getAncestorID()],faceValue*rz);
                                    }
                                }
                            }
                        }
                        itmp+=3;
                    }
                }
                
                // Faces translations
                for (UN k=0; k<dofsPerNode; k++) {
                    for (UN i=0; i<faces->getNumEntities(); i++) {
                        InterfaceEntityPtr face = faces->getEntity(i);
                        EntitySetPtr AncestorEdges = face->getAncestors();
                        if (AncestorEdges->getNumEntities()==0) {
                            for (UN j=0; j<face->getNumNodes(); j++) {
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,k),partMappings[itmp][face->getAncestorID()],1.0);
                            }
                        }
                    }
                    itmp++;
                }
                // Faces rotations
                if (buildRotations) {
                    for (UN i=0; i<faces->getNumEntities(); i++) {
                        InterfaceEntityPtr face = faces->getEntity(i);
                        EntitySetPtr AncestorEdges = face->getAncestors();
                        if (AncestorEdges->getNumEntities()==0) {
                            for (UN j=0; j<face->getNumNodes(); j++) {
                                x = nodeList->getData(0)[face->getLocalNodeID(j)];
                                y = nodeList->getData(1)[face->getLocalNodeID(j)];
                                z = nodeList->getData(2)[face->getLocalNodeID(j)];
                                
                                // Rotation 1
                                rx = y;
                                ry = -x;
                                rz = 0;
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp][face->getAncestorID()],rx);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp][face->getAncestorID()],ry);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp][face->getAncestorID()],rz);
                                
                                // Rotation 2
                                rx = -z;
                                ry = 0;
                                rz = x;
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp+1][face->getAncestorID()],rx);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp+1][face->getAncestorID()],ry);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp+1][face->getAncestorID()],rz);
                                
                                // Rotation 3
                                rx = 0;
                                ry = z;
                                rz = -y;
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp+2][face->getAncestorID()],rx);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp+2][face->getAncestorID()],ry);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp+2][face->getAncestorID()],rz);
                            }
                        }
                    }
                    itmp+=3;
                }
                break;
            }
            case 2:
            {
                FROSCH_ASSERT(0!=0,"Only options 1 and 3 are implemented so far...");
                break;
            }
            case 3:
            {
                SCVecPtr edgeValues;
                SCVecPtr faceValues;
                EntitySetPtr vertexAncestorsFace;
                
                // Vertices translations
                if (dimension==2) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        for (UN i=0; i<edges->getNumEntities(); i++) {
                            InterfaceEntityPtr edge = edges->getEntity(i);
                            edgeValues = SCVecPtr(edge->getNumNodes(),0.0);
                            EntitySetPtr AncestorVertices = edge->getAncestors();
                            for (UN ii=0; ii<AncestorVertices->getNumEntities(); ii++) {
                                InterfaceEntityPtr AncestorVertex = AncestorVertices->getEntity(ii);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,k),partMappings[itmp][AncestorVertex->getAncestorID()],1.0);
                                for (UN j=0; j<edge->getNumNodes(); j++) {
                                    // compute distance
                                    x = nodeList->getData(0)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(0)[edge->getLocalNodeID(j)];
                                    y = nodeList->getData(1)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(1)[edge->getLocalNodeID(j)];
                                    edgeValues[j] += 1.0/sqrt(x*x+y*y);
                                }
                            }
                            for (UN j=0; j<edge->getNumNodes(); j++) {
                                for (UN ii=0; ii<AncestorVertices->getNumEntities(); ii++) {
                                    InterfaceEntityPtr AncestorVertex = AncestorVertices->getEntity(ii);
                                    x = nodeList->getData(0)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(0)[edge->getLocalNodeID(j)];
                                    y = nodeList->getData(1)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(1)[edge->getLocalNodeID(j)];
                                    edgeValue = (1.0/sqrt(x*x+y*y))/(edgeValues[j]);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,k),partMappings[itmp][AncestorVertex->getAncestorID()],edgeValue);
                                }
                            }
                        }
                        itmp++;
                    }
                } else if (dimension==3) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        for (UN i=0; i<faces->getNumEntities(); i++) {
                            vertexAncestorsFace.reset(new EntitySet<SC,LO,GO,NO>(VertexType));
                            InterfaceEntityPtr face = faces->getEntity(i);
                            faceValues = SCVecPtr(face->getNumNodes(),0.0);
                            EntitySetPtr AncestorEdges = face->getAncestors();
                            for (UN ii=0; ii<AncestorEdges->getNumEntities(); ii++) {
                                InterfaceEntityPtr AncestorEdge = AncestorEdges->getEntity(ii);
                                edgeValues = SCVecPtr(AncestorEdge->getNumNodes(),0.0);
                                EntitySetPtr AncestorVertices = AncestorEdge->getAncestors();
                                for (UN iii=0; iii<AncestorVertices->getNumEntities(); iii++) {
                                    InterfaceEntityPtr AncestorVertex = AncestorVertices->getEntity(iii);
                                    vertexAncestorsFace->addEntity(AncestorVertex);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,k),partMappings[itmp][AncestorVertex->getAncestorID()],1.0);
                                    for (UN j=0; j<AncestorEdge->getNumNodes(); j++) {
                                        // compute distance
                                        x = nodeList->getData(0)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(0)[AncestorEdge->getLocalNodeID(j)];
                                        y = nodeList->getData(1)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(1)[AncestorEdge->getLocalNodeID(j)];
                                        z = nodeList->getData(2)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(2)[AncestorEdge->getLocalNodeID(j)];
                                        edgeValues[j] += 1.0/sqrt(x*x+y*y+z*z);
                                    }
                                }
                                for (UN j=0; j<AncestorEdge->getNumNodes(); j++) {
                                    for (UN iii=0; iii<AncestorVertices->getNumEntities(); iii++) {
                                        InterfaceEntityPtr AncestorVertex = AncestorVertices->getEntity(iii);
                                        x = nodeList->getData(0)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(0)[AncestorEdge->getLocalNodeID(j)];
                                        y = nodeList->getData(1)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(1)[AncestorEdge->getLocalNodeID(j)];
                                        z = nodeList->getData(2)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(2)[AncestorEdge->getLocalNodeID(j)];
                                        edgeValue = (1.0/sqrt(x*x+y*y+z*z))/(edgeValues[j]);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,k),partMappings[itmp][AncestorVertex->getAncestorID()],edgeValue);
                                    }
                                }
                            }
                            vertexAncestorsFace->sortUnique();
                            for (UN ii=0; ii<vertexAncestorsFace->getNumEntities(); ii++) {
                                InterfaceEntityPtr AncestorVertex = vertexAncestorsFace->getEntity(ii);
                                for (UN j=0; j<face->getNumNodes(); j++) {
                                    // compute distance
                                    x = nodeList->getData(0)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(0)[face->getLocalNodeID(j)];
                                    y = nodeList->getData(1)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(1)[face->getLocalNodeID(j)];
                                    z = nodeList->getData(2)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(2)[face->getLocalNodeID(j)];
                                    faceValues[j] += 1.0/sqrt(x*x+y*y+z*z);
                                }
                            }
                            for (UN j=0; j<face->getNumNodes(); j++) {
                                for (UN ii=0; ii<vertexAncestorsFace->getNumEntities(); ii++) {
                                    InterfaceEntityPtr AncestorVertex = vertexAncestorsFace->getEntity(ii);
                                    x = nodeList->getData(0)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(0)[face->getLocalNodeID(j)];
                                    y = nodeList->getData(1)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(1)[face->getLocalNodeID(j)];
                                    z = nodeList->getData(2)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(2)[face->getLocalNodeID(j)];
                                    faceValue = (1.0/sqrt(x*x+y*y+z*z))/(faceValues[j]);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,k),partMappings[itmp][AncestorVertex->getAncestorID()],faceValue);
                                }
                            }
                        }
                        itmp++;
                    }
                }
                
                // Vertices rotations
                if (buildRotations) {
                    if (dimension==2) {
                        for (UN i=0; i<edges->getNumEntities(); i++) {
                            InterfaceEntityPtr edge = edges->getEntity(i);
                            edgeValues = SCVecPtr(edge->getNumNodes(),0.0);
                            EntitySetPtr AncestorVertices = edge->getAncestors();
                            for (UN ii=0; ii<AncestorVertices->getNumEntities(); ii++) {
                                InterfaceEntityPtr AncestorVertex = AncestorVertices->getEntity(ii);
                                
                                x = nodeList->getData(0)[AncestorVertex->getLocalNodeID(0)];
                                y = nodeList->getData(1)[AncestorVertex->getLocalNodeID(0)];
                                rx = -y;
                                ry = x;
                                this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,0),partMappings[itmp][AncestorVertex->getAncestorID()],rx);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,1),partMappings[itmp][AncestorVertex->getAncestorID()],ry);
                                
                                for (UN j=0; j<edge->getNumNodes(); j++) {
                                    // compute distance
                                    x = nodeList->getData(0)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(0)[edge->getLocalNodeID(j)];
                                    y = nodeList->getData(1)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(1)[edge->getLocalNodeID(j)];
                                    edgeValues[j] += 1.0/sqrt(x*x+y*y);
                                }
                            }
                            for (UN j=0; j<edge->getNumNodes(); j++) {
                                for (UN ii=0; ii<AncestorVertices->getNumEntities(); ii++) {
                                    InterfaceEntityPtr AncestorVertex = AncestorVertices->getEntity(ii);
                                    x = nodeList->getData(0)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(0)[edge->getLocalNodeID(j)];
                                    y = nodeList->getData(1)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(1)[edge->getLocalNodeID(j)];
                                    edgeValue = (1.0/sqrt(x*x+y*y))/(edgeValues[j]);
                                    
                                    x = nodeList->getData(0)[edge->getLocalNodeID(j)];
                                    y = nodeList->getData(1)[edge->getLocalNodeID(j)];
                                    rx = -y;
                                    ry = x;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,0),partMappings[itmp][AncestorVertex->getAncestorID()],edgeValue*rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,1),partMappings[itmp][AncestorVertex->getAncestorID()],edgeValue*ry);
                                }
                            }
                        }
                        itmp++;
                    } else if (dimension==3) {
                        for (UN i=0; i<faces->getNumEntities(); i++) {
                            vertexAncestorsFace.reset(new EntitySet<SC,LO,GO,NO>(VertexType));
                            InterfaceEntityPtr face = faces->getEntity(i);
                            faceValues = SCVecPtr(face->getNumNodes(),0.0);
                            EntitySetPtr AncestorEdges = face->getAncestors();
                            for (UN ii=0; ii<AncestorEdges->getNumEntities(); ii++) {
                                InterfaceEntityPtr AncestorEdge = AncestorEdges->getEntity(ii);
                                edgeValues = SCVecPtr(AncestorEdge->getNumNodes(),0.0);
                                EntitySetPtr AncestorVertices = AncestorEdge->getAncestors();
                                for (UN iii=0; iii<AncestorVertices->getNumEntities(); iii++) {
                                    InterfaceEntityPtr AncestorVertex = AncestorVertices->getEntity(iii);
                                    vertexAncestorsFace->addEntity(AncestorVertex);
                                    
                                    x = nodeList->getData(0)[AncestorVertex->getLocalNodeID(0)];
                                    y = nodeList->getData(1)[AncestorVertex->getLocalNodeID(0)];
                                    z = nodeList->getData(2)[AncestorVertex->getLocalNodeID(0)];
                                    
                                    // Rotation 1
                                    rx = y;
                                    ry = -x;
                                    rz = 0;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,0),partMappings[itmp][AncestorVertex->getAncestorID()],rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,1),partMappings[itmp][AncestorVertex->getAncestorID()],ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,2),partMappings[itmp][AncestorVertex->getAncestorID()],rz);
                                    
                                    // Rotation 2
                                    rx = -z;
                                    ry = 0;
                                    rz = x;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,0),partMappings[itmp+1][AncestorVertex->getAncestorID()],rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,1),partMappings[itmp+1][AncestorVertex->getAncestorID()],ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,2),partMappings[itmp+1][AncestorVertex->getAncestorID()],rz);
                                    
                                    // Rotation 3
                                    rx = 0;
                                    ry = z;
                                    rz = -y;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,0),partMappings[itmp+2][AncestorVertex->getAncestorID()],rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,1),partMappings[itmp+2][AncestorVertex->getAncestorID()],ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorVertex->getGammaDofID(0,2),partMappings[itmp+2][AncestorVertex->getAncestorID()],rz);
                                    
                                    for (UN j=0; j<AncestorEdge->getNumNodes(); j++) {
                                        // compute distance
                                        x = nodeList->getData(0)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(0)[AncestorEdge->getLocalNodeID(j)];
                                        y = nodeList->getData(1)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(1)[AncestorEdge->getLocalNodeID(j)];
                                        z = nodeList->getData(2)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(2)[AncestorEdge->getLocalNodeID(j)];
                                        edgeValues[j] += 1.0/sqrt(x*x+y*y+z*z);
                                    }
                                }
                                for (UN j=0; j<AncestorEdge->getNumNodes(); j++) {
                                    for (UN iii=0; iii<AncestorVertices->getNumEntities(); iii++) {
                                        InterfaceEntityPtr AncestorVertex = AncestorVertices->getEntity(iii);
                                        x = nodeList->getData(0)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(0)[AncestorEdge->getLocalNodeID(j)];
                                        y = nodeList->getData(1)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(1)[AncestorEdge->getLocalNodeID(j)];
                                        z = nodeList->getData(2)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(2)[AncestorEdge->getLocalNodeID(j)];
                                        edgeValue = (1.0/sqrt(x*x+y*y+z*z))/(edgeValues[j]);
                                        
                                        x = nodeList->getData(0)[AncestorEdge->getLocalNodeID(j)];
                                        y = nodeList->getData(1)[AncestorEdge->getLocalNodeID(j)];
                                        z = nodeList->getData(2)[AncestorEdge->getLocalNodeID(j)];
                                        
                                        // Rotation 1
                                        rx = y;
                                        ry = -x;
                                        rz = 0;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,0),partMappings[itmp][AncestorVertex->getAncestorID()],edgeValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,1),partMappings[itmp][AncestorVertex->getAncestorID()],edgeValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,2),partMappings[itmp][AncestorVertex->getAncestorID()],edgeValue*rz);
                                        
                                        // Rotation 2
                                        rx = -z;
                                        ry = 0;
                                        rz = x;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,0),partMappings[itmp+1][AncestorVertex->getAncestorID()],edgeValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,1),partMappings[itmp+1][AncestorVertex->getAncestorID()],edgeValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,2),partMappings[itmp+1][AncestorVertex->getAncestorID()],edgeValue*rz);
                                        
                                        // Rotation 3
                                        rx = 0;
                                        ry = z;
                                        rz = -y;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,0),partMappings[itmp+2][AncestorVertex->getAncestorID()],edgeValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,1),partMappings[itmp+2][AncestorVertex->getAncestorID()],edgeValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,2),partMappings[itmp+2][AncestorVertex->getAncestorID()],edgeValue*rz);
                                    }
                                }
                            }
                            vertexAncestorsFace->sortUnique();
                            for (UN ii=0; ii<vertexAncestorsFace->getNumEntities(); ii++) {
                                InterfaceEntityPtr AncestorVertex = vertexAncestorsFace->getEntity(ii);
                                for (UN j=0; j<face->getNumNodes(); j++) {
                                    // compute distance
                                    x = nodeList->getData(0)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(0)[face->getLocalNodeID(j)];
                                    y = nodeList->getData(1)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(1)[face->getLocalNodeID(j)];
                                    z = nodeList->getData(2)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(2)[face->getLocalNodeID(j)];
                                    faceValues[j] += 1.0/sqrt(x*x+y*y+z*z);
                                }
                            }
                            for (UN j=0; j<face->getNumNodes(); j++) {
                                for (UN ii=0; ii<vertexAncestorsFace->getNumEntities(); ii++) {
                                    InterfaceEntityPtr AncestorVertex = vertexAncestorsFace->getEntity(ii);
                                    x = nodeList->getData(0)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(0)[face->getLocalNodeID(j)];
                                    y = nodeList->getData(1)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(1)[face->getLocalNodeID(j)];
                                    z = nodeList->getData(2)[AncestorVertex->getLocalNodeID(0)] - nodeList->getData(2)[face->getLocalNodeID(j)];
                                    faceValue = (1.0/sqrt(x*x+y*y+z*z))/(faceValues[j]);
                                    
                                    x = nodeList->getData(0)[face->getLocalNodeID(j)];
                                    y = nodeList->getData(1)[face->getLocalNodeID(j)];
                                    z = nodeList->getData(2)[face->getLocalNodeID(j)];
                                    
                                    // Rotation 1
                                    rx = y;
                                    ry = -x;
                                    rz = 0;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp][vertexAncestorsFace->getEntity(ii)->getAncestorID()],faceValue*rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp][vertexAncestorsFace->getEntity(ii)->getAncestorID()],faceValue*ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp][vertexAncestorsFace->getEntity(ii)->getAncestorID()],faceValue*rz);
                                    
                                    // Rotation 2
                                    rx = -z;
                                    ry = 0;
                                    rz = x;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp+1][vertexAncestorsFace->getEntity(ii)->getAncestorID()],faceValue*rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp+1][vertexAncestorsFace->getEntity(ii)->getAncestorID()],faceValue*ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp+1][vertexAncestorsFace->getEntity(ii)->getAncestorID()],faceValue*rz);
                                    
                                    // Rotation 3
                                    rx = 0;
                                    ry = z;
                                    rz = -y;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp+2][vertexAncestorsFace->getEntity(ii)->getAncestorID()],faceValue*rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp+2][vertexAncestorsFace->getEntity(ii)->getAncestorID()],faceValue*ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp+2][vertexAncestorsFace->getEntity(ii)->getAncestorID()],faceValue*rz);
                                }
                            }
                        }
                        itmp+=3;
                    }
                }
                
                // Edges translations
                if (dimension==2) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        for (UN i=0; i<edges->getNumEntities(); i++) {
                            InterfaceEntityPtr edge = edges->getEntity(i);
                            EntitySetPtr AncestorVertices = edge->getAncestors();
                            if (AncestorVertices->getNumEntities()==0) {
                                for (UN j=0; j<edge->getNumNodes(); j++) {
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,k),partMappings[itmp][edge->getAncestorID()],1.0);
                                }
                            }
                        }
                        itmp++;
                    }
                } else if (dimension==3) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        for (UN i=0; i<faces->getNumEntities(); i++) {
                            InterfaceEntityPtr face = faces->getEntity(i);
                            faceValues = SCVecPtr(face->getNumNodes(),0.0);
                            EntitySetPtr AncestorEdges = face->getAncestors();
                            for (UN ii=0; ii<AncestorEdges->getNumEntities(); ii++) {
                                InterfaceEntityPtr AncestorEdge = AncestorEdges->getEntity(ii);
                                EntitySetPtr AncestorVertices = AncestorEdge->getAncestors();
                                if (AncestorVertices->getNumEntities()==0) {
                                    for (UN j=0; j<AncestorEdge->getNumNodes(); j++) {
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,k),partMappings[itmp][AncestorEdge->getAncestorID()],1.0);
                                        for (UN jj=0; jj<face->getNumNodes(); jj++) {
                                            x = nodeList->getData(0)[AncestorEdge->getLocalNodeID(j)] - nodeList->getData(0)[face->getLocalNodeID(jj)];
                                            y = nodeList->getData(1)[AncestorEdge->getLocalNodeID(j)] - nodeList->getData(1)[face->getLocalNodeID(jj)];
                                            z = nodeList->getData(2)[AncestorEdge->getLocalNodeID(j)] - nodeList->getData(2)[face->getLocalNodeID(jj)];
                                            faceValues[jj] += 1.0/sqrt(x*x+y*y+z*z);
                                        }
                                    }
                                }
                            }
                            for (UN ii=0; ii<AncestorEdges->getNumEntities(); ii++) {
                                InterfaceEntityPtr AncestorEdge = AncestorEdges->getEntity(ii);
                                EntitySetPtr AncestorVertices = AncestorEdge->getAncestors();
                                if (AncestorVertices->getNumEntities()==0) {
                                    for (UN jj=0; jj<face->getNumNodes(); jj++) {
                                        faceValue = 0.0;
                                        for (UN j=0; j<AncestorEdge->getNumNodes(); j++) {
                                            x = nodeList->getData(0)[AncestorEdge->getLocalNodeID(j)] - nodeList->getData(0)[face->getLocalNodeID(jj)];
                                            y = nodeList->getData(1)[AncestorEdge->getLocalNodeID(j)] - nodeList->getData(1)[face->getLocalNodeID(jj)];
                                            z = nodeList->getData(2)[AncestorEdge->getLocalNodeID(j)] - nodeList->getData(2)[face->getLocalNodeID(jj)];
                                            faceValue += (1.0/sqrt(x*x+y*y+z*z))/(faceValues[j]);
                                        }
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(jj,k),partMappings[itmp][AncestorEdge->getAncestorID()],faceValue);
                                    }
                                }
                            }
                        }
                        itmp++;
                    }
                }
                
                // Edges rotations
                if (buildRotations) {
                    if (dimension==2) {
                        for (UN i=0; i<edges->getNumEntities(); i++) {
                            InterfaceEntityPtr edge = edges->getEntity(i);
                            EntitySetPtr AncestorVertices = edge->getAncestors();
                            if (AncestorVertices->getNumEntities()==0) {
                                for (UN j=0; j<edge->getNumNodes(); j++) {
                                    x = nodeList->getData(0)[edge->getLocalNodeID(j)];
                                    y = nodeList->getData(1)[edge->getLocalNodeID(j)];
                                    rx = -y;
                                    ry = x;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,0),partMappings[itmp][edge->getAncestorID()],rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,1),partMappings[itmp][edge->getAncestorID()],ry);
                                }
                            }
                        }
                        itmp++;
                    } else if (dimension==3) {
                        for (UN i=0; i<faces->getNumEntities(); i++) {
                            InterfaceEntityPtr face = faces->getEntity(i);
                            faceValues = SCVecPtr(face->getNumNodes(),0.0);
                            EntitySetPtr AncestorEdges = face->getAncestors();
                            for (UN ii=0; ii<AncestorEdges->getNumEntities(); ii++) {
                                InterfaceEntityPtr AncestorEdge = AncestorEdges->getEntity(ii);
                                EntitySetPtr AncestorVertices = AncestorEdge->getAncestors();
                                if (AncestorVertices->getNumEntities()==0) {
                                    for (UN j=0; j<AncestorEdge->getNumNodes(); j++) {
                                        x = nodeList->getData(0)[AncestorEdge->getLocalNodeID(j)];
                                        y = nodeList->getData(1)[AncestorEdge->getLocalNodeID(j)];
                                        z = nodeList->getData(2)[AncestorEdge->getLocalNodeID(j)];
                                        
                                        // Rotation 1
                                        rx = y;
                                        ry = -x;
                                        rz = 0;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,0),partMappings[itmp][AncestorEdge->getAncestorID()],rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,1),partMappings[itmp][AncestorEdge->getAncestorID()],ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,2),partMappings[itmp][AncestorEdge->getAncestorID()],rz);
                                        
                                        // Rotation 2
                                        rx = -z;
                                        ry = 0;
                                        rz = x;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,0),partMappings[itmp+1][AncestorEdge->getAncestorID()],rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,1),partMappings[itmp+1][AncestorEdge->getAncestorID()],ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,2),partMappings[itmp+1][AncestorEdge->getAncestorID()],rz);
                                        
                                        // Rotation 3
                                        rx = 0;
                                        ry = z;
                                        rz = -y;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,0),partMappings[itmp+2][AncestorEdge->getAncestorID()],rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,1),partMappings[itmp+2][AncestorEdge->getAncestorID()],ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(AncestorEdge->getGammaDofID(j,2),partMappings[itmp+2][AncestorEdge->getAncestorID()],rz);
                                        
                                        for (UN jj=0; jj<face->getNumNodes(); jj++) {
                                            x = nodeList->getData(0)[AncestorEdge->getLocalNodeID(j)] - nodeList->getData(0)[face->getLocalNodeID(jj)];
                                            y = nodeList->getData(1)[AncestorEdge->getLocalNodeID(j)] - nodeList->getData(1)[face->getLocalNodeID(jj)];
                                            z = nodeList->getData(2)[AncestorEdge->getLocalNodeID(j)] - nodeList->getData(2)[face->getLocalNodeID(jj)];
                                            faceValues[jj] += 1.0/sqrt(x*x+y*y+z*z);
                                        }
                                    }
                                }
                            }
                            for (UN ii=0; ii<AncestorEdges->getNumEntities(); ii++) {
                                InterfaceEntityPtr AncestorEdge = AncestorEdges->getEntity(ii);
                                EntitySetPtr AncestorVertices = AncestorEdge->getAncestors();
                                if (AncestorVertices->getNumEntities()==0) {
                                    for (UN jj=0; jj<face->getNumNodes(); jj++) {
                                        faceValue = 0.0;
                                        for (UN j=0; j<AncestorEdge->getNumNodes(); j++) {
                                            x = nodeList->getData(0)[AncestorEdge->getLocalNodeID(j)] - nodeList->getData(0)[face->getLocalNodeID(jj)];
                                            y = nodeList->getData(1)[AncestorEdge->getLocalNodeID(j)] - nodeList->getData(1)[face->getLocalNodeID(jj)];
                                            z = nodeList->getData(2)[AncestorEdge->getLocalNodeID(j)] - nodeList->getData(2)[face->getLocalNodeID(jj)];
                                            faceValue += (1.0/sqrt(x*x+y*y+z*z))/(faceValues[j]);
                                        }
                                        x = nodeList->getData(0)[face->getLocalNodeID(jj)];
                                        y = nodeList->getData(1)[face->getLocalNodeID(jj)];
                                        z = nodeList->getData(2)[face->getLocalNodeID(jj)];
                                        
                                        // Rotation 1
                                        rx = y;
                                        ry = -x;
                                        rz = 0;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(jj,0),partMappings[itmp][AncestorEdge->getAncestorID()],faceValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(jj,1),partMappings[itmp][AncestorEdge->getAncestorID()],faceValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(jj,2),partMappings[itmp][AncestorEdge->getAncestorID()],faceValue*rz);
                                        
                                        // Rotation 2
                                        rx = -z;
                                        ry = 0;
                                        rz = x;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(jj,0),partMappings[itmp+1][AncestorEdge->getAncestorID()],faceValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(jj,1),partMappings[itmp+1][AncestorEdge->getAncestorID()],faceValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(jj,2),partMappings[itmp+1][AncestorEdge->getAncestorID()],faceValue*rz);
                                        
                                        // Rotation 3
                                        rx = 0;
                                        ry = z;
                                        rz = -y;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(jj,0),partMappings[itmp+2][AncestorEdge->getAncestorID()],faceValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(jj,1),partMappings[itmp+2][AncestorEdge->getAncestorID()],faceValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(jj,2),partMappings[itmp+2][AncestorEdge->getAncestorID()],faceValue*rz);
                                    }
                                }
                            }
                        }
                        itmp+=3;
                    }
                }
                
                // Faces translations
                for (UN k=0; k<dofsPerNode; k++) {
                    for (UN i=0; i<faces->getNumEntities(); i++) {
                        InterfaceEntityPtr face = faces->getEntity(i);
                        EntitySetPtr AncestorEdges = face->getAncestors();
                        if (AncestorEdges->getNumEntities()==0) {
                            for (UN j=0; j<face->getNumNodes(); j++) {
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,k),partMappings[itmp][face->getAncestorID()],1.0);
                            }
                        }
                    }
                    itmp++;
                }
                
                // Faces rotations
                if (buildRotations) {
                    for (UN i=0; i<faces->getNumEntities(); i++) {
                        InterfaceEntityPtr face = faces->getEntity(i);
                        EntitySetPtr AncestorEdges = face->getAncestors();
                        if (AncestorEdges->getNumEntities()==0) {
                            for (UN j=0; j<face->getNumNodes(); j++) {
                                x = nodeList->getData(0)[face->getLocalNodeID(j)];
                                y = nodeList->getData(1)[face->getLocalNodeID(j)];
                                z = nodeList->getData(2)[face->getLocalNodeID(j)];
                                
                                // Rotation 1
                                rx = y;
                                ry = -x;
                                rz = 0;
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp][face->getAncestorID()],rx);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp][face->getAncestorID()],ry);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp][face->getAncestorID()],rz);
                                
                                // Rotation 2
                                rx = -z;
                                ry = 0;
                                rz = x;
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp+1][face->getAncestorID()],rx);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp+1][face->getAncestorID()],ry);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp+1][face->getAncestorID()],rz);
                                
                                // Rotation 3
                                rx = 0;
                                ry = z;
                                rz = -y;
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp+2][face->getAncestorID()],rx);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp+2][face->getAncestorID()],ry);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp+2][face->getAncestorID()],rz);
                            }
                        }
                    }
                    itmp+=3;
                }
                break;
            }
            default:
            {
                FROSCH_ASSERT(0!=0,"Only options 1 and 3 are implemented so far...");
                break;
            }
        }
        
        
        return 0;
    }
    
}

#endif
