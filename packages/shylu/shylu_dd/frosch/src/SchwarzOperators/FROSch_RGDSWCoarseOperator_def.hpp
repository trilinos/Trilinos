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
                                                                MapPtr &nodesMap,
                                                                MapPtrVecPtr &dofsMaps,
                                                                GOVecPtr &myGlobalDirichletBoundaryDofs,
                                                                SCVecPtr2D &localNodeList)
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
        
        if (useRotations && (localNodeList.size()==0) && this->Verbose_) {
            FROSCH_ASSERT(option==1,"Only option 1 can be constructed without a valid node list.");
            useRotations = false;
            if (this->Verbose_) std::cout << "\nWarning: Rotations cannot be used!\n";
        }
        
        this->DofsMaps_[blockId] = dofsMaps;
        this->DofsPerNode_[blockId] = dofsPerNode;
        
        Teuchos::Array<GO> globalDirichletBoundaryDofs(myGlobalDirichletBoundaryDofs()); // Here, we do a copy. Maybe, this is not necessary
        sortunique(globalDirichletBoundaryDofs);
        
        this->DDInterface_.reset(new DDInterface<SC,LO,GO,NO>(dimension,dofsPerNode,nodesMap));
        this->DDInterface_->resetGlobalDofs(dofsMaps);
        this->DDInterface_->removeDirichletNodes(globalDirichletBoundaryDofs);
        this->DDInterface_->divideUnconnectedEntities(this->K_);
        
        EntitySetPtr vertices,edges,faces,interface,interior,parentVertices,parentEdges,parentFaces;
        MapPtr parentVerticesMap,parentEdgesMap,parentFacesMap;
        
        interface = this->DDInterface_->getInterface();
        interior = this->DDInterface_->getInterior();
        
        this->IndicesGammaDofs_[blockId] = LOVecPtr(dofsPerNode*interface->getEntity(0)->getNumNodes());
        this->IndicesIDofs_[blockId] = LOVecPtr(dofsPerNode*interior->getEntity(0)->getNumNodes());
        for (UN k=0; k<dofsPerNode; k++) {
            for (UN i=0; i<interface->getEntity(0)->getNumNodes(); i++) {
                this->IndicesGammaDofs_[blockId][dofsPerNode*i+k] = interface->getEntity(0)->getLocalDofID(i,k);
            }
            for (UN i=0; i<interior->getEntity(0)->getNumNodes(); i++) {
                this->IndicesIDofs_[blockId][dofsPerNode*i+k] = interior->getEntity(0)->getLocalDofID(i,k);
            }
        }
        
        if (useForCoarseSpace) {
            this->DDInterface_->findParents();
            
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
            parentVertices = this->DDInterface_->getParentVertices();
            parentVertices->buildEntityMap(nodesMap);
            
            parentEdges = this->DDInterface_->getParentEdges();
            parentEdges->buildEntityMap(nodesMap);
            
            parentFaces = this->DDInterface_->getParentFaces();
            parentFaces->buildEntityMap(nodesMap);
            
            UN ii=0;
            for (UN i=0; i<dofsPerNode; i++) {
                mapVector[ii] = parentVertices->getEntityMap();
                ii++;
            }
            if (useRotations) {
                for (UN i=0; i<dofsPerNode-1+((dimension==3)&&(dofsPerNode==3)); i++) {
                    mapVector[ii] = parentVertices->getEntityMap();
                    ii++;
                }
            }
            for (UN i=0; i<dofsPerNode; i++) {
                mapVector[ii] = parentEdges->getEntityMap();
                ii++;
            }
            if (useRotations) {
                for (UN i=0; i<dofsPerNode-1+((dimension==3)&&(dofsPerNode==3)); i++) {
                    mapVector[ii] = parentEdges->getEntityMap();
                    ii++;
                }
            }
            for (UN i=0; i<dofsPerNode; i++) {
                mapVector[ii] = parentFaces->getEntityMap();
                ii++;
            }
            if (useRotations) {
                for (UN i=0; i<dofsPerNode-1+((dimension==3)&&(dofsPerNode==3)); i++) {
                    mapVector[ii] = parentFaces->getEntityMap();
                    ii++;
                }
            }
            
            LOVec numEntitiesGlobal(3);
            numEntitiesGlobal[0] = parentVertices->getEntityMap()->getMaxAllGlobalIndex();
            if (parentVertices->getEntityMap()->lib()==Xpetra::UseEpetra || parentVertices->getEntityMap()->getGlobalNumElements()>0) {
                numEntitiesGlobal[0] += 1;
            }
            numEntitiesGlobal[1] = parentEdges->getEntityMap()->getMaxAllGlobalIndex();
            if (parentEdges->getEntityMap()->lib()==Xpetra::UseEpetra || parentEdges->getEntityMap()->getGlobalNumElements()>0) {
                numEntitiesGlobal[1] += 1;
            }
            numEntitiesGlobal[2] = parentFaces->getEntityMap()->getMaxAllGlobalIndex();
            if (parentFaces->getEntityMap()->lib()==Xpetra::UseEpetra || parentFaces->getEntityMap()->getGlobalNumElements()>0) {
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
            this->BlockCoarseMaps_[blockId] = AssembleMaps(mapVector,partMappings);
            
            ////////////////////
            // Build PhiGamma //
            ////////////////////
            phiGammaReducedGDSW(blockId,option,useRotations,dimension,dofsPerNode,localNodeList,partMappings,vertices,edges,faces);            
        }
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int RGDSWCoarseOperator<SC,LO,GO,NO>::phiGammaReducedGDSW(UN blockId,
                                                              int option,
                                                              bool buildRotations,
                                                              UN dimension,
                                                              UN dofsPerNode,
                                                              SCVecPtr2D &localNodeList,
                                                              LOVecPtr2D &partMappings,
                                                              EntitySetPtr &vertices,
                                                              EntitySetPtr &edges,
                                                              EntitySetPtr &faces)
    {
        if (buildRotations || (option == 3) ) {
            FROSCH_ASSERT(localNodeList[0].size()==dimension,"dimension of the localNodeList is wrong.");
        }
        
        //Epetra_SerialComm serialComm;
        MapPtr serialGammaMap = Xpetra::MapFactory<LO,GO,NO>::Build(this->BlockCoarseMaps_[blockId]->lib(),this->IndicesGammaDofs_[blockId].size(),0,this->SerialComm_);
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
                LOVec vertexParentsFace(0);
                
                // Vertices translations
                if (dimension==2) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        for (UN i=0; i<edges->getNumEntities(); i++) {
                            InterfaceEntityPtr edge = edges->getEntity(i);
                            EntitySetPtr parentVertices = edge->getParents();
                            edgeValue = 1.0/SC(parentVertices->getNumEntities());
                            for (UN ii=0; ii<parentVertices->getNumEntities(); ii++) {
                                InterfaceEntityPtr parentVertex = parentVertices->getEntity(ii);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,k),partMappings[itmp][parentVertex->getParentID()],1.0);
                                for (UN j=0; j<edge->getNumNodes(); j++) {
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,k),partMappings[itmp][parentVertex->getParentID()],edgeValue);
                                }
                            }
                        }
                        itmp++;
                    }
                } else if (dimension==3) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        for (UN i=0; i<faces->getNumEntities(); i++) {
                            vertexParentsFace.resize(0);
                            InterfaceEntityPtr face = faces->getEntity(i);
                            EntitySetPtr parentEdges = face->getParents();
                            for (UN ii=0; ii<parentEdges->getNumEntities(); ii++) {
                                InterfaceEntityPtr parentEdge = parentEdges->getEntity(ii);
                                EntitySetPtr parentVertices = parentEdge->getParents();
                                edgeValue = 1.0/SC(parentVertices->getNumEntities());
                                for (UN iii=0; iii<parentVertices->getNumEntities(); iii++) {
                                    InterfaceEntityPtr parentVertex = parentVertices->getEntity(iii);
                                    vertexParentsFace.push_back(parentVertex->getParentID());
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,k),partMappings[itmp][parentVertex->getParentID()],1.0);
                                    for (UN j=0; j<parentEdge->getNumNodes(); j++) {
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,k),partMappings[itmp][parentVertex->getParentID()],edgeValue);
                                    }
                                }
                            }
                            sortunique(vertexParentsFace);
                            faceValue = 1.0/SC(vertexParentsFace.size());
                            for (UN ii=0; ii<vertexParentsFace.size(); ii++) {
                                for (UN j=0; j<face->getNumNodes(); j++) {
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,k),partMappings[itmp][vertexParentsFace[ii]],faceValue);
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
                            EntitySetPtr parentVertices = edge->getParents();
                            edgeValue = 1.0/SC(parentVertices->getNumEntities());
                            for (UN ii=0; ii<parentVertices->getNumEntities(); ii++) {
                                InterfaceEntityPtr parentVertex = parentVertices->getEntity(ii);
                                
                                x = localNodeList[parentVertex->getLocalNodeID(0)][0];
                                y = localNodeList[parentVertex->getLocalNodeID(0)][1];
                                rx = -y;
                                ry = x;
                                this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,0),partMappings[itmp][parentVertex->getParentID()],rx);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,1),partMappings[itmp][parentVertex->getParentID()],ry);
                                
                                for (UN j=0; j<edge->getNumNodes(); j++) {
                                    x = localNodeList[edge->getLocalNodeID(j)][0];
                                    y = localNodeList[edge->getLocalNodeID(j)][1];
                                    rx = -y;
                                    ry = x;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,0),partMappings[itmp][parentVertex->getParentID()],edgeValue*rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,1),partMappings[itmp][parentVertex->getParentID()],edgeValue*ry);
                                    
                                }
                            }
                        }
                        itmp++;
                    } else if (dimension==3) {
                        for (UN i=0; i<faces->getNumEntities(); i++) {
                            vertexParentsFace.resize(0);
                            InterfaceEntityPtr face = faces->getEntity(i);
                            EntitySetPtr parentEdges = face->getParents();
                            for (UN ii=0; ii<parentEdges->getNumEntities(); ii++) {
                                InterfaceEntityPtr parentEdge = parentEdges->getEntity(ii);
                                EntitySetPtr parentVertices = parentEdge->getParents();
                                edgeValue = 1.0/SC(parentVertices->getNumEntities());
                                for (UN iii=0; iii<parentVertices->getNumEntities(); iii++) {
                                    InterfaceEntityPtr parentVertex = parentVertices->getEntity(iii);
                                    vertexParentsFace.push_back(parentVertex->getParentID());
                                    
                                    x = localNodeList[parentVertex->getLocalNodeID(0)][0];
                                    y = localNodeList[parentVertex->getLocalNodeID(0)][1];
                                    z = localNodeList[parentVertex->getLocalNodeID(0)][2];
                                    
                                    // Rotation 1
                                    rx = y;
                                    ry = -x;
                                    rz = 0;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,0),partMappings[itmp][parentVertex->getParentID()],rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,1),partMappings[itmp][parentVertex->getParentID()],ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,2),partMappings[itmp][parentVertex->getParentID()],rz);
                                    
                                    // Rotation 2
                                    rx = -z;
                                    ry = 0;
                                    rz = x;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,0),partMappings[itmp+1][parentVertex->getParentID()],rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,1),partMappings[itmp+1][parentVertex->getParentID()],ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,2),partMappings[itmp+1][parentVertex->getParentID()],rz);
                                    
                                    // Rotation 3
                                    rx = 0;
                                    ry = z;
                                    rz = -y;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,0),partMappings[itmp+2][parentVertex->getParentID()],rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,1),partMappings[itmp+2][parentVertex->getParentID()],ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,2),partMappings[itmp+2][parentVertex->getParentID()],rz);
                                    
                                    for (UN j=0; j<parentEdge->getNumNodes(); j++) {
                                        x = localNodeList[parentEdge->getLocalNodeID(j)][0];
                                        y = localNodeList[parentEdge->getLocalNodeID(j)][1];
                                        z = localNodeList[parentEdge->getLocalNodeID(j)][2];
                                        
                                        // Rotation 1
                                        rx = y;
                                        ry = -x;
                                        rz = 0;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,0),partMappings[itmp][parentVertex->getParentID()],edgeValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,1),partMappings[itmp][parentVertex->getParentID()],edgeValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,2),partMappings[itmp][parentVertex->getParentID()],edgeValue*rz);
                                        
                                        // Rotation 2
                                        rx = -z;
                                        ry = 0;
                                        rz = x;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,0),partMappings[itmp+1][parentVertex->getParentID()],edgeValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,1),partMappings[itmp+1][parentVertex->getParentID()],edgeValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,2),partMappings[itmp+1][parentVertex->getParentID()],edgeValue*rz);
                                        
                                        // Rotation 3
                                        rx = 0;
                                        ry = z;
                                        rz = -y;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,0),partMappings[itmp+2][parentVertex->getParentID()],edgeValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,1),partMappings[itmp+2][parentVertex->getParentID()],edgeValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,2),partMappings[itmp+2][parentVertex->getParentID()],edgeValue*rz);
                                    }
                                    
                                }
                            }
                            sortunique(vertexParentsFace);
                            faceValue = 1.0/SC(vertexParentsFace.size());
                            for (UN ii=0; ii<vertexParentsFace.size(); ii++) {
                                for (UN j=0; j<face->getNumNodes(); j++) {
                                    x = localNodeList[face->getLocalNodeID(j)][0];
                                    y = localNodeList[face->getLocalNodeID(j)][1];
                                    z = localNodeList[face->getLocalNodeID(j)][2];
                                    
                                    // Rotation 1
                                    rx = y;
                                    ry = -x;
                                    rz = 0;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp][vertexParentsFace[ii]],faceValue*rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp][vertexParentsFace[ii]],faceValue*ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp][vertexParentsFace[ii]],faceValue*rz);
                                    
                                    // Rotation 2
                                    rx = -z;
                                    ry = 0;
                                    rz = x;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp+1][vertexParentsFace[ii]],faceValue*rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp+1][vertexParentsFace[ii]],faceValue*ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp+1][vertexParentsFace[ii]],faceValue*rz);
                                    
                                    // Rotation 3
                                    rx = 0;
                                    ry = z;
                                    rz = -y;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp+2][vertexParentsFace[ii]],faceValue*rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp+2][vertexParentsFace[ii]],faceValue*ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp+2][vertexParentsFace[ii]],faceValue*rz);
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
                            EntitySetPtr parentVertices = edge->getParents();
                            if (parentVertices->getNumEntities()==0) {
                                for (UN j=0; j<edge->getNumNodes(); j++) {
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,k),partMappings[itmp][edge->getParentID()],1.0);
                                }
                            }
                        }
                        itmp++;
                    }
                } else if (dimension==3) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        for (UN i=0; i<faces->getNumEntities(); i++) {
                            InterfaceEntityPtr face = faces->getEntity(i);
                            EntitySetPtr parentEdges = face->getParents();
                            faceValue = 1.0/SC(parentEdges->getNumEntities());
                            for (UN ii=0; ii<parentEdges->getNumEntities(); ii++) {
                                InterfaceEntityPtr parentEdge = parentEdges->getEntity(ii);
                                EntitySetPtr parentVertices = parentEdge->getParents();
                                if (parentVertices->getNumEntities()==0) {
                                    for (UN j=0; j<parentEdge->getNumNodes(); j++) {
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,k),partMappings[itmp][parentEdge->getParentID()],1.0);
                                    }
                                    for (UN j=0; j<face->getNumNodes(); j++) {
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,k),partMappings[itmp][parentEdge->getParentID()],faceValue);
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
                            EntitySetPtr parentVertices = edge->getParents();
                            if (parentVertices->getNumEntities()==0) {
                                for (UN j=0; j<edge->getNumNodes(); j++) {
                                    x = localNodeList[edge->getLocalNodeID(j)][0];
                                    y = localNodeList[edge->getLocalNodeID(j)][1];
                                    rx = -y;
                                    ry = x;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,0),partMappings[itmp][edge->getParentID()],rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,1),partMappings[itmp][edge->getParentID()],ry);
                                }
                            }
                        }
                        itmp++;
                    } else if (dimension==3) {
                        for (UN i=0; i<faces->getNumEntities(); i++) {
                            InterfaceEntityPtr face = faces->getEntity(i);
                            EntitySetPtr parentEdges = face->getParents();
                            faceValue = 1.0/SC(parentEdges->getNumEntities());
                            for (UN ii=0; ii<parentEdges->getNumEntities(); ii++) {
                                InterfaceEntityPtr parentEdge = parentEdges->getEntity(ii);
                                EntitySetPtr parentVertices = parentEdge->getParents();
                                if (parentVertices->getNumEntities()==0) {
                                    for (UN j=0; j<parentEdge->getNumNodes(); j++) {
                                        x = localNodeList[parentEdge->getLocalNodeID(j)][0];
                                        y = localNodeList[parentEdge->getLocalNodeID(j)][1];
                                        z = localNodeList[parentEdge->getLocalNodeID(j)][2];
                                        
                                        // Rotation 1
                                        rx = y;
                                        ry = -x;
                                        rz = 0;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,0),partMappings[itmp][parentEdge->getParentID()],rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,1),partMappings[itmp][parentEdge->getParentID()],ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,2),partMappings[itmp][parentEdge->getParentID()],rz);
                                        
                                        // Rotation 2
                                        rx = -z;
                                        ry = 0;
                                        rz = x;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,0),partMappings[itmp+1][parentEdge->getParentID()],rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,1),partMappings[itmp+1][parentEdge->getParentID()],ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,2),partMappings[itmp+1][parentEdge->getParentID()],rz);
                                        
                                        // Rotation 3
                                        rx = 0;
                                        ry = z;
                                        rz = -y;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,0),partMappings[itmp+2][parentEdge->getParentID()],rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,1),partMappings[itmp+2][parentEdge->getParentID()],ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,2),partMappings[itmp+2][parentEdge->getParentID()],rz);
                                    }
                                    for (UN j=0; j<face->getNumNodes(); j++) {
                                        x = localNodeList[face->getLocalNodeID(j)][0];
                                        y = localNodeList[face->getLocalNodeID(j)][1];
                                        z = localNodeList[face->getLocalNodeID(j)][2];
                                        
                                        // Rotation 1
                                        rx = y;
                                        ry = -x;
                                        rz = 0;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp][parentEdge->getParentID()],faceValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp][parentEdge->getParentID()],faceValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp][parentEdge->getParentID()],faceValue*rz);
                                        
                                        // Rotation 2
                                        rx = -z;
                                        ry = 0;
                                        rz = x;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp+1][parentEdge->getParentID()],faceValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp+1][parentEdge->getParentID()],faceValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp+1][parentEdge->getParentID()],faceValue*rz);
                                        
                                        // Rotation 3
                                        rx = 0;
                                        ry = z;
                                        rz = -y;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp+2][parentEdge->getParentID()],faceValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp+2][parentEdge->getParentID()],faceValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp+2][parentEdge->getParentID()],faceValue*rz);
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
                        EntitySetPtr parentEdges = face->getParents();
                        if (parentEdges->getNumEntities()==0) {
                            for (UN j=0; j<face->getNumNodes(); j++) {
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,k),partMappings[itmp][face->getParentID()],1.0);
                            }
                        }
                    }
                    itmp++;
                }
                // Faces rotations
                if (buildRotations) {
                    for (UN i=0; i<faces->getNumEntities(); i++) {
                        InterfaceEntityPtr face = faces->getEntity(i);
                        EntitySetPtr parentEdges = face->getParents();
                        if (parentEdges->getNumEntities()==0) {
                            for (UN j=0; j<face->getNumNodes(); j++) {
                                x = localNodeList[face->getLocalNodeID(j)][0];
                                y = localNodeList[face->getLocalNodeID(j)][1];
                                z = localNodeList[face->getLocalNodeID(j)][2];
                                
                                // Rotation 1
                                rx = y;
                                ry = -x;
                                rz = 0;
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp][face->getParentID()],rx);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp][face->getParentID()],ry);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp][face->getParentID()],rz);
                                
                                // Rotation 2
                                rx = -z;
                                ry = 0;
                                rz = x;
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp+1][face->getParentID()],rx);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp+1][face->getParentID()],ry);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp+1][face->getParentID()],rz);
                                
                                // Rotation 3
                                rx = 0;
                                ry = z;
                                rz = -y;
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp+2][face->getParentID()],rx);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp+2][face->getParentID()],ry);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp+2][face->getParentID()],rz);
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
                EntitySetPtr vertexParentsFace;
                
                // Vertices translations
                if (dimension==2) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        for (UN i=0; i<edges->getNumEntities(); i++) {
                            InterfaceEntityPtr edge = edges->getEntity(i);
                            edgeValues = SCVecPtr(edge->getNumNodes(),0.0);
                            EntitySetPtr parentVertices = edge->getParents();
                            for (UN ii=0; ii<parentVertices->getNumEntities(); ii++) {
                                InterfaceEntityPtr parentVertex = parentVertices->getEntity(ii);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,k),partMappings[itmp][parentVertex->getParentID()],1.0);
                                for (UN j=0; j<edge->getNumNodes(); j++) {
                                    // compute distance
                                    x = localNodeList[parentVertex->getLocalNodeID(0)][0] - localNodeList[edge->getLocalNodeID(j)][0];
                                    y = localNodeList[parentVertex->getLocalNodeID(0)][1] - localNodeList[edge->getLocalNodeID(j)][1];
                                    edgeValues[j] += 1.0/sqrt(x*x+y*y);
                                }
                            }
                            for (UN j=0; j<edge->getNumNodes(); j++) {
                                for (UN ii=0; ii<parentVertices->getNumEntities(); ii++) {
                                    InterfaceEntityPtr parentVertex = parentVertices->getEntity(ii);
                                    x = localNodeList[parentVertex->getLocalNodeID(0)][0] - localNodeList[edge->getLocalNodeID(j)][0];
                                    y = localNodeList[parentVertex->getLocalNodeID(0)][1] - localNodeList[edge->getLocalNodeID(j)][1];
                                    edgeValue = (1.0/sqrt(x*x+y*y))/(edgeValues[j]);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,k),partMappings[itmp][parentVertex->getParentID()],edgeValue);
                                }
                            }
                        }
                        itmp++;
                    }
                } else if (dimension==3) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        for (UN i=0; i<faces->getNumEntities(); i++) {
                            vertexParentsFace.reset(new EntitySet<SC,LO,GO,NO>(VertexType));
                            InterfaceEntityPtr face = faces->getEntity(i);
                            faceValues = SCVecPtr(face->getNumNodes(),0.0);
                            EntitySetPtr parentEdges = face->getParents();
                            for (UN ii=0; ii<parentEdges->getNumEntities(); ii++) {
                                InterfaceEntityPtr parentEdge = parentEdges->getEntity(ii);
                                edgeValues = SCVecPtr(parentEdge->getNumNodes(),0.0);
                                EntitySetPtr parentVertices = parentEdge->getParents();
                                for (UN iii=0; iii<parentVertices->getNumEntities(); iii++) {
                                    InterfaceEntityPtr parentVertex = parentVertices->getEntity(iii);
                                    vertexParentsFace->addEntity(parentVertex);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,k),partMappings[itmp][parentVertex->getParentID()],1.0);
                                    for (UN j=0; j<parentEdge->getNumNodes(); j++) {
                                        // compute distance
                                        x = localNodeList[parentVertex->getLocalNodeID(0)][0] - localNodeList[parentEdge->getLocalNodeID(j)][0];
                                        y = localNodeList[parentVertex->getLocalNodeID(0)][1] - localNodeList[parentEdge->getLocalNodeID(j)][1];
                                        z = localNodeList[parentVertex->getLocalNodeID(0)][2] - localNodeList[parentEdge->getLocalNodeID(j)][2];
                                        edgeValues[j] += 1.0/sqrt(x*x+y*y+z*z);
                                    }
                                }
                                for (UN j=0; j<parentEdge->getNumNodes(); j++) {
                                    for (UN iii=0; iii<parentVertices->getNumEntities(); iii++) {
                                        InterfaceEntityPtr parentVertex = parentVertices->getEntity(iii);
                                        x = localNodeList[parentVertex->getLocalNodeID(0)][0] - localNodeList[parentEdge->getLocalNodeID(j)][0];
                                        y = localNodeList[parentVertex->getLocalNodeID(0)][1] - localNodeList[parentEdge->getLocalNodeID(j)][1];
                                        z = localNodeList[parentVertex->getLocalNodeID(0)][2] - localNodeList[parentEdge->getLocalNodeID(j)][2];
                                        edgeValue = (1.0/sqrt(x*x+y*y+z*z))/(edgeValues[j]);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,k),partMappings[itmp][parentVertex->getParentID()],edgeValue);
                                    }
                                }
                            }
                            vertexParentsFace->sortUnique();
                            for (UN ii=0; ii<vertexParentsFace->getNumEntities(); ii++) {
                                InterfaceEntityPtr parentVertex = vertexParentsFace->getEntity(ii);
                                for (UN j=0; j<face->getNumNodes(); j++) {
                                    // compute distance
                                    x = localNodeList[parentVertex->getLocalNodeID(0)][0] - localNodeList[face->getLocalNodeID(j)][0];
                                    y = localNodeList[parentVertex->getLocalNodeID(0)][1] - localNodeList[face->getLocalNodeID(j)][1];
                                    z = localNodeList[parentVertex->getLocalNodeID(0)][2] - localNodeList[face->getLocalNodeID(j)][2];
                                    faceValues[j] += 1.0/sqrt(x*x+y*y+z*z);
                                }
                            }
                            for (UN j=0; j<face->getNumNodes(); j++) {
                                for (UN ii=0; ii<vertexParentsFace->getNumEntities(); ii++) {
                                    InterfaceEntityPtr parentVertex = vertexParentsFace->getEntity(ii);
                                    x = localNodeList[parentVertex->getLocalNodeID(0)][0] - localNodeList[face->getLocalNodeID(j)][0];
                                    y = localNodeList[parentVertex->getLocalNodeID(0)][1] - localNodeList[face->getLocalNodeID(j)][1];
                                    z = localNodeList[parentVertex->getLocalNodeID(0)][2] - localNodeList[face->getLocalNodeID(j)][2];
                                    faceValue = (1.0/sqrt(x*x+y*y+z*z))/(faceValues[j]);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,k),partMappings[itmp][parentVertex->getParentID()],faceValue);
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
                            EntitySetPtr parentVertices = edge->getParents();
                            for (UN ii=0; ii<parentVertices->getNumEntities(); ii++) {
                                InterfaceEntityPtr parentVertex = parentVertices->getEntity(ii);
                                
                                x = localNodeList[parentVertex->getLocalNodeID(0)][0];
                                y = localNodeList[parentVertex->getLocalNodeID(0)][1];
                                rx = -y;
                                ry = x;
                                this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,0),partMappings[itmp][parentVertex->getParentID()],rx);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,1),partMappings[itmp][parentVertex->getParentID()],ry);
                                
                                for (UN j=0; j<edge->getNumNodes(); j++) {
                                    // compute distance
                                    x = localNodeList[parentVertex->getLocalNodeID(0)][0] - localNodeList[edge->getLocalNodeID(j)][0];
                                    y = localNodeList[parentVertex->getLocalNodeID(0)][1] - localNodeList[edge->getLocalNodeID(j)][1];
                                    edgeValues[j] += 1.0/sqrt(x*x+y*y);
                                }
                            }
                            for (UN j=0; j<edge->getNumNodes(); j++) {
                                for (UN ii=0; ii<parentVertices->getNumEntities(); ii++) {
                                    InterfaceEntityPtr parentVertex = parentVertices->getEntity(ii);
                                    x = localNodeList[parentVertex->getLocalNodeID(0)][0] - localNodeList[edge->getLocalNodeID(j)][0];
                                    y = localNodeList[parentVertex->getLocalNodeID(0)][1] - localNodeList[edge->getLocalNodeID(j)][1];
                                    edgeValue = (1.0/sqrt(x*x+y*y))/(edgeValues[j]);
                                    
                                    x = localNodeList[edge->getLocalNodeID(j)][0];
                                    y = localNodeList[edge->getLocalNodeID(j)][1];
                                    rx = -y;
                                    ry = x;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,0),partMappings[itmp][parentVertex->getParentID()],edgeValue*rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,1),partMappings[itmp][parentVertex->getParentID()],edgeValue*ry);
                                }
                            }
                        }
                        itmp++;
                    } else if (dimension==3) {
                        for (UN i=0; i<faces->getNumEntities(); i++) {
                            vertexParentsFace.reset(new EntitySet<SC,LO,GO,NO>(VertexType));
                            InterfaceEntityPtr face = faces->getEntity(i);
                            faceValues = SCVecPtr(face->getNumNodes(),0.0);
                            EntitySetPtr parentEdges = face->getParents();
                            for (UN ii=0; ii<parentEdges->getNumEntities(); ii++) {
                                InterfaceEntityPtr parentEdge = parentEdges->getEntity(ii);
                                edgeValues = SCVecPtr(parentEdge->getNumNodes(),0.0);
                                EntitySetPtr parentVertices = parentEdge->getParents();
                                for (UN iii=0; iii<parentVertices->getNumEntities(); iii++) {
                                    InterfaceEntityPtr parentVertex = parentVertices->getEntity(iii);
                                    vertexParentsFace->addEntity(parentVertex);
                                    
                                    x = localNodeList[parentVertex->getLocalNodeID(0)][0];
                                    y = localNodeList[parentVertex->getLocalNodeID(0)][1];
                                    z = localNodeList[parentVertex->getLocalNodeID(0)][2];
                                    
                                    // Rotation 1
                                    rx = y;
                                    ry = -x;
                                    rz = 0;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,0),partMappings[itmp][parentVertex->getParentID()],rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,1),partMappings[itmp][parentVertex->getParentID()],ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,2),partMappings[itmp][parentVertex->getParentID()],rz);
                                    
                                    // Rotation 2
                                    rx = -z;
                                    ry = 0;
                                    rz = x;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,0),partMappings[itmp+1][parentVertex->getParentID()],rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,1),partMappings[itmp+1][parentVertex->getParentID()],ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,2),partMappings[itmp+1][parentVertex->getParentID()],rz);
                                    
                                    // Rotation 3
                                    rx = 0;
                                    ry = z;
                                    rz = -y;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,0),partMappings[itmp+2][parentVertex->getParentID()],rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,1),partMappings[itmp+2][parentVertex->getParentID()],ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(parentVertex->getGammaDofID(0,2),partMappings[itmp+2][parentVertex->getParentID()],rz);
                                    
                                    for (UN j=0; j<parentEdge->getNumNodes(); j++) {
                                        // compute distance
                                        x = localNodeList[parentVertex->getLocalNodeID(0)][0] - localNodeList[parentEdge->getLocalNodeID(j)][0];
                                        y = localNodeList[parentVertex->getLocalNodeID(0)][1] - localNodeList[parentEdge->getLocalNodeID(j)][1];
                                        z = localNodeList[parentVertex->getLocalNodeID(0)][2] - localNodeList[parentEdge->getLocalNodeID(j)][2];
                                        edgeValues[j] += 1.0/sqrt(x*x+y*y+z*z);
                                    }
                                }
                                for (UN j=0; j<parentEdge->getNumNodes(); j++) {
                                    for (UN iii=0; iii<parentVertices->getNumEntities(); iii++) {
                                        InterfaceEntityPtr parentVertex = parentVertices->getEntity(iii);
                                        x = localNodeList[parentVertex->getLocalNodeID(0)][0] - localNodeList[parentEdge->getLocalNodeID(j)][0];
                                        y = localNodeList[parentVertex->getLocalNodeID(0)][1] - localNodeList[parentEdge->getLocalNodeID(j)][1];
                                        z = localNodeList[parentVertex->getLocalNodeID(0)][2] - localNodeList[parentEdge->getLocalNodeID(j)][2];
                                        edgeValue = (1.0/sqrt(x*x+y*y+z*z))/(edgeValues[j]);
                                        
                                        x = localNodeList[parentEdge->getLocalNodeID(j)][0];
                                        y = localNodeList[parentEdge->getLocalNodeID(j)][1];
                                        z = localNodeList[parentEdge->getLocalNodeID(j)][2];
                                        
                                        // Rotation 1
                                        rx = y;
                                        ry = -x;
                                        rz = 0;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,0),partMappings[itmp][parentVertex->getParentID()],edgeValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,1),partMappings[itmp][parentVertex->getParentID()],edgeValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,2),partMappings[itmp][parentVertex->getParentID()],edgeValue*rz);
                                        
                                        // Rotation 2
                                        rx = -z;
                                        ry = 0;
                                        rz = x;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,0),partMappings[itmp+1][parentVertex->getParentID()],edgeValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,1),partMappings[itmp+1][parentVertex->getParentID()],edgeValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,2),partMappings[itmp+1][parentVertex->getParentID()],edgeValue*rz);
                                        
                                        // Rotation 3
                                        rx = 0;
                                        ry = z;
                                        rz = -y;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,0),partMappings[itmp+2][parentVertex->getParentID()],edgeValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,1),partMappings[itmp+2][parentVertex->getParentID()],edgeValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,2),partMappings[itmp+2][parentVertex->getParentID()],edgeValue*rz);
                                    }
                                }
                            }
                            vertexParentsFace->sortUnique();
                            for (UN ii=0; ii<vertexParentsFace->getNumEntities(); ii++) {
                                InterfaceEntityPtr parentVertex = vertexParentsFace->getEntity(ii);
                                for (UN j=0; j<face->getNumNodes(); j++) {
                                    // compute distance
                                    x = localNodeList[parentVertex->getLocalNodeID(0)][0] - localNodeList[face->getLocalNodeID(j)][0];
                                    y = localNodeList[parentVertex->getLocalNodeID(0)][1] - localNodeList[face->getLocalNodeID(j)][1];
                                    z = localNodeList[parentVertex->getLocalNodeID(0)][2] - localNodeList[face->getLocalNodeID(j)][2];
                                    faceValues[j] += 1.0/sqrt(x*x+y*y+z*z);
                                }
                            }
                            for (UN j=0; j<face->getNumNodes(); j++) {
                                for (UN ii=0; ii<vertexParentsFace->getNumEntities(); ii++) {
                                    InterfaceEntityPtr parentVertex = vertexParentsFace->getEntity(ii);
                                    x = localNodeList[parentVertex->getLocalNodeID(0)][0] - localNodeList[face->getLocalNodeID(j)][0];
                                    y = localNodeList[parentVertex->getLocalNodeID(0)][1] - localNodeList[face->getLocalNodeID(j)][1];
                                    z = localNodeList[parentVertex->getLocalNodeID(0)][2] - localNodeList[face->getLocalNodeID(j)][2];
                                    faceValue = (1.0/sqrt(x*x+y*y+z*z))/(faceValues[j]);
                                    
                                    x = localNodeList[face->getLocalNodeID(j)][0];
                                    y = localNodeList[face->getLocalNodeID(j)][1];
                                    z = localNodeList[face->getLocalNodeID(j)][2];
                                    
                                    // Rotation 1
                                    rx = y;
                                    ry = -x;
                                    rz = 0;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp][vertexParentsFace->getEntity(ii)->getParentID()],faceValue*rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp][vertexParentsFace->getEntity(ii)->getParentID()],faceValue*ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp][vertexParentsFace->getEntity(ii)->getParentID()],faceValue*rz);
                                    
                                    // Rotation 2
                                    rx = -z;
                                    ry = 0;
                                    rz = x;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp+1][vertexParentsFace->getEntity(ii)->getParentID()],faceValue*rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp+1][vertexParentsFace->getEntity(ii)->getParentID()],faceValue*ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp+1][vertexParentsFace->getEntity(ii)->getParentID()],faceValue*rz);
                                    
                                    // Rotation 3
                                    rx = 0;
                                    ry = z;
                                    rz = -y;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp+2][vertexParentsFace->getEntity(ii)->getParentID()],faceValue*rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp+2][vertexParentsFace->getEntity(ii)->getParentID()],faceValue*ry);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp+2][vertexParentsFace->getEntity(ii)->getParentID()],faceValue*rz);
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
                            EntitySetPtr parentVertices = edge->getParents();
                            if (parentVertices->getNumEntities()==0) {
                                for (UN j=0; j<edge->getNumNodes(); j++) {
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,k),partMappings[itmp][edge->getParentID()],1.0);
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
                            EntitySetPtr parentEdges = face->getParents();
                            for (UN ii=0; ii<parentEdges->getNumEntities(); ii++) {
                                InterfaceEntityPtr parentEdge = parentEdges->getEntity(ii);
                                EntitySetPtr parentVertices = parentEdge->getParents();
                                if (parentVertices->getNumEntities()==0) {
                                    for (UN j=0; j<parentEdge->getNumNodes(); j++) {
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,k),partMappings[itmp][parentEdge->getParentID()],1.0);
                                        for (UN jj=0; jj<face->getNumNodes(); jj++) {
                                            x = localNodeList[parentEdge->getLocalNodeID(j)][0] - localNodeList[face->getLocalNodeID(jj)][0];
                                            y = localNodeList[parentEdge->getLocalNodeID(j)][1] - localNodeList[face->getLocalNodeID(jj)][1];
                                            z = localNodeList[parentEdge->getLocalNodeID(j)][2] - localNodeList[face->getLocalNodeID(jj)][2];
                                            faceValues[jj] += 1.0/sqrt(x*x+y*y+z*z);
                                        }
                                    }
                                }
                            }
                            for (UN ii=0; ii<parentEdges->getNumEntities(); ii++) {
                                InterfaceEntityPtr parentEdge = parentEdges->getEntity(ii);
                                EntitySetPtr parentVertices = parentEdge->getParents();
                                if (parentVertices->getNumEntities()==0) {
                                    for (UN jj=0; jj<face->getNumNodes(); jj++) {
                                        faceValue = 0.0;
                                        for (UN j=0; j<parentEdge->getNumNodes(); j++) {
                                            x = localNodeList[parentEdge->getLocalNodeID(j)][0] - localNodeList[face->getLocalNodeID(jj)][0];
                                            y = localNodeList[parentEdge->getLocalNodeID(j)][1] - localNodeList[face->getLocalNodeID(jj)][1];
                                            z = localNodeList[parentEdge->getLocalNodeID(j)][2] - localNodeList[face->getLocalNodeID(jj)][2];
                                            faceValue += (1.0/sqrt(x*x+y*y+z*z))/(faceValues[j]);
                                        }
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(jj,k),partMappings[itmp][parentEdge->getParentID()],faceValue);
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
                            EntitySetPtr parentVertices = edge->getParents();
                            if (parentVertices->getNumEntities()==0) {
                                for (UN j=0; j<edge->getNumNodes(); j++) {
                                    x = localNodeList[edge->getLocalNodeID(j)][0];
                                    y = localNodeList[edge->getLocalNodeID(j)][1];
                                    rx = -y;
                                    ry = x;
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,0),partMappings[itmp][edge->getParentID()],rx);
                                    this->MVPhiGamma_[blockId]->replaceLocalValue(edge->getGammaDofID(j,1),partMappings[itmp][edge->getParentID()],ry);
                                }
                            }
                        }
                        itmp++;
                    } else if (dimension==3) {
                        for (UN i=0; i<faces->getNumEntities(); i++) {
                            InterfaceEntityPtr face = faces->getEntity(i);
                            faceValues = SCVecPtr(face->getNumNodes(),0.0);
                            EntitySetPtr parentEdges = face->getParents();
                            for (UN ii=0; ii<parentEdges->getNumEntities(); ii++) {
                                InterfaceEntityPtr parentEdge = parentEdges->getEntity(ii);
                                EntitySetPtr parentVertices = parentEdge->getParents();
                                if (parentVertices->getNumEntities()==0) {
                                    for (UN j=0; j<parentEdge->getNumNodes(); j++) {
                                        x = localNodeList[parentEdge->getLocalNodeID(j)][0];
                                        y = localNodeList[parentEdge->getLocalNodeID(j)][1];
                                        z = localNodeList[parentEdge->getLocalNodeID(j)][2];
                                        
                                        // Rotation 1
                                        rx = y;
                                        ry = -x;
                                        rz = 0;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,0),partMappings[itmp][parentEdge->getParentID()],rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,1),partMappings[itmp][parentEdge->getParentID()],ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,2),partMappings[itmp][parentEdge->getParentID()],rz);
                                        
                                        // Rotation 2
                                        rx = -z;
                                        ry = 0;
                                        rz = x;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,0),partMappings[itmp+1][parentEdge->getParentID()],rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,1),partMappings[itmp+1][parentEdge->getParentID()],ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,2),partMappings[itmp+1][parentEdge->getParentID()],rz);
                                        
                                        // Rotation 3
                                        rx = 0;
                                        ry = z;
                                        rz = -y;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,0),partMappings[itmp+2][parentEdge->getParentID()],rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,1),partMappings[itmp+2][parentEdge->getParentID()],ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(parentEdge->getGammaDofID(j,2),partMappings[itmp+2][parentEdge->getParentID()],rz);
                                        
                                        for (UN jj=0; jj<face->getNumNodes(); jj++) {
                                            x = localNodeList[parentEdge->getLocalNodeID(j)][0] - localNodeList[face->getLocalNodeID(jj)][0];
                                            y = localNodeList[parentEdge->getLocalNodeID(j)][1] - localNodeList[face->getLocalNodeID(jj)][1];
                                            z = localNodeList[parentEdge->getLocalNodeID(j)][2] - localNodeList[face->getLocalNodeID(jj)][2];
                                            faceValues[jj] += 1.0/sqrt(x*x+y*y+z*z);
                                        }
                                    }
                                }
                            }
                            for (UN ii=0; ii<parentEdges->getNumEntities(); ii++) {
                                InterfaceEntityPtr parentEdge = parentEdges->getEntity(ii);
                                EntitySetPtr parentVertices = parentEdge->getParents();
                                if (parentVertices->getNumEntities()==0) {
                                    for (UN jj=0; jj<face->getNumNodes(); jj++) {
                                        faceValue = 0.0;
                                        for (UN j=0; j<parentEdge->getNumNodes(); j++) {
                                            x = localNodeList[parentEdge->getLocalNodeID(j)][0] - localNodeList[face->getLocalNodeID(jj)][0];
                                            y = localNodeList[parentEdge->getLocalNodeID(j)][1] - localNodeList[face->getLocalNodeID(jj)][1];
                                            z = localNodeList[parentEdge->getLocalNodeID(j)][2] - localNodeList[face->getLocalNodeID(jj)][2];
                                            faceValue += (1.0/sqrt(x*x+y*y+z*z))/(faceValues[j]);
                                        }
                                        x = localNodeList[face->getLocalNodeID(jj)][0];
                                        y = localNodeList[face->getLocalNodeID(jj)][1];
                                        z = localNodeList[face->getLocalNodeID(jj)][2];
                                        
                                        // Rotation 1
                                        rx = y;
                                        ry = -x;
                                        rz = 0;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(jj,0),partMappings[itmp][parentEdge->getParentID()],faceValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(jj,1),partMappings[itmp][parentEdge->getParentID()],faceValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(jj,2),partMappings[itmp][parentEdge->getParentID()],faceValue*rz);
                                        
                                        // Rotation 2
                                        rx = -z;
                                        ry = 0;
                                        rz = x;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(jj,0),partMappings[itmp+1][parentEdge->getParentID()],faceValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(jj,1),partMappings[itmp+1][parentEdge->getParentID()],faceValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(jj,2),partMappings[itmp+1][parentEdge->getParentID()],faceValue*rz);
                                        
                                        // Rotation 3
                                        rx = 0;
                                        ry = z;
                                        rz = -y;
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(jj,0),partMappings[itmp+2][parentEdge->getParentID()],faceValue*rx);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(jj,1),partMappings[itmp+2][parentEdge->getParentID()],faceValue*ry);
                                        this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(jj,2),partMappings[itmp+2][parentEdge->getParentID()],faceValue*rz);
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
                        EntitySetPtr parentEdges = face->getParents();
                        if (parentEdges->getNumEntities()==0) {
                            for (UN j=0; j<face->getNumNodes(); j++) {
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,k),partMappings[itmp][face->getParentID()],1.0);
                            }
                        }
                    }
                    itmp++;
                }
                
                // Faces rotations
                if (buildRotations) {
                    for (UN i=0; i<faces->getNumEntities(); i++) {
                        InterfaceEntityPtr face = faces->getEntity(i);
                        EntitySetPtr parentEdges = face->getParents();
                        if (parentEdges->getNumEntities()==0) {
                            for (UN j=0; j<face->getNumNodes(); j++) {
                                x = localNodeList[face->getLocalNodeID(j)][0];
                                y = localNodeList[face->getLocalNodeID(j)][1];
                                z = localNodeList[face->getLocalNodeID(j)][2];
                                
                                // Rotation 1
                                rx = y;
                                ry = -x;
                                rz = 0;
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp][face->getParentID()],rx);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp][face->getParentID()],ry);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp][face->getParentID()],rz);
                                
                                // Rotation 2
                                rx = -z;
                                ry = 0;
                                rz = x;
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp+1][face->getParentID()],rx);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp+1][face->getParentID()],ry);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp+1][face->getParentID()],rz);
                                
                                // Rotation 3
                                rx = 0;
                                ry = z;
                                rz = -y;
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,0),partMappings[itmp+2][face->getParentID()],rx);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,1),partMappings[itmp+2][face->getParentID()],ry);
                                this->MVPhiGamma_[blockId]->replaceLocalValue(face->getGammaDofID(j,2),partMappings[itmp+2][face->getParentID()],rz);
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
