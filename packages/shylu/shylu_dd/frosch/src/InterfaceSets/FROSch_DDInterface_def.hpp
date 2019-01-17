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
    
    template <class SC,class LO,class GO,class NO>
    DDInterface<SC,LO,GO,NO>::DDInterface(UN dimension,
                                          UN dofsPerNode,
                                          MapPtr localToGlobalMap) :
    MpiComm_ (localToGlobalMap->getComm()),
    Dimension_ (dimension),
    DofsPerNode_ (dofsPerNode),
    NumMyNodes_ (localToGlobalMap->getNodeNumElements()),
    Vertices_ (new EntitySet<SC,LO,GO,NO>(VertexType)),
    Edges_ (new EntitySet<SC,LO,GO,NO>(EdgeType)),
    Faces_ (new EntitySet<SC,LO,GO,NO>(FaceType)),
    Interface_ (new EntitySet<SC,LO,GO,NO>(InterfaceType)),
    Interior_ (new EntitySet<SC,LO,GO,NO>(InteriorType)),
    CoarseNodes_ (new EntitySet<SC,LO,GO,NO>(DefaultType)),
    EntityVector_ (),
    NodesMap_ (localToGlobalMap),
    UniqueNodesMap_ ()
    {
        FROSCH_ASSERT(((Dimension_==2)||(Dimension_==3)),"Only dimension 2 and 3 are available");
        
        UniqueNodesMap_ = BuildUniqueMap<LO,GO,NO>(NodesMap_);
        
        GOVecVecPtr componentsSubdomains;
        GOVecVec componentsSubdomainsUnique;
        
        communicateLocalComponents(componentsSubdomains,componentsSubdomainsUnique);
        identifyLocalComponents(componentsSubdomains,componentsSubdomainsUnique);
    }
    
    template <class SC,class LO,class GO,class NO>
    DDInterface<SC,LO,GO,NO>::~DDInterface()
    {
        
    } // Do we need sth here?
    
    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::resetGlobalDofs(MapPtrVecPtr dofsMaps)
    {
        // Vertices
        for (UN i=0; i<Vertices_->getNumEntities(); i++) {
            for (UN j=0; j<Vertices_->getEntity(i)->getNumNodes(); j++) {
                LO localID = Vertices_->getEntity(i)->getLocalNodeID(j);
                UNVecPtr dofIDs(DofsPerNode_);
                GOVecPtr dofsGlobal(DofsPerNode_);
                for (UN k=0; k<DofsPerNode_; k++) {
                    dofIDs[k] = k;
                    dofsGlobal[k] = dofsMaps[k]->getGlobalElement(localID);
                }
                Vertices_->getEntity(i)->resetGlobalDofs(j,DofsPerNode_,&(dofIDs[0]),&(dofsGlobal[0]));
            }
        }
        
        // ShortEdges
        for (UN i=0; i<ShortEdges_->getNumEntities(); i++) {
            for (UN j=0; j<ShortEdges_->getEntity(i)->getNumNodes(); j++) {
                LO localID = ShortEdges_->getEntity(i)->getLocalNodeID(j);
                UNVecPtr dofIDs(DofsPerNode_);
                GOVecPtr dofsGlobal(DofsPerNode_);
                for (UN k=0; k<DofsPerNode_; k++) {
                    dofIDs[k] = k;
                    dofsGlobal[k] = dofsMaps[k]->getGlobalElement(localID);
                }
                ShortEdges_->getEntity(i)->resetGlobalDofs(j,DofsPerNode_,&(dofIDs[0]),&(dofsGlobal[0]));
            }
        }
        
        // StraightEdges
        for (UN i=0; i<StraightEdges_->getNumEntities(); i++) {
            for (UN j=0; j<StraightEdges_->getEntity(i)->getNumNodes(); j++) {
                LO localID = StraightEdges_->getEntity(i)->getLocalNodeID(j);
                UNVecPtr dofIDs(DofsPerNode_);
                GOVecPtr dofsGlobal(DofsPerNode_);
                for (UN k=0; k<DofsPerNode_; k++) {
                    dofIDs[k] = k;
                    dofsGlobal[k] = dofsMaps[k]->getGlobalElement(localID);
                }
                StraightEdges_->getEntity(i)->resetGlobalDofs(j,DofsPerNode_,&(dofIDs[0]),&(dofsGlobal[0]));
            }
        }
        
        // Edges
        for (UN i=0; i<Edges_->getNumEntities(); i++) {
            for (UN j=0; j<Edges_->getEntity(i)->getNumNodes(); j++) {
                LO localID = Edges_->getEntity(i)->getLocalNodeID(j);
                UNVecPtr dofIDs(DofsPerNode_);
                GOVecPtr dofsGlobal(DofsPerNode_);
                for (UN k=0; k<DofsPerNode_; k++) {
                    dofIDs[k] = k;
                    dofsGlobal[k] = dofsMaps[k]->getGlobalElement(localID);
                }
                Edges_->getEntity(i)->resetGlobalDofs(j,DofsPerNode_,&(dofIDs[0]),&(dofsGlobal[0]));
            }
        }
        
        // Faces
        for (UN i=0; i<Faces_->getNumEntities(); i++) {
            for (UN j=0; j<Faces_->getEntity(i)->getNumNodes(); j++) {
                LO localID = Faces_->getEntity(i)->getLocalNodeID(j);
                UNVecPtr dofIDs(DofsPerNode_);
                GOVecPtr dofsGlobal(DofsPerNode_);
                for (UN k=0; k<DofsPerNode_; k++) {
                    dofIDs[k] = k;
                    dofsGlobal[k] = dofsMaps[k]->getGlobalElement(localID);
                }
                Faces_->getEntity(i)->resetGlobalDofs(j,DofsPerNode_,&(dofIDs[0]),&(dofsGlobal[0]));
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
        // Vertices
        for (UN i=0; i<Vertices_->getNumEntities(); i++) {
            UN length = Vertices_->getEntity(i)->getNumNodes();
            for (UN j=0; j<length; j++) {
                UN itmp = length-1-j;
                UN k = 0;
                while (k<DofsPerNode_) {
                    GO dofGlobal = Vertices_->getEntity(i)->getGlobalDofID(itmp,k);
                    if (std::binary_search(dirichletBoundaryDofs.begin(),dirichletBoundaryDofs.end(),dofGlobal)) {
                        Vertices_->getEntity(i)->removeNode(itmp);
                        break;
                    }
                    k++;
                }
            }
        }
        Vertices_->removeEmptyEntities();
        
        // ShortEdges
        for (UN i=0; i<ShortEdges_->getNumEntities(); i++) {
            UN length = ShortEdges_->getEntity(i)->getNumNodes();
            for (UN j=0; j<length; j++) {
                UN itmp = length-1-j;
                UN k = 0;
                while (k<DofsPerNode_) {
                    GO dofGlobal = ShortEdges_->getEntity(i)->getGlobalDofID(itmp,k);
                    if (std::binary_search(dirichletBoundaryDofs.begin(),dirichletBoundaryDofs.end(),dofGlobal)) {
                        ShortEdges_->getEntity(i)->removeNode(itmp);
                        break;
                    }
                    k++;
                }
            }
        }
        ShortEdges_->removeEmptyEntities();
        
        // StraightEdges
        for (UN i=0; i<StraightEdges_->getNumEntities(); i++) {
            UN length = StraightEdges_->getEntity(i)->getNumNodes();
            for (UN j=0; j<length; j++) {
                UN itmp = length-1-j;
                UN k = 0;
                while (k<DofsPerNode_) {
                    GO dofGlobal = StraightEdges_->getEntity(i)->getGlobalDofID(itmp,k);
                    if (std::binary_search(dirichletBoundaryDofs.begin(),dirichletBoundaryDofs.end(),dofGlobal)) {
                        StraightEdges_->getEntity(i)->removeNode(itmp);
                        break;
                    }
                    k++;
                }
            }
        }
        StraightEdges_->removeEmptyEntities();
        
        // Edges
        for (UN i=0; i<Edges_->getNumEntities(); i++) {
            UN length = Edges_->getEntity(i)->getNumNodes();
            for (UN j=0; j<length; j++) {
                UN itmp = length-1-j;
                UN k = 0;
                while (k<DofsPerNode_) {
                    GO dofGlobal = Edges_->getEntity(i)->getGlobalDofID(itmp,k);
                    if (std::binary_search(dirichletBoundaryDofs.begin(),dirichletBoundaryDofs.end(),dofGlobal)) {
                        Edges_->getEntity(i)->removeNode(itmp);
                        break;
                    }
                    k++;
                }
            }
        }
        Edges_->removeEmptyEntities();
        
        // Faces
        for (UN i=0; i<Faces_->getNumEntities(); i++) {
            UN length = Faces_->getEntity(i)->getNumNodes();
            for (UN j=0; j<length; j++) {
                UN itmp = length-1-j;
                UN k = 0;
                while (k<DofsPerNode_) {
                    GO dofGlobal = Faces_->getEntity(i)->getGlobalDofID(itmp,k);
                    if (std::binary_search(dirichletBoundaryDofs.begin(),dirichletBoundaryDofs.end(),dofGlobal)) {
                        Faces_->getEntity(i)->removeNode(itmp);
                        break;
                    }
                    k++;
                }
            }
        }
        Faces_->removeEmptyEntities();
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::divideUnconnectedEntities(CrsMatrixPtr matrix)
    {
#ifdef INTERFACE_OUTPUT
        GOVecPtr indicesGammaDofs(DofsPerNode_*Interface_->getEntity(0)->getNumNodes());
        for (UN k=0; k<DofsPerNode_; k++) {
            for (UN i=0; i<Interface_->getEntity(0)->getNumNodes(); i++) {
                indicesGammaDofs[Interface_->getEntity(0)->getGammaDofID(i,k)] = Interface_->getEntity(0)->getGlobalDofID(i,k);
            }
        }
        MapPtr map = Xpetra::MapFactory<LO,GO,NO>::Build(matrix->getRowMap()->lib(),-1,indicesGammaDofs(),0,MpiComm_);
        matrix = FROSch::ExtractLocalSubdomainMatrix(matrix,map,1.0);
        LO numSeparateEdges = Edges_->divideUnconnectedEntities(matrix);
        LO numSeparateFaces = Faces_->divideUnconnectedEntities(matrix);
        
        if (MpiComm_->getRank()==0) {
            std::cout << "\n\
            --------------------------------------------\n\
            # separate edges:     --- " << numSeparateEdges << "\n\
            # separate faces:     --- " << numSeparateFaces << "\n\
            --------------------------------------------\n";
        }
#else
        if (MpiComm_->getRank()==0) std::cout << "WARNING: divideUnconnectedEntities has not been really tested yet...\n";
        
        GOVecPtr indicesGammaDofs(DofsPerNode_*Interface_->getEntity(0)->getNumNodes());
        for (UN k=0; k<DofsPerNode_; k++) {
            for (UN i=0; i<Interface_->getEntity(0)->getNumNodes(); i++) {
                indicesGammaDofs[Interface_->getEntity(0)->getGammaDofID(i,k)] = Interface_->getEntity(0)->getGlobalDofID(i,k);
            }
        }
        MapPtr map = Xpetra::MapFactory<LO,GO,NO>::Build(matrix->getRowMap()->lib(),-1,indicesGammaDofs(),0,MpiComm_);
        matrix = FROSch::ExtractLocalSubdomainMatrix(matrix,map,1.0);
        Edges_->divideUnconnectedEntities(matrix,MpiComm_->getRank());
        Faces_->divideUnconnectedEntities(matrix,MpiComm_->getRank());
#endif
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::sortEntities()
    {
        // Edges_
        InterfaceEntityPtrVecPtr extractedVertices = Edges_->sortOutNodes();
        for (UN i=0; i<extractedVertices.size(); i++) {
            Vertices_->addEntity(extractedVertices[i]);
        }
        
        InterfaceEntityPtrVecPtr extractedShortEdges = Edges_->sortOutShortEntities();
        for (UN i=0; i<extractedShortEdges.size(); i++) {
            ShortEdges_->addEntity(extractedShortEdges[i]);
        }
        
        // Faces_
        extractedVertices = Faces_->sortOutNodes();
        for (UN i=0; i<extractedVertices.size(); i++) {
            Vertices_->addEntity(extractedVertices[i]);
        }
        
        InterfaceEntityPtrVecPtr extractedShortFaces = Faces_->sortOutShortEntities();
        for (UN i=0; i<extractedShortEdges.size(); i++) {
            ShortFaces_->addEntity(extractedShortFaces[i]);
        }
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::sortEntities(MultiVectorPtr nodeList)
    {
        if (nodeList.is_null()) {
            this->sortEntities();
        } else {
            // Edges_
            InterfaceEntityPtrVecPtr extractedVertices = Edges_->sortOutNodes();
            for (UN i=0; i<extractedVertices.size(); i++) {
                Vertices_->addEntity(extractedVertices[i]);
            }
            
            InterfaceEntityPtrVecPtr extractedShortEdges = Edges_->sortOutShortEntities();
            for (UN i=0; i<extractedShortEdges.size(); i++) {
                ShortEdges_->addEntity(extractedShortEdges[i]);
            }
            
            InterfaceEntityPtrVecPtr extractedStraightEdges = Edges_->sortOutStraightEntities(Dimension_,nodeList);
            for (UN i=0; i<extractedStraightEdges.size(); i++) {
                StraightEdges_->addEntity(extractedStraightEdges[i]);
            }
            
            // Faces_
            extractedVertices = Faces_->sortOutNodes();
            for (UN i=0; i<extractedVertices.size(); i++) {
                Vertices_->addEntity(extractedVertices[i]);
            }
            
            InterfaceEntityPtrVecPtr extractedShortFaces = Faces_->sortOutShortEntities();
            for (UN i=0; i<extractedShortFaces.size(); i++) {
                ShortFaces_->addEntity(extractedShortFaces[i]);
            }
            
            InterfaceEntityPtrVecPtr extractedStraightFaces = Faces_->sortOutStraightEntities(Dimension_,nodeList);
            for (UN i=0; i<extractedStraightFaces.size(); i++) {
                StraightFaces_->addEntity(extractedStraightFaces[i]);
            }
        }
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::buildEntityHierarchy()
    {
        for (UN i=0; i<EntityVector.size(); i++) {
            for (UN j=i+1; j<EntityVector.size(); j++) {
                EntityVector[i]->findAncestors(EntityVector[j]);
            }
        }
        
        UN itmp = 0;
        for (UN i=0; i<EntityVector.size(); i++) {
            for (UN j=0; j<EntityVector[i]->getNumEntities(); j++) {
                tmpAncestors = EntityVector[i]->getEntity(j)->getAncestors();
                if (tmpAncestors->getNumEntities() == 0) {
                    CoarseNodes_->addEntity(EntityVector[i]->getEntity(j));
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
    typename DDInterface<SC,LO,GO,NO>::EntitySetConstPtr & DDInterface<SC,LO,GO,NO>::getCoarseNodes() const
    {
        return CoarseNodes_;
    }
    
    template <class SC,class LO,class GO,class NO>
    typename DDInterface<SC,LO,GO,NO>::EntitySetPtrConstVecPtr & DDInterface<SC,LO,GO,NO>::getEntityVector() const
    {
        return EntityVector_;
    }
    
    template <class SC,class LO,class GO,class NO>
    typename DDInterface<SC,LO,GO,NO>::ConstMapPtr DDInterface<SC,LO,GO,NO>::getNodesMap() const
    {
        return NodesMap_.getConst();
    }

    // Issue with OpenMP and offset maps: When starting with OffsetMap (GID not starting at 0) in constructor of an Xpetra::Matrix with underlyinglib=tpetra that is not square (fillComplete called with this OffsetMap as RangeMap and a continuous DomainMap starting at GID 0). Interface will be identified with offset removed if "FROSCH_OFFSET_MAPS" is defined. Otherwise the previous (standard) communicateLocalComponents method will be used.
#ifdef FROSCH_OFFSET_MAPS
    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::communicateLocalComponents(GOVecVecPtr &componentsSubdomains,
                                                             GOVecVec &componentsSubdomainsUnique,
                                                             UN priorDofsPerNode)
    {
        

        if (UniqueNodesMap_->getMinAllGlobalIndex() > 0 || priorDofsPerNode > 0) {
            GO minGID = UniqueNodesMap_->getMinAllGlobalIndex();
            ConstGOVecView ElementList;
            MapPtr tmpUniqueNodesMap;
            {
                GOVec indices(UniqueNodesMap_->getNodeNumElements());
                ElementList = UniqueNodesMap_->getNodeElementList();
                for (UN i=0; i<indices.size(); i++) {
                    indices[i] = ElementList[i] - minGID - (GO) i*priorDofsPerNode;
                }
                tmpUniqueNodesMap = Xpetra::MapFactory<LO,GO,NO>::Build(UniqueNodesMap_->lib(),-1,indices(),0,UniqueNodesMap_->getComm());
            }

            MapPtr tmpNodesMap;
            {
                GOVec indices(NodesMap_->getNodeNumElements());
                ElementList = NodesMap_->getNodeElementList();
                for (UN i=0; i<indices.size(); i++) {
                    indices[i] = ElementList[i] - minGID - (GO) i*priorDofsPerNode;
                }
                tmpNodesMap = Xpetra::MapFactory<LO,GO,NO>::Build(NodesMap_->lib(),-1,indices(),0,NodesMap_->getComm());
            }

            Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > commMat = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(tmpNodesMap,10);
            Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > commMatTmp = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(tmpUniqueNodesMap,10);
#ifdef Tpetra_issue_1752
            // AH 10/10/2017: Can we get away with using just one importer/exporter after the Tpetra issue is fixed?
            Teuchos::RCP<Xpetra::Import<LO,GO,NO> > commImporter = Xpetra::ImportFactory<LO,GO,NO>::Build(tmpUniqueNodesMap,tmpNodesMap);
            Teuchos::RCP<Xpetra::Export<LO,GO,NO> > commExporter = Xpetra::ExportFactory<LO,GO,NO>::Build(tmpNodesMap,tmpUniqueNodesMap);
#else
            Teuchos::RCP<Xpetra::Export<LO,GO,NO> > commExporter = Xpetra::ExportFactory<LO,GO,NO>::Build(tmpNodesMap,tmpUniqueNodesMap);
#endif

            Teuchos::Array<SC> one(1,Teuchos::ScalarTraits<SC>::one());
            Teuchos::Array<GO> myPID(1,tmpUniqueNodesMap->getComm()->getRank());
            for (int i=0; i<NumMyNodes_; i++) {
                commMat->insertGlobalValues(tmpNodesMap->getGlobalElement(i),myPID(),one());
            }

            Teuchos::RCP<Xpetra::Map<LO,GO,NO> > rangeMap = Xpetra::MapFactory<LO,GO,NO>::Build(tmpNodesMap->lib(),-1,myPID(),0,tmpNodesMap->getComm());
            
            commMat->fillComplete(tmpNodesMap,rangeMap);
            
            commMatTmp->doExport(*commMat,*commExporter,Xpetra::INSERT);
            
            commMatTmp->fillComplete(tmpUniqueNodesMap,rangeMap);

            
            commMat = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(tmpNodesMap,10);
#ifdef Tpetra_issue_1752
            commMat->doImport(*commMatTmp,*commImporter,Xpetra::INSERT);
#else
            commMat->doImport(*commMatTmp,*commExporter,Xpetra::INSERT);
#endif
            
            componentsSubdomains = GOVecVecPtr(NumMyNodes_);
            componentsSubdomainsUnique = GOVecVec(NumMyNodes_);
            Teuchos::ArrayView<const GO> indices2;
            Teuchos::ArrayView<const SC> values2;
            for (LO i=0; i<NumMyNodes_; i++) {
                commMat->getGlobalRowView(tmpNodesMap->getGlobalElement(i),indices2,values2);
                componentsSubdomains[i].resize(indices2.size());
                componentsSubdomainsUnique[i].resize(indices2.size());
                for (LO j=0; j<indices2.size(); j++) {
                    componentsSubdomains[i][j] = indices2[j];
                    componentsSubdomainsUnique[i][j] = indices2[j];
                }
                sortunique(componentsSubdomains[i]);
                sortunique(componentsSubdomainsUnique[i]);
            }
            
            sortunique(componentsSubdomainsUnique);


        } else {
            Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > commMat = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(NodesMap_,10);
            Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > commMatTmp = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(UniqueNodesMap_,10);
    #ifdef Tpetra_issue_1752
            // AH 10/10/2017: Can we get away with using just one importer/exporter after the Tpetra issue is fixed?
            Teuchos::RCP<Xpetra::Import<LO,GO,NO> > commImporter = Xpetra::ImportFactory<LO,GO,NO>::Build(UniqueNodesMap_,NodesMap_);
            Teuchos::RCP<Xpetra::Export<LO,GO,NO> > commExporter = Xpetra::ExportFactory<LO,GO,NO>::Build(NodesMap_,UniqueNodesMap_);
    #else
            Teuchos::RCP<Xpetra::Export<LO,GO,NO> > commExporter = Xpetra::ExportFactory<LO,GO,NO>::Build(NodesMap_,UniqueNodesMap_);
    #endif
            
            Teuchos::Array<SC> one(1,Teuchos::ScalarTraits<SC>::one());
            Teuchos::Array<GO> myPID(1,UniqueNodesMap_->getComm()->getRank());
            for (int i=0; i<NumMyNodes_; i++) {
                commMat->insertGlobalValues(NodesMap_->getGlobalElement(i),myPID(),one());
            }
            Teuchos::RCP<Xpetra::Map<LO,GO,NO> > rangeMap = Xpetra::MapFactory<LO,GO,NO>::Build(NodesMap_->lib(),-1,myPID(),0,NodesMap_->getComm());
            
            commMat->fillComplete(NodesMap_,rangeMap);
            
            commMatTmp->doExport(*commMat,*commExporter,Xpetra::INSERT);
            
            commMatTmp->fillComplete(UniqueNodesMap_,rangeMap);
            
            commMat = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(NodesMap_,10);
    #ifdef Tpetra_issue_1752
            commMat->doImport(*commMatTmp,*commImporter,Xpetra::INSERT);
    #else
            commMat->doImport(*commMatTmp,*commExporter,Xpetra::INSERT);
    #endif
            
            componentsSubdomains = GOVecVecPtr(NumMyNodes_);
            componentsSubdomainsUnique = GOVecVec(NumMyNodes_);
            
            Teuchos::ArrayView<const GO> indices2;
            Teuchos::ArrayView<const SC> values2;
            for (LO i=0; i<NumMyNodes_; i++) {
                commMat->getGlobalRowView(NodesMap_->getGlobalElement(i),indices2,values2);
                componentsSubdomains[i].resize(indices2.size());
                componentsSubdomainsUnique[i].resize(indices2.size());
                for (LO j=0; j<indices2.size(); j++) {
                    componentsSubdomains[i][j] = indices2[j];
                    componentsSubdomainsUnique[i][j] = indices2[j];
                }
                sortunique(componentsSubdomains[i]);
                sortunique(componentsSubdomainsUnique[i]);
            }
            
            sortunique(componentsSubdomainsUnique);
        }
        return 0;
    }
#else
    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::communicateLocalComponents(GOVecVecPtr &componentsSubdomains,
                                                             GOVecVec &componentsSubdomainsUnique)
    {
        Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > commMat = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(NodesMap_,10);
        Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > commMatTmp = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(UniqueNodesMap_,10);
#ifdef Tpetra_issue_1752
        // AH 10/10/2017: Can we get away with using just one importer/exporter after the Tpetra issue is fixed?
        Teuchos::RCP<Xpetra::Import<LO,GO,NO> > commImporter = Xpetra::ImportFactory<LO,GO,NO>::Build(UniqueNodesMap_,NodesMap_);
        Teuchos::RCP<Xpetra::Export<LO,GO,NO> > commExporter = Xpetra::ExportFactory<LO,GO,NO>::Build(NodesMap_,UniqueNodesMap_);
#else
        Teuchos::RCP<Xpetra::Export<LO,GO,NO> > commExporter = Xpetra::ExportFactory<LO,GO,NO>::Build(NodesMap_,UniqueNodesMap_);
#endif

        Teuchos::Array<SC> one(1,Teuchos::ScalarTraits<SC>::one());
        Teuchos::Array<GO> myPID(1,UniqueNodesMap_->getComm()->getRank());
        for (int i=0; i<NumMyNodes_; i++) {
            commMat->insertGlobalValues(NodesMap_->getGlobalElement(i),myPID(),one());
        }
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > rangeMap = Xpetra::MapFactory<LO,GO,NO>::Build(NodesMap_->lib(),-1,myPID(),0,NodesMap_->getComm());

        commMat->fillComplete(NodesMap_,rangeMap);

        commMatTmp->doExport(*commMat,*commExporter,Xpetra::INSERT);

        commMatTmp->fillComplete(UniqueNodesMap_,rangeMap);

        commMat = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(NodesMap_,10);
#ifdef Tpetra_issue_1752
        commMat->doImport(*commMatTmp,*commImporter,Xpetra::INSERT);
#else
        commMat->doImport(*commMatTmp,*commExporter,Xpetra::INSERT);
#endif

        componentsSubdomains = GOVecVecPtr(NumMyNodes_);
        componentsSubdomainsUnique = GOVecVec(NumMyNodes_);

        Teuchos::ArrayView<const GO> indices2;
        Teuchos::ArrayView<const SC> values2;
        for (LO i=0; i<NumMyNodes_; i++) {
            commMat->getGlobalRowView(NodesMap_->getGlobalElement(i),indices2,values2);
            componentsSubdomains[i].resize(indices2.size());
            componentsSubdomainsUnique[i].resize(indices2.size());
            for (LO j=0; j<indices2.size(); j++) {
                componentsSubdomains[i][j] = indices2[j];
                componentsSubdomainsUnique[i][j] = indices2[j];
            }
            sortunique(componentsSubdomains[i]);
            sortunique(componentsSubdomainsUnique[i]);
        }
        
        sortunique(componentsSubdomainsUnique);
        
        return 0;
    }
#endif

    
    template <class SC,class LO,class GO,class NO>
    int DDInterface<SC,LO,GO,NO>::identifyLocalComponents(GOVecVecPtr &componentsSubdomains,
                                                          GOVecVec &componentsSubdomainsUnique)
    {
        // Hier herausfinden, ob Ecke, Kante oder Fläche
        UNVecPtr componentsMultiplicity(componentsSubdomainsUnique.size());
        GOVecVecPtr components(componentsSubdomainsUnique.size());
        GOVecVecPtr componentsGamma(componentsSubdomainsUnique.size());
        UN maxMultiplicity = 0;
        for (UN i=0; i<componentsSubdomainsUnique.size(); i++) {
            componentsMultiplicity[i] = componentsSubdomainsUnique[i].size();
            maxMultiplicity = std::max(maxMultiplicity,componentsMultiplicity[i]);
        }
        EntityVector_ = EntitySetPtrVecPtr(maxMultiplicity);
        for (UN i=0; i<maxMultiplicity; i++) {
            EntityVector_[i].reset(new EntitySet<SC,LO,GO,NO>(DefaultType));
        }
        
        typename GOVecVec::iterator classIterator;
        LOVecPtr localComponentIndices(NumMyNodes_);
        for (int i=0; i<NumMyNodes_; i++) {
            classIterator = std::lower_bound(componentsSubdomainsUnique.begin(),componentsSubdomainsUnique.end(),componentsSubdomains[i]);
            localComponentIndices[i] = classIterator - componentsSubdomainsUnique.begin();
        }
        
        LO tmp1 = 0;
        GO *tmp2 = NULL;
        Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> > interior(new InterfaceEntity<SC,LO,GO,NO>(InteriorType,DofsPerNode_,tmp1,tmp2));
        Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> > interface(new InterfaceEntity<SC,LO,GO,NO>(InterfaceType,DofsPerNode_,tmp1,tmp2));
        for (LO i=0; i<NumMyNodes_; i++) {
            if (componentsMultiplicity[localComponentIndices[i]] == 1) {
                LO nodeIDI = volume->getNumNodes();
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
                LO nodeIDGamma = surface->getNumNodes();
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
                componentsGamma[localComponentIndices[i]].push_back(surface->getNumNodes()-1);
            }
        }
        Interior_->addEntity(interior);
        Interface_->addEntity(interface);
        
        for (UN i=0; i<componentsSubdomainsUnique.size(); i++) {
            Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> > tmpEntity(new InterfaceEntity<SC,LO,GO,NO>(VertexType,DofsPerNode_,componentsMultiplicity[i],&(componentsSubdomainsUnique[i][0])));
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
            EntityVector_[componentsMultiplicity[i]]->addEntity(tmpEntity);
        }

        return 0;
    }
    
}

#endif
