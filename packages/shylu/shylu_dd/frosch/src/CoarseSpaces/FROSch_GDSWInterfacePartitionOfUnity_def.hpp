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

#ifndef _FROSCH_GDSWINTERFACEPARTITIONOFUNITY_DEF_HPP
#define _FROSCH_GDSWINTERFACEPARTITIONOFUNITY_DEF_HPP

#include <FROSch_GDSWInterfacePartitionOfUnity_decl.hpp>

namespace FROSch {
    
    template <class SC,class LO,class GO,class NO>
    GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::GDSWInterfacePartitionOfUnity(CommPtr mpiComm,
                                                                              CommPtr serialComm,
                                                                              UN dimension,
                                                                              UN dofsPerNode,
                                                                              MapPtr nodesMap,
                                                                              MapPtrVecPtr dofsMaps,
                                                                              ParameterListPtr parameterList) :
    InterfacePartitionOfUnity<SC,LO,GO,NO> (mpiComm,serialComm,dimension,dofsPerNode,nodesMap,dofsMaps,sublist(parameterList,"GDSW")),
    UseVertices_ (false),
    UseShortEdges_ (false),
    UseStraightEdges_ (false),
    UseEdges_ (false),
    UseFaces_ (false),
    Vertices_ (),
    ShortEdges_ (),
    StraightEdges_ (),
    Edges_ (),
    Faces_ ()
    {
        if (!this->ParameterList_->get("Type","Full").compare("Full")) {
            UseVertices_ = true;
            UseShortEdges_ = true;
            UseStraightEdges_ = true;
            UseEdges_ = true;
            UseFaces_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("Vertices")) {
            UseVertices_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("ShortEdges")) {
            UseShortEdges_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("StraightEdges")) {
            UseStraightEdges_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("Edges")) {
            UseEdges_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("AllEdges")) {
            UseFaces_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("NoVertices")) {
            UseShortEdges_ = true;
            UseStraightEdges_ = true;
            UseEdges_ = true;
            UseFaces_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("NoShortEdges")) {
            UseVertices_ = true;
            UseStraightEdges_ = true;
            UseEdges_ = true;
            UseFaces_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("NoStraightEdges")) {
            UseVertices_ = true;
            UseShortEdges_ = true;
            UseEdges_ = true;
            UseFaces_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("NoEdges")) {
            UseVertices_ = true;
            UseShortEdges_ = true;
            UseStraightEdges_ = true;
            UseFaces_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("NoAllEdges")) {
            UseVertices_ = true;
            UseFaces_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("NoFaces")) {
            UseVertices_ = true;
            UseShortEdges_ = true;
            UseStraightEdges_ = true;
            UseEdges_ = true;
        } else if (!this->ParameterList_->get("Type","Full").compare("Custom")) {
            UseVertices_ = this->ParameterList_->sublist("Custom").get("Vertices",false);
            UseShortEdges_ = this->ParameterList_->sublist("Custom").get("Vertices",false);
            UseStraightEdges_ = this->ParameterList_->sublist("Custom").get("Vertices",false);
            UseEdges_ = this->ParameterList_->sublist("Custom").get("Vertices",false);
            UseFaces_ = this->ParameterList_->sublist("Custom").get("Vertices",false);
        }
        this->LocalPartitionOfUnity_ = MultiVectorPtrVecPtr(5);
        this->PartitionOfUnityMaps_ = MapPtrVecPtr(5);
    }
    
    template <class SC,class LO,class GO,class NO>
    GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::~GDSWInterfacePartitionOfUnity()
    {
        
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::removeDirichletNodes(GOVecView dirichletBoundaryDofs,
                                                                         MultiVectorPtr nodeList)
    {
        if (!dirichletBoundaryDofs.is_null()) {
            GOVec tmpDirichletBoundaryDofs(dirichletBoundaryDofs());
            sortunique(tmpDirichletBoundaryDofs);            
            this->DDInterface_->removeDirichletNodes(tmpDirichletBoundaryDofs());
            this->DDInterface_->sortEntities(nodeList);
        }
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::sortInterface(CrsMatrixPtr matrix,
                                                                  MultiVectorPtr nodeList)
    {
        this->DDInterface_->divideUnconnectedEntities(matrix);
        this->DDInterface_->sortEntities(nodeList);
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>::computePartitionOfUnity()
    {
        // Interface
        UN dofsPerNode = this->DDInterface_->getInterface()->getEntity(0)->getDofsPerNode();
        UN numInterfaceDofs = dofsPerNode*this->DDInterface_->getInterface()->getEntity(0)->getNumNodes();
        
        
        // Maps
        LO numVerticesGlobal;
        if (UseVertices_) {
            Vertices_ = this->DDInterface_->getVertices();
            Vertices_->buildEntityMap(this->DDInterface_->getNodesMap());
            this->PartitionOfUnityMaps_[0] = Vertices_->getEntityMap();
            
            numVerticesGlobal = Vertices_->getEntityMap()->getMaxAllGlobalIndex();
            if (Vertices_->getEntityMap()->lib()==Xpetra::UseEpetra || Vertices_->getEntityMap()->getGlobalNumElements()>0) {
                numVerticesGlobal += 1;
            }
        } else {
            numVerticesGlobal = -1;
        }
        if (numVerticesGlobal < 0) numVerticesGlobal = 0;
        
        LO numShortEdgesGlobal;
        if (UseShortEdges_) {
            ShortEdges_ = this->DDInterface_->getShortEdges();
            ShortEdges_->buildEntityMap(this->DDInterface_->getNodesMap());
            this->PartitionOfUnityMaps_[1] = ShortEdges_->getEntityMap();
            
            numShortEdgesGlobal = ShortEdges_->getEntityMap()->getMaxAllGlobalIndex();
            if (ShortEdges_->getEntityMap()->lib()==Xpetra::UseEpetra || ShortEdges_->getEntityMap()->getGlobalNumElements()>0) {
                numShortEdgesGlobal += 1;
            }
        } else {
            numShortEdgesGlobal = -1;
        }
        if (numShortEdgesGlobal < 0) numShortEdgesGlobal = 0;
        
        LO numStraightEdgesGlobal;
        if (UseStraightEdges_) {
            StraightEdges_ = this->DDInterface_->getStraightEdges();
            StraightEdges_->buildEntityMap(this->DDInterface_->getNodesMap());
            this->PartitionOfUnityMaps_[2] = StraightEdges_->getEntityMap();
            
            numStraightEdgesGlobal = StraightEdges_->getEntityMap()->getMaxAllGlobalIndex();
            if (StraightEdges_->getEntityMap()->lib()==Xpetra::UseEpetra || StraightEdges_->getEntityMap()->getGlobalNumElements()>0) {
                numStraightEdgesGlobal += 1;
            }
        } else {
            numStraightEdgesGlobal = -1;
        }
        if (numStraightEdgesGlobal < 0) numStraightEdgesGlobal = 0;
        
        LO numEdgesGlobal;
        if (UseEdges_) {
            Edges_ = this->DDInterface_->getEdges();
            Edges_->buildEntityMap(this->DDInterface_->getNodesMap());
            this->PartitionOfUnityMaps_[3] = Edges_->getEntityMap();
            
            numEdgesGlobal = Edges_->getEntityMap()->getMaxAllGlobalIndex();
            if (Edges_->getEntityMap()->lib()==Xpetra::UseEpetra || Edges_->getEntityMap()->getGlobalNumElements()>0) {
                numEdgesGlobal += 1;
            }
        } else {
            numEdgesGlobal = -1;
        }
        if (numEdgesGlobal < 0) numEdgesGlobal = 0;
        
        LO numFacesGlobal;
        if (UseFaces_) {
            Faces_ = this->DDInterface_->getFaces();
            Faces_->buildEntityMap(this->DDInterface_->getNodesMap());
            this->PartitionOfUnityMaps_[4] = Faces_->getEntityMap();
            
            numFacesGlobal = Faces_->getEntityMap()->getMaxAllGlobalIndex();
            if (Faces_->getEntityMap()->lib()==Xpetra::UseEpetra || Faces_->getEntityMap()->getGlobalNumElements()>0) {
                numFacesGlobal += 1;
            }
        } else {
            numFacesGlobal = -1;
        }
        if (numFacesGlobal < 0) numFacesGlobal = 0;
        
        if (this->Verbose_) {
            std::cout << "-----------------------------------------------\n";
            std::cout << " GDSW Interface Partition Of Unity (GDSW IPOU) \n";
            std::cout << "-----------------------------------------------\n";
            std::cout << std::boolalpha;
            std::cout << " \tUse Vertices\t\t--\t" << UseVertices_ << "\n";
            std::cout << " \tUse Short Edges\t\t--\t" << UseShortEdges_ << "\n";
            std::cout << " \tUse Straight Edges\t--\t" << UseStraightEdges_ << "\n";
            std::cout << " \tUse Edges\t\t--\t" << UseEdges_ << "\n";
            std::cout << " \tUse Faces\t\t--\t" << UseFaces_ << "\n";
            std::cout << std::noboolalpha;
            std::cout << "-----------------------------------------------\n";
            std::cout << " \t# Vertices\t\t--\t" << numVerticesGlobal << "\n";
            std::cout << " \t# Short Edges\t\t--\t" << numShortEdgesGlobal << "\n";
            std::cout << " \t# Straight Edges\t--\t" << numStraightEdgesGlobal << "\n";
            std::cout << " \t# Edges\t\t\t--\t" << numEdgesGlobal << "\n";
            std::cout << " \t# Faces\t\t\t--\t" << numFacesGlobal << "\n";
            std::cout << "-----------------------------------------------\n";
        }
        
        // Build Partition Of Unity Vectors
        MapPtr serialInterfaceMap = Xpetra::MapFactory<LO,GO,NO>::Build(this->DDInterface_->getNodesMap()->lib(),numInterfaceDofs,0,this->SerialComm_);
        
        if (UseVertices_ && Vertices_->getNumEntities()>0) {
            MultiVectorPtr tmpVector = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(serialInterfaceMap,Vertices_->getNumEntities());
            
            for (UN i=0; i<Vertices_->getNumEntities(); i++) {
                for (UN j=0; j<Vertices_->getEntity(i)->getNumNodes(); j++) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        tmpVector->replaceLocalValue(Vertices_->getEntity(i)->getGammaDofID(j,k),i,1.0);
                    }
                }
            }
            
            this->LocalPartitionOfUnity_[0] = tmpVector;
        }
        
        if (UseShortEdges_ && ShortEdges_->getNumEntities()>0) {
            MultiVectorPtr tmpVector = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(serialInterfaceMap,ShortEdges_->getNumEntities());
            
            for (UN i=0; i<ShortEdges_->getNumEntities(); i++) {
                for (UN j=0; j<ShortEdges_->getEntity(i)->getNumNodes(); j++) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        tmpVector->replaceLocalValue(ShortEdges_->getEntity(i)->getGammaDofID(j,k),i,1.0);
                    }
                }
            }
            
            this->LocalPartitionOfUnity_[1] = tmpVector;
        }
        
        if (UseStraightEdges_ && StraightEdges_->getNumEntities()>0) {
            MultiVectorPtr tmpVector = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(serialInterfaceMap,StraightEdges_->getNumEntities());
            
            for (UN i=0; i<StraightEdges_->getNumEntities(); i++) {
                for (UN j=0; j<StraightEdges_->getEntity(i)->getNumNodes(); j++) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        tmpVector->replaceLocalValue(StraightEdges_->getEntity(i)->getGammaDofID(j,k),i,1.0);
                    }
                }
            }
            this->LocalPartitionOfUnity_[2] = tmpVector;
        }
        
        if (UseEdges_ && Edges_->getNumEntities()>0) {
            MultiVectorPtr tmpVector = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(serialInterfaceMap,Edges_->getNumEntities());
            
            for (UN i=0; i<Edges_->getNumEntities(); i++) {
                for (UN j=0; j<Edges_->getEntity(i)->getNumNodes(); j++) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        tmpVector->replaceLocalValue(Edges_->getEntity(i)->getGammaDofID(j,k),i,1.0);
                    }
                }
            }
            
            this->LocalPartitionOfUnity_[3] = tmpVector;
        }
        
        if (UseFaces_ && Faces_->getNumEntities()>0) {
            MultiVectorPtr tmpVector = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(serialInterfaceMap,Faces_->getNumEntities());
            
            for (UN i=0; i<Faces_->getNumEntities(); i++) {
                for (UN j=0; j<Faces_->getEntity(i)->getNumNodes(); j++) {
                    for (UN k=0; k<dofsPerNode; k++) {
                        tmpVector->replaceLocalValue(Faces_->getEntity(i)->getGammaDofID(j,k),i,1.0);
                    }
                }
            }
            
            this->LocalPartitionOfUnity_[4] = tmpVector;
        }
        
        return 0;
    }
    
}

#endif
