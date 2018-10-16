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

#ifndef _FROSCH_INTERFACEENTITY_DEF_HPP
#define _FROSCH_INTERFACEENTITY_DEF_HPP

#include <FROSch_InterfaceEntity_decl.hpp>

namespace FROSch {
    
    template <class SC,class LO,class GO>
    bool Node<SC,LO,GO>::operator< (const Node &n) const
    {
        return NodeIDGlobal_<n.NodeIDGlobal_;
    }
    
    template <class SC,class LO,class GO>
    bool Node<SC,LO,GO>::operator== (const Node &n) const
    {
        return NodeIDGlobal_==n.NodeIDGlobal_;
    }
    
    
    template <class SC,class LO,class GO,class NO>
    InterfaceEntity<SC,LO,GO,NO>::InterfaceEntity(EntityType type,
                                                  UN dofsPerNode,
                                                  UN multiplicity,
                                                  GO *subdomains) :
    Type_ (type),
    NodeVector_ (0),
    SubdomainsVector_ (multiplicity),
    Ancestors_ (),
    DofsPerNode_ (dofsPerNode),
    Multiplicity_ (multiplicity),
    UniqueID_ (-1),
    LocalID_ (-1),
    AncestorID_ (-1)
    {
        for (UN i=0; i<multiplicity; i++) {
            SubdomainsVector_[i] = subdomains[i];
        }
        sortunique(SubdomainsVector_);
        if ((Type_==ShortEdgeType)||(Type_==StraightEdgeType)||(Type_==EdgeType)) {
            Ancestors_.reset(new EntitySet<SC,LO,GO,NO>(VertexType));
        } else if (Type_==FaceType) {
            Ancestors_.reset(new EntitySet<SC,LO,GO,NO>(EdgeType));
        }
    }
    
    template <class SC,class LO,class GO,class NO>
    InterfaceEntity<SC,LO,GO,NO>::~InterfaceEntity()
    {
        
    } // Do we need sth here?
    
    template <class SC,class LO,class GO,class NO>
    int InterfaceEntity<SC,LO,GO,NO>::addNode(LO nodeIDGamma,
                                              LO nodeIDLocal,
                                              GO nodeIDGlobal,
                                              UN nDofs,
                                              const LOVecPtr dofsGamma,
                                              const LOVecPtr dofsLocal,
                                              const GOVecPtr dofsGlobal)
    {
        FROSCH_ASSERT(nDofs<=DofsPerNode_,"nDofs>NumDofs_.");
        
        FROSCH_ASSERT(dofsGamma.size()==nDofs,"dofIDs.size()!=nDofs");
        FROSCH_ASSERT(dofsLocal.size()==nDofs,"dofIDs.size()!=nDofs");
        FROSCH_ASSERT(dofsGlobal.size()==nDofs,"dofIDs.size()!=nDofs");

        Node<SC,LO,GO> node;
        
        node.NodeIDGamma_ = nodeIDGamma;
        node.NodeIDLocal_ = nodeIDLocal;
        node.NodeIDGlobal_ = nodeIDGlobal;
        
        node.DofsGamma_.deepCopy(dofsGamma());
        node.DofsLocal_.deepCopy(dofsLocal());
        node.DofsGlobal_.deepCopy(dofsGlobal());
        
        NodeVector_.push_back(node);
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int InterfaceEntity<SC,LO,GO,NO>::addNode(const NodePtr &node)
    {
        return addNode(*node);
    }
    
    template <class SC,class LO,class GO,class NO>
    int InterfaceEntity<SC,LO,GO,NO>::addNode(const Node<SC,LO,GO> &node)
    {
        FROSCH_ASSERT(node.DofsGamma_.size()<=DofsPerNode_,"node.DofsGamma_ is too large.");
        FROSCH_ASSERT(node.DofsLocal_.size()<=DofsPerNode_,"node.DofsLocal_ is too large.");
        FROSCH_ASSERT(node.DofsGlobal_.size()<=DofsPerNode_,"node.DofsGlobal_ is too large.");
        
        NodeVector_.push_back(node);
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int InterfaceEntity<SC,LO,GO,NO>::resetGlobalDofs(UN iD,
                                                      UN nDofs,
                                                      UN *dofIDs,
                                                      GO *dofsGlobal)
    {
        FROSCH_ASSERT(iD<getNumNodes(),"iD=>getNumNodes()");
        FROSCH_ASSERT(nDofs<=DofsPerNode_,"nDofs>DofsPerNode_.");
        
        for (unsigned i=0; i<nDofs; i++) {
            if ((0<=dofIDs[i])&&(dofIDs[i]<=DofsPerNode_)) {
                NodeVector_[iD].DofsGlobal_[dofIDs[i]] = dofsGlobal[dofIDs[i]];
            } else {
                FROSCH_ASSERT(0!=0,"dofIDs[i] is out of range.");
            }
        }
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int InterfaceEntity<SC,LO,GO,NO>::removeNode(UN iD)
    {
        FROSCH_ASSERT(iD<getNumNodes(),"iD=>getNumNodes()");
        NodeVector_.erase(NodeVector_.begin()+iD);
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int InterfaceEntity<SC,LO,GO,NO>::sortByGlobalID()
    {
        sortunique(NodeVector_);
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int InterfaceEntity<SC,LO,GO,NO>::setUniqueID(GO uniqueID)
    {
        UniqueID_ = uniqueID;
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int InterfaceEntity<SC,LO,GO,NO>::setLocalID(LO localID)
    {
        LocalID_ = localID;
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int InterfaceEntity<SC,LO,GO,NO>::setAncestorID(LO AncestorID)
    {
        AncestorID_ = AncestorID;
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int InterfaceEntity<SC,LO,GO,NO>::setUniqueIDToFirstGlobalID()
    {
        UniqueID_ = NodeVector_[0].NodeIDGlobal_;
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int InterfaceEntity<SC,LO,GO,NO>::resetEntityType(EntityType type)
    {
        Type_ = type;
        if ((Type_==ShortEdgeType)||(Type_==StraightEdgeType)||(Type_==EdgeType)) {
            Ancestors_.reset(new EntitySet<SC,LO,GO,NO>(VertexType));
        } else if (Type_==FaceType) {
            Ancestors_.reset(new EntitySet<SC,LO,GO,NO>(EdgeType));
        }
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int InterfaceEntity<SC,LO,GO,NO>::findAncestors(EntitySetPtr entitySet)
    {
        if (Type_ == VertexType) {
            FROSCH_ASSERT(0!=0,"There are no Ancestors to vertices.")
        } else if ((Type_ == ShortEdgeType) || (Type_ == StraightEdgeType) || (Type_ == EdgeType)) {
            FROSCH_ASSERT(entitySet->getEntityType()==VertexType,"entitySet has the wrong type.")
        } else if (Type_ == FaceType) {
            FROSCH_ASSERT((entitySet->getEntityType()==ShortEdgeType)||(entitySet->getEntityType()==StraightEdgeType)||(entitySet->getEntityType()==EdgeType),"entitySet has the wrong type.")
        } else if (Type_ == SurfaceType) {
            FROSCH_ASSERT(0!=0,"SurfaceType not yet implemented...")
        } else if (Type_ == VolumeType) {
            FROSCH_ASSERT(0!=0,"VolumeType not yet implemented...")
        }
        
        //
        EntitySetPtr Ancestors(new EntitySet<SC,LO,GO,NO>(*entitySet));
        GOVec tmpVector;
        for (UN i=0; i<Multiplicity_; i++) {
            UN length = Ancestors->getNumEntities();
            for (UN j=0; j<length; j++) {
                tmpVector = Ancestors->getEntity(length-1-j)->getSubdomainsVector();
                if (!std::binary_search(tmpVector.begin(),tmpVector.end(),SubdomainsVector_[i])) {
                    Ancestors->removeEntity(length-1-j);
                }
            }
        }
        
        for (UN i=0; i<Ancestors->getNumEntities(); i++) {
            Ancestors_->addEntity(Ancestors->getEntity(i));
        }
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    typename InterfaceEntity<SC,LO,GO,NO>::InterfaceEntityPtr InterfaceEntity<SC,LO,GO,NO>::divideEntity(CrsMatrixPtr matrix,int pID)
    {
        InterfaceEntityPtr entity(new InterfaceEntity<SC,LO,GO,NO>(Type_,DofsPerNode_,Multiplicity_,&(SubdomainsVector_[0])));
        
        if (getNumNodes()>=2) {
            sortByGlobalID();
            
            GOVecPtr mapVector(getNumNodes());
            for (UN i=0; i<getNumNodes(); i++) {
                mapVector[i] = NodeVector_[i].DofsGamma_[0];
            }
            CrsMatrixPtr localMatrix,mat1,mat2,mat3;
            BuildSubmatrices(matrix,mapVector(),localMatrix,mat1,mat2,mat3);
            
            VectorPtr iterationVector = Xpetra::VectorFactory<SC,LO,GO,NO>::Build(localMatrix->getRowMap());
            iterationVector->getDataNonConst(0)[0] = 1.0;
            for (UN i=0; i<getNumNodes()-1; i++) {
                localMatrix->apply(*iterationVector,*iterationVector);
            }
            
            for (UN i=0; i<getNumNodes(); i++) {
                if (fabs(iterationVector->getData(0)[i])<1.0e-10) {
                    entity->addNode(getNode(i));
                    removeNode(i);
                }
            }
        }
        return entity;
    }
    
    /////////////////
    // Get Methods //
    /////////////////
    
    template <class SC,class LO,class GO,class NO>
    EntityType InterfaceEntity<SC,LO,GO,NO>::getEntityType() const
    {
        return Type_;
    }
    
    template <class SC,class LO,class GO,class NO>
    typename InterfaceEntity<SC,LO,GO,NO>::UN InterfaceEntity<SC,LO,GO,NO>::getDofsPerNode() const
    {
        return DofsPerNode_;
    }
    
    template <class SC,class LO,class GO,class NO>
    typename InterfaceEntity<SC,LO,GO,NO>::UN InterfaceEntity<SC,LO,GO,NO>::getMultiplicity() const
    {
        return Multiplicity_;
    }
    
    template <class SC,class LO,class GO,class NO>
    GO InterfaceEntity<SC,LO,GO,NO>::getUniqueID() const
    {
        return UniqueID_;
    }
    
    template <class SC,class LO,class GO,class NO>
    LO InterfaceEntity<SC,LO,GO,NO>::getLocalID() const
    {
        return LocalID_;
    }
    
    template <class SC,class LO,class GO,class NO>
    LO InterfaceEntity<SC,LO,GO,NO>::getAncestorID() const
    {
        return AncestorID_;
    }
    
    template <class SC,class LO,class GO,class NO>
    const Node<SC,LO,GO>& InterfaceEntity<SC,LO,GO,NO>::getNode(UN iDNode) const
    {
        return NodeVector_[iDNode];
    }
    
    template <class SC,class LO,class GO,class NO>
    LO InterfaceEntity<SC,LO,GO,NO>::getGammaNodeID(UN iDNode) const
    {
        return NodeVector_[iDNode].NodeIDGamma_;
    }
    
    template <class SC,class LO,class GO,class NO>
    LO InterfaceEntity<SC,LO,GO,NO>::getLocalNodeID(UN iDNode) const
    {
        return NodeVector_[iDNode].NodeIDLocal_;
    }
    
    template <class SC,class LO,class GO,class NO>
    GO InterfaceEntity<SC,LO,GO,NO>::getGlobalNodeID(UN iDNode) const
    {
        return NodeVector_[iDNode].NodeIDGlobal_;
    }
    
    template <class SC,class LO,class GO,class NO>
    LO InterfaceEntity<SC,LO,GO,NO>::getGammaDofID(UN iDNode, UN iDDof) const
    {
        return NodeVector_[iDNode].DofsGamma_[iDDof];
    }
    
    template <class SC,class LO,class GO,class NO>
    LO InterfaceEntity<SC,LO,GO,NO>::getLocalDofID(UN iDNode, UN iDDof) const
    {
        return NodeVector_[iDNode].DofsLocal_[iDDof];
    }
    
    template <class SC,class LO,class GO,class NO>
    GO InterfaceEntity<SC,LO,GO,NO>::getGlobalDofID(UN iDNode, UN iDDof) const
    {
        return NodeVector_[iDNode].DofsGlobal_[iDDof];
    }
    
    template <class SC,class LO,class GO,class NO>
    const typename InterfaceEntity<SC,LO,GO,NO>::GOVec & InterfaceEntity<SC,LO,GO,NO>::getSubdomainsVector() const
    {
        return SubdomainsVector_;
    }
    
    template <class SC,class LO,class GO,class NO>
    typename InterfaceEntity<SC,LO,GO,NO>::UN InterfaceEntity<SC,LO,GO,NO>::getNumNodes() const
    {
        return NodeVector_.size();
    }
    
    template <class SC,class LO,class GO,class NO>
    const typename InterfaceEntity<SC,LO,GO,NO>::EntitySetPtr InterfaceEntity<SC,LO,GO,NO>::getAncestors() const
    {
        return Ancestors_;
    }
    
    
    template <class SC,class LO,class GO,class NO>
    bool compareInterfaceEntities(Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> > iEa, Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> > iEb)
    {
        return iEa->getUniqueID()<iEb->getUniqueID();
    }
    
    template <class SC,class LO,class GO,class NO>
    bool equalInterfaceEntities(Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> > iEa, Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> > iEb)
    {
        return iEa->getUniqueID()==iEb->getUniqueID();
    }
    
}

#endif
