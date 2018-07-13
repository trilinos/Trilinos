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

#ifndef _FROSCH_INTERFACEENTITY_DECL_HPP
#define _FROSCH_INTERFACEENTITY_DECL_HPP

#define FROSCH_ASSERT(A,S) if(!(A)) { std::cerr<<"Assertion failed. "<<S<<std::endl; std::cout.flush(); throw std::out_of_range("Assertion.");};

#include <Xpetra_VectorFactory_fwd.hpp>

#include <FROSch_ExtractSubmatrices_def.hpp>
#include <FROSch_Tools_def.hpp>


namespace FROSch {
    
    template <class SC = Xpetra::Operator<>::scalar_type,
    class LO = typename Xpetra::Operator<SC>::local_ordinal_type,
    class GO = typename Xpetra::Operator<SC,LO>::global_ordinal_type,
    class NO = typename Xpetra::Operator<SC,LO,GO>::node_type>
    class EntitySet;
    
    enum EntityType {VertexType,ShortEdgeType,StraightEdgeType,EdgeType,FaceType,SurfaceType,VolumeType};
    
    template <class SC = Xpetra::Operator<>::scalar_type,
    class LO = typename Xpetra::Operator<SC>::local_ordinal_type,
    class GO = typename Xpetra::Operator<SC,LO>::global_ordinal_type>
    struct Node {
        LO NodeIDGamma_;
        LO NodeIDLocal_;
        GO NodeIDGlobal_;
        
        Teuchos::ArrayRCP<LO> DofsGamma_;
        Teuchos::ArrayRCP<LO> DofsLocal_;
        Teuchos::ArrayRCP<GO> DofsGlobal_;
        
        bool operator< (const Node &n) const;
        
        bool operator== (const Node &n) const;
    };
    
    template <class SC = Xpetra::Operator<>::scalar_type,
    class LO = typename Xpetra::Operator<SC>::local_ordinal_type,
    class GO = typename Xpetra::Operator<SC,LO>::global_ordinal_type,
    class NO = typename Xpetra::Operator<SC,LO,GO>::node_type>
    class InterfaceEntity {
        
    public:
        
        typedef Xpetra::Matrix<SC,LO,GO,NO> CrsMatrix;
        typedef Teuchos::RCP<CrsMatrix> CrsMatrixPtr;
        
        typedef Xpetra::Vector<SC,LO,GO,NO> Vector;
        typedef Teuchos::RCP<Vector> VectorPtr;
        
        typedef Teuchos::RCP<EntitySet<SC,LO,GO,NO> > EntitySetPtr;
        
        typedef Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> > InterfaceEntityPtr;
        
        typedef Teuchos::Array<Node<SC,LO,GO> > NodeVec;
        typedef Teuchos::RCP<Node<SC,LO,GO> > NodePtr;
        typedef Teuchos::Array<NodePtr> NodePtrVec;
        
        typedef unsigned UN;
        
        typedef Teuchos::ArrayRCP<LO> LOVecPtr;
        
        typedef Teuchos::Array<GO> GOVec;
        typedef Teuchos::ArrayRCP<GO> GOVecPtr;
        
        
        InterfaceEntity(EntityType type,
                        UN dofsPerNode,
                        UN multiplicity,
                        GO *subdomains);
        
        ~InterfaceEntity();
        
        int addNode(LO nodeIDGamma,
                    LO nodeIDLocal,
                    GO nodeIDGlobal,
                    UN nDofs,
                    const LOVecPtr dofsGamma,
                    const LOVecPtr dofsLocal,
                    const GOVecPtr dofsGlobal);
        
        int addNode(const NodePtr &node);
        
        int addNode(const Node<SC,LO,GO> &node);
        
        int resetGlobalDofs(UN iD,
                            UN nDofs,
                            UN *dofIDs,
                            GO *dofsGlobal);
        
        int removeNode(UN iD);
        
        int sortByGlobalID();
        
        int setUniqueID(GO uniqueID);
        
        int setLocalID(LO localID);
        
        int setAncestorID(LO AncestorID);
        
        int setUniqueIDToFirstGlobalID();
        
        int resetEntityType(EntityType type);
        
        int findAncestors(EntitySetPtr entitySet);
        
        InterfaceEntityPtr divideEntity(CrsMatrixPtr matrix, int pID);
        
        /////////////////
        // Get Methods //
        /////////////////
        
        EntityType getEntityType() const;
        
        UN getDofsPerNode() const;
        
        UN getMultiplicity() const;
        
        GO getUniqueID() const;
        
        LO getLocalID() const;
        
        LO getAncestorID() const;
        
        const Node<SC,LO,GO>& getNode(UN iDNode) const;
        
        LO getGammaNodeID(UN iDNode) const;
        
        LO getLocalNodeID(UN iDNode) const;
        
        GO getGlobalNodeID(UN iDNode) const;
        
        LO getGammaDofID(UN iDNode, UN iDDof) const;
        
        LO getLocalDofID(UN iDNode, UN iDDof) const;
        
        GO getGlobalDofID(UN iDNode, UN iDDof) const;
        
        const GOVec & getSubdomainsVector() const;
        
        UN getNumNodes() const;
        
        const EntitySetPtr getAncestors() const;
        
    protected:
        
        EntityType Type_;
        
        NodeVec NodeVector_;
        
        GOVec SubdomainsVector_;
        
        EntitySetPtr Ancestors_;
        
        UN DofsPerNode_;
        UN Multiplicity_;
        GO UniqueID_;
        LO LocalID_;
        LO AncestorID_;
    };
    
    template <class SC,class LO,class GO,class NO>
    bool compareInterfaceEntities(Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> > iEa, Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> > iEb);
    
    template <class SC,class LO,class GO,class NO>
    bool equalInterfaceEntities(Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> > iEa, Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> > iEb);
    
}

#endif
