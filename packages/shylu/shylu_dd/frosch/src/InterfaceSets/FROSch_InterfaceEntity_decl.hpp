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

#include <Xpetra_VectorFactory_fwd.hpp>

#include <FROSch_ExtractSubmatrices_def.hpp>
#include <FROSch_Tools_def.hpp>


namespace FROSch {
    
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = KokkosClassic::DefaultNode::DefaultNodeType>
    class EntitySet;

    enum EntityType {DefaultType,VertexType,EdgeType,FaceType,InteriorType,InterfaceType};
    enum EntityFlag {DefaultFlag,StraightFlag,ShortFlag,NodeFlag};
    enum DistanceFunction {ConstantDistanceFunction,InverseEuclideanDistanceFunction};

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal>
    struct Node {
        LO NodeIDGamma_;
        LO NodeIDLocal_;
        GO NodeIDGlobal_;

        ArrayRCP<LO> DofsGamma_;
        ArrayRCP<LO> DofsLocal_;
        ArrayRCP<GO> DofsGlobal_;

        bool operator< (const Node &n) const;

        bool operator== (const Node &n) const;
    };

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = KokkosClassic::DefaultNode::DefaultNodeType>
    class InterfaceEntity {

    protected:

        using XMatrix               = Matrix<SC,LO,GO,NO>;
        using XMatrixPtr            = RCP<XMatrix>;
        using ConstXMatrixPtr       = RCP<const XMatrix>;

        using XVector               = Vector<SC,LO,GO,NO>;
        using XVectorPtr            = RCP<XVector>;

        using XMultiVector          = MultiVector<SC,LO,GO,NO>;
        using XMultiVectorPtr       = RCP<XMultiVector>;
        using ConstXMultiVectorPtr  = RCP<const XMultiVector>;

        using EntitySetPtr          = RCP<EntitySet<SC,LO,GO,NO> >;

        using InterfaceEntityPtr    = RCP<InterfaceEntity<SC,LO,GO,NO> >;

        using NodeVec               = Array<Node<SC,LO,GO> >;
        using NodePtr               = RCP<Node<SC,LO,GO> >;
        using NodePtrVec            = Array<NodePtr>;

        using UN                    = unsigned;

        using IntVec                = Array<int>;

        using LOVecPtr              = ArrayRCP<LO>;

        using GOVec                 = Array<GO>;
        using GOVecPtr              = ArrayRCP<GO>;

        using SCVecPtr              = ArrayRCP<SC>;
        using SCVecPtrVec           = Array<SCVecPtr>;

    public:

        InterfaceEntity(EntityType type,
                        UN dofsPerNode,
                        UN multiplicity,
                        const int *subdomains,
                        EntityFlag flag = DefaultFlag);

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

        int setCoarseNodeID(LO coarseNodeID);

        int setUniqueIDToFirstGlobalID();

        int resetEntityType(EntityType type);

        int resetEntityFlag(EntityFlag flag);

        int findAncestorsInSet(EntitySetPtr entitySet);

        int clearAncestors();

        int addOffspring(InterfaceEntityPtr interfaceEntity);

        int clearOffspring();

        EntitySetPtr findCoarseNodes();

        int clearCoarseNodes();

        int computeDistancesToCoarseNodes(UN dimension,
                                          ConstXMultiVectorPtr &nodeList = null,
                                          DistanceFunction distanceFunction = ConstantDistanceFunction);

        InterfaceEntityPtr divideEntity(ConstXMatrixPtr matrix,
                                        int pID);

        /////////////////
        // Get Methods //
        /////////////////

        EntityType getEntityType() const;

        EntityFlag getEntityFlag() const;

        UN getDofsPerNode() const;

        UN getMultiplicity() const;

        GO getUniqueID() const;

        LO getLocalID() const;

        LO getCoarseNodeID() const;

        const Node<SC,LO,GO>& getNode(UN iDNode) const;

        LO getGammaNodeID(UN iDNode) const;

        LO getLocalNodeID(UN iDNode) const;

        GO getGlobalNodeID(UN iDNode) const;

        LO getGammaDofID(UN iDNode, UN iDDof) const;

        LO getLocalDofID(UN iDNode, UN iDDof) const;

        GO getGlobalDofID(UN iDNode, UN iDDof) const;

        const IntVec & getSubdomainsVector() const;

        UN getNumNodes() const;

        const EntitySetPtr getAncestors() const;

        const EntitySetPtr getOffspring() const;

        const EntitySetPtr getCoarseNodes() const;

        SC getDistanceToCoarseNode(UN iDNode,
                                   UN iDCoarseNode) const;

    protected:

        EntityType Type_;

        EntityFlag Flag_;

        NodeVec NodeVector_;

        IntVec SubdomainsVector_;

        EntitySetPtr Ancestors_;
        EntitySetPtr Offspring_;
        EntitySetPtr CoarseNodes_;

        SCVecPtrVec DistancesVector_; // AH 08/08/2019 TODO: make a MultiVector out of this

        UN DofsPerNode_;
        UN Multiplicity_;
        GO UniqueID_;
        LO LocalID_;
        LO CoarseNodeID_;
    };

    template <class SC,class LO,class GO,class NO>
    bool compareInterfaceEntities(RCP<InterfaceEntity<SC,LO,GO,NO> > iEa,
                                  RCP<InterfaceEntity<SC,LO,GO,NO> > iEb);

    template <class SC,class LO,class GO,class NO>
    bool equalInterfaceEntities(RCP<InterfaceEntity<SC,LO,GO,NO> > iEa,
                                RCP<InterfaceEntity<SC,LO,GO,NO> > iEb);

}

#endif
