// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
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
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class InterfaceEntity {

    protected:

        using XMatrix               = Xpetra::Matrix<SC,LO,GO,NO>;
        using XMatrixPtr            = RCP<XMatrix>;
        using ConstXMatrixPtr       = RCP<const XMatrix>;

        using XVector               = Xpetra::Vector<SC,LO,GO,NO>;
        using XVectorPtr            = RCP<XVector>;

        using XMultiVector          = Xpetra::MultiVector<SC,LO,GO,NO>;
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

        int setRootID(LO rootID);
        
        int setLeafID(LO leafID);

        int setUniqueIDToFirstGlobalID();

        int resetEntityType(EntityType type);

        int resetEntityFlag(EntityFlag flag);

        int findAncestorsInSet(EntitySetPtr entitySet);

        int clearAncestors();

        int addOffspring(InterfaceEntityPtr interfaceEntity);

        int clearOffspring();

        EntitySetPtr findRoots();

        int clearRoots();

        int computeDistancesToRoots(UN dimension,
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

        LO getRootID() const;
        
        LO getLeafID() const;

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

        const EntitySetPtr getRoots() const;

        SC getDistanceToRoot(UN iDNode,
                             UN iDRoot) const;

    protected:

        EntityType Type_ = DefaultType;

        EntityFlag Flag_ = DefaultFlag;

        NodeVec NodeVector_ = NodeVec(0);

        IntVec SubdomainsVector_ = IntVec(0);

        EntitySetPtr Ancestors_;
        EntitySetPtr Offspring_;
        EntitySetPtr Roots_;

        SCVecPtrVec DistancesVector_ = SCVecPtrVec(0); // AH 08/08/2019 TODO: make a MultiVector out of this

        UN DofsPerNode_ = 1;
        UN Multiplicity_ = 1;
        GO UniqueID_ = -1;
        LO LocalID_ = -1;
        LO RootID_ = -1;
        LO LeafID_ = -1;
    };

    template <class SC,class LO,class GO,class NO>
    bool compareInterfaceEntities(RCP<InterfaceEntity<SC,LO,GO,NO> > iEa,
                                  RCP<InterfaceEntity<SC,LO,GO,NO> > iEb);

    template <class SC,class LO,class GO,class NO>
    bool equalInterfaceEntities(RCP<InterfaceEntity<SC,LO,GO,NO> > iEa,
                                RCP<InterfaceEntity<SC,LO,GO,NO> > iEb);

}

#endif
