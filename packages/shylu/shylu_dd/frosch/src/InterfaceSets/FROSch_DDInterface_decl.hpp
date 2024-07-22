// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_DDINTERFACE_DECL_HPP
#define _FROSCH_DDINTERFACE_DECL_HPP

//#define INTERFACE_OUTPUT

#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_ExportFactory_fwd.hpp>
#include <Xpetra_CrsGraphFactory.hpp>

#include <FROSch_EntitySet_def.hpp>
#include <FROSch_InterfaceEntity_decl.hpp>

#include <FROSch_ExtractSubmatrices_def.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    enum CommunicationStrategy {CommCrsMatrix,CommCrsGraph,CreateOneToOneMap};

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class DDInterface {

    protected:

        using CommPtr                   = RCP<const Comm<int> >;

        using XMap                      = Map<LO,GO,NO>;
        using XMapPtr                   = RCP<XMap>;
        using ConstXMapPtr              = RCP<const XMap>;
        using XMapPtrVecPtr             = ArrayRCP<XMapPtr>;
        using ConstXMapPtrVecPtr        = ArrayRCP<ConstXMapPtr>;

        using XMatrix                   = Matrix<SC,LO,GO,NO>;
        using XMatrixPtr                = RCP<XMatrix>;
        using ConstXMatrixPtr           = RCP<const XMatrix>;

        using XCrsGraph                 = CrsGraph<LO,GO,NO>;
        using XCrsGraphPtr              = RCP<XCrsGraph>;

        using XMultiVector              = MultiVector<SC,LO,GO,NO>;
        using XMultiVectorPtr           = RCP<XMultiVector>;
        using ConstXMultiVectorPtr      = RCP<const XMultiVector>;

        using XImport                   = Import<LO,GO,NO>;
        using XImportPtr                = RCP<XImport>;

        using XExport                   = Export<LO,GO,NO>;
        using XExportPtr                = RCP<XExport>;

        using EntitySetPtr              = RCP<EntitySet<SC,LO,GO,NO> >;
        using EntitySetConstPtr         = const EntitySetPtr;
        using EntitySetPtrVecPtr        = ArrayRCP<EntitySetPtr>;
        using EntitySetPtrConstVecPtr   = const EntitySetPtrVecPtr;

        using EntityFlagVecPtr          = ArrayRCP<EntityFlag>;

        using InterfaceEntityPtr        = RCP<InterfaceEntity<SC,LO,GO,NO> >;
        using InterfaceEntityPtrVecPtr  = ArrayRCP<InterfaceEntityPtr>;

        using UN                        = unsigned;
        using ConstUN                   = const UN;
        using UNVecPtr                  = ArrayRCP<UN>;

        using IntVec                    = Array<int>;
        using IntVecVec                 = Array<IntVec>;
        using IntVecVecPtr              = ArrayRCP<IntVec>;

        using LOVec                     = Array<LO>;
        using LOVecPtr                  = ArrayRCP<LO>;

        using GOVec                     = Array<GO>;
        using ConstGOVecView            = ArrayView<const GO>;
        using GOVecPtr                  = ArrayRCP<GO>;
        using GOVecView                 = ArrayView<GO>;
        using GOVecVec                  = Array<GOVec>;
        using GOVecVecPtr               = ArrayRCP<GOVec>;

        using SCVec                     = Array<SC>;
        using SCVecPtr                  = ArrayRCP<SC>;

    public:

        DDInterface(UN dimension,
                    UN dofsPerNode,
                    ConstXMapPtr localToGlobalMap,
                    Verbosity verbosity = All,
                    UN levelID = 1,
                    CommunicationStrategy commStrategy = CommCrsGraph);

        ~DDInterface();

        int resetGlobalDofs(ConstXMapPtrVecPtr dofsMaps);

        int removeDirichletNodes(GOVecView dirichletBoundaryDofs);

        int divideUnconnectedEntities(ConstXMatrixPtr matrix);

        int flagEntities(ConstXMultiVectorPtr nodeList = null);

        int removeEmptyEntities();

        int sortVerticesEdgesFaces(ConstXMultiVectorPtr nodeList = null);

        int buildEntityMaps(bool buildVerticesMap = true,
                            bool buildShortEdgesMap = true,
                            bool buildStraightEdgesMap = true,
                            bool buildEdgesMap = true,
                            bool buildFacesMap = true,
                            bool buildRootsMap = false,
                            bool buildLeafsMap = false);

        int buildEntityHierarchy();

        int computeDistancesToRoots(UN dimension,
                                    ConstXMultiVectorPtr &nodeList = null,
                                    DistanceFunction distanceFunction = ConstantDistanceFunction);

        //! This function extracts those entities which are to be used to build a connectivity graph on the subdomain
        //! level. By default, we identify all entities with multiplicity 2. Afterwards, the corresponding entities can
        //! be obtained using the function getConnectivityEntities().
        //! If short or straight edges should be omitted, the function flagEntities() has to be called in advance.
        int identifyConnectivityEntities(UNVecPtr multiplicities = null,
                                         EntityFlagVecPtr flags = null);

        UN getDimension() const;

        UN getDofsPerNode() const;

        LO getNumMyNodes() const;

        //
        // Remove the references below?
        //

        EntitySetConstPtr & getVertices() const;

        EntitySetConstPtr & getShortEdges() const;

        EntitySetConstPtr & getStraightEdges() const;

        EntitySetConstPtr & getEdges() const;

        EntitySetConstPtr & getFaces() const;

        EntitySetConstPtr & getInterface() const;

        EntitySetConstPtr & getInterior() const;

        EntitySetConstPtr & getRoots() const;

        EntitySetConstPtr & getLeafs() const;

        EntitySetPtrConstVecPtr & getEntitySetVector() const;

        GOVec getNumEnt() const;
        //! This function returns those entities which are to be used to build a connectivity graph on the subdomain
        //! level. They have to identified first using the function identifyConnectivityEntities().
        EntitySetConstPtr & getConnectivityEntities() const;

        ConstXMapPtr getNodesMap() const;


    protected:

        int communicateLocalComponents(IntVecVecPtr &componentsSubdomains,
                                       IntVecVec &componentsSubdomainsUnique);

        int identifyLocalComponents(IntVecVecPtr &componentsSubdomains,
                                    IntVecVec &componentsSubdomainsUnique);


        CommPtr MpiComm_;

        UN Dimension_ = 3;
        UN DofsPerNode_ = 1;
        LO NumMyNodes_ = 0;

        EntitySetPtr Vertices_ = EntitySetPtr(new EntitySet<SC,LO,GO,NO>(VertexType));
        EntitySetPtr ShortEdges_ = EntitySetPtr(new EntitySet<SC,LO,GO,NO>(EdgeType));
        EntitySetPtr StraightEdges_ = EntitySetPtr(new EntitySet<SC,LO,GO,NO>(EdgeType));
        EntitySetPtr Edges_ = EntitySetPtr(new EntitySet<SC,LO,GO,NO>(EdgeType));
        EntitySetPtr Faces_ = EntitySetPtr(new EntitySet<SC,LO,GO,NO>(FaceType));
        EntitySetPtr Interface_ = EntitySetPtr(new EntitySet<SC,LO,GO,NO>(InterfaceType));
        EntitySetPtr Interior_ = EntitySetPtr(new EntitySet<SC,LO,GO,NO>(InteriorType));
        EntitySetPtr Roots_ = EntitySetPtr(new EntitySet<SC,LO,GO,NO>(DefaultType));
        EntitySetPtr Leafs_ = EntitySetPtr(new EntitySet<SC,LO,GO,NO>(DefaultType));
        EntitySetPtr ConnectivityEntities_ = EntitySetPtr(new EntitySet<SC,LO,GO,NO>(DefaultType));
        EntitySetPtrVecPtr EntitySetVector_;

        ConstXMapPtr NodesMap_;
        ConstXMapPtr UniqueNodesMap_;

        CommunicationStrategy CommStrategy_ = CommCrsGraph;

        bool Verbose_ = false;

        Verbosity Verbosity_ = All;

        ConstUN LevelID_ = 1;
        GOVec NumEntity_;
    };

}

#endif
