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

#ifndef _FROSCH_DDINTERFACE_DECL_HPP
#define _FROSCH_DDINTERFACE_DECL_HPP

#define FROSCH_ASSERT(A,S) if(!(A)) { std::cerr<<"Assertion failed. "<<S<<std::endl; std::cout.flush(); throw std::out_of_range("Assertion.");};

//#define INTERFACE_OUTPUT
//#define FROSCH_OFFSET_MAPS

#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_ExportFactory_fwd.hpp>

#include <FROSch_EntitySet_def.hpp>
#include <FROSch_InterfaceEntity_decl.hpp>

#include <FROSch_ExtractSubmatrices_def.hpp>


namespace FROSch {
    
    template <class SC = Xpetra::Operator<>::scalar_type,
    class LO = typename Xpetra::Operator<SC>::local_ordinal_type,
    class GO = typename Xpetra::Operator<SC, LO>::global_ordinal_type,
    class NO = typename Xpetra::Operator<SC, LO, GO>::node_type>
    class DDInterface {
        
    public:
        
        typedef Teuchos::RCP<const Teuchos::Comm<int> > CommPtr;
        
        typedef Xpetra::Map<LO,GO,NO> Map;
        typedef Teuchos::RCP<Map> MapPtr;
        typedef Teuchos::RCP<const Map> ConstMapPtr;
        typedef Teuchos::ArrayRCP<MapPtr> MapPtrVecPtr;
        
        typedef Xpetra::Matrix<SC,LO,GO,NO> CrsMatrix;
        typedef Teuchos::RCP<CrsMatrix> CrsMatrixPtr;
        
        typedef Xpetra::MultiVector<SC,LO,GO,NO> MultiVector;
        typedef Teuchos::RCP<MultiVector> MultiVectorPtr;
        
        typedef Teuchos::RCP<EntitySet<SC,LO,GO,NO> > EntitySetPtr;
        typedef const EntitySetPtr EntitySetConstPtr;
        typedef Teuchos::ArrayRCP<EntitySetPtr> EntitySetPtrVecPtr;
        typedef const EntitySetPtrVecPtr EntitySetPtrConstVecPtr;
        
        typedef Teuchos::ArrayRCP<EntityFlag> EntityFlagVecPtr;
        
        typedef Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> > InterfaceEntityPtr;
        typedef Teuchos::ArrayRCP<InterfaceEntityPtr> InterfaceEntityPtrVecPtr;
        
        typedef unsigned UN;
        typedef Teuchos::ArrayRCP<UN> UNVecPtr;
        
        typedef Teuchos::ArrayRCP<LO> LOVecPtr;
        
        typedef Teuchos::Array<GO> GOVec;
        typedef Teuchos::ArrayView<const GO> ConstGOVecView;
        typedef Teuchos::ArrayRCP<GO> GOVecPtr;
        typedef Teuchos::ArrayView<GO> GOVecView;
        typedef Teuchos::Array<GOVec> GOVecVec;
        typedef Teuchos::ArrayRCP<GOVec> GOVecVecPtr;
        
        typedef Teuchos::ArrayRCP<SC> SCVecPtr;
        
        
        DDInterface(UN dimension,
                    UN dofsPerNode,
                    MapPtr localToGlobalMap);
        
        ~DDInterface();
        
        int resetGlobalDofs(MapPtrVecPtr dofsMaps);
        
        int removeDirichletNodes(GOVecView dirichletBoundaryDofs);
        
        int divideUnconnectedEntities(CrsMatrixPtr matrix);
        
        int flagEntities(MultiVectorPtr nodeList = Teuchos::null);
        
        int removeEmptyEntities();
        
        int sortVerticesEdgesFaces(MultiVectorPtr nodeList = Teuchos::null);
        
        int buildEntityHierarchy();
        
        int computeDistancesToCoarseNodes(UN dimension,
                                          MultiVectorPtr &nodeList = Teuchos::null,
                                          DistanceFunction distanceFunction = ConstantDistanceFunction);
        
        //! This function extracts those entities which are to be used to build a connectivity graph on the subdomain
        //! level. By default, we identify all entities with multiplicity 2. Afterwards, the corresponding entities can
        //! be obtained using the function getConnectivityEntities().
        //! If short or straight edges should be omitted, the function flagEntities() has to be called in advance.
        int identifyConnectivityEntities(UNVecPtr multiplicities = Teuchos::null,
                                         EntityFlagVecPtr flags = Teuchos::null);
        
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
        
        EntitySetConstPtr & getCoarseNodes() const;
        
        EntitySetPtrConstVecPtr & getEntitySetVector() const;
        
        //! This function returns those entities which are to be used to build a connectivity graph on the subdomain
        //! level. They have to identified first using the function identifyConnectivityEntities().
        EntitySetConstPtr & getConnectivityEntities() const;
        
        ConstMapPtr getNodesMap() const;
        
        
    protected:
#ifdef FROSCH_OFFSET_MAPS
        int communicateLocalComponents(GOVecVecPtr &componentsSubdomains,
                                       GOVecVec &componentsSubdomainsUnique,
                                       UN priorDofsPerNode = 0);
#else
        int communicateLocalComponents(GOVecVecPtr &componentsSubdomains,
                                       GOVecVec &componentsSubdomainsUnique);
#endif
        int identifyLocalComponents(GOVecVecPtr &componentsSubdomains,
                                    GOVecVec &componentsSubdomainsUnique);
        
        
        CommPtr MpiComm_;
        
        UN Dimension_;
        UN DofsPerNode_;
        LO NumMyNodes_;
        
        EntitySetPtr Vertices_;
        EntitySetPtr ShortEdges_;
        EntitySetPtr StraightEdges_;
        EntitySetPtr Edges_;
        EntitySetPtr Faces_;
        EntitySetPtr Interface_;
        EntitySetPtr Interior_;
        EntitySetPtr CoarseNodes_;
        EntitySetPtr ConnectivityEntities_;
        EntitySetPtrVecPtr EntitySetVector_;
        
        MapPtr NodesMap_;
        MapPtr UniqueNodesMap_;
    };
    
}

#endif
