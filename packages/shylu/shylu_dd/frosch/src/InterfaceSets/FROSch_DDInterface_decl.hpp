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

#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_ExportFactory_fwd.hpp>

#include <FROSch_EntitySet_def.hpp>
#include <FROSch_InterfaceEntity_decl.hpp>

#include <FROSch_ExtractSubmatrices_def.hpp>

// TODO
// -> "Parent" -> "Anchestor"

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
        typedef Teuchos::ArrayRCP<MapPtr> MapPtrVecPtr;
        
        typedef Xpetra::Matrix<SC,LO,GO,NO> CrsMatrix;
        typedef Teuchos::RCP<CrsMatrix> CrsMatrixPtr;
        
        typedef Teuchos::RCP<EntitySet<SC,LO,GO,NO> > EntitySetPtr;
        
        typedef Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> > InterfaceEntityPtr;
        typedef Teuchos::ArrayRCP<InterfaceEntityPtr> InterfaceEntityPtrVecPtr;
        
        typedef unsigned UN;
        typedef Teuchos::ArrayRCP<UN> UNVecPtr;
        
        typedef Teuchos::ArrayRCP<LO> LOVecPtr;
        
        typedef Teuchos::Array<GO> GOVec;
        typedef Teuchos::ArrayRCP<GO> GOVecPtr;
        typedef Teuchos::ArrayView<GO> GOVecView;
        typedef Teuchos::Array<GOVec> GOVecVec;
        typedef Teuchos::ArrayRCP<GOVec> GOVecVecPtr;
        
        typedef Teuchos::ArrayRCP<SC> SCVecPtr;
        typedef Teuchos::ArrayRCP<SCVecPtr> SCVecPtr2D;
        
        
        DDInterface(UN dimension,
                    UN dofsPerNode,
                    MapPtr localToGlobalMap);
        
        ~DDInterface();
        
        int resetGlobalDofs(MapPtrVecPtr dofsMaps);
        
        int removeDirichletNodes(GOVecView myGlobalDirichletBoundaryDofs);
        
        int divideUnconnectedEntities(CrsMatrixPtr matrix);
        
        int sortEntities();
        
        int sortEntities(SCVecPtr2D localNodeList);
        
        int findParents();
        
        EntitySetPtr & getVertices();
        
        EntitySetPtr & getShortEdges();
        
        EntitySetPtr & getStraightEdges();
        
        EntitySetPtr & getEdges();
        
        EntitySetPtr & getFaces();
        
        EntitySetPtr & getInterface();
        
        EntitySetPtr & getInterior();
        
        EntitySetPtr & getParentVertices();
        
        EntitySetPtr & getParentEdges();
        
        EntitySetPtr & getParentFaces();
        
    protected:
        
        int communicateLocalComponents(GOVecVecPtr &componentsSubdomains,
                                       GOVecVec &componentsSubdomainsUnique);
        
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
        EntitySetPtr ParentVertices_;
        EntitySetPtr ParentEdges_;
        EntitySetPtr ParentFaces_;
        
        MapPtr LocalToGlobalNodesMap_;
        MapPtr LocalToGlobalNodesUniqueMap_;
    };
    
}

#endif
