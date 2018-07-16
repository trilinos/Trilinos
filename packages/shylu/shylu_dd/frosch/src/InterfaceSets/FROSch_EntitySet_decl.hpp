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

#ifndef _FROSCH_ENTITYSET_DECL_HPP
#define _FROSCH_ENTITYSET_DECL_HPP

#define FROSCH_ASSERT(A,S) if(!(A)) { std::cerr<<"Assertion failed. "<<S<<std::endl; std::cout.flush(); throw std::out_of_range("Assertion.");};

#include <FROSch_InterfaceEntity_def.hpp>

#include <FROSch_Tools_def.hpp>

namespace FROSch {
    
    template <class SC,class LO,class GO,class NO>
    class EntitySet {
        
    public:
        
        typedef Xpetra::Map<LO,GO,NO> Map;
        typedef Teuchos::RCP<Map> MapPtr;
        typedef Teuchos::RCP<const Map> ConstMapPtr;
        
        typedef Xpetra::Matrix<SC,LO,GO,NO> CrsMatrix;
        typedef Teuchos::RCP<CrsMatrix> CrsMatrixPtr;
        
        typedef Xpetra::MultiVector<SC,LO,GO,NO> MultiVector;
        typedef Teuchos::RCP<MultiVector> MultiVectorPtr;
        
        typedef Teuchos::RCP<EntitySet<SC,LO,GO,NO> > EntitySetPtr;
        
        typedef Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> > InterfaceEntityPtr;
        typedef Teuchos::Array<InterfaceEntityPtr> InterfaceEntityPtrVec;
        typedef Teuchos::ArrayRCP<InterfaceEntityPtr> InterfaceEntityPtrVecPtr;
        
        typedef unsigned UN;
        
        typedef Teuchos::Array<GO> GOVec;
        
        typedef Teuchos::Array<SC> SCVec;
        typedef Teuchos::ArrayRCP<SC> SCVecPtr;
        
        
        EntitySet(EntityType type);
        
        EntitySet(const EntitySet &entitySet);
        
        ~EntitySet();
        
        int addEntity(InterfaceEntityPtr entity);
        
        int buildEntityMap(ConstMapPtr localToGlobalNodesMap);
        
        int findAncestors(EntitySetPtr entitySet);
        
        int divideUnconnectedEntities(CrsMatrixPtr matrix, int pID);
        
        InterfaceEntityPtrVecPtr sortOutVertices();
        
        InterfaceEntityPtrVecPtr sortOutShortEdges();
        
        InterfaceEntityPtrVecPtr sortOutStraightEdges(UN dimension,
                                                      MultiVectorPtr &nodeList);
        
        int removeEntity(UN iD);
        
        int removeEmptyEntities();
        
        int sortUnique();
        
        bool checkForVertices();
        
        bool checkForShortEdges();
        
        bool checkForStraightEdges(UN dimension,
                                   MultiVectorPtr &nodeList);
        
        bool checkForEmptyEntities();
        
        /////////////////
        // Set Methods //
        /////////////////
        
        int setUniqueIDToFirstGlobalNodeID();
        
        int resetEntityType(EntityType type);
        
        /////////////////
        // Get Methods //
        /////////////////
        
        EntityType getEntityType() const;
        
        UN getNumEntities() const;
        
        const InterfaceEntityPtrVec & getEntityVector() const;
        
        const InterfaceEntityPtr getEntity(UN iD) const;
        
        const MapPtr getEntityMap() const;
        
        const SCVecPtr getDirection(UN dimension,
                                    MultiVectorPtr &nodeList,
                                    UN iD) const;
        
    protected:
        
        ///////////////
        // Variables //
        ///////////////
        
        EntityType Type_;
        
        InterfaceEntityPtrVec EntityVector_;
        
        bool EntityMapIsUpToDate_;
        
        MapPtr EntityMap_;
    };
    
}

#endif
