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

#include <FROSch_InterfaceEntity_def.hpp>

#include <FROSch_Tools_def.hpp>


namespace FROSch {
    
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,
              class LO,
              class GO,
              class NO>
    class EntitySet {

    protected:

        using XMap                          = Map<LO,GO,NO>;
        using XMapPtr                       = RCP<XMap>;
        using ConstXMapPtr                  = RCP<const XMap>;

        using XMatrix                       = Matrix<SC,LO,GO,NO>;
        using XMatrixPtr                    = RCP<XMatrix>;
        using ConstXMatrixPtr               = RCP<const XMatrix>;

        using XMultiVector                  = MultiVector<SC,LO,GO,NO>;
        using XMultiVectorPtr               = RCP<XMultiVector>;
        using ConstXMultiVectorPtr          = RCP<const XMultiVector>;

        using EntitySetPtr                  = RCP<EntitySet<SC,LO,GO,NO> >;

        using InterfaceEntityPtr            = RCP<InterfaceEntity<SC,LO,GO,NO> >;
        using InterfaceEntityPtrVec         = Array<InterfaceEntityPtr>;
        using InterfaceEntityPtrVecPtr      = ArrayRCP<InterfaceEntityPtr>;

        using UN                            = unsigned;

        using GOVec                         = Array<GO>;
        using GOVecView                     = ArrayView<GO>;

        using SCVec                         = Array<SC>;
        using SCVecPtr                      = ArrayRCP<SC>;

    public:

        EntitySet(EntityType type);

        EntitySet(const EntitySet &entitySet);

        ~EntitySet();

        int addEntity(InterfaceEntityPtr entity);

        int addEntitySet(EntitySetPtr entitySet);

        EntitySetPtr deepCopy();
        
        int buildEntityMap(ConstXMapPtr localToGlobalNodesMap);

        int findAncestorsInSet(EntitySetPtr entitySet);

        EntitySetPtr findRoots();
        
        EntitySetPtr findLeafs();
        
        int clearAncestors();

        int clearOffspring();

        int clearRoots();
        
        int clearLeafs();

        int computeDistancesToRoots(UN dimension,
                                    ConstXMultiVectorPtr &nodeList = null,
                                    DistanceFunction distanceFunction = ConstantDistanceFunction);

        int divideUnconnectedEntities(ConstXMatrixPtr matrix,
                                      int pID);

        int flagNodes();

        int flagShortEntities();

        int flagStraightEntities(UN dimension,
                                 ConstXMultiVectorPtr &nodeList);

        EntitySetPtr sortOutEntities(EntityFlag flag);

        int removeEntity(UN iD);

        int removeNodesWithDofs(GOVecView dirichletBoundaryDofs);
        
        int removeEmptyEntities();

        int sortUnique();

        bool checkForVertices();

        bool checkForShortEdges();

        bool checkForStraightEdges(UN dimension,
                                   ConstXMultiVectorPtr &nodeList);

        bool checkForEmptyEntities();

        /////////////////
        // Set Methods //
        /////////////////

        int setUniqueIDToFirstGlobalNodeID();

        int setRootID();
        
        int setLeafID();

        int resetEntityType(EntityType type);

        /////////////////
        // Get Methods //
        /////////////////

        EntityType getEntityType() const;

        UN getNumEntities() const;

        const InterfaceEntityPtrVec & getEntityVector() const;

        const InterfaceEntityPtr getEntity(UN iD) const;

        const XMapPtr getEntityMap() const;

        const SCVecPtr getDirection(UN dimension,
                                    ConstXMultiVectorPtr &nodeList,
                                    UN iD) const;

    protected:

        ///////////////
        // Variables //
        ///////////////

        EntityType Type_ = DefaultType;

        InterfaceEntityPtrVec EntityVector_ = InterfaceEntityPtrVec(0);

        bool EntityMapIsUpToDate_ = false;

        XMapPtr EntityMap_;
    };

}

#endif
