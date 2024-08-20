// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_GDSWPRECONDITIONER_DEF_HPP
#define _FROSCH_GDSWPRECONDITIONER_DEF_HPP

#include <FROSch_GDSWPreconditioner_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    GDSWPreconditioner<SC,LO,GO,NO>::GDSWPreconditioner(ConstXMatrixPtr k,
                                                        ParameterListPtr parameterList) :
    AlgebraicOverlappingPreconditioner<SC,LO,GO,NO> (k,parameterList)
    {
        FROSCH_DETAILTIMER_START_LEVELID(gDSWPreconditionerTime,"GDSWPreconditioner::GDSWPreconditioner");
        // Set the LevelID in the sublist
        parameterList->sublist("GDSWCoarseOperator").set("Level ID",this->LevelID_);
        CoarseOperator_.reset(new GDSWCoarseOperator<SC,LO,GO,NO>(k,sublist(parameterList,"GDSWCoarseOperator")));
        this->SumOperator_->addOperator(CoarseOperator_);
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(bool useDefaultParameters)
    {
        ConstXMapPtr repeatedMap = BuildRepeatedMap(this->K_->getCrsGraph());
        return initialize(repeatedMap,useDefaultParameters);
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(ConstXMapPtr repeatedMap,
                                                    bool useDefaultParameters)
    {
        if (useDefaultParameters) {
            return initialize(3,1,repeatedMap);
        } else {
            DofOrdering dofOrdering = NodeWise;
            if (!this->ParameterList_->get("DofOrdering","NodeWise").compare("NodeWise")) {
                dofOrdering = NodeWise;
            } else if (!this->ParameterList_->get("DofOrdering","NodeWise").compare("DimensionWise")) {
                dofOrdering = DimensionWise;
            } else if (!this->ParameterList_->get("DofOrdering","NodeWise").compare("Custom")) {
                dofOrdering = Custom;
            } else {
                FROSCH_ASSERT(false,"ERROR: Specify a valid DofOrdering.");
            }

            return initialize(this->ParameterList_->get("Dimension",1),this->ParameterList_->get("DofsPerNode",1),dofOrdering,this->ParameterList_->get("Overlap",1),repeatedMap);
        }
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(GOVecPtr &dirichletBoundaryDofs,
                                                    bool useDefaultParameters)
    {
        ConstXMapPtr repeatedMap = BuildRepeatedMap(this->K_->getCrsGraph());
        return initialize(repeatedMap,dirichletBoundaryDofs,useDefaultParameters);
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(ConstXMapPtr repeatedMap,
                                                    GOVecPtr &dirichletBoundaryDofs,
                                                    bool useDefaultParameters)
    {
        if (useDefaultParameters) {
            return initialize(3,1,repeatedMap,dirichletBoundaryDofs);
        } else {
            DofOrdering dofOrdering = NodeWise;
            if (!this->ParameterList_->get("DofOrdering","NodeWise").compare("NodeWise")) {
                dofOrdering = NodeWise;
            } else if (!this->ParameterList_->get("DofOrdering","NodeWise").compare("DimensionWise")) {
                dofOrdering = DimensionWise;
            } else if (!this->ParameterList_->get("DofOrdering","NodeWise").compare("Custom")) {
                dofOrdering = Custom;
            } else {
                FROSCH_ASSERT(false,"ERROR: Specify a valid DofOrdering.");
            }

            return initialize(this->ParameterList_->get("Dimension",1),this->ParameterList_->get("DofsPerNode",1),dofOrdering,this->ParameterList_->get("Overlap",1),repeatedMap,dirichletBoundaryDofs);
        }
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                    int overlap)
    {
        ConstXMapPtr repeatedMap = BuildRepeatedMap(this->K_->getCrsGraph());
        return initialize(dimension,overlap,repeatedMap);
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                    int overlap,
                                                    ConstXMapPtr repeatedMap)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"GDSWPreconditioner::initialize");
        int ret = 0;
        if (0>this->OverlappingOperator_->initialize(overlap,repeatedMap)) ret -= 1;
        if (0>CoarseOperator_->initialize(dimension,repeatedMap)) ret -= 10;

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                    int overlap,
                                                    ConstXMapPtr repeatedMap,
                                                    GOVecPtr &dirichletBoundaryDofs)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"GDSWPreconditioner::initialize");
        int ret = 0;
        if (0>this->OverlappingOperator_->initialize(overlap,repeatedMap)) ret -= 1;
        if (0>CoarseOperator_->initialize(dimension,repeatedMap,dirichletBoundaryDofs)) ret -= 10;

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UN dofsPerNode,
                                                    DofOrdering dofOrdering,
                                                    int overlap,
                                                    ConstXMapPtr repeatedMap)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"GDSWPreconditioner::initialize");
        FROSCH_ASSERT(dofOrdering == NodeWise || dofOrdering == DimensionWise,"ERROR: Specify a valid DofOrdering.");
        int ret = 0;
        if (0>this->OverlappingOperator_->initialize(overlap,repeatedMap)) ret -= 1;
        ConstXMapPtr repeatedNodesMap;
        ConstXMapPtrVecPtr repeatedDofMaps;
        if (0>BuildDofMaps(repeatedMap,dofsPerNode,dofOrdering,repeatedNodesMap,repeatedDofMaps)) ret -= 100;
        if (0>CoarseOperator_->initialize(dimension,dofsPerNode,repeatedNodesMap,repeatedDofMaps)) ret -=10;

        return ret;
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UN dofsPerNode,
                                                    DofOrdering dofOrdering,
                                                    int overlap,
                                                    ConstXMapPtr repeatedMap,
                                                    GOVecPtr &dirichletBoundaryDofs)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"GDSWPreconditioner::initialize");
        FROSCH_ASSERT(dofOrdering == NodeWise || dofOrdering == DimensionWise,"ERROR: Specify a valid DofOrdering.");
        int ret = 0;
        if (0>this->OverlappingOperator_->initialize(overlap,repeatedMap)) ret -= 1;
        ConstXMapPtr repeatedNodesMap;
        ConstXMapPtrVecPtr repeatedDofMaps;
        if (0>BuildDofMaps(repeatedMap,dofsPerNode,dofOrdering,repeatedNodesMap,repeatedDofMaps)) ret -= 100;
        if (0>CoarseOperator_->initialize(dimension,dofsPerNode,repeatedNodesMap,repeatedDofMaps,dirichletBoundaryDofs)) ret -=10;

        return ret;
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UN dofsPerNode,
                                                    DofOrdering dofOrdering,
                                                    int overlap,
                                                    ConstXMapPtr repeatedMap,
                                                    ConstXMultiVectorPtr &nodeList)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"GDSWPreconditioner::initialize");
        FROSCH_ASSERT(dofOrdering == NodeWise || dofOrdering == DimensionWise,"ERROR: Specify a valid DofOrdering.");
        int ret = 0;
        if (0>this->OverlappingOperator_->initialize(overlap,repeatedMap)) ret -= 1;
        ConstXMapPtr repeatedNodesMap;
        ConstXMapPtrVecPtr repeatedDofMaps;
        if (0>BuildDofMaps(repeatedMap,dofsPerNode,dofOrdering,repeatedNodesMap,repeatedDofMaps)) ret -= 100;
        if (0>CoarseOperator_->initialize(dimension,dofsPerNode,repeatedNodesMap,repeatedDofMaps,nodeList)) ret -=10;

        return ret;
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UN dofsPerNode,
                                                    DofOrdering dofOrdering,
                                                    int overlap,
                                                    ConstXMapPtr repeatedMap,
                                                    GOVecPtr &dirichletBoundaryDofs,
                                                    ConstXMultiVectorPtr &nodeList)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"GDSWPreconditioner::initialize");
        FROSCH_ASSERT(dofOrdering == NodeWise || dofOrdering == DimensionWise,"ERROR: Specify a valid DofOrdering.");
        int ret = 0;
        if (0>this->OverlappingOperator_->initialize(overlap,repeatedMap)) ret -= 1;

        ConstXMapPtr repeatedNodesMap;
        ConstXMapPtrVecPtr repeatedDofMaps;
        if (0>BuildDofMaps(repeatedMap,dofsPerNode,dofOrdering,repeatedNodesMap,repeatedDofMaps)) ret -= 100;
        if (0>CoarseOperator_->initialize(dimension,dofsPerNode,repeatedNodesMap,repeatedDofMaps,dirichletBoundaryDofs,nodeList)) ret -=10;

        return ret;
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_LEVELID(computeTime,"GDSWPreconditioner::compute");
        int ret = 0;
        if (0>this->OverlappingOperator_->compute()) ret -= 1;
        if (0>CoarseOperator_->compute()) ret -= 10;
        return ret;
    }

    template <class SC,class LO,class GO,class NO>
    void GDSWPreconditioner<SC,LO,GO,NO>::describe(FancyOStream &out,
                                                   const EVerbosityLevel verbLevel) const
    {
        FROSCH_ASSERT(false,"describe() has to be implemented properly...");
    }

    template <class SC,class LO,class GO,class NO>
    string GDSWPreconditioner<SC,LO,GO,NO>::description() const
    {
        return "GDSW Preconditioner";
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::resetMatrix(ConstXMatrixPtr &k)
    {
        FROSCH_DETAILTIMER_START_LEVELID(resetMatrixTime,"GDSWPreconditioner::resetMatrix");
        this->K_ = k;
        this->OverlappingOperator_->resetMatrix(this->K_);
        CoarseOperator_->resetMatrix(this->K_);
        return 0;
    }
}

#endif
