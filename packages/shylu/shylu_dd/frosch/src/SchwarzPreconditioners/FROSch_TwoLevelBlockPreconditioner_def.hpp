// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_TWOLEVELBLOCKPRECONDITIONER_DEF_HPP
#define _FROSCH_TWOLEVELBLOCKPRECONDITIONER_DEF_HPP
#include <FROSch_TwoLevelBlockPreconditioner_decl.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    TwoLevelBlockPreconditioner<SC,LO,GO,NO>::TwoLevelBlockPreconditioner(ConstXMatrixPtr k,
                                                                          ParameterListPtr parameterList) :
    OneLevelPreconditioner<SC,LO,GO,NO> (k,parameterList)
    {
        FROSCH_DETAILTIMER_START_LEVELID(twoLevelBlockPreconditionerTime,"TwoLevelBlockPreconditioner::TwoLevelBlockPreconditioner");
        if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("IPOUHarmonicCoarseOperator")) {
            // Set the LevelID in the sublist
            parameterList->sublist("IPOUHarmonicCoarseOperator").set("Level ID",this->LevelID_);
            //                FROSCH_ASSERT(false,"not implemented for block.");
            // this->ParameterList_->sublist("IPOUHarmonicCoarseOperator").sublist("InterfacePartitionOfUnity").set("Test Unconnected Interface",true);
            CoarseOperator_ = IPOUHarmonicCoarseOperatorPtr(new IPOUHarmonicCoarseOperator<SC,LO,GO,NO>(k,sublist(parameterList,"IPOUHarmonicCoarseOperator")));
        } else if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("GDSWCoarseOperator")) {
            // Set the LevelID in the sublist
            parameterList->sublist("GDSWCoarseOperator").set("Level ID",this->LevelID_);
            // this->ParameterList_->sublist("GDSWCoarseOperator").set("Test Unconnected Interface",true);
            CoarseOperator_ = GDSWCoarseOperatorPtr(new GDSWCoarseOperator<SC,LO,GO,NO>(k,sublist(parameterList,"GDSWCoarseOperator")));
        } else if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("RGDSWCoarseOperator")) {
            // Set the LevelID in the sublist
            parameterList->sublist("RGDSWCoarseOperator").set("Level ID",this->LevelID_);
            // this->ParameterList_->sublist("RGDSWCoarseOperator").set("Test Unconnected Interface",true);
            CoarseOperator_ = RGDSWCoarseOperatorPtr(new RGDSWCoarseOperator<SC,LO,GO,NO>(k,sublist(parameterList,"RGDSWCoarseOperator")));
        } else {
            FROSCH_ASSERT(false,"CoarseOperator Type unkown.");
        } // TODO: Add ability to disable individual levels
        if (this->UseMultiplicative_) {
            this->MultiplicativeOperator_->addOperator(CoarseOperator_);
        }
        else{
            this->SumOperator_->addOperator(CoarseOperator_);
        }
    }


    template <class SC,class LO,class GO,class NO>
    int TwoLevelBlockPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                             UNVecPtr dofsPerNodeVec,
                                                             DofOrderingVecPtr dofOrderingVec,
                                                             int overlap,
                                                             ConstXMapPtrVecPtr repeatedMapVec,
                                                             ConstXMultiVectorPtrVecPtr nullSpaceBasisVec,
                                                             ConstXMultiVectorPtrVecPtr nodeListVec,
                                                             ConstXMapPtrVecPtr2D dofsMapsVec,
                                                             GOVecPtr2D dirichletBoundaryDofsVec)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"TwoLevelBlockPreconditioner::initialize");

        ////////////
        // Checks //
        ////////////
        UN nmbBlocks = dofsPerNodeVec.size();
        for (UN i = 0; i < dofOrderingVec.size(); i++ ) {
            DofOrdering dofOrdering = dofOrderingVec[i];
            FROSCH_ASSERT(dofOrdering == NodeWise || dofOrdering == DimensionWise || dofOrdering == Custom,"ERROR: Specify a valid DofOrdering.");
        }
        int ret = 0;

        //////////
        // Maps //
        //////////
        this->buildRepeatedDoFMapPerBlock(repeatedMapVec,dofsMapsVec,dofsPerNodeVec);

        // Build dofsMaps and repeatedNodesMap
        ConstXMapPtrVecPtr repeatedNodesMapVec;
        if (dofsMapsVec.is_null()) {
            FROSCH_DETAILTIMER_START_LEVELID(buildDofMapsTime,"BuildDofMaps");
            if (0>BuildDofMapsVec(repeatedMapVec,dofsPerNodeVec,dofOrderingVec,repeatedNodesMapVec,dofsMapsVec)) ret -= 100; // Todo: Rückgabewerte
        } else {
            FROSCH_ASSERT(dofsMapsVec.size()==dofsPerNodeVec.size(),"dofsMapsVec.size()!=dofsPerNodeVec.size()");
            for (UN j=0; j<dofsMapsVec.size(); j++) {
                FROSCH_ASSERT(dofsMapsVec[j].size()==dofsPerNodeVec[j],"dofsMapsVec[block].size()!=dofsPerNodeVec[block]");
                for (UN i=0; i<dofsMapsVec[j].size(); i++) {
                    FROSCH_ASSERT(!dofsMapsVec[j][i].is_null(),"dofsMapsVec[block][i].is_null()");
                }
            }
            repeatedNodesMapVec = BuildNodeMapsFromDofMaps(dofsMapsVec,dofsPerNodeVec,dofOrderingVec);
        }
        ConstXMapPtr repeatedMap = MergeMaps(repeatedMapVec);

        //////////////////////////
        // Communicate nodeList //
        //////////////////////////
        if (!nodeListVec.is_null()) {
            if (this->MpiComm_->getRank() == 0) std::cout<<"Node List is NOT null\n";
            FROSCH_DETAILTIMER_START_LEVELID(communicateNodeListTime,"Communicate Node List");
            for (UN i=0; i<nodeListVec.size(); i++) {
                if (!nodeListVec[i]->getMap()->isSameAs(*repeatedNodesMapVec[i])) {
                    RCP<MultiVector<SC,LO,GO,NO> > tmpNodeList = MultiVectorFactory<SC,LO,GO,NO>::Build(repeatedNodesMapVec[i],nodeListVec[i]->getNumVectors());
                    RCP<Import<LO,GO,NO> > scatter = ImportFactory<LO,GO,NO>::Build(nodeListVec[i]->getMap(),repeatedNodesMapVec[i]);
                    tmpNodeList->doImport(*nodeListVec[i],*scatter,INSERT);
                    nodeListVec[i] = tmpNodeList.getConst();
                }
            }
        } else {
            if (this->MpiComm_->getRank() == 0) std::cout<<"Node List is null\n";
            nodeListVec.resize(nmbBlocks);
        }

        /////////////////////////////////////
        // Determine dirichletBoundaryDofs //
        /////////////////////////////////////
        this->determineDirichletBoundaryDofs(dirichletBoundaryDofsVec,
                                             repeatedMapVec,
                                             repeatedMap,
                                             dofsMapsVec);

        ////////////////////////////////////
        // Initialize OverlappingOperator //
        ////////////////////////////////////
        if (!this->ParameterList_->get("OverlappingOperator Type","AlgebraicOverlappingOperator").compare("AlgebraicOverlappingOperator")) {
            AlgebraicOverlappingOperatorPtr algebraicOverlappigOperator = rcp_static_cast<AlgebraicOverlappingOperator<SC,LO,GO,NO> >(this->OverlappingOperator_);
            if (0>algebraicOverlappigOperator->initialize(overlap,repeatedMap)) ret -= 1;
        } else {
            FROSCH_ASSERT(false,"OverlappingOperator Type unkown.");
        }

        ///////////////////////////////
        // Initialize CoarseOperator //
        ///////////////////////////////
        if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("IPOUHarmonicCoarseOperator")) {
            // Build Null Space
            if (!this->ParameterList_->get("Null Space Type","Laplace").compare("Laplace")) {
                nullSpaceBasisVec.resize(nmbBlocks);
                for (UN i = 0; i < nmbBlocks; i++) {                
                    nullSpaceBasisVec[i] = BuildNullSpace<SC,LO,GO,NO>(dimension,NullSpaceType::Laplace,repeatedMap,dofsPerNodeVec[i],dofsMapsVec[i]);
                }
            } else if (!this->ParameterList_->get("Null Space Type","Laplace").compare("Linear Elasticity")) {
                nullSpaceBasisVec.resize(repeatedMapVec.size());
                for (int i = 0; i<repeatedMapVec.size(); i++) {
                    nullSpaceBasisVec[i] = BuildNullSpace(dimension,NullSpaceType::Elasticity,repeatedMap,dofsPerNodeVec[i],dofsMapsVec[i],nodeListVec[i]);
                }
            } else if (!this->ParameterList_->get("Null Space Type","Laplace").compare("Stokes")) {
                nullSpaceBasisVec.resize(2);
                nullSpaceBasisVec[0] = BuildNullSpace<SC,LO,GO,NO>(dimension,NullSpaceType::Laplace,repeatedMap,dofsPerNodeVec[0],dofsMapsVec[0]);
                nullSpaceBasisVec[1] = BuildNullSpace<SC,LO,GO,NO>(dimension,NullSpaceType::Laplace,repeatedMap,dofsPerNodeVec[1],dofsMapsVec[1]);
            } else if (!this->ParameterList_->get("Null Space Type","Laplace").compare("Input")) {
                FROSCH_ASSERT(!nullSpaceBasisVec.is_null(),"Null Space Type is 'Input', but nullSpaceBasis.is_null().");
            } else {
                FROSCH_ASSERT(false,"Null Space Type unknown.");
            }

            IPOUHarmonicCoarseOperatorPtr iPOUHarmonicCoarseOperator = rcp_static_cast<IPOUHarmonicCoarseOperator<SC,LO,GO,NO> >(CoarseOperator_);
            if (0>iPOUHarmonicCoarseOperator->initialize(dimension,dofsPerNodeVec,repeatedNodesMapVec,dofsMapsVec,nullSpaceBasisVec,nodeListVec,dirichletBoundaryDofsVec)) ret -=10;
        } else if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("GDSWCoarseOperator")) {
            this->ParameterList_->sublist("GDSWCoarseOperator").sublist("CoarseSolver").sublist("MueLu").set("Dimension",(int)dimension);
            GDSWCoarseOperatorPtr gDSWCoarseOperator = rcp_static_cast<GDSWCoarseOperator<SC,LO,GO,NO> >(CoarseOperator_);
            if (0>gDSWCoarseOperator->initialize(dimension,dofsPerNodeVec,repeatedNodesMapVec,dofsMapsVec,dirichletBoundaryDofsVec,nodeListVec)) ret -=10;
        } else if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("RGDSWCoarseOperator")) {
            this->ParameterList_->sublist("RGDSWCoarseOperator").sublist("CoarseSolver").sublist("MueLu").set("Dimension",(int)dimension);
            RGDSWCoarseOperatorPtr rGDSWCoarseOperator = rcp_static_cast<RGDSWCoarseOperator<SC,LO,GO,NO> >(CoarseOperator_);
            if (0>rGDSWCoarseOperator->initialize(dimension,dofsPerNodeVec,repeatedNodesMapVec,dofsMapsVec,dirichletBoundaryDofsVec,nodeListVec)) ret -=10;
        } else {
            FROSCH_ASSERT(false,"CoarseOperator Type unkown.");
        }
        return ret;
    }

    template <class SC,class LO,class GO,class NO>
    int TwoLevelBlockPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                             UNVecPtr dofsPerNodeVec,
                                                             DofOrderingVecPtr dofOrderingVec,
                                                             int overlap,
                                                             ConstXMultiVectorPtrVecPtr nodeListVec,
                                                             ConstXMapPtrVecPtr repeatedMapVec,
                                                             ConstXMultiVectorPtrVecPtr nullSpaceBasisVec,
                                                             ConstXMapPtrVecPtr2D dofsMapsVec,
                                                             GOVecPtr2D dirichletBoundaryDofsVec)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"TwoLevelBlockPreconditioner::initialize");

        ////////////
        // Checks //
        ////////////
        UN nmbBlocks = dofsPerNodeVec.size();
        for (UN i = 0; i < dofOrderingVec.size(); i++ ) {
            DofOrdering dofOrdering = dofOrderingVec[i];
            FROSCH_ASSERT(dofOrdering == NodeWise || dofOrdering == DimensionWise || dofOrdering == Custom,"ERROR: Specify a valid DofOrdering.");
        }
        int ret = 0;

        //////////
        // Maps //
        //////////
        this->buildRepeatedDoFMapPerBlock(repeatedMapVec,dofsMapsVec,dofsPerNodeVec);

        // Build dofsMaps and repeatedNodesMap
        ConstXMapPtrVecPtr repeatedNodesMapVec;
        if (dofsMapsVec.is_null()) {
            FROSCH_DETAILTIMER_START_LEVELID(buildDofMapsTime,"BuildDofMaps");
            if (0>BuildDofMapsVec(repeatedMapVec,dofsPerNodeVec,dofOrderingVec,repeatedNodesMapVec,dofsMapsVec)) ret -= 100; // Todo: Rückgabewerte
        } else {
            FROSCH_ASSERT(dofsMapsVec.size()==dofsPerNodeVec.size(),"dofsMapsVec.size()!=dofsPerNodeVec.size()");
            for (UN j=0; j<dofsMapsVec.size(); j++) {
                FROSCH_ASSERT(dofsMapsVec[j].size()==dofsPerNodeVec[j],"dofsMapsVec[block].size()!=dofsPerNodeVec[block]");
                for (UN i=0; i<dofsMapsVec[j].size(); i++) {
                    FROSCH_ASSERT(!dofsMapsVec[j][i].is_null(),"dofsMapsVec[block][i].is_null()");
                }
            }
            repeatedNodesMapVec = BuildNodeMapsFromDofMaps( dofsMapsVec, dofsPerNodeVec, dofOrderingVec );
        }

        ConstXMapPtr repeatedMap;
        if (this->ParameterList_->get("Continuous Blocks",false)) {
//            repeatedMap = MergeMapsCont( repeatedMapVec );
            FROSCH_ASSERT(false,"Implement MergeMapsCont.");
        } else {
            repeatedMap = MergeMaps(repeatedMapVec);
        }

        //////////////////////////
        // Communicate nodeList //
        //////////////////////////
        if (!nodeListVec.is_null()) {
            FROSCH_DETAILTIMER_START_LEVELID(communicateNodeListTime,"Communicate Node List");
            for (UN i=0; i<nodeListVec.size(); i++) {
                ConstXMapPtr nodeListVecMap_i = nodeListVec[i]->getMap();
                if (!nodeListVecMap_i->isSameAs(*repeatedNodesMapVec[i])) {
                    RCP<MultiVector<SC,LO,GO,NO> > tmpNodeList = MultiVectorFactory<SC,LO,GO,NO>::Build(repeatedNodesMapVec[i],nodeListVec[i]->getNumVectors());
                    RCP<Import<LO,GO,NO> > scatter = ImportFactory<LO,GO,NO>::Build(nodeListVecMap_i,repeatedNodesMapVec[i]);
                    tmpNodeList->doImport(*nodeListVec[i],*scatter,INSERT);
                    nodeListVec[i] = tmpNodeList.getConst();
                }
            }
        } else {
            nodeListVec.resize(nmbBlocks);
        }

        /////////////////////////////////////
        // Determine dirichletBoundaryDofs //
        /////////////////////////////////////
        this->determineDirichletBoundaryDofs(dirichletBoundaryDofsVec,
                                             repeatedMapVec,
                                             repeatedMap,
                                             dofsMapsVec);

        ////////////////////////////////////
        // Initialize OverlappingOperator //
        ////////////////////////////////////
        if (!this->ParameterList_->get("OverlappingOperator Type","AlgebraicOverlappingOperator").compare("AlgebraicOverlappingOperator")) {
            AlgebraicOverlappingOperatorPtr algebraicOverlappigOperator = rcp_static_cast<AlgebraicOverlappingOperator<SC,LO,GO,NO> >(this->OverlappingOperator_);
            if (0>algebraicOverlappigOperator->initialize(overlap,repeatedMap)) ret -= 1;
        } else {
            FROSCH_ASSERT(false,"OverlappingOperator Type unkown.");
        }

        ///////////////////////////////
        // Initialize CoarseOperator //
        ///////////////////////////////

        if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("IPOUHarmonicCoarseOperator")) {
            // Build Null Space
            if (!this->ParameterList_->get("Null Space Type","Laplace").compare("Laplace")) {
                nullSpaceBasisVec.resize(nmbBlocks);
                for (UN i = 0; i < nmbBlocks; i++) {
                    // Todo: For now we have to use repeatedMap as the map for the vectors, such that the offsets are consistent with dofsMapsVec[i]. This would not be the case for repeatedMapVec[i]. However, it would be sufficient to use repeatedMapVec[i] with an offset instead
                    nullSpaceBasisVec[i] = BuildNullSpace<SC,LO,GO,NO>(dimension,NullSpaceType::Laplace,repeatedMap,dofsPerNodeVec[i],dofsMapsVec[i]);
                }
            } else if (!this->ParameterList_->get("Null Space Type","Laplace").compare("Linear Elasticity")) {
                nullSpaceBasisVec.resize(repeatedMapVec.size());
                for (int i = 0; i<repeatedMapVec.size(); i++) {
                    nullSpaceBasisVec[i] = BuildNullSpace(dimension,NullSpaceType::Elasticity,repeatedMap,dofsPerNodeVec[i],dofsMapsVec[i],nodeListVec[i]);
                }
            } else if (!this->ParameterList_->get("Null Space Type","Laplace").compare("Stokes")) {
                nullSpaceBasisVec.resize(2);
                nullSpaceBasisVec[0] = BuildNullSpace<SC,LO,GO,NO>(dimension,NullSpaceType::Laplace,repeatedMap,dofsPerNodeVec[0],dofsMapsVec[0]);
                nullSpaceBasisVec[1] = BuildNullSpace<SC,LO,GO,NO>(dimension,NullSpaceType::Laplace,repeatedMap,dofsPerNodeVec[1],dofsMapsVec[1]);
            } else if (!this->ParameterList_->get("Null Space Type","Laplace").compare("Input")) {
                FROSCH_ASSERT(!nullSpaceBasisVec.is_null(),"Null Space Type is 'Input', but nullSpaceBasis.is_null().");
            } else {
                FROSCH_ASSERT(false,"Null Space Type unknown.");
            }
            IPOUHarmonicCoarseOperatorPtr iPOUHarmonicCoarseOperator = rcp_static_cast<IPOUHarmonicCoarseOperator<SC,LO,GO,NO> >(CoarseOperator_);
            if (0>iPOUHarmonicCoarseOperator->initialize(dimension,dofsPerNodeVec,repeatedNodesMapVec,dofsMapsVec,nullSpaceBasisVec,nodeListVec,dirichletBoundaryDofsVec)) ret -=10;
        } else if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("GDSWCoarseOperator")) {
            this->ParameterList_->sublist("GDSWCoarseOperator").sublist("CoarseSolver").sublist("MueLu").set("Dimension",(int)dimension);
            GDSWCoarseOperatorPtr gDSWCoarseOperator = rcp_static_cast<GDSWCoarseOperator<SC,LO,GO,NO> >(CoarseOperator_);
            if (0>gDSWCoarseOperator->initialize(dimension,dofsPerNodeVec,repeatedNodesMapVec,dofsMapsVec,dirichletBoundaryDofsVec,nodeListVec)) ret -=10;
        }
        else if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("RGDSWCoarseOperator")) {
            this->ParameterList_->sublist("RGDSWCoarseOperator").sublist("CoarseSolver").sublist("MueLu").set("Dimension",(int)dimension);
            RGDSWCoarseOperatorPtr rGDSWCoarseOperator = rcp_static_cast<RGDSWCoarseOperator<SC,LO,GO,NO> >(CoarseOperator_);
            if (0>rGDSWCoarseOperator->initialize(dimension,dofsPerNodeVec,repeatedNodesMapVec,dofsMapsVec,dirichletBoundaryDofsVec,nodeListVec)) ret -=10;
        }
        else {
            FROSCH_ASSERT(false,"CoarseOperator Type unkown.");
        }
        return ret;
    }

    template <class SC,class LO,class GO,class NO>
    int TwoLevelBlockPreconditioner<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_LEVELID(computeTime,"TwoLevelBlockPreconditioner::compute");
        int ret = 0;
        if (0>this->OverlappingOperator_->compute()) ret -= 1;
        if (0>CoarseOperator_->compute()) ret -= 10;
        return ret;
    }

    template <class SC,class LO,class GO,class NO>
    void TwoLevelBlockPreconditioner<SC,LO,GO,NO>::describe(FancyOStream &out,
                                                            const EVerbosityLevel verbLevel) const
    {
        FROSCH_ASSERT(false,"describe() has to be implemented properly...");
    }

    template <class SC,class LO,class GO,class NO>
    string TwoLevelBlockPreconditioner<SC,LO,GO,NO>::description() const
    {
        return "GDSW Preconditioner";
    }

    template <class SC,class LO,class GO,class NO>
    int TwoLevelBlockPreconditioner<SC,LO,GO,NO>::resetMatrix(ConstXMatrixPtr &k)
    {
        FROSCH_DETAILTIMER_START_LEVELID(resetMatrixTime,"TwoLevelBlockPreconditioner::resetMatrix");
        this->K_ = k;
        this->OverlappingOperator_->resetMatrix(this->K_);
        CoarseOperator_->resetMatrix(this->K_);
        if (this->UseMultiplicative_) this->MultiplicativeOperator_->resetMatrix(this->K_);
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int TwoLevelBlockPreconditioner<SC,LO,GO,NO>::preApplyCoarse(XMultiVectorPtr &x,
                                                                 XMultiVectorPtr &y)
    {
        FROSCH_DETAILTIMER_START_LEVELID(preApplyCoarseTime,"TwoLevelBlockPreconditioner::preApplyCoarse");
        if (this->UseMultiplicative_) {
            this->MultiplicativeOperator_->preApplyCoarse(*x,*y);
        }
        else{
            FROSCH_ASSERT(false,"preApplyCoarse(XMultiVectorPtr &x) only implemented for MultiplicativeOperator.")
        }
        return 0;
    }

    template <class SC, class LO, class GO, class NO>
    void TwoLevelBlockPreconditioner<SC, LO, GO, NO>::buildRepeatedDoFMapPerBlock(
            ConstXMapPtrVecPtr &repeatedMapVec,
            const ConstXMapPtrVecPtr2D &dofsMapsVec,
            const UNVecPtr &dofsPerNodeVec)
    {
        // We need to know the number of blocks.
        const UN numBlocks = dofsPerNodeVec.size();
        FROSCH_ASSERT(numBlocks > 0, "No blocks specified: Set up dofsPerNodeVec to fix it.");
        for (UN i = 0; i < numBlocks; ++i) {
            FROSCH_ASSERT(dofsPerNodeVec[i] > 0, "Block " << i << " has zero degrees of freedom.");
        }

        // Make sure the block vectors' sizes match.
        if (!repeatedMapVec.is_null()) {
            FROSCH_ASSERT(repeatedMapVec.size() == numBlocks, "repeatedMapVec.size() != numBlocks");
        }
        if (!dofsMapsVec.is_null()) {
            FROSCH_ASSERT(dofsMapsVec.size() == numBlocks, "dofsMapsVec.size() != numBlocks");
        }

        // If the repeatedMapVec is not passed and also not the dofsMapsVec, we need to generate the repeated maps
        // based on the unique map of the system matrix K_. However, this is not supported for block systems (with
        // more than one block). So for true block systems, either the dofsMapsVec or the repeatedMapVec need to be
        // provided by the caller.
        if ((repeatedMapVec.is_null() || repeatedMapVec[0].is_null()) &&
            (dofsMapsVec.is_null() || dofsMapsVec[0].is_null())) {
            FROSCH_ASSERT(
                    numBlocks == 1, "If neither the repeatedMaps per block (repeatedMapVec) "
                                    "nor the dofsMaps per block (dofsMapsVec) are provided, "
                                    "the repeated maps must be generated from the unique map "
                                    "from the system matrix K. "
                                    "However, this is not supported for true block systems (with more than 1 block).");
            repeatedMapVec.resize(numBlocks);
            repeatedMapVec[0] = BuildRepeatedMap(this->K_);
            return;
        }

        // We only resize repeatedMapVec. dofsMapsVec is not taken care of in this function.
        if (repeatedMapVec.is_null()) {
            repeatedMapVec.resize(numBlocks);
        }

        // The indidual blocks of the repeated map may still be null.
        // If they are, the DoFs maps should not be null.
        // Otherwise, either the code above is wrong, or
        // the caller passed "half initialized" block maps (e.g.,
        // where the first entry is not null, but later entries are null).
        for (UN i = 0; i < numBlocks; ++i) {
            if (repeatedMapVec[i].is_null()) {
                FROSCH_ASSERT(
                        !dofsMapsVec.is_null() && !dofsMapsVec[i].is_null(),
                        "dofsMapVec[i] should not be null at thi point.");
                FROSCH_ASSERT(
                        dofsMapsVec[i].size() == dofsPerNodeVec[i],
                        "dofsMapsVec[block].size() != dofsPerNodeVec[block]");

                for (UN k = 0; k < dofsMapsVec[i].size(); k++) {
                    FROSCH_ASSERT(!dofsMapsVec[i][k].is_null(), "dofsMapsVec[block][component].is_null()");
                }
            }
        }

        // Unless the merged, repeated maps already exist (repeatedMapVec[i]),
        // merge the DoF maps to generate them.
        GO blockOffset = 0;
        for (UN i = 0; i < numBlocks; ++i) {
            if (repeatedMapVec[i].is_null()) {
                // Merge the DoFs maps (x,y,...) of length n to a single map of length d*n.
                ConstXMapPtrVecPtr2D dofsMapsOfBlock(1);
                dofsMapsOfBlock[0] = dofsMapsVec[i];
                UNVecPtr blockDofsPerNode(1, dofsPerNodeVec[i]);
                ConstXMapPtr mergedMapOfBlock = AssembleSubdomainMap(1, dofsMapsOfBlock, blockDofsPerNode);

                // The merged map uses the numbering from the global system.
                // repeatedMapVec, however, must use an independent number for each block:
                // repeatedMapVec[0] = {0, 1, 2, 3}  -->  repeatedMapVec[0] = {0, 1, 2, 3}
                // repeatedMapVec[1] = {4, 5, 6}     -->  repeatedMapVec[1] = {0, 1, 2}
                if (blockOffset == 0) {
                    repeatedMapVec[i] = mergedMapOfBlock;
                } else {
                    Array<GO> blockIDs(mergedMapOfBlock->getLocalElementList());
                    for (auto &id : blockIDs) {
                        id -= blockOffset;
                    }
                    repeatedMapVec[i] = MapFactory<LO, GO, NO>::Build(
                            mergedMapOfBlock->lib(), Teuchos::OrdinalTraits<GO>::invalid(), blockIDs(), 0,
                            mergedMapOfBlock->getComm());
                }
            }
            blockOffset += repeatedMapVec[i]->getMaxAllGlobalIndex() + 1;
        }
    }

    template <class SC, class LO, class GO, class NO>
    void TwoLevelBlockPreconditioner<SC, LO, GO, NO>::determineDirichletBoundaryDofs(
            GOVecPtr2D &dirichletBoundaryDofsVec,
            ConstXMapPtrVecPtr repeatedMapVec,
            ConstXMapPtr repeatedMap,
            ConstXMapPtrVecPtr2D dofsMapsVec)
    {

        const bool boundaryDofsProvided = !dirichletBoundaryDofsVec.is_null();
        ParameterListPtr paramList_DirichletBnd;
        bool detectBoundaryDofs = false;
        string source;
        if (dirichletBoundaryDofsVec.is_null()) {
            FROSCH_DETAILTIMER_START_LEVELID(determineDirichletRowsTime, "Determine Dirichlet Rows");

            // Extract the parameter list that stores if and how to detect Dirichlet boundary rows.
            ParameterListPtr paramList_coarseOp;
            if (this->ParameterList_->get("CoarseOperator Type", "IPOUHarmonicCoarseOperator") ==
                "IPOUHarmonicCoarseOperator") {
                paramList_coarseOp = sublist(this->ParameterList_, "IPOUHarmonicCoarseOperator");
            } else if (
                    this->ParameterList_->get("CoarseOperator Type", "IPOUHarmonicCoarseOperator") ==
                    "GDSWCoarseOperator") {
                paramList_coarseOp = sublist(this->ParameterList_, "GDSWCoarseOperator");
            } else if (
                    this->ParameterList_->get("CoarseOperator Type", "IPOUHarmonicCoarseOperator") ==
                    "RGDSWCoarseOperator") {
                paramList_coarseOp = sublist(this->ParameterList_, "RGDSWCoarseOperator");
            }
            paramList_DirichletBnd = sublist(paramList_coarseOp, "Dirichlet boundary rows in global matrix");
            detectBoundaryDofs = paramList_DirichletBnd->get("Find Dirichlet boundary rows", true);
            // Coarse operators expect one boundary list per block, even with detection disabled.
            dirichletBoundaryDofsVec.resize(repeatedMapVec.size());
            if (detectBoundaryDofs) {
                source = getDirichletBoundarySource(paramList_DirichletBnd);
                // Initialize for each block
                for (UN j = 0; j < dirichletBoundaryDofsVec.size(); ++j) {
                    dirichletBoundaryDofsVec[j] = GOVecPtr(repeatedMapVec[j]->getLocalNumElements());
                }

                // The Dirichlet boundary rows are determined by counting the nonzeros in a row of the system matrix.
                // If there is only one nonzero, it is considered as a Dirichlet boundary row.
                // If only nonzero values are stored explicitly, the graph of the matrix can be used.
                // Otherwise, the matrix has to be used to compare the values against a tolerance of zero.
                // If the tolerance is chosen as larger than zero, small values (in magnitude) can be ignored (treated
                // as zeros).
                GOVecPtr dirichletBoundaryDofs;
                if (source == "matrix") {
                    typename ScalarTraits<SC>::magnitudeType tol =
                            getDirichletBoundaryMatrixTolerance(paramList_DirichletBnd);
                    dirichletBoundaryDofs = FindOneEntryOnlyRowsGlobal(this->K_.getConst(), repeatedMap, tol);
                } else if (source == "graph") {
                    dirichletBoundaryDofs = FindOneEntryOnlyRowsGlobal(this->K_->getCrsGraph(), repeatedMap);
                } else {
                    FROSCH_ASSERT(
                            false, "Source type to determine the Dirichlet boundary is unknown." +
                                           source);
                }

                // Assign the found Dirichlet boundary degrees of freedom to their respective block.
                LOVecPtr counterSub(repeatedMapVec.size(), 0);
                for (UN i = 0; i < dirichletBoundaryDofs.size(); i++) {
                    LO subNumber = -1;
                    for (UN j = dofsMapsVec.size(); j > 0; j--) {
                        for (UN k = 0; k < dofsMapsVec[j - 1].size(); k++) {
                            if (dirichletBoundaryDofs[i] <= dofsMapsVec[j - 1][k]->getMaxAllGlobalIndex()) {
                                subNumber = j - 1;
                            }
                        }
                    }
                    dirichletBoundaryDofsVec[subNumber][counterSub[subNumber]] = dirichletBoundaryDofs[i];
                    counterSub[subNumber]++;
                }
                for (UN i = 0; i < dirichletBoundaryDofsVec.size(); i++) {
                    dirichletBoundaryDofsVec[i].resize(counterSub[i]);
                }
            }
        }

        Array<long long> boundaryCounts(repeatedMapVec.size(),0);
        if (!dirichletBoundaryDofsVec.is_null()) {
            for (UN i=0; i<dirichletBoundaryDofsVec.size(); ++i) {
                boundaryCounts[i] = static_cast<long long>(dirichletBoundaryDofsVec[i].size());
            }
        }
        PrintDirichletBoundaryStatistics(this->MpiComm_,this->LevelID_,paramList_DirichletBnd,
                                         source,boundaryDofsProvided,detectBoundaryDofs,boundaryCounts);

    }
}

#endif
