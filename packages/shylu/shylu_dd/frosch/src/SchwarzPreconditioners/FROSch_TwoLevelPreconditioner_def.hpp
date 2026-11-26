// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_TWOLEVELPRECONDITIONER_DEF_HPP
#define _FROSCH_TWOLEVELPRECONDITIONER_DEF_HPP

#include <FROSch_TwoLevelPreconditioner_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    TwoLevelPreconditioner<SC,LO,GO,NO>::TwoLevelPreconditioner(ConstXMatrixPtr k,
                                                                ParameterListPtr parameterList) :
    OneLevelPreconditioner<SC,LO,GO,NO> (k,parameterList)
    {
        FROSCH_DETAILTIMER_START_LEVELID(twoLevelPreconditionerTime,"TwoLevelPreconditioner::TwoLevelPreconditioner::");
        if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("IPOUHarmonicCoarseOperator")) {
            // Set the LevelID in the sublist
            parameterList->sublist("IPOUHarmonicCoarseOperator").set("Level ID",this->LevelID_);
            CoarseOperator_ = IPOUHarmonicCoarseOperatorPtr(new IPOUHarmonicCoarseOperator<SC,LO,GO,NO>(k,sublist(parameterList,"IPOUHarmonicCoarseOperator")));
        } else if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("GDSWCoarseOperator")) {
            // Set the LevelID in the sublist
            parameterList->sublist("GDSWCoarseOperator").set("Level ID",this->LevelID_);
            CoarseOperator_ = GDSWCoarseOperatorPtr(new GDSWCoarseOperator<SC,LO,GO,NO>(k,sublist(parameterList,"GDSWCoarseOperator")));
        } else if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("RGDSWCoarseOperator")) {
            // Set the LevelID in the sublist
            parameterList->sublist("RGDSWCoarseOperator").set("Level ID",this->LevelID_);
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
    int TwoLevelPreconditioner<SC,LO,GO,NO>::initialize(bool useDefaultParameters)
    {
        if (useDefaultParameters) {
            return initialize(3,1,1);
        } else {
            return initialize(this->ParameterList_->get("Dimension",1),this->ParameterList_->get("DofsPerNode",1),this->ParameterList_->get("Overlap",1));
        }
    }

    template <class SC,class LO,class GO,class NO>
    int TwoLevelPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                        int overlap,
                                                        UN dofsPerNode,
                                                        DofOrdering dofOrdering)
    {
        return initialize(dimension,dofsPerNode,overlap,null,null,dofOrdering,null,null,null);
    }

    template <class SC,class LO,class GO,class NO>
    int TwoLevelPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                        int overlap,
                                                        ConstXMapPtr repeatedMap,
                                                        UN dofsPerNode,
                                                        DofOrdering dofOrdering,
                                                        ConstXMultiVectorPtr nodeList)
    {
        return initialize(dimension,dofsPerNode,overlap,null,nodeList,dofOrdering,repeatedMap,null,null);
    }

    template <class SC,class LO,class GO,class NO>
    int TwoLevelPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                        UN dofsPerNode,
                                                        int overlap,
                                                        ConstXMultiVectorPtr nullSpaceBasis,
                                                        ConstXMultiVectorPtr nodeList,
                                                        DofOrdering dofOrdering,
                                                        ConstXMapPtr repeatedMap,
                                                        ConstXMapPtrVecPtr dofsMaps,
                                                        GOVecPtr dirichletBoundaryDofs)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"TwoLevelPreconditioner::initialize");
        ////////////
        // Checks //
        ////////////
        FROSCH_ASSERT(dofOrdering == NodeWise || dofOrdering == DimensionWise || dofOrdering == Custom,"ERROR: Specify a valid DofOrdering.");
        int ret = 0;

        //////////
        // Maps //
        //////////
        if (repeatedMap.is_null()) {
            FROSCH_DETAILTIMER_START_LEVELID(buildRepeatedMapTime, "BuildRepeatedMap");
            if (dofsMaps.is_null()) {
                repeatedMap = BuildRepeatedMap(this->K_->getCrsGraph());
            } else {
                FROSCH_ASSERT(dofsMaps.size() == dofsPerNode, "dofsMaps.size()!=dofsPerNode");
                FROSCH_ASSERT(dofsMaps.size() > 0, "dofsMaps is empty.");
                for (UN i = 0; i < dofsMaps.size(); i++) {
                    FROSCH_ASSERT(!dofsMaps[i].is_null(), "dofsMaps[i].is_null()");
                }
                // Assemble one block without changing the matrix global IDs.
                ArrayRCP<ConstXMapPtrVecPtr> blockDofsMaps(1);
                blockDofsMaps[0] = dofsMaps;
                ArrayRCP<UN> blockDofsPerNode(1, dofsPerNode);
                repeatedMap = AssembleSubdomainMap(1, blockDofsMaps, blockDofsPerNode);
            }
        }

        // Build dofsMaps and repeatedNodesMap
        ConstXMapPtr repeatedNodesMap;
        if (dofsMaps.is_null()) {
            FROSCH_DETAILTIMER_START_LEVELID(buildDofMapsTime, "BuildDofMaps");
            if (0 > BuildDofMaps(repeatedMap, dofsPerNode, dofOrdering, repeatedNodesMap, dofsMaps))
                ret -= 100; // Todo: Rückgabewerte
        } else {
            FROSCH_ASSERT(dofsMaps.size() == dofsPerNode, "dofsMaps.size()!=dofsPerNode");
            for (UN i = 0; i < dofsMaps.size(); i++) {
                FROSCH_ASSERT(!dofsMaps[i].is_null(), "dofsMaps[i].is_null()");
            }
            if (repeatedNodesMap.is_null()) {
                repeatedNodesMap = dofsMaps[0];
            }
        }

        //////////////////////////
        // Communicate nodeList //
        //////////////////////////
        if (!nodeList.is_null()) {
            FROSCH_DETAILTIMER_START_LEVELID(communicateNodeListTime,"Communicate Node List");
            ConstXMapPtr nodeListMap = nodeList->getMap();
            if (!nodeListMap->isSameAs(*repeatedNodesMap)) {
                RCP<MultiVector<SC,LO,GO,NO> > tmpNodeList = MultiVectorFactory<SC,LO,GO,NO>::Build(repeatedNodesMap,nodeList->getNumVectors());
                RCP<Import<LO,GO,NO> > scatter = ImportFactory<LO,GO,NO>::Build(nodeListMap,repeatedNodesMap);
                tmpNodeList->doImport(*nodeList,*scatter,INSERT);
                nodeList = tmpNodeList.getConst();
            }
        }

        /////////////////////////////////////
        // Determine dirichletBoundaryDofs //
        /////////////////////////////////////
        const bool boundaryDofsProvided = !dirichletBoundaryDofs.is_null();
        ParameterListPtr paramList_DirichletBnd;
        bool detectBoundaryDofs = false;
        string source;
        if (not boundaryDofsProvided) {
            FROSCH_DETAILTIMER_START_LEVELID(determineDirichletRowsTime,"Determine Dirichlet Rows");

            // Extract the parameter list that contains information on how to detect the boundary dofs.
            {
                ParameterListPtr paramList_coarseOp;
                if (this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator") == "IPOUHarmonicCoarseOperator") {
                    paramList_coarseOp = sublist(this->ParameterList_, "IPOUHarmonicCoarseOperator");
                } else if (this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator") == "GDSWCoarseOperator") {
                    paramList_coarseOp = sublist(this->ParameterList_, "GDSWCoarseOperator");
                } else if (this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator") == "RGDSWCoarseOperator") {
                    paramList_coarseOp = sublist(this->ParameterList_, "RGDSWCoarseOperator");
                }
                paramList_DirichletBnd = sublist(paramList_coarseOp, "Dirichlet boundary rows in global matrix");
            }

            detectBoundaryDofs = paramList_DirichletBnd->get("Find Dirichlet boundary rows", false);
            if (detectBoundaryDofs) {
                source = getDirichletBoundarySource(paramList_DirichletBnd);
                if (source == "matrix") {
                    typename ScalarTraits<SC>::magnitudeType tol = getDirichletBoundaryMatrixTolerance(paramList_DirichletBnd);
                    dirichletBoundaryDofs = FindOneEntryOnlyRowsGlobal(this->K_.getConst(),     repeatedMap, tol);
                } else if (source == "graph") {
                    dirichletBoundaryDofs = FindOneEntryOnlyRowsGlobal(this->K_->getCrsGraph(), repeatedMap);
                } else {
                    FROSCH_ASSERT(false,"Source type to determine the Dirichlet boundary is unknown." + source);
                }
            }
        }
        Teuchos::Array<long long> boundaryCounts(1, static_cast<long long>(dirichletBoundaryDofs.size()));
        PrintDirichletBoundaryStatistics(
                this->MpiComm_, this->LevelID_, paramList_DirichletBnd, source, boundaryDofsProvided,
                detectBoundaryDofs, boundaryCounts);

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
                nullSpaceBasis = BuildNullSpace<SC,LO,GO,NO>(dimension,NullSpaceType::Laplace,repeatedMap,dofsPerNode,dofsMaps);
            } else if (!this->ParameterList_->get("Null Space Type","Laplace").compare("Linear Elasticity")) {
                nullSpaceBasis = BuildNullSpace(dimension,NullSpaceType::Elasticity,repeatedMap,dofsPerNode,dofsMaps,nodeList);
            } else if (!this->ParameterList_->get("Null Space Type","Laplace").compare("Input")) {
                FROSCH_ASSERT(!nullSpaceBasis.is_null(),"Null Space Type is 'Input', but nullSpaceBasis.is_null().");
                ConstXMapPtr nullSpaceBasisMap = nullSpaceBasis->getMap();
                if (!nullSpaceBasisMap->isSameAs(*repeatedMap)) {
                    FROSCH_DETAILTIMER_START_LEVELID(communicateNullSpaceBasis,"Communicate Null Space");
                    RCP<MultiVector<SC,LO,GO,NO> > tmpNullSpaceBasis = MultiVectorFactory<SC,LO,GO,NO>::Build(repeatedMap,nullSpaceBasis->getNumVectors());
                    RCP<Import<LO,GO,NO> > scatter = ImportFactory<LO,GO,NO>::Build(nullSpaceBasisMap,repeatedMap);
                    tmpNullSpaceBasis->doImport(*nullSpaceBasis,*scatter,INSERT);
                    nullSpaceBasis = tmpNullSpaceBasis.getConst();
                }
            } else {
                FROSCH_ASSERT(false,"Null Space Type unknown.");
            }
            IPOUHarmonicCoarseOperatorPtr iPOUHarmonicCoarseOperator = rcp_static_cast<IPOUHarmonicCoarseOperator<SC,LO,GO,NO> >(CoarseOperator_);
            if (0>iPOUHarmonicCoarseOperator->initialize(dimension,dofsPerNode,repeatedNodesMap,dofsMaps,nullSpaceBasis,nodeList,dirichletBoundaryDofs)) ret -=10;
        } else if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("GDSWCoarseOperator")) {
            GDSWCoarseOperatorPtr gDSWCoarseOperator = rcp_static_cast<GDSWCoarseOperator<SC,LO,GO,NO> >(CoarseOperator_);
            if (0>gDSWCoarseOperator->initialize(dimension,dofsPerNode,repeatedNodesMap,dofsMaps,dirichletBoundaryDofs,nodeList)) ret -=10;
        } else if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("RGDSWCoarseOperator")) {
            RGDSWCoarseOperatorPtr rGDSWCoarseOperator = rcp_static_cast<RGDSWCoarseOperator<SC,LO,GO,NO> >(CoarseOperator_);
            if (0>rGDSWCoarseOperator->initialize(dimension,dofsPerNode,repeatedNodesMap,dofsMaps,dirichletBoundaryDofs,nodeList)) ret -=10;
        } else {
            FROSCH_ASSERT(false,"CoarseOperator Type unkown.");
        }
        return ret;
    }

    template <class SC,class LO,class GO,class NO>
    int TwoLevelPreconditioner<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_LEVELID(computeTime,"TwoLevelPreconditioner::compute");
        int ret = 0;
        if (0>this->OverlappingOperator_->compute()) ret -= 1;
        if (0>CoarseOperator_->compute()) ret -= 10;
        return ret;
    }

    template <class SC,class LO,class GO,class NO>
    void TwoLevelPreconditioner<SC,LO,GO,NO>::describe(FancyOStream &out,
                                                       const EVerbosityLevel verbLevel) const
    {
        FROSCH_ASSERT(false,"describe() has to be implemented properly...");
    }

    template <class SC,class LO,class GO,class NO>
    string TwoLevelPreconditioner<SC,LO,GO,NO>::description() const
    {
        return "GDSW Preconditioner";
    }

    template <class SC,class LO,class GO,class NO>
    int TwoLevelPreconditioner<SC,LO,GO,NO>::resetMatrix(ConstXMatrixPtr &k)
    {
        FROSCH_DETAILTIMER_START_LEVELID(resetMatrixTime,"TwoLevelPreconditioner::resetMatrix");
        this->K_ = k;
        this->OverlappingOperator_->resetMatrix(this->K_);
        CoarseOperator_->resetMatrix(this->K_);
        if (this->UseMultiplicative_) this->MultiplicativeOperator_->resetMatrix(this->K_);
        return 0;
    }
}

#endif
