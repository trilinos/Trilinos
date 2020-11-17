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
            FROSCH_DETAILTIMER_START_LEVELID(buildRepeatedMapTime,"BuildRepeatedMap");
            repeatedMap = BuildRepeatedMap(this->K_->getCrsGraph()); // Todo: Achtung, die UniqueMap könnte unsinnig verteilt sein. Falls es eine repeatedMap gibt, sollte dann die uniqueMap neu gebaut werden können. In diesem Fall, sollte man das aber basierend auf der repeatedNodesMap tun
        }
        // Build dofsMaps and repeatedNodesMap
        ConstXMapPtr repeatedNodesMap;
        if (dofsMaps.is_null()) {
            FROSCH_DETAILTIMER_START_LEVELID(buildDofMapsTime,"BuildDofMaps");
            if (0>BuildDofMaps(repeatedMap,dofsPerNode,dofOrdering,repeatedNodesMap,dofsMaps)) ret -= 100; // Todo: Rückgabewerte
        } else {
            FROSCH_ASSERT(dofsMaps.size()==dofsPerNode,"dofsMaps.size()!=dofsPerNode");
            for (UN i=0; i<dofsMaps.size(); i++) {
                FROSCH_ASSERT(!dofsMaps[i].is_null(),"dofsMaps[i].is_null()");
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
        if (dirichletBoundaryDofs.is_null()) {
            FROSCH_DETAILTIMER_START_LEVELID(determineDirichletRowsTime,"Determine Dirichlet Rows");
#ifdef FindOneEntryOnlyRowsGlobal_Matrix
            GOVecPtr dirichletBoundaryDofs = FindOneEntryOnlyRowsGlobal(this->K_.getConst(),repeatedMap);
#else
            GOVecPtr dirichletBoundaryDofs = FindOneEntryOnlyRowsGlobal(this->K_->getCrsGraph(),repeatedMap);
#endif
        }

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
                nullSpaceBasis = BuildNullSpace<SC,LO,GO,NO>(dimension,LaplaceNullSpace,repeatedMap,dofsPerNode,dofsMaps);
            } else if (!this->ParameterList_->get("Null Space Type","Laplace").compare("Linear Elasticity")) {
                nullSpaceBasis = BuildNullSpace(dimension,LinearElasticityNullSpace,repeatedMap,dofsPerNode,dofsMaps,nodeList);
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
