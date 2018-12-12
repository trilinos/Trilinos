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
// Questions? Contact   Alexander Heinlein (alexander.heinlein@uni-koeln.de)
//                      Christian Hochmuth (c.hochmuth@uni-koeln.de)
//
// ************************************************************************
//@HEADER

#ifndef _FROSCH_TWOLEVELBLOCKPRECONDITIONER_DEF_HPP
#define _FROSCH_TWOLEVELBLOCKPRECONDITIONER_DEF_HPP

#include <FROSch_TwoLevelBlockPreconditioner_decl.hpp>
using namespace Teuchos;
namespace FROSch {
    
    template <class SC,class LO,class GO,class NO>
    TwoLevelBlockPreconditioner<SC,LO,GO,NO>::TwoLevelBlockPreconditioner(CrsMatrixPtr k,
                                                                ParameterListPtr parameterList) :
    OneLevelPreconditioner<SC,LO,GO,NO> (k,parameterList),
    CoarseOperator_ ()
    {
        if (this->ParameterList_->get("TwoLevel",true)) {            
            if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("IPOUHarmonicCoarseOperator")) {
//                FROSCH_ASSERT(false,"not implemented for block.");
                this->ParameterList_->sublist("IPOUHarmonicCoarseOperator").sublist("InterfacePartitionOfUnity").set("Test Unconnected Interface",false);
                CoarseOperator_ = IPOUHarmonicCoarseOperatorPtr(new IPOUHarmonicCoarseOperator<SC,LO,GO,NO>(k,sublist(parameterList,"IPOUHarmonicCoarseOperator")));
            } else if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("GDSWCoarseOperator")) {
                this->ParameterList_->sublist("GDSWCoarseOperator").set("Test Unconnected Interface",false);
                CoarseOperator_ = GDSWCoarseOperatorPtr(new GDSWCoarseOperator<SC,LO,GO,NO>(k,sublist(parameterList,"GDSWCoarseOperator")));
            } else if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("RGDSWCoarseOperator")) {
                this->ParameterList_->sublist("RGDSWCoarseOperator").set("Test Unconnected Interface",false);
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
    }
    
    
    template <class SC,class LO,class GO,class NO>
    int TwoLevelBlockPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                             UNVecPtr dofsPerNodeVec,
                                                             DofOrderingVecPtr dofOrderingVec,
                                                             int overlap,
                                                             MapPtrVecPtr repeatedMapVec,
                                                             MultiVectorPtrVecPtr nullSpaceBasisVec,
                                                             MultiVectorPtrVecPtr nodeListVec,
                                                             MapPtrVecPtr2D dofsMapsVec,
                                                             GOVecPtr2D dirichletBoundaryDofsVec)
    {
        ////////////
        // Checks //
        ////////////
        UN nmbBlocks = dofsPerNodeVec.size();
        for (UN i = 0; i < dofOrderingVec.size(); i++ ) {
            DofOrdering dofOrdering = dofOrderingVec[i];
            FROSCH_ASSERT(dofOrdering == NodeWise || dofOrdering == DimensionWise || dofOrdering == Custom,"ERROR: Specify a valid DofOrdering.");
        }
        int ret = 0;
//        //////////
//        // Maps //
//        //////////
//        if (repeatedMapVec.is_null()) {
//            ConstMapPtr tmpMap =  this->K_->getRowMap();
//            MapPtrVecPtr subMapVec = BuildSubMaps(tmpMap,blockMaxGIDVec);// Todo: Achtung, die UniqueMap könnte unsinnig verteilt sein. Falls es eine repeatedMap gibt, sollte dann die uniqueMap neu gebaut werden können. In diesem Fall, sollte man das aber basierend auf der repeatedNodesMap tun
//            repeatedMapVec = BuildRepeatedSubMaps(this->K_,subMapVec);
//        
//        }
        
        // Build dofsMaps and repeatedNodesMap
        MapPtrVecPtr repeatedNodesMapVec;
        if (dofsMapsVec.is_null()) {
            if (0>BuildDofMapsVec(repeatedMapVec,dofsPerNodeVec,dofOrderingVec,repeatedNodesMapVec,dofsMapsVec)) ret -= 100; // Todo: Rückgabewerte
            } else {
            FROSCH_ASSERT(dofsMapsVec.size()==dofsPerNodeVec.size(),"dofsMapsVec.size()!=dofsPerNodeVec.size()");
            for (UN j=0; j<dofsMapsVec.size(); j++) {
                FROSCH_ASSERT(dofsMapsVec[j].size()==dofsPerNodeVec[j],"dofsMapsVec[block].size()!=dofsPerNodeVec[block]");
                for (UN i=0; i<dofsMapsVec.size(); i++) {
                    FROSCH_ASSERT(!dofsMapsVec[j][i].is_null(),"dofsMapsVec[block][i].is_null()");
                }
            }
        }
        
        //////////////////////////
        // Communicate nodeList //
        //////////////////////////
        if (!nodeListVec.is_null()) {
            for (UN i=0; i<nodeListVec.size(); i++) {
                if (!nodeListVec[i]->getMap()->isSameAs(*repeatedNodesMapVec[i])) {
                    Teuchos::RCP<Xpetra::MultiVector<SC,LO,GO,NO> > tmpNodeList = nodeListVec[i];
                    nodeListVec[i] = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(repeatedNodesMapVec[i],tmpNodeList->getNumVectors());
                    Teuchos::RCP<Xpetra::Import<LO,GO,NO> > scatter = Xpetra::ImportFactory<LO,GO,NO>::Build(tmpNodeList->getMap(),repeatedNodesMapVec[i]);
                    nodeListVec[i]->doImport(*tmpNodeList,*scatter,Xpetra::INSERT);
                }
            }
        }
        else{
            nodeListVec.resize(nmbBlocks);
        }
        //////////////////////////////////////////
        // Determine dirichletBoundaryDofs //
        //////////////////////////////////////////
        MapPtr repeatedMap = MergeMaps(repeatedMapVec);
        if (dirichletBoundaryDofsVec.is_null()) {
            dirichletBoundaryDofsVec.resize(repeatedMapVec.size());
            LOVecPtr counterSub(repeatedMapVec.size(),0);
            for (UN j=0; j<dirichletBoundaryDofsVec.size(); j++) {
                dirichletBoundaryDofsVec[j] = GOVecPtr(repeatedMapVec[j]->getNodeNumElements());
            }
            GOVecPtr dirichletBoundaryDofs = FindOneEntryOnlyRowsGlobal(this->K_,repeatedMap);
            for (UN i=0; i<dirichletBoundaryDofs.size(); i++) {
                LO subNumber = -1;
                for (UN j = dofsMapsVec.size(); j > 0 ; j--) {
                    for (UN k=0; k<dofsMapsVec[j-1].size(); k++) {
                        if ( dirichletBoundaryDofs[i] <= dofsMapsVec[j-1][k]->getMaxAllGlobalIndex() ) {
                            subNumber = j-1;
                        }
                    }
                }
                dirichletBoundaryDofsVec[subNumber][counterSub[subNumber]] = dirichletBoundaryDofs[i];
                counterSub[subNumber]++;
            }
            
            //dirichletBoundaryDofsVec = GOVecPtr2D(repeatedMapVec.size());
            for (UN i=0; i<dirichletBoundaryDofsVec.size(); i++) {
                dirichletBoundaryDofsVec[i].resize(counterSub[i]);
            }
            
        }
        ////////////////////////////////////
        // Initialize OverlappingOperator //
        ////////////////////////////////////
        if (!this->ParameterList_->get("OverlappingOperator Type","AlgebraicOverlappingOperator").compare("AlgebraicOverlappingOperator")) {
            AlgebraicOverlappingOperatorPtr algebraicOverlappigOperator = Teuchos::rcp_static_cast<AlgebraicOverlappingOperator<SC,LO,GO,NO> >(this->OverlappingOperator_);
            if (0>algebraicOverlappigOperator->initialize(overlap,repeatedMap)) ret -= 1;
        } else {
            FROSCH_ASSERT(false,"OverlappingOperator Type unkown.");
        }
        ///////////////////////////////
        // Initialize CoarseOperator //
        ///////////////////////////////
        if (this->ParameterList_->get("TwoLevel",true)) {
            if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("IPOUHarmonicCoarseOperator")) {
                this->ParameterList_->sublist("IPOUHarmonicCoarseOperator").sublist("CoarseSolver").sublist("MueLu").set("Dimension",(int)dimension);
                // Build Null Space
                if (!this->ParameterList_->get("Null Space Type","Stokes").compare("Stokes")) {
                    nullSpaceBasisVec.resize(2);
                    nullSpaceBasisVec[0] = BuildNullSpace<SC,LO,GO,NO>(dimension,LaplaceNullSpace,repeatedMapVec[0],dofsPerNodeVec[0],dofsMapsVec[0]);
                    nullSpaceBasisVec[1] = BuildNullSpace<SC,LO,GO,NO>(dimension,LaplaceNullSpace,repeatedMapVec[1],dofsPerNodeVec[1],dofsMapsVec[1]);
                } else if (!this->ParameterList_->get("Null Space Type","Stokes").compare("Input")) {
                    FROSCH_ASSERT(!nullSpaceBasisVec.is_null(),"Null Space Type is 'Input', but nullSpaceBasis.is_null().");
                } else {
                    FROSCH_ASSERT(false,"Null Space Type unknown.");
                }
                IPOUHarmonicCoarseOperatorPtr iPOUHarmonicCoarseOperator = Teuchos::rcp_static_cast<IPOUHarmonicCoarseOperator<SC,LO,GO,NO> >(CoarseOperator_);
                if (0>iPOUHarmonicCoarseOperator->initialize(dimension,dofsPerNodeVec,repeatedNodesMapVec,dofsMapsVec,nullSpaceBasisVec,nodeListVec,dirichletBoundaryDofsVec)) ret -=10;
            } else if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("GDSWCoarseOperator")) {
                this->ParameterList_->sublist("GDSWCoarseOperator").sublist("CoarseSolver").sublist("MueLu").set("Dimension",(int)dimension);
                GDSWCoarseOperatorPtr gDSWCoarseOperator = Teuchos::rcp_static_cast<GDSWCoarseOperator<SC,LO,GO,NO> >(CoarseOperator_);
                if (0>gDSWCoarseOperator->initialize(dimension,dofsPerNodeVec,repeatedNodesMapVec,dofsMapsVec,dirichletBoundaryDofsVec,nodeListVec)) ret -=10;
            }
            else if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("RGDSWCoarseOperator")) {
                this->ParameterList_->sublist("RGDSWCoarseOperator").sublist("CoarseSolver").sublist("MueLu").set("Dimension",(int)dimension);
                RGDSWCoarseOperatorPtr rGDSWCoarseOperator = Teuchos::rcp_static_cast<RGDSWCoarseOperator<SC,LO,GO,NO> >(CoarseOperator_);
                if (0>rGDSWCoarseOperator->initialize(dimension,dofsPerNodeVec,repeatedNodesMapVec,dofsMapsVec,dirichletBoundaryDofsVec,nodeListVec)) ret -=10;
            }
            else {
                FROSCH_ASSERT(false,"CoarseOperator Type unkown.");
            }
        }
        return ret;
    }
    
    template <class SC,class LO,class GO,class NO>
    int TwoLevelBlockPreconditioner<SC,LO,GO,NO>::compute()
    {
        int ret = 0;
        if (0>this->OverlappingOperator_->compute()) ret -= 1;
        if (this->ParameterList_->get("TwoLevel",true)) {
            if (0>CoarseOperator_->compute()) ret -= 10;
        }
        return ret;
    }
    
    template <class SC,class LO,class GO,class NO>
    void TwoLevelBlockPreconditioner<SC,LO,GO,NO>::describe(Teuchos::FancyOStream &out,
                                                   const Teuchos::EVerbosityLevel verbLevel) const
    {
        FROSCH_ASSERT(false,"describe() has be implemented properly...");
    }
    
    template <class SC,class LO,class GO,class NO>
    std::string TwoLevelBlockPreconditioner<SC,LO,GO,NO>::description() const
    {
        return "GDSW Preconditioner";
    }
    
    template <class SC,class LO,class GO,class NO>
    int TwoLevelBlockPreconditioner<SC,LO,GO,NO>::resetMatrix(CrsMatrixPtr &k)
    {
        this->K_ = k;
        this->OverlappingOperator_->resetMatrix(this->K_);
        if (this->ParameterList_->get("TwoLevel",true)) {
            CoarseOperator_->resetMatrix(this->K_);
            if (this->UseMultiplicative_) this->MultiplicativeOperator_->resetMatrix(this->K_);
        }
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int TwoLevelBlockPreconditioner<SC,LO,GO,NO>::preApplyCoarse(MultiVectorPtr &x,MultiVectorPtr &y)
    {
        if (this->UseMultiplicative_) {
            this->MultiplicativeOperator_->preApplyCoarse(*x,*y);
        }
        else{
            FROSCH_ASSERT(false,"preApplyCoarse(MultiVectorPtr &x) only implemented for MultiplicativeOperator.")
        }
        return 0;
    }
    
}

#endif
