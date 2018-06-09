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
    
    template <class SC,class LO,class GO,class NO>
    TwoLevelPreconditioner<SC,LO,GO,NO>::TwoLevelPreconditioner(CrsMatrixPtr k,
                                                                ParameterListPtr parameterList) :
    OneLevelPreconditioner<SC,LO,GO,NO> (k,parameterList),
    CoarseOperator_ ()
    {
        if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("IPOUHarmonicCoarseOperator")) {
            CoarseOperator_ = IPOUHarmonicCoarseOperatorPtr(new IPOUHarmonicCoarseOperator<SC,LO,GO,NO>(k,sublist(parameterList,"IPOUHarmonicCoarseOperator")));
        } else if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("GDSWCoarseOperator")) {
            CoarseOperator_ = GDSWCoarseOperatorPtr(new GDSWCoarseOperator<SC,LO,GO,NO>(k,sublist(parameterList,"GDSWCoarseOperator")));
        } else if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("RGDSWCoarseOperator")) {
            CoarseOperator_ = RGDSWCoarseOperatorPtr(new RGDSWCoarseOperator<SC,LO,GO,NO>(k,sublist(parameterList,"RGDSWCoarseOperator")));
        } else {
            FROSCH_ASSERT(0!=0,"CoarseOperator Type unkown.");
        } // Todo: Möglichkeit die einzelnen Level auszuschalten
        this->SumOperator_->addOperator(CoarseOperator_);
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
        return initialize(dimension,dofsPerNode,overlap,Teuchos::null,Teuchos::null,dofOrdering,Teuchos::null,Teuchos::null,Teuchos::null);
    }
    
    template <class SC,class LO,class GO,class NO>
    int TwoLevelPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                        int overlap,
                                                        MapPtr repeatedMap,
                                                        UN dofsPerNode,
                                                        DofOrdering dofOrdering,
                                                        MultiVectorPtr nodeList)
    {
        return initialize(dimension,dofsPerNode,overlap,Teuchos::null,nodeList,dofOrdering,repeatedMap,Teuchos::null,Teuchos::null);
    }
    
    template <class SC,class LO,class GO,class NO>
    int TwoLevelPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                        UN dofsPerNode,
                                                        int overlap,
                                                        MultiVectorPtr nullSpaceBasis,
                                                        MultiVectorPtr nodeList,
                                                        DofOrdering dofOrdering,
                                                        MapPtr repeatedMap,
                                                        MapPtrVecPtr dofsMaps,
                                                        GOVecPtr dirichletBoundaryDofs)
    {
        ////////////
        // Checks //
        ////////////
        if (this->Verbose_) std::cout << "HIER CHECK EINFUEGEN!!!\n";
        
        FROSCH_ASSERT(dofOrdering == NodeWise || dofOrdering == DimensionWise || dofOrdering == Custom,"ERROR: Specify a valid DofOrdering.");
        int ret = 0;

        //////////
        // Maps //
        //////////
        if (repeatedMap.is_null()) {
            repeatedMap = BuildRepeatedMap(this->K_); // Todo: Achtung, die UniqueMap könnte unsinnig verteilt sein. Falls es eine repeatedMap gibt, sollte dann die uniqueMap neu gebaut werden können. In diesem Fall, sollte man das aber basierend auf der repeatedNodesMap tun
        }
        // Build dofsMaps and repeatedNodesMap
        MapPtr repeatedNodesMap;
        if (dofsMaps.is_null()) {
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
            if (!nodeList->getMap()->isSameAs(*repeatedNodesMap)) {
                Teuchos::RCP<Xpetra::MultiVector<SC,LO,GO,NO> > tmpNodeList = nodeList;
                nodeList = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(repeatedNodesMap,tmpNodeList->getNumVectors());
                Teuchos::RCP<Xpetra::Import<LO,GO,NO> > scatter = Xpetra::ImportFactory<LO,GO,NO>::Build(tmpNodeList->getMap(),repeatedNodesMap);
                nodeList->doImport(*tmpNodeList,*scatter,Xpetra::INSERT);
            }
        }
        
        //////////////////////////////////////////
        // Determine dirichletBoundaryDofs //
        //////////////////////////////////////////
        if (dirichletBoundaryDofs.is_null()) {
            dirichletBoundaryDofs = FindOneEntryOnlyRowsGlobal(this->K_,repeatedMap);            
        }
        
        ////////////////////////////////////
        // Initialize OverlappingOperator //
        ////////////////////////////////////
        if (!this->ParameterList_->get("OverlappingOperator Type","AlgebraicOverlappingOperator").compare("AlgebraicOverlappingOperator")) {
            AlgebraicOverlappingOperatorPtr algebraicOverlappigOperator = Teuchos::rcp_static_cast<AlgebraicOverlappingOperator<SC,LO,GO,NO> >(this->OverlappingOperator_);
            if (0>algebraicOverlappigOperator->initialize(overlap,repeatedMap)) ret -= 1;
        } else {
            FROSCH_ASSERT(0!=0,"OverlappingOperator Type unkown.");
        }
        
        ////////////////////////////////////
        // Initialize OverlappingOperator //
        ////////////////////////////////////
        if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("IPOUHarmonicCoarseOperator")) {
            // Build Null Space
            if (!this->ParameterList_->get("Null Space Type","Laplace").compare("Laplace")) {
                nullSpaceBasis = BuildNullSpace<SC,LO,GO,NO>(dimension,LaplaceNullSpace,repeatedMap,dofsPerNode,dofsMaps);
            } else if (!this->ParameterList_->get("Null Space Type","Laplace").compare("Linear Elasticity")) {
                nullSpaceBasis = BuildNullSpace(dimension,LinearElasticityNullSpace,repeatedMap,dofsPerNode,dofsMaps,nodeList);
            } else if (!this->ParameterList_->get("Null Space Type","Laplace").compare("Input")) {
                FROSCH_ASSERT(!nullSpaceBasis.is_null(),"Null Space Type is 'Input', but nullSpaceBasis.is_null().");
            } else {
                FROSCH_ASSERT(0!=0,"Null Space Type unknown.");
            }
            IPOUHarmonicCoarseOperatorPtr iPOUHarmonicCoarseOperator = Teuchos::rcp_static_cast<IPOUHarmonicCoarseOperator<SC,LO,GO,NO> >(CoarseOperator_);
            if (0>iPOUHarmonicCoarseOperator->initialize(dimension,dofsPerNode,repeatedNodesMap,dofsMaps,nullSpaceBasis,nodeList,dirichletBoundaryDofs)) ret -=10;
        } else if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("GDSWCoarseOperator")) {
            GDSWCoarseOperatorPtr gDSWCoarseOperator = Teuchos::rcp_static_cast<GDSWCoarseOperator<SC,LO,GO,NO> >(CoarseOperator_);
            if (0>gDSWCoarseOperator->initialize(dimension,dofsPerNode,repeatedNodesMap,dofsMaps,dirichletBoundaryDofs,nodeList)) ret -=10; 
        } else if (!this->ParameterList_->get("CoarseOperator Type","IPOUHarmonicCoarseOperator").compare("RGDSWCoarseOperator")) {
            RGDSWCoarseOperatorPtr rGDSWCoarseOperator = Teuchos::rcp_static_cast<RGDSWCoarseOperator<SC,LO,GO,NO> >(CoarseOperator_);
            if (0>rGDSWCoarseOperator->initialize(dimension,dofsPerNode,repeatedNodesMap,dofsMaps,dirichletBoundaryDofs,nodeList)) ret -=10;
        } else {
            FROSCH_ASSERT(0!=0,"CoarseOperator Type unkown.");
        }

        return ret;
    }
    
    template <class SC,class LO,class GO,class NO>
    int TwoLevelPreconditioner<SC,LO,GO,NO>::compute()
    {
        int ret = 0;
        if (0>this->OverlappingOperator_->compute()) ret -= 1;
        if (0>CoarseOperator_->compute()) ret -= 10;
        return ret;
    }
    
    template <class SC,class LO,class GO,class NO>
    void TwoLevelPreconditioner<SC,LO,GO,NO>::describe(Teuchos::FancyOStream &out,
                                                   const Teuchos::EVerbosityLevel verbLevel) const
    {
        FROSCH_ASSERT(0!=0,"describe() has be implemented properly...");
    }
    
    template <class SC,class LO,class GO,class NO>
    std::string TwoLevelPreconditioner<SC,LO,GO,NO>::description() const
    {
        return "GDSW Preconditioner";
    }
    
}

#endif
