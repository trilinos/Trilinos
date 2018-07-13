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

#ifndef _FROSCH_GDSWPRECONDITIONER_DEF_HPP
#define _FROSCH_GDSWPRECONDITIONER_DEF_HPP

#include <FROSch_GDSWPreconditioner_decl.hpp>

namespace FROSch {
    
    template <class SC,class LO,class GO,class NO>
    GDSWPreconditioner<SC,LO,GO,NO>::GDSWPreconditioner(CrsMatrixPtr k,
                                                        ParameterListPtr parameterList) :
    AlgebraicOverlappingPreconditioner<SC,LO,GO,NO> (k,parameterList),
    CoarseLevelOperator_ (new GDSWCoarseOperator<SC,LO,GO,NO>(k,sublist(parameterList,"GDSWOperator")))
    {
        this->SumOperator_->addOperator(CoarseLevelOperator_);
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(bool useDefaultParameters)
    {
        MapPtr repeatedMap = BuildRepeatedMap(this->K_);
        return initialize(repeatedMap,useDefaultParameters);
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(MapPtr repeatedMap,
                                                    bool useDefaultParameters)
    {
        if (useDefaultParameters) {
            return initialize(3,1,repeatedMap);
        } else {
            DofOrdering dofOrdering;
            if (!this->ParameterList_->get("DofOrdering","NodeOrdering").compare("NodeWise")) {
                dofOrdering = NodeWise;
            } else if (!this->ParameterList_->get("DofOrdering","NodeOrdering").compare("DimensionWise")) {
                dofOrdering = DimensionWise;
            } else if (!this->ParameterList_->get("DofOrdering","NodeOrdering").compare("Custom")) {
                dofOrdering = Custom;
            } else {
                FROSCH_ASSERT(0!=0,"ERROR: Specify a valid DofOrdering.");
            }
            
            return initialize(this->ParameterList_->get("Dimension",1),this->ParameterList_->get("DofsPerNode",1),dofOrdering,this->ParameterList_->get("Overlap",1),repeatedMap);
        }
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(GOVecPtr &dirichletBoundaryDofs,
                                                    bool useDefaultParameters)
    {
        MapPtr repeatedMap = BuildRepeatedMap(this->K_);
        return initialize(repeatedMap,dirichletBoundaryDofs,useDefaultParameters);
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(MapPtr repeatedMap,
                                                    GOVecPtr &dirichletBoundaryDofs,
                                                    bool useDefaultParameters)
    {
        if (useDefaultParameters) {
            return initialize(3,1,repeatedMap,dirichletBoundaryDofs);
        } else {
            DofOrdering dofOrdering;
            if (!this->ParameterList_->get("DofOrdering","NodeOrdering").compare("NodeWise")) {
                dofOrdering = NodeWise;
            } else if (!this->ParameterList_->get("DofOrdering","NodeOrdering").compare("DimensionWise")) {
                dofOrdering = DimensionWise;
            } else if (!this->ParameterList_->get("DofOrdering","NodeOrdering").compare("Custom")) {
                dofOrdering = Custom;
            } else {
                FROSCH_ASSERT(0!=0,"ERROR: Specify a valid DofOrdering.");
            }
            
            return initialize(this->ParameterList_->get("Dimension",1),this->ParameterList_->get("DofsPerNode",1),dofOrdering,this->ParameterList_->get("Overlap",1),repeatedMap,dirichletBoundaryDofs);
        }
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                    int overlap)
    {
        MapPtr repeatedMap = BuildRepeatedMap(this->K_);
        return initialize(dimension,overlap,repeatedMap);
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                    int overlap,
                                                    MapPtr repeatedMap)
    {
        int ret = 0;
        if (0>this->FirstLevelOperator_->initialize(overlap,repeatedMap)) ret -= 1;
        if (0>CoarseLevelOperator_->initialize(dimension,repeatedMap)) ret -= 10; 
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                    int overlap,
                                                    MapPtr repeatedMap,
                                                    GOVecPtr &dirichletBoundaryDofs)
    {
        int ret = 0;
        if (0>this->FirstLevelOperator_->initialize(overlap,repeatedMap)) ret -= 1;
        if (0>CoarseLevelOperator_->initialize(dimension,repeatedMap,dirichletBoundaryDofs)) ret -= 10;
        
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UN dofsPerNode,
                                                    DofOrdering dofOrdering,
                                                    int overlap,
                                                    MapPtr repeatedMap)
    {
        FROSCH_ASSERT(dofOrdering == NodeWise || dofOrdering == DimensionWise,"ERROR: Specify a valid DofOrdering.");
        int ret = 0;
        if (0>this->FirstLevelOperator_->initialize(overlap,repeatedMap)) ret -= 1;
        MapPtr repeatedNodesMap;
        MapPtrVecPtr repeatedDofMaps;
        if (0>BuildDofMaps(repeatedMap,dofsPerNode,dofOrdering,repeatedNodesMap,repeatedDofMaps)) ret -= 100;
        if (0>CoarseLevelOperator_->initialize(dimension,dofsPerNode,repeatedNodesMap,repeatedDofMaps)) ret -=10;
        
        return ret;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UN dofsPerNode,
                                                    DofOrdering dofOrdering,
                                                    int overlap,
                                                    MapPtr repeatedMap,
                                                    GOVecPtr &dirichletBoundaryDofs)
    {
        FROSCH_ASSERT(dofOrdering == NodeWise || dofOrdering == DimensionWise,"ERROR: Specify a valid DofOrdering.");
        int ret = 0;
        if (0>this->FirstLevelOperator_->initialize(overlap,repeatedMap)) ret -= 1;
        MapPtr repeatedNodesMap;
        MapPtrVecPtr repeatedDofMaps;
        if (0>BuildDofMaps(repeatedMap,dofsPerNode,dofOrdering,repeatedNodesMap,repeatedDofMaps)) ret -= 100;
        if (0>CoarseLevelOperator_->initialize(dimension,dofsPerNode,repeatedNodesMap,repeatedDofMaps,dirichletBoundaryDofs)) ret -=10;
        
        return ret;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UN dofsPerNode,
                                                    DofOrdering dofOrdering,
                                                    int overlap,
                                                    MapPtr repeatedMap,
                                                    MultiVectorPtr &nodeList)
    {
        FROSCH_ASSERT(dofOrdering == NodeWise || dofOrdering == DimensionWise,"ERROR: Specify a valid DofOrdering.");
        int ret = 0;
        if (0>this->FirstLevelOperator_->initialize(overlap,repeatedMap)) ret -= 1;
        MapPtr repeatedNodesMap;
        MapPtrVecPtr repeatedDofMaps;
        if (0>BuildDofMaps(repeatedMap,dofsPerNode,dofOrdering,repeatedNodesMap,repeatedDofMaps)) ret -= 100;
        if (0>CoarseLevelOperator_->initialize(dimension,dofsPerNode,repeatedNodesMap,repeatedDofMaps,nodeList)) ret -=10;
        
        return ret;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UN dofsPerNode,
                                                    DofOrdering dofOrdering,
                                                    int overlap,
                                                    MapPtr repeatedMap,
                                                    GOVecPtr &dirichletBoundaryDofs,
                                                    MultiVectorPtr &nodeList)
    {
        FROSCH_ASSERT(dofOrdering == NodeWise || dofOrdering == DimensionWise,"ERROR: Specify a valid DofOrdering.");
        int ret = 0;
        if (0>this->FirstLevelOperator_->initialize(overlap,repeatedMap)) ret -= 1;
        
        MapPtr repeatedNodesMap;
        MapPtrVecPtr repeatedDofMaps;
        if (0>BuildDofMaps(repeatedMap,dofsPerNode,dofOrdering,repeatedNodesMap,repeatedDofMaps)) ret -= 100;
        if (0>CoarseLevelOperator_->initialize(dimension,dofsPerNode,repeatedNodesMap,repeatedDofMaps,dirichletBoundaryDofs,nodeList)) ret -=10;
        
        return ret;
    }
    
    template <class SC,class LO,class GO,class NO>
    int GDSWPreconditioner<SC,LO,GO,NO>::compute()
    {
        int ret = 0;
        if (0>this->FirstLevelOperator_->compute()) ret -= 1;
        if (0>CoarseLevelOperator_->compute()) ret -= 10;
        return ret;
    }
    
    template <class SC,class LO,class GO,class NO>
    void GDSWPreconditioner<SC,LO,GO,NO>::describe(Teuchos::FancyOStream &out,
                                                   const Teuchos::EVerbosityLevel verbLevel) const
    {
        FROSCH_ASSERT(0!=0,"describe() has be implemented properly...");
    }
    
    template <class SC,class LO,class GO,class NO>
    std::string GDSWPreconditioner<SC,LO,GO,NO>::description() const
    {
        return "GDSW Preconditioner";
    }
    
}

#endif
