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

#ifndef _FROSCH_ALGEBRAICOVERLAPPINGOPERATOR_DEF_HPP
#define _FROSCH_ALGEBRAICOVERLAPPINGOPERATOR_DEF_HPP

#include <FROSch_AlgebraicOverlappingOperator_decl.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    AlgebraicOverlappingOperator<SC,LO,GO,NO>::AlgebraicOverlappingOperator(ConstXMatrixPtr k,
                                                                            ParameterListPtr parameterList) :
    OverlappingOperator<SC,LO,GO,NO> (k,parameterList),
    AddingLayersStrategy_ ()
    {
        FROSCH_TIMER_START_LEVELID(algebraicOverlappingOperatorTime,"AlgebraicOverlappingOperator::AlgebraicOverlappingOperator");
        if (!this->ParameterList_->get("Adding Layers Strategy","CrsGraph").compare("CrsGraph")) {
            AddingLayersStrategy_ = LayersFromGraph;
        } else if (!this->ParameterList_->get("Adding Layers Strategy","CrsGraph").compare("CrsMatrix")) {
            AddingLayersStrategy_ = LayersFromMatrix;
        } else if (!this->ParameterList_->get("Adding Layers Strategy","CrsGraph").compare("Old")) {
            AddingLayersStrategy_ = LayersOld;
        } else {
            FROSCH_ASSERT(false,"FROSch::AlgebraicOverlappingOperator : ERROR: Specify a valid strategy for adding layers.");
        }
    }

    template <class SC,class LO,class GO,class NO>
    int AlgebraicOverlappingOperator<SC,LO,GO,NO>::initialize(int overlap,
                                                              ConstXMapPtr repeatedMap)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"AlgebraicOverlappingOperator::initialize");
        if (this->Verbose_) {
            std::cout << "\n\
+------------------------------+\n\
| AlgebraicOverlappingOperator |\n\
+------------------------------+\n";
        }

        if (repeatedMap.is_null()) repeatedMap = this->K_->getRangeMap();
        this->buildOverlappingMatrices(overlap,repeatedMap);
        this->initializeOverlappingOperator();

        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0; // RETURN VALUE!!!
    }

    template <class SC,class LO,class GO,class NO>
    int AlgebraicOverlappingOperator<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_LEVELID(computeTime,"AlgebraicOverlappingOperator::compute");
        FROSCH_ASSERT(this->IsInitialized_,"ERROR: AlgebraicOverlappingOperator has to be initialized before calling compute()");
        this->computeOverlappingOperator();
        return 0; // RETURN VALUE!!!
    }

    template <class SC,class LO,class GO,class NO>
    void AlgebraicOverlappingOperator<SC,LO,GO,NO>::describe(FancyOStream &out,
                                                             const EVerbosityLevel verbLevel) const
    {
        FROSCH_ASSERT(false,"describe() has to be implemented properly...");
    }

    template <class SC,class LO,class GO,class NO>
    std::string AlgebraicOverlappingOperator<SC,LO,GO,NO>::description() const
    {
        return "Algebraic Overlapping Operator";
    }

    template <class SC,class LO,class GO,class NO>
    int AlgebraicOverlappingOperator<SC,LO,GO,NO>::buildOverlappingMatrices(int overlap,
                                                                            ConstXMapPtr repeatedMap)
    {
        FROSCH_TIMER_START_LEVELID(buildOverlappingMatricesTime,"AlgebraicOverlappingOperator::buildOverlappingMatrices");
        // ====================================================================================
        // AH 08/09/2019: This is just temporary. Implement this properly in all the classes
        Verbosity verbosity = All;
        if (!this->ParameterList_->get("Verbosity","All").compare("None")) {
            verbosity = None;
        } else if (!this->ParameterList_->get("Verbosity","All").compare("All")) {
            verbosity = All;
        } else {
            FROSCH_ASSERT(false,"FROSch::AlgebraicOverlappingOperator : ERROR: Specify a valid verbosity level.");
        }
        // ====================================================================================

        this->OverlappingMap_ = repeatedMap;
        this->OverlappingMatrix_ = this->K_;

        LO local,sum,min,max;
        SC avg;
        if (verbosity==All) {
            local = (LO) std::max((LO) this->OverlappingMap_->getNodeNumElements(),(LO) 0);
            reduceAll(*this->MpiComm_,REDUCE_SUM,local,ptr(&sum));
            avg = std::max(sum/double(this->MpiComm_->getSize()),0.0);
            reduceAll(*this->MpiComm_,REDUCE_MIN,local,ptr(&min));
            reduceAll(*this->MpiComm_,REDUCE_MAX,local,ptr(&max));

            if (this->Verbose_) {
            std::cout << "\n\
    ------------------------------------------------------------------------------\n\
     Overlapping subdomains statistics\n\
    ------------------------------------------------------------------------------\n\
      Layer " << 0 << ":        avg / min / max             ---  " << avg << " / " << min << " / " << max << "\n";
            }
        }

        // Adding Layers of Elements to the overlapping subdomains
        ConstXCrsGraphPtr overlappingGraph = this->OverlappingMatrix_->getCrsGraph();
        for (int i=0; i<overlap; i++) {
            switch (AddingLayersStrategy_) {
                case LayersFromGraph:
                    ExtendOverlapByOneLayer(overlappingGraph,this->OverlappingMap_,overlappingGraph,this->OverlappingMap_);
                    break;

                case LayersFromMatrix:
                    ExtendOverlapByOneLayer(this->OverlappingMatrix_,this->OverlappingMap_,this->OverlappingMatrix_,this->OverlappingMap_);
                    break;

                case LayersOld:
                    ExtendOverlapByOneLayer_Old(this->OverlappingMatrix_,this->OverlappingMap_,this->OverlappingMatrix_,this->OverlappingMap_);
                    break;

                default:
                    FROSCH_ASSERT(false,"FROSch::AlgebraicOverlappingOperator : ERROR: Specify a valid strategy for adding layers.");
                    break;
            }
            if (verbosity==All) {
                local = (LO) std::max((LO) this->OverlappingMap_->getNodeNumElements(),(LO) 0);
                reduceAll(*this->MpiComm_,REDUCE_SUM,local,ptr(&sum));
                avg = std::max(sum/double(this->MpiComm_->getSize()),0.0);
                reduceAll(*this->MpiComm_,REDUCE_MIN,local,ptr(&min));
                reduceAll(*this->MpiComm_,REDUCE_MAX,local,ptr(&max));

                if (this->Verbose_) {
                    std::cout << "\
      layer " << i+1 << ":        avg / min / max             ---  " << avg << " / " << min << " / " << max << "\n";
                }
            }
        }
        if (this->Verbose_ && verbosity==All) {
            std::cout << "\
    ------------------------------------------------------------------------------\n";
        }
        
        // AH 08/28/2019 TODO: It seems that ExtendOverlapByOneLayer_Old is currently the fastest method because the map is sorted. This seems to be better for the direct solver. (At least Klu)
        if (this->ParameterList_->get("Sort Overlapping Map",true)) {
            switch (AddingLayersStrategy_) {
                case LayersFromGraph:
                    this->OverlappingMap_ = SortMapByGlobalIndex(this->OverlappingMap_);
                    break;
                    
                case LayersFromMatrix:
                    this->OverlappingMap_ = SortMapByGlobalIndex(this->OverlappingMap_);
                    break;
                    
                case LayersOld:
                    if (this->Verbose_) std::cout << "FROSch::AlgebraicOverlappingOperator : The overlapping map is already sorted" << std::endl;
                    break;
                    
                default:
                    FROSCH_ASSERT(false,"FROSch::AlgebraicOverlappingOperator : ERROR: Specify a valid strategy for adding layers.");
                    break;
            }
        }
        
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int AlgebraicOverlappingOperator<SC,LO,GO,NO>::updateLocalOverlappingMatrices()
    {
        if (this->IsComputed_) { // already computed once and we want to recycle the information. That is why we reset OverlappingMatrix_ to K_, because K_ has been reset at this point
            this->OverlappingMatrix_ = this->K_;
        }
        this->OverlappingMatrix_ = ExtractLocalSubdomainMatrix(this->OverlappingMatrix_,this->OverlappingMap_);
        return 0;
    }
}

#endif
