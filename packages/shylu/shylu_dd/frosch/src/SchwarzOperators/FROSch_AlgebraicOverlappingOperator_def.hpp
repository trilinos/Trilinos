// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_ALGEBRAICOVERLAPPINGOPERATOR_DEF_HPP
#define _FROSCH_ALGEBRAICOVERLAPPINGOPERATOR_DEF_HPP

#include <FROSch_AlgebraicOverlappingOperator_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    AlgebraicOverlappingOperator<SC,LO,GO,NO>::AlgebraicOverlappingOperator(ConstXMatrixPtr k,
                                                                            ParameterListPtr parameterList) :
    OverlappingOperator<SC,LO,GO,NO> (k,parameterList)
    {
        FROSCH_DETAILTIMER_START_LEVELID(algebraicOverlappingOperatorTime,"AlgebraicOverlappingOperator::AlgebraicOverlappingOperator");
        if (!this->ParameterList_->get("Adding Layers Strategy","CrsGraph").compare("CrsGraph")) {
            AddingLayersStrategy_ = LayersFromGraph;
        } else if (!this->ParameterList_->get("Adding Layers Strategy","CrsGraph").compare("CrsMatrix")) {
            AddingLayersStrategy_ = LayersFromMatrix;
        } else if (!this->ParameterList_->get("Adding Layers Strategy","CrsGraph").compare("Old")) {
            AddingLayersStrategy_ = LayersOld;
        } else {
            FROSCH_ASSERT(false,"FROSch::AlgebraicOverlappingOperator: Specify a valid strategy for adding layers.");
        }
    }

    template <class SC,class LO,class GO,class NO>
    int AlgebraicOverlappingOperator<SC,LO,GO,NO>::initialize(int overlap,
                                                              ConstXMapPtr repeatedMap)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"AlgebraicOverlappingOperator::initialize");

        if (this->Verbose_) {
            cout
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << setw(89) << "-----------------------------------------------------------------------------------------"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << "| "
            << left << setw(74) << "AlgebraicOverlappingOperator " << right << setw(8) << "(Level " << setw(2) << this->LevelID_ << ")"
            << " |"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << setw(89) << "========================================================================================="
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << "| " << left << setw(41) << "Adding layers strategy" << right
            << " | " << setw(41) << this->ParameterList_->get("Adding Layers Strategy","CrsGraph")
            << " |"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << "| " << left << setw(41) << "Combine mode in overlap" << right
            << " | " << setw(41) << this->ParameterList_->get("Combine Values in Overlap","Restricted")
            << " |"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << "| " << left << setw(41) << "Solver type" << right
            << " | " << setw(41) << this->ParameterList_->sublist("Solver").get("SolverType","Amesos2")
            << " |"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << "| " << left << setw(41) << "Solver" << right
            << " | " << setw(41) << this->ParameterList_->sublist("Solver").get("Solver","Klu")
            << " |"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << "| " << left << setw(41) << "Reuse symbolic factorization" << right
            << " | " << setw(41) << this->ParameterList_->get("Reuse: Symbolic Factorization",true)
            << " |"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << setw(89) << "-----------------------------------------------------------------------------------------"
            << endl;
        }

        if (repeatedMap.is_null()) repeatedMap = this->K_->getRangeMap();
        this->buildOverlappingMatrices(overlap,repeatedMap);
        this->initializeOverlappingOperator();
        this->updateLocalOverlappingMatrices_Symbolic();
        bool reuseSymbolicFactorization = this->ParameterList_->get("Reuse: Symbolic Factorization",true);
        if (this->ExtractLocalSubdomainMatrix_Symbolic_Done_ && reuseSymbolicFactorization) {
            // if reuseSymbolicFactorization=false, we call initializeSubdomainSolver is called during compute
            this->initializeSubdomainSolver(this->localSubdomainMatrix_);
        }

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
    string AlgebraicOverlappingOperator<SC,LO,GO,NO>::description() const
    {
        return "Algebraic Overlapping Operator";
    }

    template <class SC,class LO,class GO,class NO>
    int AlgebraicOverlappingOperator<SC,LO,GO,NO>::buildOverlappingMatrices(int overlap,
                                                                            ConstXMapPtr repeatedMap)
    {
        FROSCH_DETAILTIMER_START_LEVELID(buildOverlappingMatricesTime,"AlgebraicOverlappingOperator::buildOverlappingMatrices");
        // ====================================================================================
        // AH 08/09/2019: This is just temporary. Implement this properly in all the classes
        Verbosity verbosity = All;
        if (!this->ParameterList_->get("Verbosity","All").compare("None")) {
            verbosity = None;
        } else if (!this->ParameterList_->get("Verbosity","All").compare("All")) {
            verbosity = All;
        } else {
            FROSCH_ASSERT(false,"FROSch::AlgebraicOverlappingOperator: Specify a valid verbosity level.");
        }
        // ====================================================================================

        this->OverlappingMap_ = repeatedMap;
        this->OverlappingMatrix_ = this->K_;

        GO global = 0, sum = 0;
        LO local,minVal,maxVal;
        SC avg;
        if (verbosity==All) {
            FROSCH_DETAILTIMER_START_LEVELID(printStatisticsTime,"print statistics");

            global = this->OverlappingMap_->getMaxAllGlobalIndex();
            if (this->OverlappingMap_->lib()==UseEpetra || this->OverlappingMap_->getGlobalNumElements()>0) {
                global += 1;
            }

            local = (LO) max((LO) this->OverlappingMap_->getLocalNumElements(),(LO) 0);
            reduceAll(*this->MpiComm_,REDUCE_SUM,GO(local),ptr(&sum));
            avg = max(sum/double(this->MpiComm_->getSize()),0.0);
            reduceAll(*this->MpiComm_,REDUCE_MIN,local,ptr(&minVal));
            reduceAll(*this->MpiComm_,REDUCE_MAX,local,ptr(&maxVal));

            if (this->Verbose_) {
                cout
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << setw(89) << "-----------------------------------------------------------------------------------------"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| "
                << left << setw(74) << "> Overlapping Subdomains Statistics " << right << setw(8) << "(Level " << setw(2) << this->LevelID_ << ")"
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << setw(89) << "========================================================================================="
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| " << left << setw(20) << " " << right
                << " | " << setw(10) << "total"
                << " | " << setw(10) << "avg"
                << " | " << setw(10) << "min"
                << " | " << setw(10) << "max"
                << " | " << setw(10) << "global sum"
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << setw(89) << "-----------------------------------------------------------------------------------------"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| " << left << setw(20) << "Layer 0" << right
                << " | " << setw(10) << global
                << " | " << setw(10) << setprecision(5) << avg
                << " | " << setw(10) << minVal
                << " | " << setw(10) << maxVal
                << " | " << setw(10) << sum
                << " |";
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
                    FROSCH_ASSERT(false,"FROSch::AlgebraicOverlappingOperator: Specify a valid strategy for adding layers.");
            }
            if (verbosity==All) {
                FROSCH_DETAILTIMER_START_LEVELID(printStatisticsTime,"print statistics");
                local = (LO) max((LO) this->OverlappingMap_->getLocalNumElements(),(LO) 0);
                reduceAll(*this->MpiComm_,REDUCE_SUM,GO(local),ptr(&sum));
                avg = max(sum/double(this->MpiComm_->getSize()),0.0);
                reduceAll(*this->MpiComm_,REDUCE_MIN,local,ptr(&minVal));
                reduceAll(*this->MpiComm_,REDUCE_MAX,local,ptr(&maxVal));

                if (this->Verbose_) {
                    cout
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << "| " << left << "Layer " << setw(14) << i+1 << right
                    << " | " << setw(10) << global
                    << " | " << setw(10) << setprecision(5) << avg
                    << " | " << setw(10) << minVal
                    << " | " << setw(10) << maxVal
                    << " | " << setw(10) << sum
                    << " |";
                }
            }
        }

        if (verbosity==All) {
            FROSCH_DETAILTIMER_START_LEVELID(printStatisticsTime,"print statistics");
            if (this->Verbose_) {
                cout
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << setw(89) << "-----------------------------------------------------------------------------------------"
                << endl;
            }
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
                    if (this->Verbose_) cout << "FROSch::AlgebraicOverlappingOperator : The overlapping map is already sorted" << endl;
                    break;

                default:
                    FROSCH_ASSERT(false,"FROSch::AlgebraicOverlappingOperator: Specify a valid strategy for adding layers.");
            }
        }

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int AlgebraicOverlappingOperator<SC,LO,GO,NO>::updateLocalOverlappingMatrices()
    {
        FROSCH_DETAILTIMER_START_LEVELID(updateLocalOverlappingMatricesTime,"AlgebraicOverlappingOperator::updateLocalOverlappingMatrices");
        if (this->ExtractLocalSubdomainMatrix_Symbolic_Done_) {
            // using original K_ as input
            ExtractLocalSubdomainMatrix_Compute(this->subdomainScatter_, this->K_, this->subdomainMatrix_, this->localSubdomainMatrix_);
            this->OverlappingMatrix_ = this->localSubdomainMatrix_.getConst();
        } else {
            if (this->IsComputed_) {
                // already computed once and we want to recycle the information. That is why we reset OverlappingMatrix_ to K_, because K_ has been reset at this point
                this->OverlappingMatrix_ = this->K_;
            }
            this->OverlappingMatrix_ = ExtractLocalSubdomainMatrix(this->OverlappingMatrix_, this->OverlappingMap_);
        }
        return 0;
    }


    template <class SC,class LO,class GO,class NO>
    int AlgebraicOverlappingOperator<SC,LO,GO,NO>::updateLocalOverlappingMatrices_Symbolic()
    {
        FROSCH_DETAILTIMER_START_LEVELID(updateLocalOverlappingMatrices_SymbolicTime, "AlgebraicOverlappingOperator::updateLocalOverlappingMatrices_Symbolic");
        this->extractLocalSubdomainMatrix_Symbolic();
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    void AlgebraicOverlappingOperator<SC,LO,GO,NO>::extractLocalSubdomainMatrix_Symbolic()
    {
        if (this->OverlappingMap_->lib() == UseTpetra) {
            FROSCH_DETAILTIMER_START_LEVELID(AlgebraicOverlappin_extractLocalSubdomainMatrix_SymbolicTime,"AlgebraicOverlappinOperator::extractLocalSubdomainMatrix_Symbolic");
            // buid sudomain matrix
            this->subdomainMatrix_ = MatrixFactory<SC,LO,GO,NO>::Build(this->OverlappingMap_, this->OverlappingMatrix_->getGlobalMaxNumRowEntries());
            RCP<Import<LO,GO,NO> > scatter = ImportFactory<LO,GO,NO>::Build(this->OverlappingMatrix_->getRowMap(), this->OverlappingMap_);
            this->subdomainMatrix_->doImport(*(this->OverlappingMatrix_), *scatter, ADD);

            // Used to Map original K_ to overlapping suubdomainMatrix
            this->subdomainScatter_ = ImportFactory<LO,GO,NO>::Build(this->K_->getRowMap(), this->OverlappingMap_);

            // build local subdomain matrix
            RCP<const Comm<LO> > SerialComm = rcp(new MpiComm<LO>(MPI_COMM_SELF));
            RCP<Map<LO,GO,NO> > localSubdomainMap = MapFactory<LO,GO,NO>::Build(this->OverlappingMap_->lib(), this->OverlappingMap_->getLocalNumElements(), 0, SerialComm);
            this->localSubdomainMatrix_ = MatrixFactory<SC,LO,GO,NO>::Build(localSubdomainMap, localSubdomainMap, this->OverlappingMatrix_->getGlobalMaxNumRowEntries());

            // fill in column indexes
            ExtractLocalSubdomainMatrix_Symbolic(this->subdomainMatrix_, // input
                                                 this->localSubdomainMatrix_);   // output

            // turn flag on
            this->ExtractLocalSubdomainMatrix_Symbolic_Done_ = true;
        }
    }
}

#endif
