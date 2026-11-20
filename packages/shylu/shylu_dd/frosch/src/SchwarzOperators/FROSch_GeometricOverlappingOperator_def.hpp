// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_GEOMETRICOVERLAPPINGOPERATOR_DEF_HPP
#define _FROSCH_GEOMETRICOVERLAPPINGOPERATOR_DEF_HPP

#include <FROSch_GeometricOverlappingOperator_decl.hpp>
#include <ostream>

namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;



    template <class SC,class LO,class GO,class NO>
    GeometricOverlappingOperator<SC,LO,GO,NO>::GeometricOverlappingOperator(
        ConstXMatrixPtr  k,
        int              overlap,
        GraphPtr         dualGraph,
        ParameterListPtr parameterList) 
    : OverlappingOperator<SC,LO,GO,NO> (k,parameterList)
    , DualGraph_(dualGraph)
    {
        FROSCH_DETAILTIMER_START_LEVELID(
            optimizedSchwarzOperatorTime, 
            "GeometricOverlappingOperator::GeometricOverlappingOperator");

        this->buildOverlappingMap(overlap);
    }



    template <class SC,class LO,class GO,class NO>
    int
    GeometricOverlappingOperator<SC,LO,GO,NO>::initialize()
    {
        FROSCH_TIMER_START_LEVELID(
            initializeTime, 
            "GeometricOverlappingOperator::initialize");

        FROSCH_ASSERT(false, "GeometricOverlappingOperator cannot be built without input parameters.");
    }



    template <class SC,class LO,class GO,class NO>
    int
    GeometricOverlappingOperator<SC,LO,GO,NO>::initialize(
        XMapPtr overlappingMap
    )
    {
        FROSCH_TIMER_START_LEVELID(
            initializeTime, 
            "GeometricOverlappingOperator::initialize");

        if (this->Verbose_) {
            cout
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << setw(89) << "-----------------------------------------------------------------------------------------"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << "| "
            << left << setw(74) << "GeometricOverlappingOperator " << right << setw(8) << "(Level " << setw(2) << this->LevelID_ << ")"
            << " |"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << setw(89) << "========================================================================================="
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

        this->OverlappingMatrix_ = this->K_;
        this->OverlappingMap_    = overlappingMap;
        OverlappingGraph_        = this->OverlappingMatrix_->getCrsGraph();
        
        this->initializeOverlappingOperator();

        this->updateLocalOverlappingMatrices_Symbolic();

        this->IsInitialized_ = true;
        this->IsComputed_    = false;

        return 0; // RETURN VALUE!!!
    }



    template <class SC,class LO,class GO,class NO>
    int 
    GeometricOverlappingOperator<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_LEVELID(
            computeTime,
            "GeometricOverlappingOperator::compute");

        FROSCH_ASSERT(this->IsInitialized_,
                      "ERROR: GeometricOverlappingOperator has to be initialized before calling compute()");

        this->computeOverlappingOperator();

        return 0; // RETURN VALUE!!!
    }



    template <class SC,class LO,class GO,class NO>
    int 
    GeometricOverlappingOperator<SC,LO,GO,NO>::compute(ConstXMatrixPtr neumannMatrix, ConstXMatrixPtr robinMatrix)
    {
        FROSCH_TIMER_START_LEVELID(
            computeTime,
            "GeometricOverlappingOperator::compute");

        FROSCH_ASSERT(this->IsInitialized_,
                      "ERROR: GeometricOverlappingOperator has to be initialized before calling compute()");

        NeumannMatrix_ = neumannMatrix;
        RobinMatrix_   = robinMatrix;
        
        this->computeOverlappingOperator();

        return 0; // RETURN VALUE!!!
    }



    template <class SC, class LO, class GO, class NO>
    int
    GeometricOverlappingOperator<SC, LO, GO, NO>::communicateOverlappingTriangulation(
      XMultiVectorPtr                     nodeList,
      XMultiVectorTemplatePtr<long long>  elementList,
      XMultiVectorTemplatePtr<long long>  auxillaryList,
      XMultiVectorPtr                    &nodeListOverlapping,
      XMultiVectorTemplatePtr<long long> &elementListOverlapping,
      XMultiVectorTemplatePtr<long long> &auxillaryListOverlapping)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"GeometricOverlappingOperator::communicateTriangulation");

        // get information about the system
        long long nodesPerCell = elementList->getNumVectors();

        // Create the output lists:
        elementListOverlapping =
          Xpetra::MultiVectorFactory<long long,LO,GO,NO>::Build(this->OverlappingElementMap_, elementList->getNumVectors());
        auxillaryListOverlapping = 
          Xpetra::MultiVectorFactory<long long,LO,GO,NO>::Build(this->OverlappingElementMap_, auxillaryList->getNumVectors());

        // Create the communication plan:
        Teuchos::RCP<Xpetra::Import<LO,GO,NO>> elementImporter 
          = Xpetra::ImportFactory<LO,GO,NO>::Build(DualGraph_->getRowMap(), this->OverlappingElementMap_);

        // communicate:
        elementListOverlapping->doImport(*elementList, *elementImporter, Xpetra::INSERT);
        auxillaryListOverlapping->doImport(*auxillaryList, *elementImporter, Xpetra::INSERT);

        long long numLocalElemtents = elementListOverlapping->getLocalLength();

        // extract the new local to global vertex map from elementListOverlapping
        Teuchos::Array<long long> array(numLocalElemtents * nodesPerCell);
        for (unsigned int i = 0; i < nodesPerCell; ++i ) {
            auto data = elementListOverlapping->getData(i);
            for (unsigned int j = 0; j < numLocalElemtents; ++j) {
              array[(j * nodesPerCell) + i] = data[j];
            }
        }

        // Since, we just added all global vertex indices from the elementList,
        // there are now a lot of duplicates. Let's remove those.
        FROSch::sortunique(array);

        // Create the output node list
        XMapPtr overlappingNodeMap = 
          Xpetra::MapFactory<LO,GO,NO>::Build(DualGraph_->getMap()->lib(), DualGraph_->getMap()->getMaxGlobalIndex()+1, array(), 0, DualGraph_->getMap()->getComm());

        nodeListOverlapping = 
          Xpetra::MultiVectorFactory<double,LO,GO,NO>::Build(overlappingNodeMap,nodeList->getNumVectors());

        // Create the communication plan for the node list
        RCP<Xpetra::Import<LO,GO,NO> > vertex_importer = 
          Xpetra::ImportFactory<LO,GO,NO>::Build(nodeList->getMap(), overlappingNodeMap);

        // communicate the node list:
        nodeListOverlapping->doImport(*nodeList,*vertex_importer,Xpetra::INSERT);

        // use the overlappingNodeMap to replace the globalDofIndices with localDofIndices
        // in the elementListOverlapping
        for (unsigned int i = 0; i < nodesPerCell; ++i ) {
            auto data = elementListOverlapping->getDataNonConst(i);
            for (unsigned int j = 0; j < numLocalElemtents; ++j) {
                data[j] = overlappingNodeMap->getLocalElement(data[j]);
            }
        }

        return 0;
    }



    template <class SC,class LO,class GO,class NO>
    void 
    GeometricOverlappingOperator<SC,LO,GO,NO>::describe(
        FancyOStream          &/*out*/,
        const EVerbosityLevel  /*verbLevel*/) const
    {
        FROSCH_ASSERT(false,"GeometricOverlappingOperator::describe() is not available.");
    }



    template <class SC,class LO,class GO,class NO>
    string 
    GeometricOverlappingOperator<SC,LO,GO,NO>::description() const
    {
        return "Optimized Schwarz Operator";
    }



    template <class SC,class LO,class GO,class NO>
    int 
    GeometricOverlappingOperator<SC,LO,GO,NO>::buildOverlappingMap(int overlap)
    {
        FROSCH_DETAILTIMER_START_LEVELID(
            buildOverlappingMatricesTime,
            "GeometricOverlappingOperator::buildOverlappingMatrices");

        // ====================================================================================
        // AH 08/09/2019: This is just temporary. Implement this properly in all the classes
        Verbosity verbosity = All;
        if (!this->ParameterList_->get("Verbosity","All").compare("None")) {
            verbosity = None;
        } else if (!this->ParameterList_->get("Verbosity","All").compare("All")) {
            verbosity = All;
        } else {
            FROSCH_ASSERT(false,"FROSch::GeometricOverlappingOperator: Specify a valid verbosity level.");
        }
        // ====================================================================================

        OverlappingGraph_      = DualGraph_;
        OverlappingElementMap_ = DualGraph_->getColMap();

        GO global = 0, sum = 0;
        LO local,minVal,maxVal;
        SC avg;
        if (verbosity==All) {
            FROSCH_DETAILTIMER_START_LEVELID(printStatisticsTime,"print statistics");

            global = this->OverlappingElementMap_->getMaxAllGlobalIndex();
            if (this->OverlappingElementMap_->lib()==UseEpetra || this->OverlappingElementMap_->getGlobalNumElements()>0) {
                global += 1;
            }

            local = (LO) max((LO) this->OverlappingElementMap_->getLocalNumElements(),(LO) 0);
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
        for (int i = 0; i < overlap; i++) {
            ExtendOverlapByOneLayer(
                OverlappingGraph_,
                OverlappingElementMap_,
                OverlappingGraph_,
                OverlappingElementMap_);

            if (verbosity==All) {
                FROSCH_DETAILTIMER_START_LEVELID(printStatisticsTime,"print statistics");
                local = (LO) max((LO) this->OverlappingElementMap_->getLocalNumElements(),(LO) 0);
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

        // TODO: MOVE!
        // AH 08/28/2019 TODO: It seems that ExtendOverlapByOneLayer_Old 
        // is currently the fastest method because the map is sorted. 
        // This seems to be better for the direct solver. (At least Klu)
        //if (this->ParameterList_->get("Sort Overlapping Map", true)) {
        //    this->OverlappingMap_ = SortMapByGlobalIndex(this->OverlappingMap_);
        //}

        return 0;
    }



    template <class SC,class LO,class GO,class NO>
    int 
    GeometricOverlappingOperator<SC,LO,GO,NO>::updateLocalOverlappingMatrices()
    {
        FROSCH_DETAILTIMER_START_LEVELID(
            updateLocalOverlappingMatricesTime,
            "GeometricOverlappingOperator::updateLocalOverlappingMatrices");

        if ( this->ExtractLocalSubdomainMatrix_Symbolic_Done_ ) {
            // using original K_ as input
            ExtractLocalSubdomainMatrix_Compute(this->K_, this->subdomainMatrix_, this->localSubdomainMatrix_);
            this->OverlappingMatrix_ = this->localSubdomainMatrix_.getConst();
        } else {
            if ( this->IsComputed_ ) {
                // already computed once and we want to recycle the information. 
                // That is why we reset OverlappingMatrix_ to K_, because K_ has 
                // been reset at this point
                this->OverlappingMatrix_ = this->K_;
            }
            this->OverlappingMatrix_ = ExtractLocalSubdomainMatrix(this->OverlappingMatrix_, this->OverlappingMap_);
        }

        auto out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));
        

        Teuchos::RCP<XMatrix> overlappingMatrix = null;
        MatrixMatrix<SC,LO,GO,NO>::TwoMatrixAdd(*NeumannMatrix_, false, 1.0, *RobinMatrix_, false, 1.0, overlappingMatrix, *out);
        overlappingMatrix->fillComplete();
        this->OverlappingMatrix_ = overlappingMatrix;

        return 0;
    }



    template <class SC,class LO,class GO,class NO>
    int 
    GeometricOverlappingOperator<SC,LO,GO,NO>::updateLocalOverlappingMatrices_Symbolic()
    {
        FROSCH_DETAILTIMER_START_LEVELID(
            updateLocalOverlappingMatrices_SymbolicTime, 
            "GeometricOverlappingOperator::updateLocalOverlappingMatrices_Symbolic");

        this->extractLocalSubdomainMatrix_Symbolic();

        return 0;
    }



    template <class SC,class LO,class GO,class NO>
    void 
    GeometricOverlappingOperator<SC,LO,GO,NO>::extractLocalSubdomainMatrix_Symbolic()
    {
        if ( this->OverlappingMap_->lib() == UseTpetra ) {
            FROSCH_DETAILTIMER_START_LEVELID(
                extractLocalSubdomainMatrix_SymbolicTime,
                "GeometricOverlappingOperator::extractLocalSubdomainMatrix_Symbolic");

            // build sudomain matrix
            this->subdomainMatrix_ = 
              MatrixFactory<SC,LO,GO,NO>::Build(this->OverlappingMap_, this->OverlappingMatrix_->getGlobalMaxNumRowEntries());

            RCP<Import<LO,GO,NO> > scatter = 
              ImportFactory<LO,GO,NO>::Build(this->OverlappingMatrix_->getRowMap(), this->OverlappingMap_);

            this->subdomainMatrix_->doImport(*(this->OverlappingMatrix_), *scatter, ADD);

            // build local subdomain matrix
            RCP<const Comm<LO> > SerialComm = rcp(new MpiComm<LO>(MPI_COMM_SELF));
            RCP<Map<LO,GO,NO> > localSubdomainMap = 
              MapFactory<LO,GO,NO>::Build(
                  this->OverlappingMap_->lib(), 
                  this->OverlappingMap_->getLocalNumElements(), 0, SerialComm);

            this->localSubdomainMatrix_ = 
              MatrixFactory<SC,LO,GO,NO>::Build(
                  localSubdomainMap, localSubdomainMap, 
                  this->OverlappingMatrix_->getGlobalMaxNumRowEntries());

            // fill in column indexes
            ExtractLocalSubdomainMatrix_Symbolic(this->subdomainMatrix_, // input
                                                 this->localSubdomainMatrix_);   // output

            // turn flag on
            this->ExtractLocalSubdomainMatrix_Symbolic_Done_ = true;
        }
    }



} //namespace FROSch

#endif // _FROSCH_GEOMETRICOVERLAPPINGOPERATOR_DECL_HPP
