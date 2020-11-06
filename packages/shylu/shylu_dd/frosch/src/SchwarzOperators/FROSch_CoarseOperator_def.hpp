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

#ifndef _FROSCH_COARSEOPERATOR_DEF_HPP
#define _FROSCH_COARSEOPERATOR_DEF_HPP

#include <FROSch_CoarseOperator_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template<class SC,class LO,class GO,class NO>
    CoarseOperator<SC,LO,GO,NO>::CoarseOperator(ConstXMatrixPtr k,
                                                ParameterListPtr parameterList) :
    SchwarzOperator<SC,LO,GO,NO> (k,parameterList),
    CoarseSpace_ (new CoarseSpace<SC,LO,GO,NO>(this->MpiComm_,this->SerialComm_)),
    DistributionList_ (sublist(parameterList,"Distribution"))
    {
        FROSCH_DETAILTIMER_START_LEVELID(coarseOperatorTime,"CoarseOperator::CoarseOperator");
    }

    template<class SC,class LO,class GO,class NO>
    CoarseOperator<SC,LO,GO,NO>::~CoarseOperator()
    {
        CoarseSolver_.reset();
    }

    template <class SC,class LO,class GO,class NO>
    int CoarseOperator<SC,LO,GO,NO>::compute()
    {
        FROSCH_TIMER_START_LEVELID(computeTime,"CoarseOperator::compute");
        FROSCH_ASSERT(this->IsInitialized_,"FROSch::CoarseOperator : ERROR: CoarseOperator has to be initialized before calling compute()");
        // This is not optimal yet... Some work could be moved to Initialize
        //if (this->Verbose_) cout << "FROSch::CoarseOperator : WARNING: Some of the operations could probably be moved from initialize() to Compute().\n";

        bool reuseCoarseBasis = this->ParameterList_->get("Reuse: Coarse Basis",true);
        bool reuseCoarseMatrix = this->ParameterList_->get("Reuse: Coarse Matrix",false);
        if (!this->IsComputed_) {
            reuseCoarseBasis = false;
            reuseCoarseMatrix = false;
        }

        if (!reuseCoarseBasis) {
            if (this->IsComputed_ && this->Verbose_) cout << "FROSch::CoarseOperator : Recomputing the Coarse Basis" << endl;
            clearCoarseSpace(); // AH 12/11/2018: If we do not clear the coarse space, we will always append just append the coarse space
            XMapPtr subdomainMap = this->computeCoarseSpace(CoarseSpace_); // AH 12/11/2018: This map could be overlapping, repeated, or unique. This depends on the specific coarse operator
            if (CoarseSpace_->hasUnassembledMaps()) { // If there is no unassembled basis, the current Phi_ should already be correct
                CoarseSpace_->assembleCoarseSpace();
                FROSCH_ASSERT(CoarseSpace_->hasAssembledBasis(),"FROSch::CoarseOperator : !CoarseSpace_->hasAssembledBasis()");
                CoarseSpace_->buildGlobalBasisMatrix(this->K_->getRowMap(),this->K_->getRangeMap(),subdomainMap,this->ParameterList_->get("Threshold Phi",1.e-8));
                FROSCH_ASSERT(CoarseSpace_->hasGlobalBasisMatrix(),"FROSch::CoarseOperator : !CoarseSpace_->hasGlobalBasisMatrix()");
                Phi_ = CoarseSpace_->getGlobalBasisMatrix();
            }
        }
        if (!reuseCoarseMatrix) {
            if (this->IsComputed_ && this->Verbose_) cout << "FROSch::CoarseOperator : Recomputing the Coarse Matrix" << endl;
            this->setUpCoarseOperator();
        }
        this->IsComputed_ = true;

        // Store current Phi in ParameterList_
        if ( this->ParameterList_->get("Store Phi",false) ) {
            FROSCH_NOTIFICATION("FROSch::CoarseOperator",this->Verbose_,"Storing current Phi in Parameterlist.");
            this->ParameterList_->set("RCP(Phi)", Phi_);
        }

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int CoarseOperator<SC,LO,GO,NO>::clearCoarseSpace()
    {
        return CoarseSpace_->clearCoarseSpace();
    }

    template<class SC,class LO,class GO,class NO>
    void CoarseOperator<SC,LO,GO,NO>::apply(const XMultiVector &x,
                                            XMultiVector &y,
                                            bool usePreconditionerOnly,
                                            ETransp mode,
                                            SC alpha,
                                            SC beta) const
    {
        FROSCH_TIMER_START_LEVELID(applyTime,"CoarseOperator::apply");
        static int i = 0;
        if (!Phi_.is_null() && this->IsComputed_) {
            if (XTmp_.is_null()) XTmp_ = MultiVectorFactory<SC,LO,GO,NO>::Build(x.getMap(),x.getNumVectors());
            *XTmp_ = x;
            if (!usePreconditionerOnly && mode == NO_TRANS) {
                this->K_->apply(x,*XTmp_,mode,ScalarTraits<SC>::one(),ScalarTraits<SC>::zero());
            }
            if (XCoarseSolve_.is_null()) XCoarseSolve_ = MultiVectorFactory<SC,LO,GO,NO>::Build(GatheringMaps_[GatheringMaps_.size()-1],x.getNumVectors());
            else XCoarseSolve_->replaceMap(GatheringMaps_[GatheringMaps_.size()-1]); // The map is replaced in applyCoarseSolve(). If we do not build it from scratch, we should at least replace the map here. This may be important since the maps live on different communicators.
            if (YCoarseSolve_.is_null()) YCoarseSolve_ = MultiVectorFactory<SC,LO,GO,NO>::Build(GatheringMaps_[GatheringMaps_.size()-1],y.getNumVectors());
            applyPhiT(*XTmp_,*XCoarseSolve_);
            applyCoarseSolve(*XCoarseSolve_,*YCoarseSolve_,mode);
            applyPhi(*YCoarseSolve_,*XTmp_);
            if (!usePreconditionerOnly && mode != NO_TRANS) {
                this->K_->apply(*XTmp_,*XTmp_,mode,ScalarTraits<SC>::one(),ScalarTraits<SC>::zero());
            }
            y.update(alpha,*XTmp_,beta);
        } else {
            if (i==0) {
                FROSCH_WARNING("FROSch::CoarseOperator",(this->Verbose_ && Phi_.is_null()),"Coarse Basis is empty => The CoarseOperator will just act as the identity...");
                FROSCH_WARNING("FROSch::CoarseOperator",(this->Verbose_ && !this->IsComputed_),"CoarseOperator has not been computed yet => The CoarseOperator will just act as the identity...");
                i++;
            }
            y.update(ScalarTraits<SC>::one(),x,ScalarTraits<SC>::zero());
        }
    }

    template<class SC,class LO,class GO,class NO>
    void CoarseOperator<SC,LO,GO,NO>::applyPhiT(const XMultiVector& x,
                                                XMultiVector& y) const
    {
        FROSCH_DETAILTIMER_START_LEVELID(applyPhiTTime,"CoarseOperator::applyPhiT");
        // AH 08/22/2019 TODO: We cannot ger rid of the Build() calls because of "XCoarse_ = XCoarseSolveTmp_;". This is basically caused by the whole Gathering Map strategy. As soon as we have replaced this, we can get rid of the Build() calls
        XCoarse_ = MultiVectorFactory<SC,LO,GO,NO>::Build(CoarseSpace_->getBasisMapUnique(),x.getNumVectors()); // AH 08/22/2019 TODO: Can we get rid of this? If possible, we should remove the whole GatheringMaps idea and replace it by some smart all-to-all MPI communication
        {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
            FROSCH_DETAILTIMER_START_LEVELID(applyTime,"apply");
#endif
            Phi_->apply(x,*XCoarse_,TRANS);
        }
        for (UN j=0; j<GatheringMaps_.size(); j++) {
            XCoarseSolveTmp_ = MultiVectorFactory<SC,LO,GO,NO>::Build(GatheringMaps_[j],x.getNumVectors()); // AH 08/22/2019 TODO: Can we get rid of this? If possible, we should remove the whole GatheringMaps idea and replace it by some smart all-to-all MPI communication
            {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                FROSCH_DETAILTIMER_START_LEVELID(applyTime,"doExport");
#endif
                XCoarseSolveTmp_->doExport(*XCoarse_,*CoarseSolveExporters_[j],ADD);
            }
            XCoarse_ = XCoarseSolveTmp_;
        }
        y = *XCoarseSolveTmp_;
    }

    template<class SC,class LO,class GO,class NO>
    void CoarseOperator<SC,LO,GO,NO>::applyCoarseSolve(XMultiVector& x,
                                                       XMultiVector& y,
                                                       ETransp mode) const
    {
        FROSCH_DETAILTIMER_START_LEVELID(applyCoarseSolveTime,"CoarseOperator::applyCoarseSolve");
        if (OnCoarseSolveComm_) {
            x.replaceMap(CoarseSolveMap_);
            if (YTmp_.is_null()) YTmp_ = MultiVectorFactory<SC,LO,GO,NO>::Build(CoarseSolveMap_,x.getNumVectors());
            else YTmp_->replaceMap(CoarseSolveMap_); // The map is replaced later in this function. If we do not build it from scratch, we should at least replace the map here. This may be important since the maps live on different communicators.
            CoarseSolver_->apply(x,*YTmp_,mode);
        } else {
            if (YTmp_.is_null()) YTmp_ = MultiVectorFactory<SC,LO,GO,NO>::Build(CoarseSolveMap_,x.getNumVectors());
            else YTmp_->replaceMap(CoarseSolveMap_); // The map is replaced later in this function. If we do not build it from scratch, we should at least replace the map here. This may be important since the maps live on different communicators.
        }
        YTmp_->replaceMap(GatheringMaps_[GatheringMaps_.size()-1]);
        y = *YTmp_;
    }

    template<class SC,class LO,class GO,class NO>
    void CoarseOperator<SC,LO,GO,NO>::applyPhi(const XMultiVector& x,
                                               XMultiVector& y) const
    {
        FROSCH_DETAILTIMER_START_LEVELID(applyPhiTime,"CoarseOperator::applyPhi");
        // AH 08/22/2019 TODO: We have the same issue here as in applyPhiT()
        YCoarseSolveTmp_ = MultiVectorFactory<SC,LO,GO,NO>::Build(x.getMap(),x.getNumVectors());
        *YCoarseSolveTmp_ = x;
        for (int j=GatheringMaps_.size()-1; j>0; j--) {
            YCoarse_ = MultiVectorFactory<SC,LO,GO,NO>::Build(GatheringMaps_[j-1],x.getNumVectors());
            {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                FROSCH_DETAILTIMER_START_LEVELID(applyTime,"doImport");
#endif
#ifdef FROSCH_COARSEOPERATOR_EXPORT_AND_IMPORT
                YCoarse_->doImport(*YCoarseSolveTmp_,*CoarseSolveImporters_[j],INSERT);
#else
                YCoarse_->doImport(*YCoarseSolveTmp_,*CoarseSolveExporters_[j],INSERT);
#endif
            }
            YCoarseSolveTmp_ = YCoarse_;
        }
        YCoarse_ = MultiVectorFactory<SC,LO,GO,NO>::Build(CoarseSpace_->getBasisMapUnique(),x.getNumVectors());
        {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
            FROSCH_DETAILTIMER_START_LEVELID(applyTime,"doImport");
#endif
#ifdef FROSCH_COARSEOPERATOR_EXPORT_AND_IMPORT
            YCoarse_->doImport(*YCoarseSolveTmp_,*CoarseSolveImporters_[0],INSERT);
#else
            YCoarse_->doImport(*YCoarseSolveTmp_,*CoarseSolveExporters_[0],INSERT);
#endif
        }
        {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
            FROSCH_DETAILTIMER_START_LEVELID(applyTime,"apply");
#endif
            Phi_->apply(*YCoarse_,y,NO_TRANS);
        }
    }

    template<class SC,class LO,class GO,class NO>
    typename CoarseOperator<SC,LO,GO,NO>::CoarseSpacePtr CoarseOperator<SC,LO,GO,NO>::getCoarseSpace() const
    {
        return CoarseSpace_;
    }

    template<class SC,class LO,class GO,class NO>
    int CoarseOperator<SC,LO,GO,NO>::setUpCoarseOperator()
    {
        FROSCH_DETAILTIMER_START_LEVELID(setUpCoarseOperatorTime,"CoarseOperator::setUpCoarseOperator");
        if (!Phi_.is_null()) {
            // Build CoarseMatrix_
            XMatrixPtr k0 = buildCoarseMatrix();

            //------------------------------------------------------------------------------------------------------------------------
            // Communicate coarse matrix
            FROSCH_DETAILTIMER_START_LEVELID(communicateCoarseMatrixTime,"communicate coarse matrix");
            if (!DistributionList_->get("Type","linear").compare("linear")) {
                XMatrixPtr tmpCoarseMatrix = MatrixFactory<SC,LO,GO,NO>::Build(GatheringMaps_[0]);
                {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                    FROSCH_DETAILTIMER_START_LEVELID(coarseMatrixExportTime,"Export Coarse Matrix");
#endif
                    tmpCoarseMatrix->doExport(*k0,*CoarseSolveExporters_[0],INSERT);
                }

                for (UN j=1; j<GatheringMaps_.size(); j++) {
                    tmpCoarseMatrix->fillComplete();
                    k0 = tmpCoarseMatrix;
                    tmpCoarseMatrix = MatrixFactory<SC,LO,GO,NO>::Build(GatheringMaps_[j]);
                    {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                        FROSCH_DETAILTIMER_START_LEVELID(coarseMatrixExportTime,"Export Coarse Matrix");
#endif
                        tmpCoarseMatrix->doExport(*k0,*CoarseSolveExporters_[j],INSERT);
                    }
                }
                k0 = tmpCoarseMatrix;

            } else if (!DistributionList_->get("Type","linear").compare("ZoltanDual")) {
                //ZoltanDual
                CoarseSolveExporters_[0] = Xpetra::ExportFactory<LO,GO,NO>::Build(k0->getMap(),GatheringMaps_[0]);
                XMatrixPtr tmpCoarseMatrix = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(GatheringMaps_[0],k0->getGlobalMaxNumRowEntries());
                tmpCoarseMatrix->doExport(*k0,*CoarseSolveExporters_[0],Xpetra::INSERT);

                for (UN j=1; j<GatheringMaps_.size(); j++) {
                    tmpCoarseMatrix->fillComplete();
                    k0 = tmpCoarseMatrix;
                    CoarseSolveExporters_[j] = Xpetra::ExportFactory<LO,GO,NO>::Build(GatheringMaps_[j-1],GatheringMaps_[j]);
                    tmpCoarseMatrix = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(GatheringMaps_[j],k0->getGlobalMaxNumRowEntries());
                    tmpCoarseMatrix->doExport(*k0,*CoarseSolveExporters_[j],Xpetra::INSERT);
                }

                tmpCoarseMatrix->fillComplete();
                k0 = tmpCoarseMatrix;

            } else if (!DistributionList_->get("Type","linear").compare("Zoltan2")) {
#ifdef HAVE_SHYLU_DDFROSCH_ZOLTAN2
                GatheringMaps_[0] = rcp_const_cast<XMap> (BuildUniqueMap(k0->getRowMap()));
                CoarseSolveExporters_[0] = ExportFactory<LO,GO,NO>::Build(CoarseSpace_->getBasisMapUnique(),GatheringMaps_[0]);

                if (NumProcsCoarseSolve_ < this->MpiComm_->getSize()) {
                    XMatrixPtr k0Unique = MatrixFactory<SC,LO,GO,NO>::Build(GatheringMaps_[0]);
                    k0Unique->doExport(*k0,*CoarseSolveExporters_[0],INSERT);
                    k0Unique->fillComplete(GatheringMaps_[0],GatheringMaps_[0]);

                    if (NumProcsCoarseSolve_<this->MpiComm_->getSize()) {
                        ParameterListPtr tmpList = sublist(DistributionList_,"Zoltan2 Parameter");
                        tmpList->set("num_global_parts",NumProcsCoarseSolve_);
                        FROSch::RepartionMatrixZoltan2(k0Unique,tmpList);
                    }

                    k0 = k0Unique;
                    GatheringMaps_[0] = k0->getRowMap();
                    CoarseSolveExporters_[0] = ExportFactory<LO,GO,NO>::Build(CoarseSpace_->getBasisMapUnique(),GatheringMaps_[0]);

                    if (GatheringMaps_[0]->getNodeNumElements()>0) {
                        OnCoarseSolveComm_=true;
                    }
                    CoarseSolveComm_ = this->MpiComm_->split(!OnCoarseSolveComm_,this->MpiComm_->getRank());
                    CoarseSolveMap_ = MapFactory<LO,GO,NO>::Build(CoarseSpace_->getBasisMapUnique()->lib(),-1,GatheringMaps_[0]->getNodeElementList(),0,CoarseSolveComm_);
                }
#else
                ThrowErrorMissingPackage("FROSch::CoarseOperator","Zoltan2");
#endif
                //------------------------------------------------------------------------------------------------------------------------
            } else {
                FROSCH_ASSERT(false,"Distribution Type unknown!");
            }
            FROSCH_DETAILTIMER_STOP(communicateCoarseMatrixTime);

            //------------------------------------------------------------------------------------------------------------------------
            // Matrix to the new communicator
            if (OnCoarseSolveComm_) {
                FROSCH_DETAILTIMER_START_LEVELID(replicateCoarseMatrixOnCoarseCommTime,"replicate coarse matrix on coarse comm");
                LO numRows = k0->getNodeNumRows();
                ArrayRCP<size_t> elemsPerRow(numRows);
                LO numDiagonalsAdded = 0;
                if (k0->isFillComplete()) {
                    ConstLOVecView indices;
                    ConstSCVecView values;
                    for (LO i = 0; i < numRows; i++) {
                        size_t numEntries;
                        numEntries = k0->getNumEntriesInLocalRow(i);
                        if (numEntries == 0) {
                            //Always add the diagonal for empty rows
                            numEntries = 1;
                        }
                        elemsPerRow[i] = numEntries;
                    }
                    CoarseMatrix_ = MatrixFactory<SC,LO,GO,NO>::Build(CoarseSolveMap_,elemsPerRow);
                    for (LO i = 0; i < numRows; i++) {
                        GO globalRow = CoarseSolveMap_->getGlobalElement(i);
                        k0->getLocalRowView(i,indices,values);
                        if (indices.size()>0) {
                            GOVec indicesGlob(indices.size());
                            for (UN j=0; j<indices.size(); j++) {
                                indicesGlob[j] = k0->getColMap()->getGlobalElement(indices[j]);
                            }
                            CoarseMatrix_->insertGlobalValues(globalRow,indicesGlob(),values);
                        } else { // Add diagonal unit for zero rows // Todo: Do you we need to sort the coarse matrix "NodeWise"?
                            GOVec indicesGlob(1,CoarseSolveMap_->getGlobalElement(i));
                            SCVec values(1,ScalarTraits<SC>::one());
                            CoarseMatrix_->insertGlobalValues(globalRow,indicesGlob(),values());
                            numDiagonalsAdded++;
                        }
                    }
                    CoarseMatrix_->fillComplete(CoarseSolveMap_,CoarseSolveMap_); //RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout)); CoarseMatrix_->describe(*fancy,VERB_EXTREME);
                } else {
                    ConstGOVecView indices;
                    ConstSCVecView values;
                    for (LO i = 0; i < numRows; i++) {
                        GO globalRow = CoarseSolveMap_->getGlobalElement(i);
                        size_t numEntries;
                        numEntries = k0->getNumEntriesInGlobalRow(globalRow);
                        if (numEntries == 0) {
                            //Always add the diagonal for empty rows
                            numEntries = 1;
                        }
                        elemsPerRow[i] = numEntries;
                    }
                    CoarseMatrix_ = MatrixFactory<SC,LO,GO,NO>::Build(CoarseSolveMap_,elemsPerRow);
                    for (LO i = 0; i < numRows; i++) {
                        GO globalRow = CoarseSolveMap_->getGlobalElement(i);
                        k0->getGlobalRowView(globalRow,indices,values);
                        if (indices.size()>0) {
                            CoarseMatrix_->insertGlobalValues(globalRow,indices,values);
                        } else { // Add diagonal unit for zero rows // Todo: Do you we need to sort the coarse matrix "NodeWise"?
                            GOVec indices(1,globalRow);
                            SCVec values(1,ScalarTraits<SC>::one());
                            CoarseMatrix_->insertGlobalValues(globalRow,indices(),values());
                            numDiagonalsAdded++;
                        }
                    }
                    CoarseMatrix_->fillComplete(CoarseSolveMap_,CoarseSolveMap_); //RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout)); CoarseMatrix_->describe(*fancy,VERB_EXTREME);
                }
                FROSCH_DETAILTIMER_STOP(replicateCoarseMatrixOnCoarseCommTime);

                FROSCH_DETAILTIMER_START_LEVELID(printStatisticsTime,"print statistics");
                // Statistics on adding diagonal entries
                GOVec globalVec(5);
                LOVec localVec(5);
                LOVec sumVec(5);
                SCVec avgVec(5);
                LOVec minVec(5);
                LOVec maxVec(5);

                globalVec[0] = CoarseMatrix_->getGlobalNumRows();
                localVec[0] = CoarseMatrix_->getNodeNumRows();
                reduceAll(*CoarseSolveComm_,REDUCE_SUM,localVec[0],ptr(&sumVec[0]));
                avgVec[0] = max(sumVec[0]/double(CoarseSolveComm_->getSize()),0.0);
                reduceAll(*CoarseSolveComm_,REDUCE_MIN,localVec[0],ptr(&minVec[0]));
                reduceAll(*CoarseSolveComm_,REDUCE_MAX,localVec[0],ptr(&maxVec[0]));

                globalVec[1] = CoarseMatrix_->getGlobalNumEntries();
                localVec[1] = CoarseMatrix_->getNodeNumEntries();
                reduceAll(*CoarseSolveComm_,REDUCE_SUM,localVec[1],ptr(&sumVec[1]));
                avgVec[1] = max(sumVec[1]/double(CoarseSolveComm_->getSize()),0.0);
                reduceAll(*CoarseSolveComm_,REDUCE_MIN,localVec[1],ptr(&minVec[1]));
                reduceAll(*CoarseSolveComm_,REDUCE_MAX,localVec[1],ptr(&maxVec[1]));

                globalVec[2] = double(globalVec[1])/double(globalVec[0]);
                localVec[2] = double(localVec[1])/double(localVec[0]);
                reduceAll(*CoarseSolveComm_,REDUCE_SUM,localVec[2],ptr(&sumVec[2]));
                avgVec[2] = max(sumVec[2]/double(CoarseSolveComm_->getSize()),0.0);
                reduceAll(*CoarseSolveComm_,REDUCE_MIN,localVec[2],ptr(&minVec[2]));
                reduceAll(*CoarseSolveComm_,REDUCE_MAX,localVec[2],ptr(&maxVec[2]));

                localVec[3] = CoarseMatrix_->getNodeMaxNumRowEntries();
                reduceAll(*CoarseSolveComm_,REDUCE_SUM,localVec[3],ptr(&sumVec[3]));
                avgVec[3] = max(sumVec[3]/double(CoarseSolveComm_->getSize()),0.0);
                reduceAll(*CoarseSolveComm_,REDUCE_MIN,localVec[3],ptr(&minVec[3]));
                reduceAll(*CoarseSolveComm_,REDUCE_MAX,localVec[3],ptr(&maxVec[3]));

                localVec[4] = numDiagonalsAdded;
                reduceAll(*CoarseSolveComm_,REDUCE_SUM,localVec[4],ptr(&sumVec[4]));
                avgVec[4] = max(sumVec[4]/double(CoarseSolveComm_->getSize()),0.0);
                reduceAll(*CoarseSolveComm_,REDUCE_MIN,localVec[4],ptr(&minVec[4]));
                reduceAll(*CoarseSolveComm_,REDUCE_MAX,localVec[4],ptr(&maxVec[4]));

                if (CoarseSolveComm_->getRank() == 0) {
                    cout
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << setw(89) << "-----------------------------------------------------------------------------------------"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << "| "
                    << left << setw(74) << "> Coarse Problem Statistics (Coarse Communicator) " << right << setw(8) << "(Level " << setw(2) << this->LevelID_ << ")"
                    << " |"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << setw(89) << "========================================================================================="
                    // << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    // << "| " << left << setw(41) << "Dimension of the coarse problem" << right
                    // << " | " << setw(41) << dimCoarseProblem
                    // << " |"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << "| " << left << setw(41) << "Number of ranks on the coarse comm" << right
                    << " | " << setw(41) << NumProcsCoarseSolve_
                    << " |"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << "| " << left << setw(41) << "Solver type" << right
                    << " | " << setw(41) << this->ParameterList_->sublist("CoarseSolver").get("SolverType","Amesos")
                    << " |"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << "| " << left << setw(41) << "Solver" << right
                    << " | " << setw(41) << this->ParameterList_->sublist("CoarseSolver").get("Solver","Mumps")
                    << " |"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << "| " << left << setw(41) << "Reuse symbolic factorization" << right
                    << " | " << setw(41) << this->ParameterList_->get("Reuse: Coarse Matrix Symbolic Factorization",true)
                    << " |"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << "| " << left << setw(41) << "Reuse coarse basis" << right
                    << " | " << setw(41) << this->ParameterList_->get("Reuse: Coarse Basis",true)
                    << " |"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << "| " << left << setw(41) << "Reuse coarse matrix" << right
                    << " | " << setw(41) << this->ParameterList_->get("Reuse: Coarse Matrix",false)
                    << " |"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << setw(89) << "-----------------------------------------------------------------------------------------"
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
                    << "| " << left << setw(20) << "Number of rows" << right
                    << " | " << setw(10) << globalVec[0]
                    << " | " << setw(10) << setprecision(5) << avgVec[0]
                    << " | " << setw(10) << minVec[0]
                    << " | " << setw(10) << maxVec[0]
                    << " | " << setw(10) << sumVec[0]
                    << " |"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << "| " << left << setw(20) << "Entries" << right
                    << " | " << setw(10) << globalVec[1]
                    << " | " << setw(10) << setprecision(5) << avgVec[1]
                    << " | " << setw(10) << minVec[1]
                    << " | " << setw(10) << maxVec[1]
                    << " | " << setw(10) << sumVec[1]
                    << " |"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << "| " << left << setw(20) << "Avg entries per row" << right
                    << " | " << setw(10) << globalVec[2]
                    << " | " << setw(10) << setprecision(5) << avgVec[2]
                    << " | " << setw(10) << minVec[2]
                    << " | " << setw(10) << maxVec[2]
                    << " | " << setw(10) << sumVec[2]
                    << " |"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << "| " << left << setw(20) << "Max entries per row" << right
                    << " | " << setw(10) << " "
                    << " | " << setw(10) << setprecision(5) << avgVec[3]
                    << " | " << setw(10) << minVec[3]
                    << " | " << setw(10) << maxVec[3]
                    << " | " << setw(10) << " "
                    << " |"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << "| " << left << setw(20) << "Unit diagonals added" << right
                    << " | " << setw(10) << sumVec[4]
                    << " | " << setw(10) << setprecision(5) << avgVec[4]
                    << " | " << setw(10) << minVec[4]
                    << " | " << setw(10) << maxVec[4]
                    << " | " << setw(10) << sumVec[4]
                    << " |"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << setw(89) << "-----------------------------------------------------------------------------------------"
                    << endl;
                }
                FROSCH_DETAILTIMER_STOP(printStatisticsTime);

                bool reuseCoarseMatrixSymbolicFactorization = this->ParameterList_->get("Reuse: Coarse Matrix Symbolic Factorization",true);
                if (!this->IsComputed_) {
                    reuseCoarseMatrixSymbolicFactorization = false;
                }
                if (!reuseCoarseMatrixSymbolicFactorization) {
                    if (this->IsComputed_ && this->Verbose_) cout << "FROSch::CoarseOperator : Recomputing the Symbolic Factorization of the coarse matrix" << endl;
                    CoarseSolver_.reset(new SubdomainSolver<SC,LO,GO,NO>(CoarseMatrix_,
                                                                         sublist(this->ParameterList_,"CoarseSolver"),
                                                                         string("CoarseSolver (Level ") + to_string(this->LevelID_) + string(")")));
                    CoarseSolver_->initialize();
                } else {
                    FROSCH_ASSERT(!CoarseSolver_.is_null(),"FROSch::CoarseOperator : ERROR: CoarseSolver_.is_null()");
                    CoarseSolver_->resetMatrix(CoarseMatrix_.getConst(),true);
                }
                CoarseSolver_->compute();
            }
        } else {
            FROSCH_WARNING("FROSch::CoarseOperator",this->Verbose_,"No coarse basis has been set up. Neglecting CoarseOperator.");
        }
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    typename CoarseOperator<SC,LO,GO,NO>::XMatrixPtr CoarseOperator<SC,LO,GO,NO>::buildCoarseMatrix()
    {
        FROSCH_DETAILTIMER_START_LEVELID(buildCoarseMatrixTime,"CoarseOperator::buildCoarseMatrix");
        XMatrixPtr k0;
        if (this->ParameterList_->get("Use Triple MatrixMultiply",false)) {
            k0 = MatrixFactory<SC,LO,GO,NO>::Build(CoarseSpace_->getBasisMapUnique(),as<LO>(0));
            TripleMatrixMultiply<SC,LO,GO,NO>::MultiplyRAP(*Phi_,true,*this->K_,false,*Phi_,false,*k0);
        } else {
            RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout)); //Phi_->describe(*fancy,VERB_EXTREME);
            XMatrixPtr tmp = MatrixMatrix<SC,LO,GO,NO>::Multiply(*this->K_,false,*Phi_,false,*fancy);
            k0 = MatrixMatrix<SC,LO,GO,NO>::Multiply(*Phi_,true,*tmp,false,*fancy); //k0->describe(*fancy,VERB_EXTREME);
        }
        return k0;
    }

    template<class SC,class LO,class GO,class NO>
    int CoarseOperator<SC,LO,GO,NO>::buildCoarseSolveMap(ConstXMapPtr coarseMapUnique)
    {
        FROSCH_DETAILTIMER_START_LEVELID(buildCoarseSolveMapTime,"CoarseOperator::buildCoarseSolveMap");
        NumProcsCoarseSolve_ = DistributionList_->get("NumProcs",1);
        double factor = DistributionList_->get("Factor",0.0);

        switch (NumProcsCoarseSolve_) {
            case -1:
                FROSCH_ASSERT(false,"We do not know the size of the matrix yet. Therefore, we cannot use the formula NumProcsCoarseSolve_ = int(0.5*(1+max(k0->getGlobalNumRows()/10000,k0->getGlobalNumEntries()/100000)));");
                //NumProcsCoarseSolve_ = int(0.5*(1+max(k0->getGlobalNumRows()/10000,k0->getGlobalNumEntries()/100000)));
                break;

            case 0:
                NumProcsCoarseSolve_ = this->MpiComm_->getSize();
                break;

            default:
                if (NumProcsCoarseSolve_>this->MpiComm_->getSize()) NumProcsCoarseSolve_ = this->MpiComm_->getSize();
                if (fabs(factor) > 1.0e-12) NumProcsCoarseSolve_ = int(NumProcsCoarseSolve_/factor);
                if (NumProcsCoarseSolve_<1) NumProcsCoarseSolve_ = 1;
                break;
        }

        if (!DistributionList_->get("Type","linear").compare("linear")) {

            int gatheringSteps = DistributionList_->get("GatheringSteps",1);
            GatheringMaps_.resize(gatheringSteps);
            CoarseSolveExporters_.resize(gatheringSteps);
#ifdef FROSCH_COARSEOPERATOR_EXPORT_AND_IMPORT
            CoarseSolveImporters_.resize(gatheringSteps);
#endif

            LO numProcsGatheringStep = this->MpiComm_->getSize();
            GO numGlobalIndices = coarseMapUnique->getMaxAllGlobalIndex()+1;
            int numMyRows;
            double gatheringFactor = pow(double(this->MpiComm_->getSize())/double(NumProcsCoarseSolve_),1.0/double(gatheringSteps));

            for (int i=0; i<gatheringSteps-1; i++) {
                numMyRows = 0;
                numProcsGatheringStep = LO(numProcsGatheringStep/gatheringFactor);
                //if (this->Verbose_) cout << i << " " << numProcsGatheringStep << " " << numGlobalIndices << endl;
                if (this->MpiComm_->getRank()%(this->MpiComm_->getSize()/numProcsGatheringStep) == 0 && this->MpiComm_->getRank()/(this->MpiComm_->getSize()/numProcsGatheringStep) < numProcsGatheringStep) {
                    if (this->MpiComm_->getRank()==0) {
                        numMyRows = numGlobalIndices - (numGlobalIndices/numProcsGatheringStep)*(numProcsGatheringStep-1);
                    } else {
                        numMyRows = numGlobalIndices/numProcsGatheringStep;
                    }
                }
                {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                    FROSCH_DETAILTIMER_START_LEVELID(gatheringMapsTime,"Gathering Maps");
#endif
                    GatheringMaps_[i] = MapFactory<LO,GO,NO>::Build(coarseMapUnique->lib(),-1,numMyRows,0,this->MpiComm_);
                }
            }

            numMyRows = 0;
            if (this->MpiComm_->getRank()%(this->MpiComm_->getSize()/NumProcsCoarseSolve_) == 0 && this->MpiComm_->getRank()/(this->MpiComm_->getSize()/NumProcsCoarseSolve_) < NumProcsCoarseSolve_) {
                if (this->MpiComm_->getRank()==0) {
                    numMyRows = numGlobalIndices - (numGlobalIndices/NumProcsCoarseSolve_)*(NumProcsCoarseSolve_-1);
                } else {
                    numMyRows = numGlobalIndices/NumProcsCoarseSolve_;
                }
            }
            {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                FROSCH_DETAILTIMER_START_LEVELID(gatheringMapsTime,"Gathering Maps");
#endif
                GatheringMaps_[gatheringSteps-1] = MapFactory<LO,GO,NO>::Build(coarseMapUnique->lib(),-1,numMyRows,0,this->MpiComm_);
            }
            //cout << *GatheringMaps_->at(gatheringSteps-1);

            //------------------------------------------------------------------------------------------------------------------------
            // Use a separate Communicator for the coarse problem
            if (GatheringMaps_[GatheringMaps_.size()-1]->getNodeNumElements()>0) {
                OnCoarseSolveComm_=true;
            }
            {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                FROSCH_DETAILTIMER_START_LEVELID(commSplitTime,"Coarse Communicator Split");
#endif
                CoarseSolveComm_ = this->MpiComm_->split(!OnCoarseSolveComm_,this->MpiComm_->getRank());
            }
            {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                FROSCH_DETAILTIMER_START_LEVELID(coarseCommMapTime,"Coarse Communicator Map");
#endif
                CoarseSolveMap_ = MapFactory<LO,GO,NO>::Build(coarseMapUnique->lib(),-1,GatheringMaps_[GatheringMaps_.size()-1]->getNodeElementList(),0,CoarseSolveComm_);
            }

            // Possibly change the Send type for this Exporter
            ParameterListPtr gatheringCommunicationList = sublist(DistributionList_,"Gathering Communication");
            // Set communication type "Alltoall" if not specified differently
            if (!gatheringCommunicationList->isParameter("Send type")) gatheringCommunicationList->set("Send type","Send");

            // Create Import and Export objects
            {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                FROSCH_DETAILTIMER_START_LEVELID(coarseSolveExportersTime,"Build Exporters");
#endif
                CoarseSolveExporters_[0] = ExportFactory<LO,GO,NO>::Build(coarseMapUnique,GatheringMaps_[0]);
                CoarseSolveExporters_[0]->setDistributorParameters(gatheringCommunicationList); // Set the parameter list for the communication of the exporter
            }
#ifdef FROSCH_COARSEOPERATOR_EXPORT_AND_IMPORT
            {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                FROSCH_DETAILTIMER_START_LEVELID(coarseSolveImportersTime,"Build Importers");
#endif
                CoarseSolveImporters_[0] = ImportFactory<LO,GO,NO>::Build(GatheringMaps_[0],coarseMapUnique);
                CoarseSolveImporters_[0]->setDistributorParameters(gatheringCommunicationList); // Set the parameter list for the communication of the exporter
            }
#endif

            for (UN j=1; j<GatheringMaps_.size(); j++) {
                {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                    FROSCH_DETAILTIMER_START_LEVELID(coarseSolveExportersTime,"Build Exporters");
#endif
                    CoarseSolveExporters_[j] = ExportFactory<LO,GO,NO>::Build(GatheringMaps_[j-1],GatheringMaps_[j]);
                    CoarseSolveExporters_[j]->setDistributorParameters(gatheringCommunicationList); // Set the parameter list for the communication of the exporter
                }
#ifdef FROSCH_COARSEOPERATOR_EXPORT_AND_IMPORT
                {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                    FROSCH_DETAILTIMER_START_LEVELID(coarseSolveImportersTime,"Build Importers");
#endif
                    CoarseSolveImporters_[j] = ImportFactory<LO,GO,NO>::Build(GatheringMaps_[j],GatheringMaps_[j-1]);
                    CoarseSolveImporters_[j]->setDistributorParameters(gatheringCommunicationList); // Set the parameter list for the communication of the exporter
                }
#endif
            }
        } else if (!DistributionList_->get("Type","linear").compare("ZoltanDual")) {
            //ZoltanDual provides a partition of the coarse problem with Zoltan2 inlcuding the
            //build of a Repeated map suited for the next level
            //GatheringSteps to communicate Matrix
            int gatheringSteps = DistributionList_->get("GatheringSteps",1);
            GatheringMaps_.resize(gatheringSteps);
            CoarseSolveExporters_.resize(gatheringSteps);

            double gatheringFactor = pow(double(this->MpiComm_->getSize())/double(NumProcsCoarseSolve_),1.0/double(gatheringSteps));
            LO numProcsGatheringStep = this->MpiComm_->getSize();
            GO numGlobalIndices = CoarseMap_->getMaxAllGlobalIndex();
            GO numMyRows;
            numMyRows = 0;

            if (this->MpiComm_->getRank()%(this->MpiComm_->getSize()/NumProcsCoarseSolve_) == 0 && this->MpiComm_->getRank()/(this->MpiComm_->getSize()/NumProcsCoarseSolve_) < NumProcsCoarseSolve_) {
                if (this->MpiComm_->getRank()==0) {
                    numMyRows = numGlobalIndices - (numGlobalIndices/NumProcsCoarseSolve_)*(NumProcsCoarseSolve_-1);
                } else {
                    numMyRows = numGlobalIndices/NumProcsCoarseSolve_;
                }
            }

            XMapPtr tmpCoarseMap = Xpetra::MapFactory<LO,GO,NO>::Build(CoarseMap_->lib(),-1,numMyRows,0,this->MpiComm_);
            if (tmpCoarseMap->getNodeNumElements()>0) {
                OnCoarseSolveComm_=true;
            }
            CoarseSolveComm_ = this->MpiComm_->split(!OnCoarseSolveComm_,this->MpiComm_->getRank());

            //Gathering Steps for RepeatedMap#################################################
            //-> Have to test that
            RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout));

            int MLgatheringSteps = DistributionList_->get("MLGatheringSteps",2);
            MLGatheringMaps_.resize(MLgatheringSteps);
            MLCoarseSolveExporters_.resize(MLgatheringSteps);
            double MLgatheringFactor = pow(double(this->MpiComm_->getSize())/double(NumProcsCoarseSolve_),1.0/double(MLgatheringSteps));
            GO MLnumGlobalIndices = SubdomainConnectGraph_->getRowMap()->getMaxAllGlobalIndex()+1;
            GO MLnumMyRows;

            MLGatheringMaps_[0] =  Xpetra::MapFactory<LO,GO,NO>::Build(this->K_->getMap()->lib(),-1,1,0,this->K_->getMap()->getComm());
            for (int i=1; i<MLgatheringSteps-1; i++) {
                MLnumMyRows = 0;
                numProcsGatheringStep = LO(numProcsGatheringStep/MLgatheringFactor);
                if (this->MpiComm_->getRank()%(this->MpiComm_->getSize()/numProcsGatheringStep) == 0 && this->MpiComm_->getRank()/(this->MpiComm_->getSize()/numProcsGatheringStep) < numProcsGatheringStep) {
                    if (this->MpiComm_->getRank()==0) {
                        MLnumMyRows = MLnumGlobalIndices - (MLnumGlobalIndices/numProcsGatheringStep)*(numProcsGatheringStep-1);
                    } else {
                        MLnumMyRows = MLnumGlobalIndices/numProcsGatheringStep;
                    }
                }
                MLGatheringMaps_[i] = Xpetra::MapFactory<LO,GO,NO>::Build(CoarseMap_->lib(),-1,MLnumMyRows,0,this->MpiComm_);

            }

            MLnumMyRows = 0;
            if (this->MpiComm_->getRank()%(this->MpiComm_->getSize()/NumProcsCoarseSolve_) == 0 && this->MpiComm_->getRank()/(this->MpiComm_->getSize()/NumProcsCoarseSolve_) < NumProcsCoarseSolve_) {
                if (this->MpiComm_->getRank()==0) {
                    MLnumMyRows = MLnumGlobalIndices - (MLnumGlobalIndices/NumProcsCoarseSolve_)*(NumProcsCoarseSolve_-1);
                } else {
                    MLnumMyRows = MLnumGlobalIndices/NumProcsCoarseSolve_;
                }
            }
            MLGatheringMaps_[MLgatheringSteps-1] = Xpetra::MapFactory<LO,GO,NO>::Build(CoarseMap_->lib(),-1,MLnumMyRows,0,this->MpiComm_);
            for (UN j=1; j<MLGatheringMaps_.size(); j++) {
                MLCoarseSolveExporters_[j] = Xpetra::ExportFactory<LO,GO,NO>::Build(MLGatheringMaps_[j-1],MLGatheringMaps_[j]);
            }
            int nSubs = this->MpiComm_->getSize();
            GOVec RowsCoarseSolve;
            if (OnCoarseSolveComm_) {
                int start = (nSubs*(CoarseSolveComm_->getRank()))/NumProcsCoarseSolve_;
                int end = (nSubs*(CoarseSolveComm_->getRank()+1))/NumProcsCoarseSolve_;
                RowsCoarseSolve.resize(end-start);
                for (int i = 0; i<end-start; i++) {
                    RowsCoarseSolve[i] = start+i;
                }
            }
            Teuchos::ArrayView< const GO > CList = MLGatheringMaps_[MLgatheringSteps-1]->getNodeElementList();
            //MLCoarseMap_ =  Xpetra::MapFactory<LO,GO,NO>::Build(CoarseMap_->lib(),-1,CList,0,CoarseSolveComm_);->Should work but does not WHY?!?!
            MLCoarseMap_ = Xpetra::MapFactory<LO,GO,NO>::Build(CoarseMap_->lib(),MLGatheringMaps_[MLgatheringSteps-1]->getGlobalNumElements(),0,CoarseSolveComm_);

            //#####################################################################
            // Build Repeated Map Zoltan2
            // build ElementNodeList_ to have adjacent entities to one subdomain
            this->buildElementNodeList();
            // Connectivity Graph on the CoarseSolveComm_
            this->buildCoarseGraph();
            //Build Repeatd Map on CoarseComm------------
            //Initialize Maps...
            ConstXMapPtr UniqueMap;
            XMapPtr UniqueMapAll;
            XMapPtr tmpRepMap;
            ConstXMapPtr ConstRepMap;
            GOVec uniEle;

            if (OnCoarseSolveComm_) {
                //Coarse DofsMaps so far only one Block will work
                ConstXMapPtrVecPtr2D CoarseDofsMaps(1);
                FROSch::BuildRepMapZoltan(SubdomainConnectGraph_,ElementNodeList_, DistributionList_,MLCoarseMap_->getComm(),CoarseSolveRepeatedMap_);
                ConstRepMap = CoarseSolveRepeatedMap_;
                ConstXMapPtrVecPtr NodesMapVector(1);
                //MapVector for next Level
                //So far only one Block is allowed ; needs to be adapetd fpr Block Ops
                ConstXMapPtrVecPtr RepMapVector(1);
                //Create DofMaps according to counting the interface entities
                ConstXMapPtrVecPtr DMap(CoarseDofsPerNode_);
                ConstXMapPtrVecPtr DMapRep(CoarseDofsPerNode_);

                tmpRepMap = this->BuildRepeatedMapCoarseLevel(ConstRepMap,CoarseDofsPerNode_,DMapRep,PartitionType_);
                RepMapCoarse_ = tmpRepMap;
                RepMapVector[0] = tmpRepMap;

                NodesMapVector[0] = ConstRepMap;
                //Pass Repeated Map Vector on to the next Level

                //Create uniqueMap following the repeatedMap
                //Create uniqueNodeMap so that dof belonging to one node are on the same process
                UniqueMap = FROSch::BuildUniqueMap<LO,GO,NO>(CoarseSolveRepeatedMap_);
                UniqueMapAll = this->BuildRepeatedMapCoarseLevel(UniqueMap,CoarseDofsPerNode_,DMap,PartitionType_);

                uniEle = UniqueMapAll->getNodeElementList();
                //Set DofOderingVec and DofsPerNodeVec to ParameterList for the next Level
                //Create Here DofsMaps for the next Level->DofOrdering will become redundant
                Teuchos::ArrayRCP<DofOrdering> dofOrderings(1);
                dofOrderings[0] = Custom; //special Ordering for Coarse Level
                Teuchos::ArrayRCP<UN> dofsPerNodeVector(1);
                dofsPerNodeVector[0] = CoarseDofsPerNode_;
                CoarseDofsMaps[0] = DMapRep;

                sublist(this->ParameterList_,"CoarseSolver")->set("Repeated Map Vector",RepMapVector);
                sublist(this->ParameterList_,"CoarseSolver")->set("Dofs Maps Vector",CoarseDofsMaps);
                sublist(this->ParameterList_,"CoarseSolver")->set("DofOrdering Vector",dofOrderings);
                sublist(this->ParameterList_,"CoarseSolver")->set("DofsPerNode Vector",dofsPerNodeVector);
                sublist(this->ParameterList_,"CoarseSolver")->set("Nodes Map Vector",NodesMapVector);
            }

            Teuchos::RCP<Xpetra::Map<LO,GO,NO> > tmpMap = Xpetra::MapFactory<LO,GO,NO>::Build(CoarseMap_->lib(),-1,uniEle,0,this->MpiComm_);

            for (int i=0; i<gatheringSteps-1; i++) {
                numMyRows = 0;
                numProcsGatheringStep = LO(numProcsGatheringStep/gatheringFactor);
                if (this->MpiComm_->getRank()%(this->MpiComm_->getSize()/numProcsGatheringStep) == 0 && this->MpiComm_->getRank()/(this->MpiComm_->getSize()/numProcsGatheringStep) < numProcsGatheringStep) {
                    if (this->MpiComm_->getRank()==0) {
                        numMyRows = numGlobalIndices - (numGlobalIndices/numProcsGatheringStep)*(numProcsGatheringStep-1);
                    } else {
                        numMyRows = numGlobalIndices/numProcsGatheringStep;
                    }
                }
                GatheringMaps_[i] = Xpetra::MapFactory<LO,GO,NO>::Build(CoarseMap_->lib(),-1,numMyRows,0,this->MpiComm_);
            }
            GatheringMaps_[gatheringSteps-1] = tmpMap;
            CoarseSolveMap_ = Xpetra::MapFactory<LO,GO,NO>::Build(CoarseMap_->lib(),-1,tmpMap->getNodeElementList(),0,CoarseSolveComm_);
        } else if (!DistributionList_->get("Type","linear").compare("Zoltan2")) {
#ifdef HAVE_SHYLU_DDFROSCH_ZOLTAN2
            GatheringMaps_.resize(1);
            CoarseSolveExporters_.resize(1);
#ifdef FROSCH_COARSEOPERATOR_EXPORT_AND_IMPORT
            CoarseSolveImporters_.resize(1);
#endif
#else
            ThrowErrorMissingPackage("FROSch::CoarseOperator","Zoltan2");
#endif
        } else {
            FROSCH_ASSERT(false,"FROSch::CoarseOperator : ERROR: Distribution type unknown.");
        }

        return 0;
    }

}

#endif
