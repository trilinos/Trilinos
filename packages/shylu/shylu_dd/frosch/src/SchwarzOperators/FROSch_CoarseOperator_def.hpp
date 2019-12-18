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

    using namespace Teuchos;
    using namespace Xpetra;

    template<class SC,class LO,class GO,class NO>
    CoarseOperator<SC,LO,GO,NO>::CoarseOperator(ConstXMatrixPtr k,
                                                ParameterListPtr parameterList) :
    SchwarzOperator<SC,LO,GO,NO> (k,parameterList),
    CoarseSpace_ (new CoarseSpace<SC,LO,GO,NO>(this->MpiComm_,this->SerialComm_)),
    DistributionList_ (sublist(parameterList,"Distribution"))
    {
        FROSCH_TIMER_START_LEVELID(coarseOperatorTime,"CoarseOperator::CoarseOperator");
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
        //if (this->Verbose_) std::cout << "FROSch::CoarseOperator : WARNING: Some of the operations could probably be moved from initialize() to Compute().\n";

        bool reuseCoarseBasis = this->ParameterList_->get("Reuse: Coarse Basis",true);
        bool reuseCoarseMatrix = this->ParameterList_->get("Reuse: Coarse Matrix",false);
        if (!this->IsComputed_) {
            reuseCoarseBasis = false;
            reuseCoarseMatrix = false;
        }

        if (!reuseCoarseBasis) {
            if (this->IsComputed_ && this->Verbose_) std::cout << "FROSch::CoarseOperator : Recomputing the Coarse Basis" << std::endl;
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
        if ( this->ParameterList_->get("Set Phi to PList", false ) ){
            if (this->Verbose_)
                std::cout << "\t### Setting Phi (RCP<Xpetra::Matrix>) to ParameterList.\n";
            
            this->ParameterList_->set("Phi Pointer", Phi_);
        }
        if (!reuseCoarseMatrix) {
            if (this->IsComputed_ && this->Verbose_) std::cout << "FROSch::CoarseOperator : Recomputing the Coarse Matrix" << std::endl;
            this->setUpCoarseOperator();
        }
        this->IsComputed_ = true;
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
        FROSCH_TIMER_START_LEVELID(applyPhiTTime,"CoarseOperator::applyPhiT");
        // AH 08/22/2019 TODO: We cannot ger rid of the Build() calls because of "XCoarse_ = XCoarseSolveTmp_;". This is basically caused by the whole Gathering Map strategy. As soon as we have replaced this, we can get rid of the Build() calls
        XCoarse_ = MultiVectorFactory<SC,LO,GO,NO>::Build(CoarseSpace_->getBasisMapUnique(),x.getNumVectors()); // AH 08/22/2019 TODO: Can we get rid of this? If possible, we should remove the whole GatheringMaps idea and replace it by some smart all-to-all MPI communication
        {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
            FROSCH_TIMER_START_LEVELID(applyTime,"apply");
#endif
            Phi_->apply(x,*XCoarse_,TRANS);
        }
        for (UN j=0; j<GatheringMaps_.size(); j++) {
            XCoarseSolveTmp_ = MultiVectorFactory<SC,LO,GO,NO>::Build(GatheringMaps_[j],x.getNumVectors()); // AH 08/22/2019 TODO: Can we get rid of this? If possible, we should remove the whole GatheringMaps idea and replace it by some smart all-to-all MPI communication
            {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                FROSCH_TIMER_START_LEVELID(applyTime,"doExport");
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
        FROSCH_TIMER_START_LEVELID(applyCoarseSolveTime,"CoarseOperator::applyCoarseSolve");
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
        FROSCH_TIMER_START_LEVELID(applyPhiTime,"CoarseOperator::applyPhi");
        // AH 08/22/2019 TODO: We have the same issue here as in applyPhiT()
        YCoarseSolveTmp_ = MultiVectorFactory<SC,LO,GO,NO>::Build(x.getMap(),x.getNumVectors());
        *YCoarseSolveTmp_ = x;
        for (int j=GatheringMaps_.size()-1; j>0; j--) {
            YCoarse_ = MultiVectorFactory<SC,LO,GO,NO>::Build(GatheringMaps_[j-1],x.getNumVectors());
            {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                FROSCH_TIMER_START_LEVELID(applyTime,"doImport");
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
            FROSCH_TIMER_START_LEVELID(applyTime,"doImport");
#endif
#ifdef FROSCH_COARSEOPERATOR_EXPORT_AND_IMPORT
            YCoarse_->doImport(*YCoarseSolveTmp_,*CoarseSolveImporters_[0],INSERT);
#else
            YCoarse_->doImport(*YCoarseSolveTmp_,*CoarseSolveExporters_[0],INSERT);
#endif
        }
        {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
            FROSCH_TIMER_START_LEVELID(applyTime,"apply");
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
        FROSCH_TIMER_START_LEVELID(setUpCoarseOperatorTime,"CoarseOperator::setUpCoarseOperator");
        if (!Phi_.is_null()) {
            // Build CoarseMatrix_
            XMatrixPtr k0 = buildCoarseMatrix();

            //------------------------------------------------------------------------------------------------------------------------
            // Communicate coarse matrix
            if (!DistributionList_->get("Type","linear").compare("linear")) {
                XMatrixPtr tmpCoarseMatrix = MatrixFactory<SC,LO,GO,NO>::Build(GatheringMaps_[0]);
                {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                    FROSCH_TIMER_START_LEVELID(coarseMatrixExportTime,"Export Coarse Matrix");
#endif
                    tmpCoarseMatrix->doExport(*k0,*CoarseSolveExporters_[0],INSERT);
                }

                for (UN j=1; j<GatheringMaps_.size(); j++) {
                    tmpCoarseMatrix->fillComplete();
                    k0 = tmpCoarseMatrix;
                    tmpCoarseMatrix = MatrixFactory<SC,LO,GO,NO>::Build(GatheringMaps_[j]);
                    {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                        FROSCH_TIMER_START_LEVELID(coarseMatrixExportTime,"Export Coarse Matrix");
#endif
                        tmpCoarseMatrix->doExport(*k0,*CoarseSolveExporters_[j],INSERT);
                    }
                }
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

            //------------------------------------------------------------------------------------------------------------------------
            // Matrix to the new communicator
            if (OnCoarseSolveComm_) {
                LO numRows = k0->getNodeNumRows();
                ArrayRCP<size_t> elemsPerRow(numRows);
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
                        }
                    }
                    CoarseMatrix_->fillComplete(CoarseSolveMap_,CoarseSolveMap_); //RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(std::cout)); CoarseMatrix_->describe(*fancy,VERB_EXTREME);
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
                        }
                    }
                    CoarseMatrix_->fillComplete(CoarseSolveMap_,CoarseSolveMap_); //RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(std::cout)); CoarseMatrix_->describe(*fancy,VERB_EXTREME);
                }
                
                bool reuseCoarseMatrixSymbolicFactorization = this->ParameterList_->get("Reuse: Coarse Matrix Symbolic Factorization",true);
                if (!this->IsComputed_) {
                    reuseCoarseMatrixSymbolicFactorization = false;
                }
                if (!reuseCoarseMatrixSymbolicFactorization) {
                    if (this->IsComputed_ && this->Verbose_) std::cout << "FROSch::CoarseOperator : Recomputing the Symbolic Factorization of the coarse matrix" << std::endl;
                    CoarseSolver_.reset(new SubdomainSolver<SC,LO,GO,NO>(CoarseMatrix_,sublist(this->ParameterList_,"CoarseSolver")));
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
        FROSCH_TIMER_START_LEVELID(buildCoarseMatrixTime,"CoarseOperator::buildCoarseMatrix");
        XMatrixPtr k0;
        if (this->ParameterList_->get("Use Triple MatrixMultiply",false)) {
            k0 = MatrixFactory<SC,LO,GO,NO>::Build(CoarseSpace_->getBasisMapUnique(),as<LO>(0));
            TripleMatrixMultiply<SC,LO,GO,NO>::MultiplyRAP(*Phi_,true,*this->K_,false,*Phi_,false,*k0);
        } else {
            RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(std::cout));
            XMatrixPtr tmp = MatrixMatrix<SC,LO,GO,NO>::Multiply(*this->K_,false,*Phi_,false,*fancy);
            k0 = MatrixMatrix<SC,LO,GO,NO>::Multiply(*Phi_,true,*tmp,false,*fancy);
        }
        return k0;
    }

    template<class SC,class LO,class GO,class NO>
    int CoarseOperator<SC,LO,GO,NO>::buildCoarseSolveMap(ConstXMapPtr coarseMapUnique)
    {
        FROSCH_TIMER_START_LEVELID(buildCoarseSolveMapTime,"CoarseOperator::buildCoarseSolveMap");
        NumProcsCoarseSolve_ = DistributionList_->get("NumProcs",1);
        double factor = DistributionList_->get("Factor",0.0);
        
        switch (NumProcsCoarseSolve_) {
            case -1:
                FROSCH_ASSERT(false,"We do not know the size of the matrix yet. Therefore, we cannot use the formula NumProcsCoarseSolve_ = int(0.5*(1+std::max(k0->getGlobalNumRows()/10000,k0->getGlobalNumEntries()/100000)));");
                //NumProcsCoarseSolve_ = int(0.5*(1+std::max(k0->getGlobalNumRows()/10000,k0->getGlobalNumEntries()/100000)));
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
                //if (this->Verbose_) std::cout << i << " " << numProcsGatheringStep << " " << numGlobalIndices << std::endl;
                if (this->MpiComm_->getRank()%(this->MpiComm_->getSize()/numProcsGatheringStep) == 0 && this->MpiComm_->getRank()/(this->MpiComm_->getSize()/numProcsGatheringStep) < numProcsGatheringStep) {
                    if (this->MpiComm_->getRank()==0) {
                        numMyRows = numGlobalIndices - (numGlobalIndices/numProcsGatheringStep)*(numProcsGatheringStep-1);
                    } else {
                        numMyRows = numGlobalIndices/numProcsGatheringStep;
                    }
                }
                {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                    FROSCH_TIMER_START_LEVELID(gatheringMapsTime,"Gathering Maps");
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
                FROSCH_TIMER_START_LEVELID(gatheringMapsTime,"Gathering Maps");
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
                FROSCH_TIMER_START_LEVELID(commSplitTime,"Coarse Communicator Split");
#endif
                CoarseSolveComm_ = this->MpiComm_->split(!OnCoarseSolveComm_,this->MpiComm_->getRank());
            }
            {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                FROSCH_TIMER_START_LEVELID(coarseCommMapTime,"Coarse Communicator Map");
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
                FROSCH_TIMER_START_LEVELID(coarseSolveExportersTime,"Build Exporters");
#endif
                CoarseSolveExporters_[0] = ExportFactory<LO,GO,NO>::Build(coarseMapUnique,GatheringMaps_[0]);
                CoarseSolveExporters_[0]->setDistributorParameters(gatheringCommunicationList); // Set the parameter list for the communication of the exporter
            }
#ifdef FROSCH_COARSEOPERATOR_EXPORT_AND_IMPORT
            {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                FROSCH_TIMER_START_LEVELID(coarseSolveImportersTime,"Build Importers");
#endif
                CoarseSolveImporters_[0] = ImportFactory<LO,GO,NO>::Build(GatheringMaps_[0],coarseMapUnique);
                CoarseSolveImporters_[0]->setDistributorParameters(gatheringCommunicationList); // Set the parameter list for the communication of the exporter
            }
#endif
            
            for (UN j=1; j<GatheringMaps_.size(); j++) {
                {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                    FROSCH_TIMER_START_LEVELID(coarseSolveExportersTime,"Build Exporters");
#endif
                    CoarseSolveExporters_[j] = ExportFactory<LO,GO,NO>::Build(GatheringMaps_[j-1],GatheringMaps_[j]);
                    CoarseSolveExporters_[j]->setDistributorParameters(gatheringCommunicationList); // Set the parameter list for the communication of the exporter
                }
#ifdef FROSCH_COARSEOPERATOR_EXPORT_AND_IMPORT
                {
#ifdef FROSCH_COARSEOPERATOR_DETAIL_TIMERS
                    FROSCH_TIMER_START_LEVELID(coarseSolveImportersTime,"Build Importers");
#endif
                    CoarseSolveImporters_[j] = ImportFactory<LO,GO,NO>::Build(GatheringMaps_[j],GatheringMaps_[j-1]);
                    CoarseSolveImporters_[j]->setDistributorParameters(gatheringCommunicationList); // Set the parameter list for the communication of the exporter
                }
#endif
            }
        } else if(!DistributionList_->get("Type","linear").compare("Zoltan2")) {
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

        if (this->Verbose_) {
            std::cout << "\n\
    ------------------------------------------------------------------------------\n\
     Coarse problem statistics\n\
    ------------------------------------------------------------------------------\n\
      Dimension of the coarse problem             --- " << coarseMapUnique->getMaxAllGlobalIndex()+1 << "\n\
      Number of processes                         --- " << NumProcsCoarseSolve_ << "\n\
    ------------------------------------------------------------------------------\n";
        }

        return 0;
    }
    
}

#endif
