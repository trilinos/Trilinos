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
    CoarseSolveComm_ (),
    OnCoarseSolveComm_ (false),
    NumProcsCoarseSolve_ (0),
    CoarseSpace_ (new CoarseSpace<SC,LO,GO,NO>()),
    Phi_ (),
    CoarseMatrix_ (),
    XTmp_ (),
    XCoarse_ (),
    XCoarseSolve_ (),
    XCoarseSolveTmp_ (),
    YTmp_ (),
    YCoarse_ (),
    YCoarseSolve_ (),
    YCoarseSolveTmp_ (),
    GatheringMaps_ (0),
    CoarseSolveMap_ (),
    CoarseSolveRepeatedMap_ (),
    CoarseSolver_ (),
    DistributionList_ (sublist(parameterList,"Distribution")),
    CoarseSolveExporters_ (0)
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
        FROSCH_ASSERT(this->IsInitialized_,"ERROR: CoarseOperator has to be initialized before calling compute()");
        // This is not optimal yet... Some work could be moved to Initialize
        if (this->Verbose_) std::cout << "FROSch::CoarseOperator : WARNING: Some of the operations could probably be moved from initialize() to Compute().\n";

        if (!this->ParameterList_->get("Recycling","none").compare("basis") && this->IsComputed_) {
            this->setUpCoarseOperator();

            this->IsComputed_ = true;
        } else if(!this->ParameterList_->get("Recycling","none").compare("all") && this->IsComputed_) {
            // Maybe use some advanced settings in the future
        } else {
            clearCoarseSpace(); // AH 12/11/2018: If we do not clear the coarse space, we will always append just append the coarse space
            XMapPtr subdomainMap = this->computeCoarseSpace(CoarseSpace_); // AH 12/11/2018: This map could be overlapping, repeated, or unique. This depends on the specific coarse operator
            CoarseSpace_->assembleCoarseSpace();
            CoarseSpace_->buildGlobalBasisMatrix(this->K_->getRangeMap(),subdomainMap,this->ParameterList_->get("Threshold Phi",1.e-8));
            Phi_ = CoarseSpace_->getGlobalBasisMatrix();
            this->setUpCoarseOperator();

            this->IsComputed_ = true;
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
        if (this->IsComputed_) {
            if (XTmp_.is_null()) XTmp_ = MultiVectorFactory<SC,LO,GO,NO>::Build(x.getMap(),x.getNumVectors());
            *XTmp_ = x;
            if (!usePreconditionerOnly && mode == NO_TRANS) {
                this->K_->apply(x,*XTmp_,mode,ScalarTraits<SC>::one(),ScalarTraits<SC>::zero());
            }
            if (XCoarseSolve_.is_null()) XCoarseSolve_ = MultiVectorFactory<SC,LO,GO,NO>::Build(GatheringMaps_[GatheringMaps_.size()-1],x.getNumVectors());
            if (YCoarseSolve_.is_null()) YCoarseSolve_ = MultiVectorFactory<SC,LO,GO,NO>::Build(GatheringMaps_[GatheringMaps_.size()-1],y.getNumVectors());
            applyPhiT(*XTmp_,*XCoarseSolve_);
            applyCoarseSolve(*XCoarseSolve_,*YCoarseSolve_,mode);
            applyPhi(*YCoarseSolve_,*XTmp_);
            if (!usePreconditionerOnly && mode != NO_TRANS) {
                this->K_->apply(*XTmp_,*XTmp_,mode,ScalarTraits<SC>::one(),ScalarTraits<SC>::zero());
            }
            y.update(alpha,*XTmp_,beta);
        } else {
            if (i==1) {
                if (this->Verbose_) std::cout << "WARNING: CoarseOperator has not been computed yet => It will just act as the identity...\n";
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
        XCoarse_ = MultiVectorFactory<SC,LO,GO,NO>::Build(CoarseSpace_->getBasisMap(),x.getNumVectors()); // AH 08/22/2019 TODO: Can we get rid of this? If possible, we should remove the whole GatheringMaps idea and replace it by some smart all-to-all MPI communication
        Phi_->apply(x,*XCoarse_,TRANS);
        for (UN j=0; j<GatheringMaps_.size(); j++) {
            XCoarseSolveTmp_ = MultiVectorFactory<SC,LO,GO,NO>::Build(GatheringMaps_[j],x.getNumVectors()); // AH 08/22/2019 TODO: Can we get rid of this? If possible, we should remove the whole GatheringMaps idea and replace it by some smart all-to-all MPI communication
            XCoarseSolveTmp_->doExport(*XCoarse_,*CoarseSolveExporters_[j],ADD);
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
            CoarseSolver_->apply(x,*YTmp_,mode);
        } else {
            if (YTmp_.is_null()) YTmp_ = MultiVectorFactory<SC,LO,GO,NO>::Build(CoarseSolveMap_,x.getNumVectors());
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
            YCoarse_->doImport(*YCoarseSolveTmp_,*CoarseSolveExporters_[j],INSERT);
            YCoarseSolveTmp_ = YCoarse_;
        }
        YCoarse_ = MultiVectorFactory<SC,LO,GO,NO>::Build(CoarseSpace_->getBasisMap(),x.getNumVectors());
        YCoarse_->doImport(*YCoarseSolveTmp_,*CoarseSolveExporters_[0],INSERT);
        Phi_->apply(*YCoarse_,y,NO_TRANS);
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
        // Build CoarseMatrix_
        XMatrixPtr k0 = buildCoarseMatrix();

        // Build Map for the coarse solver
        buildCoarseSolveMap(k0);

        //------------------------------------------------------------------------------------------------------------------------
        // Communicate coarse matrix
        if (!DistributionList_->get("Type","linear").compare("linear")) {
            CoarseSolveExporters_[0] = ExportFactory<LO,GO,NO>::Build(CoarseSpace_->getBasisMap(),GatheringMaps_[0]);

            XMatrixPtr tmpCoarseMatrix = MatrixFactory<SC,LO,GO,NO>::Build(GatheringMaps_[0],k0->getGlobalMaxNumRowEntries());

            tmpCoarseMatrix->doExport(*k0,*CoarseSolveExporters_[0],INSERT);

            for (UN j=1; j<GatheringMaps_.size(); j++) {
                tmpCoarseMatrix->fillComplete();
                k0 = tmpCoarseMatrix;
                CoarseSolveExporters_[j] = ExportFactory<LO,GO,NO>::Build(GatheringMaps_[j-1],GatheringMaps_[j]);
                tmpCoarseMatrix = MatrixFactory<SC,LO,GO,NO>::Build(GatheringMaps_[j],k0->getGlobalMaxNumRowEntries());

                tmpCoarseMatrix->doExport(*k0,*CoarseSolveExporters_[j],INSERT);
            }
            //------------------------------------------------------------------------------------------------------------------------
            // Matrix to the new communicator
            if (OnCoarseSolveComm_) {
                LO numRows = tmpCoarseMatrix->getNodeNumRows();
                ArrayRCP<size_t> elemsPerRow(numRows);
                ConstGOVecView indices;
                ConstSCVecView values;
                for (LO i = 0; i < numRows; i++) {
                  GO globalRow = CoarseSolveMap_->getGlobalElement(i);
                  size_t numEntries = tmpCoarseMatrix->getNumEntriesInGlobalRow(globalRow);
                  if(numEntries == 0)
                  {
                    //Always add the diagonal for empty rows
                    numEntries = 1;
                  }
                  elemsPerRow[i] = numEntries;
                }
                CoarseMatrix_ = MatrixFactory<SC,LO,GO,NO>::Build(CoarseSolveMap_, elemsPerRow, StaticProfile);
                for (LO i = 0; i < numRows; i++) {
                    GO globalRow = CoarseSolveMap_->getGlobalElement(i);
                    tmpCoarseMatrix->getGlobalRowView(globalRow,indices,values);
                    if (indices.size()>0) {
                        CoarseMatrix_->insertGlobalValues(globalRow,indices,values);
                    } else { // Add diagonal unit for zero rows // Todo: Do you we need to sort the coarse matrix "NodeWise"?
                        GOVec indices(1,globalRow);
                        SCVec values(1,ScalarTraits<SC>::one());
                        CoarseMatrix_->insertGlobalValues(globalRow,indices(),values());
                    }

                }

                CoarseMatrix_->fillComplete(CoarseSolveMap_,CoarseSolveMap_); //RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(std::cout)); CoarseMatrix_->describe(*fancy,VERB_EXTREME);
                CoarseSolver_.reset(new SubdomainSolver<SC,LO,GO,NO>(CoarseMatrix_,sublist(this->ParameterList_,"CoarseSolver")));
                CoarseSolver_->initialize();

                CoarseSolver_->compute();

            }
#ifdef HAVE_SHYLU_DDFROSCH_ZOLTAN2
        } else if (!DistributionList_->get("Type","linear").compare("Zoltan2")) {
            //------------------------------------------------------------------------------------------------------------------------
            //coarse matrix already communicated with Zoltan2. Communicate to CoarseSolveComm.
            //------------------------------------------------------------------------------------------------------------------------

            // Matrix to the new communicator
            if (OnCoarseSolveComm_) {
                CoarseMatrix_ = MatrixFactory<SC,LO,GO,NO>::Build(CoarseSolveMap_,k0->getGlobalMaxNumRowEntries());
                ConstLOVecView indices;
                ConstSCVecView values;
                for (UN i=0; i<k0->getNodeNumRows(); i++) {
                    // different sorted maps: CoarseSolveMap_ and k0
                    LO locRow = k0->getRowMap()->getLocalElement(CoarseSolveMap_->getGlobalElement(i));
                    k0->getLocalRowView(locRow,indices,values);
                    if (indices.size()>0) {
                        GOVec indicesGlob(indices.size());
                        for (UN j=0; j<indices.size(); j++) {
                            indicesGlob[j] = k0->getColMap()->getGlobalElement(indices[j]);
                        }
                        CoarseMatrix_->insertGlobalValues(CoarseSolveMap_->getGlobalElement(i),indicesGlob(),values);
                    } else { // Add diagonal unit for zero rows // Todo: Do you we need to sort the coarse matrix "NodeWise"?
                        GOVec indices(1,CoarseSolveMap_->getGlobalElement(i));
                        SCVec values(1,ScalarTraits<SC>::one());
                        CoarseMatrix_->insertGlobalValues(CoarseSolveMap_->getGlobalElement(i),indices(),values());
                    }

                }

                CoarseMatrix_->fillComplete(CoarseSolveMap_,CoarseSolveMap_); //RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(std::cout)); CoarseMatrix_->describe(*fancy,VERB_EXTREME);

                if (!this->ParameterList_->sublist("CoarseSolver").get("SolverType","Amesos").compare("MueLu")) {
                    CoarseSolver_.reset(new SubdomainSolver<SC,LO,GO,NO>(CoarseMatrix_,sublist(this->ParameterList_,"CoarseSolver")));
                }
                else{
                    CoarseSolver_.reset(new SubdomainSolver<SC,LO,GO,NO>(CoarseMatrix_,sublist(this->ParameterList_,"CoarseSolver")));
                }

                CoarseSolver_->initialize();

                CoarseSolver_->compute();

            }
            //------------------------------------------------------------------------------------------------------------------------
#endif
        } else {
            FROSCH_ASSERT(false,"Distribution Type unknown!");
        }



        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    typename CoarseOperator<SC,LO,GO,NO>::XMatrixPtr CoarseOperator<SC,LO,GO,NO>::buildCoarseMatrix()
    {
        FROSCH_TIMER_START_LEVELID(buildCoarseMatrixTime,"CoarseOperator::buildCoarseMatrix");
        XMatrixPtr k0 = MatrixFactory<SC,LO,GO,NO>::Build(CoarseSpace_->getBasisMap(),CoarseSpace_->getBasisMap()->getNodeNumElements());

        if (this->ParameterList_->get("Use Triple MatrixMultiply",false)) {
            TripleMatrixMultiply<SC,LO,GO,NO>::MultiplyRAP(*Phi_,true,*this->K_,false,*Phi_,false,*k0);
        }
        else{
            XMatrixPtr tmp = MatrixFactory<SC,LO,GO,NO>::Build(this->K_->getRowMap(),50);
            MatrixMatrix<SC,LO,GO,NO>::Multiply(*this->K_,false,*Phi_,false,*tmp);
            MatrixMatrix<SC,LO,GO,NO>::Multiply(*Phi_,true,*tmp,false,*k0);
        }
        return k0;
    }

    template<class SC,class LO,class GO,class NO>
    int CoarseOperator<SC,LO,GO,NO>::buildCoarseSolveMap(XMatrixPtr &k0)
    {
        FROSCH_TIMER_START_LEVELID(buildCoarseSolveMapTime,"CoarseOperator::buildCoarseSolveMap");
        NumProcsCoarseSolve_ = DistributionList_->get("NumProcs",0);
        double fac = DistributionList_->get("Factor",1.0);

        // Redistribute Matrix
        if (NumProcsCoarseSolve_==0) {
            NumProcsCoarseSolve_ = this->MpiComm_->getSize();//Phi->DomainMap().Comm().getSize();
        } else if (NumProcsCoarseSolve_==1) {
            NumProcsCoarseSolve_ = 1;
        } else if (NumProcsCoarseSolve_==-1) {
            NumProcsCoarseSolve_ = int(1+std::max(k0->getGlobalNumRows()/10000,k0->getGlobalNumEntries()/100000));
        } else if (NumProcsCoarseSolve_>1) {

        } else if (NumProcsCoarseSolve_<-1) {
            NumProcsCoarseSolve_ = round(pow(1.0*this->MpiComm_->getSize(), 1./(-NumProcsCoarseSolve_)));
        } else {
            FROSCH_ASSERT(false,"This should never happen...");
        }

        NumProcsCoarseSolve_ = (LO)  NumProcsCoarseSolve_ * fac;
        if (NumProcsCoarseSolve_<1) {
            NumProcsCoarseSolve_ = 1;
        }

        if (NumProcsCoarseSolve_ >= this->MpiComm_->getSize() && DistributionList_->get("Type","linear").compare("Zoltan2")) {
            GatheringMaps_.resize(1);
            CoarseSolveExporters_.resize(1);
            GatheringMaps_[0] = BuildUniqueMap<LO,GO,NO>(Phi_->getColMap()); // DO WE NEED THIS IN ANY CASE???
            return 0;
        }
        //cout << DistributionList_->get("Type","linear") << std::endl;
        if (!DistributionList_->get("Type","linear").compare("linear")) {

            int gatheringSteps = DistributionList_->get("GatheringSteps",1);
            GatheringMaps_.resize(gatheringSteps);
            CoarseSolveExporters_.resize(gatheringSteps);

            LO numProcsGatheringStep = this->MpiComm_->getSize();
            GO numGlobalIndices = CoarseSpace_->getBasisMap()->getMaxAllGlobalIndex()+1;
            GO numMyRows;
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
                GatheringMaps_[i] = MapFactory<LO,GO,NO>::Build(CoarseSpace_->getBasisMap()->lib(),-1,numMyRows,0,this->MpiComm_);
            }

            numMyRows = 0;
            if (this->MpiComm_->getRank()%(this->MpiComm_->getSize()/NumProcsCoarseSolve_) == 0 && this->MpiComm_->getRank()/(this->MpiComm_->getSize()/NumProcsCoarseSolve_) < NumProcsCoarseSolve_) {
                if (this->MpiComm_->getRank()==0) {
                    numMyRows = numGlobalIndices - (numGlobalIndices/NumProcsCoarseSolve_)*(NumProcsCoarseSolve_-1);
                } else {
                    numMyRows = numGlobalIndices/NumProcsCoarseSolve_;
                }
            }
            GatheringMaps_[gatheringSteps-1] = MapFactory<LO,GO,NO>::Build(CoarseSpace_->getBasisMap()->lib(),-1,numMyRows,0,this->MpiComm_);
            //cout << *GatheringMaps_->at(gatheringSteps-1);

            //------------------------------------------------------------------------------------------------------------------------
            // Use a separate Communicator for the coarse problem
            ConstXMapPtr tmpCoarseMap = GatheringMaps_[GatheringMaps_.size()-1];

            if (tmpCoarseMap->getNodeNumElements()>0) {
                OnCoarseSolveComm_=true;
            }
            CoarseSolveComm_ = this->MpiComm_->split(!OnCoarseSolveComm_,this->MpiComm_->getRank());
            CoarseSolveMap_ = MapFactory<LO,GO,NO>::Build(CoarseSpace_->getBasisMap()->lib(),-1,tmpCoarseMap->getNodeElementList(),0,CoarseSolveComm_);

#ifdef HAVE_SHYLU_DDFROSCH_ZOLTAN2
        } else if(!DistributionList_->get("Type","linear").compare("Zoltan2")){

            RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(std::cout));

            GatheringMaps_.resize(1);
            CoarseSolveExporters_.resize(1);

            GatheringMaps_[0] = rcp_const_cast<XMap> (BuildUniqueMap(k0->getRowMap()));
            //
            CoarseSolveExporters_[0] = ExportFactory<LO,GO,NO>::Build(CoarseSpace_->getBasisMap(),GatheringMaps_[0]);

            XMatrixPtr k0Unique = MatrixFactory<SC,LO,GO,NO>::Build(GatheringMaps_[0],k0->getGlobalMaxNumRowEntries());

            k0Unique->doExport(*k0,*CoarseSolveExporters_[0],INSERT);
            k0Unique->fillComplete(GatheringMaps_[0],GatheringMaps_[0]);
            if (NumProcsCoarseSolve_<this->MpiComm_->getSize()) {
                ParameterListPtr tmpList = sublist(DistributionList_,"Zoltan2 Parameter");
                tmpList->set("num_global_parts", NumProcsCoarseSolve_);
                FROSch::RepartionMatrixZoltan2(k0Unique,tmpList);
            }

            k0 = k0Unique;

            GatheringMaps_[0] = k0->getRowMap();
            CoarseSolveExporters_[0] = ExportFactory<LO,GO,NO>::Build(CoarseSpace_->getBasisMap(),GatheringMaps_[0]);

            ConstXMapPtr tmpCoarseMap = GatheringMaps_[0];

            if (tmpCoarseMap->getNodeNumElements()>0) {
                OnCoarseSolveComm_=true;
            }

            GOVec elementList(tmpCoarseMap->getNodeElementList());
            CoarseSolveComm_ = this->MpiComm_->split(!OnCoarseSolveComm_,this->MpiComm_->getRank());
            CoarseSolveMap_ = MapFactory<LO,GO,NO>::Build(CoarseSpace_->getBasisMap()->lib(),-1,elementList,0,CoarseSolveComm_);
#endif
        } else {
            FROSCH_ASSERT(false,"Distribution type not defined...");
        }

        if (this->Verbose_) {
            std::cout << "\n\
    ------------------------------------------------------------------------------\n\
     Coarse problem statistics\n\
    ------------------------------------------------------------------------------\n\
      dimension of the coarse problem             --- " << CoarseSpace_->getBasisMap()->getMaxAllGlobalIndex()+1 << "\n\
      number of processes                         --- " << NumProcsCoarseSolve_ << "\n\
    ------------------------------------------------------------------------------\n";
        }

        return 0;
    }

}

#endif
