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
    
    template<class SC,class LO,class GO,class NO>
    CoarseOperator<SC,LO,GO,NO>::CoarseOperator(CrsMatrixPtr k,
                                                ParameterListPtr parameterList) :
    SchwarzOperator<SC,LO,GO,NO> (k,parameterList),
    CoarseSolveComm_ (),
    OnCoarseSolveComm_ (false),
    NumProcsCoarseSolve_ (0),
    Phi_ (),
    CoarseMatrix_ (),
    CoarseMap_ (),
    GatheringMaps_ (0),
    CoarseSolveMap_ (),
    CoarseSolveRepeatedMap_ (),
    CoarseSolver_ (),
    DistributionList_ (sublist(parameterList,"Distribution")),
    CoarseSolveExporters_ (0)
    {
        
    }
    
    template<class SC,class LO,class GO,class NO>
    CoarseOperator<SC,LO,GO,NO>::~CoarseOperator()
    {
        CoarseSolver_.reset();
    }
    
    template<class SC,class LO,class GO,class NO>
    void CoarseOperator<SC,LO,GO,NO>::apply(const MultiVector &x,
                                            MultiVector &y,
                                            bool usePreconditionerOnly,
                                            Teuchos::ETransp mode,
                                            SC alpha,
                                            SC beta) const
    {
        static int i = 0;
        
        if (this->IsComputed_) {
            MultiVectorPtr xTmp = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(x.getMap(),x.getNumVectors());
            *xTmp = x;
            
            if (!usePreconditionerOnly && mode == Teuchos::NO_TRANS) {
                this->K_->apply(x,*xTmp,mode,1.0,0.0);
            }
            
            MultiVectorPtr xCoarseSolve = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(GatheringMaps_[GatheringMaps_.size()-1],x.getNumVectors());
            MultiVectorPtr yCoarseSolve = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(GatheringMaps_[GatheringMaps_.size()-1],y.getNumVectors());
            
            applyPhiT(*xTmp,*xCoarseSolve); 
            applyCoarseSolve(*xCoarseSolve,*yCoarseSolve,mode);
            applyPhi(*yCoarseSolve,*xTmp);
            
            if (!usePreconditionerOnly && mode != Teuchos::NO_TRANS) {
                this->K_->apply(*xTmp,*xTmp,mode,1.0,0.0);
            }
            y.update(alpha,*xTmp,beta);
        } else {
            if (i==1) {
                if (this->Verbose_) std::cout << "WARNING: CoarseOperator has not been computed yet => It will just act as the identity...\n";
                i++;
            }
            y.update(1.0,x,0.0);
        }
    }
    
    template<class SC,class LO,class GO,class NO>
    void CoarseOperator<SC,LO,GO,NO>::applyPhiT(MultiVector& x,
                                                MultiVector& y) const
    {
        MultiVectorPtr xCoarse = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(CoarseMap_,x.getNumVectors());
        
        Phi_->apply(x,*xCoarse,Teuchos::TRANS);
        
        MultiVectorPtr xCoarseSolveTmp;
        for (UN j=0; j<GatheringMaps_.size(); j++) {
            xCoarseSolveTmp = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(GatheringMaps_[j],x.getNumVectors());
            xCoarseSolveTmp->doExport(*xCoarse,*CoarseSolveExporters_[j],Xpetra::ADD);
            xCoarse = xCoarseSolveTmp;
        }
        y = *xCoarseSolveTmp;
    }
    
    template<class SC,class LO,class GO,class NO>
    void CoarseOperator<SC,LO,GO,NO>::applyCoarseSolve(MultiVector& x,
                                                       MultiVector& y,
                                                       Teuchos::ETransp mode) const
    {
        MultiVectorPtr yTmp;
        if (OnCoarseSolveComm_) {
            x.replaceMap(CoarseSolveMap_);
            yTmp = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(CoarseSolveMap_,x.getNumVectors());
            CoarseSolver_->apply(x,*yTmp,mode);
        } else {
            yTmp = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(CoarseSolveMap_,x.getNumVectors());
        }
        yTmp->replaceMap(GatheringMaps_[GatheringMaps_.size()-1]);
        y = *yTmp;
    }
    
    template<class SC,class LO,class GO,class NO>
    void CoarseOperator<SC,LO,GO,NO>::applyPhi(MultiVector& x,
                                               MultiVector& y) const
    {
        MultiVectorPtr yCoarseSolveTmp = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(x.getMap(),x.getNumVectors());
        *yCoarseSolveTmp = x;
        
        MultiVectorPtr yCoarse;
        for (int j=GatheringMaps_.size()-1; j>0; j--) {
            yCoarse = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(GatheringMaps_[j-1],x.getNumVectors());
            yCoarse->doImport(*yCoarseSolveTmp,*CoarseSolveExporters_[j],Xpetra::INSERT);
            yCoarseSolveTmp = yCoarse;
        }
        yCoarse = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(CoarseMap_,x.getNumVectors());
        
        yCoarse->doImport(*yCoarseSolveTmp,*CoarseSolveExporters_[0],Xpetra::INSERT);
        
        Phi_->apply(*yCoarse,y,Teuchos::NO_TRANS);
    }
    
    template<class SC,class LO,class GO,class NO>
    int CoarseOperator<SC,LO,GO,NO>::setUpCoarseOperator()
    {
        // Build CoarseMatrix_
        CrsMatrixPtr k0 = buildCoarseMatrix();
        
        // Build CoarseMap_
        buildCoarseSolveMap(k0);
        
        //------------------------------------------------------------------------------------------------------------------------
        // Communicate coarse matrix
        CoarseSolveExporters_[0] = Xpetra::ExportFactory<LO,GO,NO>::Build(CoarseMap_,GatheringMaps_[0]);
        CrsMatrixPtr tmpCoarseMatrix = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(GatheringMaps_[0],k0->getGlobalMaxNumRowEntries());
        tmpCoarseMatrix->doExport(*k0,*CoarseSolveExporters_[0],Xpetra::INSERT);
        
        for (UN j=1; j<GatheringMaps_.size(); j++) {
            tmpCoarseMatrix->fillComplete();
            k0 = tmpCoarseMatrix;
            CoarseSolveExporters_[j] = Xpetra::ExportFactory<LO,GO,NO>::Build(GatheringMaps_[j-1],GatheringMaps_[j]);
            tmpCoarseMatrix = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(GatheringMaps_[j],k0->getGlobalMaxNumRowEntries());
            
            tmpCoarseMatrix->doExport(*k0,*CoarseSolveExporters_[j],Xpetra::INSERT);
        }
        
        //------------------------------------------------------------------------------------------------------------------------
        // Matrix to the new communicator
        if (OnCoarseSolveComm_) {
            CoarseMatrix_ = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(CoarseSolveMap_,k0->getGlobalMaxNumRowEntries());
            ConstGOVecView indices;
            ConstSCVecView values;
            for (UN i=0; i<tmpCoarseMatrix->getNodeNumRows(); i++) {
                tmpCoarseMatrix->getGlobalRowView(CoarseSolveMap_->getGlobalElement(i),indices,values);
                if (indices.size()>0) {
                    CoarseMatrix_->insertGlobalValues(CoarseSolveMap_->getGlobalElement(i),indices,values);
                } else { // Add diagonal unit for zero rows // Todo: Do you we need to sort the coarse matrix "NodeWise"?
                    GOVec indices(1,CoarseSolveMap_->getGlobalElement(i));
                    SCVec values(1,1.0);
                    CoarseMatrix_->insertGlobalValues(CoarseSolveMap_->getGlobalElement(i),indices(),values());
                }
                
            }
            CoarseMatrix_->fillComplete(); //Teuchos::RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)); CoarseMatrix_->describe(*fancy,Teuchos::VERB_EXTREME);
            
            CoarseSolver_.reset(new SubdomainSolver<SC,LO,GO,NO>(CoarseMatrix_,sublist(this->ParameterList_,"CoarseSolver")));
            CoarseSolver_->initialize();
            CoarseSolver_->compute();
        }
        //------------------------------------------------------------------------------------------------------------------------
        
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    typename CoarseOperator<SC,LO,GO,NO>::CrsMatrixPtr CoarseOperator<SC,LO,GO,NO>::buildCoarseMatrix()
    {
        CoarseMap_ = Xpetra::MapFactory<LO,GO,NO>::Build(Phi_->getDomainMap(),1);
        CrsMatrixPtr k0 = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(CoarseMap_,CoarseMap_->getNodeNumElements());
        CrsMatrixPtr tmp = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(this->K_->getRowMap(),50);
        
        Xpetra::MatrixMatrix<SC,LO,GO,NO>::Multiply(*this->K_,false,*Phi_,false,*tmp);
        Xpetra::MatrixMatrix<SC,LO,GO,NO>::Multiply(*Phi_,true,*tmp,false,*k0);
        
        return k0;
    }
    
    /*
    CoarseMap_.reset(new Epetra_Map(Phi_->DomainMap()));
    CrsMatrixPtr k0(new Epetra_CrsMatrix(Copy,*CoarseMap_,CoarseMap_->getNodeNumElements())); // schoener machen mit ColMap
    Epetra_CrsMatrix tmp(Copy,K_->getRowMap(),50);
    EpetraExt::MatrixMatrix::Multiply(*K_,false,*Phi_,false,tmp);
    EpetraExt::MatrixMatrix::Multiply(*Phi_,true,tmp,false,*k0);
     */
    
    template<class SC,class LO,class GO,class NO>
    int CoarseOperator<SC,LO,GO,NO>::buildCoarseSolveMap(CrsMatrixPtr k0)
    {
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
            FROSCH_ASSERT(0!=0,"This should never happen...");
        }
        
        NumProcsCoarseSolve_ = (LO) ((double) NumProcsCoarseSolve_)*fac;
        //cout << NumProcsCoarseSolve_ << std::endl;
        if (NumProcsCoarseSolve_<1) {
            NumProcsCoarseSolve_ = 1;
        }
        
        if (NumProcsCoarseSolve_ >= this->MpiComm_->getSize()) {
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
            GO numGlobalIndices = CoarseMap_->getMaxAllGlobalIndex()+1;
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
                GatheringMaps_[i] = Xpetra::MapFactory<LO,GO,NO>::Build(CoarseMap_->lib(),-1,numMyRows,0,this->MpiComm_);
            }
            
            numMyRows = 0;
            if (this->MpiComm_->getRank()%(this->MpiComm_->getSize()/NumProcsCoarseSolve_) == 0 && this->MpiComm_->getRank()/(this->MpiComm_->getSize()/NumProcsCoarseSolve_) < NumProcsCoarseSolve_) {
                if (this->MpiComm_->getRank()==0) {
                    numMyRows = numGlobalIndices - (numGlobalIndices/NumProcsCoarseSolve_)*(NumProcsCoarseSolve_-1);
                } else {
                    numMyRows = numGlobalIndices/NumProcsCoarseSolve_;
                }
            }
            GatheringMaps_[gatheringSteps-1] = Xpetra::MapFactory<LO,GO,NO>::Build(CoarseMap_->lib(),-1,numMyRows,0,this->MpiComm_);
            //cout << *GatheringMaps_->at(gatheringSteps-1);
            
            //------------------------------------------------------------------------------------------------------------------------
            // Use a separate Communicator for the coarse problem
            MapPtr tmpCoarseMap = GatheringMaps_[GatheringMaps_.size()-1];
            
            if (tmpCoarseMap->getNodeNumElements()>0) {
                OnCoarseSolveComm_=true;
            }
            CoarseSolveComm_ = this->MpiComm_->split(!OnCoarseSolveComm_,this->MpiComm_->getRank());
            CoarseSolveMap_ = Xpetra::MapFactory<LO,GO,NO>::Build(CoarseMap_->lib(),-1,tmpCoarseMap->getNodeElementList(),0,CoarseSolveComm_);
            
        } else {
            FROSCH_ASSERT(0!=0,"Distribution type not defined...");
        }
        
        return 0;
    }
    
}

#endif
