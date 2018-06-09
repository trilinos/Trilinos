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

#ifndef _FROSCH_HARMONICCOARSEOPERATOR_DEF_HPP
#define _FROSCH_HARMONICCOARSEOPERATOR_DEF_HPP

#include <FROSch_HarmonicCoarseOperator_decl.hpp>

namespace FROSch {
    
    template <class SC,class LO,class GO,class NO>
    HarmonicCoarseOperator<SC,LO,GO,NO>::HarmonicCoarseOperator(CrsMatrixPtr k,
                                                                ParameterListPtr parameterList) :
    CoarseOperator<SC,LO,GO,NO> (k,parameterList),
    ExtensionSolver_ (),
    MVPhiGamma_ (0),
    BlockCoarseMaps_ (0),
    Dimensions_ (0),
    DofsPerNode_ (0),
    GammaDofs_ (0),
    IDofs_ (0),
    DofsMaps_ (0),
    NumberOfBlocks_ (0)
    {
        
    }
    
    template <class SC,class LO,class GO,class NO>
    int HarmonicCoarseOperator<SC,LO,GO,NO>::compute()
    {
        // This is not optimal yet... Some work could be moved to Initialize
        if (this->Verbose_) {
            cerr << "WARNING: Some of the operations could be moved from initialize() to Compute().\n";
        }
        this->computeHarmonicExtensions();
        this->setUpCoarseOperator();
        this->IsComputed_ = true;
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int HarmonicCoarseOperator<SC,LO,GO,NO>::computeHarmonicExtensions()
    {
        // Build repeatedMap for the local saddle point problem
        MapPtr repeatedMap = assembleRepeatedMap(); // Todo:: Eigentlich gehört das in Initialize !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        // Build local saddle point problem
        CrsMatrixPtr repeatedMatrix = FROSch::ExtractLocalSubdomainMatrix(this->K_,repeatedMap);

        // Extract submatrices
        GOVec indicesGammaDofsAll(0);
        GOVec indicesIDofsAll(0);
        LO tmp = 0;

        for (UN i=0; i<NumberOfBlocks_; i++) {
            for (UN j=0; j<GammaDofs_[i].size(); j++) {
                indicesGammaDofsAll.push_back(tmp+GammaDofs_[i][j]);
            }
            for (UN j=0; j<IDofs_[i].size(); j++) {
                indicesIDofsAll.push_back(tmp+IDofs_[i][j]);
            }
            tmp += GammaDofs_[i].size()+IDofs_[i].size(); // Ist das mit tmp korrekt?
        }

        CrsMatrixPtr kII;
        CrsMatrixPtr kIGamma;
        CrsMatrixPtr kGammaI;
        CrsMatrixPtr kGammaGamma;

        FROSch::BuildSubmatrices(repeatedMatrix,indicesIDofsAll(),kII,kIGamma,kGammaI,kGammaGamma);

        // Assemble coarse map
        MapPtr coarseMap = assembleCoarseMap(); // Todo:: Eigentlich gehört das in Initialize !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        // Build the saddle point harmonic extensions
        computeAndFillPhi(repeatedMatrix,repeatedMap,coarseMap,indicesGammaDofsAll(),indicesIDofsAll(),kII,kIGamma);

        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    typename HarmonicCoarseOperator<SC,LO,GO,NO>::MapPtr HarmonicCoarseOperator<SC,LO,GO,NO>::assembleRepeatedMap()
    {
        FROSCH_ASSERT(DofsMaps_.size()==NumberOfBlocks_,"DofsMaps_.size()!=NumberOfBlocks_");
        FROSCH_ASSERT(DofsPerNode_.size()==NumberOfBlocks_,"DofsPerNode_.size()!=NumberOfBlocks_");
        
        GOVec mapVector(0);
        for (UN i=0; i<NumberOfBlocks_; i++) {
            FROSCH_ASSERT(DofsMaps_[i].size()==DofsPerNode_[i],"DofsMaps_[i].size()!=DofsPerNode_[i]");
            UN numMyElements = DofsMaps_[i][0]->getNodeNumElements();
            for (UN j=1; j<DofsPerNode_[i]; j++) {
                FROSCH_ASSERT(DofsMaps_[i][j]->getNodeNumElements()==(unsigned) numMyElements,"DofsMaps_[i][j]->getNodeNumElements()==numMyElements");
            }
            for (UN j=0; j<numMyElements; j++) {
                for (UN k=0; k<DofsPerNode_[i]; k++) {
                    mapVector.push_back(DofsMaps_[i][k]->getGlobalElement(j));
                }
            }
        }
        return Xpetra::MapFactory<LO,GO,NO>::Build(DofsMaps_[0][0]->lib(),-1,mapVector(),0,this->MpiComm_);
    }
    
    template <class SC,class LO,class GO,class NO>
    typename HarmonicCoarseOperator<SC,LO,GO,NO>::MapPtr HarmonicCoarseOperator<SC,LO,GO,NO>::assembleCoarseMap()
    {
        GOVec mapVector(0);
        GO tmp = 0;
        for (UN i=0; i<NumberOfBlocks_; i++) {
            if (!BlockCoarseMaps_[i].is_null()) {
                for (UN j=0; j<BlockCoarseMaps_[i]->getNodeNumElements(); j++) {
                    mapVector.push_back(BlockCoarseMaps_[i]->getGlobalElement(j)+tmp);
                }
                tmp += BlockCoarseMaps_[i]->getMaxAllGlobalIndex()+1;
            }
        }
        return Xpetra::MapFactory<LO,GO,NO>::Build(DofsMaps_[0][0]->lib(),-1,mapVector(),0,this->MpiComm_);
    }
    
    template <class SC,class LO,class GO,class NO>
        int  HarmonicCoarseOperator<SC,LO,GO,NO>::addZeroCoarseSpaceBlock(MapPtr dofsMap)
    	{
        	// Das könnte man noch ändern
        	this->GammaDofs_->resize(this->GammaDofs_.size()+1);
        	this->IDofs_->resize(this->IDofs_.size()+1);
        	this->BlockCoarseMaps_->resize(this->BlockCoarseMaps_.size()+1);
        	this->MVPhiGamma_->resize(this->MVPhiGamma_.size()+1);
        	this->DofsMaps_->resize(this->DofsMaps_.size()+1);
        	this->DofsPerNode_->resize(this->DofsPerNode_.size()+1);

        	this->NumberOfBlocks_++;

        	/////
    		int blockId = this->NumberOfBlocks_-1;

    		// Process the parameter list
    		std::stringstream blockIdStringstream;
    		blockIdStringstream << blockId+1;
    		std::string blockIdString = blockIdStringstream.str();
    		Teuchos::RCP<Teuchos::ParameterList> coarseSpaceList = sublist(sublist(this->ParameterList_,"Blocks"),blockIdString.c_str());

    		bool useForCoarseSpace = coarseSpaceList->get("Use For Coarse Space",true);

    		this->GammaDofs_[blockId] = LOVecPtr(0);

    		if (useForCoarseSpace) {
    			//Epetra_SerialComm serialComm;
    			MapPtr serialGammaMap = Xpetra::MapFactory<LO,GO,NO>::Build(dofsMap->lib(),dofsMap->getNodeNumElements(),0,this->SerialComm_);
    			this->MVPhiGamma_[blockId] = Xpetra::MultiVectorFactory<LO,GO,NO>::Build(serialGammaMap,dofsMap->getNodeNumElements());
    		}

    		for (int i=0; i<dofsMap->getNodeNumElements(); i++) {
    			this->GammaDofs_[blockId]->push_back(i);

    			if (useForCoarseSpace) {
    				this->MVPhiGamma_[blockId]->replaceLocalValue(i,i,1.0);
    			}
    		}

    		this->IDofs_[blockId] = LOVecPtr(0);

    		if (useForCoarseSpace) {
    			this->BlockCoarseMaps_[blockId] = Xpetra::MapFactory<LO,GO,NO>::Build(dofsMap->lib(),-1,this->GammaDofs_[blockId](),0,this->MpiComm_);
    		}

    		this->DofsMaps_[blockId] = MapPtrVecPtr(0);
    		this->DofsMaps_[blockId].push_back(dofsMap);

    		this->DofsPerNode_[blockId] = 1;


    		return 0;
    	}

    template <class SC,class LO,class GO,class NO>
    int HarmonicCoarseOperator<SC,LO,GO,NO>::computeAndFillPhi(CrsMatrixPtr repeatedMatrix,
                                                           MapPtr repeatedMap,
                                                           MapPtr coarseMap,
                                                           GOVecView indicesGammaDofsAll,
                                                           GOVecView indicesIDofsAll,
                                                           CrsMatrixPtr kII,
                                                           CrsMatrixPtr kIGamma)
    {
        this->Phi_ = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(this->K_->getRangeMap(),coarseMap,coarseMap->getNodeNumElements()); // Nonzeroes abhängig von dim/dofs!!!
        MultiVectorPtr mVtmp = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(kII->getRowMap(),coarseMap->getNodeNumElements());
        MultiVectorPtr mVPhiI = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(kII->getRowMap(),coarseMap->getNodeNumElements());
        
        //Build mVPhiGamma
        MultiVectorPtr mVPhiGamma = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(kIGamma->getDomainMap(),coarseMap->getNodeNumElements());
        LO jj=0;
        LO kk=0;
        for (UN i=0; i<NumberOfBlocks_; i++) {
            UN j = 0;
            UN k = 0;
            //if (this->Verbose_) std::cout << MVPhiGamma_[i]->MyLength() << std::endl;
            if (!MVPhiGamma_[i].is_null()) {
                for (j=0; j<MVPhiGamma_[i]->getNumVectors(); j++) {
                    for (k=0; k<MVPhiGamma_[i]->getLocalLength(); k++) {
                        //if (this->Verbose_) std::cout << j << " " << k << " " <<  (*(*MVPhiGamma_[i])(j))[k] << std::endl;
                        mVPhiGamma->replaceLocalValue(k+kk,j+jj,MVPhiGamma_[i]->getData(j)[k]);
                    }
                }
            } else { // Das ist für den Fall, dass keine Basisfunktionen für einen Block gebaut werden sollen
                //mVPhiGamma->replaceLocalValue(k+kk,j+jj,1.0);
                k=GammaDofs_[i].size();
            }
            jj += j;
            kk += k;
        }
        
        LO iD;
        SC valueTmp;
        LOVec indices;
        SCVec values;
        // IST DAS LANGSAM??????? TESTEN!!!!!
        for (UN i=0; i<mVPhiGamma->getLocalLength(); i++) {
            indices.resize(0);
            values.resize(0);
            for (UN j=0; j<mVPhiGamma->getNumVectors(); j++) {
                valueTmp=mVPhiGamma->getData(j)[i];
                if (fabs(valueTmp) > 1.0e-8) {
                    indices.push_back(j);
                    values.push_back(valueTmp);
                }
            }
            
            // CHECK, IF OK?!?!?
            iD = this->K_->getRowMap()->getLocalElement(repeatedMap->getGlobalElement(indicesGammaDofsAll[i]));
            //cout << ProcessorMapRepeatedDofs->GID(IndicesGamma[0][i]) << " " << ID << std::endl;
            
            //ID = IndicesGamma[0][i];
            if (iD!=-1) {
                this->Phi_->insertLocalValues(iD,indices(),values());
            }
            
        }
//        Teuchos::RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)); this->Phi_->describe(*fancy,Teuchos::VERB_EXTREME);
        // Hier Multiplikation kIGamma*PhiGamma
        kIGamma->apply(*mVPhiGamma,*mVtmp);
        
        mVtmp->scale(-1.0);
        
        // Jetzt der solver für kII
        ExtensionSolver_.reset(new SubdomainSolver<SC,LO,GO,NO>(kII,sublist(this->ParameterList_,"ExtensionSolver")));
        // DAS MÜSSEN WIR NOCH ÄNDERN -> initialize, compute, apply...
        ExtensionSolver_->initialize();
        ExtensionSolver_->compute();
        ExtensionSolver_->apply(*mVtmp,*mVPhiI);
        
        // Now we have to insert the values to the global matrix Phi
        // IST DAS LANGSAM??????? TESTEN!!!!!
        for (UN i=0; i<mVPhiI->getLocalLength(); i++) {
            indices.resize(0);
            values.resize(0);
            for (UN j=0; j<mVPhiI->getNumVectors(); j++) {
                valueTmp=mVPhiI->getData(j)[i]; //if (this->Verbose_) std::cout << i << " " << this->K_->getRowMap()->getLocalElement(repeatedMap->getGlobalElement(indicesIDofsAll[i])) << " " << j << " " << valueTmp << std::endl;
                if (fabs(valueTmp) > 1.0e-8) {
                    indices.push_back(j);
                    values.push_back(valueTmp);
                }
            }
            
            // CHECK, IF OK?!?!?
            iD = this->K_->getRowMap()->getLocalElement(repeatedMap->getGlobalElement(indicesIDofsAll[i]));
            //ID = IndicesI[0][i];
            if (iD!=-1) {
                this->Phi_->insertLocalValues(iD,indices(),values());
            }
        }
        this->Phi_->fillComplete(coarseMap,this->Phi_->getRowMap());
//        Teuchos::RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)); this->Phi_->describe(*fancy,Teuchos::VERB_EXTREME);
        return 0;
    }
    
}

#endif
