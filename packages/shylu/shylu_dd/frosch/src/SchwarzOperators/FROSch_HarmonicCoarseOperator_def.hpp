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
    typename HarmonicCoarseOperator<SC,LO,GO,NO>::MapPtr HarmonicCoarseOperator<SC,LO,GO,NO>::computeCoarseSpace(CoarseSpacePtr coarseSpace)
    {
        MapPtr repeatedMap = assembleSubdomainMap();
        
        // Build local saddle point problem
        CrsMatrixPtr repeatedMatrix = FROSch::ExtractLocalSubdomainMatrix(this->K_,repeatedMap); // AH 12/11/2018: Should this be in initalize?

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
            tmp += GammaDofs_[i].size()+IDofs_[i].size();
        }

        CrsMatrixPtr kII;
        CrsMatrixPtr kIGamma;
        CrsMatrixPtr kGammaI;
        CrsMatrixPtr kGammaGamma;

        FROSch::BuildSubmatrices(repeatedMatrix,indicesIDofsAll(),kII,kIGamma,kGammaI,kGammaGamma);

        // Assemble coarse map
        MapPtr coarseMap = assembleCoarseMap(); // AH 12/11/2018: Should this be in initalize?

        // Build the saddle point harmonic extensions
        MultiVectorPtr localCoarseSpaceBasis = computeExtensions(repeatedMatrix->getRowMap(),coarseMap,indicesGammaDofsAll(),indicesIDofsAll(),kII,kIGamma);

        coarseSpace->addSubspace(coarseMap,localCoarseSpaceBasis);
        
        return repeatedMap;
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
    typename HarmonicCoarseOperator<SC,LO,GO,NO>::MapPtr HarmonicCoarseOperator<SC,LO,GO,NO>::assembleSubdomainMap()
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
    typename HarmonicCoarseOperator<SC,LO,GO,NO>::MultiVectorPtr HarmonicCoarseOperator<SC,LO,GO,NO>::computeExtensions(ConstMapPtr localMap,
                                                                                                                        ConstMapPtr coarseMap,
                                                                                                                        GOVecView indicesGammaDofsAll,
                                                                                                                        GOVecView indicesIDofsAll,
                                                                                                                        CrsMatrixPtr kII,
                                                                                                                        CrsMatrixPtr kIGamma)
    {
        //this->Phi_ = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(this->K_->getRangeMap(),coarseMap,coarseMap->getNodeNumElements()); // Nonzeroes abhängig von dim/dofs!!!
        MultiVectorPtr mVPhi = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(localMap,coarseMap->getNodeNumElements());
        MultiVectorPtr mVtmp = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(kII->getRowMap(),coarseMap->getNodeNumElements());
        MultiVectorPtr mVPhiI = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(kII->getRowMap(),coarseMap->getNodeNumElements());

        //Build mVPhiGamma
        MultiVectorPtr mVPhiGamma = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(kIGamma->getDomainMap(),coarseMap->getNodeNumElements());
        LO jj=0;
        LO kk=0;
        UNVec numLocalBlockRows(NumberOfBlocks_);
        for (UN i=0; i<NumberOfBlocks_; i++) {
            UN j = 0;
            UN k = 0;
            //if (this->Verbose_) std::cout << MVPhiGamma_[i]->MyLength() << std::endl;
            if (!MVPhiGamma_[i].is_null()) {
                numLocalBlockRows[i] = MVPhiGamma_[i]->getNumVectors();
                for (j=0; j<MVPhiGamma_[i]->getNumVectors(); j++) {
                    for (k=0; k<MVPhiGamma_[i]->getLocalLength(); k++) {
                        //if (this->Verbose_) std::cout << j << " " << k << " " <<  (*(*MVPhiGamma_[i])(j))[k] << std::endl;
                        mVPhiGamma->replaceLocalValue(k+kk,j+jj,MVPhiGamma_[i]->getData(j)[k]);
                        mVPhi->replaceLocalValue(indicesGammaDofsAll[k+kk],j+jj,MVPhiGamma_[i]->getData(j)[k]);
                    }
                }
            } else { // Das ist für den Fall, dass keine Basisfunktionen für einen Block gebaut werden sollen
                //mVPhiGamma->replaceLocalValue(k+kk,j+jj,1.0);
                k=GammaDofs_[i].size();
            }
            jj += j;
            kk += k;
        }
        // Teuchos::RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)); this->Phi_->describe(*fancy,Teuchos::VERB_EXTREME);
        // Hier Multiplikation kIGamma*PhiGamma
        kIGamma->apply(*mVPhiGamma,*mVtmp);
        
        mVtmp->scale(-1.0);
        
        // Jetzt der solver für kII
        ExtensionSolver_.reset(new SubdomainSolver<SC,LO,GO,NO>(kII,sublist(this->ParameterList_,"ExtensionSolver")));
        ExtensionSolver_->initialize();
        ExtensionSolver_->compute();
        ExtensionSolver_->apply(*mVtmp,*mVPhiI);

        GOVec priorIndex(NumberOfBlocks_,0);
        GOVec postIndex(NumberOfBlocks_,0);
        
        GOVec2D excludeCols(NumberOfBlocks_);
        
        LOVec bound(NumberOfBlocks_+1,(LO)0);
        for (UN i=0; i<NumberOfBlocks_; i++) {
            bound[i+1] = bound[i] + this->IDofs_[i].size();
        }
        
        LO itmp = 0;
        for (UN i=0; i<NumberOfBlocks_; i++) {
            std::stringstream blockIdStringstream;
            blockIdStringstream << i+1;
            std::string blockIdString = blockIdStringstream.str();
            Teuchos::RCP<Teuchos::ParameterList> coarseSpaceList = sublist(sublist(this->ParameterList_,"Blocks"),blockIdString.c_str());
            std::string excludeBlocksString = coarseSpaceList->get("Exclude","0");
            UNVec excludeBlocks = FROSch::GetIndicesFromString(excludeBlocksString,(UN)0);
            sortunique(excludeBlocks);
            for (UN j=0; j<excludeBlocks.size(); j++) {
                excludeBlocks[j] -= 1;
            }
            UNVec extensionBlocks(0);
            for (UN j=0; j<NumberOfBlocks_; j++) {
                typename UNVec::iterator it = std::find(excludeBlocks.begin(),excludeBlocks.end(),j);
                if (it == excludeBlocks.end()) {
                    extensionBlocks.push_back(j);
                }
            }
            
            for (UN j=0; j<numLocalBlockRows[i]; j++) {
                for (UN ii=0; ii<extensionBlocks.size(); ii++) {
                    for (LO k=bound[extensionBlocks[ii]]; k<bound[extensionBlocks[ii]+1]; k++) {
                        if (k>=indicesIDofsAll.size()) {
                            std::cout << k << " " << indicesIDofsAll.size() << std::endl;
                        }
                        indicesIDofsAll[k];
                        mVPhiI->getData(itmp)[k];
                        mVPhi->replaceLocalValue(indicesIDofsAll[k],itmp,mVPhiI->getData(itmp)[k]);
                    }
                }
                itmp++;
            }
        }
        return mVPhi;
    }
    
}

#endif
