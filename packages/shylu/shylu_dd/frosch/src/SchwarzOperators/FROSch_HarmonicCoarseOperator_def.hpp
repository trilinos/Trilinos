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
    InterfaceCoarseSpaces_ (0),
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
            for (UN j=0; j<InterfaceCoarseSpaces_[i]->getBasisMap()->getNodeNumElements(); j++) {
                    mapVector.push_back(InterfaceCoarseSpaces_[i]->getBasisMap()->getGlobalElement(j)+tmp);
                }
                tmp += InterfaceCoarseSpaces_[i]->getBasisMap()->getMaxAllGlobalIndex()+1;
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
        GammaDofs_->resize(GammaDofs_.size()+1);
        IDofs_->resize(IDofs_.size()+1);
        InterfaceCoarseSpaces_->resize(InterfaceCoarseSpaces_.size()+1);
        DofsMaps_->resize(DofsMaps_.size()+1);
        DofsPerNode_->resize(DofsPerNode_.size()+1);

        NumberOfBlocks_++;

        /////
        int blockId = NumberOfBlocks_-1;

        // Process the parameter list
        std::stringstream blockIdStringstream;
        blockIdStringstream << blockId+1;
        std::string blockIdString = blockIdStringstream.str();
        Teuchos::RCP<Teuchos::ParameterList> coarseSpaceList = sublist(sublist(this->ParameterList_,"Blocks"),blockIdString.c_str());

        bool useForCoarseSpace = coarseSpaceList->get("Use For Coarse Space",true);

        GammaDofs_[blockId] = LOVecPtr(0);

        MultiVectorPtr mVPhiGamma;
        MapPtr blockCoarseMap;
        if (useForCoarseSpace) {
            InterfaceCoarseSpaces_[blockId].reset(new CoarseSpace<SC,LO,GO,NO>());
            
            //Epetra_SerialComm serialComm;
            MapPtr serialGammaMap = Xpetra::MapFactory<LO,GO,NO>::Build(dofsMap->lib(),dofsMap->getNodeNumElements(),0,this->SerialComm_);
            mVPhiGamma = Xpetra::MultiVectorFactory<LO,GO,NO>::Build(serialGammaMap,dofsMap->getNodeNumElements());
        }

        for (int i=0; i<dofsMap->getNodeNumElements(); i++) {
            GammaDofs_[blockId]->push_back(i);

            if (useForCoarseSpace) {
                mVPhiGamma->replaceLocalValue(i,i,1.0);
            }
        }

        IDofs_[blockId] = LOVecPtr(0);

        if (useForCoarseSpace) {
            blockCoarseMap = Xpetra::MapFactory<LO,GO,NO>::Build(dofsMap->lib(),-1,GammaDofs_[blockId](),0,this->MpiComm_);
            
            InterfaceCoarseSpaces_[blockId]->addSubspace(blockCoarseMap,mVPhiGamma);
            InterfaceCoarseSpaces_[blockId]->assembleCoarseSpace();
        }

        DofsMaps_[blockId] = MapPtrVecPtr(0);
        DofsMaps_[blockId].push_back(dofsMap);

        DofsPerNode_[blockId] = 1;

        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    int HarmonicCoarseOperator<SC,LO,GO,NO>::computeVolumeFunctions(UN blockId,
                                                                    UN dimension,
                                                                    MapPtr nodesMap,
                                                                    MultiVectorPtr nodeList,
                                                                    EntitySetPtr interior)
    {
        // Process the parameter list
        std::stringstream blockIdStringstream;
        blockIdStringstream << blockId+1;
        std::string blockIdString = blockIdStringstream.str();
        Teuchos::RCP<Teuchos::ParameterList> coarseSpaceList = sublist(sublist(this->ParameterList_,"Blocks"),blockIdString.c_str());
        
        bool useForCoarseSpace = coarseSpaceList->get("Use For Coarse Space",true);
        bool useRotations = coarseSpaceList->get("Rotations",true);
        
        this->GammaDofs_[blockId] = LOVecPtr(this->DofsPerNode_[blockId]*interior->getEntity(0)->getNumNodes());
        this->IDofs_[blockId] = LOVecPtr(0);
        for (UN k=0; k<this->DofsPerNode_[blockId]; k++) {
            for (UN i=0; i<interior->getEntity(0)->getNumNodes(); i++) {
                this->GammaDofs_[blockId][this->DofsPerNode_[blockId]*i+k] = interior->getEntity(0)->getLocalDofID(i,k);
            }
        }
        
        if (useForCoarseSpace) {            
            InterfaceCoarseSpaces_[blockId].reset(new CoarseSpace<SC,LO,GO,NO>());
            
            interior->buildEntityMap(nodesMap);
            
            MultiVectorPtrVecPtr translations = computeTranslations(blockId,interior);
            for (UN i=0; i<translations.size(); i++) {
                this->InterfaceCoarseSpaces_[blockId]->addSubspace(interior->getEntityMap(),translations[i]);
            }
            if (useRotations) {
                MultiVectorPtrVecPtr rotations = computeRotations(blockId,dimension,nodeList,interior);
                for (UN i=0; i<rotations.size(); i++) {
                    this->InterfaceCoarseSpaces_[blockId]->addSubspace(interior->getEntityMap(),rotations[i]);
                }
            }
            
            InterfaceCoarseSpaces_[blockId]->assembleCoarseSpace();
            
            // Count entities
            GO numEntitiesGlobal = interior->getEntityMap()->getMaxAllGlobalIndex();
            if (interior->getEntityMap()->lib()==Xpetra::UseEpetra || interior->getEntityMap()->getGlobalNumElements()>0) {
                numEntitiesGlobal += 1;
            }
            
            if (this->MpiComm_->getRank() == 0) {
                std::cout << "\n\
                --------------------------------------------\n\
                # volumes:               --- " << numEntitiesGlobal << "\n\
                --------------------------------------------\n\
                Coarse space:\n\
                --------------------------------------------\n\
                volume: translations      --- " << 1 << "\n\
                volume: rotations         --- " << useRotations << "\n\
                --------------------------------------------\n";
            }
        }
        return 0;
    }
    
    template <class SC,class LO,class GO,class NO>
    typename HarmonicCoarseOperator<SC,LO,GO,NO>::MultiVectorPtrVecPtr HarmonicCoarseOperator<SC,LO,GO,NO>::computeTranslations(UN blockId,
                                                                                                                                EntitySetPtr entitySet)
    {
        MultiVectorPtrVecPtr translations(this->DofsPerNode_[blockId]);
        MapPtr serialGammaMap = Xpetra::MapFactory<LO,GO,NO>::Build(this->K_->getRangeMap()->lib(),this->GammaDofs_[blockId].size(),0,this->SerialComm_);
        for (UN i=0; i<this->DofsPerNode_[blockId]; i++) {
            if (entitySet->getNumEntities()>0) {
                translations[i] = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(serialGammaMap,entitySet->getNumEntities());
            } else {
                translations[i] = Teuchos::null;
            }
        }
        
        for (UN k=0; k<this->DofsPerNode_[blockId]; k++) {
            for (UN i=0; i<entitySet->getNumEntities(); i++) {
                for (UN j=0; j<entitySet->getEntity(i)->getNumNodes(); j++) {
                    translations[k]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,k),i,1.0);
                }
            }
        }
        return translations;
    }
    
    template <class SC,class LO,class GO,class NO>
    typename HarmonicCoarseOperator<SC,LO,GO,NO>::MultiVectorPtrVecPtr HarmonicCoarseOperator<SC,LO,GO,NO>::computeRotations(UN blockId,
                                                                                                                             UN dimension,
                                                                                                                             MultiVectorPtr nodeList,
                                                                                                                             EntitySetPtr entitySet)
    {
        FROSCH_ASSERT(nodeList->getNumVectors()==dimension,"dimension of the nodeList is wrong.");
        FROSCH_ASSERT(dimension==this->DofsPerNode_[blockId],"dimension!=this->DofsPerNode_[blockId]");
        
        UN rotationsPerEntity = 0;
        switch (dimension) {
            case 1:
                return Teuchos::null;
                break;
            case 2:
                rotationsPerEntity = 1;
                break;
            case 3:
                rotationsPerEntity = 3;
                break;
            default:
                FROSCH_ASSERT(false,"The dimension is neither 2 nor 3!");
                break;
        }
        
        MultiVectorPtrVecPtr rotations(rotationsPerEntity);
        MapPtr serialGammaMap = Xpetra::MapFactory<LO,GO,NO>::Build(this->K_->getRangeMap()->lib(),this->GammaDofs_[blockId].size(),0,this->SerialComm_);
        for (UN i=0; i<rotationsPerEntity; i++) {
            if (entitySet->getNumEntities()>0) {
                rotations[i] = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(serialGammaMap,entitySet->getNumEntities());
            } else {
                rotations[i] = Teuchos::null;
            }
        }
        
        SC x,y,z,rx,ry,rz;
        for (UN i=0; i<entitySet->getNumEntities(); i++) {
            for (UN j=0; j<entitySet->getEntity(i)->getNumNodes(); j++) {
                x = nodeList->getData(0)[entitySet->getEntity(i)->getLocalNodeID(j)];
                y = nodeList->getData(1)[entitySet->getEntity(i)->getLocalNodeID(j)];
                
                // Rotation 1
                rx = y;
                ry = -x;
                rz = 0;
                rotations[0]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,0),i,rx);
                rotations[0]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,1),i,ry);
                if (dimension == 3) {
                    z = nodeList->getData(2)[entitySet->getEntity(i)->getLocalNodeID(j)];
                    
                    rotations[0]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,2),i,rz);
                    
                    // Rotation 2
                    rx = -z;
                    ry = 0;
                    rz = x;
                    rotations[1]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,0),i,rx);
                    rotations[1]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,1),i,ry);
                    rotations[1]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,2),i,rz);
                    
                    // Rotation 3
                    rx = 0;
                    ry = z;
                    rz = -y;
                    rotations[2]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,0),i,rx);
                    rotations[2]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,1),i,ry);
                    rotations[2]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,2),i,rz);
                }
            }
        }
        return rotations;
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
            if (InterfaceCoarseSpaces_[i]->hasAssembledBasis()) {
                numLocalBlockRows[i] = InterfaceCoarseSpaces_[i]->getLocalBasis()->getNumVectors();
                for (j=0; j<numLocalBlockRows[i]; j++) {
                    for (k=0; k<InterfaceCoarseSpaces_[i]->getLocalBasis()->getLocalLength(); k++) {
                        mVPhiGamma->replaceLocalValue(k+kk,j+jj,InterfaceCoarseSpaces_[i]->getLocalBasis()->getData(j)[k]);
                        mVPhi->replaceLocalValue(indicesGammaDofsAll[k+kk],j+jj,InterfaceCoarseSpaces_[i]->getLocalBasis()->getData(j)[k]);
                    }
                }
            } else { // Das ist für den Fall, dass keine Basisfunktionen für einen Block gebaut werden sollen
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
