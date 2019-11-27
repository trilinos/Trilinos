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

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    HarmonicCoarseOperator<SC,LO,GO,NO>::HarmonicCoarseOperator(ConstXMatrixPtr k,
                                                                ParameterListPtr parameterList) :
    CoarseOperator<SC,LO,GO,NO> (k,parameterList),
    ExtensionSolver_ (),
    InterfaceCoarseSpaces_ (0),
    AssembledInterfaceCoarseSpace_ (new CoarseSpace<SC,LO,GO,NO>(this->MpiComm_,this->SerialComm_)),
    Dimensions_ (0),
    DofsPerNode_ (0),
    GammaDofs_ (0),
    IDofs_ (0),
    DofsMaps_ (0),
    NumberOfBlocks_ (0)
    {
        FROSCH_TIMER_START_LEVELID(harmonicCoarseOperatorTime,"HarmonicCoarseOperator::HarmonicCoarseOperator");
    }

    template <class SC,class LO,class GO,class NO>
    typename HarmonicCoarseOperator<SC,LO,GO,NO>::XMapPtr HarmonicCoarseOperator<SC,LO,GO,NO>::computeCoarseSpace(CoarseSpacePtr coarseSpace)
    {
        FROSCH_TIMER_START_LEVELID(computeCoarseSpaceTime,"HarmonicCoarseOperator::computeCoarseSpace");
        XMapPtr repeatedMap = AssembleSubdomainMap(NumberOfBlocks_,DofsMaps_,DofsPerNode_);

        // Build local saddle point problem
        ConstXMatrixPtr repeatedMatrix = ExtractLocalSubdomainMatrix(this->K_.getConst(),repeatedMap.getConst()); // AH 12/11/2018: Should this be in initalize?

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

        XMatrixPtr kII;
        XMatrixPtr kIGamma;
        XMatrixPtr kGammaI;
        XMatrixPtr kGammaGamma;

        BuildSubmatrices(repeatedMatrix,indicesIDofsAll(),kII,kIGamma,kGammaI,kGammaGamma);

        // Build the saddle point harmonic extensions
        XMultiVectorPtr localCoarseSpaceBasis;
        if (AssembledInterfaceCoarseSpace_->getBasisMap()->getNodeNumElements()) {
            localCoarseSpaceBasis = computeExtensions(repeatedMatrix->getRowMap(),AssembledInterfaceCoarseSpace_->getBasisMap(),indicesGammaDofsAll(),indicesIDofsAll(),kII,kIGamma);
            
            coarseSpace->addSubspace(AssembledInterfaceCoarseSpace_->getBasisMap(),AssembledInterfaceCoarseSpace_->getBasisMapUnique(),localCoarseSpaceBasis);
        } else {
            if (this->Verbose_) std::cout << "FROSch::HarmonicCoarseOperator : WARNING: The Coarse Space is empty. No extensions are computed" << std::endl;
        }

        return repeatedMap;
    }

    template <class SC,class LO,class GO,class NO>
    int HarmonicCoarseOperator<SC,LO,GO,NO>::assembleInterfaceCoarseSpace()
    {
        FROSCH_TIMER_START_LEVELID(assembleInterfaceCoarseSpaceTime,"HarmonicCoarseOperator::assembleInterfaceCoarseSpace");
        LO ii=0;
        for (UN i=0; i<NumberOfBlocks_; i++) {
            if (!InterfaceCoarseSpaces_[i].is_null()) {
                if (InterfaceCoarseSpaces_[i]->hasBasisMap()) {
                    FROSCH_ASSERT(InterfaceCoarseSpaces_[i]->hasBasisMapUnique(),"FROSch::HarmonicCoarseOperator : ERROR: !InterfaceCoarseSpaces_[i]->hasAssembledBasis()");
                    this->AssembledInterfaceCoarseSpace_->addSubspace(InterfaceCoarseSpaces_[i]->getBasisMap(),InterfaceCoarseSpaces_[i]->getBasisMapUnique(),InterfaceCoarseSpaces_[i]->getAssembledBasis(),ii);
                }
            }
            ii += InterfaceCoarseSpaces_[i]->getAssembledBasis()->getLocalLength();
            InterfaceCoarseSpaces_[i].reset();
        }
        return this->AssembledInterfaceCoarseSpace_->assembleCoarseSpace();
    }

    template <class SC,class LO,class GO,class NO>
    int HarmonicCoarseOperator<SC,LO,GO,NO>::addZeroCoarseSpaceBlock(ConstXMapPtr dofsMap)
    {
        FROSCH_TIMER_START_LEVELID(addZeroCoarseSpaceBlockTime,"HarmonicCoarseOperator::addZeroCoarseSpaceBlock");
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
        RCP<ParameterList> coarseSpaceList = sublist(sublist(this->ParameterList_,"Blocks"),blockIdString.c_str());

        bool useForCoarseSpace = coarseSpaceList->get("Use For Coarse Space",true);

        GammaDofs_[blockId] = LOVecPtr(0);

        XMultiVectorPtr mVPhiGamma;
        XMapPtr blockCoarseMap;
        if (useForCoarseSpace) {
            InterfaceCoarseSpaces_[blockId].reset(new CoarseSpace<SC,LO,GO,NO>(this->MpiComm_,this->SerialComm_));

            //Epetra_SerialComm serialComm;
            XMapPtr serialGammaMap = MapFactory<LO,GO,NO>::Build(dofsMap->lib(),dofsMap->getNodeNumElements(),0,this->SerialComm_);
            mVPhiGamma = MultiVectorFactory<LO,GO,NO>::Build(serialGammaMap,dofsMap->getNodeNumElements());
        }

        for (int i=0; i<dofsMap->getNodeNumElements(); i++) {
            GammaDofs_[blockId]->push_back(i);

            if (useForCoarseSpace) {
                mVPhiGamma->replaceLocalValue(i,i,ScalarTraits<SC>::one());
            }
        }

        IDofs_[blockId] = LOVecPtr(0);

        if (useForCoarseSpace) {
            blockCoarseMap = MapFactory<LO,GO,NO>::Build(dofsMap->lib(),-1,GammaDofs_[blockId](),0,this->MpiComm_);

            InterfaceCoarseSpaces_[blockId]->addSubspace(blockCoarseMap,mVPhiGamma);
            InterfaceCoarseSpaces_[blockId]->assembleCoarseSpace();
        }

        DofsMaps_[blockId] = XMapPtrVecPtr(0);
        DofsMaps_[blockId].push_back(dofsMap);

        DofsPerNode_[blockId] = 1;

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int HarmonicCoarseOperator<SC,LO,GO,NO>::computeVolumeFunctions(UN blockId,
                                                                    UN dimension,
                                                                    ConstXMapPtr nodesMap,
                                                                    ConstXMultiVectorPtr nodeList,
                                                                    EntitySetConstPtr interior)
    {
        FROSCH_TIMER_START_LEVELID(computeVolumeFunctionsTime,"HarmonicCoarseOperator::computeVolumeFunctions");
        // Process the parameter list
        std::stringstream blockIdStringstream;
        blockIdStringstream << blockId+1;
        std::string blockIdString = blockIdStringstream.str();
        RCP<ParameterList> coarseSpaceList = sublist(sublist(this->ParameterList_,"Blocks"),blockIdString.c_str());

        bool useForCoarseSpace = coarseSpaceList->get("Use For Coarse Space",true);
        bool useRotations = coarseSpaceList->get("Rotations",true);
        if (useRotations && nodeList.is_null()) {
            useRotations = false;
            if (this->Verbose_) std::cout << "FROSch::HarmonicCoarseOperator : WARNING: Rotations cannot be used" << std::endl;
        }

        this->GammaDofs_[blockId] = LOVecPtr(this->DofsPerNode_[blockId]*interior->getEntity(0)->getNumNodes());
        this->IDofs_[blockId] = LOVecPtr(0);
        for (UN k=0; k<this->DofsPerNode_[blockId]; k++) {
            for (UN i=0; i<interior->getEntity(0)->getNumNodes(); i++) {
                this->GammaDofs_[blockId][this->DofsPerNode_[blockId]*i+k] = interior->getEntity(0)->getLocalDofID(i,k);
            }
        }

        if (useForCoarseSpace) {
            InterfaceCoarseSpaces_[blockId].reset(new CoarseSpace<SC,LO,GO,NO>(this->MpiComm_,this->SerialComm_));

            interior->buildEntityMap(nodesMap);

            XMultiVectorPtrVecPtr translations = computeTranslations(blockId,interior);
            for (UN i=0; i<translations.size(); i++) {
                this->InterfaceCoarseSpaces_[blockId]->addSubspace(interior->getEntityMap(),null,translations[i]);
            }
            if (useRotations) {
                XMultiVectorPtrVecPtr rotations = computeRotations(blockId,dimension,nodeList,interior);
                for (UN i=0; i<rotations.size(); i++) {
                    this->InterfaceCoarseSpaces_[blockId]->addSubspace(interior->getEntityMap(),null,rotations[i]);
                }
            }

            InterfaceCoarseSpaces_[blockId]->assembleCoarseSpace();

            // Count entities
            GO numEntitiesGlobal = interior->getEntityMap()->getMaxAllGlobalIndex();
            if (interior->getEntityMap()->lib()==UseEpetra || interior->getEntityMap()->getGlobalNumElements()>0) {
                numEntitiesGlobal += 1;
            }

            if (this->MpiComm_->getRank() == 0) {
                std::cout << std::boolalpha << "\n\
    ------------------------------------------------------------------------------\n\
     GDSW coarse space\n\
    ------------------------------------------------------------------------------\n\
      Volumes: translations                       --- " << useForCoarseSpace << "\n\
      Volumes: rotations                          --- " << useRotations << "\n\
    ------------------------------------------------------------------------------\n" << std::noboolalpha;
            }
        }
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    typename HarmonicCoarseOperator<SC,LO,GO,NO>::XMultiVectorPtrVecPtr HarmonicCoarseOperator<SC,LO,GO,NO>::computeTranslations(UN blockId,
                                                                                                                                 EntitySetConstPtr entitySet)
    {
        FROSCH_TIMER_START_LEVELID(computeTranslationsTime,"HarmonicCoarseOperator::computeTranslations");
        XMultiVectorPtrVecPtr translations(this->DofsPerNode_[blockId]);
        XMapPtr serialGammaMap = MapFactory<LO,GO,NO>::Build(this->K_->getRangeMap()->lib(),this->GammaDofs_[blockId].size(),0,this->SerialComm_);
        for (UN i=0; i<this->DofsPerNode_[blockId]; i++) {
            if (entitySet->getNumEntities()>0) {
                translations[i] = MultiVectorFactory<SC,LO,GO,NO>::Build(serialGammaMap,entitySet->getNumEntities());
            } else {
                translations[i] = null;
            }
        }

        for (UN k=0; k<this->DofsPerNode_[blockId]; k++) {
            for (UN i=0; i<entitySet->getNumEntities(); i++) {
                for (UN j=0; j<entitySet->getEntity(i)->getNumNodes(); j++) {
                    translations[k]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,k),i,ScalarTraits<SC>::one());
                }
            }
        }
        return translations;
    }

    template <class SC,class LO,class GO,class NO>
    typename HarmonicCoarseOperator<SC,LO,GO,NO>::XMultiVectorPtrVecPtr HarmonicCoarseOperator<SC,LO,GO,NO>::computeRotations(UN blockId,
                                                                                                                              UN dimension,
                                                                                                                              ConstXMultiVectorPtr nodeList,
                                                                                                                              EntitySetConstPtr entitySet,
                                                                                                                              UN discardRotations)
    {
        FROSCH_TIMER_START_LEVELID(computeRotationsTime,"HarmonicCoarseOperator::computeRotations");
        FROSCH_ASSERT(nodeList->getNumVectors()==dimension,"FROSch::HarmonicCoarseOperator : ERROR: Dimension of the nodeList is wrong.");
        FROSCH_ASSERT(dimension==this->DofsPerNode_[blockId],"FROSch::HarmonicCoarseOperator : ERROR: Dimension!=this->DofsPerNode_[blockId]");

        UN rotationsPerEntity = 0;
        switch (dimension) {
            case 1:
                return null;
                break;
            case 2:
                rotationsPerEntity = 1;
                break;
            case 3:
                rotationsPerEntity = 3;
                break;
            default:
                FROSCH_ASSERT(false,"FROSch::HarmonicCoarseOperator : ERROR: The dimension is neither 2 nor 3!");
                break;
        }

        XMultiVectorPtrVecPtr rotations(rotationsPerEntity);
        XMapPtr serialGammaMap = MapFactory<LO,GO,NO>::Build(this->K_->getRangeMap()->lib(),this->GammaDofs_[blockId].size(),0,this->SerialComm_);
        for (UN i=0; i<rotationsPerEntity; i++) {
            if (entitySet->getNumEntities()>0) {
                rotations[i] = MultiVectorFactory<SC,LO,GO,NO>::Build(serialGammaMap,entitySet->getNumEntities());
            } else {
                rotations[i] = null;
            }
        }

        SC x,y,z,rx,ry,rz;
        SCVec rx0(3,0.0),ry0(3,0.0),rz0(3,0.0),errx(3,0.0),erry(3,0.0),errz(3,0.0);
        for (UN i=0; i<entitySet->getNumEntities(); i++) {
            // Compute values for the first node to check if rotation is constant
            x = nodeList->getData(0)[entitySet->getEntity(i)->getLocalNodeID(0)];
            y = nodeList->getData(1)[entitySet->getEntity(i)->getLocalNodeID(0)];
            
            rx0[0] = y;
            ry0[0] = -x;
            if (dimension == 3) {
                z = nodeList->getData(2)[entitySet->getEntity(i)->getLocalNodeID(0)];
                
                rz0[0] = 0;
                
                rx0[1] = -z;
                ry0[1] = 0;
                rz0[1] = x;
                
                rx0[2] = 0;
                ry0[2] = z;
                rz0[2] = -y;
            }
            for (UN j=0; j<entitySet->getEntity(i)->getNumNodes(); j++) {
                x = nodeList->getData(0)[entitySet->getEntity(i)->getLocalNodeID(j)];
                y = nodeList->getData(1)[entitySet->getEntity(i)->getLocalNodeID(j)];

                ////////////////
                // Rotation 1 //
                ////////////////
                rx = y; 
                ry = -x;
                rotations[0]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,0),i,rx);
                rotations[0]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,1),i,ry);
                
                // Compute difference (Euclidean norm of error to constant function)
                errx[0] += (rx0[0]-rx)*(rx0[0]-rx);
                erry[0] += (ry0[0]-ry)*(ry0[0]-ry);
                
                if (dimension == 3) {
                    z = nodeList->getData(2)[entitySet->getEntity(i)->getLocalNodeID(j)];
                    
                    rz = 0;
                    rotations[0]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,2),i,rz);
                    
                    // Compute difference (Euclidean norm of error to constant function)
                    errz[0] += (rz0[0]-rz)*(rz0[0]-rz);

                    ////////////////
                    // Rotation 2 //
                    ////////////////
                    rx = -z;
                    ry = 0;
                    rz = x;
                    rotations[1]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,0),i,rx);
                    rotations[1]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,1),i,ry);
                    rotations[1]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,2),i,rz);
                    
                    // Compute difference (Euclidean norm of error to constant function)
                    errx[1] += (rx0[1]-rx)*(rx0[1]-rx);
                    erry[1] += (ry0[1]-ry)*(ry0[1]-ry);
                    errz[1] += (rz0[1]-rz)*(rz0[1]-rz);

                    ////////////////
                    // Rotation 3 //
                    ////////////////
                    rx = 0;
                    ry = z;
                    rz = -y;
                    rotations[2]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,0),i,rx);
                    rotations[2]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,1),i,ry);
                    rotations[2]->replaceLocalValue(entitySet->getEntity(i)->getGammaDofID(j,2),i,rz);
                    
                    // Compute difference (Euclidean norm of error to constant function)
                    errx[2] += (rx0[2]-rx)*(rx0[2]-rx);
                    erry[2] += (ry0[2]-ry)*(ry0[2]-ry);
                    errz[2] += (rz0[2]-rz)*(rz0[2]-rz);
                }
            }
            
            // If error to constant function is almost zero => scale rotation with zero
            SC err;
            UN numZeroRotations = 0;
            switch (dimension) {
                case 2:
                    err = errx[0]+erry[0];
                    err = ScalarTraits<SC>::squareroot(err);
                    if (std::fabs(err)<1.0e-12) {
                        FROSCH_ASSERT(false,"FROSch::HarmonicCoarseOperator : ERROR: In 2D, no rotation can be constant!");
                        rotations[0]->getVectorNonConst(i)->scale(ScalarTraits<SC>::zero());
                        numZeroRotations++;
                    }
                    break;
                case 3:
                    for (UN j=0; j<3; j++) {
                        err = errx[j]+erry[j]+errz[j];
                        err = ScalarTraits<SC>::squareroot(err);
                        if (std::fabs(err)<1.0e-12) {
                            rotations[j]->getVectorNonConst(i)->scale(ScalarTraits<SC>::zero());
                            numZeroRotations++;
                        }
                    }
                    break;
                default:
                    FROSCH_ASSERT(false,"FROSch::HarmonicCoarseOperator : ERROR: The dimension is neither 1 nor 2 nor 3!");
                    break;
            }
            // If necessary, discard additional rotations
            UN rotationsToDiscard = discardRotations - numZeroRotations;
            if (rotationsToDiscard<0) {
                if (this->Verbose_) std::cout << "FROSch::HarmonicCoarseOperator : WARNING: More rotations have been discarded than expected." << std::endl;
            } else if (rotationsToDiscard>0) {
                UN it=0;
                UN rotationsDiscarded=0;
                while (rotationsDiscarded<rotationsToDiscard) {
                    if (rotations[it]->getVector(i)->norm2()>1.0e-12) {
                        rotations[it]->getVectorNonConst(i)->scale(ScalarTraits<SC>::zero());
                        rotationsDiscarded++;
                    }
                    it++;
                }
            }
        }
        return rotations;
    }

    template <class SC,class LO,class GO,class NO>
    typename HarmonicCoarseOperator<SC,LO,GO,NO>::XMultiVectorPtr HarmonicCoarseOperator<SC,LO,GO,NO>::computeExtensions(ConstXMapPtr localMap,
                                                                                                                         ConstXMapPtr coarseMap,
                                                                                                                         GOVecView indicesGammaDofsAll,
                                                                                                                         GOVecView indicesIDofsAll,
                                                                                                                         XMatrixPtr kII,
                                                                                                                         XMatrixPtr kIGamma)
    {
        FROSCH_TIMER_START_LEVELID(computeExtensionsTime,"HarmonicCoarseOperator::computeExtensions");
        //this->Phi_ = MatrixFactory<SC,LO,GO,NO>::Build(this->K_->getRangeMap(),coarseMap,coarseMap->getNodeNumElements()); // Nonzeroes abhängig von dim/dofs!!!
        XMultiVectorPtr mVPhi = MultiVectorFactory<SC,LO,GO,NO>::Build(localMap,coarseMap->getNodeNumElements());
        XMultiVectorPtr mVtmp = MultiVectorFactory<SC,LO,GO,NO>::Build(kII->getRowMap(),coarseMap->getNodeNumElements());
        XMultiVectorPtr mVPhiI = MultiVectorFactory<SC,LO,GO,NO>::Build(kII->getRowMap(),coarseMap->getNodeNumElements());

        //Build mVPhiGamma
        XMultiVectorPtr mVPhiGamma = MultiVectorFactory<SC,LO,GO,NO>::Build(kIGamma->getDomainMap(),coarseMap->getNodeNumElements());
        if (AssembledInterfaceCoarseSpace_->hasAssembledBasis()) {
            for (UN i=0; i<AssembledInterfaceCoarseSpace_->getAssembledBasis()->getNumVectors(); i++) {
                ConstSCVecPtr AssembledInterfaceCoarseSpaceData = AssembledInterfaceCoarseSpace_->getAssembledBasis()->getData(i);
                for (UN j=0; j<AssembledInterfaceCoarseSpace_->getAssembledBasis()->getLocalLength(); j++) {
                    mVPhiGamma->replaceLocalValue(j,i,AssembledInterfaceCoarseSpaceData[j]);
                    mVPhi->replaceLocalValue(indicesGammaDofsAll[j],i,AssembledInterfaceCoarseSpaceData[j]);
                }
            }
        }
//        LO jj=0;
//        LO kk=0;
//        UNVec numLocalBlockRows(NumberOfBlocks_);
//        for (UN i=0; i<NumberOfBlocks_; i++) {
//            UN j = 0;
//            UN k = 0;
//            if (InterfaceCoarseSpaces_[i]->hasAssembledBasis()) {
//                numLocalBlockRows[i] = InterfaceCoarseSpaces_[i]->getAssembledBasis()->getNumVectors();
//                for (j=0; j<numLocalBlockRows[i]; j++) {
//                    for (k=0; k<InterfaceCoarseSpaces_[i]->getAssembledBasis()->getLocalLength(); k++) {
//                        mVPhiGamma->replaceLocalValue(k+kk,j+jj,InterfaceCoarseSpaces_[i]->getAssembledBasis()->getData(j)[k]);
//                        mVPhi->replaceLocalValue(indicesGammaDofsAll[k+kk],j+jj,InterfaceCoarseSpaces_[i]->getAssembledBasis()->getData(j)[k]);
//                    }
//                }
//            } else { // Das ist für den Fall, dass keine Basisfunktionen für einen Block gebaut werden sollen
//                k=GammaDofs_[i].size();
//            }
//            jj += j;
//            kk += k;
//        }
        // RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(std::cout)); this->Phi_->describe(*fancy,VERB_EXTREME);
        // Hier Multiplikation kIGamma*PhiGamma
        kIGamma->apply(*mVPhiGamma,*mVtmp);

        mVtmp->scale(-ScalarTraits<SC>::one());

        // Jetzt der solver für kII
        if (indicesIDofsAll.size()>0) {
            ExtensionSolver_.reset(new SubdomainSolver<SC,LO,GO,NO>(kII,sublist(this->ParameterList_,"ExtensionSolver")));
            ExtensionSolver_->initialize();
            ExtensionSolver_->compute();
            ExtensionSolver_->apply(*mVtmp,*mVPhiI);
        }

        GOVec priorIndex(NumberOfBlocks_,0);
        GOVec postIndex(NumberOfBlocks_,0);

        GOVec2D excludeCols(NumberOfBlocks_);

        LOVec bound(NumberOfBlocks_+1,(LO)0);
        for (UN i=0; i<NumberOfBlocks_; i++) {
            bound[i+1] = bound[i] + this->IDofs_[i].size();
        }

        LO itmp = 0;
        ConstUNVecView numLocalBlockRows = AssembledInterfaceCoarseSpace_->getLocalSubspaceSizes();
        FROSCH_ASSERT(numLocalBlockRows.size()==NumberOfBlocks_,"FROSch::HarmonicCoarseOperator : ERROR: numLocalBlockRows.size()!=NumberOfBlocks_");
        for (UN i=0; i<NumberOfBlocks_; i++) {
            std::stringstream blockIdStringstream;
            blockIdStringstream << i+1;
            std::string blockIdString = blockIdStringstream.str();
            RCP<ParameterList> coarseSpaceList = sublist(sublist(this->ParameterList_,"Blocks"),blockIdString.c_str());
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
                ConstSCVecPtr mVPhiIData = mVPhiI->getData(itmp);
                for (UN ii=0; ii<extensionBlocks.size(); ii++) {
                    for (LO k=bound[extensionBlocks[ii]]; k<bound[extensionBlocks[ii]+1]; k++) {
                        mVPhi->replaceLocalValue(indicesIDofsAll[k],itmp,mVPhiIData[k]);
                    }
                }
                itmp++;
            }
        }
        return mVPhi;
    }

}

#endif
