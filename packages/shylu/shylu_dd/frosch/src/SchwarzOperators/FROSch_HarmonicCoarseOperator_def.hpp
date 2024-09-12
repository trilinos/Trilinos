// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_HARMONICCOARSEOPERATOR_DEF_HPP
#define _FROSCH_HARMONICCOARSEOPERATOR_DEF_HPP

#include <FROSch_HarmonicCoarseOperator_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    HarmonicCoarseOperator<SC,LO,GO,NO>::HarmonicCoarseOperator(ConstXMatrixPtr k,
                                                                ParameterListPtr parameterList) :
    CoarseOperator<SC,LO,GO,NO> (k,parameterList),
    AssembledInterfaceCoarseSpace_ (new CoarseSpace<SC,LO,GO,NO>(this->MpiComm_,this->SerialComm_))
    {
        FROSCH_DETAILTIMER_START_LEVELID(harmonicCoarseOperatorTime,"HarmonicCoarseOperator::HarmonicCoarseOperator");
    }

    template <class SC,class LO,class GO,class NO>
    typename HarmonicCoarseOperator<SC,LO,GO,NO>::ConstXMapPtr HarmonicCoarseOperator<SC,LO,GO,NO>::computeCoarseSpace(CoarseSpacePtr coarseSpace)
    {
        FROSCH_DETAILTIMER_START_LEVELID(computeCoarseSpaceTime,"HarmonicCoarseOperator::computeCoarseSpace");

        // Build local saddle point problem
        ConstXMapPtr repeatedMap;
        ConstXMatrixPtr repeatedMatrix;
        if (this->coarseExtractLocalSubdomainMatrix_Symbolic_Done_) {
            ExtractLocalSubdomainMatrix_Compute(this->coarseScatter_,
                                                this->K_.getConst(),
                                                this->coarseSubdomainMatrix_,
                                                this->coarseLocalSubdomainMatrix_);
            repeatedMap = this->coarseSubdomainMatrix_->getRowMap();
            repeatedMatrix = this->coarseLocalSubdomainMatrix_;
        } else {
            repeatedMap = AssembleSubdomainMap(NumberOfBlocks_,DofsMaps_,DofsPerNode_);
            repeatedMatrix = ExtractLocalSubdomainMatrix(this->K_.getConst(),repeatedMap.getConst());
        }

        // Remove coupling blocks
        if (this->ParameterList_->get("Extensions: Remove Coupling",false)) {
            FROSCH_ASSERT(this->ParameterList_->isParameter("Extensions: Coupling IDs to Remove"),"FROSch::HarmonicCoarseOperator: Coupling IDs to remove are not specified (\"Extensions: Coupling IDs to Remove\").");
            //FROSCH_ASSERT(this->ParameterList_->isType<decltype(couplingIDsToRemove)>("Extensions: Coupling IDs to Remove"),"FROSch::HarmonicCoarseOperator: \"Extensions: Coupling IDs to Remove\" is not of type Teuchos::TwoDArray<int>.");

            TwoDArray<int> couplingIDsToRemove = this->ParameterList_->get("Extensions: Coupling IDs to Remove",TwoDArray<int>(0,0,0));
            repeatedMatrix = removeCouplingBetweenDofs(repeatedMatrix,repeatedMap,couplingIDsToRemove);
        }

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

        BuildSubmatrices(repeatedMatrix.getConst(),indicesIDofsAll(),kII,kIGamma,kGammaI,kGammaGamma);

        //Detect linear dependencies
        if (!this->ParameterList_->get("Skip DetectLinearDependencies",false)) {
            LOVecPtr linearDependentVectors = detectLinearDependencies(indicesGammaDofsAll(),this->K_->getRowMap(),this->K_->getRangeMap(),repeatedMap,this->ParameterList_->get("Phi: Dropping Threshold",1.e-8),this->ParameterList_->get("Phi: Orthogonalization Threshold",1.e-12));
            // cout << this->MpiComm_->getRank() << " " << linearDependentVectors.size() << endl;
            AssembledInterfaceCoarseSpace_->zeroOutBasisVectors(linearDependentVectors());
        }

        // Build the saddle point harmonic extensions
        XMultiVectorPtr localCoarseSpaceBasis;
        if (AssembledInterfaceCoarseSpace_->getBasisMap()->getLocalNumElements()) {
            localCoarseSpaceBasis = computeExtensions(repeatedMatrix->getRowMap(),indicesGammaDofsAll(),indicesIDofsAll(),kII,kIGamma);

            coarseSpace->addSubspace(AssembledInterfaceCoarseSpace_->getBasisMap(),AssembledInterfaceCoarseSpace_->getBasisMapUnique(),localCoarseSpaceBasis);
        } else {
            FROSCH_NOTIFICATION("FROSch::HarmonicCoarseOperator",this->Verbose_,"The Coarse Space is empty. No extensions are computed.");
        }

        return repeatedMap;
    }

    template <class SC,class LO,class GO,class NO>
    int HarmonicCoarseOperator<SC,LO,GO,NO>::assembleInterfaceCoarseSpace()
    {
        FROSCH_DETAILTIMER_START_LEVELID(assembleInterfaceCoarseSpaceTime,"HarmonicCoarseOperator::assembleInterfaceCoarseSpace");
        LO ii=0;
        for (UN i=0; i<NumberOfBlocks_; i++) {
            if (!InterfaceCoarseSpaces_[i].is_null()) {
                if (InterfaceCoarseSpaces_[i]->hasBasisMap()) {
                    FROSCH_ASSERT(InterfaceCoarseSpaces_[i]->hasBasisMapUnique(),"FROSch::HarmonicCoarseOperator: !InterfaceCoarseSpaces_[i]->hasAssembledBasis()");
                    this->AssembledInterfaceCoarseSpace_->addSubspace(InterfaceCoarseSpaces_[i]->getBasisMap(),InterfaceCoarseSpaces_[i]->getBasisMapUnique(),InterfaceCoarseSpaces_[i]->getAssembledBasis(),ii);
                    ii += InterfaceCoarseSpaces_[i]->getAssembledBasis()->getLocalLength();
                }
                InterfaceCoarseSpaces_[i].reset();
            }
        }
        return this->AssembledInterfaceCoarseSpace_->assembleCoarseSpace();
    }

    template <class SC,class LO,class GO,class NO>
    int HarmonicCoarseOperator<SC,LO,GO,NO>::addZeroCoarseSpaceBlock(ConstXMapPtr dofsMap)
    {
        FROSCH_DETAILTIMER_START_LEVELID(addZeroCoarseSpaceBlockTime,"HarmonicCoarseOperator::addZeroCoarseSpaceBlock");
        /////
        int blockId = NumberOfBlocks_-1;

        // Process the parameter list
        stringstream blockIdStringstream;
        blockIdStringstream << blockId+1;
        string blockIdString = blockIdStringstream.str();
        RCP<ParameterList> coarseSpaceList = sublist(sublist(this->ParameterList_,"Blocks"),blockIdString.c_str());
        bool useForCoarseSpace = coarseSpaceList->get("Use For Coarse Space",true);

	if (useForCoarseSpace) {
            // Das könnte man noch ändern
            NumberOfBlocks_++;

            GammaDofs_[blockId] = LOVecPtr(0);

            XMultiVectorPtr mVPhiGamma;
            XMapPtr blockCoarseMap;
            if (useForCoarseSpace) {
                InterfaceCoarseSpaces_[blockId].reset(new CoarseSpace<SC,LO,GO,NO>(this->MpiComm_,this->SerialComm_));

                //Epetra_SerialComm serialComm;
                XMapPtr serialGammaMap = MapFactory<LO,GO,NO>::Build(dofsMap->lib(),dofsMap->getLocalNumElements(),0,this->SerialComm_);
                mVPhiGamma = MultiVectorFactory<LO,GO,NO>::Build(serialGammaMap,dofsMap->getLocalNumElements());
            }

            for (int i=0; i<dofsMap->getLocalNumElements(); i++) {
                GammaDofs_[blockId]->push_back(i);

                if (useForCoarseSpace) {
                    mVPhiGamma->replaceLocalValue(i,i,ScalarTraits<SC>::one());
                }
            }

            GammaDofs_->resize(GammaDofs_.size()+1);
            IDofs_->resize(IDofs_.size()+1);
            InterfaceCoarseSpaces_->resize(InterfaceCoarseSpaces_.size()+1);
            DofsMaps_->resize(DofsMaps_.size()+1);
            DofsPerNode_->resize(DofsPerNode_.size()+1);

            IDofs_[blockId] = LOVecPtr(0);

            if (useForCoarseSpace) {
                blockCoarseMap = MapFactory<LO,GO,NO>::Build(dofsMap->lib(),-1,GammaDofs_[blockId](),0,this->MpiComm_);

                InterfaceCoarseSpaces_[blockId]->addSubspace(blockCoarseMap,mVPhiGamma);
                InterfaceCoarseSpaces_[blockId]->assembleCoarseSpace();
            }

            DofsMaps_[blockId] = XMapPtrVecPtr(0);
            DofsMaps_[blockId].push_back(dofsMap);

            DofsPerNode_[blockId] = 1;
        }
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int HarmonicCoarseOperator<SC,LO,GO,NO>::computeVolumeFunctions(UN blockId,
                                                                    UN dimension,
                                                                    ConstXMapPtr nodesMap,
                                                                    ConstXMultiVectorPtr nodeList,
                                                                    EntitySetConstPtr interior)
    {
        FROSCH_DETAILTIMER_START_LEVELID(computeVolumeFunctionsTime,"HarmonicCoarseOperator::computeVolumeFunctions");
        // Process the parameter list
        stringstream blockIdStringstream;
        blockIdStringstream << blockId+1;
        string blockIdString = blockIdStringstream.str();
        RCP<ParameterList> coarseSpaceList = sublist(sublist(this->ParameterList_,"Blocks"),blockIdString.c_str());

        bool useForCoarseSpace = coarseSpaceList->get("Use For Coarse Space",true);
        bool useRotations = coarseSpaceList->get("Rotations",true);
        if (useRotations && nodeList.is_null()) {
            useRotations = false;
            FROSCH_WARNING("FROSch::HarmonicCoarseOperator",this->Verbose_,"Rotations cannot be used since nodeList.is_null().");
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

            if (this->Verbose_) {
                cout
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << setw(89) << "-----------------------------------------------------------------------------------------"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| "
                << left << setw(74) << "> RGDSW coarse space " << right << setw(8) << "(Level " << setw(2) << this->LevelID_ << ")"
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << setw(89) << "========================================================================================="
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| " << left << setw(19) << "Volumes " << " | " << setw(19) << "Translations " << right
                << " | " << setw(41) << boolalpha << useForCoarseSpace << noboolalpha
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| " << left << setw(19) << "Volumes " << " | " << setw(19) << "Rotations " << right
                << " | " << setw(41) << boolalpha << useRotations << noboolalpha
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << setw(89) << "-----------------------------------------------------------------------------------------"
                << endl;
            }
        }
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    typename HarmonicCoarseOperator<SC,LO,GO,NO>::XMultiVectorPtrVecPtr HarmonicCoarseOperator<SC,LO,GO,NO>::computeTranslations(UN blockId,
                                                                                                                                 EntitySetConstPtr entitySet)
    {
        FROSCH_DETAILTIMER_START_LEVELID(computeTranslationsTime,"HarmonicCoarseOperator::computeTranslations");
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
    typename HarmonicCoarseOperator<SC,LO,GO,NO>::XMapPtr HarmonicCoarseOperator<SC,LO,GO,NO>::assembleCoarseMap()
    {
        FROSCH_DETAILTIMER_START_LEVELID(assembleCoarseMapTime,"HarmonicCoarseOperator::assembleCoarseMap");
        GOVec mapVector(0);
        GO tmp = 0;
        for (UN i=0; i<NumberOfBlocks_; i++) {
            if (!InterfaceCoarseSpaces_[i].is_null()) {
                if (InterfaceCoarseSpaces_[i]->hasBasisMap()) {
                    for (UN j=0; j<InterfaceCoarseSpaces_[i]->getBasisMap()->getLocalNumElements(); j++) {
                        mapVector.push_back(InterfaceCoarseSpaces_[i]->getBasisMap()->getGlobalElement(j)+tmp);
                    }
                    if (InterfaceCoarseSpaces_[i]->getBasisMap()->getMaxAllGlobalIndex()>=0) {
                        tmp += InterfaceCoarseSpaces_[i]->getBasisMap()->getMaxAllGlobalIndex()+1;
                    }
                }
            }
        }
        return MapFactory<LO,GO,NO>::Build(DofsMaps_[0][0]->lib(),-1,mapVector(),0,this->MpiComm_);
    }

    template <class SC,class LO,class GO,class NO>
    typename HarmonicCoarseOperator<SC,LO,GO,NO>::XMultiVectorPtrVecPtr HarmonicCoarseOperator<SC,LO,GO,NO>::computeRotations(UN blockId,
                                                                                                                              UN dimension,
                                                                                                                              ConstXMultiVectorPtr nodeList,
                                                                                                                              EntitySetConstPtr entitySet,
                                                                                                                              UN discardRotations)
    {
        FROSCH_DETAILTIMER_START_LEVELID(computeRotationsTime,"HarmonicCoarseOperator::computeRotations");
        FROSCH_ASSERT(nodeList->getNumVectors()==dimension,"FROSch::HarmonicCoarseOperator: Dimension of the nodeList is wrong.");
        FROSCH_ASSERT(dimension==this->DofsPerNode_[blockId],"FROSch::HarmonicCoarseOperator: Dimension!=this->DofsPerNode_[blockId]");

        UN rotationsPerEntity = 0;
        switch (dimension) {
            case 1:
                return null;
            case 2:
                rotationsPerEntity = 1;
                break;
            case 3:
                rotationsPerEntity = 3;
                break;
            default:
                FROSCH_ASSERT(false,"FROSch::HarmonicCoarseOperator: The dimension is neither 2 nor 3!");
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
                    if (fabs(err)<1.0e-12) {
                        FROSCH_ASSERT(false,"FROSch::HarmonicCoarseOperator: In 2D, no rotation can be constant!");
                        rotations[0]->getVectorNonConst(i)->scale(ScalarTraits<SC>::zero());
                        numZeroRotations++;
                    }
                    break;
                case 3:
                    for (UN j=0; j<3; j++) {
                        err = errx[j]+erry[j]+errz[j];
                        err = ScalarTraits<SC>::squareroot(err);
                        if (fabs(err)<1.0e-12) {
                            rotations[j]->getVectorNonConst(i)->scale(ScalarTraits<SC>::zero());
                            numZeroRotations++;
                        }
                    }
                    break;
                default:
                    FROSCH_ASSERT(false,"FROSch::HarmonicCoarseOperator: The dimension is neither 1 nor 2 nor 3!");
            }
            // If necessary, discard additional rotations
            UN rotationsToDiscard = discardRotations - numZeroRotations;
            if (rotationsToDiscard<0) {
                FROSCH_WARNING("FROSch::HarmonicCoarseOperator",this->Verbose_,"More rotations have been discarded than expected.");
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
    typename HarmonicCoarseOperator<SC,LO,GO,NO>::LOVecPtr HarmonicCoarseOperator<SC,LO,GO,NO>::detectLinearDependencies(GOVecView indicesGammaDofsAll,
                                                                                                                         ConstXMapPtr rowMap,
                                                                                                                         ConstXMapPtr rangeMap,
                                                                                                                         ConstXMapPtr repeatedMap,
                                                                                                                         SC tresholdDropping,
                                                                                                                         SC tresholdOrthogonalization)
    {
        FROSCH_DETAILTIMER_START_LEVELID(detectLinearDependenciesTime,"HarmonicCoarseOperator::detectLinearDependencies");
        LOVecPtr linearDependentVectors(AssembledInterfaceCoarseSpace_->getBasisMap()->getLocalNumElements()); //if (this->Verbose_) cout << AssembledInterfaceCoarseSpace_->getAssembledBasis()->getNumVectors() << " " << AssembledInterfaceCoarseSpace_->getAssembledBasis()->getLocalLength() << " " << indicesGammaDofsAll.size() << endl;

        auto asembledBasis = AssembledInterfaceCoarseSpace_->getAssembledBasis();
        auto basisMap = AssembledInterfaceCoarseSpace_->getBasisMap();
        auto basisMapUnique = AssembledInterfaceCoarseSpace_->getBasisMapUnique();
        UN numMVRows = asembledBasis->getLocalLength();
        UN numCols = asembledBasis->getNumVectors();
        if (numMVRows > 0 && numCols > 0) {

            XMatrixPtr phiGamma;
#if defined(HAVE_XPETRA_KOKKOS_REFACTOR) && defined(HAVE_XPETRA_TPETRA)
            if (rowMap->lib() == UseTpetra)
            {
                using crsmat_type  = typename Matrix<SC,LO,GO,NO>::local_matrix_type;
                using graph_type   = typename crsmat_type::StaticCrsGraphType;
                using rowptr_type  = typename graph_type::row_map_type::non_const_type;
                using indices_type = typename graph_type::entries_type::non_const_type;
                using values_type  = typename crsmat_type::values_type::non_const_type;
                using local_mv_type = typename MultiVector<SC,LO,GO,NO>::dual_view_type::t_dev_const_um;
                using local_map_type = typename Map<LO,GO,NO>::local_map_type;

                UN numRows = rowMap->getLocalNumElements();

                auto localRowMap = rowMap->getLocalMap();
                auto localRepeatedMap = repeatedMap->getLocalMap();
                auto localBasisMap = basisMap->getLocalMap();

                auto localMVBasis = asembledBasis->getDeviceLocalView(Xpetra::Access::ReadOnly);

                // Array for scaling the columns of PhiGamma (1/norm(PhiGamma(:,i)))
                using execution_space = typename Map<LO,GO,NO>::local_map_type::execution_space;
                using SCView = Kokkos::View<SC*, execution_space>;
                SCView scale(Kokkos::ViewAllocateWithoutInitializing("FROSch::HarmonicCoarseOperator::detectLinearDependencies::scale"), numCols);
                detectLinearDependenciesFunctor<GOVecView, SCView, local_map_type, local_mv_type, rowptr_type, indices_type, values_type>
                    scale_functor(numMVRows, numCols, scale, localMVBasis);
                Kokkos::RangePolicy<ScaleTag, execution_space> policy_scale (0, numCols);
                Kokkos::parallel_for(
                    "FROSch_HarmonicCoarseOperator::detectLinearDependencies::scale", policy_scale, scale_functor);

                // count nnz
                rowptr_type Rowptr (Kokkos::ViewAllocateWithoutInitializing("Rowptr"), numRows+1);
                Kokkos::deep_copy(Rowptr, 0);
                detectLinearDependenciesFunctor<GOVecView, SCView, local_map_type, local_mv_type, rowptr_type, indices_type, values_type>
                    count_functor(numMVRows, numCols, localMVBasis, tresholdDropping, indicesGammaDofsAll,
                                  localRowMap, localRepeatedMap, Rowptr);
                Kokkos::RangePolicy<CountNnzTag, execution_space> policy_count (0, numMVRows);
                Kokkos::parallel_for(
                    "FROSch_HarmonicCoarseOperator::detectLinearDependencies::countNnz", policy_count, count_functor);
                Kokkos::fence();

                // get Nnz
                UN nnz = 0;
                Kokkos::RangePolicy<TotalNnzTag, execution_space> policy_nnz (0, 1+numRows);
                Kokkos::parallel_reduce(
                    "FROSch_HarmonicCoarseOperator::detectLinearDependencies::reduceNnz", policy_nnz, count_functor, nnz);
                Kokkos::fence();

                // make it into offsets
                KokkosKernels::Impl::kk_inclusive_parallel_prefix_sum<rowptr_type, execution_space>
                    (1+numRows, Rowptr);
                Kokkos::fence();

                // allocate local matrix
                indices_type Indices (Kokkos::ViewAllocateWithoutInitializing("Indices"), nnz);
                values_type  Values  (Kokkos::ViewAllocateWithoutInitializing("Values"),  nnz);

                // fill in all the nz entries
                detectLinearDependenciesFunctor<GOVecView, SCView, local_map_type, local_mv_type, rowptr_type, indices_type, values_type>
                    fill_functor(numMVRows, numCols, scale, localMVBasis, tresholdDropping, indicesGammaDofsAll,
                                 localRowMap, localRepeatedMap, Rowptr, Indices, Values);
                Kokkos::RangePolicy<FillNzEntriesTag, execution_space> policy_fill (0, numMVRows);
                Kokkos::parallel_for(
                    "FROSch_HarmonicCoarseOperator::detectLinearDependencies::fillNz", policy_fill, fill_functor);
                Kokkos::fence();

                Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp (new Teuchos::ParameterList());
                params->set("sorted", false);
                params->set("No Nonlocal Changes", true);

                graph_type crsgraph (Indices, Rowptr);
                crsmat_type crsmat = crsmat_type ("CrsMatrix", numRows, Values, crsgraph);

                phiGamma = MatrixFactory<SC,LO,GO,NO>::Build(crsmat, rowMap, basisMap, basisMapUnique, rangeMap,
                                                             params);
            } else
#endif
            {
                // Array for scaling the columns of PhiGamma (1/norm(PhiGamma(:,i)))
                SCVec scale(numCols, 0.0);
                for (UN j = 0; j < numCols; j++) {
                    ConstSCVecPtr assembledInterfaceCoarseSpaceData = asembledBasis->getData(j);
                    for (UN i = 0; i < numMVRows; i++) {
                        scale[j] += assembledInterfaceCoarseSpaceData[i]*assembledInterfaceCoarseSpaceData[i];
                    }
                    scale[j] = 1.0/sqrt(scale[j]);
                }

                //Construct matrix phiGamma
                phiGamma = MatrixFactory<SC,LO,GO,NO>::Build(rowMap,basisMap->getLocalNumElements());
                LO iD;
                SC valueTmp;
                for (UN i=0; i<asembledBasis->getLocalLength(); i++) {
                    iD = repeatedMap->getGlobalElement(indicesGammaDofsAll[i]);
                    if (rowMap->getLocalElement(iD)!=-1) { // This should prevent duplicate entries on the interface
                        GOVec indices;
                        SCVec values;
                        for (UN j=0; j<asembledBasis->getNumVectors(); j++) {
                            valueTmp=asembledBasis->getData(j)[i];
                            if (fabs(valueTmp)>tresholdDropping) {
                                indices.push_back(basisMap->getGlobalElement(j));
                                values.push_back(valueTmp*scale[j]);
                            }
                        }
                        phiGamma->insertGlobalValues(iD,indices(),values());
                    }
                }
                RCP<ParameterList> fillCompleteParams(new ParameterList);
                fillCompleteParams->set("No Nonlocal Changes", true);
                phiGamma->fillComplete(basisMapUnique,rangeMap,fillCompleteParams);
            }

            //Compute Phi^T * Phi
            RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout));
            XMatrixPtr phiTPhi = MatrixMatrix<SC,LO,GO,NO>::Multiply(*phiGamma,true,*phiGamma,false,*fancy); //phiGamma->describe(*fancy,VERB_EXTREME); phiTPhi->describe(*fancy,VERB_EXTREME); //AssembledInterfaceCoarseSpace_->getBasisMap()->describe(*fancy,VERB_EXTREME);

            // Extract local part of the matrix
            ConstXMatrixPtr repeatedPhiTPhi = ExtractLocalSubdomainMatrix(phiTPhi.getConst(),AssembledInterfaceCoarseSpace_->getBasisMap());
            //Communicate matrix to the repeated map
            // repeatedPhiTPhi = MatrixFactory<SC,LO,GO,NO>::Build(AssembledInterfaceCoarseSpace_->getBasisMap());
            // XExportPtr exporter = ExportFactory<LO,GO,NO>::Build(rowMap,repeatedMap);
            // repeatedPhiTPhi->doExport(*phiTPhi,*exporter,INSERT);

            UN numRows = repeatedPhiTPhi->getRowMap()->getLocalNumElements();
            TSerialDenseMatrixPtr denseRepeatedPhiTPhi(new SerialDenseMatrix<LO,SC>(numRows,repeatedPhiTPhi->getColMap()->getLocalNumElements()));
            for (UN i=0; i<numRows; i++) {
                ConstLOVecView indices;
                ConstSCVecView values;
                repeatedPhiTPhi->getLocalRowView(i,indices,values);
                for (UN j=0; j<indices.size(); j++) {
                    (*denseRepeatedPhiTPhi)(i,indices[j]) = values[j];
                }
            }
            // if (this->MpiComm_->getRank()==3) {
            //     for (LO i=0; i<denseRepeatedPhiTPhi->numRows(); i++) {
            //         for (LO j=0; j<denseRepeatedPhiTPhi->numCols(); j++) {
            //             cout << (*denseRepeatedPhiTPhi)(i,j) << " ";
            //         }
            //         cout << endl;
            //     }
            // }

            //Compute local QR factorization
            TSerialQRDenseSolverPtr qRSolver(new SerialQRDenseSolver<LO,SC>());
            qRSolver->setMatrix(denseRepeatedPhiTPhi);
            qRSolver->factor();
            qRSolver->formQ();
            qRSolver->formR();

            //Find rows of R approx. zero
            TSerialDenseMatrixPtr r = qRSolver->getR();
            LO tmp = 0;
            for (LO i=0; i<r->numRows(); i++) {
                // SC normRow = 0.0;
                // for (LO j=0; j<r->numCols(); j++) {
                //     normRow += (*r)(i,j)*(*r)(i,j);
                // }
                // if (sqrt(normRow)<treshold) {
                //     //cout << this->MpiComm_->getRank() << " " << i << " " << AssembledInterfaceCoarseSpace_->getBasisMap()->getGlobalElement(i) << " " << sqrt(normRow) << std::endl;
                //     linearDependentVectors[tmp] = i;
                //     tmp++;
                // }
                if (fabs((*r)(i,i))<tresholdOrthogonalization) {
                    //cout << this->MpiComm_->getRank() << " " << i << " " << AssembledInterfaceCoarseSpace_->getBasisMap()->getGlobalElement(i) << " " << sqrt(normRow) << std::endl;
                    linearDependentVectors[tmp] = i;
                    tmp++;
                }
            }
            linearDependentVectors.resize(tmp);
        }

        FROSCH_DETAILTIMER_START_LEVELID(printStatisticsTime,"print statistics");
        // Statistics on linear dependencies
        GO global = AssembledInterfaceCoarseSpace_->getBasisMap()->getMaxAllGlobalIndex();
        if (AssembledInterfaceCoarseSpace_->getBasisMap()->lib()==UseEpetra || AssembledInterfaceCoarseSpace_->getBasisMap()->getGlobalNumElements()>0) {
            global += 1;
        }
        LOVec localVec(3);
        LOVec sumVec(3);
        SCVec avgVec(3);
        LOVec minVec(3);
        LOVec maxVec(3);

        LO numLocalBasisFunctions = AssembledInterfaceCoarseSpace_->getBasisMap()->getLocalNumElements();
        LO numLocalLinearDependencies = linearDependentVectors.size();
        LO numLocalLinearIndependencies = numLocalBasisFunctions-numLocalLinearDependencies;

        sumVec[0] = AssembledInterfaceCoarseSpace_->getBasisMap()->getGlobalNumElements();
        avgVec[0] = max(sumVec[0]/double(this->MpiComm_->getSize()),0.0);
        reduceAll(*this->MpiComm_,REDUCE_MIN,numLocalBasisFunctions,ptr(&minVec[0]));
        reduceAll(*this->MpiComm_,REDUCE_MAX,numLocalBasisFunctions,ptr(&maxVec[0]));

        reduceAll(*this->MpiComm_,REDUCE_SUM,numLocalLinearDependencies,ptr(&sumVec[1]));
        avgVec[1] = max(sumVec[1]/double(this->MpiComm_->getSize()),0.0);
        reduceAll(*this->MpiComm_,REDUCE_MIN,numLocalLinearDependencies,ptr(&minVec[1]));
        reduceAll(*this->MpiComm_,REDUCE_MAX,numLocalLinearDependencies,ptr(&maxVec[1]));

        reduceAll(*this->MpiComm_,REDUCE_SUM,numLocalLinearIndependencies,ptr(&sumVec[2]));
        avgVec[2] = max(sumVec[2]/double(this->MpiComm_->getSize()),0.0);
        reduceAll(*this->MpiComm_,REDUCE_MIN,numLocalLinearIndependencies,ptr(&minVec[2]));
        reduceAll(*this->MpiComm_,REDUCE_MAX,numLocalLinearIndependencies,ptr(&maxVec[2]));

        if (this->Verbose_) {
            cout
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << setw(89) << "-----------------------------------------------------------------------------------------"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << "| "
            << left << setw(74) << "> Local Linear Dependencies of Coarse Basis Functions Statistics " << right << setw(8) << "(Level " << setw(2) << this->LevelID_ << ")"
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
            << "| " << left << setw(20) << "Basis functions" << right
            << " | " << setw(10) << global
            << " | " << setw(10) << setprecision(5) << avgVec[0]
            << " | " << setw(10) << minVec[0]
            << " | " << setw(10) << maxVec[0]
            << " | " << setw(10) << sumVec[0]
            << " |"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << "| " << left << setw(20) << "Dependent" << right
            << " | " << setw(10) << " "
            << " | " << setw(10) << setprecision(5) << avgVec[1]
            << " | " << setw(10) << minVec[1]
            << " | " << setw(10) << maxVec[1]
            << " | " << setw(10) << sumVec[1]
            << " |"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << "| " << left << setw(20) << "Independent" << right
            << " | " << setw(10) << " "
            << " | " << setw(10) << setprecision(5) << avgVec[2]
            << " | " << setw(10) << minVec[2]
            << " | " << setw(10) << maxVec[2]
            << " | " << setw(10) << sumVec[2]
            << " |"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << setw(89) << "-----------------------------------------------------------------------------------------"
            << endl;
        }
        FROSCH_DETAILTIMER_STOP(printStatisticsTime);

        return linearDependentVectors;
    }

    template <class SC,class LO,class GO,class NO>
    typename HarmonicCoarseOperator<SC,LO,GO,NO>::XMultiVectorPtr HarmonicCoarseOperator<SC,LO,GO,NO>::computeExtensions(ConstXMapPtr localMap,
                                                                                                                         GOVecView indicesGammaDofsAll,
                                                                                                                         GOVecView indicesIDofsAll,
                                                                                                                         XMatrixPtr kII,
                                                                                                                         XMatrixPtr kIGamma)
    {
        FROSCH_DETAILTIMER_START_LEVELID(computeExtensionsTime,"HarmonicCoarseOperator::computeExtensions");

        if (this->Verbose_) {
            cout
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << setw(89) << "-----------------------------------------------------------------------------------------"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << "| "
            << left << setw(74) << "> Extensions solver" << right << setw(8) << "(Level " << setw(2) << this->LevelID_ << ")"
            << " |"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << setw(89) << "========================================================================================="
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << "| " << left << setw(41) << "Solver type" << right
            << " | " << setw(41) << this->ParameterList_->sublist("ExtensionSolver").get("SolverType","Amesos2")
            << " |"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << "| " << left << setw(41) << "Solver" << right
            << " | " << setw(41) << this->ParameterList_->sublist("ExtensionSolver").get("Solver","Klu")
            << " |"
            << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
            << setw(89) << "-----------------------------------------------------------------------------------------"
            << endl;
        }

        //this->Phi_ = MatrixFactory<SC,LO,GO,NO>::Build(this->K_->getRangeMap(),AssembledInterfaceCoarseSpace_->getBasisMap(),AssembledInterfaceCoarseSpace_->getBasisMap()->getLocalNumElements()); // Nonzeroes abhängig von dim/dofs!!!
        XMultiVectorPtr mVPhi = MultiVectorFactory<SC,LO,GO,NO>::Build(localMap,AssembledInterfaceCoarseSpace_->getBasisMap()->getLocalNumElements());
        XMultiVectorPtr mVtmp = MultiVectorFactory<SC,LO,GO,NO>::Build(kII->getRowMap(),AssembledInterfaceCoarseSpace_->getBasisMap()->getLocalNumElements());
        XMultiVectorPtr mVPhiI = MultiVectorFactory<SC,LO,GO,NO>::Build(kII->getRowMap(),AssembledInterfaceCoarseSpace_->getBasisMap()->getLocalNumElements());

        //Build mVPhiGamma
        XMultiVectorPtr mVPhiGamma = MultiVectorFactory<SC,LO,GO,NO>::Build(kIGamma->getDomainMap(),AssembledInterfaceCoarseSpace_->getBasisMap()->getLocalNumElements());
        if (AssembledInterfaceCoarseSpace_->hasAssembledBasis()) {
            for (UN i=0; i<AssembledInterfaceCoarseSpace_->getAssembledBasis()->getNumVectors(); i++) {
                ConstSCVecPtr assembledInterfaceCoarseSpaceData = AssembledInterfaceCoarseSpace_->getAssembledBasis()->getData(i);
                for (UN j=0; j<AssembledInterfaceCoarseSpace_->getAssembledBasis()->getLocalLength(); j++) {
                    mVPhiGamma->replaceLocalValue(j,i,assembledInterfaceCoarseSpaceData[j]);
                    mVPhi->replaceLocalValue(indicesGammaDofsAll[j],i,assembledInterfaceCoarseSpaceData[j]);
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
        // RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout)); this->Phi_->describe(*fancy,VERB_EXTREME);
        // Hier Multiplikation kIGamma*PhiGamma
        kIGamma->apply(*mVPhiGamma,*mVtmp);

        mVtmp->scale(-ScalarTraits<SC>::one());

        // Jetzt der solver für kII
        if (indicesIDofsAll.size()>0) {
            if (this->coarseExtractLocalSubdomainMatrix_Symbolic_Done_) {
                ExtensionSolver_->updateMatrix(kII, true);
            } else {
                ExtensionSolver_ = SolverFactory<SC,LO,GO,NO>::Build(kII,
                                                                     sublist(this->ParameterList_,"ExtensionSolver"),
                                                                     string("ExtensionSolver (Level ") + to_string(this->LevelID_) + string(")"));
                ExtensionSolver_->initialize();
            }
            ExtensionSolver_->compute();
            ExtensionSolver_->apply(*mVtmp,*mVPhiI);
        }

        GOVec2D excludeCols(NumberOfBlocks_);

        LOVec bound(NumberOfBlocks_+1,(LO)0);
        for (UN i=0; i<NumberOfBlocks_; i++) {
            bound[i+1] = bound[i] + this->IDofs_[i].size();
        }

        LO itmp = 0;
        ConstUNVecView numLocalBlockColumns = AssembledInterfaceCoarseSpace_->getLocalSubspaceSizes();
        FROSCH_ASSERT(numLocalBlockColumns.size()==NumberOfBlocks_,"FROSch::HarmonicCoarseOperator: numLocalBlockColumns.size()!=NumberOfBlocks_("+to_string(numLocalBlockColumns.size())+", "+to_string(NumberOfBlocks_)+") ");
        for (UN i=0; i<NumberOfBlocks_; i++) {
            stringstream blockIdStringstream;
            blockIdStringstream << i+1;
            string blockIdString = blockIdStringstream.str();
            RCP<ParameterList> coarseSpaceList = sublist(sublist(this->ParameterList_,"Blocks"),blockIdString.c_str());
            string excludeBlocksString = coarseSpaceList->get("Exclude","0");
            UNVec excludeBlocks = FROSch::GetIndicesFromString(excludeBlocksString,(UN)0);
            sortunique(excludeBlocks);
            for (UN j=0; j<excludeBlocks.size(); j++) {
                excludeBlocks[j] -= 1;
            }
            UNVec extensionBlocks(0);
            for (UN j=0; j<NumberOfBlocks_; j++) {
                typename UNVec::iterator it = find(excludeBlocks.begin(),excludeBlocks.end(),j);
                if (it == excludeBlocks.end()) {
                    extensionBlocks.push_back(j);
                }
            }

            #if defined(HAVE_XPETRA_TPETRA)
            if (mVPhi->getMap()->lib() == UseTpetra) {
                using XMap            = typename SchwarzOperator<SC,LO,GO,NO>::XMap;
                using execution_space = typename XMap::local_map_type::execution_space;

                // explicitly copying indicesIDofsAll in Teuchos::Array to Kokkos::View on "device"
                using GOIndViewHost = Kokkos::View<GO*, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
                using GOIndView     = Kokkos::View<GO*, execution_space>;
                UN numIndices = indicesIDofsAll.size();
                GOIndViewHost indicesIDofsAllHostData(indicesIDofsAll.getRawPtr(), numIndices);
                GOIndView     indicesIDofsAllData (Kokkos::ViewAllocateWithoutInitializing("indicesIDofsAllData"), numIndices);
                Kokkos::deep_copy(indicesIDofsAllData, indicesIDofsAllHostData);

                using xTMVector    = Xpetra::TpetraMultiVector<SC,LO,GO,NO>;
                using TMVector     = Tpetra::MultiVector<SC,LO,GO,NO>;
                using SCView       = typename TMVector::dual_view_type::t_dev;
                using ConstSCView  = typename TMVector::dual_view_type::t_dev::const_type;
                // Xpetra wrapper for Tpetra MV
                auto mVPhiIXTpetraMVector = rcp_dynamic_cast<const xTMVector>(mVPhiI, true);
                auto mVPhiXTpetraMVector = rcp_dynamic_cast<       xTMVector>(mVPhi, true);
                // Tpetra MV
                auto mVPhiITpetraMVector = mVPhiIXTpetraMVector->getTpetra_MultiVector();
                auto mVPhiTpetraMVector = mVPhiXTpetraMVector->getTpetra_MultiVector();
                // Kokkos-Kernels Views
                auto mvPhiIView = mVPhiITpetraMVector->getLocalViewDevice(Tpetra::Access::ReadOnly);
                auto mvPhiView = mVPhiTpetraMVector->getLocalViewDevice(Tpetra::Access::ReadWrite);
                auto mvPhiICols = Tpetra::getMultiVectorWhichVectors(*mVPhiITpetraMVector);
                auto mvPhiCols = Tpetra::getMultiVectorWhichVectors(*mVPhiTpetraMVector);
                for (UN j=0; j<numLocalBlockColumns[i]; j++) {
                    int col_in = mVPhiITpetraMVector->isConstantStride() ? itmp : mvPhiICols[itmp];
                    int col_out = mVPhiTpetraMVector->isConstantStride() ? itmp : mvPhiCols[itmp];
                    CopyPhiViewFunctor<GOIndView, ConstSCView, SCView> functor(col_in, indicesIDofsAllData, mvPhiIView, col_out, mvPhiView);
                    for (UN ii=0; ii<extensionBlocks.size(); ii++) {
                        Kokkos::RangePolicy<execution_space> policy (bound[extensionBlocks[ii]], bound[extensionBlocks[ii]+1]);
                        Kokkos::parallel_for(
                            "FROSch_HarmonicCoarseOperator::fillPhiData", policy, functor);
                    }
                    itmp++;
                }
                Kokkos::fence();
            } else
            #endif
            {
                for (UN j=0; j<numLocalBlockColumns[i]; j++) {
                    ConstSCVecPtr mVPhiIData = mVPhiI->getData(itmp);
                    for (UN ii=0; ii<extensionBlocks.size(); ii++) {
                        for (LO k=bound[extensionBlocks[ii]]; k<bound[extensionBlocks[ii]+1]; k++) {
                            mVPhi->replaceLocalValue(indicesIDofsAll[k],itmp,mVPhiIData[k]);
                        }
                    }
                    itmp++;
                }
            }
        }
        return mVPhi;
    }

    // Should we move this to Tools?
    template <class SC,class LO,class GO,class NO>
    typename HarmonicCoarseOperator<SC,LO,GO,NO>::ConstXMatrixPtr HarmonicCoarseOperator<SC,LO,GO,NO>::removeCouplingBetweenDofs(ConstXMatrixPtr matrix,
                                                                                                                                 ConstXMapPtr map,
                                                                                                                                 TwoDArray<int> &couplingIDsToRemove)
    {
        FROSCH_DETAILTIMER_START_LEVELID(removeCouplingBetweenDofsTime,"HarmonicCoarseOperator::removeCouplingBetweenDofs");

        ConstUNVecView numLocalBlockRows = AssembledInterfaceCoarseSpace_->getLocalSubspaceSizes();
        FROSCH_ASSERT(numLocalBlockRows.size()==NumberOfBlocks_,"FROSch::HarmonicCoarseOperator: numLocalBlockRows.size()!=NumberOfBlocks_");

        // Array<RCP<MultiVector<bool,LO,GO,NO> > > mask(NumberOfBlocks_);
        // for (UN i=0; i<NumberOfBlocks_; i++) {
        //     mask[i] = MultiVectorFactory<bool,LO,GO,NO>::Build(matrix->getRowMap(),DofsPerNode_[i]);
        //     mask[i]->putScalar(false);
        // }
        Array<Array<bool> > mask(NumberOfBlocks_);
        for (UN i=0; i<NumberOfBlocks_; i++) {
            mask[i] = Array<bool>(DofsPerNode_[i]*map->getLocalNumElements(),false);
        }

        FROSCH_ASSERT(couplingIDsToRemove.getNumCols()==4,"FROSch::HarmonicCoarseOperator: couplingIDsToRemove.getNumCols()!=4");
        for (UN i=0; i<couplingIDsToRemove.getNumRows() ; i++) {
            FROSCH_ASSERT(couplingIDsToRemove[i][0]>=0,"FROSch::HarmonicCoarseOperator: couplingIDsToRemove[i][0]<0");
            FROSCH_ASSERT(couplingIDsToRemove[i][1]>=0,"FROSch::HarmonicCoarseOperator: couplingIDsToRemove[i][1]<0");
            FROSCH_ASSERT(couplingIDsToRemove[i][2]>=0,"FROSch::HarmonicCoarseOperator: couplingIDsToRemove[i][2]<0");
            FROSCH_ASSERT(couplingIDsToRemove[i][3]>=0,"FROSch::HarmonicCoarseOperator: couplingIDsToRemove[i][3]<0");

            int blockID1 = couplingIDsToRemove[i][0];
            int dofID1 = couplingIDsToRemove[i][1];

            int blockID2 = couplingIDsToRemove[i][2];
            int dofID2 = couplingIDsToRemove[i][3];

            Array<bool>& tmpMask = mask[blockID1];
            for (UN j=0; j<DofsMaps_[blockID2][dofID2]->getLocalNumElements(); j++) {
                GO globalIndex = DofsMaps_[blockID2][dofID2]->getGlobalElement(j);
                FROSCH_ASSERT(globalIndex>=0,"FROSch::HarmonicCoarseOperator: globalIndex<0");
                LO localIndex = map->getLocalElement(globalIndex);
                FROSCH_ASSERT(localIndex>=0,"FROSch::HarmonicCoarseOperator: localIndex<0");
                //mask[blockID1]->replaceLocalValue(localIndex,dofID1,true);
                tmpMask[dofID1*map->getLocalNumElements()+localIndex] = true;
            }
        }

        XMatrixPtr reducedMatrix = MatrixFactory<SC,LO,GO,NO>::Build(matrix->getRowMap(),matrix->getLocalMaxNumRowEntries());
        for (UN i=0; i<NumberOfBlocks_; i++) {
            Array<bool>& tmpMask = mask[i];
            for (UN j=0; j<DofsPerNode_[i]; j++) {
                //ConstBoolVecPtr maskData = mask[i]->getData(j);
                for (UN k=0; k<DofsMaps_[i][j]->getLocalNumElements(); k++) {
                    GO globalIndex = DofsMaps_[i][j]->getGlobalElement(k);
                    FROSCH_ASSERT(globalIndex>=0,"FROSch::HarmonicCoarseOperator: globalIndex<0");
                    LO localIndex = map->getLocalElement(globalIndex);
                    FROSCH_ASSERT(localIndex>=0,"FROSch::HarmonicCoarseOperator: localIndex<0");

                    ArrayView<const LO> indices;
                    ArrayView<const SC> values;
                    matrix->getLocalRowView(localIndex,indices,values);

                    LO size = indices.size();
                    if (size>0) {
                        Array<GO> indicesGlobal;
                        Array<SC> valuesGlobal;
                        for (LO l=0; l<size; l++) {
                            //if (!maskData[indices[l]]) {
                            if (!tmpMask[j*map->getLocalNumElements()+indices[l]]) {
                                indicesGlobal.push_back(matrix->getRowMap()->getGlobalElement(indices[l]));
                                valuesGlobal.push_back(values[l]);
                            }
                        }
                        reducedMatrix->insertGlobalValues(matrix->getRowMap()->getGlobalElement(localIndex),indicesGlobal(),valuesGlobal());
                    }
                }
            }
        }
        reducedMatrix->fillComplete();

        // for (UN i=0; i<NumberOfBlocks_; i++) {
        //     if (this->Verbose_) {RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout)); matrix->describe(*fancy,VERB_EXTREME);}
        //     if (this->Verbose_) {RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout)); reducedMatrix->describe(*fancy,VERB_EXTREME);}
        // }

        return reducedMatrix.getConst();
    }

    template <class SC,class LO,class GO, class NO>
    int HarmonicCoarseOperator<SC,LO,GO,NO>::buildElementNodeList()
    {
        //get elements belonging to one subdomain
        FROSCH_DETAILTIMER_START_LEVELID(buildElementNodeListTime,"CoarseOperator::buildElementNodeList");

        Teuchos::ArrayView<const GO> elements =  KRowMap_->getLocalElementList();
        UN maxNumElements = -1;
        UN numElementsLocal = elements.size();

        reduceAll(*this->MpiComm_,Teuchos::REDUCE_MAX,numElementsLocal,Teuchos::ptr(&maxNumElements));


        GraphPtr elemGraph = Xpetra::CrsGraphFactory<LO,GO,NO>::Build(this->MLGatheringMaps_[0],maxNumElements);
        Teuchos::ArrayView<const GO> myGlobals = this->SubdomainConnectGraph_->getRowMap()->getLocalElementList();

        for (size_t i = 0; i < this->SubdomainConnectGraph_->getRowMap()->getLocalNumElements(); i++) {
            elemGraph->insertGlobalIndices(myGlobals[i],elements);
        }
        elemGraph->fillComplete();

        GraphPtr tmpElemGraph = Xpetra::CrsGraphFactory<LO,GO,NO>::Build(this->MLGatheringMaps_[1],maxNumElements);
        GraphPtr elemSGraph;
        //communicate ElemGrapg to CoarseComm_
        tmpElemGraph->doExport(*elemGraph,*this->MLCoarseSolveExporters_[1],Xpetra::INSERT);
        UN gathered = 0;
        for (int i  = 2;i<this->MLGatheringMaps_.size();i++) {
            gathered = 1;
            tmpElemGraph->fillComplete();
            elemSGraph = tmpElemGraph;
            tmpElemGraph = Xpetra::CrsGraphFactory<LO,GO,NO>::Build(this->MLGatheringMaps_[i],maxNumElements);
            tmpElemGraph->doExport(*elemSGraph,*this->MLCoarseSolveExporters_[i],Xpetra::INSERT);
        }
        //tmpElemGraph->fillComplete();
        elemSGraph = tmpElemGraph;

        if (gathered == 0) {
            elemSGraph = tmpElemGraph;
        }
        this->ElementNodeList_ =Xpetra::CrsGraphFactory<LO,GO,NO>::Build(this->MLCoarseMap_,maxNumElements);

        if (this->OnCoarseSolveComm_) {
            const size_t numMyElementS = elemSGraph->getMap()->getLocalNumElements();
            Teuchos::ArrayView<const GO> va;
            for (UN i = 0; i < numMyElementS; i++) {
                GO kg = this->MLGatheringMaps_[this->MLGatheringMaps_.size()-1]->getGlobalElement(i);
                tmpElemGraph->getGlobalRowView(kg,va);
                Teuchos::Array<GO> vva(va);
                this->ElementNodeList_->insertGlobalIndices(kg,vva());//mal va nehmen
            }
            this->ElementNodeList_->fillComplete();

        }
        return 0;
    }

    template<class SC,class LO, class GO, class NO>
    int HarmonicCoarseOperator<SC,LO,GO,NO>::buildGlobalGraph(Teuchos::RCP<DDInterface<SC,LO,GO,NO> > theDDInterface)
    {
        FROSCH_DETAILTIMER_START_LEVELID(buildGlobalGraphTime,"HarmonicCoarseOperator::buildGlobalGraph");
        std::map<GO,int> rep;
        Teuchos::Array<GO> entries;
        IntVec2D conn;
        InterfaceEntityPtrVec connVec;
        int connrank;
        //get connected subdomains
        theDDInterface->identifyConnectivityEntities();
        EntitySetConstPtr connect=  theDDInterface->getConnectivityEntities();
        connect->buildEntityMap(theDDInterface->getNodesMap());
        connVec = connect->getEntityVector();
        connrank = connect->getEntityMap()->getComm()->getRank();

        GO connVecSize = connVec.size();
        conn.resize(connVecSize);
        {
            if (connVecSize>0) {
                for (GO i = 0;i<connVecSize;i++) {
                    conn[i] = connVec[i]->getSubdomainsVector();
                    for (int j = 0; j<conn[i].size(); j++) rep.insert(std::pair<GO,int>(conn.at(i).at(j),connrank));
                }
                for (auto& x: rep) {
                    entries.push_back(x.first);
                }
            }
        }

        const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > graphMap = Xpetra::MapFactory<LO,GO,NO>::Build(this->K_->getMap()->lib(),INVALID,1,0,this->K_->getMap()->getComm());

        //UN maxNumElements = -1;
        //get the maximum number of neighbors for a subdomain
        MaxNumNeigh_ = -1;
        UN numElementsLocal = entries.size();
        reduceAll(*this->MpiComm_,Teuchos::REDUCE_MAX,numElementsLocal,Teuchos::ptr(&MaxNumNeigh_));
        this->SubdomainConnectGraph_ = Xpetra::CrsGraphFactory<LO,GO,NO>::Build(graphMap,MaxNumNeigh_);
        this->SubdomainConnectGraph_->insertGlobalIndices(graphMap->getComm()->getRank(),entries());

        return 0;
    }

    template <class SC,class LO, class GO,class NO>
    int HarmonicCoarseOperator<SC,LO,GO,NO>::buildCoarseGraph()
    {
        //bring graph to CoarseSolveComm_
        FROSCH_DETAILTIMER_START_LEVELID(buildCoarseGraphTime,"CoarseOperator::buildCoarseGraph");

        //create temporary graphs to perform communication (possibly with gathering steps)
        GraphPtr tmpGraph =  Xpetra::CrsGraphFactory<LO,GO,NO>::Build(this->MLGatheringMaps_[1],MaxNumNeigh_);;
        GraphPtr tmpGraphGathering;

        tmpGraph->doExport(*this->SubdomainConnectGraph_,*this->MLCoarseSolveExporters_[1],Xpetra::INSERT);

        for (int i  = 2;i<this->MLGatheringMaps_.size();i++) {
            tmpGraph->fillComplete();
            tmpGraphGathering = tmpGraph;
            tmpGraph = Xpetra::CrsGraphFactory<LO,GO,NO>::Build(this->MLGatheringMaps_[i],MaxNumNeigh_);
            tmpGraph->doExport(*tmpGraphGathering,*this->MLCoarseSolveExporters_[i],Xpetra::INSERT);
        }
        const size_t numMyElementS = this->MLGatheringMaps_[this->MLGatheringMaps_.size()-1]->getLocalNumElements();
        this->SubdomainConnectGraph_= Xpetra::CrsGraphFactory<LO,GO,NO>::Build(this->MLCoarseMap_,MaxNumNeigh_);

        if (this->OnCoarseSolveComm_) {
            for (size_t k = 0; k<numMyElementS; k++) {
                Teuchos::ArrayView<const LO> in;
                Teuchos::ArrayView<const GO> vals_graph;
                GO kg = this->MLGatheringMaps_[this->MLGatheringMaps_.size()-1]->getGlobalElement(k);
                tmpGraph->getGlobalRowView(kg,vals_graph);
                Teuchos::Array<GO> vals(vals_graph);
                this->SubdomainConnectGraph_->insertGlobalIndices(kg,vals());
            }
            this->SubdomainConnectGraph_->fillComplete();
        }
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    void HarmonicCoarseOperator<SC,LO,GO,NO>::extractLocalSubdomainMatrix_Symbolic()
    {
        if (this->K_->getRowMap()->lib() == UseTpetra) {
            FROSCH_DETAILTIMER_START_LEVELID(extractLocalSubdomainMatrix_SymbolicTime,"HarmonicCoarseOperatorCoarseOperator::extractLocalSubdomainMatrix_Symbolic");
            XMapPtr repeatedMap = AssembleSubdomainMap(this->NumberOfBlocks_, this->DofsMaps_, this->DofsPerNode_);

            // buid sudomain matrix
            this->coarseSubdomainMatrix_ = MatrixFactory<SC,LO,GO,NO>::Build(repeatedMap, repeatedMap, this->K_->getGlobalMaxNumRowEntries());
            this->coarseScatter_ = ImportFactory<LO,GO,NO>::Build(this->K_->getRowMap(), repeatedMap);
            this->coarseSubdomainMatrix_->doImport(*(this->K_.getConst()), *(this->coarseScatter_),ADD);

            // build local subdomain matrix
            RCP<const Comm<LO> > SerialComm = rcp(new MpiComm<LO>(MPI_COMM_SELF));
            RCP<Map<LO,GO,NO> > localSubdomainMap = MapFactory<LO,GO,NO>::Build(repeatedMap->lib(), repeatedMap->getLocalNumElements(), 0, SerialComm);
            this->coarseLocalSubdomainMatrix_ = MatrixFactory<SC,LO,GO,NO>::Build(localSubdomainMap, localSubdomainMap, this->K_->getGlobalMaxNumRowEntries());

            // fill in column indexes
            ExtractLocalSubdomainMatrix_Symbolic(this->coarseSubdomainMatrix_, // input
                                                 this->coarseLocalSubdomainMatrix_); // output

            // turn flag on
            this->coarseExtractLocalSubdomainMatrix_Symbolic_Done_ = true;


            // -------------------------------------------------------
            // symbolic for computeExtensions to factor kII
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

            // extract local submatrices
            XMatrixPtr kII;
            XMatrixPtr kIGamma;
            XMatrixPtr kGammaI;
            XMatrixPtr kGammaGamma;
            auto repeatedMatrix = this->coarseLocalSubdomainMatrix_;
            BuildSubmatrices(repeatedMatrix.getConst(),indicesIDofsAll(),kII,kIGamma,kGammaI,kGammaGamma);

            // perform symbolic on kII
            ExtensionSolver_ = SolverFactory<SC,LO,GO,NO>::Build(kII,
                                                                 sublist(this->ParameterList_,"ExtensionSolver"),
                                                                 string("ExtensionSolver (Level ") + to_string(this->LevelID_) + string(")"));
            ExtensionSolver_->initialize();
        }
    }
    
}

#endif
