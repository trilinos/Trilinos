// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_GDSWCOARSEOPERATOR_DEF_HPP
#define _FROSCH_GDSWCOARSEOPERATOR_DEF_HPP

#include <FROSch_GDSWCoarseOperator_decl.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    GDSWCoarseOperator<SC,LO,GO,NO>::GDSWCoarseOperator(ConstXMatrixPtr k,
                                                        ParameterListPtr parameterList) :
    HarmonicCoarseOperator<SC,LO,GO,NO> (k,parameterList)
    {
        FROSCH_DETAILTIMER_START_LEVELID(gDSWCoarseOperatorTime,"GDSWCoarseOperator::GDSWCoarseOperator");
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                    ConstXMapPtr repeatedMap)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"GDSWCoarseOperator::initialize");
        buildCoarseSpace(dimension,repeatedMap);
        this->assembleInterfaceCoarseSpace();
        this->buildCoarseSolveMap(this->AssembledInterfaceCoarseSpace_->getBasisMapUnique());
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                    ConstXMapPtr repeatedMap,
                                                    GOVecPtr dirichletBoundaryDofs)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"GDSWCoarseOperator::initialize");
        buildCoarseSpace(dimension,repeatedMap,dirichletBoundaryDofs);
        this->assembleInterfaceCoarseSpace();
        this->buildCoarseSolveMap(this->AssembledInterfaceCoarseSpace_->getBasisMapUnique());
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UN dofsPerNode,
                                                    ConstXMapPtr repeatedNodesMap,
                                                    ConstXMapPtrVecPtr repeatedDofMaps)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"GDSWCoarseOperator::initialize");
        buildCoarseSpace(dimension,dofsPerNode,repeatedNodesMap,repeatedDofMaps);
        this->assembleInterfaceCoarseSpace();
        this->buildCoarseSolveMap(this->AssembledInterfaceCoarseSpace_->getBasisMapUnique());
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UN dofsPerNode,
                                                    ConstXMapPtr repeatedNodesMap,
                                                    ConstXMapPtrVecPtr repeatedDofMaps,
                                                    GOVecPtr dirichletBoundaryDofs)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"GDSWCoarseOperator::initialize");
        buildCoarseSpace(dimension,dofsPerNode,repeatedNodesMap,repeatedDofMaps,dirichletBoundaryDofs);
        this->assembleInterfaceCoarseSpace();
        this->buildCoarseSolveMap(this->AssembledInterfaceCoarseSpace_->getBasisMapUnique());
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UN dofsPerNode,
                                                    ConstXMapPtr repeatedNodesMap,
                                                    ConstXMapPtrVecPtr repeatedDofMaps,
                                                    ConstXMultiVectorPtr nodeList)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"GDSWCoarseOperator::initialize");
        buildCoarseSpace(dimension,dofsPerNode,repeatedNodesMap,repeatedDofMaps,nodeList);
        this->assembleInterfaceCoarseSpace();
        this->buildCoarseSolveMap(this->AssembledInterfaceCoarseSpace_->getBasisMapUnique());
        this->extractLocalSubdomainMatrix_Symbolic();
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UN dofsPerNode,
                                                    ConstXMapPtr repeatedNodesMap,
                                                    ConstXMapPtrVecPtr repeatedDofMaps,
                                                    GOVecPtr dirichletBoundaryDofs,
                                                    ConstXMultiVectorPtr nodeList)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"GDSWCoarseOperator::initialize");
        buildCoarseSpace(dimension,dofsPerNode,repeatedNodesMap,repeatedDofMaps,dirichletBoundaryDofs,nodeList);
        this->assembleInterfaceCoarseSpace();
        this->buildCoarseSolveMap(this->AssembledInterfaceCoarseSpace_->getBasisMapUnique());
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                    UNVecPtr dofsPerNodeVec,
                                                    ConstXMapPtrVecPtr repeatedNodesMapVec,
                                                    ConstXMapPtrVecPtr2D repeatedDofMapsVec,
                                                    GOVecPtr2D dirichletBoundaryDofsVec,
                                                    ConstXMultiVectorPtrVecPtr nodeListVec)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"GDSWCoarseOperator::initialize");
        buildCoarseSpace(dimension,dofsPerNodeVec,repeatedNodesMapVec,repeatedDofMapsVec,dirichletBoundaryDofsVec,nodeListVec);
        this->assembleInterfaceCoarseSpace();
        this->buildCoarseSolveMap(this->AssembledInterfaceCoarseSpace_->getBasisMapUnique());
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    void GDSWCoarseOperator<SC,LO,GO,NO>::describe(FancyOStream &out,
                                                   const EVerbosityLevel verbLevel) const
    {
        FROSCH_ASSERT(false,"describe() has to be implemented properly...");
    }

    template <class SC,class LO,class GO,class NO>
    string GDSWCoarseOperator<SC,LO,GO,NO>::description() const
    {
        return "GDSW Coarse Operator";
    }

    template<class SC,class LO,class GO,class NO>
    typename GDSWCoarseOperator<SC,LO,GO,NO>::XMapPtr GDSWCoarseOperator<SC,LO,GO,NO>::BuildRepeatedMapCoarseLevel(ConstXMapPtr &nodesMap,
                                                UN dofsPerNode,
                                                ConstXMapPtrVecPtr dofsMaps,
                                                UN partitionType)
    {
      FROSCH_ASSERT(false,"For GDSWCoarseOperator the ZoltanDual Option is not implemented!");
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                          ConstXMapPtr nodesMap)
    {
        ConstXMapPtrVecPtr dofsMaps(1);
        dofsMaps[0] = nodesMap;
        buildCoarseSpace(dimension,1,nodesMap,dofsMaps);

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                          ConstXMapPtr nodesMap,
                                                          GOVecPtr dirichletBoundaryDofs)
    {
        ConstXMapPtrVecPtr dofsMaps(1);
        dofsMaps[0] = nodesMap;
        buildCoarseSpace(dimension,1,nodesMap,dofsMaps,dirichletBoundaryDofs);

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                          UN dofsPerNode,
                                                          ConstXMapPtr nodesMap,
                                                          ConstXMapPtrVecPtr dofsMaps)
    {
/*
#ifdef FindOneEntryOnlyRowsGlobal_Matrix
        GOVecPtr dirichletBoundaryDofs = FindOneEntryOnlyRowsGlobal(this->K_.getConst(),nodesMap);
#else
        GOVecPtr dirichletBoundaryDofs = FindOneEntryOnlyRowsGlobal(this->K_->getCrsGraph(),nodesMap);
#end
 */
        FROSCH_WARNING("FROSch::GDSWCoarseOperator",this->Verbose_,"We do not have the right map (repeatedMap) to use FindOneEntryOnlyRowsGlobal. A variant that uses the row map could be implemented?! => We use dirichletBoundaryDofs = null for now.");
        GOVecPtr dirichletBoundaryDofs = null;
        buildCoarseSpace(dimension,dofsPerNode,nodesMap,dofsMaps,dirichletBoundaryDofs);

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                          UN dofsPerNode,
                                                          ConstXMapPtr nodesMap,
                                                          ConstXMapPtrVecPtr dofsMaps,
                                                          GOVecPtr dirichletBoundaryDofs)
    {
        ConstXMultiVectorPtr nodeList;
        buildCoarseSpace(dimension,dofsPerNode,nodesMap,dofsMaps,dirichletBoundaryDofs,nodeList);

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                          UN dofsPerNode,
                                                          ConstXMapPtr nodesMap,
                                                          ConstXMapPtrVecPtr dofsMaps,
                                                          ConstXMultiVectorPtr nodeList)
    {
/*
#ifdef FindOneEntryOnlyRowsGlobal_Matrix
        GOVecPtr dirichletBoundaryDofs = FindOneEntryOnlyRowsGlobal(this->K_.getConst(),nodesMap);
#else
        GOVecPtr dirichletBoundaryDofs = FindOneEntryOnlyRowsGlobal(this->K_->getCrsGraph(),nodesMap);
#end
 */
        FROSCH_WARNING("FROSch::GDSWCoarseOperator",this->Verbose_,"We do not have the right map (repeatedMap) to use FindOneEntryOnlyRowsGlobal. A variant that uses the row map could be implemented?! => We use dirichletBoundaryDofs = null for now.");
        GOVecPtr dirichletBoundaryDofs = null;
        buildCoarseSpace(dimension,dofsPerNode,nodesMap,dofsMaps,dirichletBoundaryDofs,nodeList);

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                          UN dofsPerNode,
                                                          ConstXMapPtr nodesMap,
                                                          ConstXMapPtrVecPtr dofsMaps,
                                                          GOVecPtr dirichletBoundaryDofs,
                                                          ConstXMultiVectorPtr nodeList)
    {
        FROSCH_DETAILTIMER_START_LEVELID(buildCoarseSpaceTime,"GDSWCoarseOperator::buildCoarseSpace");
        FROSCH_ASSERT(dofsMaps.size()==dofsPerNode,"dofsMaps.size()!=dofsPerNode");

        // Das könnte man noch ändern
        // TODO: DAS SOLLTE ALLES IN EINE FUNKTION IN HARMONICCOARSEOPERATOR
        resetCoarseSpaceBlock(this->NumberOfBlocks_,dimension,dofsPerNode,nodesMap,dofsMaps,dirichletBoundaryDofs,nodeList);

        return 0;
    }



    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                          UNVecPtr dofsPerNodeVec,
                                                          ConstXMapPtrVecPtr repeatedNodesMapVec,
                                                          ConstXMapPtrVecPtr2D repeatedDofMapsVec,
                                                          GOVecPtr2D dirichletBoundaryDofsVec,
                                                          ConstXMultiVectorPtrVecPtr nodeListVec)
    {
        FROSCH_DETAILTIMER_START_LEVELID(buildCoarseSpaceTime,"GDSWCoarseOperator::buildCoarseSpace");
        // Das könnte man noch ändern
        // TODO: DAS SOLLTE ALLES IN EINE FUNKTION IN HARMONICCOARSEOPERATOR
        for (UN i=0; i<repeatedNodesMapVec.size(); i++) {
            resetCoarseSpaceBlock(this->NumberOfBlocks_,dimension,dofsPerNodeVec[i],repeatedNodesMapVec[i],repeatedDofMapsVec[i],dirichletBoundaryDofsVec[i],nodeListVec[i]);
        }
        return 0;
    }


    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::resetCoarseSpaceBlock(UN blockId,
                                                               UN dimension,
                                                               UN dofsPerNode,
                                                               ConstXMapPtr nodesMap,
                                                               ConstXMapPtrVecPtr dofsMaps,
                                                               GOVecPtr dirichletBoundaryDofs,
                                                               ConstXMultiVectorPtr nodeList)
    {
        FROSCH_DETAILTIMER_START_LEVELID(resetCoarseSpaceBlockTime,"GDSWCoarseOperator::resetCoarseSpaceBlock");
        FROSCH_ASSERT(dofsMaps.size()==dofsPerNode,"dofsMaps.size()!=dofsPerNode");
        FROSCH_ASSERT(blockId<=this->NumberOfBlocks_,"Block does not exist yet and can therefore not be reset("+to_string(blockId)+" <= "+to_string(this->NumberOfBlocks_)+". ");

        if (!this->DistributionList_->get("Type","linear").compare("ZoltanDual")) {
            FROSCH_ASSERT(false,"RGDSWCoarseOperator:: Distribution Type ZoltanDual only works for IPOUHarmonicCoarseOperator");
        }

        // Process the parameter list
        stringstream blockIdStringstream;
        blockIdStringstream << blockId+1;
        string blockIdString = blockIdStringstream.str();
        RCP<ParameterList> coarseSpaceList = sublist(sublist(this->ParameterList_,"Blocks"),blockIdString.c_str());

        CommunicationStrategy communicationStrategy = CreateOneToOneMap;
        if (!coarseSpaceList->get("Interface Communication Strategy","CreateOneToOneMap").compare("CrsMatrix")) {
            communicationStrategy = CommCrsMatrix;
        } else if (!coarseSpaceList->get("Interface Communication Strategy","CreateOneToOneMap").compare("CrsGraph")) {
            communicationStrategy = CommCrsGraph;
        } else if (!coarseSpaceList->get("Interface Communication Strategy","CreateOneToOneMap").compare("CreateOneToOneMap")) {
            communicationStrategy = CreateOneToOneMap;
        } else {
            FROSCH_ASSERT(false,"FROSch::GDSWCoarseOperator: Specify a valid communication strategy for the identification of the interface components.");
        }

        Verbosity verbosity = All;
        if (!coarseSpaceList->get("Verbosity","All").compare("None")) {
            verbosity = None;
        } else if (!coarseSpaceList->get("Verbosity","All").compare("All")) {
            verbosity = All;
        } else {
            FROSCH_ASSERT(false,"FROSch::GDSWCoarseOperator: Specify a valid verbosity level.");
        }

        bool useForCoarseSpace = coarseSpaceList->get("Use For Coarse Space",true);

        bool useVertexTranslations = coarseSpaceList->sublist("Custom").get("Vertices: translations",true);

        bool useShortEdgeTranslations = coarseSpaceList->sublist("Custom").get("ShortEdges: translations",true);
        bool useShortEdgeRotations = coarseSpaceList->sublist("Custom").get("ShortEdges: rotations",true);

        bool useStraightEdgeTranslations = coarseSpaceList->sublist("Custom").get("StraightEdges: translations",true);
        bool useStraightEdgeRotations = coarseSpaceList->sublist("Custom").get("StraightEdges: rotations",true);

        bool useEdgeTranslations = coarseSpaceList->sublist("Custom").get("Edges: translations",true);
        bool useEdgeRotations = coarseSpaceList->sublist("Custom").get("Edges: rotations",true);

        bool useFaceTranslations = coarseSpaceList->sublist("Custom").get("Faces: translations",true);
        bool useFaceRotations = coarseSpaceList->sublist("Custom").get("Faces: rotations",true);

        bool useRotations = coarseSpaceList->get("Rotations",true);
        if (useRotations && nodeList.is_null()) {
            useRotations = false;
            FROSCH_WARNING("FROSch::GDSWCoarseOperator",this->Verbose_,"Rotations cannot be used since nodeList.is_null().");
        }
        if (!useRotations) {
            useShortEdgeRotations = false;
            useStraightEdgeRotations = false;
            useEdgeRotations = false;
            useFaceRotations = false;
        }

        if (useForCoarseSpace) {
            this->NumberOfBlocks_++;

            this->GammaDofs_.resize(this->GammaDofs_.size()+1);
            this->IDofs_.resize(this->IDofs_.size()+1);
            this->InterfaceCoarseSpaces_.resize(this->InterfaceCoarseSpaces_.size()+1);
            this->DofsMaps_.resize(this->DofsMaps_.size()+1);
            this->DofsPerNode_.resize(this->DofsPerNode_.size()+1);

            this->DofsMaps_[blockId] = dofsMaps;
            this->DofsPerNode_[blockId] = dofsPerNode;

            Array<GO> tmpDirichletBoundaryDofs(dirichletBoundaryDofs()); // Here, we do a copy. Maybe, this is not necessary
            sortunique(tmpDirichletBoundaryDofs);

            DDInterface_.reset(new DDInterface<SC,LO,GO,NO>(dimension,this->DofsPerNode_[blockId],nodesMap.getConst(),verbosity,this->LevelID_,communicationStrategy));
            DDInterface_->resetGlobalDofs(dofsMaps);
            DDInterface_->removeDirichletNodes(tmpDirichletBoundaryDofs());

            if (useVertexTranslations||useShortEdgeTranslations||useShortEdgeRotations||useStraightEdgeTranslations||useStraightEdgeRotations||useEdgeTranslations||useEdgeRotations||useFaceTranslations||useFaceRotations) {
                EntitySetPtr interface = this->DDInterface_->getInterface();
                EntitySetPtr interior = this->DDInterface_->getInterior();

                if (this->Verbose_) {
                    cout
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << setw(89) << "-----------------------------------------------------------------------------------------"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << "| "
                    << left << setw(74) << "GDSWCoarseOperator " << right << setw(8) << "(Level " << setw(2) << this->LevelID_ << ")"
                    << " |"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << setw(89) << "========================================================================================="
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << "| " << left << setw(41) << "Block" << right
                    << " | " << setw(41) << blockId
                    << " |"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << "| " << left << setw(41) << "Spatial dimensions" << right
                    << " | " << setw(41) << dimension
                    << " |"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << "| " << left << setw(41) << "Number of degrees of freedom per node" << right
                    << " | " << setw(41) << dofsPerNode
                    << " |"
                    << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                    << setw(89) << "-----------------------------------------------------------------------------------------"
                    << endl;
                }

                // Check for interface
                if (interface->getEntity(0)->getNumNodes()==0) {
                    FROSCH_NOTIFICATION("FROSch::GDSWCoarseOperator",this->Verbose_,"No interface found => Volume functions will be used instead.");
                    this->computeVolumeFunctions(blockId,dimension,nodesMap,nodeList,interior);
                } else {
                    this->GammaDofs_[blockId] = LOVecPtr(this->DofsPerNode_[blockId]*interface->getEntity(0)->getNumNodes());
                    this->IDofs_[blockId] = LOVecPtr(this->DofsPerNode_[blockId]*interior->getEntity(0)->getNumNodes());
                    for (UN k=0; k<this->DofsPerNode_[blockId]; k++) {
                        for (UN i=0; i<interface->getEntity(0)->getNumNodes(); i++) {
                            this->GammaDofs_[blockId][this->DofsPerNode_[blockId]*i+k] = interface->getEntity(0)->getLocalDofID(i,k);
                        }
                        for (UN i=0; i<interior->getEntity(0)->getNumNodes(); i++) {
                            this->IDofs_[blockId][this->DofsPerNode_[blockId]*i+k] = interior->getEntity(0)->getLocalDofID(i,k);
                        }
                    }

                    this->InterfaceCoarseSpaces_[blockId].reset(new CoarseSpace<SC,LO,GO,NO>(this->MpiComm_,this->SerialComm_));

                    if (this->ParameterList_->get("Test Unconnected Interface",true)) {
                        DDInterface_->divideUnconnectedEntities(this->K_);
                    }

                    DDInterface_->sortVerticesEdgesFaces(nodeList);

                    EntitySetPtr interface = DDInterface_->getInterface();
                    EntitySetPtr interior = DDInterface_->getInterior();

                    ////////////////////////////////
                    // Build Processor Map Coarse //
                    ////////////////////////////////
                    DDInterface_->buildEntityMaps(useVertexTranslations,
                                                  useShortEdgeTranslations||useShortEdgeRotations,
                                                  useStraightEdgeTranslations || useStraightEdgeRotations,
                                                  useEdgeTranslations || useEdgeRotations,
                                                  useFaceTranslations || useFaceRotations,
                                                  false);

                    // Vertices
                    if (useVertexTranslations) {
                        XMultiVectorPtrVecPtr translations = this->computeTranslations(blockId,DDInterface_->getVertices());
                        ConstXMapPtr verticesEntityMap = DDInterface_->getVertices()->getEntityMap();
                        for (UN i=0; i<translations.size(); i++) {
                            this->InterfaceCoarseSpaces_[blockId]->addSubspace(verticesEntityMap,null,translations[i]);
                        }
                    }
                    // ShortEdges
                    if (useShortEdgeTranslations) {
                        XMultiVectorPtrVecPtr translations = this->computeTranslations(blockId,DDInterface_->getShortEdges());
                        ConstXMapPtr shortEdgesEntityMap = DDInterface_->getShortEdges()->getEntityMap();
                        for (UN i=0; i<translations.size(); i++) {
                            this->InterfaceCoarseSpaces_[blockId]->addSubspace(shortEdgesEntityMap,null,translations[i]);
                        }
                    }
                    if (useShortEdgeRotations) {
                        XMultiVectorPtrVecPtr rotations = this->computeRotations(blockId,dimension,nodeList,DDInterface_->getShortEdges(),(dimension==3));
                        ConstXMapPtr shortEdgesEntityMap = DDInterface_->getShortEdges()->getEntityMap();
                        for (UN i=0; i<rotations.size(); i++) {
                            this->InterfaceCoarseSpaces_[blockId]->addSubspace(shortEdgesEntityMap,null,rotations[i]);
                        }
                    }
                    // StraightEdges
                    if (useStraightEdgeTranslations) {
                        XMultiVectorPtrVecPtr translations = this->computeTranslations(blockId,DDInterface_->getStraightEdges());
                        ConstXMapPtr straightEdgesEntityMap = DDInterface_->getStraightEdges()->getEntityMap();
                        for (UN i=0; i<translations.size(); i++) {
                            this->InterfaceCoarseSpaces_[blockId]->addSubspace(straightEdgesEntityMap,null,translations[i]);
                        }
                    }
                    if (useStraightEdgeRotations) {
                        XMultiVectorPtrVecPtr rotations = this->computeRotations(blockId,dimension,nodeList,DDInterface_->getStraightEdges(),(dimension==3));
                        ConstXMapPtr straightEdgesEntityMap = DDInterface_->getStraightEdges()->getEntityMap();
                        for (UN i=0; i<rotations.size(); i++) {
                            this->InterfaceCoarseSpaces_[blockId]->addSubspace(straightEdgesEntityMap,null,rotations[i]);
                        }
                    }
                    // Edges
                    if (useEdgeTranslations) {
                        XMultiVectorPtrVecPtr translations = this->computeTranslations(blockId,DDInterface_->getEdges());
                        ConstXMapPtr edgesEntityMap = DDInterface_->getEdges()->getEntityMap();
                        for (UN i=0; i<translations.size(); i++) {
                            this->InterfaceCoarseSpaces_[blockId]->addSubspace(edgesEntityMap,null,translations[i]);
                        }
                    }
                    if (useEdgeRotations) {
                        XMultiVectorPtrVecPtr rotations = this->computeRotations(blockId,dimension,nodeList,DDInterface_->getEdges());
                        ConstXMapPtr edgesEntityMap = DDInterface_->getEdges()->getEntityMap();
                        for (UN i=0; i<rotations.size(); i++) {
                            this->InterfaceCoarseSpaces_[blockId]->addSubspace(edgesEntityMap,null,rotations[i]);
                        }
                    }
                    // Faces
                    if (useFaceTranslations) {
                        XMultiVectorPtrVecPtr translations = this->computeTranslations(blockId,DDInterface_->getFaces());
                        ConstXMapPtr facesEntityMap = DDInterface_->getFaces()->getEntityMap();
                        for (UN i=0; i<translations.size(); i++) {
                            this->InterfaceCoarseSpaces_[blockId]->addSubspace(facesEntityMap,null,translations[i]);
                        }
                    }
                    if (useFaceRotations) {
                        XMultiVectorPtrVecPtr rotations = this->computeRotations(blockId,dimension,nodeList,DDInterface_->getFaces());
                        ConstXMapPtr facesEntityMap = DDInterface_->getFaces()->getEntityMap();
                        for (UN i=0; i<rotations.size(); i++) {
                            this->InterfaceCoarseSpaces_[blockId]->addSubspace(facesEntityMap,null,rotations[i]);
                        }
                    }

                    this->InterfaceCoarseSpaces_[blockId]->assembleCoarseSpace();

                    if (this->Verbose_) {
                        cout
                        << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                        << setw(89) << "-----------------------------------------------------------------------------------------"
                        << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                        << "| "
                        << left << setw(74) << "> GDSW coarse space " << right << setw(8) << "(Level " << setw(2) << this->LevelID_ << ")"
                        << " |"
                        << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                        << setw(89) << "========================================================================================="
                        << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                        << "| " << left << setw(19) << "Vertices " << " | " << setw(19) << "Translations " << right
                        << " | " << setw(41) << boolalpha << useVertexTranslations << noboolalpha
                        << " |"
                        << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                        << "| " << left << setw(19) << "ShortEdges " << " | " << setw(19) << "Translations " << right
                        << " | " << setw(41) << boolalpha << useShortEdgeTranslations << noboolalpha
                        << " |"
                        << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                        << "| " << left << setw(19) << "ShortEdges " << " | " << setw(19) << "Rotations " << right
                        << " | " << setw(41) << boolalpha << useShortEdgeRotations << noboolalpha
                        << " |"
                        << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                        << "| " << left << setw(19) << "StraightEdges " << " | " << setw(19) << "Translations " << right
                        << " | " << setw(41) << boolalpha << useStraightEdgeTranslations << noboolalpha
                        << " |"
                        << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                        << "| " << left << setw(19) << "StraightEdges " << " | " << setw(19) << "Rotations " << right
                        << " | " << setw(41) << boolalpha << useStraightEdgeRotations << noboolalpha
                        << " |"
                        << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                        << "| " << left << setw(19) << "Edges " << " | " << setw(19) << "Translations " << right
                        << " | " << setw(41) << boolalpha << useEdgeTranslations << noboolalpha
                        << " |"
                        << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                        << "| " << left << setw(19) << "Edges " << " | " << setw(19) << "Rotations " << right
                        << " | " << setw(41) << boolalpha << useEdgeRotations << noboolalpha
                        << " |"
                        << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                        << "| " << left << setw(19) << "Faces " << " | " << setw(19) << "Translations " << right
                        << " | " << setw(41) << boolalpha << useFaceTranslations << noboolalpha
                        << " |"
                        << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                        << "| " << left << setw(19) << "Faces " << " | " << setw(19) << "Rotations " << right
                        << " | " << setw(41) << boolalpha << useFaceRotations << noboolalpha
                        << " |"
                        << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                        << setw(89) << "-----------------------------------------------------------------------------------------"
                        << endl;
                    }
                }
            }
        }
        return 0;
    }
}

#endif
