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

#ifndef _FROSCH_GDSWCOARSEOPERATOR_DEF_HPP
#define _FROSCH_GDSWCOARSEOPERATOR_DEF_HPP

#include <FROSch_GDSWCoarseOperator_decl.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    GDSWCoarseOperator<SC,LO,GO,NO>::GDSWCoarseOperator(ConstXMatrixPtr k,
                                                        ParameterListPtr parameterList) :
    HarmonicCoarseOperator<SC,LO,GO,NO> (k,parameterList),
    DDInterface_ ()
    {
        FROSCH_TIMER_START_LEVELID(gDSWCoarseOperatorTime,"GDSWCoarseOperator::GDSWCoarseOperator");
    }

    template <class SC,class LO,class GO,class NO>
    int GDSWCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                    ConstXMapPtr repeatedMap)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"GDSWCoarseOperator::initialize");
        buildCoarseSpace(dimension,repeatedMap);
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
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }


    template <class SC,class LO,class GO,class NO>
    void GDSWCoarseOperator<SC,LO,GO,NO>::describe(FancyOStream &out,
                                                   const EVerbosityLevel verbLevel) const
    {
        FROSCH_ASSERT(false,"describe() has be implemented properly...");
    }

    template <class SC,class LO,class GO,class NO>
    std::string GDSWCoarseOperator<SC,LO,GO,NO>::description() const
    {
        return "GDSW Coarse Operator";
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
        GOVecPtr dirichletBoundaryDofs = FindOneEntryOnlyRowsGlobal(this->K_.getConst(),nodesMap);
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
        GOVecPtr dirichletBoundaryDofs = FindOneEntryOnlyRowsGlobal(this->K_.getConst(),nodesMap);
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
        FROSCH_TIMER_START_LEVELID(buildCoarseSpaceTime,"GDSWCoarseOperator::buildCoarseSpace");
        FROSCH_ASSERT(dofsMaps.size()==dofsPerNode,"dofsMaps.size()!=dofsPerNode");

        // Das könnte man noch ändern
        // TODO: DAS SOLLTE ALLES IN EINE FUNKTION IN HARMONICCOARSEOPERATOR
        this->GammaDofs_.resize(this->GammaDofs_.size()+1);
        this->IDofs_.resize(this->IDofs_.size()+1);
        this->InterfaceCoarseSpaces_.resize(this->InterfaceCoarseSpaces_.size()+1);
        this->DofsMaps_.resize(this->DofsMaps_.size()+1);
        this->DofsPerNode_.resize(this->DofsPerNode_.size()+1);
        this->NumberOfBlocks_++;

        resetCoarseSpaceBlock(this->NumberOfBlocks_-1,dimension,dofsPerNode,nodesMap,dofsMaps,dirichletBoundaryDofs,nodeList);

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
        FROSCH_TIMER_START_LEVELID(buildCoarseSpaceTime,"GDSWCoarseOperator::buildCoarseSpace");
        // Das könnte man noch ändern
        // TODO: DAS SOLLTE ALLES IN EINE FUNKTION IN HARMONICCOARSEOPERATOR
        for (UN i=0; i<repeatedNodesMapVec.size(); i++) {
            this->GammaDofs_.resize(this->GammaDofs_.size()+1);
            this->IDofs_.resize(this->IDofs_.size()+1);
            this->InterfaceCoarseSpaces_.resize(this->InterfaceCoarseSpaces_.size()+1);
            this->DofsMaps_.resize(this->DofsMaps_.size()+1);
            this->DofsPerNode_.resize(this->DofsPerNode_.size()+1);
            this->NumberOfBlocks_++;
            resetCoarseSpaceBlock(this->NumberOfBlocks_-1,dimension,dofsPerNodeVec[i],repeatedNodesMapVec[i],repeatedDofMapsVec[i],dirichletBoundaryDofsVec[i],nodeListVec[i]);
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
        FROSCH_TIMER_START_LEVELID(resetCoarseSpaceBlockTime,"GDSWCoarseOperator::resetCoarseSpaceBlock");
        FROSCH_ASSERT(dofsMaps.size()==dofsPerNode,"dofsMaps.size()!=dofsPerNode");
        FROSCH_ASSERT(blockId<this->NumberOfBlocks_,"Block does not exist yet and can therefore not be reset.");

        if (this->Verbose_) {
            std::cout << "\n\
+--------------------+\n\
| GDSWCoarseOperator |\n\
|  Block " << blockId << "           |\n\
+--------------------+\n";
        }


        // Process the parameter list
        std::stringstream blockIdStringstream;
        blockIdStringstream << blockId+1;
        std::string blockIdString = blockIdStringstream.str();
        RCP<ParameterList> coarseSpaceList = sublist(sublist(this->ParameterList_,"Blocks"),blockIdString.c_str());

        CommunicationStrategy communicationStrategy = CreateOneToOneMap;
        if (!coarseSpaceList->get("Interface Communication Strategy","CreateOneToOneMap").compare("CrsMatrix")) {
            communicationStrategy = CommCrsMatrix;
        } else if (!coarseSpaceList->get("Interface Communication Strategy","CreateOneToOneMap").compare("CrsGraph")) {
            communicationStrategy = CommCrsGraph;
        } else if (!coarseSpaceList->get("Interface Communication Strategy","CreateOneToOneMap").compare("CreateOneToOneMap")) {
            communicationStrategy = CreateOneToOneMap;
        } else {
            FROSCH_ASSERT(false,"FROSch::GDSWCoarseOperator : ERROR: Specify a valid communication strategy for the identification of the interface components.");
        }

        Verbosity verbosity = All;
        if (!coarseSpaceList->get("Verbosity","All").compare("None")) {
            verbosity = None;
        } else if (!coarseSpaceList->get("Verbosity","All").compare("All")) {
            verbosity = All;
        } else {
            FROSCH_ASSERT(false,"FROSch::GDSWCoarseOperator : ERROR: Specify a valid verbosity level.");
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
            if (this->Verbose_) std::cout << "FROSch::GDSWCoarseOperator : WARNING: Rotations cannot be used" << std::endl;
        }
        if (!useRotations) {
            useShortEdgeRotations = false;
            useStraightEdgeRotations = false;
            useEdgeRotations = false;
            useFaceRotations = false;
        }

        this->DofsMaps_[blockId] = dofsMaps;
        this->DofsPerNode_[blockId] = dofsPerNode;

        Array<GO> tmpDirichletBoundaryDofs(dirichletBoundaryDofs()); // Here, we do a copy. Maybe, this is not necessary
        sortunique(tmpDirichletBoundaryDofs);

        DDInterface_.reset(new DDInterface<SC,LO,GO,NO>(dimension,this->DofsPerNode_[blockId],nodesMap.getConst(),verbosity,this->LevelID_,communicationStrategy));
        DDInterface_->resetGlobalDofs(dofsMaps);
        DDInterface_->removeDirichletNodes(tmpDirichletBoundaryDofs());
        if (this->ParameterList_->get("Test Unconnected Interface",true)) {
            DDInterface_->divideUnconnectedEntities(this->K_);
        }

        DDInterface_->sortVerticesEdgesFaces(nodeList);

        EntitySetPtr interface = DDInterface_->getInterface();
        EntitySetPtr interior = DDInterface_->getInterior();

        // Check for interface
        if (this->DofsPerNode_[blockId]*interface->getEntity(0)->getNumNodes()==0) {
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

            this->InterfaceCoarseSpaces_[blockId].reset(new CoarseSpace<SC,LO,GO,NO>());

            if (useForCoarseSpace && (useVertexTranslations||useShortEdgeTranslations||useShortEdgeRotations||useStraightEdgeTranslations||useStraightEdgeRotations||useEdgeTranslations||useEdgeRotations||useFaceTranslations||useFaceRotations)) {

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
                    for (UN i=0; i<translations.size(); i++) {
                        this->InterfaceCoarseSpaces_[blockId]->addSubspace(DDInterface_->getVertices()->getEntityMap(),translations[i]);
                    }
                }
                // ShortEdges
                if (useShortEdgeTranslations) {
                    XMultiVectorPtrVecPtr translations = this->computeTranslations(blockId,DDInterface_->getShortEdges());
                    for (UN i=0; i<translations.size(); i++) {
                        this->InterfaceCoarseSpaces_[blockId]->addSubspace(DDInterface_->getShortEdges()->getEntityMap(),translations[i]);
                    }
                }
                if (useShortEdgeRotations) {
                    XMultiVectorPtrVecPtr rotations = this->computeRotations(blockId,dimension,nodeList,DDInterface_->getShortEdges());
                    for (UN i=0; i<rotations.size(); i++) {
                        this->InterfaceCoarseSpaces_[blockId]->addSubspace(DDInterface_->getShortEdges()->getEntityMap(),rotations[i]);
                    }
                }
                // StraightEdges
                if (useStraightEdgeTranslations) {
                    XMultiVectorPtrVecPtr translations = this->computeTranslations(blockId,DDInterface_->getStraightEdges());
                    for (UN i=0; i<translations.size(); i++) {                        this->InterfaceCoarseSpaces_[blockId]->addSubspace(DDInterface_->getStraightEdges()->getEntityMap(),translations[i]);
                    }
                }
                if (useStraightEdgeRotations) {
                    XMultiVectorPtrVecPtr rotations = this->computeRotations(blockId,dimension,nodeList,DDInterface_->getStraightEdges());
                    for (UN i=0; i<rotations.size(); i++) {                        this->InterfaceCoarseSpaces_[blockId]->addSubspace(DDInterface_->getStraightEdges()->getEntityMap(),rotations[i]);
                    }
                }
                // Edges
                if (useEdgeTranslations) {
                    XMultiVectorPtrVecPtr translations = this->computeTranslations(blockId,DDInterface_->getEdges());
                    for (UN i=0; i<translations.size(); i++) {
                        this->InterfaceCoarseSpaces_[blockId]->addSubspace(DDInterface_->getEdges()->getEntityMap(),translations[i]);
                    }
                }
                if (useEdgeRotations) {
                    XMultiVectorPtrVecPtr rotations = this->computeRotations(blockId,dimension,nodeList,DDInterface_->getEdges());
                    for (UN i=0; i<rotations.size(); i++) {
                        this->InterfaceCoarseSpaces_[blockId]->addSubspace(DDInterface_->getEdges()->getEntityMap(),rotations[i]);
                    }
                }
                // Faces
                if (useFaceTranslations) {
                    XMultiVectorPtrVecPtr translations = this->computeTranslations(blockId,DDInterface_->getFaces());
                    for (UN i=0; i<translations.size(); i++) {
                        this->InterfaceCoarseSpaces_[blockId]->addSubspace(DDInterface_->getFaces()->getEntityMap(),translations[i]);
                    }
                }
                if (useFaceRotations) {
                    XMultiVectorPtrVecPtr rotations = this->computeRotations(blockId,dimension,nodeList,DDInterface_->getFaces());
                    for (UN i=0; i<rotations.size(); i++) {
                        this->InterfaceCoarseSpaces_[blockId]->addSubspace(DDInterface_->getFaces()->getEntityMap(),rotations[i]);
                    }
                }

                this->InterfaceCoarseSpaces_[blockId]->assembleCoarseSpace();

                if (this->Verbose_) {
                    std::cout << std::boolalpha << "\n\
    ------------------------------------------------------------------------------\n\
     GDSW coarse space\n\
    ------------------------------------------------------------------------------\n\
      Vertices: translations                      --- " << useVertexTranslations << "\n\
      ShortEdges: translations                    --- " << useShortEdgeTranslations << "\n\
      ShortEdges: rotations                       --- " << useShortEdgeRotations << "\n\
      StraightEdges: translations                 --- " << useStraightEdgeTranslations << "\n\
      StraightEdges: rotations                    --- " << useStraightEdgeRotations << "\n\
      Edges: translations                         --- " << useEdgeTranslations << "\n\
      Edges: rotations                            --- " << useEdgeRotations << "\n\
      Faces: translations                         --- " << useFaceTranslations << "\n\
      Faces: rotations                            --- " << useFaceRotations << "\n\
    ------------------------------------------------------------------------------\n" << std::noboolalpha;
                }
            }
        }
        return 0;
    }
}

#endif
