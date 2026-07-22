// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_IPOUHARMONICCOARSEOPERATOR_DEF_HPP
#define _FROSCH_IPOUHARMONICCOARSEOPERATOR_DEF_HPP

#include <FROSch_IPOUHarmonicCoarseOperator_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    IPOUHarmonicCoarseOperator<SC,LO,GO,NO>::IPOUHarmonicCoarseOperator(ConstXMatrixPtr k,
                                                                         ParameterListPtr parameterList) :
    HarmonicCoarseOperator<SC,LO,GO,NO> (k,parameterList)
    {
        FROSCH_DETAILTIMER_START_LEVELID(iPOUHarmonicCoarseOperatorTime,"IPOUHarmonicCoarseOperator::IPOUHarmonicCoarseOperator");
    }

    template <class SC,class LO,class GO,class NO>
    int IPOUHarmonicCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                            UN dofsPerNode,
                                                            ConstXMapPtr nodesMap,
                                                            ConstXMapPtrVecPtr dofsMaps,
                                                            ConstXMultiVectorPtr nullSpaceBasis,
                                                            ConstXMultiVectorPtr nodeList,
                                                            GOVecPtr dirichletBoundaryDofs)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"IPOUHarmonicCoarseOperator::initialize");
        int ret = buildCoarseSpace(dimension,dofsPerNode,nodesMap,dofsMaps,nullSpaceBasis,dirichletBoundaryDofs,nodeList);
        this->CoarseMap_ = this->assembleCoarseMap();
        this->assembleInterfaceCoarseSpace();
        this->buildCoarseSolveMap(this->AssembledInterfaceCoarseSpace_->getBasisMapUnique());
        this->extractLocalSubdomainMatrix_Symbolic();
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return ret;
    }

    template <class SC,class LO,class GO,class NO>
    int IPOUHarmonicCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                            UNVecPtr dofsPerNodeVec,
                                                            ConstXMapPtrVecPtr repeatedNodesMapVec,
                                                            ConstXMapPtrVecPtr2D repeatedDofMapsVec,
                                                            ConstXMultiVectorPtrVecPtr nullSpaceBasisVec,
                                                            ConstXMultiVectorPtrVecPtr nodeListVec,
                                                            GOVecPtr2D dirichletBoundaryDofsVec)
    {
        FROSCH_TIMER_START_LEVELID(initializeTime,"IPOUHarmonicCoarseOperator::initialize");
        buildCoarseSpace(dimension,dofsPerNodeVec,repeatedNodesMapVec,repeatedDofMapsVec,nullSpaceBasisVec,dirichletBoundaryDofsVec,nodeListVec);
        this->CoarseMap_ = this->assembleCoarseMap();
        this->assembleInterfaceCoarseSpace();
        this->buildCoarseSolveMap(this->AssembledInterfaceCoarseSpace_->getBasisMapUnique());
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    void IPOUHarmonicCoarseOperator<SC,LO,GO,NO>::describe(FancyOStream &out,
                                                           const EVerbosityLevel verbLevel) const
    {
        FROSCH_ASSERT(false,"describe() has to be implemented properly...");
    }

    template <class SC,class LO,class GO,class NO>
    string IPOUHarmonicCoarseOperator<SC,LO,GO,NO>::description() const
    {
        return "Interface Partition of Unity Coarse Operator";
    }

    template <class SC,class LO,class GO,class NO>
    typename IPOUHarmonicCoarseOperator<SC,LO,GO,NO>::XMapPtr IPOUHarmonicCoarseOperator<SC,LO,GO,NO>::BuildRepeatedMapCoarseLevel(ConstXMapPtr &nodesMap,
                                                                                                                                   UN dofsPerNode,
                                                                                                                                   ConstXMapPtrVecPtr dofsMaps,
                                                                                                                                   UN partition)
    {
        //FROSCH_ASSERT(numVert+numEdg+numFac != nodesMap->getGlobalNumElements(),"ERROR: Map does not match number of Entities");
        const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();
        UN blockId = 1; //This is not implemented for  Block Variant yet
        stringstream blockIdStringstream;
        blockIdStringstream << blockId+1;
        string blockIdString = blockIdStringstream.str();
        RCP<ParameterList> coarseSpaceList = sublist(sublist(this->ParameterList_,"Blocks"),blockIdString.c_str());

        Teuchos::Array<GO> nodeEle = nodesMap->getLocalElementList();

        Teuchos::Array<GO> dofEle(nodeEle.size()*dofsPerNode);
        Teuchos::Array<GO> dmapEle(nodeEle.size());
        //GDSW Type CoarseOperator
        if (partition == 0) {
            for (unsigned j = 0; j<dofsPerNode; j++) {
                for ( unsigned i = 0; i<nodeEle.size(); i++) {
                    //vertices
                    if (nodeEle[i]<NumEnt_[0]) {
                        dofEle[i*dofsPerNode+j] = nodeEle[i]+j*NumEnt_[0];
                        dmapEle[i] = nodeEle[i]+j*NumEnt_[0];
                    }
                    //shortEdges
                    else if (nodeEle[i]>=NumEnt_[0] && nodeEle[i]<NumEnt_[1]+NumEnt_[0]) {
                        dofEle[i*dofsPerNode+j] = nodeEle[i]+j*NumEnt_[1]+(dofsPerNode-1)*NumEnt_[0];
                        dmapEle[i] = nodeEle[i]+j*NumEnt_[1]+(dofsPerNode-1)*NumEnt_[0];
                    }
                    //straightEdges
                    else if (nodeEle[i]<NumEnt_[2]+NumEnt_[0]+NumEnt_[1] && nodeEle[i]>=NumEnt_[1]+NumEnt_[0]) {
                        dofEle[i*dofsPerNode+j] = nodeEle[i]+j*NumEnt_[2]+(NumEnt_[0]*(dofsPerNode-1))+(NumEnt_[1]*(dofsPerNode-1));
                        dmapEle[i] = nodeEle[i]+j*NumEnt_[2]+(NumEnt_[0]*(dofsPerNode-1))+(NumEnt_[1]*(dofsPerNode-1));
                    }
                    //edges
                    else if (nodeEle[i]<NumEnt_[3]+NumEnt_[2]+NumEnt_[0]+NumEnt_[1] && nodeEle[i]>=NumEnt_[1]+NumEnt_[0]+NumEnt_[2]) {
                        dofEle[i*dofsPerNode+j] = nodeEle[i]+j*NumEnt_[3]+(NumEnt_[0]*(dofsPerNode-1))+(NumEnt_[1]*(dofsPerNode-1))+(NumEnt_[2]*(dofsPerNode-1));
                        dmapEle[i] = nodeEle[i]+j*NumEnt_[3]+(NumEnt_[0]*(dofsPerNode-1))+(NumEnt_[1]*(dofsPerNode-1))+(NumEnt_[2]*(dofsPerNode-1));
                    }
                    //faces
                    else if (nodeEle[i]<NumEnt_[3]+NumEnt_[2]+NumEnt_[0]+NumEnt_[1]+NumEnt_[4] && nodeEle[i]>=NumEnt_[1]+NumEnt_[0]+NumEnt_[2]+NumEnt_[3]) {
                        dofEle[i*dofsPerNode+j] = nodeEle[i]+j*NumEnt_[4]+(NumEnt_[0]*(dofsPerNode-1))+(NumEnt_[1]*(dofsPerNode-1))+(NumEnt_[2]*(dofsPerNode-1))+(NumEnt_[3]*(dofsPerNode-1));
                        dmapEle[i] = nodeEle[i]+j*NumEnt_[4]+(NumEnt_[0]*(dofsPerNode-1))+(NumEnt_[1]*(dofsPerNode-1))+(NumEnt_[2]*(dofsPerNode-1))+(NumEnt_[3]*(dofsPerNode-1));
                    } else {
                        if (nodesMap->getComm()->getRank() == 0) std::cout<<"This should never happen\n";
                    }
                }
                dofsMaps[j] =   Xpetra::MapFactory<LO,GO,NO>::Build(Xpetra::UseTpetra,INVALID,dmapEle,0,nodesMap->getComm());
            }
        } //RGDSW type CoarseOperator
        else if (partition == 1) {
            for (unsigned j = 0;j<dofsPerNode;j++) {
                for (unsigned i = 0; i<nodeEle.size(); i++) {
                    //roots
                    if (nodeEle[i]<NumEnt_[5]) {
                        dofEle[i*dofsPerNode+j] = nodeEle[i]+j*NumEnt_[5];
                        dmapEle[i] = nodeEle[i]+j*NumEnt_[5];
                    }
                    else if (NumEnt_[6] !=-1 &&nodeEle[i]>=NumEnt_[5] && nodeEle[i]<NumEnt_[6]+NumEnt_[5]) {
                        dofEle[i*dofsPerNode+j] = nodeEle[i]+j*NumEnt_[6]+(dofsPerNode-1)*NumEnt_[5];
                        dmapEle[i] = nodeEle[i]+j*NumEnt_[6]+(dofsPerNode-1)*NumEnt_[5];
                    }
                }
                dofsMaps[j] =   Xpetra::MapFactory<LO,GO,NO>::Build(Xpetra::UseTpetra,INVALID,dmapEle,0,nodesMap->getComm());
            }
        } else if (partition == 2) {
            FROSCH_ASSERT(false,"GDSWStar is not implemented yet!");
        }
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > tmpMap =   Xpetra::MapFactory<LO,GO,NO>::Build(Xpetra::UseTpetra,INVALID,dofEle,0,nodesMap->getComm());
        return tmpMap;
    }

    template <class SC,class LO,class GO,class NO>
    int  IPOUHarmonicCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                                   UN dofsPerNode,
                                                                   ConstXMapPtr nodesMap,
                                                                   ConstXMapPtrVecPtr dofsMaps,
                                                                   ConstXMultiVectorPtr nullSpaceBasis,
                                                                   GOVecPtr dirichletBoundaryDofs,
                                                                   ConstXMultiVectorPtr nodeList)
    {
        FROSCH_DETAILTIMER_START_LEVELID(buildCoarseSpaceTime,"IPOUHarmonicCoarseOperator::buildCoarseSpace");
        FROSCH_ASSERT(dofsMaps.size()==dofsPerNode,"dofsMaps.size()!=dofsPerNode");

        // Das könnte man noch ändern
        // Todo: Check the lengths of the vectors against NumberOfBlocks_
        return resetCoarseSpaceBlock(this->NumberOfBlocks_,dimension,dofsPerNode,nodesMap,dofsMaps,nullSpaceBasis,dirichletBoundaryDofs,nodeList);
    }

    template <class SC,class LO,class GO,class NO>
    int IPOUHarmonicCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                                  UNVecPtr dofsPerNodeVec,
                                                                  ConstXMapPtrVecPtr repeatedNodesMapVec,
                                                                  ConstXMapPtrVecPtr2D repeatedDofMapsVec,
                                                                  ConstXMultiVectorPtrVecPtr nullSpaceBasisVec,
                                                                  GOVecPtr2D dirichletBoundaryDofsVec,
                                                                  ConstXMultiVectorPtrVecPtr nodeListVec)
    {
        FROSCH_DETAILTIMER_START_LEVELID(buildCoarseSpaceTime,"IPOUHarmonicCoarseOperator::buildCoarseSpace");

        UN TotalNumberOfBlocks = repeatedNodesMapVec.size();

        FROSCH_ASSERT(dofsPerNodeVec.size()==TotalNumberOfBlocks,"dofsPerNodeVec.size()!=TotalNumberOfBlocks");
        FROSCH_ASSERT(repeatedDofMapsVec.size()==TotalNumberOfBlocks,"repeatedDofMapsVec.size()!=TotalNumberOfBlocks");
        FROSCH_ASSERT(nullSpaceBasisVec.size()==TotalNumberOfBlocks,"nullSpaceBasisVec.size()!=TotalNumberOfBlocks");
        FROSCH_ASSERT(dirichletBoundaryDofsVec.size()==TotalNumberOfBlocks,"dirichletBoundaryDofsVec.size()!=TotalNumberOfBlocks");
        FROSCH_ASSERT(nodeListVec.size()==TotalNumberOfBlocks,"nodeListVec.size()!=TotalNumberOfBlocks");

        // Todo: Move this to a function in HarmonicCoarseOperator at some point
        for (UN i=0; i<TotalNumberOfBlocks; i++) {
            FROSCH_ASSERT(!repeatedNodesMapVec[i].is_null(),"repeatedNodesMapVec[i].is_null()");
            FROSCH_ASSERT(!repeatedDofMapsVec[i].is_null(),"repeatedDofMapsVec[i].is_null()");
            FROSCH_ASSERT(!nullSpaceBasisVec[i].is_null(),"nullSpaceBasisVec[i].is_null()");

            resetCoarseSpaceBlock(this->NumberOfBlocks_,
                                  dimension,
                                  dofsPerNodeVec[i],
                                  repeatedNodesMapVec[i],
                                  repeatedDofMapsVec[i],
                                  nullSpaceBasisVec[i],
                                  dirichletBoundaryDofsVec[i],
                                  nodeListVec[i]);
        }
        return 0;
    }


    template <class SC,class LO,class GO,class NO>
    int IPOUHarmonicCoarseOperator<SC,LO,GO,NO>::resetCoarseSpaceBlock(UN blockId,
                                                                       UN dimension,
                                                                       UN dofsPerNode,
                                                                       ConstXMapPtr nodesMap,
                                                                       ConstXMapPtrVecPtr dofsMaps,
                                                                       ConstXMultiVectorPtr nullSpaceBasis,
                                                                       GOVecPtr dirichletBoundaryDofs,
                                                                       ConstXMultiVectorPtr nodeList)
    {
        FROSCH_DETAILTIMER_START_LEVELID(resetCoarseSpaceBlockTime,"IPOUHarmonicCoarseOperator::resetCoarseSpaceBlock");
        FROSCH_ASSERT(dofsMaps.size()==dofsPerNode,"dofsMaps.size()!=dofsPerNode");
        FROSCH_ASSERT(blockId<=this->NumberOfBlocks_,"Block does not exist yet and can therefore not be reset ("+to_string(blockId)+", "+to_string(this->NumberOfBlocks_)+").");

        // Process the parameter list
        stringstream blockIdStringstream;
        blockIdStringstream << blockId+1;
        string blockIdString = blockIdStringstream.str();
        RCP<ParameterList> coarseSpaceList = sublist(sublist(this->ParameterList_,"Blocks"),blockIdString.c_str());

        Verbosity verbosity = All;
        if (!coarseSpaceList->get("Verbosity","All").compare("None")) {
            verbosity = None;
        } else if (!coarseSpaceList->get("Verbosity","All").compare("All")) {
            verbosity = All;
        } else {
            FROSCH_ASSERT(false,"FROSch::IPOUHarmonicCoarseOperator: Specify a valid verbosity level.");
        }

        bool useForCoarseSpace = coarseSpaceList->get("Use For Coarse Space",true);

        if (useForCoarseSpace) {
            this->GammaDofs_.resize(this->GammaDofs_.size()+1);
            this->IDofs_.resize(this->IDofs_.size()+1);
            this->InterfaceCoarseSpaces_.resize(this->InterfaceCoarseSpaces_.size()+1);
            this->DofsMaps_.resize(this->DofsMaps_.size()+1);
            this->DofsPerNode_.resize(this->DofsPerNode_.size()+1);
            this->NumberOfBlocks_++;

            if (this->Verbose_) {
                cout
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << setw(89) << "-----------------------------------------------------------------------------------------"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| "
                << left << setw(74) << "IPOUHarmonicCoarseOperator " << right << setw(8) << "(Level " << setw(2) << this->LevelID_ << ")"
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
                << "| " << left << setw(41) << "Interface partition of unity type" << right
                << " | " << setw(41) << coarseSpaceList->sublist("InterfacePartitionOfUnity").get("Type","GDSW")
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << "| " << left << setw(41) << "Dimension of the null space" << right
                << " | " << setw(41) << nullSpaceBasis->getNumVectors()
                << " |"
                << "\n" << setw(FROSCH_OUTPUT_INDENT) << " "
                << setw(89) << "-----------------------------------------------------------------------------------------"
                << endl;
            }

            this->DofsMaps_[blockId] = dofsMaps;
            this->DofsPerNode_[blockId] = dofsPerNode;

            // Compute Interface Partition of Unity
            // AH: Can we get rid of the PartitionType_?
            InterfacePartitionOfUnityPtr interfacePartitionOfUnity;
            if (!coarseSpaceList->sublist("InterfacePartitionOfUnity").get("Type","GDSW").compare("GDSW")) {
                coarseSpaceList->sublist("InterfacePartitionOfUnity").sublist("GDSW").set("Test Unconnected Interface",this->ParameterList_->get("Test Unconnected Interface",true));
                interfacePartitionOfUnity = InterfacePartitionOfUnityPtr(new GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>(this->MpiComm_,this->SerialComm_,dimension,this->DofsPerNode_[blockId],nodesMap,this->DofsMaps_[blockId],sublist(sublist(coarseSpaceList,"InterfacePartitionOfUnity"),"GDSW"),verbosity,this->LevelID_));
                this->PartitionType_ = 0;
            } else if (!coarseSpaceList->sublist("InterfacePartitionOfUnity").get("Type","GDSW").compare("GDSWStar")) {
                coarseSpaceList->sublist("InterfacePartitionOfUnity").sublist("GDSWStar").set("Test Unconnected Interface",this->ParameterList_->get("Test Unconnected Interface",true));
                interfacePartitionOfUnity = InterfacePartitionOfUnityPtr(new GDSWStarInterfacePartitionOfUnity<SC,LO,GO,NO>(this->MpiComm_,this->SerialComm_,dimension,this->DofsPerNode_[blockId],nodesMap,this->DofsMaps_[blockId],sublist(sublist(coarseSpaceList,"InterfacePartitionOfUnity"),"GDSWStar"),verbosity,this->LevelID_));
                this->PartitionType_ = 2;
            } else if (!coarseSpaceList->sublist("InterfacePartitionOfUnity").get("Type","GDSW").compare("RGDSW")) {
                coarseSpaceList->sublist("InterfacePartitionOfUnity").sublist("RGDSW").set("Test Unconnected Interface",this->ParameterList_->get("Test Unconnected Interface",true));
                interfacePartitionOfUnity = InterfacePartitionOfUnityPtr(new RGDSWInterfacePartitionOfUnity<SC,LO,GO,NO>(this->MpiComm_,this->SerialComm_,dimension,this->DofsPerNode_[blockId],nodesMap,this->DofsMaps_[blockId],sublist(sublist(coarseSpaceList,"InterfacePartitionOfUnity"),"RGDSW"),verbosity,this->LevelID_));
                this->PartitionType_ = 1;
            } else {
                FROSCH_ASSERT(false,"InterfacePartitionOfUnity Type is unknown.");
            }

            // Extract the interface and the interior from the DDInterface stored in the Interface Partition of Unity object
            InterfaceEntityPtr interface = interfacePartitionOfUnity->getDDInterface()->getInterface()->getEntity(0);
            InterfaceEntityPtr interior = interfacePartitionOfUnity->getDDInterface()->getInterior()->getEntity(0);

            // Check whether the interface is empty. If so, we use a ConstantPartitionOfUnity instead, because we assume that the subdomains are decoupled.
            if (interface->getNumNodes()==0) {
                FROSCH_NOTIFICATION("FROSch::IPOUHarmonicCoarseOperator",this->Verbose_,"No interface found => A Constant Partition of Unity will be used instead.");
                coarseSpaceList->sublist("InterfacePartitionOfUnity").sublist("RGDSW").set("Test Unconnected Interface",this->ParameterList_->get("Test Unconnected Interface",true));
                PartitionOfUnity_ = PartitionOfUnityPtr(new ConstantPartitionOfUnity<SC,LO,GO,NO>(this->MpiComm_,
                                                                                                  this->SerialComm_,
                                                                                                  dimension,
                                                                                                  this->DofsPerNode_[blockId],
                                                                                                  nodesMap,
                                                                                                  this->DofsMaps_[blockId],
                                                                                                  sublist(sublist(coarseSpaceList,"PartitionOfUnity"),"Constant"),
                                                                                                  verbosity,
                                                                                                  this->LevelID_,
                                                                                                  interfacePartitionOfUnity->getDDInterfaceNonConst()));

                if (this->ParameterList_->get("Remove Dirichlet Nodes",true)) PartitionOfUnity_->removeDirichletNodes(dirichletBoundaryDofs());

                interface = interior;

                // Construct Interface and Interior index sets
                this->GammaDofs_[blockId] = LOVecPtr(this->DofsPerNode_[blockId]*interface->getNumNodes());
                this->IDofs_[blockId] = LOVecPtr(0);
                for (UN k=0; k<this->DofsPerNode_[blockId]; k++) {
                    for (UN i=0; i<interface->getNumNodes(); i++) {
                        this->GammaDofs_[blockId][interface->getGammaDofID(i,k)] = interface->getLocalDofID(i,k);
                    }
                }

                PartitionOfUnity_->computePartitionOfUnity(nodeList);
            } else {
                if (this->ParameterList_->get("Remove Dirichlet Nodes",true)) interfacePartitionOfUnity->removeDirichletNodes(dirichletBoundaryDofs(),nodeList);
                interfacePartitionOfUnity->sortInterface(this->K_,nodeList);

                // Construct Interface and Interior index sets
                this->GammaDofs_[blockId] = LOVecPtr(this->DofsPerNode_[blockId]*interface->getNumNodes());
                this->IDofs_[blockId] = LOVecPtr(this->DofsPerNode_[blockId]*interior->getNumNodes());
                for (UN k=0; k<this->DofsPerNode_[blockId]; k++) {
                    for (UN i=0; i<interface->getNumNodes(); i++) {
                        this->GammaDofs_[blockId][interface->getGammaDofID(i,k)] = interface->getLocalDofID(i,k);
                    }
                    for (UN i=0; i<interior->getNumNodes(); i++) {
                        this->IDofs_[blockId][interior->getGammaDofID(i,k)] = interior->getLocalDofID(i,k);
                    }
                }

                interfacePartitionOfUnity->computePartitionOfUnity(nodeList);
                PartitionOfUnity_ = interfacePartitionOfUnity;
            }

            // Construct local Interface nullspace basis (in case the interface was empty before, it was replaced by the interior. Therfore, this should be correct as well)
            ConstXMapPtr nullSpaceBasisMap = nullSpaceBasis->getMap();
            XMapPtr serialInterfaceMap = MapFactory<LO,GO,NO>::Build(nullSpaceBasisMap->lib(),this->GammaDofs_[blockId].size(),this->GammaDofs_[blockId].size(),0,this->SerialComm_);
            XMultiVectorPtr interfaceNullspaceBasis = MultiVectorFactory<SC,LO,GO,NO>::Build(serialInterfaceMap,nullSpaceBasis->getNumVectors());
            for (UN i=0; i<nullSpaceBasis->getNumVectors(); i++) {
                SCVecPtr interfaceNullspaceBasisData = interfaceNullspaceBasis->getDataNonConst(i);
                ConstSCVecPtr nullSpaceBasisData = nullSpaceBasis->getData(i);
                for (UN k=0; k<this->DofsPerNode_[blockId]; k++) {
                    for (UN j=0; j<interface->getNumNodes(); j++) {
                        interfaceNullspaceBasisData[interface->getGammaDofID(j,k)] = nullSpaceBasisData[nullSpaceBasisMap->getLocalElement(interface->getGlobalDofID(j,k))];
                    }
                }
            }

            PartitionOfUnity_->assembledPartitionOfUnityMaps();

            // This is setup for the multilevel variant. Can this be moved inside the following if block?
            this->CoarseDofsPerNode_ = nullSpaceBasis->getNumVectors();

            this->KRowMap_ = PartitionOfUnity_->getAssembledPartitionOfUnityMap();
            NumEnt_ = interfacePartitionOfUnity->getDDInterface()->getNumEnt();

            if (!this->DistributionList_->get("Type","linear").compare("ZoltanDual")) {
                FROSCH_ASSERT(this->NumberOfBlocks_==1,"Distribution Type ZoltanDual only works for one Block");
                Teuchos::RCP<DDInterface<SC,LO,GO,NO> > theInterface =Teuchos::rcp_const_cast<DDInterface<SC,LO,GO,NO> >(interfacePartitionOfUnity->getDDInterface());
                this->buildGlobalGraph(theInterface);
                int dim = dimension;
                sublist(this->ParameterList_,"CoarseSolver")->set("Dimension",dim);
            }

            // Build local basis
            LocalPartitionOfUnityBasis_ = LocalPartitionOfUnityBasisPtr(new LocalPartitionOfUnityBasis<SC,LO,GO,NO>(this->MpiComm_,this->SerialComm_,this->DofsPerNode_[blockId],sublist(coarseSpaceList,"LocalPartitionOfUnityBasis"),interfaceNullspaceBasis.getConst(),PartitionOfUnity_->getLocalPartitionOfUnity(),PartitionOfUnity_->getPartitionOfUnityMaps())); // sublist(coarseSpaceList,"LocalPartitionOfUnityBasis") testen

            LocalPartitionOfUnityBasis_->buildLocalPartitionOfUnityBasis();

            this->InterfaceCoarseSpaces_[blockId] = LocalPartitionOfUnityBasis_->getLocalPartitionOfUnitySpace();
            FROSCH_NOTIFICATION("FROSch::IPOUHarmonicCoarseOperator",this->Verbose_,"Need to build block coarse sizes for use in MueLu nullspace. This is not performed here yet.");
            //if (this->Verbose_) {RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout)); this->MVPhiGamma_[blockId]->describe(*fancy,VERB_EXTREME);}
        }

        return 0;
    }

}

#endif
