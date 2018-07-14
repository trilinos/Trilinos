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

#ifndef _FROSCH_IPOUHARMONICCOARSEOPERATOR_DEF_HPP
#define _FROSCH_IPOUHARMONICCOARSEOPERATOR_DEF_HPP

#include <FROSch_IPOUHarmonicCoarseOperator_decl.hpp>

namespace FROSch {
    
    template <class SC,class LO,class GO,class NO>
     IPOUHarmonicCoarseOperator<SC,LO,GO,NO>:: IPOUHarmonicCoarseOperator(CrsMatrixPtr k,
                                                                          ParameterListPtr parameterList) :
    HarmonicCoarseOperator<SC,LO,GO,NO> (k,parameterList),
    InterfaceCoarseSpace_ (new CoarseSpace<SC,LO,GO,NO>()),
    InterfacePartitionOfUnity_ (),
    LocalPartitionOfUnityBasis_ ()
    {
        
    }

    template <class SC,class LO,class GO,class NO>
    int IPOUHarmonicCoarseOperator<SC,LO,GO,NO>::initialize(UN dimension,
                                                            UN dofsPerNode,
                                                            MapPtr nodesMap,
                                                            MapPtrVecPtr dofsMaps,
                                                            MultiVectorPtr nullSpaceBasis,
                                                            MultiVectorPtr nodeList,
                                                            GOVecPtr dirichletBoundaryDofs)
    {
        int ret = buildCoarseSpace(dimension,dofsPerNode,nodesMap,dofsMaps,nullSpaceBasis,dirichletBoundaryDofs,nodeList);
        this->IsInitialized_ = true;
        this->IsComputed_ = false;
        return ret;
    }
    
    template <class SC,class LO,class GO,class NO>
    void IPOUHarmonicCoarseOperator<SC,LO,GO,NO>::describe(Teuchos::FancyOStream &out,
                                                            const Teuchos::EVerbosityLevel verbLevel) const
    {
        FROSCH_ASSERT(0!=0,"describe() has be implemented properly...");
    }
    
    template <class SC,class LO,class GO,class NO>
    std::string IPOUHarmonicCoarseOperator<SC,LO,GO,NO>::description() const
    {
        return "Interface Partition of Unity Coarse Operator";
    }
    
    template <class SC,class LO,class GO,class NO>
    int  IPOUHarmonicCoarseOperator<SC,LO,GO,NO>::buildCoarseSpace(UN dimension,
                                                                   UN dofsPerNode,
                                                                   MapPtr nodesMap,
                                                                   MapPtrVecPtr dofsMaps,
                                                                   MultiVectorPtr nullSpaceBasis,
                                                                   GOVecPtr dirichletBoundaryDofs,
                                                                   MultiVectorPtr nodeList)
    {
        FROSCH_ASSERT(dofsMaps.size()==dofsPerNode,"dofsMaps.size()!=dofsPerNode");
        
        // Das könnte man noch ändern
        // LÄNGEN NOCHMAL GEGEN NumberOfBlocks_ checken!!!
        this->GammaDofs_.resize(this->GammaDofs_.size()+1);
        this->IDofs_.resize(this->IDofs_.size()+1);
        this->BlockCoarseMaps_.resize(this->BlockCoarseMaps_.size()+1);
        this->MVPhiGamma_.resize(this->MVPhiGamma_.size()+1);
        this->DofsMaps_.resize(this->DofsMaps_.size()+1);
        this->DofsPerNode_.resize(this->DofsPerNode_.size()+1);
        
        this->NumberOfBlocks_++;
        
        return resetCoarseSpaceBlock(this->NumberOfBlocks_-1,dimension,dofsPerNode,nodesMap,dofsMaps,nullSpaceBasis,dirichletBoundaryDofs,nodeList);
    }
    
    template <class SC,class LO,class GO,class NO>
    int  IPOUHarmonicCoarseOperator<SC,LO,GO,NO>::resetCoarseSpaceBlock(UN blockId,
                                                                        UN dimension,
                                                                        UN dofsPerNode,
                                                                        MapPtr nodesMap,
                                                                        MapPtrVecPtr dofsMaps,
                                                                        MultiVectorPtr nullSpaceBasis,
                                                                        GOVecPtr dirichletBoundaryDofs,
                                                                        MultiVectorPtr nodeList)
    {
        FROSCH_ASSERT(dofsMaps.size()==dofsPerNode,"dofsMaps.size()!=dofsPerNode");
        FROSCH_ASSERT(blockId<this->NumberOfBlocks_,"Block does not exist yet and can therefore not be reset.");
        
        // Process the parameter list
        std::stringstream blockIdStringstream;
        blockIdStringstream << blockId+1;
        std::string blockIdString = blockIdStringstream.str();
        Teuchos::RCP<Teuchos::ParameterList> coarseSpaceList = sublist(sublist(this->ParameterList_,"Blocks"),blockIdString.c_str());
        
        bool useForCoarseSpace = coarseSpaceList->get("Use For Coarse Space",true);
        
        if (useForCoarseSpace) {
            this->DofsMaps_[blockId] = dofsMaps;
            this->DofsPerNode_[blockId] = dofsPerNode;
            
            // Compute Interface Partition of Unity
            if (!coarseSpaceList->sublist("InterfacePartitionOfUnity").get("Type","GDSW").compare("GDSW")) {
                InterfacePartitionOfUnity_ = InterfacePartitionOfUnityPtr(new GDSWInterfacePartitionOfUnity<SC,LO,GO,NO>(this->MpiComm_,this->SerialComm_,dimension,this->DofsPerNode_[blockId],nodesMap,this->DofsMaps_[blockId],sublist(sublist(coarseSpaceList,"InterfacePartitionOfUnity"),"GDSW")));
            } else {
                FROSCH_ASSERT(0!=0,"InterfacePartitionOfUnity Type is unknown.");
            }            
            InterfacePartitionOfUnity_->removeDirichletNodes(dirichletBoundaryDofs(),nodeList);
            InterfacePartitionOfUnity_->sortInterface(this->K_,nodeList);
            InterfacePartitionOfUnity_->computePartitionOfUnity();
            
            InterfaceEntityPtr interface = InterfacePartitionOfUnity_->getDDInterface()->getInterface()->getEntity(0);
            InterfaceEntityPtr interior = InterfacePartitionOfUnity_->getDDInterface()->getInterior()->getEntity(0);
            
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
            
            // Construct local Interface nullspace basis
            MapPtr serialInterfaceMap = Xpetra::MapFactory<LO,GO,NO>::Build(nullSpaceBasis->getMap()->lib(),this->GammaDofs_[blockId].size(),this->GammaDofs_[blockId].size(),0,this->SerialComm_);
            MultiVectorPtr interfaceNullspaceBasis = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(serialInterfaceMap,nullSpaceBasis->getNumVectors());
            for (UN i=0; i<nullSpaceBasis->getNumVectors(); i++) {
                for (UN k=0; k<this->DofsPerNode_[blockId]; k++) {
                    for (UN j=0; j<interface->getNumNodes(); j++) {
                        interfaceNullspaceBasis->getDataNonConst(i)[interface->getGammaDofID(j,k)] = nullSpaceBasis->getData(i)[nullSpaceBasis->getMap()->getLocalElement(interface->getGlobalDofID(j,k))];
                    }
                }
            }
            
            // Build local basis
            LocalPartitionOfUnityBasis_ = LocalPartitionOfUnityBasisPtr(new LocalPartitionOfUnityBasis<SC,LO,GO,NO>(this->MpiComm_,this->SerialComm_,this->DofsPerNode_[blockId],sublist(coarseSpaceList,"LocalPartitionOfUnityBasis"),interfaceNullspaceBasis,InterfacePartitionOfUnity_->getLocalPartitionOfUnity(),InterfacePartitionOfUnity_->getPartitionOfUnityMaps())); // sublist(coarseSpaceList,"LocalPartitionOfUnityBasis") testen
            
            LocalPartitionOfUnityBasis_->buildLocalPartitionOfUnityBasis();
            
            InterfaceCoarseSpace_ = LocalPartitionOfUnityBasis_->getLocalPartitionOfUnitySpace();
            
            this->BlockCoarseMaps_[blockId] = InterfaceCoarseSpace_->getBasisMap();
            this->MVPhiGamma_[blockId] = InterfaceCoarseSpace_->getLocalBasis(); //if (this->Verbose_) {Teuchos::RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)); this->MVPhiGamma_[blockId]->describe(*fancy,Teuchos::VERB_EXTREME);}
        }
        
        return 0;
    }
}

#endif
