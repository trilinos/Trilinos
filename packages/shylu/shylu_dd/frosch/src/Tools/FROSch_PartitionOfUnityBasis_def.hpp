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

#ifndef _FROSCH_PARTITIONOFUNITYBASIS_DEF_hpp
#define _FROSCH_PARTITIONOFUNITYBASIS_DEF_hpp

#include <FROSch_PartitionOfUnityBasis_decl.hpp>

namespace FROSch {
    
    template<class SC,class LO,class GO,class NO>
    PartitionOfUnityBasis<SC,LO,GO,NO>::PartitionOfUnityBasis(MapPtr &nodesMap,
                                                              MapPtrVecPtr &dofsMaps) :
    DofsPerNode_ (dofsMaps.size()),
    NodesMap_ (nodesMap),
    DofsMaps_ (dofsMaps),
    PartitionOfUnity_ (),
    GlobalBasis_ (),
    Basis_ ()
    {
        
    }
    
    template<class SC,class LO,class GO,class NO>
    PartitionOfUnityBasis<SC,LO,GO,NO>::PartitionOfUnityBasis(MapPtr &nodesMap,
                                                              MapPtrVecPtr &dofsMaps,
                                                              MultiVectorPtr &partitionOfUnity,
                                                              MultiVectorPtr &globalBasis) :
    DofsPerNode_ (dofsMaps.size()),
    NodesMap_ (nodesMap),
    DofsMaps_ (dofsMaps),
    PartitionOfUnity_ (partitionOfUnity),
    GlobalBasis_ (globalBasis),
    Basis_ ()
    {
        FROSCH_ASSERT(partitionOfUnity->getMap()->isSameAs(*NodesMap_),"The Maps are not compatible.");
        PartitionOfUnity_ = partitionOfUnity;
        FROSCH_ASSERT(globalBasis->getMap()->isSameAs(*NodesMap_),"The Maps are not compatible.");
        GlobalBasis_ = globalBasis;
    }
    
    template<class SC,class LO,class GO,class NO>
    int PartitionOfUnityBasis<SC,LO,GO,NO>::addPartitionOfUnity(MultiVectorPtr &partitionOfUnity)
    {
        FROSCH_ASSERT(partitionOfUnity->getMap()->isSameAs(*NodesMap_),"The Maps are not compatible.");
        PartitionOfUnity_ = partitionOfUnity;
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    int PartitionOfUnityBasis<SC,LO,GO,NO>::addGlobalBasis(MultiVectorPtr &globalBasis)
    {
        FROSCH_ASSERT(globalBasis->getMap()->isSameAs(*NodesMap_),"The Maps are not compatible.");
        GlobalBasis_ = globalBasis;
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    int PartitionOfUnityBasis<SC,LO,GO,NO>::buildPartitionOfUnityBasis()
    {
        FROSCH_ASSERT(!PartitionOfUnity_.is_null(),"PartitionOfUnity is not set.");
        FROSCH_ASSERT(!GlobalBasis_.is_null(),"GlobalBasis is not set.");
        
        MultiVectorPtrVecPtr tmpBasis(PartitionOfUnity_->getNumVectors());
        for (UN i=0; i<PartitionOfUnity_->getNumVectors(); i++) {
            tmpBasis[i] = Xpetra::MultiVectorFactory<LO,GO,NO>::Build(NodesMap_,GlobalBasis_->getNumVectors());
            tmpBasis[i]->elementWiseMultiply(1.0,PartitionOfUnity_->getVector(i),GlobalBasis_,1.0);
            
            // Hier muss noch orthogonalisiert werden
        }
        
    }
    
    template<class SC,class LO,class GO,class NO>
    typename PartitionOfUnityBasis<SC,LO,GO,NO>::ConstMultiVectorPtr PartitionOfUnityBasis<SC,LO,GO,NO>::getPartitionOfUnity() const
    {
        return PartitionOfUnity_.getConst();
    }
    
    template<class SC,class LO,class GO,class NO>
    typename PartitionOfUnityBasis<SC,LO,GO,NO>::ConstMultiVectorPtr PartitionOfUnityBasis<SC,LO,GO,NO>::getGlobalBasis() const
    {
        return GlobalBasis_.getConst();
    }
    
    template<class SC,class LO,class GO,class NO>
    typename PartitionOfUnityBasis<SC,LO,GO,NO>::ConstMultiVectorPtr PartitionOfUnityBasis<SC,LO,GO,NO>::getBasis() const
    {
        FROSCH_ASSERT(!Basis_.is_null(),"Basis is not built yet.");
        return Basis_.getConst();
    }
}

#endif
