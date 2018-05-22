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

#include <FROSch_LocalPartitionOfUnityBasis_decl.hpp>

namespace FROSch {
    
    template<class SC,class LO,class GO,class NO>
    LocalPartitionOfUnityBasis<SC,LO,GO,NO>::LocalPartitionOfUnityBasis(MapPtr &nodesMap,
                                                              MapPtrVecPtr &dofsMaps) :
    DofsPerNode_ (dofsMaps.size()),
    NodesMap_ (nodesMap),
    DofsMaps_ (dofsMaps),
    PartitionOfUnity_ (),
    NullSpaceBasis_ (),
    LocalBasis_ ()
    {
        
    }
    
    template<class SC,class LO,class GO,class NO>
    LocalPartitionOfUnityBasis<SC,LO,GO,NO>::LocalPartitionOfUnityBasis(MapPtr &nodesMap,
                                                              MapPtrVecPtr &dofsMaps,
                                                              MultiVectorPtr &partitionOfUnity,
                                                              MultiVectorPtr &nullspacebasis) :
    DofsPerNode_ (dofsMaps.size()),
    NodesMap_ (nodesMap),
    DofsMaps_ (dofsMaps),
    PartitionOfUnity_ (partitionOfUnity),
    NullSpaceBasis_ (nullspacebasis),
    LocalBasis_ ()
    {
        FROSCH_ASSERT(partitionOfUnity->getMap()->isSameAs(*NodesMap_),"The Maps are not compatible.");
    }
    
    template<class SC,class LO,class GO,class NO>
    int LocalPartitionOfUnityBasis<SC,LO,GO,NO>::addPartitionOfUnity(MultiVectorPtr &partitionOfUnity)
    {
        FROSCH_ASSERT(partitionOfUnity->getMap()->isSameAs(*NodesMap_),"The Maps are not compatible.");
        PartitionOfUnity_ = partitionOfUnity;
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    int LocalPartitionOfUnityBasis<SC,LO,GO,NO>::addGlobalBasis(MultiVectorPtr &nullspacebasis)
    {
        NullSpaceBasis_ = nullspacebasis;
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    int LocalPartitionOfUnityBasis<SC,LO,GO,NO>::buildLocalPartitionOfUnityBasis()
    {
        FROSCH_ASSERT(!PartitionOfUnity_.is_null(),"PartitionOfUnity is not set.");
        FROSCH_ASSERT(!NullSpaceBasis_.is_null(),"GlobalBasis is not set.");
        
        LocalBasis_ = ConstMultiVectorPtrVecPtr(PartitionOfUnity_->getNumVectors());
        
        MultiVectorPtr tmpBasis;
        for (UN i=0; i<PartitionOfUnity_->getNumVectors(); i++) {
            tmpBasis = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(NullSpaceBasis_->getMap(),NullSpaceBasis_->getNumVectors());
            for (UN j=0; j<PartitionOfUnity_->getLocalLength(); j++) {
                for (UN k=0; k<NullSpaceBasis_->getNumVectors(); k++) {
                    for (UN l=0; l<DofsPerNode_; l++) {
                        tmpBasis->getDataNonConst(k)[DofsMaps_[l]->getGlobalElement(j)] = PartitionOfUnity_->getData(i)[j]*NullSpaceBasis_->getData(k)[DofsMaps_[l]->getGlobalElement(j)];
                    }
                }
            }
            
            // Orthogonalization
            LocalBasis_[i] = ModifiedGramSchmidt(tmpBasis.getConst()).getConst();
        }
        return 0;
    }
    
    template<class SC,class LO,class GO,class NO>
    typename LocalPartitionOfUnityBasis<SC,LO,GO,NO>::ConstMultiVectorPtr LocalPartitionOfUnityBasis<SC,LO,GO,NO>::getPartitionOfUnity() const
    {
        FROSCH_ASSERT(!PartitionOfUnity_.is_null(),"Basis is not built yet.");
        return PartitionOfUnity_.getConst();
    }
    
    template<class SC,class LO,class GO,class NO>
    typename LocalPartitionOfUnityBasis<SC,LO,GO,NO>::ConstMultiVectorPtr LocalPartitionOfUnityBasis<SC,LO,GO,NO>::getGlobalBasis() const
    {
        FROSCH_ASSERT(!NullSpaceBasis_.is_null(),"Basis is not built yet.");
        return NullSpaceBasis_.getConst();
    }
    
    template<class SC,class LO,class GO,class NO>
    typename LocalPartitionOfUnityBasis<SC,LO,GO,NO>::ConstMultiVectorPtrVecPtr LocalPartitionOfUnityBasis<SC,LO,GO,NO>::getBasis() const
    {
        FROSCH_ASSERT(!LocalBasis_.is_null(),"Basis is not built yet.");
        return LocalBasis_;
    }
}

#endif
