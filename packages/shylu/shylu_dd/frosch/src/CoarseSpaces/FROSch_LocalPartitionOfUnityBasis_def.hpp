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

    using namespace Teuchos;
    using namespace Xpetra;

    template<class SC,class LO,class GO,class NO>
    LocalPartitionOfUnityBasis<SC,LO,GO,NO>::LocalPartitionOfUnityBasis(CommPtr mpiComm,
                                                                        CommPtr serialComm,
                                                                        UN dofsPerNode,
                                                                        ParameterListPtr parameterList,
                                                                        ConstXMultiVectorPtr nullSpaceBasis,
                                                                        ConstXMultiVectorPtrVecPtr partitionOfUnity,
                                                                        ConstXMapPtrVecPtr partitionOfUnityMaps) :
    MpiComm_ (mpiComm),
    SerialComm_ (serialComm),
    DofsPerNode_ (dofsPerNode),
    ParameterList_ (parameterList),
    PartitionOfUnity_ (partitionOfUnity),
    NullspaceBasis_ (nullSpaceBasis),
    PartitionOfUnityMaps_ (partitionOfUnityMaps)
    {

    }

    template<class SC,class LO,class GO,class NO>
    int LocalPartitionOfUnityBasis<SC,LO,GO,NO>::addPartitionOfUnity(ConstXMultiVectorPtrVecPtr partitionOfUnity,
                                                                     ConstXMapPtrVecPtr partitionOfUnityMaps)
    {
        PartitionOfUnity_ = partitionOfUnity;
        PartitionOfUnityMaps_ = partitionOfUnityMaps;
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int LocalPartitionOfUnityBasis<SC,LO,GO,NO>::addGlobalBasis(ConstXMultiVectorPtr nullSpaceBasis)
    {
        NullspaceBasis_ = nullSpaceBasis;
        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    int LocalPartitionOfUnityBasis<SC,LO,GO,NO>::buildLocalPartitionOfUnityBasis()
    {
        FROSCH_ASSERT(!NullspaceBasis_.is_null(),"Nullspace Basis is not set.");
        FROSCH_ASSERT(!PartitionOfUnity_.is_null(),"Partition Of Unity is not set.");
        FROSCH_ASSERT(!PartitionOfUnityMaps_.is_null(),"Partition Of Unity Map is not set.");

        LocalPartitionOfUnitySpace_ = CoarseSpacePtr(new CoarseSpace<SC,LO,GO,NO>(this->MpiComm_,this->SerialComm_));

        XMultiVectorPtrVecPtr2D tmpBasis(PartitionOfUnity_.size());
        ConstXMapPtr nullspaceBasisMap = NullspaceBasis_->getMap();
        for (UN i=0; i<PartitionOfUnity_.size(); i++) {
            if (!PartitionOfUnity_[i].is_null()) {
                FROSCH_ASSERT(PartitionOfUnityMaps_[i]->getNodeNumElements()>0,"PartitionOfUnityMaps_[i]->getNodeNumElements()==0");
                tmpBasis[i] = XMultiVectorPtrVecPtr(PartitionOfUnity_[i]->getNumVectors());
                for (UN j=0; j<PartitionOfUnity_[i]->getNumVectors(); j++) {
                    XMultiVectorPtr tmpBasisJ = MultiVectorFactory<SC,LO,GO,NO>::Build(nullspaceBasisMap,NullspaceBasis_->getNumVectors());
                    tmpBasisJ->elementWiseMultiply(ScalarTraits<SC>::one(),*PartitionOfUnity_[i]->getVector(j),*NullspaceBasis_,ScalarTraits<SC>::one());
                    tmpBasis[i][j] = tmpBasisJ;
                }
            } else {
                FROSCH_ASSERT(PartitionOfUnityMaps_[i]->getNodeNumElements()==0,"PartitionOfUnityMaps_[i]->getNodeNumElements()!=0");
            }
        }

        // Kann man das schöner machen?
        for (UN i=0; i<PartitionOfUnity_.size(); i++) {
            if (!PartitionOfUnityMaps_[i].is_null()) {
                if (!PartitionOfUnity_[i].is_null()) {
                    ConstXMapPtr partitionOfUnityMap_i = PartitionOfUnity_[i]->getMap();
                    for (UN j=0; j<NullspaceBasis_->getNumVectors(); j++) {
                        XMultiVectorPtr entityBasis = MultiVectorFactory<SC,LO,GO,NO >::Build(partitionOfUnityMap_i,PartitionOfUnity_[i]->getNumVectors());
                        entityBasis->scale(ScalarTraits<SC>::zero());
                        for (UN k=0; k<PartitionOfUnity_[i]->getNumVectors(); k++) {
                            if (j<tmpBasis[i][k]->getNumVectors()) {
                                entityBasis->getDataNonConst(k).deepCopy(tmpBasis[i][k]->getData(j)()); // Here, we copy data. Do we need to do this?
                            }
                        }
                        LocalPartitionOfUnitySpace_->addSubspace(PartitionOfUnityMaps_[i],null,entityBasis);
                    }
                } else {
                    for (UN j=0; j<NullspaceBasis_->getNumVectors(); j++) {
                        LocalPartitionOfUnitySpace_->addSubspace(PartitionOfUnityMaps_[i]);
                    }
                }
            } else {
                FROSCH_WARNING("FROSch::LocalPartitionOfUnityBasis",this->MpiComm_->getRank()==0,"PartitionOfUnityMaps_[i].is_null()");
            }
        }

        LocalPartitionOfUnitySpace_->assembleCoarseSpace();

        return 0;
    }

    template<class SC,class LO,class GO,class NO>
    typename LocalPartitionOfUnityBasis<SC,LO,GO,NO>::XMultiVectorPtrVecPtr LocalPartitionOfUnityBasis<SC,LO,GO,NO>::getPartitionOfUnity() const
    {
        FROSCH_ASSERT(!PartitionOfUnity_.is_null(),"Partition Of Unity is not set.");
        return PartitionOfUnity_;
    }

    template<class SC,class LO,class GO,class NO>
    typename LocalPartitionOfUnityBasis<SC,LO,GO,NO>::XMultiVectorPtr LocalPartitionOfUnityBasis<SC,LO,GO,NO>::getNullspaceBasis() const
    {
        FROSCH_ASSERT(!NullspaceBasis_.is_null(),"Nullspace Basis is not set.");
        return NullspaceBasis_;
    }

    template<class SC,class LO,class GO,class NO>
    typename LocalPartitionOfUnityBasis<SC,LO,GO,NO>::CoarseSpacePtr LocalPartitionOfUnityBasis<SC,LO,GO,NO>::getLocalPartitionOfUnitySpace() const
    {
        FROSCH_ASSERT(!LocalPartitionOfUnitySpace_.is_null(),"Local Partition Of Unity Space is not built yet.");
        return LocalPartitionOfUnitySpace_;
    }
}

#endif
