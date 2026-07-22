// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_PARTITIONOFUNITYBASIS_DEF_HPP
#define _FROSCH_PARTITIONOFUNITYBASIS_DEF_HPP

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
                FROSCH_ASSERT(PartitionOfUnityMaps_[i]->getLocalNumElements()>0,"PartitionOfUnityMaps_[i]->getLocalNumElements()==0");
                tmpBasis[i] = XMultiVectorPtrVecPtr(PartitionOfUnity_[i]->getNumVectors());
                for (UN j=0; j<PartitionOfUnity_[i]->getNumVectors(); j++) {
                    XMultiVectorPtr tmpBasisJ = MultiVectorFactory<SC,LO,GO,NO>::Build(nullspaceBasisMap,NullspaceBasis_->getNumVectors());
                    tmpBasisJ->elementWiseMultiply(ScalarTraits<SC>::one(),*PartitionOfUnity_[i]->getVector(j),*NullspaceBasis_,ScalarTraits<SC>::one());
                    tmpBasis[i][j] = tmpBasisJ;
                }
            } else {
                FROSCH_ASSERT(PartitionOfUnityMaps_[i]->getLocalNumElements()==0,"PartitionOfUnityMaps_[i]->getLocalNumElements()!=0");
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
