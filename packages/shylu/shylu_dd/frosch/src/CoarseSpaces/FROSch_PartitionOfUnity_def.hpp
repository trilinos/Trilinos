// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_PARTITIONOFUNITY_DEF_HPP
#define _FROSCH_PARTITIONOFUNITY_DEF_HPP

#include <FROSch_PartitionOfUnity_decl.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    PartitionOfUnity<SC,LO,GO,NO>::PartitionOfUnity(CommPtr mpiComm,
                                                    CommPtr serialComm,
                                                    UN dofsPerNode,
                                                    ConstXMapPtr nodesMap,
                                                    ConstXMapPtrVecPtr dofsMaps,
                                                    ParameterListPtr parameterList,
                                                    Verbosity verbosity,
                                                    UN levelID) :
    MpiComm_ (mpiComm),
    SerialComm_ (serialComm),
    ParameterList_ (parameterList),
    Verbose_ (MpiComm_->getRank() == 0),
    Verbosity_ (verbosity),
    LevelID_ (levelID)
    {

    }

    template <class SC,class LO,class GO,class NO>
    PartitionOfUnity<SC,LO,GO,NO>::~PartitionOfUnity()
    {

    }

    template <class SC,class LO,class GO,class NO>
    int PartitionOfUnity<SC,LO,GO,NO>::assembledPartitionOfUnityMaps()
    {
        if (!AssmbledPartitionOfUnityMap_.is_null()) {
            FROSCH_NOTIFICATION("FROSch::PartitionOfUnity",Verbosity_,"AssmbledPartitionOfUnityMap_ has already been assembled previously.");
        }
        LOVecPtr2D partMappings;
        AssmbledPartitionOfUnityMap_ = AssembleMaps(PartitionOfUnityMaps_(),partMappings);
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    typename PartitionOfUnity<SC,LO,GO,NO>::ConstXMultiVectorPtrVecPtr PartitionOfUnity<SC,LO,GO,NO>::getLocalPartitionOfUnity() const
    {
        return LocalPartitionOfUnity_;
    }

    template <class SC,class LO,class GO,class NO>
    typename PartitionOfUnity<SC,LO,GO,NO>::ConstXMapPtrVecPtr PartitionOfUnity<SC,LO,GO,NO>::getPartitionOfUnityMaps() const
    {
        return PartitionOfUnityMaps_;
    }

    template <class SC,class LO,class GO,class NO>
    typename PartitionOfUnity<SC,LO,GO,NO>::ConstXMapPtr PartitionOfUnity<SC,LO,GO,NO>::getAssembledPartitionOfUnityMap() const
    {
        return AssmbledPartitionOfUnityMap_;
    }
}

#endif
