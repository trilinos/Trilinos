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

#ifndef _FROSCH_CONSTANTPARTITIONOFUNITY_DECL_HPP
#define _FROSCH_CONSTANTPARTITIONOFUNITY_DECL_HPP

#include <FROSch_PartitionOfUnity_def.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = KokkosClassic::DefaultNode::DefaultNodeType>
    class ConstantPartitionOfUnity : public PartitionOfUnity<SC,LO,GO,NO> {

    protected:

        using CommPtr                       = typename PartitionOfUnity<SC,LO,GO,NO>::CommPtr;

        using XMapPtr                       = typename PartitionOfUnity<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr                  = typename PartitionOfUnity<SC,LO,GO,NO>::ConstXMapPtr;
        using XMapPtrVecPtr                 = typename PartitionOfUnity<SC,LO,GO,NO>::XMapPtrVecPtr;
        using ConstXMapPtrVecPtr            = typename PartitionOfUnity<SC,LO,GO,NO>::ConstXMapPtrVecPtr;

        using XMultiVectorPtr               = typename PartitionOfUnity<SC,LO,GO,NO>::XMultiVectorPtr;
        using ConstXMultiVectorPtr          = typename PartitionOfUnity<SC,LO,GO,NO>::ConstXMultiVectorPtr;
        using XMultiVectorPtrVecPtr         = typename PartitionOfUnity<SC,LO,GO,NO>::XMultiVectorPtrVecPtr;
        using ConstXMultiVectorPtrVecPtr    = typename PartitionOfUnity<SC,LO,GO,NO>::ConstXMultiVectorPtrVecPtr;

        using ParameterListPtr              = typename PartitionOfUnity<SC,LO,GO,NO>::ParameterListPtr;

        using DDInterfacePtr                = typename PartitionOfUnity<SC,LO,GO,NO>::DDInterfacePtr;

        using EntitySetPtr                  = typename PartitionOfUnity<SC,LO,GO,NO>::EntitySetPtr;

        using UN                            = typename PartitionOfUnity<SC,LO,GO,NO>::UN;

        using LOVec                         = typename PartitionOfUnity<SC,LO,GO,NO>::LOVec;

        using GOVec                         = typename PartitionOfUnity<SC,LO,GO,NO>::GOVec;
        using GOVecView                     = typename PartitionOfUnity<SC,LO,GO,NO>::GOVecView;

        using SCVec                         = typename PartitionOfUnity<SC,LO,GO,NO>::SCVec;

    public:

        ConstantPartitionOfUnity(CommPtr mpiComm,
                                 CommPtr serialComm,
                                 UN dimension,
                                 UN dofsPerNode,
                                 ConstXMapPtr nodesMap,
                                 ConstXMapPtrVecPtr dofsMaps,
                                 ParameterListPtr parameterList,
                                 Verbosity verbosity = All,
                                 UN levelID = 1,
                                 DDInterfacePtr ddInterface = null);

        virtual ~ConstantPartitionOfUnity();

        virtual int removeDirichletNodes(GOVecView dirichletBoundaryDofs,
                                         ConstXMultiVectorPtr nodeList = null);

        virtual int computePartitionOfUnity(ConstXMultiVectorPtr nodeList = null);

    protected:

        DDInterfacePtr DDInterface_;

        bool UseVolumes_ = false;

        EntitySetPtr Volumes_;
    };

}

#endif
