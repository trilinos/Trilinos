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

#ifndef _FROSCH_PARTITIONOFUNITYBASIS_DECL_hpp
#define _FROSCH_PARTITIONOFUNITYBASIS_DECL_hpp

#include <Xpetra_Operator.hpp>
#include <Xpetra_MapFactory_fwd.hpp>

#include <Teuchos_ScalarTraits.hpp>
#include <FROSch_CoarseSpace_def.hpp>

#include "FROSch_Tools_def.hpp"


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = KokkosClassic::DefaultNode::DefaultNodeType>
    class LocalPartitionOfUnityBasis {

    protected:

        using CommPtr                       = RCP<const Comm<int> >;

        using XMap                          = Map<LO,GO,NO>;
        using XMapPtr                       = RCP<XMap>;
        using ConstXMapPtr                  = RCP<const XMap>;
        using XMapPtrVecPtr                 = ArrayRCP<XMapPtr>;
        using ConstXMapPtrVecPtr            = ArrayRCP<ConstXMapPtr>;

        using XMultiVector                  = MultiVector<SC,LO,GO,NO>;
        using ConstXMultiVector             = const MultiVector<SC,LO,GO,NO>;
        using XMultiVectorPtr               = RCP<XMultiVector>;
        using ConstXMultiVectorPtr          = RCP<ConstXMultiVector>;
        using XMultiVectorPtrVecPtr         = ArrayRCP<XMultiVectorPtr>;
        using ConstXMultiVectorPtrVecPtr    = ArrayRCP<ConstXMultiVectorPtr>;
        using XMultiVectorPtrVecPtr2D       = ArrayRCP<XMultiVectorPtrVecPtr>;

        using ParameterListPtr              = RCP<ParameterList>;

        using CoarseSpacePtr                = RCP<CoarseSpace<SC,LO,GO,NO> >;

        using UN                            = unsigned;
        using UNVecPtr                      = ArrayRCP<UN>;

        using LOVecPtr                      = ArrayRCP<LO>;
        using LOVecPtr2D                    = ArrayRCP<LOVecPtr>;

        using BoolVecPtr                    = ArrayRCP<bool>;
        using BoolVecPtr2D                  = ArrayRCP<BoolVecPtr>;

    public:

        LocalPartitionOfUnityBasis(CommPtr mpiComm,
                                   CommPtr serialComm,
                                   UN dofsPerNode,
                                   ParameterListPtr parameterList,
                                   ConstXMultiVectorPtr nullSpaceBasis = XMultiVectorPtr(),
                                   ConstXMultiVectorPtrVecPtr partitionOfUnity = XMultiVectorPtrVecPtr(),
                                   ConstXMapPtrVecPtr partitionOfUnityMaps = XMapPtrVecPtr());

//        virtual ~LocalPartitionOfUnityBasis();

        int addPartitionOfUnity(ConstXMultiVectorPtrVecPtr partitionOfUnity,
                                ConstXMapPtrVecPtr partitionOfUnityMaps);

        int addGlobalBasis(ConstXMultiVectorPtr nullSpaceBasis);

        int buildLocalPartitionOfUnityBasis();

        XMultiVectorPtrVecPtr getPartitionOfUnity() const;

        XMultiVectorPtr getNullspaceBasis() const;

        CoarseSpacePtr getLocalPartitionOfUnitySpace() const;

    protected:

        CommPtr MpiComm_;
        CommPtr SerialComm_;

        UN DofsPerNode_ = 1;

        ParameterListPtr ParameterList_;

        CoarseSpacePtr LocalPartitionOfUnitySpace_;

        ConstXMultiVectorPtrVecPtr PartitionOfUnity_;
        ConstXMultiVectorPtr NullspaceBasis_;

        ConstXMapPtrVecPtr PartitionOfUnityMaps_;

    };

}

#endif
