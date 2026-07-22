// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_PARTITIONOFUNITYBASIS_DECL_HPP
#define _FROSCH_PARTITIONOFUNITYBASIS_DECL_HPP

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
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
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
