// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_COARSESPACE_DECL_HPP
#define _FROSCH_COARSESPACE_DECL_HPP

//#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_Access.hpp>

#include<KokkosKernels_Utils.hpp>

#include <FROSch_Tools_def.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class CoarseSpace {

    protected:

        using CommPtr                   = RCP<const Comm<int> >;

        using XMap                      = Map<LO,GO,NO>;
        using XMapPtr                   = RCP<XMap>;
        using ConstXMapPtr              = RCP<const XMap>;
        using XMapPtrVec                = Array<XMapPtr>;
        using ConstXMapPtrVec           = Array<ConstXMapPtr>;

        using XMatrix                   = Matrix<SC,LO,GO,NO>;
        using XMatrixPtr                = RCP<XMatrix>;

        using XMultiVector              = MultiVector<SC,LO,GO,NO>;
        using XMultiVectorPtr           = RCP<XMultiVector>;
        using ConstXMultiVectorPtr      = RCP<const XMultiVector>;
        using ConstXMultiVectorPtrVec   = Array<ConstXMultiVectorPtr>;

        using ParameterListPtr          = RCP<ParameterList>;

        using UN                        = unsigned;
        using UNVec                     = Array<UN>;
        using ConstUNVecView            = ArrayView<const UN>;

        using LOVec                     = Array<LO>;
        using LOVecPtr                  = ArrayRCP<LO>;
        using ConstLOVecView            = ArrayView<const LO>;
        using LOVecPtr2D                = ArrayRCP<LOVecPtr>;

        using GOVec                     = Array<GO>;

        using SCVec                     = Array<SC>;
        using ConstSCVecPtr             = ArrayRCP<const SC>;

    public:

        CoarseSpace(CommPtr mpiComm,
                    CommPtr serialComm);

        int addSubspace(ConstXMapPtr subspaceBasisMap,
                        ConstXMapPtr subspaceBasisMapUnique = null,
                        ConstXMultiVectorPtr subspaceBasis = null,
                        UN offset = 0);

        int assembleCoarseSpace();

        int buildGlobalBasisMatrix(ConstXMapPtr rowMap,
                                   ConstXMapPtr rangeMap,
                                   ConstXMapPtr repeatedMap,
                                   SC tresholdDropping);

        int clearCoarseSpace();

        int zeroOutBasisVectors(ConstLOVecView zeros);

        bool hasUnassembledMaps() const;

        bool hasBasisMap() const;

        ConstXMapPtr getBasisMap() const;

        bool hasBasisMapUnique() const;

        ConstXMapPtr getBasisMapUnique() const;

        bool hasAssembledBasis() const;

        ConstXMultiVectorPtr getAssembledBasis() const;

        ConstUNVecView getLocalSubspaceSizes() const;

        bool hasGlobalBasisMatrix() const;

        XMatrixPtr getGlobalBasisMatrix() const;

    protected:

        CommPtr MpiComm_;
        CommPtr SerialComm_;

        ConstXMapPtrVec UnassembledBasesMaps_ = ConstXMapPtrVec(0);
        ConstXMapPtrVec UnassembledBasesMapsUnique_ = ConstXMapPtrVec(0);

        ConstXMultiVectorPtrVec UnassembledSubspaceBases_ = ConstXMultiVectorPtrVec(0);

        LOVec Offsets_ = LOVec(0);

        ConstXMapPtr AssembledBasisMap_;
        ConstXMapPtr AssembledBasisMapUnique_;

        XMultiVectorPtr AssembledBasis_;

        UNVec LocalSubspacesSizes_ = UNVec(0);

        XMatrixPtr GlobalBasisMatrix_;
    };

}

#endif
