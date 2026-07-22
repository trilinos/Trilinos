// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_SCHWARZPRECONDITIONER_DECL_HPP
#define _FROSCH_SCHWARZPRECONDITIONER_DECL_HPP

#include <Xpetra_Operator.hpp>
#include <Xpetra_Matrix_fwd.hpp>

#include <ShyLU_DDFROSch_config.h>

#include <FROSch_Types.h>

#include <FROSch_SchwarzOperators_fwd.hpp>

namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class SchwarzPreconditioner : public Operator<SC,LO,GO,NO> {

    protected:

        using CommPtr                             = RCP<const Comm<int> >;

        using XMap                                = Map<LO,GO,NO>;
        using XMapPtr                             = RCP<XMap>;
        using ConstXMapPtr                        = RCP<const XMap>;
        using XMapPtrVecPtr                       = ArrayRCP<XMapPtr>;
        using ConstXMapPtrVecPtr                  = ArrayRCP<ConstXMapPtr>;
        using XMapPtrVecPtr2D                     = ArrayRCP<XMapPtrVecPtr>;
        using ConstXMapPtrVecPtr2D                = ArrayRCP<ConstXMapPtrVecPtr>;

        using XMatrix                             = Matrix<SC,LO,GO,NO>;
        using XMatrixPtr                          = RCP<XMatrix>;
        using ConstXMatrixPtr                     = RCP<const XMatrix>;

        using XMultiVector                        = MultiVector<SC,LO,GO,NO>;
        using XMultiVectorPtr                     = RCP<XMultiVector>;
        using ConstXMultiVectorPtr                = RCP<const XMultiVector>;
        using XMultiVectorPtrVecPtr               = ArrayRCP<XMultiVectorPtr>;
        using ConstXMultiVectorPtrVecPtr          = ArrayRCP<ConstXMultiVectorPtr>;

        using ParameterListPtr                    = RCP<ParameterList>;

        using SumOperatorPtr                      = RCP<SumOperator<SC,LO,GO,NO> >;
        using MultiplicativeOperatorPtr           = RCP<MultiplicativeOperator<SC,LO,GO,NO> >;
        using OverlappingOperatorPtr              = RCP<OverlappingOperator<SC,LO,GO,NO> >;
        using AlgebraicOverlappingOperatorPtr     = RCP<AlgebraicOverlappingOperator<SC,LO,GO,NO> >;
        using CoarseOperatorPtr                   = RCP<CoarseOperator<SC,LO,GO,NO> >;
        using GDSWCoarseOperatorPtr               = RCP<GDSWCoarseOperator<SC,LO,GO,NO> >;
        using RGDSWCoarseOperatorPtr              = RCP<RGDSWCoarseOperator<SC,LO,GO,NO> >;
        using IPOUHarmonicCoarseOperatorPtr       = RCP<IPOUHarmonicCoarseOperator<SC,LO,GO,NO> >;

        using DofOrderingVecPtr                   = ArrayRCP<DofOrdering>;

        using UN                                  = unsigned;
        using ConstUN                             = const UN;
        using UNVecPtr                            = ArrayRCP<UN>;

        using LOVecPtr                            = ArrayRCP<LO>;

        using GOVec                               = Array<GO>;
        using GOVec2D                             = Array<GOVec>;
        using GOVecPtr                            = ArrayRCP<GO>;
        using GOVecPtr2D                          = ArrayRCP<GOVecPtr>;

        using SCVecPtr                            = ArrayRCP<SC>;

    public:

        SchwarzPreconditioner(ParameterListPtr parameterList,
                              CommPtr comm);

        virtual ~SchwarzPreconditioner();

        virtual int initialize(bool useDefaultParameters = true) = 0;

        virtual int compute() = 0;

        // Y = alpha * A^mode * X + beta * Y
        virtual void apply(const XMultiVector &X,
                           XMultiVector &Y,
                           ETransp mode=NO_TRANS,
                           SC alpha=ScalarTraits<SC>::one(),
                           SC beta=ScalarTraits<SC>::zero()) const = 0;

        virtual const ConstXMapPtr getDomainMap() const = 0;

        virtual const ConstXMapPtr getRangeMap() const = 0;

        virtual void describe(FancyOStream &out,
                              const EVerbosityLevel verbLevel=Describable::verbLevel_default) const = 0;

        virtual string description() const = 0;

        bool isInitialized() const;

        bool isComputed() const;

        virtual void residual(const XMultiVector & X,
                              const XMultiVector & B,
                              XMultiVector& R) const;
    protected:

        CommPtr MpiComm_;

        ParameterListPtr ParameterList_;

        bool UseTranspose_ = false;
        bool IsInitialized_ = false;
        bool IsComputed_ = false;
        bool Verbose_ = false;

        ConstUN LevelID_ = 1;
    };

}

#endif
