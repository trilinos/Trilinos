// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_SCHWARZOPERATOR_DECL_HPP
#define _FROSCH_SCHWARZOPERATOR_DECL_HPP

#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_TripleMatrixMultiply.hpp>
#include <Xpetra_Export.hpp>

#include <Teuchos_DefaultSerialComm.hpp>

#include <Teuchos_SerialQRDenseSolver.hpp>
#include <Teuchos_TwoDArray.hpp>

#include <ShyLU_DDFROSch_config.h>

#include <FROSch_DDInterface_def.hpp>
#include <FROSch_EntitySet_def.hpp>

#include <FROSch_CoarseSpace_def.hpp>
#include <FROSch_InterfacePartitionOfUnity_def.hpp>
#include <FROSch_LocalPartitionOfUnityBasis_def.hpp>

#include <FROSch_SolverFactory_def.hpp>

// TODO: Auf const überprüfen
// TODO: #ifndef überprüfen ??????


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class SchwarzOperator : public Operator<SC,LO,GO,NO> {

    protected:

        using CommPtr                           = RCP<const Comm<int> >;

        using XMap                              = Map<LO,GO,NO>;
        using XMapPtr                           = RCP<XMap>;
        using ConstXMapPtr                      = RCP<const XMap>;
        using XMapPtrVecPtr                     = ArrayRCP<XMapPtr>;
        using ConstXMapPtrVecPtr                = ArrayRCP<ConstXMapPtr>;
        using XMapPtrVecPtr2D                   = ArrayRCP<XMapPtrVecPtr>;
        using ConstXMapPtrVecPtr2D              = ArrayRCP<ConstXMapPtrVecPtr>;

        using XMatrix                           = Matrix<SC,LO,GO,NO>;
        using XMatrixPtr                        = RCP<XMatrix>;
        using ConstXMatrixPtr                   = RCP<const XMatrix>;

        using XCrsGraph                         = CrsGraph<LO,GO,NO>;
        using GraphPtr                          = RCP<XCrsGraph>;
        using ConstXCrsGraphPtr                 = RCP<const XCrsGraph>;

        using XMultiVector                      = MultiVector<SC,LO,GO,NO>;
        using XMultiVectorPtr                   = RCP<XMultiVector>;
        using ConstXMultiVectorPtr              = RCP<const XMultiVector>;
        using XMultiVectorPtrVecPtr             = ArrayRCP<XMultiVectorPtr>;
        using ConstXMultiVectorPtrVecPtr        = ArrayRCP<ConstXMultiVectorPtr>;

        using XImport                           = Import<LO,GO,NO>;
        using XImportPtr                        = RCP<XImport>;
        using XImportPtrVecPtr                  = ArrayRCP<XImportPtr>;

        using XExport                           = Export<LO,GO,NO>;
        using XExportPtr                        = RCP<XExport>;
        using XExportPtrVecPtr                  = ArrayRCP<XExportPtr>;

        using ParameterListPtr                  = RCP<ParameterList>;

        using TSerialDenseMatrixPtr             = RCP<SerialDenseMatrix<LO,SC> >;

        using TSerialQRDenseSolverPtr           = RCP<SerialQRDenseSolver<LO,SC> >;

        using DDInterfacePtr                    = RCP<DDInterface<SC,LO,GO,NO> >;

        using EntitySetPtr                      = RCP<EntitySet<SC,LO,GO,NO> >;
        using EntitySetConstPtr                 = const EntitySetPtr;
        using EntitySetPtrVecPtr                = ArrayRCP<EntitySetPtr>;
        using EntitySetPtrConstVecPtr           = const EntitySetPtrVecPtr;

        using CoarseSpacePtr                    = RCP<CoarseSpace<SC,LO,GO,NO> >;
        using CoarseSpacePtrVecPtr              = ArrayRCP<CoarseSpacePtr>;

        using InterfaceEntityPtr                = RCP<InterfaceEntity<SC,LO,GO,NO> >;
        using InterfaceEntityPtrVec             = Array<InterfaceEntityPtr>;
        using InterfaceEntityPtrVecPtr          = ArrayRCP<InterfaceEntityPtr>;

        using PartitionOfUnityPtr               = RCP<PartitionOfUnity<SC,LO,GO,NO> >;
        using InterfacePartitionOfUnityPtr      = RCP<InterfacePartitionOfUnity<SC,LO,GO,NO> >;

        using LocalPartitionOfUnityBasisPtr     = RCP<LocalPartitionOfUnityBasis<SC,LO,GO,NO> >;

        using SchwarzOperatorPtr                = RCP<SchwarzOperator<SC,LO,GO,NO> >;
        using SchwarzOperatorPtrVec             = Array<SchwarzOperatorPtr>;
        using SchwarzOperatorPtrVecPtr          = ArrayRCP<SchwarzOperatorPtr>;

        using SolverPtr                         = RCP<Solver<SC,LO,GO,NO> >;
        using SolverFactoryPtr                  = RCP<SolverFactory<SC,LO,GO,NO> >;

        using DofOrderingVecPtr                 = ArrayRCP<DofOrdering>;

        using UN                                = unsigned;
        using ConstUN                           = const UN;
        using UNVec                             = Array<UN>;
        using UNVecPtr                          = ArrayRCP<UN>;
        using ConstUNVecView                    = ArrayView<const UN>;

        using IntVec                            = Array<int>;
        using IntVec2D                          = Array<IntVec>;

        using LOVec                             = Array<LO>;
        using LOVecPtr                          = ArrayRCP<LO>;
        using LOVecView                         = ArrayView<LO>;
        using ConstLOVecView                    = ArrayView<const LO>;
        using LOVecPtr2D                        = ArrayRCP<LOVecPtr>;

        using GOVec                             = Array<GO>;
        using GOVecPtr                          = ArrayRCP<GO>;
        using GOVecView                         = ArrayView<GO>;
        using ConstGOVecView                    = ArrayView<const GO>;
        using GOVec2D                           = Array<GOVec>;
        using GOVecPtr2D                        = ArrayRCP<GOVecPtr>;

        using SCVec                             = Array<SC>;
        using SCVecPtr                          = ArrayRCP<SC>;
        using ConstSCVecPtr                     = ArrayRCP<const SC>;
        using ConstSCVecView                    = ArrayView<const SC>;

        using BoolVec                           = Array<bool>;
        using BoolVecPtr                        = ArrayRCP<bool>;
        using ConstBoolVecPtr                   = ArrayRCP<const bool>;

    public:

        SchwarzOperator(CommPtr comm);

        SchwarzOperator(ConstXMatrixPtr k,
                        ParameterListPtr parameterList);

        virtual ~SchwarzOperator();

        virtual int initialize() = 0;

        virtual int compute() = 0;

        // Y = alpha * A^mode * X + beta * Y
        virtual void apply(const XMultiVector &x,
                           XMultiVector &y,
                           ETransp mode=NO_TRANS,
                           SC alpha=ScalarTraits<SC>::one(),
                           SC beta=ScalarTraits<SC>::zero()) const;

        virtual void apply(const XMultiVector &x,
                           XMultiVector &y,
                           bool usePreconditionerOnly,
                           ETransp mode=NO_TRANS,
                           SC alpha=ScalarTraits<SC>::one(),
                           SC beta=ScalarTraits<SC>::zero()) const = 0;

        virtual const ConstXMapPtr getDomainMap() const;

        virtual const ConstXMapPtr getRangeMap() const;

        virtual void describe(FancyOStream &out,
                              const EVerbosityLevel verbLevel=Describable::verbLevel_default) const = 0;

        virtual string description() const = 0;

        bool isInitialized() const;

        bool isComputed() const;

        int resetMatrix(ConstXMatrixPtr &k);

        virtual void residual(const XMultiVector & X,
                              const XMultiVector & B,
                              XMultiVector& R) const;

    protected:

        CommPtr MpiComm_;
        CommPtr SerialComm_ = createSerialComm<int>();

        ConstXMatrixPtr K_;

        ParameterListPtr ParameterList_;

        // result of ExtractLocalSubdomainMatrix_Symbolic
        bool ExtractLocalSubdomainMatrix_Symbolic_Done_ = false;
        XMatrixPtr subdomainMatrix_;
        XMatrixPtr localSubdomainMatrix_;
        XImportPtr subdomainScatter_;

        bool Verbose_ = false;

        bool IsInitialized_ = false;
        bool IsComputed_ = false;

        ConstUN LevelID_ = 1;
    };

}

#endif
