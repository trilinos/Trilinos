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

#ifndef _FROSCH_SCHWARZOPERATOR_DECL_HPP
#define _FROSCH_SCHWARZOPERATOR_DECL_HPP

#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_TripleMatrixMultiply.hpp>
#include <Xpetra_Export.hpp>

#include <Teuchos_DefaultSerialComm.hpp>

#include <ShyLU_DDFROSch_config.h>

#include <FROSch_DDInterface_def.hpp>
#include <FROSch_EntitySet_def.hpp>

#include <FROSch_CoarseSpace_def.hpp>
#include <FROSch_InterfacePartitionOfUnity_def.hpp>
#include <FROSch_LocalPartitionOfUnityBasis_def.hpp>

#include <FROSch_SubdomainSolver_def.hpp>

// TODO: Auf const überprüfen
// TODO: #ifndef überprüfen ??????


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    class Solver;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = KokkosClassic::DefaultNode::DefaultNodeType>
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

        using XExport                           = Export<LO,GO,NO>;
        using XExportPtr                        = RCP<XExport>;
        using XExportPtrVecPtr                  = ArrayRCP<XExportPtr>;

        using ParameterListPtr                  = RCP<ParameterList>;

        using DDInterfacePtr                    = RCP<DDInterface<SC,LO,GO,NO> >;

        using EntitySetPtr                      = RCP<EntitySet<SC,LO,GO,NO> >;
        using EntitySetPtrVecPtr                = ArrayRCP<EntitySetPtr>;
        using EntitySetPtrConstVecPtr           = const EntitySetPtrVecPtr;

        using CoarseSpacePtr                    = RCP<CoarseSpace<SC,LO,GO,NO> >;
        using CoarseSpacePtrVecPtr              = ArrayRCP<CoarseSpacePtr>;

        using InterfaceEntityPtr                = RCP<InterfaceEntity<SC,LO,GO,NO> >;

        using InterfacePartitionOfUnityPtr      = RCP<InterfacePartitionOfUnity<SC,LO,GO,NO> >;

        using LocalPartitionOfUnityBasisPtr     = RCP<LocalPartitionOfUnityBasis<SC,LO,GO,NO> >;

        using SchwarzOperatorPtr                = RCP<SchwarzOperator<SC,LO,GO,NO> >;
        using SchwarzOperatorPtrVec             = Array<SchwarzOperatorPtr>;
        using SchwarzOperatorPtrVecPtr          = ArrayRCP<SchwarzOperatorPtr>;

        using SubdomainSolverPtr                = RCP<SubdomainSolver<SC,LO,GO,NO> >;

        using DofOrderingVecPtr                 = ArrayRCP<DofOrdering>;

        using UN                                = unsigned;
        using UNVec                             = Array<UN>;
        using UNVecPtr                          = ArrayRCP<UN>;

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

        virtual ConstXMapPtr getDomainMap() const;

        virtual ConstXMapPtr getRangeMap() const;

        virtual void describe(FancyOStream &out,
                              const EVerbosityLevel verbLevel=Describable::verbLevel_default) const = 0;

        virtual std::string description() const = 0;

        bool isInitialized() const;

        bool isComputed() const;

        int resetMatrix(ConstXMatrixPtr &k);

    protected:

        CommPtr MpiComm_;
        CommPtr SerialComm_;

        ConstXMatrixPtr K_;

        ParameterListPtr ParameterList_;

        bool Verbose_;

        bool IsInitialized_;
        bool IsComputed_;

        UN LevelID_;
    };

}

#endif
