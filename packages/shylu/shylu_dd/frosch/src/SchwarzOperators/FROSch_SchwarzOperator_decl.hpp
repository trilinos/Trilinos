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

#define FROSCH_ASSERT(A,S) if(!(A)) { std::cerr<<"Assertion failed. "<<S<<std::endl; std::cout.flush(); throw std::out_of_range("Assertion.");};

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

    class Solver;
    
    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = KokkosClassic::DefaultNode::DefaultNodeType>
    class SchwarzOperator : public Xpetra::Operator<SC,LO,GO,NO> {

    protected:

        using CommPtr                           = Teuchos::RCP<const Teuchos::Comm<int> >;

        using Map                               = Xpetra::Map<LO,GO,NO>;
        using MapPtr                            = Teuchos::RCP<Map>;
        using ConstMapPtr                       = Teuchos::RCP<const Map>;
        using MapPtrVecPtr                      = Teuchos::ArrayRCP<MapPtr>;
        using ConstMapPtrVecPtr                 = Teuchos::ArrayRCP<ConstMapPtr>;
        using MapPtrVecPtr2D                    = Teuchos::ArrayRCP<MapPtrVecPtr>;
        using ConstMapPtrVecPtr2D               = Teuchos::ArrayRCP<ConstMapPtrVecPtr>;

        using CrsMatrix                         = Xpetra::Matrix<SC,LO,GO,NO>;
        using CrsMatrixPtr                      = Teuchos::RCP<CrsMatrix>;
        using ConstCrsMatrixPtr                 = Teuchos::RCP<const CrsMatrix>;

        using Graph                             = Xpetra::CrsGraph<LO,GO,NO>;
        using GraphPtr                          = Teuchos::RCP<Graph>;
        using ConstGraphPtr                     = Teuchos::RCP<const Graph>;

        using MultiVector                       = Xpetra::MultiVector<SC,LO,GO,NO>;
        using MultiVectorPtr                    = Teuchos::RCP<MultiVector>;
        using ConstMultiVectorPtr               = Teuchos::RCP<const MultiVector>;
        using MultiVectorPtrVecPtr              = Teuchos::ArrayRCP<MultiVectorPtr>;
        using ConstMultiVectorPtrVecPtr         = Teuchos::ArrayRCP<ConstMultiVectorPtr>;

        using Importer                          = Xpetra::Import<LO,GO,NO>;
        using ImporterPtr                       = Teuchos::RCP<Importer>;

        using Exporter                          = Xpetra::Export<LO,GO,NO>;
        using ExporterPtr                       = Teuchos::RCP<Exporter>;
        using ExporterPtrVecPtr                 = Teuchos::ArrayRCP<ExporterPtr>;

        using ParameterList                     = Teuchos::ParameterList;
        using ParameterListPtr                  = Teuchos::RCP<Teuchos::ParameterList>;

        using DDInterfacePtr                    = Teuchos::RCP<DDInterface<SC,LO,GO,NO> >;

        using EntitySetPtr                      = Teuchos::RCP<EntitySet<SC,LO,GO,NO> >;
        using EntitySetPtrVecPtr                = Teuchos::ArrayRCP<EntitySetPtr>;
        using EntitySetPtrConstVecPtr           = const EntitySetPtrVecPtr;

        using CoarseSpacePtr                    = Teuchos::RCP<CoarseSpace<SC,LO,GO,NO> >;
        using CoarseSpacePtrVecPtr              = Teuchos::ArrayRCP<CoarseSpacePtr>;

        using InterfaceEntityPtr                = Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> >;

        using InterfacePartitionOfUnityPtr      = Teuchos::RCP<InterfacePartitionOfUnity<SC,LO,GO,NO> >;

        using LocalPartitionOfUnityBasisPtr     = Teuchos::RCP<LocalPartitionOfUnityBasis<SC,LO,GO,NO> >;

        using SchwarzOperatorPtr                = Teuchos::RCP<SchwarzOperator<SC,LO,GO,NO> >;
        using SchwarzOperatorPtrVec             = Teuchos::Array<SchwarzOperatorPtr>;
        using SchwarzOperatorPtrVecPtr          = Teuchos::ArrayRCP<SchwarzOperatorPtr>;

        using SubdomainSolverPtr                = Teuchos::RCP<SubdomainSolver<SC,LO,GO,NO> >;

        using UN                                = unsigned;
        using UNVec                             = Teuchos::Array<UN>;
        using UNVecPtr                          = Teuchos::ArrayRCP<UN>;

        using LOVec                             = Teuchos::Array<LO>;
        using LOVecPtr                          = Teuchos::ArrayRCP<LO>;
        using LOVecView                         = Teuchos::ArrayView<LO>;
        using ConstLOVecView                    = Teuchos::ArrayView<const LO>;
        using LOVecPtr2D                        = Teuchos::ArrayRCP<LOVecPtr>;

        using GOVec                             = Teuchos::Array<GO>;
        using GOVecPtr                          = Teuchos::ArrayRCP<GO>;
        using GOVecView                         = Teuchos::ArrayView<GO>;
        using ConstGOVecView                    = Teuchos::ArrayView<const GO>;
        using GOVec2D                           = Teuchos::Array<GOVec>;
        using GOVecPtr2D                        = Teuchos::ArrayRCP<GOVecPtr>;

        using SCVec                             = Teuchos::Array<SC>;
        using SCVecPtr                          = Teuchos::ArrayRCP<SC>;
        using ConstSCVecPtr                     = Teuchos::ArrayRCP<const SC>;
        using ConstSCVecView                    = Teuchos::ArrayView<const SC>;

        using BoolVec                           = Teuchos::Array<bool>;
        using BoolVecPtr                        = Teuchos::ArrayRCP<bool>;

    public:

        SchwarzOperator(CommPtr comm);

        SchwarzOperator(ConstCrsMatrixPtr k,
                        ParameterListPtr parameterList);

        virtual ~SchwarzOperator();

        virtual int initialize() = 0;

        virtual int compute() = 0;

        // Y = alpha * A^mode * X + beta * Y
        virtual void apply(const MultiVector &x,
                          MultiVector &y,
                          Teuchos::ETransp mode=Teuchos::NO_TRANS,
                          SC alpha=Teuchos::ScalarTraits<SC>::one(),
                          SC beta=Teuchos::ScalarTraits<SC>::zero()) const;

        virtual void apply(const MultiVector &x,
                          MultiVector &y,
                          bool usePreconditionerOnly,
                          Teuchos::ETransp mode=Teuchos::NO_TRANS,
                          SC alpha=Teuchos::ScalarTraits<SC>::one(),
                          SC beta=Teuchos::ScalarTraits<SC>::zero()) const = 0;

        virtual ConstMapPtr getDomainMap() const;

        virtual ConstMapPtr getRangeMap() const;

        virtual void describe(Teuchos::FancyOStream &out,
                              const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const = 0;

        virtual std::string description() const = 0;

        bool isInitialized() const;

        bool isComputed() const;

        int resetMatrix(ConstCrsMatrixPtr &k);

    protected:

        CommPtr MpiComm_;
        CommPtr SerialComm_;

        ConstCrsMatrixPtr K_;

        ParameterListPtr ParameterList_;

        bool Verbose_;

        bool IsInitialized_;
        bool IsComputed_;

    };

}

#endif
