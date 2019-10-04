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

#ifndef _FROSCH_SCHWARZPRECONDITIONER_DECL_HPP
#define _FROSCH_SCHWARZPRECONDITIONER_DECL_HPP

#include <Xpetra_Operator.hpp>
#include <Xpetra_Matrix_fwd.hpp>

#include <FROSch_SumOperator_def.hpp>
#include <FROSch_MultiplicativeOperator_def.hpp>
#include <FROSch_AlgebraicOverlappingOperator_def.hpp>
#include <FROSch_GDSWCoarseOperator_def.hpp>
#include <FROSch_RGDSWCoarseOperator_def.hpp>
#include <FROSch_IPOUHarmonicCoarseOperator_def.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = KokkosClassic::DefaultNode::DefaultNodeType>
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

        virtual ConstXMapPtr getDomainMap() const = 0;

        virtual ConstXMapPtr getRangeMap() const = 0;

        virtual void describe(FancyOStream &out,
                              const EVerbosityLevel verbLevel=Describable::verbLevel_default) const = 0;

        virtual std::string description() const = 0;

        bool isInitialized() const;

        bool isComputed() const;


    protected:

        CommPtr MpiComm_;

        ParameterListPtr ParameterList_;

        bool UseTranspose_;
        bool IsInitialized_;
        bool IsComputed_;
        bool Verbose_;

        ConstUN LevelID_;
    };

}

#endif
