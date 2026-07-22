// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_ONELEVELPRECONDITIONER_DECL_HPP
#define _FROSCH_ONELEVELPRECONDITIONER_DECL_HPP

#include <FROSch_SchwarzPreconditioner_def.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;    

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class OneLevelPreconditioner : public SchwarzPreconditioner<SC,LO,GO,NO> {

    protected:

        using XMapPtr                           = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr                      = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMapPtr;

        using XMatrixPtr                        = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr                   = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMatrixPtr;

        using XMultiVector                      = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMultiVector;

        using ParameterListPtr                  = typename SchwarzPreconditioner<SC,LO,GO,NO>::ParameterListPtr;

        using SumOperatorPtr                    = typename SchwarzPreconditioner<SC,LO,GO,NO>::SumOperatorPtr;
        using MultiplicativeOperatorPtr         = typename SchwarzPreconditioner<SC,LO,GO,NO>::MultiplicativeOperatorPtr;
        using OverlappingOperatorPtr            = typename SchwarzPreconditioner<SC,LO,GO,NO>::OverlappingOperatorPtr;
        using AlgebraicOverlappingOperatorPtr   = typename SchwarzPreconditioner<SC,LO,GO,NO>::AlgebraicOverlappingOperatorPtr;

    public:

        OneLevelPreconditioner(ConstXMatrixPtr k,
                               ParameterListPtr parameterList);

        virtual int initialize(bool useDefaultParameters = true);

        virtual int initialize(int overlap = -1,
                               bool buildRepeatedMap = false);

        virtual int initialize(int overlap,
                               ConstXMapPtr repeatedMap);

        virtual int compute();

        virtual void apply(const XMultiVector &x,
                           XMultiVector &y,
                           ETransp mode=NO_TRANS,
                           SC alpha=ScalarTraits<SC>::one(),
                           SC beta=ScalarTraits<SC>::zero()) const;

        virtual const ConstXMapPtr getDomainMap() const;

        virtual const ConstXMapPtr getRangeMap() const;

        virtual void describe(FancyOStream &out,
                              const EVerbosityLevel verbLevel=Describable::verbLevel_default) const;

        virtual string description() const;

        virtual int resetMatrix(ConstXMatrixPtr &k);

    protected:

        ConstXMatrixPtr K_;

        SumOperatorPtr SumOperator_;
        MultiplicativeOperatorPtr MultiplicativeOperator_;
        OverlappingOperatorPtr OverlappingOperator_;
        bool UseMultiplicative_ = false;
    };

}

#endif
