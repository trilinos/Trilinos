// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_MULTIPLICATIVEOPERATOR_DECL_HPP
#define _FROSCH_MULTIPLICATIVEOPERATOR_DECL_HPP

#include <FROSch_SchwarzOperator_def.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class MultiplicativeOperator : public SchwarzOperator<SC,LO,GO,NO> {

    protected:

        using CommPtr                   = typename SchwarzOperator<SC,LO,GO,NO>::CommPtr;

        using XMapPtr                   = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr              = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtr;

        using XMultiVector              = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVector;
        using XMultiVectorPtr           = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVectorPtr;

        using XMatrixPtr                = typename SchwarzOperator<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr           = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMatrixPtr;

        using SchwarzOperatorPtr        = typename SchwarzOperator<SC,LO,GO,NO>::SchwarzOperatorPtr;
        using SchwarzOperatorPtrVec     = typename SchwarzOperator<SC,LO,GO,NO>::SchwarzOperatorPtrVec;
        using SchwarzOperatorPtrVecPtr  = typename SchwarzOperator<SC,LO,GO,NO>::SchwarzOperatorPtrVecPtr;

        using UN                        = typename SchwarzOperator<SC,LO,GO,NO>::UN;

        using BoolVec                   = typename SchwarzOperator<SC,LO,GO,NO>::BoolVec;

        using ParameterListPtr          = typename SchwarzOperator<SC,LO,GO,NO>::ParameterListPtr;

    public:

        MultiplicativeOperator(ConstXMatrixPtr k,
                               ParameterListPtr parameterList);

        MultiplicativeOperator(ConstXMatrixPtr k,
                               SchwarzOperatorPtrVecPtr operators,
                               ParameterListPtr parameterList);

        ~MultiplicativeOperator();

        virtual int initialize();

        virtual int initialize(ConstXMapPtr repeatedMap);

        virtual int compute();

        void preApplyCoarse(XMultiVector &x,
                            XMultiVector &y);

        virtual void apply(const XMultiVector &x,
                           XMultiVector &y,
                           bool usePreconditionerOnly,
                           ETransp mode=NO_TRANS,
                           SC alpha=ScalarTraits<SC>::one(),
                           SC beta=ScalarTraits<SC>::zero()) const;

        virtual const ConstXMapPtr getDomainMap() const;

        virtual const ConstXMapPtr getRangeMap() const;

        virtual void describe(FancyOStream &out,
                              const EVerbosityLevel verbLevel=Describable::verbLevel_default) const;

        virtual string description() const;

        int addOperator(SchwarzOperatorPtr op);

        int addOperators(SchwarzOperatorPtrVecPtr operators);

        int resetOperator(UN iD,
                          SchwarzOperatorPtr op);

        int enableOperator(UN iD,
                           bool enable);

        UN getNumOperators();

    protected:

        SchwarzOperatorPtrVec OperatorVector_ = SchwarzOperatorPtrVec(0);

        // Temp Vectors for apply()
        mutable XMultiVectorPtr XTmp_;
        mutable XMultiVectorPtr YTmp_;

        BoolVec EnableOperators_ = BoolVec(0);
    };

}

#endif
