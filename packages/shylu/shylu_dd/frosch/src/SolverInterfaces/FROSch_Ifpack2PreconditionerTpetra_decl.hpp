// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_IFPACK2PRECONDITIONERTPETRA_DECL_HPP
#define _FROSCH_IFPACK2PRECONDITIONERTPETRA_DECL_HPP

#include <ShyLU_DDFROSch_config.h>

#include <Ifpack2_ReorderFilter_decl.hpp>
#include <Ifpack2_Details_OneLevelFactory_decl.hpp>

// FROSch
#include <FROSch_Solver_def.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class Ifpack2PreconditionerTpetra : public Solver<SC,LO,GO,NO> {

    protected:

        // Xpetra
        using ConstXMatrixPtr                   = typename Solver<SC,LO,GO,NO>::ConstXMatrixPtr;

        using XMultiVector                      = typename Solver<SC,LO,GO,NO>::XMultiVector;
        using XMultiVectorPtr                   = typename Solver<SC,LO,GO,NO>::XMultiVectorPtr;

        using XMultiVectorFactory               = typename Solver<SC,LO,GO,NO>::XMultiVectorFactory;

        // Tpetra
        using TOperator                         = Tpetra::Operator<SC,LO,GO,NO>;

        using TRowMatrix                        = Tpetra::RowMatrix<SC,LO,GO,NO>;
        using TRowMatrixPtr                     = RCP<TRowMatrix>;
        using ConstTRowMatrixPtr                = RCP<const TRowMatrix>;

        using TCrsMatrix                        = Tpetra::CrsMatrix<SC,LO,GO,NO>;
        using TCrsMatrixPtr                     = RCP<TCrsMatrix>;
        using ConstTCrsMatrixPtr                = RCP<const TCrsMatrix>;

        using TMultiVector                      = Tpetra::MultiVector<SC,LO,GO,NO>;
        using TMultiVectorPtr                   = RCP<TMultiVector>;

        // Teuchos
        using ParameterListPtr                  = typename Solver<SC,LO,GO,NO>::ParameterListPtr;

        // Ifpack2
        using Ifpack2PreconditionerPtr          = RCP<Ifpack2::Preconditioner<SC,LO,GO,NO> >;
        using TRowMatrixFilterType              = Ifpack2::ReorderFilter<TRowMatrix>;

    public:

        //! Initialize the internal solver
        virtual int initialize();

        //! Compute the internal solver
        virtual int compute();

        /*!
        \brief Computes the operator-multivector application.

        Loosely, performs \f$Y = \alpha \cdot A^{\textrm{mode}} \cdot X + \beta \cdot Y\f$. However, the details of operation
        vary according to the values of \c alpha and \c beta. Specifically
        - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y (including NaNs) are ignored.
        - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the operator, so that any values in \c X (including NaNs) are ignored.
        */
        void apply(const XMultiVector &x,
                   XMultiVector &y,
                   ETransp mode=NO_TRANS,
                   SC alpha=ScalarTraits<SC>::one(),
                   SC beta=ScalarTraits<SC>::zero()) const;

        int updateMatrix(ConstXMatrixPtr k,
                        bool reuseInitialize=false);

    protected:

        //! Constructor
        Ifpack2PreconditionerTpetra(ConstXMatrixPtr k,
                                    ParameterListPtr parameterList,
                                    string description);

        Ifpack2PreconditionerPtr Ifpack2Preconditioner_ = null;

        bool useRILUK;
        bool useZoltan2;
        bool needToApplyPerm;
        Teuchos::ArrayRCP<LO>    perm;
        Teuchos::ArrayRCP<LO>    revperm;

        friend class SolverFactory<SC,LO,GO,NO>;
    };

}

#endif
