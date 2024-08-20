// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_THYRAPRECONDITIONER_DECL_HPP
#define _FROSCH_THYRAPRECONDITIONER_DECL_HPP

#include <ShyLU_DDFROSch_config.h>

#include "Stratimikos_FROSch_decl.hpp"

#include <Xpetra_ThyraUtils.hpp>

// FROSch
#include <FROSch_Solver_def.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Thyra;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class ThyraPreconditioner : public Solver<SC,LO,GO,NO> {

    protected:

        // Xpetra
        using XMatrix                           = typename Solver<SC,LO,GO,NO>::XMatrix;
        using XMatrixPtr                        = typename Solver<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr                   = typename Solver<SC,LO,GO,NO>::ConstXMatrixPtr;

        using XMultiVector                      = typename Solver<SC,LO,GO,NO>::XMultiVector;
        using XMultiVectorPtr                   = typename Solver<SC,LO,GO,NO>::XMultiVectorPtr;
        using ConstXMultiVectorPtr              = typename Solver<SC,LO,GO,NO>::ConstXMultiVectorPtr;

        using XMultiVectorFactory               = typename Solver<SC,LO,GO,NO>::XMultiVectorFactory;

        // Tpetra
        using TCrsMatrix                        = Tpetra::CrsMatrix<SC,LO,GO,NO>;
        using TCrsMatrixPtr                     = RCP<TCrsMatrix>;
        using ConstTCrsMatrixPtr                = RCP<const TCrsMatrix>;

        using TMultiVector                      = Tpetra::MultiVector<SC,LO,GO,NO>;
        using TMultiVectorPtr                   = RCP<TMultiVector>;

        // Teuchos
        using ParameterListPtr                  = typename Solver<SC,LO,GO,NO>::ParameterListPtr;

        // Thyra
        using MultiVectorBasePtr                = RCP<MultiVectorBase<SC> >;
        using ConstMultiVectorBasePtr           = RCP<const MultiVectorBase<SC> >;
        using PreconditionerBasePtr             = RCP<PreconditionerBase<SC> >;
        using PreconditionerFactoryBasePtr      = RCP<PreconditionerFactoryBase<SC> >;

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
        ThyraPreconditioner(ConstXMatrixPtr k,
                            ParameterListPtr parameterList,
                            string description);

        // mutable MultiVectorBasePtr Y_ = null;

        PreconditionerBasePtr ThyraPreconditioner_ = null;

        friend class SolverFactory<SC,LO,GO,NO>;
    };

}

#endif
