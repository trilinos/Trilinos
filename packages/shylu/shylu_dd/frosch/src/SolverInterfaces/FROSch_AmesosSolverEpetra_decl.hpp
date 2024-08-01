// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_AMESOSSOLVEREPETRA_DECL_HPP
#define _FROSCH_AMESOSSOLVEREPETRA_DECL_HPP

#include <ShyLU_DDFROSch_config.h>

#include "Amesos_ConfigDefs.h"
#include "Amesos.h"
#include "Amesos_BaseSolver.h"

#include <Xpetra_EpetraMultiVector.hpp>

// FROSch
#include <FROSch_Solver_def.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class AmesosSolverEpetra : public Solver<SC,LO,GO,NO> {

    protected:

        // Xpetra
        using ConstXMatrixPtr                   = typename Solver<SC,LO,GO,NO>::ConstXMatrixPtr;

        using XMultiVector                      = typename Solver<SC,LO,GO,NO>::XMultiVector;
        using XMultiVectorPtr                   = typename Solver<SC,LO,GO,NO>::XMultiVectorPtr;

        using XMultiVectorFactory               = typename Solver<SC,LO,GO,NO>::XMultiVectorFactory;

        // Epetra
        using ECrsMatrix                        = Epetra_CrsMatrix;
        using ECrsMatrixPtr                     = RCP<ECrsMatrix>;
        using ConstECrsMatrixPtr                = RCP<const ECrsMatrix>;

        using EMultiVector                      = Epetra_MultiVector;
        using EMultiVectorPtr                   = RCP<EMultiVector>;

        using ELinearProblem                    = Epetra_LinearProblem;
        using ELinearProblemPtr                 = RCP<Epetra_LinearProblem>;

        // Teuchos
        using ParameterListPtr                  = typename Solver<SC,LO,GO,NO>::ParameterListPtr;

        // Amesos
        using AmesosSolverPtr                   = RCP<Amesos_BaseSolver>;

    public:

        ~AmesosSolverEpetra();

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
        AmesosSolverEpetra(ConstXMatrixPtr k,
                           ParameterListPtr parameterList,
                           string description);

        mutable XMultiVectorPtr Y_ = null;

        AmesosSolverPtr AmesosSolver_ = null;

        ELinearProblemPtr EpetraLinearProblem_ = null;

        friend class SolverFactory<SC,LO,GO,NO>;
    };

}

#endif
