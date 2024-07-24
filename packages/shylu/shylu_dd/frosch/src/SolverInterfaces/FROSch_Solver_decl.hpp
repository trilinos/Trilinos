// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_SOLVER_DECL_HPP
#define _FROSCH_SOLVER_DECL_HPP

#include <ShyLU_DDFROSch_config.h>

//Teuchos
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

// FROSch
#include <FROSch_Tools_def.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class SolverFactory;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class Solver : public Operator<SC,LO,GO,NO> {

    protected:

        // Xpetra
        using XMap                              = Map<LO,GO,NO>;
        using XMapPtr                           = RCP<XMap>;
        using ConstXMapPtr                      = RCP<const XMap>;

        using XMatrix                           = Matrix<SC,LO,GO,NO>;
        using XMatrixPtr                        = RCP<XMatrix>;
        using ConstXMatrixPtr                   = RCP<const XMatrix>;

        using XMultiVector                      = MultiVector<SC,LO,GO,NO>;
        using XMultiVectorPtr                   = RCP<XMultiVector>;
        using ConstXMultiVectorPtr              = RCP<const XMultiVector>;

        using XMultiVectorFactory               = MultiVectorFactory<SC,LO,GO,NO>;

        // Teuchos
        using ParameterListPtr                  = RCP<ParameterList>;

    public:

        //! Initialize the internal solver
        virtual int initialize() = 0;

        //! Compute the internal solver
        virtual int compute() = 0;

        /*!
        \brief Computes the operator-multivector application.

        Loosely, performs \f$Y = \alpha \cdot A^{\textrm{mode}} \cdot X + \beta \cdot Y\f$. However, the details of operation
        vary according to the values of \c alpha and \c beta. Specifically
        - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y (including NaNs) are ignored.
        - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the operator, so that any values in \c X (including NaNs) are ignored.
        */
        virtual void apply(const XMultiVector &x,
                           XMultiVector &y,
                           ETransp mode=NO_TRANS,
                           SC alpha=ScalarTraits<SC>::one(),
                           SC beta=ScalarTraits<SC>::zero()) const = 0;

        //! Get domain map
        virtual const ConstXMapPtr getDomainMap() const;

        //! Get range map
        virtual const ConstXMapPtr getRangeMap() const;

        //! Get #IsInitialized_
        bool isInitialized() const;

        //! Get #IsComputed_
        bool isComputed() const;

        virtual int updateMatrix(ConstXMatrixPtr k,
                                 bool reuseInitialize=false) = 0;

        /*!
        \brief Computes a residual using the operator
        */
        virtual void residual(const XMultiVector & X,
                              const XMultiVector & B,
                              XMultiVector& R) const;

    protected:

        //! Constructor
        Solver(ConstXMatrixPtr k,
               ParameterListPtr parameterList,
               string description);

        //! Matrix
        ConstXMatrixPtr K_;

        //! Parameter list
        ParameterListPtr ParameterList_;

        //! Description of the solver
        string Description_;

        bool IsInitialized_ = false;

        bool IsComputed_ = false;

        friend class SolverFactory<SC,LO,GO,NO>;
    };

}

#endif
