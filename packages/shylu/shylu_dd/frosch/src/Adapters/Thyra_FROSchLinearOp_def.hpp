// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _THYRA_FROSCH_LINEAR_OP_DEF_HPP
#define _THYRA_FROSCH_LINEAR_OP_DEF_HPP

#include "Thyra_FROSchLinearOp_decl.hpp"


#ifdef HAVE_SHYLU_DDFROSCH_THYRA
namespace Thyra {

    using namespace FROSch;
    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    // Constructors/initializers
    template <class SC, class LO, class GO, class NO>
    FROSchLinearOp<SC,LO,GO,NO>::FROSchLinearOp()
    {

    }

    template <class SC, class LO, class GO, class NO>
    void FROSchLinearOp<SC,LO,GO,NO>::initialize(const RCP<const VectorSpaceBase<SC> > &rangeSpace,
                                                 const RCP<const VectorSpaceBase<SC> > &domainSpace,
                                                 const RCP<Operator<SC,LO,GO,NO> > &xpetraOperator,
                                                 bool bIsTpetra)
    {
        initializeImpl(rangeSpace, domainSpace, xpetraOperator,bIsTpetra);
    }

    template <class SC, class LO, class GO, class NO>
    void FROSchLinearOp<SC,LO,GO,NO>::constInitialize(const RCP<const VectorSpaceBase<SC> > &rangeSpace,
                                                      const RCP<const VectorSpaceBase<SC> > &domainSpace,
                                                      const RCP<const Operator<SC,LO,GO,NO> > &xpetraOperator,
                                                      bool bIsTpetra)
    {
        initializeImpl(rangeSpace, domainSpace, xpetraOperator,bIsTpetra);
    }

    template <class SC, class LO, class GO, class NO>
    RCP<Operator<SC,LO,GO,NO> > FROSchLinearOp<SC,LO,GO,NO>::getXpetraOperator()
    {
        return xpetraOperator_.getNonconstObj();
    }

    template <class SC, class LO, class GO, class NO>
    RCP<const Operator<SC,LO,GO,NO> > FROSchLinearOp<SC,LO,GO,NO>::getConstXpetraOperator() const
    {
        return xpetraOperator_;
    }

    // Public Overridden functions from LinearOpBase

    template <class SC, class LO, class GO, class NO>
    RCP<const VectorSpaceBase<SC> > FROSchLinearOp<SC,LO,GO,NO>::range() const
    {
        return rangeSpace_;
    }

    template <class SC, class LO, class GO, class NO>
    RCP<const VectorSpaceBase<SC> > FROSchLinearOp<SC,LO,GO,NO>::domain() const
    {
        return domainSpace_;
    }

    // Protected Overridden functions from LinearOpBase

    template <class SC, class LO, class GO, class NO>
    bool FROSchLinearOp<SC,LO,GO,NO>::opSupportedImpl(EOpTransp M_trans) const
    {
        if (is_null(xpetraOperator_))
        return false;

        if (M_trans == NOTRANS)
        return true;

        if (M_trans == CONJ) {
            // For non-complex scalars, CONJ is always supported since it is equivalent to NO_TRANS.
            // For complex scalars, Xpetra does not support conjugation without transposition.
            return !ScalarTraits<SC>::isComplex;
        }

        return xpetraOperator_->hasTransposeApply();
    }

    template <class SC, class LO, class GO, class NO>
    void FROSchLinearOp<SC,LO,GO,NO>::applyImpl(const EOpTransp M_trans,
                                                const MultiVectorBase<SC> &X_in,
                                                const Ptr<MultiVectorBase<SC> > &Y_inout,
                                                const SC alpha,
                                                const SC beta) const
    {
        FROSCH_ASSERT(getConstXpetraOperator()!=null,"XpetraLinearOp::applyImpl: internal Operator is null.");
        RCP< const Comm<int> > comm = getConstXpetraOperator()->getRangeMap()->getComm();
        //Transform to Xpetra MultiVector
        RCP<MultiVector<SC,LO,GO,NO> > xY;

        ETransp transp = NO_TRANS;
        switch (M_trans) {
            case NOTRANS:   transp = NO_TRANS;          break;
            case TRANS:     transp = Teuchos::TRANS;    break;
            case CONJTRANS: transp = CONJ_TRANS;        break;
            default: FROSCH_ASSERT(false,"Thyra::XpetraLinearOp::apply. Unknown value for M_trans. Only NOTRANS, TRANS and CONJTRANS are supported.");
        }
        if (bIsTpetra_) {
            // Convert input vector to Xpetra
            const RCP<const Tpetra::MultiVector<SC,LO,GO,NO> > xTpMultVec = Thyra::TpetraOperatorVectorExtraction<SC,LO,GO,NO>::getConstTpetraMultiVector(rcpFromRef(X_in));
            TEUCHOS_TEST_FOR_EXCEPT(is_null(xTpMultVec));
            RCP<Tpetra::MultiVector<SC,LO,GO,NO> > tpNonConstMultVec = rcp_const_cast<Tpetra::MultiVector<SC,LO,GO,NO> >(xTpMultVec);
            TEUCHOS_TEST_FOR_EXCEPT(is_null(tpNonConstMultVec));
            const RCP<const Xpetra::MultiVector<SC,LO,GO,NO> > xX = rcp(new Xpetra::TpetraMultiVector<SC,LO,GO,NO>(tpNonConstMultVec));
            TEUCHOS_TEST_FOR_EXCEPT(is_null(xX));

            // Convert output vector to Xpetra
            const RCP<Tpetra::MultiVector<SC,LO,GO,NO> > yTpMultVec = Thyra::TpetraOperatorVectorExtraction<SC,LO,GO,NO>::getTpetraMultiVector(rcpFromPtr(Y_inout));
            TEUCHOS_TEST_FOR_EXCEPT(is_null(yTpMultVec));
            xY = rcp(new Xpetra::TpetraMultiVector<SC,LO,GO,NO>(yTpMultVec));
            TEUCHOS_TEST_FOR_EXCEPT(is_null(xY));

            // Apply operator
            xpetraOperator_->apply(*xX, *xY, transp, alpha, beta);

        } else {
            FROSCH_ASSERT(false,"There is a problem with the underlying lib in FROSchLinearOp.");
            //cout<<"Only Implemented for Epetra and Tpetra\n";
        }
    }

    // private

    template <class SC, class LO, class GO, class NO>
    template<class XpetraOperator_t>
    void FROSchLinearOp<SC,LO,GO,NO>::initializeImpl(const RCP<const VectorSpaceBase<SC> > &rangeSpace,
                                                     const RCP<const VectorSpaceBase<SC> > &domainSpace,
                                                     const RCP<XpetraOperator_t> &xpetraOperator,
                                                     bool bIsTpetra)
    {
#ifdef THYRA_DEBUG
        TEUCHOS_ASSERT(nonnull(rangeSpace));
        TEUCHOS_ASSERT(nonnull(domainSpace));
        TEUCHOS_ASSERT(nonnull(xpetraOperator));
#endif
        rangeSpace_ = rangeSpace;
        domainSpace_ = domainSpace;
        xpetraOperator_ = xpetraOperator;
        bIsTpetra_ = bIsTpetra;
    }

} // namespace Thyra

#endif

#endif  // THYRA_XPETRA_LINEAR_OP_HPP
