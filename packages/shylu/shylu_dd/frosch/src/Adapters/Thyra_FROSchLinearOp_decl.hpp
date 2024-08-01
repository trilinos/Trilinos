// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _THYRA_FROSCH_LINEAR_OP_DECL_HPP
#define _THYRA_FROSCH_LINEAR_OP_DECL_HPP

#include <ShyLU_DDFROSch_config.h>

#ifdef HAVE_SHYLU_DDFROSCH_THYRA

//Thyra
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOpBase.hpp"
#endif
#include "Thyra_SpmdMultiVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_AssertOp.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_LinearOpDefaultBase.hpp"
#include "Thyra_ScaledLinearOpBase.hpp"
#include "Thyra_RowStatLinearOpBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"

//Teuchos
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include <Teuchos_PtrDecl.hpp>
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"

//FROSch
#include <FROSch_Tools_def.hpp>

//Xpetra
#include "Xpetra_MapExtractor.hpp"
#include <Xpetra_CrsMatrixWrap.hpp>
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
#include <Xpetra_EpetraCrsMatrix.hpp>
#endif
#include <Xpetra_Parameters.hpp>
#include "Xpetra_Operator.hpp"
#include "Xpetra_ThyraUtils.hpp"

//Epetra
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
#include <Epetra_MpiComm.h>
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_RowMatrix.h"
#endif


namespace Thyra {

    using namespace FROSch;
    using namespace Teuchos;
    using namespace Thyra;
    using namespace Xpetra;

    /** \brief Concrete Thyra::LinearOpBase subclass for Operator.**/
    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class FROSchLinearOp : virtual public LinearOpDefaultBase<SC> {
        public:

        /** \name Constructors/initializers. */
        //@{

        /** \brief Construct to uninitialized. */
        FROSchLinearOp();

        /** \brief Initialize. */
        void initialize(const RCP<const VectorSpaceBase<SC> > &rangeSpace,
                        const RCP<const VectorSpaceBase<SC> > &domainSpace,
                        const RCP<Operator<SC,LO,GO,NO> > &xpetraOperator,
                        bool bIsEpetra,
                        bool bIsTpetra);

        /** \brief Initialize. */
        void constInitialize(const RCP<const VectorSpaceBase<SC> > &rangeSpace,
                             const RCP<const VectorSpaceBase<SC> > &domainSpace,
                             const RCP<const Operator<SC,LO,GO,NO> > &xpetraOperator,
                             bool bIsEpetra,
                             bool bIsTpetra);

        /** \brief Get embedded non-const Operator. */
        RCP<Operator<SC,LO,GO,NO> > getXpetraOperator();

        /** \brief Get embedded const Operator. */
        RCP<const Operator<SC,LO,GO,NO> > getConstXpetraOperator() const;

        //@}

        /** \name Public Overridden functions from LinearOpBase. */
        //@{

        /** \brief . */
        RCP<const VectorSpaceBase<SC> > range() const;

        /** \brief . */
        RCP<const VectorSpaceBase<SC> > domain() const;

        //@}

        protected:

        /** \name Protected Overridden functions from LinearOpBase. */
        //@{

        /** \brief . */
        bool opSupportedImpl(EOpTransp M_trans) const;

        /** \brief . */
        void applyImpl(const EOpTransp M_trans,
                       const MultiVectorBase<SC> &X_in,
                       const Ptr<MultiVectorBase<SC> > &Y_inout,
                       const SC alpha,
                       const SC beta
                       ) const;

        //@}

        private:

        RCP<const VectorSpaceBase<SC> >
        rangeSpace_;

        RCP<const VectorSpaceBase<SC> >
        domainSpace_;


        bool bIsEpetra_;
        bool bIsTpetra_;
        ConstNonconstObjectContainer<Operator<SC,LO,GO,NO> >
        xpetraOperator_;

        template<class XpetraOperator_t>
        void initializeImpl(const RCP<const VectorSpaceBase<SC> > &rangeSpace,
                            const RCP<const VectorSpaceBase<SC> > &domainSpace,
                            const RCP<XpetraOperator_t> &xpetraOperator,
                            bool bIsEpetra,
                            bool bIsTpetra);
    };


    /** \brief Nonmmeber constructor for XpetraLinearOp.
     *
     * \relates XpetraLinearOp
     */
    template <class SC, class LO, class GO, class NO>
    RCP<FROSchLinearOp<SC,LO,GO,NO> > fROSchLinearOp(const RCP<const VectorSpaceBase<SC> > &rangeSpace,
                                                                                   const RCP<const VectorSpaceBase<SC> > &domainSpace,
                                                                                   const RCP<Operator<SC,LO,GO,NO> > &xpetraOperator,
                                                                                   bool bIsEpetra,
                                                                                   bool bIsTpetra)
    {
        const RCP<FROSchLinearOp<SC,LO,GO,NO> > op =
        rcp(new FROSchLinearOp<SC,LO,GO,NO>);
        op->initialize(rangeSpace,domainSpace,xpetraOperator,bIsEpetra,bIsTpetra);
        return op;
    }


    /** \brief Nonmmeber constructor for XpetraLinearOp.
     *
     * \relates XpetraLinearOp
     */
    template <class SC, class LO, class GO, class NO>
    RCP<const FROSchLinearOp<SC,LO,GO,NO> > constFROSchLinearOp(const RCP<const VectorSpaceBase<SC> > &rangeSpace,
                                                                                              const RCP<const VectorSpaceBase<SC> > &domainSpace,
                                                                                              const RCP<const Operator<SC,LO,GO,NO> > &xpetraOperator,
                                                                                              bool bIsEpetra,
                                                                                              bool bIsTpetra)
    {
        const RCP<FROSchLinearOp<SC,LO,GO,NO> > op =
        rcp(new FROSchLinearOp<SC,LO,GO,NO>);
        op->constInitialize(rangeSpace, domainSpace, xpetraOperator,bIsEpetra,bIsTpetra);
        return op;
    }

}  // namespace Thyra

#endif

#endif // THYRA_XPETRA_LINEAR_OP_DECL_HPP
