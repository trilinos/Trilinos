#ifndef THYRA_FROSCH_LINEAR_OP_DECL_HPP
#define THYRA_FROSCH_LINEAR_OP_DECL_HPP

#include "Thyra_LinearOpDefaultBase.hpp"
#include "Xpetra_Operator.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"


namespace Thyra {
    
    
    /** \brief Concrete Thyra::LinearOpBase subclass for Xpetra::Operator.**/
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal=LocalOrdinal,
    class Node=KokkosClassic::DefaultNode::DefaultNodeType>
    class FROSchLinearOp
    : virtual public Thyra::LinearOpDefaultBase<Scalar>
    {
        public:
        
        /** \name Constructors/initializers. */
        //@{
        
        /** \brief Construct to uninitialized. */
        FROSchLinearOp();
        
        /** \brief Initialize. */
        void initialize(
                        const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
                        const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
                        const RCP<Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &xpetraOperator,
                        bool bIsEpetra,
                        bool bIsTpetra
                        );
        
        /** \brief Initialize. */
        void constInitialize(
                             const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
                             const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
                             const RCP<const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &xpetraOperator,
                             bool bIsEpetra,
                             bool bIsTpetra
                             );
        
        /** \brief Get embedded non-const Xpetra::Operator. */
        RCP<Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
        getXpetraOperator();
        
        /** \brief Get embedded const Xpetra::Operator. */
        RCP<const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
        getConstXpetraOperator() const;
        
        //@}
        
        /** \name Public Overridden functions from LinearOpBase. */
        //@{
        
        /** \brief . */
        RCP<const Thyra::VectorSpaceBase<Scalar> > range() const;
        
        /** \brief . */
        RCP<const Thyra::VectorSpaceBase<Scalar> > domain() const;
        
        //@}
        
        protected:
        
        /** \name Protected Overridden functions from LinearOpBase. */
        //@{
        
        /** \brief . */
        bool opSupportedImpl(Thyra::EOpTransp M_trans) const;
        
        /** \brief . */
        void applyImpl(
                       const Thyra::EOpTransp M_trans,
                       const Thyra::MultiVectorBase<Scalar> &X_in,
                       const Teuchos::Ptr<Thyra::MultiVectorBase<Scalar> > &Y_inout,
                       const Scalar alpha,
                       const Scalar beta
                       ) const;
        
        //@}
        
        private:
        
        RCP<const VectorSpaceBase<Scalar> >
        rangeSpace_;
        
        RCP<const VectorSpaceBase<Scalar> >
        domainSpace_;
        
        
        bool bIsEpetra_;
        bool bIsTpetra_;
        Teuchos::ConstNonconstObjectContainer<Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
        xpetraOperator_;
        
        template<class XpetraOperator_t>
        void initializeImpl(
                            const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
                            const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
                            const RCP<XpetraOperator_t> &xpetraOperator,
                            bool bIsEpetra,
                            bool bIsTpetra
                            );
        
    };
    
    
    /** \brief Nonmmeber constructor for XpetraLinearOp.
     *
     * \relates XpetraLinearOp
     */
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    RCP<FROSchLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    fROSchLinearOp(
                   const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
                   const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
                   const RCP<Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &xpetraOperator,
                   bool bIsEpetra,
                   bool bIsTpetra
                   )
    {
        const RCP<FROSchLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> > op =
        Teuchos::rcp(new FROSchLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>);
        op->initialize(rangeSpace, domainSpace, xpetraOperator,bIsEpetra,bIsTpetra);
        return op;
    }
    
    
    /** \brief Nonmmeber constructor for XpetraLinearOp.
     *
     * \relates XpetraLinearOp
     */
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    RCP<const FROSchLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    constFROSchLinearOp(
                        const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
                        const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
                        const RCP<const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &xpetraOperator,
                        bool bIsEpetra,
                        bool bIsTpetra
                        )
    {
        const RCP<FROSchLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> > op =
        Teuchos::rcp(new FROSchLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>);
        op->constInitialize(rangeSpace, domainSpace, xpetraOperator,bIsEpetra,bIsTpetra);
        return op;
    }
    
}  // namespace Thyra

#endif // THYRA_XPETRA_LINEAR_OP_DECL_HPP

