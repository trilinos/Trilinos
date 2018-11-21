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

#ifndef THYRA_FROSCH_LINEAR_OP_DECL_HPP
#define THYRA_FROSCH_LINEAR_OP_DECL_HPP

#include "Thyra_LinearOpDefaultBase.hpp"
#include "Xpetra_Operator.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"


namespace Thyra {
    
    using namespace Thyra;
    using namespace Xpetra;
    
    /** \brief Concrete Thyra::LinearOpBase subclass for Xpetra::Operator.**/
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal=LocalOrdinal,
    class Node=KokkosClassic::DefaultNode::DefaultNodeType>
    class FROSchLinearOp : virtual public LinearOpDefaultBase<Scalar>
    {
        public:
        
        /** \name Constructors/initializers. */
        //@{
        
        /** \brief Construct to uninitialized. */
        FROSchLinearOp();
        
        /** \brief Initialize. */
        void initialize(const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
                        const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
                        const RCP<Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &xpetraOperator,
                        bool bIsEpetra,
                        bool bIsTpetra);
        
        /** \brief Initialize. */
        void constInitialize(const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
                             const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
                             const RCP<const Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &xpetraOperator,
                             bool bIsEpetra,
                             bool bIsTpetra);
        
        /** \brief Get embedded non-const Xpetra::Operator. */
        RCP<Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
        getXpetraOperator();
        
        /** \brief Get embedded const Xpetra::Operator. */
        RCP<const Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
        getConstXpetraOperator() const;
        
        //@}
        
        /** \name Public Overridden functions from LinearOpBase. */
        //@{
        
        /** \brief . */
        RCP<const VectorSpaceBase<Scalar> > range() const;
        
        /** \brief . */
        RCP<const VectorSpaceBase<Scalar> > domain() const;
        
        //@}
        
        protected:
        
        /** \name Protected Overridden functions from LinearOpBase. */
        //@{
        
        /** \brief . */
        bool opSupportedImpl(EOpTransp M_trans) const;
        
        /** \brief . */
        void applyImpl(
                       const EOpTransp M_trans,
                       const MultiVectorBase<Scalar> &X_in,
                       const Teuchos::Ptr<MultiVectorBase<Scalar> > &Y_inout,
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
        Teuchos::ConstNonconstObjectContainer<Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
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
                        const RCP<const Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &xpetraOperator,
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

