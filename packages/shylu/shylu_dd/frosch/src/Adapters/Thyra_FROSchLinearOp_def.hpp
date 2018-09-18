#ifndef THYRA_FROSCH_LINEAR_OP_HPP
#define THYRA_FROSCH_LINEAR_OP_HPP

//Thyra

#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_SpmdMultiVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_AssertOp.hpp"
#include "Thyra_FROSchLinearOp_decl.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_EpetraLinearOpBase.hpp"
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

//FROSch
#include <FROSch_Tools_def.hpp>

//Xpetra
#include "Xpetra_MapExtractor.hpp"
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_Parameters.hpp>

//Epetra
#include <Epetra_MpiComm.h>
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_RowMatrix.h"



using namespace std;
using namespace Teuchos;
using namespace Xpetra;
using namespace FROSch;
using namespace Belos;
namespace Thyra {
    
    
    // Constructors/initializers
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    FROSchLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::FROSchLinearOp()
    {}
    
    
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void FROSchLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::initialize(
                                                                            const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
                                                                            const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
                                                                            const RCP<Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &xpetraOperator,
                                                                            bool bIsEpetra,
                                                                            bool bIsTpetra
                                                                            )
    {
        initializeImpl(rangeSpace, domainSpace, xpetraOperator,bIsEpetra,bIsTpetra);
    }
    
    
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void FROSchLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::constInitialize(
                                                                                 const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
                                                                                 const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
                                                                                 const RCP<const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &xpetraOperator,
                                                                                 bool bIsEpetra,
                                                                                 bool bIsTpetra
                                                                                 )
    {
        initializeImpl(rangeSpace, domainSpace, xpetraOperator,bIsEpetra,bIsTpetra);
    }
    
    
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    RCP<Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
    FROSchLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getXpetraOperator()
    {
        return xpetraOperator_.getNonconstObj();
    }
    
    
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    RCP<const Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
    FROSchLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getConstXpetraOperator() const
    {
        return xpetraOperator_;
    }
    
    
    // Public Overridden functions from LinearOpBase
    
    
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    RCP<const Thyra::VectorSpaceBase<Scalar> >
    FROSchLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::range() const
    {
        return rangeSpace_;
    }
    
    
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    RCP<const Thyra::VectorSpaceBase<Scalar> >
    FROSchLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::domain() const
    {
        return domainSpace_;
    }
    
    // Protected Overridden functions from LinearOpBase
    
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    bool FROSchLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::opSupportedImpl(
                                                                                 Thyra::EOpTransp M_trans) const
    {
        if (is_null(xpetraOperator_))
        return false;
        
        if (M_trans == NOTRANS)
        return true;
        
        if (M_trans == CONJ) {
            // For non-complex scalars, CONJ is always supported since it is equivalent to NO_TRANS.
            // For complex scalars, Xpetra does not support conjugation without transposition.
            return !Teuchos::ScalarTraits<Scalar>::isComplex;
        }
        
        return xpetraOperator_->hasTransposeApply();
    }
    
    
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void FROSchLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyImpl(
                                                                           const Thyra::EOpTransp M_trans,
                                                                           const Thyra::MultiVectorBase<Scalar> &X_in,
                                                                           const Teuchos::Ptr<Thyra::MultiVectorBase<Scalar> > &Y_inout,
                                                                           const Scalar alpha,
                                                                           const Scalar beta
                                                                           ) const
    {
        using Teuchos::rcpFromRef;
        using Teuchos::rcpFromPtr;
        
        const EOpTransp real_M_trans = real_trans(M_trans);

        TEUCHOS_TEST_FOR_EXCEPTION(getConstXpetraOperator() == Teuchos::null, MueLu::Exceptions::RuntimeError, "XpetraLinearOp::applyImpl: internal Xpetra::Operator is null.");
        RCP< const Teuchos::Comm<int> > comm = getConstXpetraOperator()->getRangeMap()->getComm();
        //Transform to Xpetra MultiVector
        RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > xY;
        
        Teuchos::ETransp transp;
        switch (M_trans) {
            case NOTRANS:   transp = Teuchos::NO_TRANS;   break;
            case TRANS:     transp = Teuchos::TRANS;      break;
            case CONJTRANS: transp = Teuchos::CONJ_TRANS; break;
            default: TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::NotImplemented, "Thyra::XpetraLinearOp::apply. Unknown value for M_trans. Only NOTRANS, TRANS and CONJTRANS are supported.");
        }
        //Epetra NodeType
        if(this->bIsEpetra_){
            const RCP<const VectorSpaceBase<double> > XY_domain = X_in.domain();
            
            Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > DomainM = this->xpetraOperator_->getDomainMap();
        
            Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >RangeM = this->xpetraOperator_->getRangeMap();
        
            RCP<const Xpetra::EpetraMapT<GlobalOrdinal,Node> > eDomainM = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraMapT<GlobalOrdinal,Node> >(DomainM);
        
            const Epetra_Map epetraDomain = eDomainM->getEpetra_Map();
        
            RCP<const Xpetra::EpetraMapT<GlobalOrdinal,Node> > eRangeM = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraMapT<GlobalOrdinal,Node> >(RangeM);
        
            const Epetra_Map epetraRange = eRangeM->getEpetra_Map();
            
            RCP<const Epetra_MultiVector> X;
        
            RCP<Epetra_MultiVector> Y;
       
            THYRA_FUNC_TIME_MONITOR_DIFF(
                                         "Thyra::EpetraLinearOp::euclideanApply: Convert MultiVectors", MultiVectors);
            // X
            X = Thyra::get_Epetra_MultiVector(real_M_trans==NOTRANS ? epetraDomain: epetraRange, X_in );
            RCP<Epetra_MultiVector> X_nonconst = Teuchos::rcp_const_cast<Epetra_MultiVector>(X);
            RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > xX = FROSch::ConvertToXpetra<Scalar,LocalOrdinal,GlobalOrdinal,Node>(UseEpetra,*X_nonconst,comm);
            // Y
            Y = Thyra::get_Epetra_MultiVector(real_M_trans==NOTRANS ? epetraRange: epetraDomain, *Y_inout );
            xY = FROSch::ConvertToXpetra<Scalar,LocalOrdinal,GlobalOrdinal,Node>(UseEpetra,*Y,comm);
            xpetraOperator_->apply(*xX, *xY, transp, alpha, beta);

        } //Tpetra NodeType
        else if(bIsTpetra_){
            const RCP<const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > xX =
            Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toXpetra(rcpFromRef(X_in), comm);
            xY = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toXpetra(rcpFromPtr(Y_inout), comm);
            xpetraOperator_->apply(*xX, *xY, transp, alpha, beta);
            
        }
        else{
            std::cout<<"Only Implemented for Epetra and Tpetra\n";
        }
 
       
        
        RCP<Thyra::MultiVectorBase<Scalar> >thyraX =
        Teuchos::rcp_const_cast<Thyra::MultiVectorBase<Scalar> >(Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyraMultiVector(xY));
        
        
        typedef Thyra::SpmdVectorSpaceBase<Scalar> ThySpmdVecSpaceBase;
        RCP<const ThySpmdVecSpaceBase> mpi_vs = rcp_dynamic_cast<const ThySpmdVecSpaceBase>(Teuchos::rcpFromPtr(Y_inout)->range());
        
        TEUCHOS_TEST_FOR_EXCEPTION(mpi_vs == Teuchos::null, std::logic_error, "Failed to cast Thyra::VectorSpaceBase to Thyra::SpmdVectorSpaceBase.");
        const LocalOrdinal localOffset = ( mpi_vs != Teuchos::null ? mpi_vs->localOffset() : 0 );
        const LocalOrdinal localSubDim = ( mpi_vs != Teuchos::null ? mpi_vs->localSubDim() : Teuchos::rcpFromPtr(Y_inout)->range()->dim() );
        
        RCP<Thyra::DetachedMultiVectorView<Scalar> > thyData =
        Teuchos::rcp(new Thyra::DetachedMultiVectorView<Scalar>(*Teuchos::rcpFromPtr(Y_inout),Teuchos::Range1D(localOffset,localOffset+localSubDim-1)));
        
        for(size_t j = 0; j <xY
            ->getNumVectors(); ++j) {
            Teuchos::ArrayRCP< const Scalar > xpData = xY->getData(j); // access const data from Xpetra object
            // loop over all local rows
            for(LocalOrdinal i = 0; i < localSubDim; ++i) {
                (*thyData)(i,j) = xpData[i];
            }
        }
        
        
 
    }
    
    
    // private
    
    
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    template<class XpetraOperator_t>
    void FROSchLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::initializeImpl(
                                                                                const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
                                                                                const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
                                                                                const RCP<XpetraOperator_t> &xpetraOperator,
                                                                                bool bIsEpetra,
                                                                                bool bIsTpetra
                                                                                )
    {
        
#ifdef THYRA_DEBUG
        TEUCHOS_ASSERT(nonnull(rangeSpace));
        TEUCHOS_ASSERT(nonnull(domainSpace));
        TEUCHOS_ASSERT(nonnull(xpetraOperator));
#endif
        rangeSpace_ = rangeSpace;
        domainSpace_ = domainSpace;
        xpetraOperator_ = xpetraOperator;
        bIsEpetra_ = bIsEpetra;
        bIsTpetra_ = bIsTpetra;
    }
    
    
} // namespace Thyra


#endif  // THYRA_XPETRA_LINEAR_OP_HPP
