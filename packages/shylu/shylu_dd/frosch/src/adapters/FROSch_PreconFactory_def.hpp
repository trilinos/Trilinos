#include <Xpetra_EpetraMap.hpp>


//#include <FROSch_GDSWPreconditioner_def.hpp>
#include <FROSch_TwoLevelPreconditioner_def.hpp>
#include <FROSch_TwoLevelPreconditioner_decl.hpp>

#include <FROSch_Tools_def.hpp>
#include <FROSch_TwoLevelPreconditioner_def.hpp>
//

#include "FROSch_PreconFactory_decl.hpp"

#include "Kokkos_DefaultNode.hpp"

namespace Thyra {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ParameterList;
    
    //Constructor----------------
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    FROSch_PreconFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::FROSch_PreconFactory()
    {
        paramList_ = rcp(new Teuchos::ParameterList());
    }
    //----------------------------
    
    //-----------------------------------------------------------
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    bool FROSch_PreconFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isCompatible(const LinearOpSourceBase<Scalar>& fwdOpSrc) const
    {
        std::cout<<"Only Epetra works fine->Insert CHECK\n";
        return true;
    }
    //-----------------------------------------------------------
    
    //--------------------------------------------------------------
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal , class Node>
    RCP<PreconditionerBase<Scalar> >FROSch_PreconFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::createPrec() const
    {
        return Teuchos::rcp(new DefaultPreconditioner<Scalar>);
        
    }
    //-------------------------------------------------------------
    
    //-------------------------------------------------------------
    template<class Scalar, class LocalOrdinal , class GlobalOrdinal, class Node>
    void FROSch_PreconFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::initializePrec
    (const Teuchos::RCP<const LinearOpSourceBase<Scalar> >& fwdOpSrc,
     PreconditionerBase<Scalar>* prec,
     const ESupportSolveUse supportSolveUse
     ) const{
        
        Teuchos::RCP<Teuchos::FancyOStream> fancy = fancyOStream(Teuchos::rcpFromRef(std::cout));
        
        
        using Teuchos::rcp_dynamic_cast;
        //Some Typedefs
        typedef Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                     XpMap;
        typedef Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>      XpOp;
        typedef Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>       XpThyUtils;
        typedef Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>        XpCrsMat;
        typedef Xpetra::BlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> XpBlockedCrsMat;
        typedef Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>           XpMat;
        typedef Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>      XpMultVec;
        typedef Xpetra::MultiVector<double,LocalOrdinal,GlobalOrdinal,Node>      XpMultVecDouble;
        typedef Thyra::LinearOpBase<Scalar>                                      ThyLinOpBase;
        
        TEUCHOS_ASSERT(Teuchos::nonnull(fwdOpSrc));
        TEUCHOS_ASSERT(this->isCompatible(*fwdOpSrc));
        TEUCHOS_ASSERT(prec);
        
        // Create a copy, as we may remove some things from the list
        ParameterList paramList = *paramList_;
        
        // Retrieve wrapped concrete Xpetra matrix from FwdOp
        const RCP<const ThyLinOpBase> fwdOp = fwdOpSrc->getOp();
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(fwdOp));
        
        //Check fpr Epetra
        bool bIsEpetra  = XpThyUtils::isEpetra(fwdOp);
        
        RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > K = Teuchos::null;
        
        RCP<const XpCrsMat > xpetraFwdCrsMat = XpThyUtils::toXpetra(fwdOp);
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xpetraFwdCrsMat));
        
        // FROSCH needs a non-const object as input
        RCP<XpCrsMat> xpetraFwdCrsMatNonConst = Teuchos::rcp_const_cast<XpCrsMat>(xpetraFwdCrsMat);
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xpetraFwdCrsMatNonConst));
        
        K = rcp(new Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>(xpetraFwdCrsMatNonConst));
        
        // Retrieve concrete preconditioner object
        const Teuchos::Ptr<DefaultPreconditioner<Scalar> > defaultPrec = Teuchos::ptr(dynamic_cast<DefaultPreconditioner<Scalar> *>(prec));
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(defaultPrec));
        
        // extract preconditioner operator
        RCP<ThyLinOpBase> thyra_precOp = Teuchos::null;
        thyra_precOp = rcp_dynamic_cast<Thyra::LinearOpBase<Scalar> >(defaultPrec->getNonconstUnspecifiedPrecOp(), true);
        
        if (bIsEpetra){
            std::cout<<"Korrekt\n";
        }
        else{
            std::cout<<"Did not use Epetra this might cause problem\n";
        }
        
        


    }
    //-------------------------------------------------------------
    
    //-------------------------------------------------------------
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void FROSch_PreconFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    uninitializePrec(PreconditionerBase<Scalar>* prec, RCP<const LinearOpSourceBase<Scalar> >* fwdOp, ESupportSolveUse* supportSolveUse) const {
        TEUCHOS_ASSERT(prec);
        
        // Retrieve concrete preconditioner object
        const Teuchos::Ptr<DefaultPreconditioner<Scalar> > defaultPrec = Teuchos::ptr(dynamic_cast<DefaultPreconditioner<Scalar> *>(prec));
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(defaultPrec));
        
        if (fwdOp) {
            // TODO: Implement properly instead of returning default value
            *fwdOp = Teuchos::null;
        }
        
        if (supportSolveUse) {
            // TODO: Implement properly instead of returning default value
            *supportSolveUse = Thyra::SUPPORT_SOLVE_UNSPECIFIED;
        }
        
        defaultPrec->uninitialize();
    }
    //-------------------------------------------------------------
    
    //-------------------------------------------------------------
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void FROSch_PreconFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setParameterList(RCP<ParameterList> const & paramList){
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(paramList));
        paramList_ = paramList;
    }
    //-------------------------------------------------------------
    
    //------------------------------------------------------------------
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    RCP<ParameterList> FROSch_PreconFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getNonconstParameterList(){
        return paramList_;
    }
    
    //-----------------------------------------------------------------------
    
    
    //-----------------------------------------------------------------------
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    RCP<const ParameterList>
    FROSch_PreconFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getParameterList() const {
        return paramList_;
    }
    //--------------------------------------------------------------------
    
    //--------------------------------------------------------------------
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    RCP<const ParameterList> FROSch_PreconFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getValidParameters() const {
        static RCP<const ParameterList> validPL;
        
        if (Teuchos::is_null(validPL))
            validPL = rcp(new ParameterList());
        
        return validPL;
    }
    //-----------------------------------------------------------------------
    
    //-----------------------------------------------------------------------
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    std::string FROSch_PreconFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::description() const {
        return "Thyra::FROSch_TwoLevelPreconditionerFactory";
    }
    //--------------------------------------------------------------------------
    
    //--------------------------------------------------------------------------
    template<class Scalar, class LocalOrdinal,class GlobalOrdinal, class Node>
    RCP<ParameterList> FROSch_PreconFactory<Scalar, LocalOrdinal,GlobalOrdinal,Node>::unsetParameterList(){
        RCP<ParameterList> savedParamList = paramList_;
        paramList_ = Teuchos::null;
        return savedParamList;
        
    }
    

}
