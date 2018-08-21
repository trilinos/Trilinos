#ifndef THYRA_FROSCH_TWOLEVELPRECONDITIONER_FACTORY_DECL_HPP
#define THYRA_FROSCH_TWOLEVELPRECONDITIONER_FACTORY_DECL_HPP


// Stratimikos needs Thyra, so we don't need special guards for Thyra here
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_XpetraLinearOp.hpp"
#ifdef HAVE_MUELU_TPETRA
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#endif
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"

#include "Teuchos_Ptr.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Time.hpp"

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_ThyraUtils.hpp>
//PUT FROSch includes here

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Array.hpp"

#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Version.h>
#include <Epetra_config.h>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_EpetraMap.hpp>


//#include <FROSch_GDSWPreconditioner_def.hpp>
#include <FROSch_TwoLevelPreconditioner_def.hpp>
#include <FROSch_TwoLevelPreconditioner_decl.hpp>

#include <FROSch_Tools_def.hpp>
#include <FROSch_TwoLevelPreconditioner_def.hpp>
//

#include "Thyra_PreconditionerFactoryBase.hpp"

#include "Kokkos_DefaultNode.hpp"

namespace Thyra {
    
    template <class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node=KokkosClassic::DefaultNode::DefaultNodeType>
    class FROSch_TwoLevelPreconditionerFactory:public Thyra::PreconditionerFactoryBase<Scalar>{
    public:
        
        //Constructor
        FROSch_TwoLevelPreconditionerFactory();
        
        //Overridden from PreconditionerFactory Base
        bool isCompatible(const LinearOpSourceBase<Scalar>& fwdOp) const;
        
        Teuchos::RCP<PreconditionerBase<Scalar> > createPrec() const;
        
        void initializePrec(const Teuchos::RCP<const LinearOpSourceBase<Scalar> >& fwdOpSrc,
                            PreconditionerBase<Scalar>* prec,
                            const ESupportSolveUse supportSolveUse
                            ) const;
        
        void uninitializePrec(PreconditionerBase<Scalar>* prec,
                              Teuchos::RCP<const LinearOpSourceBase<Scalar> >* fwdOp,
                              ESupportSolveUse* supportSolveUse
                              ) const;
        
        void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList);
        
        Teuchos::RCP<Teuchos::ParameterList>          unsetParameterList();

        Teuchos::RCP<Teuchos::ParameterList>          getNonconstParameterList();

        Teuchos::RCP<const Teuchos::ParameterList>    getParameterList() const;

        Teuchos::RCP<const Teuchos::ParameterList>    getValidParameters() const;
        
        std::string description() const;
    private:
        Teuchos::RCP<Teuchos::ParameterList> paramList_;


        
    };
    
#ifdef HAVE_FROSCH_EPETRA
    
    //special fpr Epetra
    template<>
    class FROSch_TwoLevelPreconditionerFactory<double,int,int, Xpetra::EpetraNode>:public PreconditionerFactoryBase<double>{
    public:
        typedef double Scalar;
        typedef int LocalOrdinal;
        typedef int GlobalOrdinal;
        typedef Xpetra::EpetraNode Node;
        
        //Constructor
        FROSch_TwoLevelPreconditionerFactory():paramList_(rcp(new ParameterList())) { };
        
        //Overridden from PreconditionerFactory Base
        bool isCompatible(const LinearOpSourceBase<Scalar>& fwdOpSrc) const{
            const RCP<const LinearOpBase<Scalar> > fwdOp = fwdOpSrc.getOp();
        
            if(Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isBlockedOperator(fwdOp)) return true;
            
            return false;
            
        };
        
        Teuchos::RCP<PreconditionerBase<Scalar> > createPrec() const {
            return Teuchos::rcp(new DefaultPreconditioner<Scalar>);
        };
        
        
        voidd initializePrec(const Teuchos::RCP<const LinearOpSourceBase<Scalar> >& fwdOp,
         PreconditionerBase<Scalar>* prec,
         const ESupportSolveUse supportSolveUse
         ) const{
            
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
            
            //PreCheck
            TEUCHOS_ASSERT(Teuchos::nonnull(fwdOpSrc));
            TEUCHOS_ASSERT(this->isCompatible(*fwdOpSrc));
            TEUCHOS_ASSERT(prec);
            
            // Create a copy, as we may remove some things from the list
            ParameterList paramList = *paramList_;
            
            // Retrieve wrapped concrete Xpetra matrix from FwdOp
            const RCP<const ThyLinOpBase> fwdOp = fwdOpSrc->getOp();
            TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(fwdOp));
            
            // Check whether it is Epetra/Tpetra
            bool bIsEpetra  = XpThyUtils::isEpetra(fwdOp);
            bool bIsTpetra  = XpThyUtils::isTpetra(fwdOp);
            bool bIsBlocked = XpThyUtils::isBlockedOperator(fwdOp);
            TEUCHOS_TEST_FOR_EXCEPT((bIsEpetra == true  && bIsTpetra == true));
            TEUCHOS_TEST_FOR_EXCEPT((bIsEpetra == bIsTpetra) && bIsBlocked == false);
            TEUCHOS_TEST_FOR_EXCEPT((bIsEpetra != bIsTpetra) && bIsBlocked == true);
            
            RCP<XpMat> A = Teuchos::null;
            if(bIsBlocked) {
                Teuchos::RCP<const Thyra::BlockedLinearOpBase<Scalar> > ThyBlockedOp =
                Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<Scalar> >(fwdOp);
                TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(ThyBlockedOp));
                
                TEUCHOS_TEST_FOR_EXCEPT(ThyBlockedOp->blockExists(0,0)==false);
                
                Teuchos::RCP<const LinearOpBase<Scalar> > b00 = ThyBlockedOp->getBlock(0,0);
                TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(b00));
                
                RCP<const XpCrsMat > xpetraFwdCrsMat00 = XpThyUtils::toXpetra(b00);
                TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xpetraFwdCrsMat00));
                
                // TwoLevelPreconditioner needs Xpetra CRS Matrix ans input object as input
                RCP<XpCrsMat> xpetraFwdCrsMatNonConst00 = Teuchos::rcp_const_cast<XpCrsMat>(xpetraFwdCrsMat00);
                TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xpetraFwdCrsMatNonConst00));
                
                // wrap the forward operator as an Xpetra::Matrix that TwoLevelPreconditioner can work with
                RCP<XpMat> A00 = rcp(new Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>(xpetraFwdCrsMatNonConst00));
                TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(A00));
                
                RCP<const XpMap> rowmap00 = A00->getRowMap();
                RCP< const Teuchos::Comm< int > > comm = rowmap00->getComm();
                
                // create a Xpetra::BlockedCrsMatrix which derives from Xpetra::Matrix that FROSCH can work with
                RCP<XpCrsMat> bMat = Teuchos::rcp(new XpCrsMat(ThyBlockedOp, comm));
                TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(bMat));
                
                // save blocked matrix
                A = bMat;
            } else {
                RCP<const XpCrsMat > xpetraFwdCrsMat = XpThyUtils::toXpetra(fwdOp);
                TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xpetraFwdCrsMat));
                
                // FROSCH needs a non-const object as input
                RCP<XpCrsMat> xpetraFwdCrsMatNonConst = Teuchos::rcp_const_cast<XpCrsMat>(xpetraFwdCrsMat);
                TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xpetraFwdCrsMatNonConst));
                
                // wrap the forward operator as an Xpetra::Matrix that MueLu can work with
                A = rcp(new Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>(xpetraFwdCrsMatNonConst));
            }
            TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(A));
            
            // Retrieve concrete preconditioner object
            const Teuchos::Ptr<DefaultPreconditioner<Scalar> > defaultPrec = Teuchos::ptr(dynamic_cast<DefaultPreconditioner<Scalar> *>(prec));
            TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(defaultPrec));
            
            // extract preconditioner operator
            RCP<ThyLinOpBase> thyra_precOp = Teuchos::null;
            thyra_precOp = rcp_dynamic_cast<Thyra::LinearOpBase<Scalar> >(defaultPrec->getNonconstUnspecifiedPrecOp(), true);
            
            //Later needs to be a utility ExtractCoordinatesFromParameterList
            //not implemented yet
            // FROSCH::Tools<SC,LO,GO,Node>::ExtractCoordinatesFromParameterList(paramList);
            
            RCP<FROSch::TwoLevelPreconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> > TwoLevelPrec (new FROSch::TwoLevelPreconditioner(A,paramList));
            
            //Initialize-> Only Works for laplce (cause defaults are used) and compute
            TwoLevelPrec->inialize();
            TwoLevelPrec->compute();
            
            //Wrap tp thyra
            RCP<ThyLinOpBase > thyraPrecOp = Teuchos::null;
            
            if(bIsBlocked) {
                TEUCHOS_TEST_FOR_EXCEPT(Teuchos::nonnull(thyraPrecOp));
                
                //create const Operator
                const RCP<FROSch::TwoLevelPreconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> > frosch_op = rcp(new FROSch::TwoLevelPreconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node>(TwoLevelPrec));
                
                RCP<const VectorSpaceBase<Scalar> > thyraRangeSpace  = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(frosch_op->getRangeMap());
                
                RCP<const VectorSpaceBase<Scalar> > thyraDomainSpace = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(frosch_op->getDomainMap());
                
                RCP <Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > xpOp = Teuchos::rcp_dynamic_cast<Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(frosch_op);
                thyraPrecOp = Thyra::xpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(thyraRangeSpace, thyraDomainSpace,xpOp);
                
                
            }
            TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraPrecOp));
            
            defaultPrec->initializeUnspecified(thyraPrecOp);
            
        };
        
        void uninitializePrec(PreconditionerBase<Scalar>* prec, RCP<const LinearOpSourceBase<Scalar> >* fwdOp, ESupportSolveUse* supportSolveUse) const {
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
        };
        
        void setParameterList(RCP<ParameterList> const & paramList){
            TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(paramList));
            paramList_ = paramList;
        };
        
        RCP<ParameterList>getNonconstParameterList(){
            return paramList_;
        };
        
        RCP<ParameterList> getParameterList() const {
            return paramList_;
        };
        
        RCP<ParameterList> getValidParameters() const {
            static RCP<const ParameterList> validPL;
            
            if (Teuchos::is_null(validPL))
                validPL = rcp(new ParameterList());
            
            return validPL;
        };
        
        std::string description() const {
            return "Thyra::FROSch_TwoLevelPreconditionerFactory";
        };
        
    private:
        Teuchos::RCP<Teuchos::ParameterList> paramList_;
    };
    
    }
#endif
}
#endif

