#ifndef THYRA_FROSCH_XPETRA_FACTORY_DECL_HPP
#define THYRA_FROSCH_XPETRA_FACTORY_DECL_HPP


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
//#include <Epetra_Version.h>
#include <Epetra_config.h>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_EpetraMap.hpp>


//#include <FROSch_GDSWPreconditioner_def.hpp>
#include <FROSch_TwoLevelPreconditioner_def.hpp>
#include "FROSch_XpetraOperator_decl.hpp"
#include "Thyra_FROSchLinearOP_def.hpp"

#include "Thyra_PreconditionerFactoryBase.hpp"

#include "Kokkos_DefaultNode.hpp"

namespace Thyra {
    
    template <class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node=KokkosClassic::DefaultNode::DefaultNodeType>
    class FROSch_XpetraFactory:public Thyra::PreconditionerFactoryBase<Scalar>{
        public:
        
        //Constructor
        FROSch_XpetraFactory();
        
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
}

    

#endif


