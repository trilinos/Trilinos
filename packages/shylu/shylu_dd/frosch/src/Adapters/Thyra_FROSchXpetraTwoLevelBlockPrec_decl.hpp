#ifndef THYRA_FROSCHXPETRATWOLEVELBLOCKPREC_DECL
#define THYRA_FROSCHXPETRATWOLEVELBLOCKPREC_DECL


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

#include <Xpetra_ThyraUtils.hpp>

#include <FROSch_TwoLevelBlockPreconditioner_def.hpp>
#include "Thyra_FROSchLinearOp_def.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"

//#include "Kokkos_DefaultNode.hpp"

namespace Thyra {
    
    template <class SC = Xpetra::Operator<>::scalar_type,
    class LO = typename Xpetra::Operator<SC>::local_ordinal_type,
    class GO = typename Xpetra::Operator<SC,LO>::global_ordinal_type,
    class NO = typename Xpetra::Operator<SC,LO,GO>::node_type>
    class Thyra_FROSchXpetraTwoLevelBlockPrec:public Thyra::PreconditionerFactoryBase<SC>{
        public:
        
 
        typedef Teuchos::RCP<const Teuchos::Comm<int> > CommPtr;
        
        typedef Xpetra::Map<LO,GO,NO> Map;
        typedef Teuchos::RCP<Map> MapPtr;
        typedef Teuchos::RCP<const Map> ConstMapPtr;
        typedef Teuchos::ArrayRCP<MapPtr> MapPtrVecPtr;
        typedef Teuchos::ArrayRCP<MapPtrVecPtr> MapPtrVecPtr2D;
        
        typedef Xpetra::Matrix<SC,LO,GO,NO> CrsMatrix;
        typedef Teuchos::RCP<CrsMatrix> CrsMatrixPtr;
        
        typedef Xpetra::MultiVector<SC,LO,GO,NO> MultiVector;
        typedef Teuchos::RCP<MultiVector> MultiVectorPtr;
        
        typedef Teuchos::RCP<Teuchos::ParameterList> ParameterListPtr;

        typedef unsigned UN;
        
        typedef Teuchos::ArrayRCP<GO> GOVecPtr;
        
        typedef Teuchos::ArrayRCP<SC> SCVecPtr;
        
        typedef Teuchos::ArrayRCP<UN> UNVecPtr;
        
        typedef Teuchos::ArrayRCP<LO> LOVecPtr;
        
        typedef Teuchos::ArrayRCP<GOVecPtr> GOVecPtr2D;
        
        typedef Teuchos::Array<GO>          GOVec;
        typedef Teuchos::Array<GOVec>       GOVec2D;
        typedef Teuchos::ArrayRCP<DofOrdering> DofOrderingVecPtr;
        
        //Constructor
        Thyra_FROSchXpetraTwoLevelBlockPrec();
        
        //Overridden from PreconditionerFactory Base
        bool isCompatible(const LinearOpSourceBase<SC>& fwdOp) const;
        
        Teuchos::RCP<PreconditionerBase<SC> > createPrec() const;
        
        void initializePrec(const Teuchos::RCP<const LinearOpSourceBase<SC> >& fwdOpSrc,
                            PreconditionerBase<SC>* prec,
                            const ESupportSolveUse supportSolveUse
                            ) const;
        
        void uninitializePrec(PreconditionerBase<SC>* prec,
                              Teuchos::RCP<const LinearOpSourceBase<SC> >* fwdOp,
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


