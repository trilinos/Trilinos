#ifndef FROSch_EPETRAOPERATOR_DEF_HPP
#define FROSch_EPETRAOPERATOR_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_EpetraMultiVector.hpp>

#include "FROSch_EpetraOp_decl.hpp"

namespace FROSch {
    
    FROSch_EpetraOperator::FROSch_EpetraOperator(RCP<Matrix> K,
                                                 ParameterListPtr paramList):
    TwoLevelPrec_(new FROSch::TwoLevelPreconditioner<SC,LO,GO,NO>(K,paramList)),
    paramList_(paramList)
    {
       
        RCP< const Teuchos::Comm< int > > Comm = K->getRowMap()->getComm();
        Comm->barrier();Comm->barrier();Comm->barrier();
        if (Comm->getRank()== 0) {
            std::cout<<"ParamList FROSchEpetraOp\n";
        }
        //sublist(paramList,"TwoLevelPreconditioner")->print(std::cout);
        RCP<FancyOStream> fancy  = fancyOStream(rcpFromRef(std::cout));
        TwoLevelPrec_->describe(*fancy,Teuchos::VERB_EXTREME);
        
        Comm->barrier();Comm->barrier();Comm->barrier();
        if (Comm->getRank()== 0) {
            std::cout<<"-------------------------------------------\n";
            //para->print();
            
        }
        
        TwoLevelPrec_->initialize();
        TwoLevelPrec_->compute();
        
    }
    
    int FROSch_EpetraOperator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
        try {
            
            RCP<const Teuchos::Comm<int> > Comm = rcp(new MpiComm<int> (MPI_COMM_WORLD));
            
            //Teuchos::RCP<ParameterList> para = TwoLevelPrec_->getParameterList();
            //Comm->barrier();Comm->barrier();Comm->barrier();
           if(Comm->getRank() == 0){
                {
                std::cout<<"Print now Parameters\n";
                //para->print();
                
            }}
            
            //TwoLevelPrec_->setParameterList(paramList_);
            //TwoLevelPrec_->description();
            
            Comm->barrier();Comm->barrier();Comm->barrier();
            RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(std::cout));
            
            //TwoLevelPrec_->describe(*fancy,Teuchos::VERB_EXTREME);
            // There is no rcpFromRef(const T&), so we need to do const_cast
            
            const Xpetra::EpetraMultiVectorT<GO,NO> eX(rcpFromRef(const_cast<Epetra_MultiVector&>(X)));
            /*Comm->barrier();Comm->barrier();Comm->barrier();
            if(Comm->getRank() == 0)std::cout<<"FROSch Epetra line 20\n";*/
            const RCP<Xpetra::MultiVector<SC,LO,GO,NO> > xX= rcpFromRef(const_cast<Xpetra::EpetraMultiVectorT<GO,NO>&> (eX));
            /*Comm->barrier();Comm->barrier();Comm->barrier();
            if(Comm->getRank() == 0)std::cout<<"FROSch Epetra line 23\n";*/
            Xpetra::EpetraMultiVectorT<GO,NO> eY(rcpFromRef(Y));
            RCP<Xpetra::MultiVector<SC,LO,GO,NO> > xY = rcpFromRef(eY);
            /*Comm->barrier();Comm->barrier();Comm->barrier();
            if(Comm->getRank() == 0)std::cout<<"FROSch Epetra line 27\n";*/
            // Generally, we assume two different vectors, but AztecOO uses a single vector
            if (X.Values() == Y.Values()) {
                /*Comm->barrier();Comm->barrier();Comm->barrier();
                if(Comm->getRank() == 0)std::cout<<"FROSch Epetra line 31\n";*/
                // X and Y point to the same memory, use an additional vector
                RCP<Xpetra::EpetraMultiVectorT<GO,NO> > tmpY = Teuchos::rcp(new Xpetra::EpetraMultiVectorT<GO,NO>(eY.getMap(), eY.getNumVectors()));
                /*Comm->barrier();Comm->barrier();Comm->barrier();
                if(Comm->getRank() == 0)std::cout<<"FROSch Epetra line 35\n";*/
                RCP<Xpetra::MultiVector<SC,LO,GO,NO> > xtmpY = tmpY;
                // InitialGuessIsZero in MueLu::Hierarchy.Iterate() does not zero out components, it
                // only assumes that user provided an already zeroed out vector
                bool usePreconOnly = true;
                tmpY->putScalar(0.0);
                // apply one V-cycle as preconditioner
                Comm->barrier();Comm->barrier();Comm->barrier();
                if(Comm->getRank() == 0)std::cout<<"FROSch Epetra Apply Two Level If\n";
                TwoLevelPrec_->apply(*xX, *xtmpY,NO_TRANS,Teuchos::ScalarTraits<SC>::one(),Teuchos::ScalarTraits<SC>::zero());
                Comm->barrier();Comm->barrier();Comm->barrier();
                if(Comm->getRank() == 0)std::cout<<"FROSch Epetra line 44\n";
                // deep copy solution from MueLu
                xY->update(1.0, *xtmpY, 0.0);
            } else {
                // X and Y point to different memory, pass the vectors through
               /* Comm->barrier();Comm->barrier();Comm->barrier();
                if(Comm->getRank() == 0)std::cout<<"FROSch Epetra line 50\n";*/
                bool usePreconOnly = true;
                xY->putScalar(0.0);
                Comm->barrier();Comm->barrier();Comm->barrier();
                if(Comm->getRank() == 0)std::cout<<"FROSch Epetra Two Level Apply else\n";
                
               TwoLevelPrec_->apply(*xX, *xY,NO_TRANS,Teuchos::ScalarTraits<SC>::one(),Teuchos::ScalarTraits<SC>::zero());
                Comm->barrier();Comm->barrier();Comm->barrier();
                if(Comm->getRank() == 0)std::cout<<"FROSch Epetra after TwoLevel Apply\n";
                
                
            }
            
        } catch (std::exception& e) {
            //TODO: error msg directly on std::cerr?
            std::cerr << "Caught an exception in FROSch::EpetraOperator::ApplyInverse():" << std::endl
            << e.what() << std::endl;
            return -1;
        }
        return 0;
    }

    const Epetra_Comm& FROSch_EpetraOperator::Comm() const {
        RCP<Matrix> A = TwoLevelPrec_->getCrsMatrix();
        
        RCP<const Xpetra::CrsMatrixWrap<SC,LO,GO,NO> > crsOp = rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<SC,LO,GO,NO> >(A);
        if (crsOp == Teuchos::null)
            std::cerr<<"Caution Bad Cast!\n";
        const RCP<const Xpetra::EpetraCrsMatrixT<GO,NO>> &tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::EpetraCrsMatrixT<GO,NO>>(crsOp->getCrsMatrix());
        if (tmp_ECrsMtx == Teuchos::null)
            std::cerr<<"Caution Bad Cast!\n";
        return tmp_ECrsMtx->getEpetra_CrsMatrixNonConst()->Comm();
    /*}
        
        RCP<const Xpetra::EpetraCrsMatrixT<GO,NO>> tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::EpetraCrsMatrixT<GO,NO> >(A);
        RCP<Epetra_CrsMatrix> epA = tmp_ECrsMtx->getEpetra_CrsMatrixNonConst();
        return epA->Comm();*/
    }
 
    const Epetra_Map& FROSch_EpetraOperator::OperatorDomainMap() const {
        RCP<Matrix> A = TwoLevelPrec_->getCrsMatrix();
        
        RCP<const Xpetra::CrsMatrixWrap<SC,LO,GO,NO> > crsOp = rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<SC,LO,GO,NO> >(A);
        if (crsOp == Teuchos::null)
            std::cerr<<"Caution Bad Cast!\n";
        const RCP<const Xpetra::EpetraCrsMatrixT<GO,NO>> &tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::EpetraCrsMatrixT<GO,NO>>(crsOp->getCrsMatrix());
        if (tmp_ECrsMtx == Teuchos::null)
            std::cerr<<"Caution Bad Cast!\n";
        return tmp_ECrsMtx->getEpetra_CrsMatrixNonConst()->DomainMap();
    }


    const Epetra_Map & FROSch_EpetraOperator::OperatorRangeMap() const {
        RCP<Matrix> A = TwoLevelPrec_->getCrsMatrix();

        
        RCP<const Xpetra::CrsMatrixWrap<SC,LO,GO,NO> > crsOp = rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<SC,LO,GO,NO> >(A);
        if (crsOp == Teuchos::null)
            std::cerr<<"Caution Bad Cast!\n";
        const RCP<const Xpetra::EpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::EpetraCrsMatrix>(crsOp->getCrsMatrix());
        if (tmp_ECrsMtx == Teuchos::null)
            std::cerr<<"Caution Bad Cast!\n";
        return tmp_ECrsMtx->getEpetra_CrsMatrixNonConst()->RangeMap();
    }

        
}
#endif
