#include "MiniEM_Utils.hpp"

namespace mini_em {

  void writeOut(const std::string & s,const Thyra::LinearOpBase<double> & op)
  {
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;
    using NT = panzer::TpetraNodeType;
    const RCP<const Thyra::TpetraLinearOp<double,int,panzer::GlobalOrdinal,NT> > tOp = rcp_dynamic_cast<const Thyra::TpetraLinearOp<double,int,panzer::GlobalOrdinal,NT> >(Teuchos::rcpFromRef(op));
    const RCP<const Thyra::EpetraLinearOp> eOp = rcp_dynamic_cast<const Thyra::EpetraLinearOp>(Teuchos::rcpFromRef(op));
    if(tOp != Teuchos::null) {
      *Teko::getOutputStream() << "Dumping matrix \'" << s << "\'" << std::endl;
      const RCP<const Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> > crsOp = rcp_dynamic_cast<const Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >(tOp->getConstTpetraOperator(),true);
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >::writeMapFile(("rowmap_"+s).c_str(),*(crsOp->getRowMap()));
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >::writeMapFile(("colmap_"+s).c_str(),*(crsOp->getColMap()));
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >::writeMapFile(("domainmap_"+s).c_str(),*(crsOp->getDomainMap()));
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >::writeMapFile(("rangemap_"+s).c_str(),*(crsOp->getRangeMap()));
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >::writeSparseFile(s.c_str(),crsOp);
    } else if (eOp != Teuchos::null) {
      *Teko::getOutputStream() << "Dumping matrix \'" << s << "\'" << std::endl;
      const RCP<const Epetra_CrsMatrix> crsOp = rcp_dynamic_cast<const Epetra_CrsMatrix>(eOp->epetra_op(),true);
      EpetraExt::BlockMapToMatrixMarketFile(("rowmap_"+s).c_str(), crsOp->RowMap());
      EpetraExt::BlockMapToMatrixMarketFile(("colmap_"+s).c_str(), crsOp->ColMap());
      EpetraExt::BlockMapToMatrixMarketFile(("domainmap_"+s).c_str(), crsOp->DomainMap());
      EpetraExt::BlockMapToMatrixMarketFile(("rangemap_"+s).c_str(), crsOp->RangeMap());
      EpetraExt::RowMatrixToMatrixMarketFile(s.c_str(), *crsOp);
    } else
      TEUCHOS_ASSERT(false);
  }


  void describeMatrix(const std::string & s,const Thyra::LinearOpBase<double> & op,Teuchos::RCP<Teuchos::FancyOStream> out)
  {
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;
    using NT = Tpetra::Map<>::node_type;
    if (out!=Teuchos::null) {
      const RCP<const Thyra::TpetraLinearOp<double,int,panzer::GlobalOrdinal,NT> > tOp = rcp_dynamic_cast<const Thyra::TpetraLinearOp<double,int,panzer::GlobalOrdinal,NT> >(Teuchos::rcpFromRef(op));
      const RCP<const Thyra::EpetraLinearOp > eOp = rcp_dynamic_cast<const Thyra::EpetraLinearOp>(Teuchos::rcpFromRef(op));
      if(tOp != Teuchos::null) {
        const RCP<const Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> > crsOp = rcp_dynamic_cast<const Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >(tOp->getConstTpetraOperator(),true);
        *out << "\nDebug: " << s << std::endl;
        crsOp->describe(*out,Teuchos::VERB_MEDIUM);
      } else if (eOp != Teuchos::null) {
        const RCP<const Epetra_CrsMatrix> crsOp = rcp_dynamic_cast<const Epetra_CrsMatrix>(eOp->epetra_op(),true);
        *out << "\nDebug: " << s << std::endl;
        // crsOp->describe(*out,Teuchos::VERB_MEDIUM);
      } else
        TEUCHOS_ASSERT(false);
    }
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > get_Tpetra_CrsMatrix(const Thyra::LinearOpBase<double> & op) {
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;
    const RCP<const Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tOp = rcp_dynamic_cast<const Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(Teuchos::rcpFromRef(op),true);
    RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsOp = rcp_dynamic_cast<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(tOp->getConstTpetraOperator(),true);
    return crsOp;
  }


  Teuchos::RCP<const Epetra_CrsMatrix> get_Epetra_CrsMatrix(const Thyra::LinearOpBase<double> & op) {
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;
    const RCP<const Thyra::EpetraLinearOp> eOp = rcp_dynamic_cast<const Thyra::EpetraLinearOp>(Teuchos::rcpFromRef(op),true);
    RCP<const Epetra_CrsMatrix> crsOp = rcp_dynamic_cast<const Epetra_CrsMatrix>(eOp->epetra_op(),true);
    return crsOp;
  }


  Teuchos::RCP<const Epetra_CrsMatrix> get_Epetra_CrsMatrix(const Thyra::DiagonalLinearOpBase<double> & op, const Epetra_Comm& comm) {
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;
    RCP<const Epetra_Map> map = Thyra::get_Epetra_Map(*op.range(),Teuchos::rcpFromRef(comm));
    int nodeNumElements = map->NumMyElements();
    RCP<Epetra_CrsMatrix> crsMatrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*map,*map,1,true));

    RCP<const Thyra::VectorBase<double> > diag = op.getDiag();
    RTOpPack::SubVectorView<double> view;
    diag->acquireDetachedView(Thyra::Range1D(),&view);

    for (int i = 0; i < nodeNumElements; i++) {
      int err = crsMatrix->InsertMyValues(i, 1, &(view[i]), &i);
      TEUCHOS_ASSERT(err==0);
    }

    diag->releaseDetachedView(&view);
    crsMatrix->FillComplete();
    return crsMatrix;
  }

}
