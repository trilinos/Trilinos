#ifndef PYTRILINOS2_MUELU_ETI
#define PYTRILINOS2_MUELU_ETI

#include <MueLu_CreateTpetraPreconditioner.hpp>

#define BINDER_MUELU_CREATETPETRAPRECONDITIONER_INSTANT_2(SCALAR,LO,GO,NO) \
  template Teuchos::RCP<MueLu::TpetraOperator<SCALAR, LO, GO, NO> > CreateTpetraPreconditioner<SCALAR, LO, GO, NO>(const Teuchos::RCP<Tpetra::Operator<SCALAR, LO, GO, NO> > &inA, Teuchos::ParameterList& inParamList);

#define BINDER_MUELU_CREATETPETRAPRECONDITIONER2_INSTANT_2(SCALAR,LO,GO,NO) \
  template Teuchos::RCP<MueLu::TpetraOperator<SCALAR, LO, GO, NO> > CreateTpetraPreconditioner2<SCALAR, LO, GO, NO>(Tpetra::Operator<SCALAR, LO, GO, NO> &inA, Teuchos::ParameterList& inParamList);

#define BINDER_MUELU_CREATETPETRAPRECONDITIONER3_INSTANT_2(SCALAR,LO,GO,NO) \
  template Teuchos::RCP<MueLu::TpetraOperator<SCALAR, LO, GO, NO> > CreateTpetraPreconditioner3<SCALAR, LO, GO, NO>(const std::shared_ptr<Tpetra::Operator<SCALAR, LO, GO, NO> > &inA, Teuchos::ParameterList& inParamList);

#define BINDER_MUELU_CREATETPETRAPRECONDITIONER4_INSTANT_2(SCALAR,LO,GO,NO) \
  template Teuchos::RCP<MueLu::TpetraOperator<SCALAR, LO, GO, NO> > CreateTpetraPreconditioner4<SCALAR, LO, GO, NO>(const Teuchos::RCP<Tpetra::CrsMatrix<SCALAR, LO, GO, NO> > &inA, Teuchos::ParameterList& inParamList);

#define BINDER_MUELU_CREATETPETRAPRECONDITIONER5_INSTANT_2(SCALAR,LO,GO,NO) \
  template Teuchos::RCP<MueLu::TpetraOperator<SCALAR, LO, GO, NO> > CreateTpetraPreconditioner5<SCALAR, LO, GO, NO>(const Teuchos::RCP<Tpetra::RowMatrix<SCALAR, LO, GO, NO> > &inA, Teuchos::ParameterList& inParamList);

namespace MueLu {

    template <typename T>
    void initiate(T) {};

  BINDER_MUELU_CREATETPETRAPRECONDITIONER_INSTANT_2(double, int, long long, Tpetra::KokkosClassic::DefaultNode::DefaultNodeType)
  BINDER_MUELU_CREATETPETRAPRECONDITIONER2_INSTANT_2(double, int, long long, Tpetra::KokkosClassic::DefaultNode::DefaultNodeType)
  BINDER_MUELU_CREATETPETRAPRECONDITIONER3_INSTANT_2(double, int, long long, Tpetra::KokkosClassic::DefaultNode::DefaultNodeType)
  BINDER_MUELU_CREATETPETRAPRECONDITIONER4_INSTANT_2(double, int, long long, Tpetra::KokkosClassic::DefaultNode::DefaultNodeType)
  BINDER_MUELU_CREATETPETRAPRECONDITIONER5_INSTANT_2(double, int, long long, Tpetra::KokkosClassic::DefaultNode::DefaultNodeType)

}

#endif // PYTRILINOS2_MUELU_ETI
