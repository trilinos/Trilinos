#ifndef PYTRILINOS2_MUELU_ETI
#define PYTRILINOS2_MUELU_ETI

#include <MueLu_CreateTpetraPreconditioner.hpp>

#define BINDER_MUELU_CREATETPETRAPRECONDITIONER_INSTANT_2(SCALAR,LO,GO,NO) \
  template Teuchos::RCP<MueLu::TpetraOperator<SCALAR, LO, GO, NO> > CreateTpetraPreconditioner<SCALAR, LO, GO, NO>(const Teuchos::RCP<Tpetra::Operator<SCALAR, LO, GO, NO> > &inA, Teuchos::ParameterList& inParamList);

namespace MueLu {

    template <typename T>
    void initiate(T) {};

  BINDER_MUELU_CREATETPETRAPRECONDITIONER_INSTANT_2(double, int, long long, Tpetra::KokkosClassic::DefaultNode::DefaultNodeType)
}

#endif // PYTRILINOS2_MUELU_ETI
