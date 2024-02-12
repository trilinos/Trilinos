
#ifndef PACKAGES_MUELU_ADAPTERS_AZTECOO_MUELU_AZTECEPETRAOPERATOR_CPP_
#define PACKAGES_MUELU_ADAPTERS_AZTECOO_MUELU_AZTECEPETRAOPERATOR_CPP_

#include "Xpetra_EpetraMultiVector.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_EpetraCrsMatrix.hpp"

#include "MueLu_config.hpp"  // for HAVE_MUELU_DEBUG
#include "MueLu_RefMaxwell.hpp"
#include "MueLu_Exceptions.hpp"

#include "MueLu_AztecEpetraOperator.hpp"

#if defined(HAVE_MUELU_SERIAL) and defined(HAVE_MUELU_EPETRA)

namespace MueLu {

int AztecEpetraOperator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
  try {
    // There is no rcpFromRef(const T&), so we need to do const_cast
    const Xpetra::EpetraMultiVectorT<GO, NO> eX(Teuchos::rcpFromRef(const_cast<Epetra_MultiVector&>(X)));
    Xpetra::EpetraMultiVectorT<GO, NO> eY(Teuchos::rcpFromRef(Y));
    // Generally, we assume two different vectors, but AztecOO uses a single vector
    if (X.Values() == Y.Values()) {
      // X and Y point to the same memory, use an additional vector
      Teuchos::RCP<Xpetra::EpetraMultiVectorT<GO, NO>> tmpY = Teuchos::rcp(new Xpetra::EpetraMultiVectorT<GO, NO>(eY.getMap(), eY.getNumVectors()));
      tmpY->putScalar(0.0);
      xOp_->apply(eX, *tmpY);
      // deep copy solution from MueLu
      eY.update(1.0, *tmpY, 0.0);
    } else {
      // X and Y point to different memory, pass the vectors through
      eY.putScalar(0.0);
      xOp_->apply(eX, eY);
    }

  } catch (std::exception& e) {
    // TODO: error msg directly on std::cerr?
    std::cerr << "Caught an exception in MueLu::AztecEpetraOperator::ApplyInverse():" << std::endl
              << e.what() << std::endl;
    return -1;
  }
  return 0;
}

const Epetra_Comm& AztecEpetraOperator::Comm() const {
  const Epetra_Map emap = Xpetra::toEpetra(xOp_->getDomainMap());
  return emap.Comm();
}

const Epetra_Map& AztecEpetraOperator::OperatorDomainMap() const {
  if (Teuchos::rcp_dynamic_cast<MueLu::RefMaxwell<SC, LO, GO, NO>>(xOp_) != Teuchos::null) {
    RCP<Xpetra::Matrix<SC, LO, GO, NO>> A                  = Teuchos::rcp_dynamic_cast<MueLu::RefMaxwell<SC, LO, GO, NO>>(xOp_)->getJacobian();
    RCP<const Xpetra::CrsMatrixWrap<SC, LO, GO, NO>> crsOp = rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<SC, LO, GO, NO>>(A);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    const RCP<const Xpetra::EpetraCrsMatrixT<GO, NO>>& tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::EpetraCrsMatrixT<GO, NO>>(crsOp->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
    return tmp_ECrsMtx->getEpetra_CrsMatrixNonConst()->DomainMap();
  }

  // TAW there is a problem with the Map version of toEeptra leading to code crashes
  //     we probably need to create a new copy of the Epetra_Map
  Teuchos::RCP<const Map> map = xOp_->getDomainMap();
  return Xpetra::toEpetra(map);
}

const Epetra_Map& AztecEpetraOperator::OperatorRangeMap() const {
  if (Teuchos::rcp_dynamic_cast<MueLu::RefMaxwell<SC, LO, GO, NO>>(xOp_) != Teuchos::null) {
    RCP<Xpetra::Matrix<SC, LO, GO, NO>> A                  = Teuchos::rcp_dynamic_cast<MueLu::RefMaxwell<SC, LO, GO, NO>>(xOp_)->getJacobian();
    RCP<const Xpetra::CrsMatrixWrap<SC, LO, GO, NO>> crsOp = rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<SC, LO, GO, NO>>(A);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    const RCP<const Xpetra::EpetraCrsMatrixT<GO, NO>>& tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::EpetraCrsMatrixT<GO, NO>>(crsOp->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
    return tmp_ECrsMtx->getEpetra_CrsMatrixNonConst()->RangeMap();
  }

  // TAW there is a problem with the Map version of toEeptra leading to code crashes
  //     we probably need to create a new copy of the Epetra_Map
  Teuchos::RCP<const Map> map = xOp_->getRangeMap();
  return Xpetra::toEpetra(map);
}

}  // namespace MueLu

#endif /*#if defined(HAVE_MUELU_SERIAL) and defined(HAVE_MUELU_EPETRA)*/

#endif /* PACKAGES_MUELU_ADAPTERS_AZTECOO_MUELU_AZTECEPETRAOPERATOR_CPP_ */
