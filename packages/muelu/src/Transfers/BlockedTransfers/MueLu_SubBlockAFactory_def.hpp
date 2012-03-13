/*
 * MueLu_SubBlockAFactory_def.hpp
 *
 *  Created on: 02.01.2012
 *      Author: tobias
 */

#ifndef MUELU_SUBBLOCKAFACTORY_DEF_HPP_
#define MUELU_SUBBLOCKAFACTORY_DEF_HPP_


#include "MueLu_SubBlockAFactory_decl.hpp"

#include <Xpetra_Operator.hpp>
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_BlockedCrsOperator.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  SubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SubBlockAFactory(Teuchos::RCP<const FactoryBase> Afact, size_t row, size_t col, LocalOrdinal blksize)
    : Afact_(Afact), row_(row), col_(col), blksize_(blksize)
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  SubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~SubBlockAFactory() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("A",Afact_.get(),this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & currentLevel) const {
    typedef Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> OOperator; //TODO
    typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrixClass; //TODO
    typedef Xpetra::CrsOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsOperatorClass; //TODO
    typedef Xpetra::BlockedCrsOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> BlockedCrsOOperator; //TODO

    RCP<OOperator> Ain = currentLevel.Get< RCP<OOperator> >("A", Afact_.get());
    RCP<BlockedCrsOOperator> bA = Teuchos::rcp_dynamic_cast<BlockedCrsOOperator>(Ain);

    TEUCHOS_TEST_FOR_EXCEPTION(bA==Teuchos::null, Exceptions::BadCast, "MueLu::SubBlockAFactory::Build: input matrix A is not of type BlockedCrsOperator! error.");
    TEUCHOS_TEST_FOR_EXCEPTION(row_>bA->Rows(), Exceptions::RuntimeError, "MueLu::SubBlockAFactory::Build: A.Rows() > rows_! error.");
    TEUCHOS_TEST_FOR_EXCEPTION(col_>bA->Cols(), Exceptions::RuntimeError, "MueLu::SubBlockAFactory::Build: A.Cols() > cols_! error.");
    TEUCHOS_TEST_FOR_EXCEPTION(row_<0, Exceptions::RuntimeError, "MueLu::SubBlockAFactory::Build: row_<0 error.");
    TEUCHOS_TEST_FOR_EXCEPTION(col_<0, Exceptions::RuntimeError, "MueLu::SubBlockAFactory::Build: col_<0 error.");

    Teuchos::RCP<CrsMatrixClass> A = bA->getMatrix(row_, col_);

    Teuchos::RCP<CrsOperatorClass> Op = Teuchos::rcp(new CrsOperatorClass(A));
    Op->SetFixedBlockSize(blksize_); // store block size information in Operator TODO: implement View mechanism as in MueMat

    currentLevel.Set("A", Teuchos::rcp_dynamic_cast<OOperator>(Op), this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetFixedBlockSize(LocalOrdinal blksize) {
    blksize_ = blksize;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  LocalOrdinal SubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetFixedBlockSize() const {
    return blksize_;
  }

} // namespace MueLu

#endif /* MUELU_SUBBLOCKAFACTORY_DEF_HPP_ */
