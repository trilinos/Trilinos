// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_CREATE_TPETRA_PRECONDITIONER_HPP
#define MUELU_CREATE_TPETRA_PRECONDITIONER_HPP

//! @file
//! @brief Various adapters that will create a MueLu preconditioner that is a Tpetra::Operator.

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Tpetra_Operator.hpp>
#include <Tpetra_RowMatrix.hpp>
#include <Xpetra_TpetraBlockCrsMatrix.hpp>
#include <Tpetra_BlockCrsMatrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include <MueLu.hpp>

#include <MueLu_Exceptions.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_MasterList.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateXpetraPreconditioner.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_HierarchyUtils.hpp>

#if defined(HAVE_MUELU_AMGX)
#include <MueLu_AMGXOperator.hpp>
#include <amgx_c.h>
#include "cuda_runtime.h"
#endif

namespace MueLu {

/*!
  @brief Helper function to create a MueLu or AMGX preconditioner that can be used by Tpetra.
  @ingroup MueLuAdapters
  Given a Tpetra::Operator, this function returns a constructed MueLu preconditioner.
  @param[in] inA Matrix
  @param[in] inParamList Parameter list
*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<MueLu::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
CreateTpetraPreconditioner(const Teuchos::RCP<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& inA,
                           Teuchos::ParameterList& inParamList) {
#include "Xpetra_UseShortNames.hpp"

  using Teuchos::ParameterList;

  using tpMultiVector         = Tpetra::MultiVector<SC, LO, GO, NO>;
  using coordMultiVector      = Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LO, GO, NO>;
  using tpCoordMultiVector    = Tpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LO, GO, NO>;
  using tpLocalOrdinalVector  = Tpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;
  using Hierarchy             = Hierarchy<SC, LO, GO, NO>;
  using crs_matrix_type       = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using block_crs_matrix_type = Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

#if defined(HAVE_MUELU_AMGX)
  std::string externalMG = "use external multigrid package";
  if (inParamList.isParameter(externalMG) && inParamList.get<std::string>(externalMG) == "amgx") {
    const RCP<crs_matrix_type> constCrsA = rcp_dynamic_cast<crs_matrix_type>(inA);
    TEUCHOS_TEST_FOR_EXCEPTION(constCrsA == Teuchos::null, Exceptions::RuntimeError, "CreateTpetraPreconditioner: failed to dynamic cast to Tpetra::CrsMatrix, which is required to be able to use AmgX.");
    return rcp(new AMGXOperator<SC, LO, GO, NO>(constCrsA, inParamList));
  }
#endif

  // Wrap A
  RCP<Matrix> A;
  RCP<block_crs_matrix_type> bcrsA = rcp_dynamic_cast<block_crs_matrix_type>(inA);
  RCP<crs_matrix_type> crsA        = rcp_dynamic_cast<crs_matrix_type>(inA);
  if (crsA != Teuchos::null)
    A = Xpetra::toXpetra(crsA);
  else if (bcrsA != Teuchos::null) {
    RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> temp = rcp(new Xpetra::TpetraBlockCrsMatrix<SC, LO, GO, NO>(bcrsA));
    TEUCHOS_TEST_FOR_EXCEPTION(temp == Teuchos::null, Exceptions::RuntimeError, "CreateTpetraPreconditioner: cast from Tpetra::BlockCrsMatrix to Xpetra::TpetraBlockCrsMatrix failed.");
    A = rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(temp));
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "CreateTpetraPreconditioner: only Tpetra CrsMatrix and BlockCrsMatrix types are supported.");
  }

  Teuchos::ParameterList& userList = inParamList.sublist("user data");
  if (userList.isParameter("Coordinates")) {
    RCP<coordMultiVector> coordinates = Teuchos::null;
    try {
      coordinates = Xpetra::toXpetra(userList.get<RCP<tpCoordMultiVector>>("Coordinates"));
    } catch (Teuchos::Exceptions::InvalidParameterType&) {
      coordinates = userList.get<RCP<coordMultiVector>>("Coordinates");
    }
    userList.set<RCP<coordMultiVector>>("Coordinates", coordinates);
  }

  if (userList.isParameter("Material")) {
    RCP<MultiVector> material = Teuchos::null;
    try {
      material = Xpetra::toXpetra(userList.get<RCP<tpMultiVector>>("Material"));
    } catch (Teuchos::Exceptions::InvalidParameterType&) {
      material = userList.get<RCP<MultiVector>>("Material");
    }
    userList.set<RCP<MultiVector>>("Material", material);
  }

  if (userList.isParameter("Nullspace")) {
    RCP<MultiVector> nullspace = Teuchos::null;
    try {
      nullspace = Xpetra::toXpetra(userList.get<RCP<tpMultiVector>>("Nullspace"));
    } catch (Teuchos::Exceptions::InvalidParameterType&) {
      nullspace = userList.get<RCP<MultiVector>>("Nullspace");
    }
    userList.set<RCP<MultiVector>>("Nullspace", nullspace);
  }

  if (userList.isParameter("BlockNumber")) {
    RCP<LocalOrdinalVector> blockNumber = Teuchos::null;
    try {
      blockNumber = Xpetra::toXpetra(userList.get<RCP<tpLocalOrdinalVector>>("BlockNumber"));
    } catch (Teuchos::Exceptions::InvalidParameterType&) {
      blockNumber = userList.get<RCP<LocalOrdinalVector>>("BlockNumber");
    }
    userList.set<RCP<LocalOrdinalVector>>("BlockNumber", blockNumber);
  }

  RCP<Hierarchy> H = MueLu::CreateXpetraPreconditioner<SC, LO, GO, NO>(A, inParamList);
  return rcp(new MueLu::TpetraOperator<SC, LO, GO, NO>(H));
}

/*!
  @brief Helper function to create a MueLu preconditioner that can be used by Tpetra.
  @ingroup MueLuAdapters

  Given a Tpetra::Operator, this function returns a constructed MueLu preconditioner.

  @param[in] inA Matrix
  @param[in] xmlFileName XML file containing MueLu options
*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<MueLu::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
CreateTpetraPreconditioner(const Teuchos::RCP<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& inA,
                           const std::string& xmlFileName) {
  Teuchos::ParameterList paramList;
  Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), *inA->getDomainMap()->getComm());
  return CreateTpetraPreconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Node>(inA, paramList);
}

/*!
  @brief Helper function to create a MueLu preconditioner that can be used by Tpetra.
  @ingroup MueLuAdapters

  Given a Tpetra::Operator, this function returns a constructed MueLu preconditioner.

  @param[in] inA Matrix
*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<MueLu::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
CreateTpetraPreconditioner(const Teuchos::RCP<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& inA) {
  Teuchos::ParameterList paramList;
  return CreateTpetraPreconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Node>(inA, paramList);
}

/*!
  @brief Helper function to reuse an existing MueLu preconditioner.
  @ingroup MueLuAdapters

  @param[in] inA Matrix
  @param[in] Op  Existing MueLu preconditioner.
*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ReuseTpetraPreconditioner(const Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& inA,
                               MueLu::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op) {
  typedef Scalar SC;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Node NO;

  typedef Xpetra::Matrix<SC, LO, GO, NO> Matrix;
  typedef MueLu ::Hierarchy<SC, LO, GO, NO> Hierarchy;

  RCP<Hierarchy> H = Op.GetHierarchy();
  RCP<Matrix> A    = Xpetra::toXpetra(inA);

  MueLu::ReuseXpetraPreconditioner<SC, LO, GO, NO>(A, H);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ReuseTpetraPreconditioner(const Teuchos::RCP<Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& inA,
                               MueLu::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op) {
  typedef Scalar SC;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Node NO;

  typedef Xpetra::Matrix<SC, LO, GO, NO> Matrix;
  typedef MueLu ::Hierarchy<SC, LO, GO, NO> Hierarchy;

  RCP<Hierarchy> H                            = Op.GetHierarchy();
  RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> temp = rcp(new Xpetra::TpetraBlockCrsMatrix<SC, LO, GO, NO>(inA));
  TEUCHOS_TEST_FOR_EXCEPTION(temp == Teuchos::null, Exceptions::RuntimeError, "ReuseTpetraPreconditioner: cast from Tpetra::BlockCrsMatrix to Xpetra::TpetraBlockCrsMatrix failed.");
  RCP<Matrix> A = rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(temp));

  MueLu::ReuseXpetraPreconditioner<SC, LO, GO, NO>(A, H);
}

}  // namespace MueLu

#endif  // ifndef MUELU_CREATE_TPETRA_PRECONDITIONER_HPP
