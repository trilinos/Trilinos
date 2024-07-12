// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_PERFUTILS_DECL_HPP
#define MUELU_PERFUTILS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_Export_fwd.hpp>
#include <Xpetra_Import_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_PerfUtils_fwd.hpp"

namespace MueLu {
// MPI helpers
#define MueLu_sumAll(rcpComm, in, out) \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out))
#define MueLu_minAll(rcpComm, in, out) \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MIN, in, Teuchos::outArg(out))
#define MueLu_maxAll(rcpComm, in, out) \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MAX, in, Teuchos::outArg(out))

template <class Scalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class PerfUtils {
#undef MUELU_PERFUTILS_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  static std::string PrintMatrixInfo(const Matrix& A, const std::string& msgTag, RCP<const Teuchos::ParameterList> params = Teuchos::null);

  static std::string PrintImporterInfo(RCP<const Import> importer, const std::string& msgTag);

  static std::string CommPattern(const Matrix& A, const std::string& msgTag, RCP<const Teuchos::ParameterList> params = Teuchos::null);

 private:
  static bool CheckMatrix(const Matrix& A);
};

}  // namespace MueLu

#define MUELU_PERFUTILS_SHORT
#endif  // MUELU_PERFUTILS_DECL_HPP
