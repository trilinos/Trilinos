// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_MAXWELL_UTILS_DECL_HPP
#define MUELU_MAXWELL_UTILS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

#include "Xpetra_Map_fwd.hpp"
#include "Xpetra_Matrix_fwd.hpp"
#include "Xpetra_VectorFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_ThresholdAFilterFactory_fwd.hpp"
#include "MueLu_RAPFactory_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {

/*!
  @brief Utility functions for Maxwell

  @ingroup MueLuAdapters
*/
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class Maxwell_Utils : public VerboseObject {
#undef MUELU_MAXWELL_UTILS_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;
  //    typedef typename Teuchos::ScalarTraits<Scalar>::coordinateType coordinateType;
  //    typedef typename Xpetra::MultiVector<coordinateType,LO,GO,NO> RealValuedMultiVector;

  //! Detect Dirichlet boundary conditions
  static void detectBoundaryConditionsSM(RCP<Matrix>& SM_Matrix,
                                         RCP<Matrix>& D0_Matrix,
                                         magnitudeType rowSumTol,
                                         bool useKokkos_,
                                         Kokkos::View<bool*, typename Node::device_type::memory_space>& BCrowsKokkos,
                                         Kokkos::View<bool*, typename Node::device_type::memory_space>& BCcolsKokkos,
                                         Kokkos::View<bool*, typename Node::device_type::memory_space>& BCdomainKokkos,
                                         int& BCedges,
                                         int& BCnodes,
                                         Teuchos::ArrayRCP<bool>& BCrows,
                                         Teuchos::ArrayRCP<bool>& BCcols,
                                         Teuchos::ArrayRCP<bool>& BCdomain,
                                         bool& allEdgesBoundary,
                                         bool& allNodesBoundary);

  //! Detect Dirichlet boundary conditions
  static void detectBoundaryConditionsSM(RCP<Matrix>& SM_Matrix,
                                         RCP<Matrix>& D0_Matrix,
                                         magnitudeType rowSumTol,
                                         Kokkos::View<bool*, typename Node::device_type::memory_space>& BCrowsKokkos,
                                         Kokkos::View<bool*, typename Node::device_type::memory_space>& BCcolsKokkos,
                                         Kokkos::View<bool*, typename Node::device_type::memory_space>& BCdomainKokkos,
                                         int& BCedges,
                                         int& BCnodes,
                                         bool& allEdgesBoundary,
                                         bool& allNodesBoundary);

  //! Remove explicit zeros
  static RCP<Matrix> removeExplicitZeros(const RCP<Matrix>& A,
                                         const magnitudeType tolerance,
                                         const bool keepDiagonal        = true,
                                         const size_t expectedNNZperRow = 0);

  static void removeExplicitZeros(Teuchos::ParameterList& parameterList,
                                  RCP<Matrix>& D0_Matrix,
                                  RCP<Matrix>& SM_Matrix,
                                  RCP<Matrix>& M1_Matrix,
                                  RCP<Matrix>& Ms_Matrix);

  static void removeExplicitZeros(Teuchos::ParameterList& parameterList,
                                  RCP<Matrix>& D0_Matrix,
                                  RCP<Matrix>& SM_Matrix) {
    RCP<Matrix> dummy;
    removeExplicitZeros(parameterList, D0_Matrix, SM_Matrix, dummy, dummy);
  }

  static void thresholdedAbs(const RCP<Matrix>& A,
                             const magnitudeType thresholded);

  //! Sets matvec params on a matrix
  static void setMatvecParams(Matrix& A, RCP<ParameterList> matvecParams);

  //! Performs an P^T AP
  static RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  PtAPWrapper(const RCP<Matrix>& A, const RCP<Matrix>& P, Teuchos::ParameterList& params, const std::string& label);
};

}  // namespace MueLu

#define MUELU_MAXWELL_UTILS_SHORT
#endif  // MUELU_MAXWELL_UTILS_DECL_HPP
