// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
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
                                         Kokkos::View<bool*, typename Node::device_type>& BCrowsKokkos,
                                         Kokkos::View<bool*, typename Node::device_type>& BCcolsKokkos,
                                         Kokkos::View<bool*, typename Node::device_type>& BCdomainKokkos,
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
                                         Kokkos::View<bool*, typename Node::device_type>& BCrowsKokkos,
                                         Kokkos::View<bool*, typename Node::device_type>& BCcolsKokkos,
                                         Kokkos::View<bool*, typename Node::device_type>& BCdomainKokkos,
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
  PtAPWrapper(const RCP<Matrix>& A, const RCP<Matrix>& P, Teuchos::ParameterList& params, std::string& label);
};

}  // namespace MueLu

#define MUELU_MAXWELL_UTILS_SHORT
#endif  // MUELU_MAXWELL_UTILS_DECL_HPP
