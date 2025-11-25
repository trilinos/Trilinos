// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_UTILITIES_DEF_HPP
#define MUELU_UTILITIES_DEF_HPP

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <iostream>

#include "MueLu_ConfigDefs.hpp"
#include "Xpetra_TpetraRowMatrix.hpp"

#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <TpetraExt_MatrixMatrix.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#include <Xpetra_TpetraOperator.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_TpetraBlockCrsMatrix.hpp>

#include <Xpetra_BlockedCrsMatrix.hpp>
//#include <Xpetra_DefaultPlatform.hpp>
#include <Xpetra_IO.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Operator.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#include <Xpetra_MatrixMatrix.hpp>

#include <MueLu_Utilities_decl.hpp>

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Transpose(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, bool /* optimizeTranspose */, const std::string& label, const Teuchos::RCP<Teuchos::ParameterList>& params) {
  std::string TorE = "tpetra";

  if (TorE == "tpetra") {
    using Helpers = Xpetra::Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    /***************************************************************/
    if (Helpers::isTpetraCrs(Op)) {
      const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tpetraOp = toTpetra(Op);

      RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> A;
      Tpetra::RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposer(rcpFromRef(tpetraOp), label);  // more than meets the eye

      {
        using Teuchos::ParameterList;
        using Teuchos::rcp;
        RCP<ParameterList> transposeParams = params.is_null() ? rcp(new ParameterList) : rcp(new ParameterList(*params));
        transposeParams->set("sort", false);
        A = transposer.createTranspose(transposeParams);
      }

      RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> AA = rcp(new Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(A));
      RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> AAA      = rcp_implicit_cast<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(AA);
      RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> AAAA        = rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(AAA));
      if (!AAAA->isFillComplete())
        AAAA->fillComplete(Op.getRangeMap(), Op.getDomainMap());

      if (Op.IsView("stridedMaps"))
        AAAA->CreateView("stridedMaps", Teuchos::rcpFromRef(Op), true /*doTranspose*/);

      return AAAA;
    } else if (Helpers::isTpetraBlockCrs(Op)) {
      using XMatrix        = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
      using XCrsMatrix     = Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
      using XCrsMatrixWrap = Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
      using BCRS           = Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
      // using CRS  = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
      const BCRS& tpetraOp = toTpetraBlock(Op);

      RCP<BCRS> At;
      {
        Tpetra::BlockCrsMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposer(rcpFromRef(tpetraOp), label);

        using Teuchos::ParameterList;
        using Teuchos::rcp;
        RCP<ParameterList> transposeParams = params.is_null() ? rcp(new ParameterList) : rcp(new ParameterList(*params));
        transposeParams->set("sort", false);
        At = transposer.createTranspose(transposeParams);
      }

      RCP<Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> AA = rcp(new Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(At));
      RCP<XCrsMatrix> AAA                                                             = rcp_implicit_cast<XCrsMatrix>(AA);
      RCP<XMatrix> AAAA                                                               = rcp(new XCrsMatrixWrap(AAA));

      if (Op.IsView("stridedMaps"))
        AAAA->CreateView("stridedMaps", Teuchos::rcpFromRef(Op), true /*doTranspose*/);

      return AAAA;
    } else {
      throw Exceptions::RuntimeError("Utilities::Transpose failed, perhaps because matrix is not a Crs matrix");
    }
  }  // if

  // Epetra case
  std::cout << "Utilities::Transpose() not implemented for Epetra" << std::endl;
  return Teuchos::null;

}  // Transpose

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    RealValuedToScalarMultiVector(RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node>> X) {
  RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Xscalar;
#if (defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE) || defined(HAVE_TPETRA_INST_COMPLEX_FLOAT))
  using range_type = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;

  // Need to cast the real-valued multivector to Scalar=complex
  if ((typeid(Scalar).name() == typeid(std::complex<double>).name()) ||
      (typeid(Scalar).name() == typeid(std::complex<float>).name())) {
    size_t numVecs  = X->getNumVectors();
    Xscalar         = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(X->getMap(), numVecs);
    auto XVec       = X->getLocalViewDevice(Tpetra::Access::ReadOnly);
    auto XVecScalar = Xscalar->getLocalViewDevice(Tpetra::Access::ReadWrite);

    Kokkos::parallel_for(
        "MueLu:Utils::RealValuedToScalarMultiVector", range_type(0, X->getLocalLength()),
        KOKKOS_LAMBDA(const size_t i) {
          for (size_t j = 0; j < numVecs; j++)
            XVecScalar(i, j) = XVec(i, j);
        });
  } else
#endif
    Xscalar = rcp_dynamic_cast<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(X);
  return Xscalar;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node>>
Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractCoordinatesFromParameterList(ParameterList& paramList) {
  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node>> coordinates = Teuchos::null;

  // check whether coordinates are contained in parameter list
  if (paramList.isParameter("Coordinates") == false)
    return coordinates;

    // define Tpetra::MultiVector type with Scalar=float only if
    // * ETI is turned off, since then the compiler will instantiate it automatically OR
    // * Tpetra is instantiated on Scalar=float
#if !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION) || defined(HAVE_TPETRA_INST_FLOAT)
  typedef Tpetra::MultiVector<float, LocalOrdinal, GlobalOrdinal, Node> tfMV;
  RCP<tfMV> floatCoords = Teuchos::null;
#endif

  // define Tpetra::MultiVector type with Scalar=double only if
  // * ETI is turned off, since then the compiler will instantiate it automatically OR
  // * Tpetra is instantiated on Scalar=double
#if !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION) || defined(HAVE_TPETRA_INST_DOUBLE)
  typedef Tpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node> tdMV;
  RCP<tdMV> doubleCoords = Teuchos::null;
  if (paramList.isType<RCP<tdMV>>("Coordinates")) {
    // Coordinates are stored as a double vector
    doubleCoords = paramList.get<RCP<tdMV>>("Coordinates");
    paramList.remove("Coordinates");
  }
#if !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION) || defined(HAVE_TPETRA_INST_FLOAT)
  else if (paramList.isType<RCP<tfMV>>("Coordinates")) {
    // check if coordinates are stored as a float vector
    floatCoords = paramList.get<RCP<tfMV>>("Coordinates");
    paramList.remove("Coordinates");
    doubleCoords = rcp(new tdMV(floatCoords->getMap(), floatCoords->getNumVectors()));
    deep_copy(*doubleCoords, *floatCoords);
  }
#endif
  // We have the coordinates in a Tpetra double vector
  if (doubleCoords != Teuchos::null) {
    // rcp(new Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(Vtpetra));
    coordinates = Xpetra::toXpetra(doubleCoords);
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(coordinates));
    TEUCHOS_TEST_FOR_EXCEPT(doubleCoords->getNumVectors() != coordinates->getNumVectors());
  }
#else
  // coordinates usually are stored as double vector
  // Tpetra is not instantiated on scalar=double
  throw Exceptions::RuntimeError("ExtractCoordinatesFromParameterList: The coordinates vector in parameter list is expected to be a Tpetra multivector with SC=double or float.");
#endif

  // check for Xpetra coordinates vector
  if (paramList.isType<decltype(coordinates)>("Coordinates")) {
    coordinates = paramList.get<decltype(coordinates)>("Coordinates");
  }

  return coordinates;
}  // ExtractCoordinatesFromParameterList

}  // namespace MueLu

#define MUELU_UTILITIES_SHORT
#endif  // MUELU_UTILITIES_DEF_HPP

//  LocalWords:  LocalOrdinal
