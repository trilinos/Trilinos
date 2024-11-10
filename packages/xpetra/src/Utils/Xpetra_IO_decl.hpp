// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_XPETRA_SUP_UTILS_XPETRA_IO_HPP_
#define PACKAGES_XPETRA_SUP_UTILS_XPETRA_IO_HPP_

#include <fstream>
#include "Xpetra_ConfigDefs.hpp"

#ifdef HAVE_XPETRA_EPETRA
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#endif

#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_MultiVectorIn.h>
#include <EpetraExt_BlockMapIn.h>
#include <Xpetra_EpetraUtils.hpp>
#include <Xpetra_EpetraMultiVector.hpp>
#include <EpetraExt_BlockMapOut.h>
#endif

#ifdef HAVE_XPETRA_TPETRA
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <TpetraExt_MatrixMatrix.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#include <Xpetra_TpetraCrsGraph.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_TpetraBlockCrsMatrix.hpp>
#include "Tpetra_Util.hpp"
#endif

#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraMap.hpp>
#endif

#include "Xpetra_Matrix.hpp"
#include "Xpetra_MatrixMatrix.hpp"
#include "Xpetra_Helpers.hpp"
#include "Xpetra_CrsGraph.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_BlockedCrsMatrix.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_StridedMap.hpp"
#include "Xpetra_StridedMapFactory.hpp"
#include "Xpetra_MapExtractor.hpp"
#include "Xpetra_MatrixFactory.hpp"

#include <Teuchos_TestForException.hpp>
#include <Teuchos_MatrixMarket_Raw_Writer.hpp>
#include <Teuchos_MatrixMarket_Raw_Reader.hpp>
#include <string>

namespace Xpetra {

#ifdef HAVE_XPETRA_EPETRA
// This non-member templated function exists so that the matrix-matrix multiply will compile if Epetra, Tpetra, and ML are enabled.
template <class SC, class LO, class GO, class NO>
RCP<Xpetra::CrsMatrixWrap<SC, LO, GO, NO>>
Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap(RCP<Epetra_CrsMatrix>& /* epAB */) {
  TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError,
                             "Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap cannot be used with Scalar != double, LocalOrdinal != int, GlobalOrdinal != int");
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

// specialization for the case of ScalarType=double and LocalOrdinal=GlobalOrdinal=int
template <>
inline RCP<Xpetra::CrsMatrixWrap<double, int, int, Xpetra::EpetraNode>> Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<double, int, int, Xpetra::EpetraNode>(RCP<Epetra_CrsMatrix>& epAB) {
  typedef double SC;
  typedef int LO;
  typedef int GO;
  typedef Xpetra::EpetraNode NO;

  RCP<Xpetra::EpetraCrsMatrixT<GO, NO>> tmpC1      = rcp(new Xpetra::EpetraCrsMatrixT<GO, NO>(epAB));
  RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> tmpC2     = Teuchos::rcp_implicit_cast<Xpetra::CrsMatrix<SC, LO, GO, NO>>(tmpC1);
  RCP<Xpetra::CrsMatrixWrap<SC, LO, GO, NO>> tmpC3 = rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(tmpC2));

  return tmpC3;
}

template <class SC, class LO, class GO, class NO>
RCP<Xpetra::MultiVector<SC, LO, GO, NO>>
Convert_Epetra_MultiVector_ToXpetra_MultiVector(RCP<Epetra_MultiVector>& epX) {
  TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError,
                             "Convert_Epetra_MultiVector_ToXpetra_MultiVector cannot be used with Scalar != double, LocalOrdinal != int, GlobalOrdinal != int");
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

// specialization for the case of ScalarType=double and LocalOrdinal=GlobalOrdinal=int
template <>
inline RCP<Xpetra::MultiVector<double, int, int, Xpetra::EpetraNode>> Convert_Epetra_MultiVector_ToXpetra_MultiVector<double, int, int, Xpetra::EpetraNode>(RCP<Epetra_MultiVector>& epX) {
  typedef double SC;
  typedef int LO;
  typedef int GO;
  typedef Xpetra::EpetraNode NO;

  RCP<Xpetra::MultiVector<SC, LO, GO, NO>> tmp = Xpetra::toXpetra<GO, NO>(epX);
  return tmp;
}

#endif

/*!
  @class IO
  @brief Xpetra utility class containing IO routines to read/write vectors, matrices etc...
  */
template <class Scalar,
          class LocalOrdinal  = int,
          class GlobalOrdinal = LocalOrdinal,
          class Node          = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class IO {
 private:
#undef XPETRA_IO_SHORT
#include "Xpetra_UseShortNames.hpp"

 public:
#ifdef HAVE_XPETRA_EPETRA
  //! Helper utility to pull out the underlying Epetra objects from an Xpetra object
  // @{
  /*static RCP<const Epetra_MultiVector>                    MV2EpetraMV(RCP<MultiVector> const vec);
    static RCP<      Epetra_MultiVector>                    MV2NonConstEpetraMV(RCP<MultiVector> vec);

    static const Epetra_MultiVector&                        MV2EpetraMV(const MultiVector& vec);
    static       Epetra_MultiVector&                        MV2NonConstEpetraMV(MultiVector& vec);

    static RCP<const Epetra_CrsMatrix>                      Op2EpetraCrs(RCP<const Matrix> Op);
    static RCP<      Epetra_CrsMatrix>                      Op2NonConstEpetraCrs(RCP<Matrix> Op);

    static const Epetra_CrsMatrix&                          Op2EpetraCrs(const Matrix& Op);
    static       Epetra_CrsMatrix&                          Op2NonConstEpetraCrs(Matrix& Op);*/

  static const Epetra_Map& Map2EpetraMap(const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& map);
  // @}
#endif

#ifdef HAVE_XPETRA_TPETRA
  //! Helper utility to pull out the underlying Tpetra objects from an Xpetra object
  // @{
  /*static RCP<const Tpetra::MultiVector<SC,LO,GO,NO> >     MV2TpetraMV(RCP<MultiVector> const vec);
    static RCP<      Tpetra::MultiVector<SC,LO,GO,NO> >     MV2NonConstTpetraMV(RCP<MultiVector> vec);
    static RCP<      Tpetra::MultiVector<SC,LO,GO,NO> >     MV2NonConstTpetraMV2(MultiVector& vec);

    static const Tpetra::MultiVector<SC,LO,GO,NO>&          MV2TpetraMV(const MultiVector& vec);
    static       Tpetra::MultiVector<SC,LO,GO,NO>&          MV2NonConstTpetraMV(MultiVector& vec);

    static RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO> >       Op2TpetraCrs(RCP<const Matrix> Op);
    static RCP<      Tpetra::CrsMatrix<SC,LO,GO,NO> >       Op2NonConstTpetraCrs(RCP<Matrix> Op);

    static const Tpetra::CrsMatrix<SC,LO,GO,NO>&            Op2TpetraCrs(const Matrix& Op);
    static       Tpetra::CrsMatrix<SC,LO,GO,NO>&            Op2NonConstTpetraCrs(Matrix& Op);

    static RCP<const Tpetra::RowMatrix<SC,LO,GO,NO> >       Op2TpetraRow(RCP<const Matrix> Op);
    static RCP<      Tpetra::RowMatrix<SC,LO,GO,NO> >       Op2NonConstTpetraRow(RCP<Matrix> Op);*/

  static const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> Map2TpetraMap(const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& map);
#endif

  //! Read/Write methods
  //@{
  /*! @brief Save map to file. */
  static void Write(const std::string& fileName, const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& M);

  /*! @brief Save multivector to file in Matrix Market format.  */
  static void Write(const std::string& fileName, const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& vec);

  // CAG:
  // The class is templated on the usual SC-LO-GO-NO.
  // Instead of instantiating the entire class using Scalar=LO or GO and then dealing with the headaches that integer-valued CrsMatrix creates, we use these two methods.

  /*! @brief Save multivector with Scalar=LocalOrdinal to file in Matrix Market format.  */
  static void WriteLOMV(const std::string& fileName, const Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>& vec);

  /*! @brief Save multivector with Scalar=GlobalOrdinal to file in Matrix Market format.  */
  static void WriteGOMV(const std::string& fileName, const Xpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>& vec);

  /*! @brief Save matrix to file in Matrix Market format. */
  static void Write(const std::string& fileName, const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, const bool& writeAllMaps = false);

  /*! @brief Save local parts of matrix to files in Matrix Market format. */
  static void WriteLocal(const std::string& fileName, const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op);

  /*! @brief Save block matrix to one file per block in Matrix Market format.

  We write one MatrixMarket file for each block of the given BlockedCrsMatrix.
  The block (row,col) indicators are included into the file name, such that filename02.m indicates
  the block in row = 0 and col = 2, for example.

  While the user can enable/disable the output of all maps of each matrix block,
  we always write the range and domain maps of each block as well of the full blocked operator
  in order to enable reconstruction of the MapExtractor objects for range and domain maps.

  @param fileName String to indicate file name of files to be written
  @param Op BlockedCrsMatrix to be written
  @param writeAllMaps Flag to control output of maps to separate files (defaults to \c false )
  */
  static void WriteBlockedCrsMatrix(const std::string& fileName, const Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, const bool& writeAllMaps = false);

  //! @brief Read matrix from file in Matrix Market or binary format.
  static Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Read(const std::string& fileName, Xpetra::UnderlyingLib lib, const RCP<const Teuchos::Comm<int>>& comm, bool binary = false);

  /*! @brief Read matrix from file in Matrix Market or binary format.

    If only rowMap is specified, then it is used for the domainMap and rangeMap, as well.
    */
  static Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  Read(const std::string& filename,
       const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowMap,
       RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> colMap          = Teuchos::null,
       const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> domainMap = Teuchos::null,
       const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rangeMap  = Teuchos::null,
       const bool callFillComplete                                               = true,
       const bool binary                                                         = false,
       const bool tolerant                                                       = false,
       const bool debug                                                          = false);

  /*! @brief Read matrix from local files in Matrix Market or binary format.

    The file name format is filename.SIZE.RANK, where SIZE is the
    size of the communicator of the rowMap and RANK is the MPI ranks
    of the calling process.

    If only rowMap is specified, then it is used for the domainMap and rangeMap, as well.
    */
  static Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> ReadLocal(const std::string& filename,
                                                                                           const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowMap,
                                                                                           RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> colMap,
                                                                                           const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> domainMap = Teuchos::null,
                                                                                           const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rangeMap  = Teuchos::null,
                                                                                           const bool callFillComplete                                               = true,
                                                                                           const bool binary                                                         = false,
                                                                                           const bool tolerant                                                       = false,
                                                                                           const bool debug                                                          = false);
  //@}

  //! @brief Read a MultiVector from file in Matrix Matrix or binary format.
  static RCP<MultiVector> ReadMultiVector(const std::string& fileName,
                                          const RCP<const Map>& map,
                                          const bool binary = false);

  //! @brief Read a MultiVector with Scalar=LocalOrdinal from file in Matrix Matrix or binary format.
  static RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> ReadMultiVectorLO(const std::string& fileName,
                                                                                                     const RCP<const Map>& map,
                                                                                                     const bool binary = false);

  static RCP<const Map> ReadMap(const std::string& fileName,
                                Xpetra::UnderlyingLib lib,
                                const RCP<const Teuchos::Comm<int>>& comm,
                                const bool binary = false);

  /*! @brief Read block matrix from one file per block in Matrix Market format.

  We read one MatrixMarket file for each block of the desired BlockedCrsMatrix.
  The block (row,col) indicators are included into the file name, such that filename02.m indicates
  the block in row = 0 and col = 2, for example.

  We also read the range and domain maps of each block as well of the full blocked operator
  in order to enable reconstruction of the MapExtractor objects for range and domain maps.

  @param fileName String to indicate file name of files to be written
  @param lib Underlying type of linear algebra package
  @param comm Communicator
  */
  static RCP<const Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> ReadBlockedCrsMatrix(const std::string& fileName, Xpetra::UnderlyingLib lib, const RCP<const Teuchos::Comm<int>>& comm);

  //! Little helper function to convert non-string types to strings
  template <class T>
  static std::string toString(const T& what);
};

#ifdef HAVE_XPETRA_EPETRA
/*!
  @class IO
  @brief Xpetra utility class containing IO routines to read/write vectors, matrices.

  Specialization for LO=GO=int

  TODO: do we need specialization for SC=double and std::complex<>???
  TODO: need specialization for <int, long long>
*/
template <class Scalar>
class IO<Scalar, int, int, EpetraNode> {
 public:
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef EpetraNode Node;

#ifdef HAVE_XPETRA_EPETRA
  //! Helper utility to pull out the underlying Epetra objects from an Xpetra object
  // @{
  static const Epetra_Map& Map2EpetraMap(const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& map) {
    RCP<const Xpetra::EpetraMapT<GlobalOrdinal, Node>> xeMap = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraMapT<GlobalOrdinal, Node>>(Teuchos::rcpFromRef(map));
    if (xeMap == Teuchos::null)
      throw Exceptions::BadCast("IO::Map2EpetraMap : Cast from Xpetra::Map to Xpetra::EpetraMap failed");
    return xeMap->getEpetra_Map();
  }
  // @}
#endif

#ifdef HAVE_XPETRA_TPETRA
  //! Helper utility to pull out the underlying Tpetra objects from an Xpetra object
  // @{
  static const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> Map2TpetraMap(const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& map) {
    const RCP<const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node>>& tmp_TMap = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node>>(rcpFromRef(map));
    if (tmp_TMap == Teuchos::null)
      throw Exceptions::BadCast("IO::Map2TpetraMap : Cast from Xpetra::Map to Xpetra::TpetraMap failed");
    return tmp_TMap->getTpetra_Map();
  }
#endif

  //! Read/Write methods
  //@{
  /*! @brief Save map to file. */
  static void Write(const std::string& fileName, const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& M) {
    RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> tmp_Map = rcpFromRef(M);
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
    const RCP<const Xpetra::EpetraMapT<GlobalOrdinal, Node>>& tmp_EMap = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraMapT<GlobalOrdinal, Node>>(tmp_Map);
    if (tmp_EMap != Teuchos::null) {
      int rv = EpetraExt::BlockMapToMatrixMarketFile(fileName.c_str(), tmp_EMap->getEpetra_Map());
      if (rv != 0)
        throw Exceptions::RuntimeError("EpetraExt::BlockMapToMatrixMarketFile() return value of " + Teuchos::toString(rv));
      return;
    }
#endif  // HAVE_XPETRA_EPETRA

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    // do nothing
#else
    const RCP<const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node>>& tmp_TMap =
        Teuchos::rcp_dynamic_cast<const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node>>(tmp_Map);
    if (tmp_TMap != Teuchos::null) {
      RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> TMap = tmp_TMap->getTpetra_Map();
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>::writeMapFile(fileName, *TMap);
      return;
    }
#endif
#endif  // HAVE_XPETRA_TPETRA
    throw Exceptions::BadCast("Could not cast to EpetraMap or TpetraMap in map writing");
  }

  static void Write(const std::string& fileName, const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& vec) {
    std::string mapfile = "map_" + fileName;
    Write(mapfile, *(vec.getMap()));

    RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmp_Vec = Teuchos::rcpFromRef(vec);
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
    const RCP<const Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node>>& tmp_EVec = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node>>(tmp_Vec);
    if (tmp_EVec != Teuchos::null) {
      int rv = EpetraExt::MultiVectorToMatrixMarketFile(fileName.c_str(), *(tmp_EVec->getEpetra_MultiVector()));
      if (rv != 0)
        throw Exceptions::RuntimeError("EpetraExt::RowMatrixToMatrixMarketFile return value of " + Teuchos::toString(rv));
      return;
    }
#endif  // HAVE_XPETRA_EPETRAEXT

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    // do nothin
#else
    const RCP<const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& tmp_TVec =
        Teuchos::rcp_dynamic_cast<const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tmp_Vec);
    if (tmp_TVec != Teuchos::null) {
      RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TVec = tmp_TVec->getTpetra_MultiVector();
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>::writeDenseFile(fileName, TVec);
      return;
    }
#endif
#endif  // HAVE_XPETRA_TPETRA

    throw Exceptions::BadCast("Could not cast to EpetraMultiVector or TpetraMultiVector in multivector writing");
  }

  static void WriteLOMV(const std::string& fileName, const Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>& vec) {
    std::string mapfile = "map_" + fileName;
    Write(mapfile, *(vec.getMap()));

    RCP<const Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> tmp_Vec = Teuchos::rcpFromRef(vec);
#ifdef HAVE_XPETRA_TPETRA
    const RCP<const Xpetra::TpetraMultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>& tmp_TVec =
        Teuchos::rcp_dynamic_cast<const Xpetra::TpetraMultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>(tmp_Vec);
    if (tmp_TVec != Teuchos::null) {
      RCP<const Tpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> TVec = tmp_TVec->getTpetra_MultiVector();
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>::writeDenseFile(fileName, TVec);
      return;
    } else
#endif  // HAVE_XPETRA_TPETRA
    {
      throw Exceptions::RuntimeError("Xpetra cannot write MV<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> when the underlying library is Epetra.");
    }

    throw Exceptions::BadCast("Could not cast to EpetraMultiVector or TpetraMultiVector in multivector writing");
  }

  static void WriteGOMV(const std::string& fileName, const Xpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>& vec) {
    std::string mapfile = "map_" + fileName;
    Write(mapfile, *(vec.getMap()));

    RCP<const Xpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> tmp_Vec = Teuchos::rcpFromRef(vec);
#ifdef HAVE_XPETRA_TPETRA
    const RCP<const Xpetra::TpetraMultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>& tmp_TVec =
        Teuchos::rcp_dynamic_cast<const Xpetra::TpetraMultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>(tmp_Vec);
    if (tmp_TVec != Teuchos::null) {
      RCP<const Tpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> TVec = tmp_TVec->getTpetra_MultiVector();
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>::writeDenseFile(fileName, TVec);
      return;
    } else
#endif  // HAVE_XPETRA_TPETRA
    {
      throw Exceptions::RuntimeError("Xpetra cannot write MV<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> when the underlying library is Epetra.");
    }

    throw Exceptions::BadCast("Could not cast to EpetraMultiVector or TpetraMultiVector in multivector writing");
  }

  static void Write(const std::string& fileName, const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, const bool& writeAllMaps = false) {
    Write("rowmap_" + fileName, *(Op.getRowMap()));
    if (!Op.getDomainMap()->isSameAs(*(Op.getRowMap())) || writeAllMaps)
      Write("domainmap_" + fileName, *(Op.getDomainMap()));
    if (!Op.getRangeMap()->isSameAs(*(Op.getRowMap())) || writeAllMaps)
      Write("rangemap_" + fileName, *(Op.getRangeMap()));
    if (!Op.getColMap()->isSameAs(*(Op.getDomainMap())) || writeAllMaps)
      Write("colmap_" + fileName, *(Op.getColMap()));

    const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>& crsOp =
        dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(Op);
    RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmp_CrsMtx = crsOp.getCrsMatrix();
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
    const RCP<const Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>>& tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>>(tmp_CrsMtx);
    if (tmp_ECrsMtx != Teuchos::null) {
      RCP<const Epetra_CrsMatrix> A = tmp_ECrsMtx->getEpetra_CrsMatrix();
      int rv                        = EpetraExt::RowMatrixToMatrixMarketFile(fileName.c_str(), *A);
      if (rv != 0)
        throw Exceptions::RuntimeError("EpetraExt::RowMatrixToMatrixMarketFile return value of " + Teuchos::toString(rv));
      return;
    }
#endif  // endif HAVE_XPETRA_EPETRA

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    // do nothin
#else
    const RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& tmp_TCrsMtx =
        Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tmp_CrsMtx);
    if (tmp_TCrsMtx != Teuchos::null) {
      RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> A = tmp_TCrsMtx->getTpetra_CrsMatrix();
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>::writeSparseFile(fileName, A);
      return;
    }
    const RCP<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& tmp_BlockCrs =
        Teuchos::rcp_dynamic_cast<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tmp_CrsMtx);
    if (tmp_BlockCrs != Teuchos::null) {
      std::ofstream outstream(fileName, std::ofstream::out);
      Teuchos::FancyOStream ofs(Teuchos::rcpFromRef(outstream));
      tmp_BlockCrs->getTpetra_BlockCrsMatrix()->describe(ofs, Teuchos::VERB_EXTREME);
      return;
    }

#endif
#endif  // HAVE_XPETRA_TPETRA

    throw Exceptions::BadCast("Could not cast to EpetraCrsMatrix or TpetraCrsMatrix in matrix writing");
  }

  static void Write(const std::string& fileName, const Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>& graph, const bool& writeAllMaps = false) {
    Write("rowmap_" + fileName, *(graph.getRowMap()));
    if (!graph.getDomainMap()->isSameAs(*(graph.getRowMap())) || writeAllMaps)
      Write("domainmap_" + fileName, *(graph.getDomainMap()));
    if (!graph.getRangeMap()->isSameAs(*(graph.getRowMap())) || writeAllMaps)
      Write("rangemap_" + fileName, *(graph.getRangeMap()));
    if (!graph.getColMap()->isSameAs(*(graph.getDomainMap())) || writeAllMaps)
      Write("colmap_" + fileName, *(graph.getColMap()));

    RCP<const Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> tmp_Graph = rcpFromRef(graph);

#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
    const RCP<const Xpetra::EpetraCrsGraphT<GlobalOrdinal, Node>>& tmp_ECrsGraph = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraCrsGraphT<GlobalOrdinal, Node>>(tmp_Graph);
    if (tmp_ECrsGraph != Teuchos::null) {
      throw Exceptions::BadCast("Writing not implemented for EpetraCrsGraphT");
    }
#endif  // endif HAVE_XPETRA_EPETRA

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    // do nothin
#else
    RCP<const Xpetra::TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>> tmp_TCrsGraph =
        Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node>>(tmp_Graph);
    if (tmp_TCrsGraph != Teuchos::null) {
      RCP<const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> G = tmp_TCrsGraph->getTpetra_CrsGraph();
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>::writeSparseGraphFile(fileName, G);
      return;
    }
#endif
#endif  // HAVE_XPETRA_TPETRA

    throw Exceptions::BadCast("Could not cast to EpetraCrsMatrix or TpetraCrsMatrix in matrix writing");
  }

  static void WriteLocal(const std::string& fileName, const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op) {
    const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>& crsOp =
        dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(Op);
    RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmp_CrsMtx = crsOp.getCrsMatrix();

    ArrayRCP<const size_t> rowptr_RCP;
    ArrayRCP<LocalOrdinal> rowptr2_RCP;
    ArrayRCP<const LocalOrdinal> colind_RCP;
    ArrayRCP<const Scalar> vals_RCP;
    tmp_CrsMtx->getAllValues(rowptr_RCP, colind_RCP, vals_RCP);

    ArrayView<const size_t> rowptr       = rowptr_RCP();
    ArrayView<const LocalOrdinal> colind = colind_RCP();
    ArrayView<const Scalar> vals         = vals_RCP();

    rowptr2_RCP.resize(rowptr.size());
    ArrayView<LocalOrdinal> rowptr2 = rowptr2_RCP();
    for (size_t j = 0; j < Teuchos::as<size_t>(rowptr.size()); j++)
      rowptr2[j] = rowptr[j];

    Teuchos::MatrixMarket::Raw::Writer<Scalar, LocalOrdinal> writer;
    writer.writeFile(fileName + "." + std::to_string(Op.getRowMap()->getComm()->getSize()) + "." + std::to_string(Op.getRowMap()->getComm()->getRank()),
                     rowptr2, colind, vals,
                     rowptr.size() - 1, Op.getColMap()->getLocalNumElements());
  }

  static void WriteBlockedCrsMatrix(const std::string& fileName, const Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, const bool& writeAllMaps = false) {
    // write all matrices with their maps
    for (size_t row = 0; row < Op.Rows(); ++row) {
      for (size_t col = 0; col < Op.Cols(); ++col) {
        auto m = Op.getMatrix(row, col);
        if (m != Teuchos::null) {  // skip empty blocks
          const bool cond = Teuchos::rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(m) == Teuchos::null;
          TEUCHOS_TEST_FOR_EXCEPTION(cond,
                                     Exceptions::BadCast,
                                     "Sub block matrix (" << row << "," << col << ") is not of type CrsMatrixWrap.");
          Write(fileName + toString(row) + toString(col) + ".m", *m, writeAllMaps);
        }
      }
    }

    // write map information of map extractors
    auto rangeMapExtractor  = Op.getRangeMapExtractor();
    auto domainMapExtractor = Op.getDomainMapExtractor();

    for (size_t row = 0; row < rangeMapExtractor->NumMaps(); ++row) {
      auto map = rangeMapExtractor->getMap(row);
      Write("subRangeMap_" + fileName + toString(row) + ".m", *map);
    }
    Write("fullRangeMap_" + fileName + ".m", *(rangeMapExtractor->getFullMap()));

    for (size_t col = 0; col < domainMapExtractor->NumMaps(); ++col) {
      auto map = domainMapExtractor->getMap(col);
      Write("subDomainMap_" + fileName + toString(col) + ".m", *map);
    }
    Write("fullDomainMap_" + fileName + ".m", *(domainMapExtractor->getFullMap()));
  }

  static Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Read(const std::string& fileName, Xpetra::UnderlyingLib lib, const RCP<const Teuchos::Comm<int>>& comm, bool binary = false) {
    if (!binary) {
      // Matrix Market file format (ASCII)
      if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
        Epetra_CrsMatrix* eA;
        const RCP<const Epetra_Comm> epcomm = Xpetra::toEpetra(comm);
        int rv                              = EpetraExt::MatrixMarketFileToCrsMatrix(fileName.c_str(), *epcomm, eA);
        if (rv != 0)
          throw Exceptions::RuntimeError("EpetraExt::MatrixMarketFileToCrsMatrix return value of " + Teuchos::toString(rv));

        RCP<Epetra_CrsMatrix> tmpA = rcp(eA);

        RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> A =
            Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tmpA);
        return A;
#else
        throw Exceptions::RuntimeError("Xpetra has not been compiled with Epetra and EpetraExt support.");
#endif
      } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
        throw Exceptions::RuntimeError("Xpetra has not been compiled with Tpetra GO=int enabled.");
#else
        typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> sparse_matrix_type;

        typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;

        bool callFillComplete = true;

        RCP<sparse_matrix_type> tA = reader_type::readSparseFile(fileName, comm, callFillComplete);

        if (tA.is_null())
          throw Exceptions::RuntimeError("The Tpetra::CrsMatrix returned from readSparseFile() is null.");

        RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmpA1 = rcp(new Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tA));
        RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmpA2       = Teuchos::rcp_implicit_cast<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tmpA1);
        RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> A              = rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tmpA2));

        return A;
#endif
#else
        throw Exceptions::RuntimeError("Xpetra has not been compiled with Tpetra support.");
#endif
      } else {
        throw Exceptions::RuntimeError("Xpetra:IO: you must specify Xpetra::UseEpetra or Xpetra::UseTpetra.");
      }
    } else {
      // Custom file format (binary)
      std::ifstream ifs(fileName.c_str(), std::ios::binary);
      TEUCHOS_TEST_FOR_EXCEPTION(!ifs.good(), Exceptions::RuntimeError, "Can not read \"" << fileName << "\"");
      int m, n, nnz;
      ifs.read(reinterpret_cast<char*>(&m), sizeof(m));
      ifs.read(reinterpret_cast<char*>(&n), sizeof(n));
      ifs.read(reinterpret_cast<char*>(&nnz), sizeof(nnz));

      int myRank = comm->getRank();

      GlobalOrdinal indexBase                                    = 0;
      RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowMap = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, m, (myRank == 0 ? m : 0), indexBase, comm), rangeMap = rowMap;
      RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> colMap = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, n, (myRank == 0 ? n : 0), indexBase, comm), domainMap = colMap;

      RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> A;

      if (myRank == 0) {
        Teuchos::Array<GlobalOrdinal> inds;
        Teuchos::Array<Scalar> vals;
        // Scan matrix to determine the exact nnz per row.
        Teuchos::ArrayRCP<size_t> numEntriesPerRow(m);
        for (int i = 0; i < m; i++) {
          int row, rownnz;
          ifs.read(reinterpret_cast<char*>(&row), sizeof(row));
          ifs.read(reinterpret_cast<char*>(&rownnz), sizeof(rownnz));
          numEntriesPerRow[i] = rownnz;
          for (int j = 0; j < rownnz; j++) {
            int index;
            ifs.read(reinterpret_cast<char*>(&index), sizeof(index));
          }
          for (int j = 0; j < rownnz; j++) {
            double value;
            ifs.read(reinterpret_cast<char*>(&value), sizeof(value));
          }
        }

        A = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap, colMap, numEntriesPerRow);

        // Now that nnz per row are known, reread and store the matrix.
        ifs.seekg(0, ifs.beg);  // rewind to beginning of file
        int junk;               // skip header info
        ifs.read(reinterpret_cast<char*>(&m), sizeof(junk));
        ifs.read(reinterpret_cast<char*>(&n), sizeof(junk));
        ifs.read(reinterpret_cast<char*>(&nnz), sizeof(junk));
        for (int i = 0; i < m; i++) {
          int row, rownnz;
          ifs.read(reinterpret_cast<char*>(&row), sizeof(row));
          ifs.read(reinterpret_cast<char*>(&rownnz), sizeof(rownnz));
          inds.resize(rownnz);
          vals.resize(rownnz);
          for (int j = 0; j < rownnz; j++) {
            int index;
            ifs.read(reinterpret_cast<char*>(&index), sizeof(index));
            inds[j] = Teuchos::as<GlobalOrdinal>(index);
          }
          for (int j = 0; j < rownnz; j++) {
            double value;
            ifs.read(reinterpret_cast<char*>(&value), sizeof(value));
            vals[j] = Teuchos::as<Scalar>(value);
          }
          A->insertGlobalValues(row, inds, vals);
        }
      }  // if (myRank == 0)
      else {
        Teuchos::ArrayRCP<size_t> numEntriesPerRow(0);
        A = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap, colMap, numEntriesPerRow);
      }

      A->fillComplete(domainMap, rangeMap);

      return A;
    }

    TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
  }

  static Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Read(const std::string& filename,
                                                                                      const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowMap,
                                                                                      RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> colMap          = Teuchos::null,
                                                                                      const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> domainMap = Teuchos::null,
                                                                                      const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rangeMap  = Teuchos::null,
                                                                                      const bool callFillComplete                                               = true,
                                                                                      const bool binary                                                         = false,
                                                                                      const bool tolerant                                                       = false,
                                                                                      const bool debug                                                          = false) {
    TEUCHOS_TEST_FOR_EXCEPTION(rowMap.is_null(), Exceptions::RuntimeError, "Utils::Read() : rowMap cannot be null");

    RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> domain = (domainMap.is_null() ? rowMap : domainMap);
    RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> range  = (rangeMap.is_null() ? rowMap : rangeMap);

    const Xpetra::UnderlyingLib lib = rowMap->lib();
    if (binary == false) {
      if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
        Epetra_CrsMatrix* eA;
        const RCP<const Epetra_Comm> epcomm = Xpetra::toEpetra(rowMap->getComm());
        const Epetra_Map& epetraRowMap      = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map2EpetraMap(*rowMap);
        const Epetra_Map& epetraDomainMap   = (domainMap.is_null() ? epetraRowMap : Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map2EpetraMap(*domainMap));
        const Epetra_Map& epetraRangeMap    = (rangeMap.is_null() ? epetraRowMap : Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map2EpetraMap(*rangeMap));
        int rv;
        if (colMap.is_null()) {
          rv = EpetraExt::MatrixMarketFileToCrsMatrix(filename.c_str(), epetraRowMap, epetraRangeMap, epetraDomainMap, eA);

        } else {
          const Epetra_Map& epetraColMap = Map2EpetraMap(*colMap);
          rv                             = EpetraExt::MatrixMarketFileToCrsMatrix(filename.c_str(), epetraRowMap, epetraColMap, epetraRangeMap, epetraDomainMap, eA);
        }

        if (rv != 0)
          throw Exceptions::RuntimeError("EpetraExt::MatrixMarketFileToCrsMatrix return value of " + Teuchos::toString(rv));

        RCP<Epetra_CrsMatrix> tmpA = rcp(eA);
        RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> A =
            Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tmpA);

        return A;
#else
        throw Exceptions::RuntimeError("Xpetra has not been compiled with Epetra and EpetraExt support.");
#endif
      } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
        throw Exceptions::RuntimeError("Xpetra has not been compiled with Tpetra GO=int support.");
#else
        typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> sparse_matrix_type;
        typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;
        typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

        const RCP<const map_type> tpetraRowMap    = Map2TpetraMap(*rowMap);
        RCP<const map_type> tpetraColMap          = (colMap.is_null() ? Teuchos::null : Map2TpetraMap(*colMap));
        const RCP<const map_type> tpetraRangeMap  = (rangeMap.is_null() ? tpetraRowMap : Map2TpetraMap(*rangeMap));
        const RCP<const map_type> tpetraDomainMap = (domainMap.is_null() ? tpetraRowMap : Map2TpetraMap(*domainMap));

        RCP<sparse_matrix_type> tA = reader_type::readSparseFile(filename, tpetraRowMap, tpetraColMap, tpetraDomainMap, tpetraRangeMap,
                                                                 callFillComplete, tolerant, debug);
        if (tA.is_null())
          throw Exceptions::RuntimeError("The Tpetra::CrsMatrix returned from readSparseFile() is null.");

        RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmpA1 = rcp(new Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tA));
        RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmpA2       = Teuchos::rcp_implicit_cast<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tmpA1);
        RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> A              = rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tmpA2));

        return A;
#endif
#else
        throw Exceptions::RuntimeError("Xpetra has not been compiled with Tpetra support.");
#endif
      } else {
        throw Exceptions::RuntimeError("Utils::Read : you must specify Xpetra::UseEpetra or Xpetra::UseTpetra.");
      }
    } else {
      // Read in on rank 0.
      auto tempA = Read(filename, lib, rowMap->getComm(), binary);

      auto A        = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap, colMap, 0);
      auto importer = Xpetra::ImportFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(tempA->getRowMap(), rowMap);
      A->doImport(*tempA, *importer, Xpetra::INSERT);
      if (callFillComplete)
        A->fillComplete(domainMap, rangeMap);

      return A;
    }

    TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
  }

  static Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> ReadLocal(const std::string& filename,
                                                                                           const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowMap,
                                                                                           RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> colMap,
                                                                                           const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> domainMap = Teuchos::null,
                                                                                           const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rangeMap  = Teuchos::null,
                                                                                           const bool callFillComplete                                               = true,
                                                                                           const bool binary                                                         = false,
                                                                                           const bool tolerant                                                       = false,
                                                                                           const bool debug                                                          = false) {
    TEUCHOS_TEST_FOR_EXCEPTION(rowMap.is_null(), Exceptions::RuntimeError, "Utils::ReadLocal() : rowMap cannot be null");
    TEUCHOS_TEST_FOR_EXCEPTION(colMap.is_null(), Exceptions::RuntimeError, "Utils::ReadLocal() : colMap cannot be null");

    using matrix_type   = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using crs_wrap_type = Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using crs_type      = Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

    RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> domain = (domainMap.is_null() ? rowMap : domainMap);
    RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> range  = (rangeMap.is_null() ? rowMap : rangeMap);

    std::string rankFilename = filename + "." + std::to_string(rowMap->getComm()->getSize()) + "." + std::to_string(rowMap->getComm()->getRank());
    RCP<matrix_type> A       = rcp(new crs_wrap_type(rowMap, colMap, 0));

    if (binary == false) {
      RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
      params->set("Parse tolerantly", tolerant);
      params->set("Debug mode", debug);

      LocalOrdinal numRows = rowMap->getLocalNumElements();
      LocalOrdinal numCols = colMap->getLocalNumElements();

      ArrayRCP<LocalOrdinal> rowptr2_RCP;
      ArrayRCP<LocalOrdinal> colind2_RCP;
      ArrayRCP<Scalar> vals2_RCP;

      Teuchos::MatrixMarket::Raw::Reader<Scalar, LocalOrdinal> reader;
      reader.readFile(rowptr2_RCP, colind2_RCP, vals2_RCP,
                      numRows, numCols,
                      rankFilename);

      RCP<crs_type> ACrs = Teuchos::rcp_dynamic_cast<crs_wrap_type>(A)->getCrsMatrix();

      ArrayRCP<size_t> rowptr_RCP;
      ArrayRCP<LocalOrdinal> colind_RCP;
      ArrayRCP<Scalar> vals_RCP;
      ACrs->allocateAllValues(colind2_RCP.size(), rowptr_RCP, colind_RCP, vals_RCP);

      rowptr_RCP.assign(rowptr2_RCP.begin(), rowptr2_RCP.end());
      colind_RCP = colind2_RCP;
      vals_RCP   = vals2_RCP;

      ACrs->setAllValues(rowptr_RCP, colind_RCP, vals_RCP);
    } else {
      // Custom file format (binary)
      std::ifstream ifs = std::ifstream(rankFilename.c_str(), std::ios::binary);
      TEUCHOS_TEST_FOR_EXCEPTION(!ifs.good(), Exceptions::RuntimeError, "Can not read \"" << filename << "\"");

      int m, n, nnz;
      ifs.read(reinterpret_cast<char*>(&m), sizeof(m));
      ifs.read(reinterpret_cast<char*>(&n), sizeof(n));
      ifs.read(reinterpret_cast<char*>(&nnz), sizeof(nnz));

      TEUCHOS_ASSERT_EQUALITY(Teuchos::as<int>(rowMap->getLocalNumElements()), m);

      Teuchos::ArrayRCP<size_t> rowptrRCP;
      Teuchos::ArrayRCP<LocalOrdinal> indicesRCP;
      Teuchos::ArrayRCP<Scalar> valuesRCP;

      RCP<crs_type> ACrs = Teuchos::rcp_dynamic_cast<crs_wrap_type>(A)->getCrsMatrix();

      ACrs->allocateAllValues(nnz, rowptrRCP, indicesRCP, valuesRCP);

      Teuchos::ArrayView<size_t> rowptr        = rowptrRCP();
      Teuchos::ArrayView<LocalOrdinal> indices = indicesRCP();
      Teuchos::ArrayView<Scalar> values        = valuesRCP();

      bool sorted = true;

      // Read in rowptr
      for (int i = 0; i < m; i++) {
        int row, rownnz;
        ifs.read(reinterpret_cast<char*>(&row), sizeof(row));
        ifs.read(reinterpret_cast<char*>(&rownnz), sizeof(rownnz));

        rowptr[row + 1] += rownnz;
        ifs.seekg(sizeof(int) * rownnz + sizeof(double) * rownnz, ifs.cur);
      }
      for (int i = 0; i < m; i++)
        rowptr[i + 1] += rowptr[i];
      TEUCHOS_ASSERT(Teuchos::as<int>(rowptr[m]) == nnz);

      // reset to where the data starts
      ifs.seekg(sizeof(int) * 3, ifs.beg);

      // read in entries
      for (int i = 0; i < m; i++) {
        int row, rownnz;
        ifs.read(reinterpret_cast<char*>(&row), sizeof(row));
        ifs.read(reinterpret_cast<char*>(&rownnz), sizeof(rownnz));
        size_t ptr = rowptr[row];
        for (int j = 0; j < rownnz; j++) {
          int index;
          ifs.read(reinterpret_cast<char*>(&index), sizeof(index));
          indices[ptr] = Teuchos::as<LocalOrdinal>(index);
          if (j > 0)
            sorted = sorted & (indices[ptr - 1] < indices[ptr]);
          ++ptr;
        }
        ptr = rowptr[row];
        for (int j = 0; j < rownnz; j++) {
          double value;
          ifs.read(reinterpret_cast<char*>(&value), sizeof(value));
          values[ptr] = Teuchos::as<Scalar>(value);
          ++ptr;
        }
        rowptr[row] += rownnz;
      }
      for (int i = m; i > 0; i--)
        rowptr[i] = rowptr[i - 1];
      rowptr[0] = 0;

#ifdef HAVE_XPETRA_TPETRA
      if (!sorted) {
        for (LocalOrdinal lclRow = 0; lclRow < m; lclRow++) {
          size_t rowBegin = rowptr[lclRow];
          size_t rowEnd   = rowptr[lclRow + 1];
          Tpetra::sort2(&indices[rowBegin], &indices[rowEnd], &values[rowBegin]);
        }
      }
#else
      TEUCHOS_ASSERT(sorted);
#endif

      ACrs->setAllValues(rowptrRCP, indicesRCP, valuesRCP);
    }

    if (callFillComplete)
      A->fillComplete(domainMap, rangeMap);
    return A;
  }

  static RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> ReadMultiVector(const std::string& fileName,
                                                                                             const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
                                                                                             const bool binary = false) {
    Xpetra::UnderlyingLib lib = map->lib();

    if (lib == Xpetra::UseEpetra) {
      // taw: Oct 9 2015: do we need a specialization for <double,int,int>??
      // TEUCHOS_TEST_FOR_EXCEPTION(true, ::Xpetra::Exceptions::BadCast, "Epetra can only be used with Scalar=double and Ordinal=int");
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
      TEUCHOS_ASSERT(!binary);
      Epetra_MultiVector* MV;
      int rv = EpetraExt::MatrixMarketFileToMultiVector(fileName.c_str(), toEpetra(map), MV);
      if (rv != 0) throw Exceptions::RuntimeError("EpetraExt::MatrixMarketFileToMultiVector failed");
      RCP<Epetra_MultiVector> MVrcp = rcp(MV);
      return Convert_Epetra_MultiVector_ToXpetra_MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(MVrcp);
#else
      throw Exceptions::RuntimeError("Xpetra has not been compiled with Epetra and EpetraExt support.");
#endif
    } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
      throw Exceptions::RuntimeError("Xpetra has not been compiled with Tpetra GO=int support.");
#else
      typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> sparse_matrix_type;
      typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;
      typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
      typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> multivector_type;

      RCP<const map_type> temp                                                = toTpetra(map);
      RCP<multivector_type> TMV                                               = reader_type::readDenseFile(fileName, map->getComm(), temp, false, false, binary);
      RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> rmv = Xpetra::toXpetra(TMV);
      return rmv;
#endif
#else
      throw Exceptions::RuntimeError("Xpetra has not been compiled with Tpetra support.");
#endif
    } else {
      throw Exceptions::RuntimeError("Utils::Read : you must specify Xpetra::UseEpetra or Xpetra::UseTpetra.");
    }

    TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
  }

  static RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> ReadMultiVectorLO(const std::string& fileName,
                                                                                                     const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
                                                                                                     const bool binary = false) {
    Xpetra::UnderlyingLib lib = map->lib();

    if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
      typedef Tpetra::CrsMatrix<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> sparse_matrix_type;
      typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;
      typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
      typedef Tpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> multivector_type;

      RCP<const map_type> temp                                                      = toTpetra(map);
      RCP<multivector_type> TMV                                                     = reader_type::readDenseFile(fileName, map->getComm(), temp, false, false, binary);
      RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> rmv = Xpetra::toXpetra(TMV);
      return rmv;
#else
      throw Exceptions::RuntimeError("Xpetra has not been compiled with Tpetra support.");
#endif
    } else {
      throw Exceptions::RuntimeError("Utils::ReadMultiVectorLO : only implemented for Tpetra");
    }

    TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
  }

  static RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> ReadMap(const std::string& fileName,
                                                                           Xpetra::UnderlyingLib lib,
                                                                           const RCP<const Teuchos::Comm<int>>& comm,
                                                                           const bool binary = false) {
    if (lib == Xpetra::UseEpetra) {
      // do we need another specialization for <double,int,int> ??
      // TEUCHOS_TEST_FOR_EXCEPTION(true, ::Xpetra::Exceptions::BadCast, "Epetra can only be used with Scalar=double and Ordinal=int");
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
      TEUCHOS_ASSERT(!binary);
      Epetra_Map* eMap;
      int rv = EpetraExt::MatrixMarketFileToMap(fileName.c_str(), *(Xpetra::toEpetra(comm)), eMap);
      if (rv != 0)
        throw Exceptions::RuntimeError("Error reading map from file " + fileName + " with EpetraExt::MatrixMarketToMap (returned " + Teuchos::toString(rv) + ")");

      RCP<Epetra_Map> eMap1 = rcp(new Epetra_Map(*eMap));
      return Xpetra::toXpetra<int, Node>(*eMap1);
#else
      throw Exceptions::RuntimeError("Xpetra has not been compiled with Epetra and EpetraExt support.");
#endif
    } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
      throw Exceptions::RuntimeError("Xpetra has not been compiled with Tpetra GO=int support.");
#else
      typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> sparse_matrix_type;
      typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;

      RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> tMap = reader_type::readMapFile(fileName, comm, false, false, binary);
      if (tMap.is_null())
        throw Exceptions::RuntimeError("The Tpetra::Map returned from readSparseFile() is null.");

      return Xpetra::toXpetra(tMap);
#endif
#else
      throw Exceptions::RuntimeError("Xpetra has not been compiled with Tpetra support.");
#endif
    } else {
      throw Exceptions::RuntimeError("Utils::Read : you must specify Xpetra::UseEpetra or Xpetra::UseTpetra.");
    }

    TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
  }

  static RCP<const Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> ReadBlockedCrsMatrix(const std::string& fileName, Xpetra::UnderlyingLib lib, const RCP<const Teuchos::Comm<int>>& comm) {
    size_t numBlocks = 2;  // TODO user parameter?

    std::vector<RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>> rangeMapVec;
    for (size_t row = 0; row < numBlocks; ++row) {
      auto map = ReadMap("subRangeMap_" + fileName + toString(row) + ".m", lib, comm);
      rangeMapVec.push_back(map);
    }
    auto fullRangeMap = ReadMap("fullRangeMap_" + fileName + ".m", lib, comm);

    std::vector<RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>> domainMapVec;
    for (size_t col = 0; col < numBlocks; ++col) {
      auto map = ReadMap("subDomainMap_" + fileName + toString(col) + ".m", lib, comm);
      domainMapVec.push_back(map);
    }
    auto fullDomainMap = ReadMap("fullDomainMap_" + fileName + ".m", lib, comm);

    /*std::vector<RCP<const XpMap> > testRgMapVec;
    for(size_t r = 0; r < numBlocks; ++r) {
      RCP<const XpMap> map = ReadMap("rangemap_" + fileName + XpIO::toString<size_t>(r) + "0.m", lib, comm);
      testRgMapVec.push_back(map);
    }
    std::vector<RCP<const XpMap> > testDoMapVec;
    for(size_t c = 0; c < numBlocks; ++c) {
      RCP<const XpMap> map = ReadMap("domainmap_" + fileName + "0" + XpIO::toString<size_t>(c) + ".m", lib, comm);
      testDoMapVec.push_back(map);
    }*/

    // create map extractors

    // range map extractor
    bool bRangeUseThyraStyleNumbering = false;
    /*
    GlobalOrdinal gMinGids = 0;
    for(size_t v = 0; v < testRgMapVec.size(); ++v) {
      gMinGids += testRgMapVec[v]->getMinAllGlobalIndex();
    }
    if ( gMinGids==0 && testRgMapVec.size() > 1 ) bRangeUseThyraStyleNumbering = true;*/
    RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>> rangeMapExtractor =
        rcp(new Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>(fullRangeMap, rangeMapVec, bRangeUseThyraStyleNumbering));

    // domain map extractor
    bool bDomainUseThyraStyleNumbering = false;
    /*gMinGids = 0;
    for(size_t v = 0; v < testDoMapVec.size(); ++v) {
      gMinGids += testDoMapVec[v]->getMinAllGlobalIndex();
    }
    if ( gMinGids==0 && testDoMapVec.size() > 1) bDomainUseThyraStyleNumbering = true;*/
    RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>> domainMapExtractor =
        rcp(new Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>(fullDomainMap, domainMapVec, bDomainUseThyraStyleNumbering));

    auto bOp = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rangeMapExtractor, domainMapExtractor, 33));

    // Read all matrices with their maps and create the BlockedCrsMatrix
    for (size_t row = 0; row < numBlocks; ++row) {
      for (size_t col = 0; col < numBlocks; ++col) {
        auto rowSubMap = ReadMap("rowmap_" + fileName + toString(row) + toString(col) + ".m", lib, comm);
        auto colSubMap = ReadMap("colmap_" + fileName + toString(row) + toString(col) + ".m", lib, comm);
        auto domSubMap = ReadMap("domainmap_" + fileName + toString(row) + toString(col) + ".m", lib, comm);
        auto ranSubMap = ReadMap("rangemap_" + fileName + toString(row) + toString(col) + ".m", lib, comm);
        auto mat       = Read(fileName + toString(row) + toString(col) + ".m", rowSubMap, colSubMap, domSubMap, ranSubMap);
        bOp->setMatrix(row, col, mat);
      }
    }

    bOp->fillComplete();

    return bOp;
  }

  //! Little helper function to convert non-string types to strings
  template <class T>
  static std::string toString(const T& what) {
    std::ostringstream buf;
    buf << what;
    return buf.str();
  }
};

#endif  // HAVE_XPETRA_EPETRA

}  // end namespace Xpetra

#define XPETRA_IO_SHORT

#endif /* PACKAGES_XPETRA_SUP_UTILS_XPETRA_IO_HPP_ */
