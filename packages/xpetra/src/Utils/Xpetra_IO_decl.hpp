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

#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <TpetraExt_MatrixMatrix.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#include <Xpetra_TpetraCrsGraph.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_TpetraBlockCrsMatrix.hpp>
#include "Tpetra_Util.hpp"

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

}  // end namespace Xpetra

#define XPETRA_IO_SHORT

#endif /* PACKAGES_XPETRA_SUP_UTILS_XPETRA_IO_HPP_ */
