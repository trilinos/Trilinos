// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_BINARYIO_DECL_HPP
#define TPETRA_BINARYIO_DECL_HPP

/// \file Tpetra_BinaryIO_decl.hpp
/// \brief Declarations of scalable binary file readers and writers for Tpetra objects.

#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"

#include <cstddef>
#include <string>

namespace Tpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class BinaryIO {
 public:
  using sparse_matrix_type = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using sparse_graph_type  = Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
  using multivector_type   = Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using vector_type        = Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using map_type           = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using trcp_tcomm_t       = Teuchos::RCP<const Teuchos::Comm<int>>;

  static void writeMapFile(const std::string& filename, const map_type& map);
  static Teuchos::RCP<const map_type> readMapFile(const std::string& filename, const trcp_tcomm_t& comm);

  static void writeDenseFile(const std::string& filename, const multivector_type& X);
  static Teuchos::RCP<multivector_type> readDenseFile(const std::string& filename, const trcp_tcomm_t& comm);
  static Teuchos::RCP<multivector_type> readDenseFile(const std::string& filename, const Teuchos::RCP<const map_type>& map);

  static void writeVectorFile(const std::string& filename, const vector_type& X);
  static Teuchos::RCP<vector_type> readVectorFile(const std::string& filename, const trcp_tcomm_t& comm);
  static Teuchos::RCP<vector_type> readVectorFile(const std::string& filename, const Teuchos::RCP<const map_type>& map);

  static void writeSparseFile(const std::string& filename, const sparse_matrix_type& A);
  static Teuchos::RCP<sparse_matrix_type> readSparseFile(const std::string& filename, const trcp_tcomm_t& comm);
  static Teuchos::RCP<sparse_matrix_type>
  readSparseFile(const std::string& filename,
                 const Teuchos::RCP<const map_type>& rowMap,
                 const Teuchos::RCP<const map_type>& domainMap = Teuchos::null,
                 const Teuchos::RCP<const map_type>& rangeMap  = Teuchos::null,
                 const bool callFillComplete                   = true);

 private:
  enum ObjectKind {
    MAP_OBJECT         = 1,
    MULTIVECTOR_OBJECT = 2,
    VECTOR_OBJECT      = 3,
    MATRIX_OBJECT      = 4
  };

  enum TypeFlags {
    TYPE_IS_ORDINAL  = 1u << 0,
    TYPE_IS_FLOATING = 1u << 1,
    TYPE_IS_COMPLEX  = 1u << 2,
    TYPE_IS_SIGNED   = 1u << 3
  };

  struct FileHeader {
    char magic[8];
    unsigned long long version;
    unsigned long long byteOrderMarker;
    unsigned long long objectKind;
    unsigned long long scalarSize;
    unsigned long long scalarFlags;
    unsigned long long localOrdinalSize;
    unsigned long long localOrdinalFlags;
    unsigned long long globalOrdinalSize;
    unsigned long long globalOrdinalFlags;
    unsigned long long objectFlags;
    unsigned long long numGlobalRows;
    unsigned long long numGlobalCols;
    unsigned long long numGlobalEntries;
    unsigned long long numVectors;
    unsigned long long rowMapOffset;
    unsigned long long domainMapOffset;
    unsigned long long rangeMapOffset;
    unsigned long long rowPtrOffset;
    unsigned long long columnIndicesOffset;
    unsigned long long valuesOffset;
  };

  struct MapSectionHeader {
    unsigned long long numGlobalElements;
    long long indexBase;
    unsigned long long mapFlags;
  };

  static const unsigned long long fileVersion     = 1ull;
  static const unsigned long long byteOrderMarker = 0x0102030405060708ull;

  static FileHeader makeBaseHeader(const unsigned long long objectKind);
  static void validateHeader(const FileHeader& header, const unsigned long long expectedObjectKind);

  template <class T>
  static unsigned long long typeFlags();

  static unsigned long long mapFlags(const map_type& map);

  static unsigned long long computeContiguousLocalCount(const unsigned long long globalCount, const int rank, const int size);
  static unsigned long long computeContiguousStart(const unsigned long long globalCount, const int rank, const int size);
  static unsigned long long checkedByteCount(const unsigned long long count, const size_t elementSize);
  static unsigned long long checkedByteOffset(const unsigned long long dataOffset,
                                              const unsigned long long globalOffset,
                                              const size_t elementSize);
  static void broadcastBytesFromRoot(char data[],
                                     const unsigned long long byteCount,
                                     const trcp_tcomm_t& comm);

  static void checkFileOpen(const bool success, const std::string& filename, const trcp_tcomm_t& comm, const char mode[]);

  static void writeHeaderToNewFile(const std::string& filename, const FileHeader& header, const trcp_tcomm_t& comm);
  static FileHeader readHeaderFromFile(const std::string& filename, const trcp_tcomm_t& comm);

  static MapSectionHeader makeMapSectionHeader(const map_type& map);
  static unsigned long long mapSectionSize();
  static unsigned long long mapSectionPayloadOffset(const unsigned long long mapSectionOffset);

  static void writeMapSection(const std::string& filename,
                              const unsigned long long mapSectionOffset,
                              const map_type& map,
                              const trcp_tcomm_t& comm);
  static Teuchos::RCP<const map_type> readMapSection(const std::string& filename,
                                                     const unsigned long long mapSectionOffset,
                                                     const trcp_tcomm_t& comm);

  static void writeMapSectionHeader(const std::string& filename,
                                    const unsigned long long mapSectionOffset,
                                    const MapSectionHeader& header,
                                    const trcp_tcomm_t& comm);
  static MapSectionHeader readMapSectionHeader(const std::string& filename,
                                               const unsigned long long mapSectionOffset,
                                               const trcp_tcomm_t& comm);

  static unsigned long long exclusiveScanUnsignedLongLong(const unsigned long long localValue, const trcp_tcomm_t& comm);

  template <class T>
  static void writeArrayCollective(const std::string& filename,
                                   const unsigned long long dataOffset,
                                   const T* data,
                                   const unsigned long long count,
                                   const unsigned long long globalOffset,
                                   const trcp_tcomm_t& comm);

  template <class T>
  static void readArrayCollective(const std::string& filename,
                                  const unsigned long long dataOffset,
                                  T* data,
                                  const unsigned long long count,
                                  const unsigned long long globalOffset,
                                  const trcp_tcomm_t& comm);
};

}  // namespace Tpetra

#endif  // TPETRA_BINARYIO_DECL_HPP
