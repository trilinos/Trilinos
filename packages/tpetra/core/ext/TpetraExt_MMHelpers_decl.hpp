// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_MMHELPERS_DECL_HPP
#define TPETRA_MMHELPERS_DECL_HPP

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_BlockCrsMatrix.hpp>
#include <Teuchos_Array.hpp>
#include <map>
#include <set>

/// \file TpetraExt_MMHelpers_decl.hpp
/// \brief Declaration of Tpetra::MMMultiMultiply and nonmember constructors.
///
/// \warning This is an implementation detail of Tpetra.  Please don't
///   rely on anything in this file.

namespace Tpetra {

/// \class CrsMatrixStruct
/// \brief Struct that holds views of the contents of a CrsMatrix.
///
/// These contents may be a mixture of local and remote rows of the
/// actual matrix.
template <class Scalar = ::Tpetra::Details::DefaultTypes::scalar_type,
          class LocalOrdinal = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
          class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
          class Node = ::Tpetra::Details::DefaultTypes::node_type>
class CrsMatrixStruct {
public:
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
  typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;

  CrsMatrixStruct ();

  virtual ~CrsMatrixStruct ();

  void deleteContents ();

  /** \brief Original row map of matrix */
  Teuchos::RCP<const map_type> origRowMap;
  /** \brief Desired row map for "imported" version of the matrix */
  Teuchos::RCP<const map_type> rowMap;
  /** \brief Col map for the original version of the matrix */
  Teuchos::RCP<const map_type> colMap;
  /** \brief Domain map for original matrix */
  Teuchos::RCP<const map_type> domainMap;
  /** \brief Colmap garnered as a result of the import */
  Teuchos::RCP<const map_type> importColMap;
  /** \brief The imported matrix */
  Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >  importMatrix;
  /** \brief The original matrix */
  Teuchos::RCP<const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >  origMatrix;

};

/// \class BlockCrsMatrixStruct
/// \brief Struct that holds views of the contents of a BlockCrsMatrix.
///
/// These contents may be a mixture of local and remote rows of the
/// actual matrix.
template <class Scalar = ::Tpetra::Details::DefaultTypes::scalar_type,
          class LocalOrdinal = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
          class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
          class Node = ::Tpetra::Details::DefaultTypes::node_type>
class BlockCrsMatrixStruct {
public:
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
  typedef BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> block_crs_matrix_type;

  BlockCrsMatrixStruct (const LocalOrdinal blocksize_);

  virtual ~BlockCrsMatrixStruct ();

  void deleteContents ();

  /** \brief Original row map of matrix */
  Teuchos::RCP<const map_type> origRowMap;
  /** \brief Desired row map for "imported" version of the matrix */
  Teuchos::RCP<const map_type> rowMap;
  /** \brief Col map for the original version of the matrix */
  Teuchos::RCP<const map_type> colMap;
  /** \brief Domain map for original matrix */
  Teuchos::RCP<const map_type> domainMap;
  /** \brief Colmap garnered as a result of the import */
  Teuchos::RCP<const map_type> importColMap;
  /** \brief The imported matrix */
  Teuchos::RCP<BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >  importMatrix;
  /** \brief The original matrix */
  Teuchos::RCP<const BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >  origMatrix;
  /** \brief The blocksize of all matrices */
  const LocalOrdinal blocksize;
};

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int
dumpCrsMatrixStruct (const CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node >& M);


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class CrsWrapper {
public:
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

  virtual ~CrsWrapper () {}
  virtual Teuchos::RCP<const map_type> getRowMap () const = 0;
  virtual bool isFillComplete () = 0;

  virtual void
  insertGlobalValues (GlobalOrdinal globalRow,
                      const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                      const Teuchos::ArrayView<const Scalar> &values) = 0;
  virtual void
  sumIntoGlobalValues (GlobalOrdinal globalRow,
                       const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                       const Teuchos::ArrayView<const Scalar> &values) = 0;
};

template <class Scalar = ::Tpetra::Details::DefaultTypes::scalar_type,
          class LocalOrdinal = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
          class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
          class Node = ::Tpetra::Details::DefaultTypes::node_type>
class CrsWrapper_CrsMatrix :
    public CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
public:
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
  typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;

  CrsWrapper_CrsMatrix (crs_matrix_type& crsmatrix);
  virtual ~CrsWrapper_CrsMatrix ();
  Teuchos::RCP<const map_type> getRowMap () const;

  bool isFillComplete ();

  void
  insertGlobalValues (GlobalOrdinal globalRow,
                      const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                      const Teuchos::ArrayView<const Scalar> &values);
  void
  sumIntoGlobalValues (GlobalOrdinal globalRow,
                       const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                       const Teuchos::ArrayView<const Scalar> &values);
private:
  crs_matrix_type& crsmat_;
};

template <class Scalar = ::Tpetra::Details::DefaultTypes::scalar_type,
          class LocalOrdinal = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
          class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
          class Node = ::Tpetra::Details::DefaultTypes::node_type>
class CrsWrapper_GraphBuilder :
    public CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
public:
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

  CrsWrapper_GraphBuilder (const Teuchos::RCP<const map_type>& map);
  virtual ~CrsWrapper_GraphBuilder ();

  Teuchos::RCP<const map_type> getRowMap () const {
    return rowmap_;
  }

  bool isFillComplete ();
  void
  insertGlobalValues (GlobalOrdinal globalRow,
                      const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                      const Teuchos::ArrayView<const Scalar> &values);
  void
  sumIntoGlobalValues (GlobalOrdinal globalRow,
                       const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                       const Teuchos::ArrayView<const Scalar> &values);

  std::map<GlobalOrdinal, std::set<GlobalOrdinal>*>& get_graph ();

  size_t get_max_row_length () {
    return max_row_length_;
  }

 private:
  std::map<GlobalOrdinal, std::set<GlobalOrdinal>*> graph_;
  const Teuchos::RCP<const map_type>& rowmap_;
  global_size_t max_row_length_;
};

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
insert_matrix_locations (CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node>& graphbuilder,
                         CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C);

} // namespace Tpetra
#endif // TPETRA_MMHELPERS_DECL_HPP

