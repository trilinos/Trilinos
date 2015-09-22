// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER


#ifndef TPETRA_MMHELPERS_DECL_HPP
#define TPETRA_MMHELPERS_DECL_HPP

#include <Tpetra_CrsMatrix.hpp>
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
template <class Scalar = Details::DefaultTypes::scalar_type,
          class LocalOrdinal = Details::DefaultTypes::local_ordinal_type,
          class GlobalOrdinal = Details::DefaultTypes::global_ordinal_type,
          class Node = Details::DefaultTypes::node_type>
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

template <class Scalar = CrsMatrix<>::scalar_type,
          class LocalOrdinal = typename CrsMatrix<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename CrsMatrix<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node = typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
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


template <class Scalar = CrsMatrix<>::scalar_type,
          class LocalOrdinal = typename CrsMatrix<Scalar>::local_ordinal_type,
          class GlobalOrdinal =
            typename CrsMatrix<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node =
            typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
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

