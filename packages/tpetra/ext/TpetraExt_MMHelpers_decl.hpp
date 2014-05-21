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

#include "Tpetra_ConfigDefs.hpp"
#include "Teuchos_Array.hpp"
#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>
#include <set>
#include <map>


/*! \file TpetraExt_MMHelpers_decl.hpp

    The declarations for the class Tpetra::MMMultiMultiply and related non-member constructors.
 */

namespace Tpetra {
#ifndef DOXYGEN_SHOULD_SKIP_THIS  
	// forward declaration
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatOps>
class CrsMatrix;

template <class LocalOrdinal, class GlobalOrdinal, class Node>
class Map;
#endif
//struct that holds views of the contents of a CrsMatrix. These
//contents may be a mixture of local and remote rows of the
//actual matrix.
template <class Scalar, 
	class LocalOrdinal=int, 
	class GlobalOrdinal=LocalOrdinal, 
	class Node=KokkosClassic::DefaultNode::DefaultNodeType, 
	class SpMatOps= typename KokkosClassic::DefaultKernels<Scalar, LocalOrdinal, Node>::SparseOps >
class CrsMatrixStruct {
public:
  CrsMatrixStruct();

  virtual ~CrsMatrixStruct();

  void deleteContents();

  /** \brief Which of the desired global rows are remote*/
  Teuchos::Array<bool> remote;

  /** \brief Original row map of matrix */
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > origRowMap;
  /** \brief Desired row map for "imported" version of the matrix */
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > rowMap;
  /** \brief Col map for the original version of the matrix */
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > colMap;
  /** \brief Domain map for original matrix */
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > domainMap;
  /** \brief Colmap garnered as a result of the import */
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > importColMap;
  /** \brief The imported matrix */
  Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> >  importMatrix;  
  /** \brief The original matrix */
  Teuchos::RCP<const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> >  origMatrix;
  
};

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatOps>
int dumpCrsMatrixStruct(const CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps >& M);

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class CrsWrapper {
 public:
  virtual ~CrsWrapper(){}

  virtual Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > getRowMap() const = 0;

  virtual bool isFillComplete() = 0;

  virtual void insertGlobalValues(GlobalOrdinal globalRow, const Teuchos::ArrayView<const GlobalOrdinal> &indices, const Teuchos::ArrayView<const Scalar> &values) = 0;

  virtual void sumIntoGlobalValues(GlobalOrdinal globalRow, const Teuchos::ArrayView<const GlobalOrdinal> &indices, const Teuchos::ArrayView<const Scalar> &values) = 0;
};

template <class Scalar, 
	class LocalOrdinal=int, 
	class GlobalOrdinal=LocalOrdinal, 
	class Node=KokkosClassic::DefaultNode::DefaultNodeType, 
	class SpMatOps= typename KokkosClassic::DefaultKernels<Scalar, LocalOrdinal, Node>::SparseOps >
class CrsWrapper_CrsMatrix : public CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node>{
 public:
  CrsWrapper_CrsMatrix(CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps >& crsmatrix);
  virtual ~CrsWrapper_CrsMatrix();

  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > getRowMap() const;

  bool isFillComplete();

  void insertGlobalValues(GlobalOrdinal globalRow, const Teuchos::ArrayView<const GlobalOrdinal> &indices, const Teuchos::ArrayView<const Scalar> &values);
  void sumIntoGlobalValues(GlobalOrdinal globalRow, const Teuchos::ArrayView<const GlobalOrdinal> &indices, const Teuchos::ArrayView<const Scalar> &values);

 private:
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& crsmat_;
};

template<class Scalar,
	class LocalOrdinal=int,
	class GlobalOrdinal=LocalOrdinal,
	class Node=KokkosClassic::DefaultNode::DefaultNodeType>
class CrsWrapper_GraphBuilder : public CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node>{
 public:
  CrsWrapper_GraphBuilder(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& map);
  virtual ~CrsWrapper_GraphBuilder();

  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > getRowMap() const {return rowmap_; }

  bool isFillComplete();

  void insertGlobalValues(GlobalOrdinal globalRow, const Teuchos::ArrayView<const GlobalOrdinal> &indices, const Teuchos::ArrayView<const Scalar> &values);
  void sumIntoGlobalValues(GlobalOrdinal globalRow, const Teuchos::ArrayView<const GlobalOrdinal> &indices, const Teuchos::ArrayView<const Scalar> &values);

  std::map<GlobalOrdinal,std::set<GlobalOrdinal>*>& get_graph();

  size_t get_max_row_length() { return max_row_length_; }

 private:
  std::map<GlobalOrdinal,std::set<GlobalOrdinal>*> graph_;
  const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rowmap_;
  global_size_t max_row_length_;
};

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatOps>
void insert_matrix_locations(CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node>& graphbuilder,
                              CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& C);



}
#endif // TPETRA_MMHELPERS_DECL_HPP

