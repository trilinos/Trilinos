//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
//                 Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER


#ifndef TPETRA_MMHELPERS_DECL_HPP
#define TPETRA_MMHELPERS_DECL_HPP

#include "Tpetra_ConfigDefs.hpp"
#include "Teuchos_Array.hpp"
#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>
#include <set>
#include <map>


/*! \file Tpetra_MMHelpers_decl.hpp 

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
	class Node=Kokkos::DefaultNode::DefaultNodeType, 
	class SpMatOps= typename Kokkos::DefaultKernels<Scalar, LocalOrdinal, Node>::SparseOps >
class CrsMatrixStruct {
public:
  CrsMatrixStruct();

  virtual ~CrsMatrixStruct();

  void deleteContents();

  //The maximum number of row entries a.k.a. the longest of the indice arrays.
  typename Array<ArrayView<const LocalOrdinal> >::size_type maxNumRowEntries;
  
  /** \brief Number of local rows that the imported version of the matrix has */
  size_t numRows;
  /** \brief Number of entries in each row of the imported version of the matrix */
  Teuchos::Array<size_t> numEntriesPerRow;
  /** \brief Indicies of entries in each row of the imported version of the matrix */
  Teuchos::Array<Teuchos::ArrayView<const LocalOrdinal> > indices;
  /** \brief Values of entries in each row of the imported version of the matrix */
  Teuchos::Array<Teuchos::ArrayView<const Scalar> > values;
  /** \brief Which of the desired global rows are remote*/
  Teuchos::Array<bool> remote;
  /** \brief number of rows in the original matrix that are remote (based on the rows that this proc actually needs*/
  global_size_t numRemote;
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
	class Node=Kokkos::DefaultNode::DefaultNodeType, 
	class SpMatOps= typename Kokkos::DefaultKernels<Scalar, LocalOrdinal, Node>::SparseOps >
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
	class Node=Kokkos::DefaultNode::DefaultNodeType>
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

