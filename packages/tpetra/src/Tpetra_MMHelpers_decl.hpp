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
#include <set>
#include <map>


/*! \file Tpetra_MMMultiply_decl.hpp 

    The declarations for the class Tpetra::MMMultiMultiply and related non-member constructors.
 */

namespace Tpetra {
namespace Tpetra {
#ifndef DOXYGEN_SHOULD_SKIP_THIS  
	// forward declaration
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class CrsMatrix;
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class Map;
#endif
//struct that holds views of the contents of a CrsMatrix. These
//contents may be a mixture of local and remote rows of the
//actual matrix.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class CrsMatrixStruct {
public:
  CrsMatrixStruct();

  virtual ~CrsMatrixStruct();

  void deleteContents();

  size_t numRows;
  size_t* numEntriesPerRow;
  LocalOrdinal** indices;
  Scalar** values;
  bool* remote;
  size_t numRemote;
  const Map<LocalOrdinal, GlobalOrdinal, Node>* origRowMap;
  const Map<LocalOrdinal, GlobalOrdinal, Node>* rowMap;
  const Map<LocalOrdinal, GlobalOrdinal, Node>* colMap;
  const Map<LocalOrdinal, GlobalOrdinal, Node>* domainMap;
  const Map<LocalOrdinal, GlobalOrdinal, Node>* importColMap;
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>* importMatrix;
};

int dumpCrsMatrixStruct(const CrsMatrixStruct& M);
/*
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class CrsWrapper {
 public:
  virtual ~CrsWrapper(){}

  virtual const Map<LocalOrdinal, GlobalOrdinal, Node>& RowMap() const = 0;

  virtual bool Filled() = 0;

  virtual int InsertGlobalValues(GlobalOrdinal GlobalRow, size_t NumEntries, Scalar* Values, GlobalOrdinal* Indices) = 0;

  virtual int SumIntoGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices) = 0;
};

class CrsWrapper_Epetra_CrsMatrix : public CrsWrapper {
 public:
  CrsWrapper_Epetra_CrsMatrix(Epetra_CrsMatrix& epetracrsmatrix);
  virtual ~CrsWrapper_Epetra_CrsMatrix();

  const Epetra_Map& RowMap() const;

  bool Filled();

  int InsertGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);
  int SumIntoGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);

 private:
  Epetra_CrsMatrix& ecrsmat_;
};

class CrsWrapper_GraphBuilder : public CrsWrapper {
 public:
  CrsWrapper_GraphBuilder(const Epetra_Map& emap);
  virtual ~CrsWrapper_GraphBuilder();

  const Epetra_Map& RowMap() const {return rowmap_; }

  bool Filled();

  int InsertGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);
  int SumIntoGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);

  std::map<int,std::set<int>*>& get_graph();

  int get_max_row_length() { return max_row_length_; }

 private:
  std::map<int,std::set<int>*> graph_;
  const Epetra_Map& rowmap_;
  int max_row_length_;
};

void insert_matrix_locations(CrsWrapper_GraphBuilder& graphbuilder,
                              Epetra_CrsMatrix& C);
*/


}
#endif // TPETRA_MMHELPERS_DECL_HPP

