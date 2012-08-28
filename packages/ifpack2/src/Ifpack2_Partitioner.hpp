//@HEADER
// ***********************************************************************
// 
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
// ***********************************************************************
//@HEADER

#ifndef IFPACK2_PARTITIONER_HPP
#define IFPACK2_PARTITIONER_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include <iostream>

namespace Ifpack2 {

//! Ifpack2::Partitioner: A class to decompose local Ifpack2::Graph objects.

/*!
 
  Class Ifpack2::Partitioner enables the decomposition of a local
  Ifpack2::Graph. It is supposed that the graph refers to
  a localized matrix (that is, a matrix that has been filtered
  through Ifpack2::LocalFilter).
  
  The overloaded operator (int i) can be used to extract the local partition
  ID of local row i.  
  
  The partitions created by Ifpack2::Partitioner derived clased 
  are non-overlapping in graph sense. This means that each row
  (or, more approriately, vertex)
  of \c G is assigned to exactly one partition.

  Partitioner can be extended using the functionalities of class
  Ifpack2_OverlappingPartitioner (itself derived from Ifpack2_Partitioner.
  This class extends the non-overlapping partitions by the required
  amount of overlap, considering local nodes only (that is, this
  overlap do \e not modify the overlap among the processes).

  Ifpack2::Partitioner is a pure virtual class. Concrete implementations
  are:
  - Ifpack2_LinearPartitioner, which allows the decomposition of the
    rows of the graph in simple consecutive chunks;
  - Ifpack2_Zoltan2Partitioner, which calls Zoltan2 to decompose the graph

  Generically, a constructor requires a Tpetra::RowGraph object. 
  
  <P>An example use of a Ifpack2::Partitioner derived class is as follows:  
  \code
#include "Ifpack2_Partitioner.hpp"
#include "Ifpack2_LinearPartitioner.hpp"
#include "Ifpack2_Graph.hpp"
#include "Ifpack2_Graph_Tpetra_CrsGraph.hpp"
...
Tpetra_CrsMatrix* A;         // A is filled
// create the wrapper from Tpetra_CrsGraph
Ifpack2_Graph* Graph = new Ifpack2_Graph_Tpetra_CrsGraph(A);

// we aim to create non-overlapping partitions only
Ifpack2_Partitioner Partitioner(Graph);

Ifpack2_Partitioner* Partitioner;
Partitioner = new Ifpack2_Graph_Tpetra_CrsGraph(&A);

// we want 16 local parts
List.set("partitioner: local parts", 16);
// and an overlap of 0 among the local parts (default option)
List.set("partitioner: overlap", 0);

// decompose the graph
Partitioner.Create(List);

// now Graph can be deleted, as Partitioner contains all the
// necessary information to use the partitions
delete Graph;

// we can get the number of parts actually created...
int NumParts = Partitioner.NumParts();

// ... and the number of rows in each of them
for (int i = 0 ; i < NumParts ; ++i) {
  cout << "rows in " << i << "=" << Partitioner.RowsInPart(i);
}  

// .. and, for non-overlapping partitions only, the partition ID 
// for each local row simply using:
for (int i = 0 ; i < A->NumMyRows() ; ++i)
  cout << "Partition[" << i <<"] = " << Partitioner(i) << endl;

\endcode
  
When overlapping partitiones are created, the user can get the 
row ID contained in each partition as follows:
\code
for (int i = 0 ; i < NumParts ; ++i) {
  for (int j = 0 ; j < Partitioner.RowsInPart(i) ; ++j) {
    cout << "Partition " << i << ", contains local row "
         << Partitioner(i,j) << endl;
  }
}  
\endcode
  
Ifpack2_Partitioner is used to create the subblocks in 
Ifpack2_BlockJacobi, Ifpack2_BlockGaussSeidel, and 
Ifpack2_BlockSymGaussSeidel.

\author Michael Heroux, SNL 9214.

\date Last modified on Nov-04.

*/  
template <class GraphType>
class Partitioner : public Teuchos::Describable {
public:
  typedef typename GraphType::local_ordinal_type LocalOrdinal;
  typedef typename GraphType::global_ordinal_type GlobalOrdinal;
  typedef typename GraphType::node_type Node;

  //! Destructor.
  virtual ~Partitioner() {};

  //! Returns the number of computed local partitions.
  virtual size_t numLocalParts() const = 0;

  //! Returns the overlapping level.
  virtual size_t overlappingLevel() const = 0;

  //! Returns the local non-overlapping partition ID of the specified row.
  /*! Returns the non-overlapping partition ID of the specified row.
   \param 
   MyRow - (In) local row number

   \return
   Local ID of non-overlapping partition for \t MyRow.
   */
  virtual LocalOrdinal operator() (LocalOrdinal MyRow) const = 0;

  //! Returns the local overlapping partition ID of the j-th node in partition i.
  virtual LocalOrdinal operator() (LocalOrdinal i, LocalOrdinal j) const = 0;

  //! Returns the number of rows contained in specified partition.
  virtual size_t numRowsInPart(const LocalOrdinal Part) const = 0;
    
  //! Copies into List the rows in the (overlapping) partition Part.
  virtual void rowsInPart(const LocalOrdinal Part, Teuchos::ArrayRCP<LocalOrdinal> &List) const = 0;
  
  //! Returns an ArrayRCP to the integer vector containing the non-overlapping partition ID of each local row.
  virtual Teuchos::ArrayView<const LocalOrdinal>  nonOverlappingPartition() const = 0;

  //! Sets all the parameters for the partitioner.
  virtual void setParameters(Teuchos::ParameterList& List) = 0;

  //! Computes the partitions. Returns 0 if successful.
  virtual void compute() = 0;

  //! Returns true if partitions have been computed successfully.
  virtual bool isComputed() const = 0;

  //! Prints basic information about the partitioning object.
  virtual std::ostream& print(std::ostream& os) const = 0;

}; // class Ifpack2::Partitioner

template <class GraphType>
inline std::ostream& operator<<(std::ostream& os, const Ifpack2::Partitioner<GraphType>& obj)
{
  return(obj.print(os));
}

} //namespace Ipack2

#endif // IFPACK2_PARTITIONER_HPP
