//@HEADER
// ***********************************************************************
// 
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
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

#ifndef TIFPACK_PARTITIONER_HPP
#define TIFPACK_PARTITIONER_HPP

#include "Tifpack_ConfigDefs.hpp"
#include "Teuchos_ParameterList.hpp"
class Tpetra_Comm;
class Tifpack_Graph;
class Tpetra_Map;
class Tpetra_BlockMap;
class Tpetra_Import;

//! Tifpack_Partitioner: A class to decompose local Tifpack_Graph's.

/*!
 
  Class Tifpack_Partitioner enables the decomposition of a local
  Tifpack_Graph's. It is supposed that the graph refers to
  a localized matrix (that is, a matrix that has been filtered
  through Tifpack_LocalFilter).
  
  The overloaded operator (int i) can be used to extract the local partition
  ID of local row i.  
  
  The partitions created by Tifpack_Partitioner derived clased 
  are non-overlapping in graph sense. This means that each row
  (or, more approriately, vertex)
  of \c G is assigned to exactly one partition.

  Partitioner can be extended using the functionalities of class
  Tifpack_OverlappingPartitioner (itself derived from Tifpack_Partitioner.
  This class extends the non-overlapping partitions by the required
  amount of overlap, considering local nodes only (that is, this
  overlap do \e not modify the overlap among the processes).

  Tifpack_Partitioner is a pure virtual class. Concrete implementations
  are:
  - Tifpack_LinearPartitioner, which allows the decomposition of the
    rows of the graph in simple consecutive chunks;
  - Tifpack_METISPartitioner, which calls METIS to decompose the graph
    (this requires the configuration option --enable-ifpack-metis);
  - Tifpack_GreedyPartitioner, a simple greedy algorith;
  - Tifpack_EquationPartitioner, which creates \c NumPDEEqns parts
    (where \c NumPDEEqns is the number of equations in the linear
    system). It is supposed that all the equations referring to the 
    same grid node are ordered consecutively. Besides, the 
    number of equations per node must be constant in the domain.

  Generically, a constructor requires an Tifpack_Graph object. 
  Tifpack_Graph is a pure virtual class. Concrete implentations are:
  - Tifpack_Graph_Tpetra_CrsGraph, a light-weight class to wrap 
    Tpetra_CrsGraph objects as Tifpack_Graph objects;
  - Tifpack_Graph_Tpetra_RowMatrix, a light-weight class to
    wrap Tpetra_RowMatrix objects as Tifpack_Graph objects.
  
  <P>An example of use is an Tifpack_Partitioner derived class is as follows:  
  \code
#include "Tifpack_Partitioner.hpp"
#include "Tifpack_LinearPartitioner.hpp"
#include "Tifpack_Graph.hpp"
#include "Tifpack_Graph_Tpetra_CrsGraph.hpp"
...
Tpetra_CrsMatrix* A;         // A is filled
// create the wrapper from Tpetra_CrsGraph
Tifpack_Graph* Graph = new Tifpack_Graph_Tpetra_CrsGraph(A);

// we aim to create non-overlapping partitions only
Tifpack_Partitioner Partitioner(Graph);

Tifpack_Partitioner* Partitioner;
Partitioner = new Tifpack_Graph_Tpetra_CrsGraph(&A);

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
  
Tifpack_Partitioner is used to create the subblocks in 
Tifpack_BlockJacobi, Tifpack_BlockGaussSeidel, and 
Tifpack_BlockSymGaussSeidel.

\author Michael Heroux, SNL 9214.

\date Last modified on Nov-04.

*/  
class Tifpack_Partitioner {

public:

  //! Destructor.
  virtual ~Tifpack_Partitioner() {};

  //! Returns the number of computed local partitions.
  virtual int NumLocalParts() const = 0;

  //! Returns the overlapping level.
  virtual int OverlappingLevel() const = 0;

  //! Returns the local non-overlapping partition ID of the specified row.
  /*! Returns the non-overlapping partition ID of the specified row.
   \param 
   MyRow - (In) local row numbe

   \return
   Local ID of non-overlapping partition for \t MyRow.
   */
  virtual int operator() (int MyRow) const = 0;

  //! Returns the local overlapping partition ID of the j-th node in partition i.
  virtual int operator() (int i, int j) const = 0;

  //! Returns the number of rows contained in specified partition.
  virtual int NumRowsInPart(const int Part) const = 0;
    
  //! Copies into List the rows in the (overlapping) partition Part.
  virtual int RowsInPart(const int Part, int* List) const = 0;
  
  //! Returns a pointer to the integer vector containing the non-overlapping partition ID of each local row.
  virtual const int* NonOverlappingPartition() const = 0;

  //! Sets all the parameters for the partitioner.
  virtual int SetParameters(Teuchos::ParameterList& List) = 0;

  //! Computes the partitions. Returns 0 if successful.
  virtual int Compute() = 0;

  //! Returns true if partitions have been computed successfully.
  virtual bool IsComputed() = 0;

  //! Prints basic information about the partitioning object.
  virtual ostream& Print(std::ostream& os) const = 0;

}; // class Tifpack_Partitioner

inline ostream& operator<<(ostream& os, const Tifpack_Partitioner& obj)
{
  return(obj.Print(os));
}

#endif // TIFPACK_PARTITIONER_HPP
