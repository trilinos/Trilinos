/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_PARTITIONER_HPP
#define IFPACK2_PARTITIONER_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include <iostream>

namespace Ifpack2 {

//! Ifpack2::Partitioner: 

/** \class Partitioner
    \brief A class to decompose local graphs.
 
  Implementations of this class can partition a local graph by rows.
  This class assumes that the graph is local.  A <i>local</i> graph is
  one for which, on every process, the column Map contains no entries
  not in the domain Map on that process.  You may use LocalFilter on
  the graph of a matrix to make a local graph; it excludes entries in
  the column Map not in the domain Map.
  
  <tt>operator () (LocalOrdinal i)</tt> returns the local partition
  index corresponding to local row i of the graph.
  
  The partitions created by Partitioner implementations are
  nonoverlapping in the graph sense. This means that each row (or,
  more approriately, vertex) of \c G is assigned to exactly one
  partition.  Partitioner can be extended using the functionalities of
  OverlappingPartitioner (itself derived from Partitioner).  This
  class extends the nonoverlapping partitions by the required amount
  of overlap, considering local nodes only (that is, this overlap do
  <i>not</i> modify the overlap among the processes).

  Partitioner is a pure virtual class. Concrete implementations include:
    - LinearPartitioner, which allows the decomposition of the
      rows of the graph in simple consecutive chunks;
    - Zoltan2Partitioner, which calls Zoltan2 to decompose the graph

  The constructor requires a Tpetra::RowGraph instance, representing
  the graph to decompose.
  
  <P>An example use of a Partitioner derived class is as follows:  
  \code
#include "Ifpack2_LinearPartitioner.hpp"
#include "Tpetra_CrsMatrix.hpp"
// ...

// The matrix A must be fill complete.
void example (Tpetra::CrsMatrix& A) {
  using std::cout;
  using std::endl;

  Ifpack2::LinearPartitioner partitioner (A.getGraph ());

  // Set up the partitioner's parameters.
  // We want 16 local partitions, 
  // and an overlap of 0 among the local partitions.
  Teuchos::ParameterList params;
  params.set ("partitioner: local parts", 16);
  params.set ("partitioner: overlap", 0);
  partitioner.setParameters (params);

  // Partition the graph.
  partitioner.compute ();

  // We can get the number of partitions created on the calling process...
  const int numParts = partitioner.numLocalParts ();

  // ... and the number of rows in each of them.
  for (int i = 0; i < numParts; ++i) {
    cout << "Number of rows in partition " << i << ": " 
         << partitioner.numRowsInPart (i) << endl;
  }  

  // For nonoverlapping partitions only, we can get
  // the partition index for each local row, by using
  // the parentheses operator:
  const size_t numLocalRows = A.getNodeNumRows ();
  for (size_t i = 0; i < numLocalRows; ++i) {
    cout << "Partition[" << i <<"] = " << partitioner(i) << endl;
  }
}
\endcode
  
When overlapping partitions are created, users can get the local
indices of the rows in each partition:
\code
const int numLocalParts = partitioner.numLocalParts ();
for (int i = 0; i < numLocalParts; ++i) {
  cout << "Local rows of Partition " << i << ": [";
  for (size_t j = 0; j < Partitioner.numRowsInPart (i); ++j) {
    cout << Partitioner(i,j) << " ";
  }
  cout << "]";
}  
\endcode
  
Partitioner is used to create the subblocks in BlockRelaxation.
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
  // See Ifpack2_OverlappingPartitioner_decl.hpp for explanation
  // of why this is an "int" instead of "LocalOrdinal"
  virtual int numLocalParts() const = 0;

  //! Returns the overlapping level.
  virtual int overlappingLevel() const = 0;

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
