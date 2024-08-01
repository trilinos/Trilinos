// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

  \section Ifpack2_Partitioner_Summary Summary

  Most Ifpack2 users will not need to create or use a Partitioner
  instance, or write a Partitioner subclass.  However, those
  implementing their own block relaxation algorithms may need to
  interact with Partitioner or write their own subclass thereof.
  Ifpack2's main application of this class is in BlockRelaxation.
  Partitioner defines the diagonal blocks of the matrix that
  BlockRelaxation uses.  BlockRelaxation creates a Partitioner
  subclass instance internally.
 
  \section Ifpack2_Partitioner_Partitions Partitions

  A Partitioner instance can partition a local graph by rows.  A
  <i>local</i> graph is one for which, on every process, the column
  Map contains no entries not in the domain Map on that process.  You
  may use LocalFilter on the graph of a matrix to make a local graph;
  it excludes entries in the column Map not in the domain Map.  This
  class assumes that the graph is local.
  
  The partitions created by Partitioner implementations are
  <i>nonoverlapping</i> in the graph sense. This means that each row
  (or, more appropriately, vertex) of the graph belongs to at most one
  partition.  Furthermore, these nonoverlapping partitions are
  <i>local</i>: partitions do not cross process boundaries.
  
  <tt>operator () (LocalOrdinal i)</tt> returns the local partition
  index corresponding to local row i of the graph.  The above implies
  that the local partition index is unique.

  \section Ifpack2_Partitioner_OverlappingPartitions Overlapping decomposition

  The OverlappingPartitioner subclass extends the nonoverlapping
  partitions by the required amount of overlap, considering local
  vertices only.  That is, this overlap does <i>not</i> modify the
  overlap among the processes.  (The mathematical definition of
  "partition" does not allow overlap, but we accept this as a useful
  extension.)

  \section Ifpack2_Partitioner_Subclasses Subclasses of Partitioner

  Partitioner is just an interface; it does not implement any
  fucntionality.  You cannot create an instance of Partitioner; you
  must instantiate a concrete implementation thereof.  Concrete
  implementations include:
    - LinearPartitioner, which partitions the graph into contiguous
      row blocks
    - Zoltan2Partitioner, which calls Zoltan2 to partition the graph

  The constructor takes a Tpetra::RowGraph instance, which is the
  graph to partition.

  \section Ifpack2_Partitioner_Example Example code

  The following example code shows how to use LinearPartitioner, a
  subclass of Partitioner.  Note that most Ifpack2 users will
  <i>not</i> need to create or interact with Partitioner instances.
  We only show this example for the sake of developers who might need
  to implement or use Preconditioner subclasses, and want an example
  of a Partitioner "in action."

  \code
#include "Ifpack2_LinearPartitioner.hpp"
#include "Tpetra_CrsMatrix.hpp"
// ...

// The matrix A must be fill complete.
void example (Tpetra::CrsMatrix<double, int>& A) {
  using std::cout;
  using std::endl;
  typedef Tpetra::CrsGraph<int> graph_type;
  typedef Ifpack2::LinearPartitioner<graph_type> partitioner_type;

  // Create the partitioner.
  partitioner_type partitioner (A.getGraph ());

  // Set up the partitioner's parameters.
  // We want 16 local partitions, 
  // and an overlap of 0 among the local partitions.
  Teuchos::ParameterList params;
  params.set ("partitioner: local parts", 16);
  params.set ("partitioner: overlap", 0);
  partitioner.setParameters (params);

  // Partition the graph.  If the structure of the 
  // graph changes, you must call compute() again,
  // but you need not call setParameters() again.
  partitioner.compute ();

  // Get the number of partitions created on the calling process.
  const int numParts = partitioner.numLocalParts ();

  // Get the number of rows in each partition.
  for (int i = 0; i < numParts; ++i) {
    cout << "Number of rows in partition " << i << ": " 
         << partitioner.numRowsInPart (i) << endl;
  }  

  // For nonoverlapping partitions only, operator()(i)
  // returns the partition index for each local row.
  const size_t numLocalRows = A.getLocalNumRows ();
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
  for (size_t j = 0; j < partitioner.numRowsInPart (i); ++j) {
    cout << partitioner(i,j) << " ";
  }
  cout << "]";
}  
\endcode
*/  
template <class GraphType>
class Partitioner : public Teuchos::Describable {
public:
  typedef typename GraphType::local_ordinal_type LocalOrdinal;
  typedef typename GraphType::global_ordinal_type GlobalOrdinal;
  typedef typename GraphType::node_type Node;

  //! Destructor.
  virtual ~Partitioner() {};

  /// \brief Number of computed local partitions.
  ///
  /// See Ifpack2_OverlappingPartitioner_decl.hpp for explanation
  /// of why this is an \c int instead of \c LocalOrdinal.
  virtual int numLocalParts () const = 0;

  //! The level of overlap.
  virtual int overlappingLevel() const = 0;

  /// \brief The local (nonoverlapping) partition index of the
  ///   specified local row.
  ///
  /// \param MyRow [in] Local index of the row.
  virtual LocalOrdinal operator() (LocalOrdinal MyRow) const = 0;

  //! The local overlapping partition index of the j-th node in partition i.
  virtual LocalOrdinal operator() (LocalOrdinal i, LocalOrdinal j) const = 0;

  //! The number of rows contained in the specified partition.
  virtual size_t numRowsInPart (const LocalOrdinal Part) const = 0;
    
  //! Copy into List the rows in the (overlapping) partition Part.
  virtual void 
  rowsInPart (const LocalOrdinal Part, 
              Teuchos::ArrayRCP<LocalOrdinal>& List) const = 0;
  
  //! The nonoverlapping partition indices of each local row.
  virtual Teuchos::ArrayView<const LocalOrdinal>
  nonOverlappingPartition () const = 0;

  //! Set all the parameters for the partitioner.
  virtual void setParameters (Teuchos::ParameterList& List) = 0;

  //! Compute the partitions.
  virtual void compute () = 0;

  //! Return true if partitions have been computed successfully.
  virtual bool isComputed () const = 0;

  //! Print basic information about the partitioning object.
  virtual std::ostream& print (std::ostream& os) const = 0;

};

// Overloaded output stream operator for Partitioner
template <class GraphType>
inline std::ostream& 
operator<< (std::ostream& os, 
            const Ifpack2::Partitioner<GraphType>& obj)
{
  return obj.print (os);
}

} // namespace Ifpack2

#endif // IFPACK2_PARTITIONER_HPP
