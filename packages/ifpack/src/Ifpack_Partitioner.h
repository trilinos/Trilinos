#ifndef IFPACK_PARTITIONER_H
#define IFPACK_PARTITIONER_H

#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Teuchos_ParameterList.hpp"
#include <vector>
class Epetra_Comm;
class Ifpack_Graph;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;

//! Ifpack_Partitioner: A class to decompose overlapping and non-overlapping Ifpack_Graph's.

/*!
 
  FIXME:

  Class Ifpack_Partitioner enables the decomposition of Ifpack_Graph's.
  The constructor of Ifpack_Partitioner requires two graphs:
  - a starting graph (G);
  - an overlapping graph (OG), required only for the construction of
    overlapping partitions.
  
  Ifpack_Partitioner operates as follows:
  - \c G will be decomposed into the specified number of parts. 
    Only the local subgraph of G is partitioned. The overloaded operator
    (int i) can be used to extract the local partition ID of local row i.
    Partitions are non-overlapping (in graph sense), as each row of \c G
    is included in exactly one partition.
  - If required, the non-overlapping partitions can be extended, so that
    they share an overlap (\c OverlapLevel). Overlapping partitions
    are created from \c OG. It is supposed that OG has been created 
    by overlapping G of \c OverlapLevel.

  Ifpack_Partitioner requires two Ifpack_Graph objects. Concrete classes
  are provided, to create Ifpack_Graph's from Epetra_CrsGraph's, and
  Epetra_CrsMatrix's.
  
  An example of use is a follows:  
  \code
#include "Ifpack_Partitioner.h"
#include "Ifpack_Graph.h"
#include "Ifpack_GraphEpetraCrs.h"
...
Epetra_CrsMatrix* A;         // A is filled
Ifpack_Graph Graph = new Ifpack_GraphEpetraCrs(A);

// we aim to create non-overlapping partitions only
Ifpack_Partitioner Partitioner(Graph,(Ifpack_Graph*)0);

Teuchos::ParameterList List;
List.set("partitioner", "linear");     // linear decomposition
List.set("overlap level", 0);          // no overlap

Partitioner.Create(List);              // decompose the graph;
delete Graph;                          // now Graph can be delete

int NumParts = Partitioner.NumParts(); // parts actually created

for (int i = 0 ; i < NumParts ; ++i) { // gets the rows in each part
  cout << "rows in " << i << "=" << Partitioner.RowsInPart(i);
}  

                                       // get the partition ID for 
				       // each local row
for (int i = 0 ; i < A->NumMyRows() ; ++i)
  cout << "Partition[" << i <<"] = " << Partitioner(i) << endl;

  \endcode
  
  The construction of overlapping partitions requires the overlapping 
  graph as well. To overlap Epetra_CrsMatrix's, function
  Ifpack_CreateOverlappingCrsMatrix() can be used. 
  To overlap Epetra_CrsGraph's, function Ifpack_CreateOverlappingCrsGraph() 
  can be employed.
  See the following code for an example.

  \code
#include "Ifpack_Partitioner.h"
#include "Ifpack_Graph.h"
#include "Ifpack_GraphEpetraCrs.h"
#include "Ifpack_Utils.h"
...
Epetra_CrsMatrix* A;         // A is filled
                             // create the overlapping matrix

int OverlapLevel = 3;			     
Ifpack_Graph OverlappingGraph = 
  new Ifpack_CreateOverlappingCrsMatrix(A, OverlapLevel)

                             // create graphs
Ifpack_Graph Graph = new Ifpack_GraphEpetraCrs(A);
Ifpack_Graph OverlappingGraph = new Ifpack_GraphEpetraCrs(OverlappingA);

// we aim to create non-overlapping partitions only
Ifpack_Partitioner Partitioner(Graph,OverlappingGraph);

Teuchos::ParameterList List;
List.set("partitioner", "linear");     // linear decomposition
List.set("overlap level", OverlapLevel);

Partitioner.Create(List);              // decompose the graph;
delete Graph;                          // now Graph can be delete

int NumParts = Partitioner.NumParts(); // parts actually created

                                       // get the rows contained
				       // in each partition
for (int i = 0 ; i < NumParts ; ++i) {
  for (int j = 0 ; j < Partitioner.RowsInPart(i) ; ++j) {
    cout << "Partition " << i << ", contains local row "
         << Partitioner(i,j) << endl;
  }
}  
  \endcode
  
  Ifpack_Partitioner is used to create the subblocks in Ifpack_BlockJacobi,
  Ifpack_OverlappingBlockJacobi, Ifpack_GaussSeidel, 
  Ifpack_OverlappingGaussSeidel. 

  \note Partition ID's are local with respect to the numbering of \c G
  and \c OG.

  \author Marzio Sala, 9214

  \date Sep-04
*/  
class Ifpack_Partitioner {

public:

  ~Ifpack_Partitioner() {};

  //! Returns the number of computed local partitions.
  virtual int NumLocalParts() const = 0;

  //! Returns the overlapping level.
  virtual int OverlappingLevel() const = 0;

  //! Returns the local non-overlapping partition ID of the specified row.
  /*! Returns the non-overlapping partition ID of the specified row.
   \param In
   MyRow - local row numbe

   \return
   Local ID of non-overlapping partition for \t MyRow.
   */
  virtual int operator() (int MyRow) const = 0;

  //! Returns the number of singletons in Graph.
  virtual int NumSingletons() const = 0;

  //! Returns a pointer to the internally-stored list of singletons (TO DO).
  virtual const int* SingletonList() const = 0;
  
  //! Returns the local overlapping partition ID of the j-th node in partition i.
  virtual int operator() (int i, int j) const = 0;

  //! Returns the number of rows contained in specified partition.
  virtual int NumRowsInPart(const int Part) const = 0;
    
  virtual int RowsInPart(const int Part, int* List) const = 0;
  
  virtual const int* NonOverlappingPartition() const = 0;

  //! Sets all the parameters for the partitioner.
  virtual int SetParameters(Teuchos::ParameterList& List) = 0;

  //! Computes the partitions. Returns 0 if successful.
  virtual int Compute() = 0;

  //! Returns true if partitions have been computed successfully.
  virtual bool IsComputed() = 0;

}; // class Ifpack_Partitioner

#endif // HAVE_IFPACK_TEUCHOS
#endif // IFPACK_PARTITIONER_H
