#ifndef IFPACK_OVERLAPPINGPARTITIONER_H
#define IFPACK_OVERLAPPINGPARTITIONER_H

#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Ifpack_Partitioner.h"
#include "Teuchos_ParameterList.hpp"
#include <vector>
class Epetra_Comm;
class Ifpack_Graph;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;

//! Ifpack_OverlappingPartitioner: A class to decompose overlapping and non-overlapping Ifpack_Graph's.

/*!

  \note Partition ID's are local with respect to the numbering of \c G
  and \c OG.

  \author Marzio Sala, 9214

  \date Sep-04
*/  
class Ifpack_OverlappingPartitioner : public Ifpack_Partitioner {

public:

  Ifpack_OverlappingPartitioner(const Ifpack_Graph* Graph);

  ~Ifpack_OverlappingPartitioner();

  //! Returns the number of computed local partitions.
  int NumLocalParts() const 
  {
    return(NumLocalParts_);
  }

  //! Returns the overlapping level.
  int OverlappingLevel() const 
  {
    return(OverlappingLevel_);
  }

  //! Returns the local non-overlapping partition ID of the specified row.
  /*! Returns the non-overlapping partition ID of the specified row.
   \param In
   MyRow - local row numbe

   \return
   Local ID of non-overlapping partition for \t MyRow.
   */
  int operator() (int MyRow) const
  {
    if ((MyRow < 0) || (MyRow > NumMyRows()))
      IFPACK_CHK_ERR(-1); // input value not valid

    return(Partition_[MyRow]);
  }

  //! Returns the number of singletons in Graph.
  int NumSingletons() const
  {
    return(NumMyRows() - NumMyNonDirichletRows());
  }

  //! Returns a pointer to the internally-stored list of singletons (TO DO).
  const int* SingletonList() const
  {
    return(0);
  }
  
  //! Returns the local overlapping partition ID of the j-th node in partition i.
  int operator() (int i, int j) const
  {
    if ((i < 0) || (i >= NumLocalParts()))
      IFPACK_CHK_ERR(-1);

    if ((j < 0) || (j > Parts_[i].size()))
      IFPACK_CHK_ERR(-2);

    return(Parts_[i][j]);
  }

  //! Returns the number of rows contained in specified partition.
  int NumRowsInPart(const int Part) const
  {
    return(Parts_[Part].size());
  }
    
  int RowsInPart(const int Part, int* List) const
  {
    for (int i = 0 ; i < NumRowsInPart(Part) ; ++i)
      List[i] = Parts_[Part][i];

    return(0);
  }
  
  const int* NonOverlappingPartition() const
  {
    return(&Partition_[0]);
  }

  //! Sets all the parameters for the partitioner.
  virtual int SetParameters(Teuchos::ParameterList& List);

  //! Sets all the parameters for the partitioner.
  virtual int SetPartitionParameters(Teuchos::ParameterList& List) = 0;

  //! Computes the partitions. Returns 0 if successful.
  virtual int Compute();

  //! Computes the partitions. Returns 0 if successful.
  virtual int ComputePartitions() = 0;

  //! Computes the partitions. Returns 0 if successful.
  virtual int ComputeOverlappingPartitions();
  
  //! Returns true if partitions have been computed successfully.
  bool IsComputed()
  {
    return(IsComputed_);
  }

protected:
   
  const int NumMyRows() const;

  const int NumMyNonzeros() const;

  const int NumMyNonDirichletRows() const;

  const int NumGlobalRows() const;

  int MaxNumEntries() const;
    
  // FIXME
  int ExtractMyRowCopy(int MyRow, int Length, int& NumIndices,
		       int* Indices);

  const Epetra_Comm& Comm() const;

  //! Number of local subgraphs
  int NumLocalParts_;
 
  //! Partition_[i] contains the ID of non-overlapping part it belongs to
  vector<int> Partition_; 

  //! Parts_[i][j] is the ID of the j-th row contained in the (overlapping) 
  // partition i
  vector<vector<int> > Parts_;
  
  //! Reference to the graph to be partitioned
  const Ifpack_Graph* Graph_;

  //! Overlapping level.
  int OverlappingLevel_;

  //! true if the partition has been computed with no errors
  bool IsComputed_;

  //! Sets the output level.
  /*! Sets the output level as follows:
   *  - 0 no output
   *  - 1 warnings
   *  - 2 normal output
   *  - 3 print information about input parameter
   *  - 5 print timing
   *  - 6 print memory allocation information
   *  - 9 verbose
   *  - 10 verbose on all processors
   */
  int verbose_;
  string PrintMsg_;
  string ErrorMsg_;

  //! Type of decomposition
  string DecompositionType_;

  // Dirichlet rows
  int NumMyNonDirichletRows_;
  vector<int> Mask_;

}; // class Ifpack_Partitioner

#endif // HAVE_IFPACK_TEUCHOS
#endif // IFPACK_OVERLAPPINGPARTITIONER_H
