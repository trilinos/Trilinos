#ifndef AMESOS_GRAPHPARTITIONER_H
#define AMESOS_GRAPHPARTITIONER_H

#include <iostream>
#include "Amesos_ConfigDefs.h"
#include <vector>

class Epetra_Comm;
#include "Epetra_CrsGraph.h"
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;
#include "Teuchos_ParameterList.hpp"

class Amesos_Partitioner {

public:

  Amesos_Partitioner(const Epetra_CrsGraph* CrsGraph,
		     const Epetra_CrsGraph* OverlappingGraph_,
		     Teuchos::ParameterList& List);

  ~Amesos_Partitioner();

  int NumLocalParts() const 
  {
    return(NumLocalParts_);
  }


  int OverlappingLevel() const 
  {
    return(OverlappingLevel_);
  }

  //! Returns the non-overlapping partition ID of row MyRow
  /*! Returns the non-overlapping partition ID of row MyRow. MyRow
   *  must be a local row of the non-overlapped graph. 
   */
  int operator() (int MyRow) const
  {
    if ((MyRow < 0) || (MyRow > NumMyRows()))
      AMESOS_CHK_ERR(-1); // input value not valid

    return(Partition_[MyRow]);
  }

  int NumSingletons() const
  {
    return(NumSingletons_);
  }

  const int* SingletonList() const
  {
    return(0);
  }
  
  int operator() (int Part, int Row) const
  {
    if ((Part < 0) || (Part >= NumLocalParts()))
      AMESOS_CHK_ERR(-1);

    if ((Row < 0) || (Row > Parts_[Part].size()))
      AMESOS_CHK_ERR(-2);

    return(Parts_[Part][Row]);
  }
    
  int operator() (int Part, int Row, bool LocalIndices) const
  {
    if (LocalIndices == false)
      return((*this)(Part,Row));

    int GID = (*this)(Part,Row);
    int LID = OverlappingGraph()->LRID(GID);
    return(LID);
  }

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

  const Epetra_CrsGraph* OverlappingGraph() const
  {
    return(OverlappingGraph_);
  }

  const int NumMyRows() const
  {
    return(Graph_->NumMyRows());
  }

  const int NumMyOverlappingRows() const
  {
    return(OverlappingGraph()->NumMyRows());
  }

  const int NumGlobalRows() const
  {
    return(Graph_->NumGlobalRows());
  }

  int Compute();

private:
   
  int SetParameters();

  int Compute1DPartition();

  int Compute2DPartition();

  int ComputeGreedyPartition();

  int ComputeMETISPartition();

  int ComputeOverlappingGraph();

  int ComputeOverlappingPartitions();

  int MaxNumEntries() const
  {
    return(Graph_->MaxNumIndices());
  }
    
  int ExtractMyRowCopy(int MyRow, int Length, int& NumIndices,
		       int* Indices)
  {
    return(Graph_->ExtractMyRowCopy(MyRow, Length, NumIndices,
				    Indices));
  }

  const Epetra_Comm& Comm() const
  {
    return(Graph_->Comm());
  }

  const int OverlappingGID(const int i) const
  {
    return(OverlappingGraph_->RowMap().GID(i));
  }

  // =========== //
  // object data //
  // =========== //

  //! Number of local subgraphs
  int NumLocalParts_;
 
  //! Partition_[i] contains the ID of non-overlapping part it belongs to
  vector<int> Partition_; 

  //! Parts_[i][j] is the ID of the j-th row contained in the (overlapping) 
  // partition i
  vector<vector<int> > Parts_;
  
  //! Reference to the graph to be partitioned
  const Epetra_CrsGraph* Graph_;

  //! Overlapping graph. It is computed only it OverlappingLevel_ > 1
  const Epetra_CrsGraph* OverlappingGraph_;
  
  /*! Map for the overlapping graph. It is computed only if OverlappingGraph_
   * > 1 */
  Epetra_BlockMap* OverlappingMap_;


  //! Overlapping level.
  int OverlappingLevel_;

  //! true if the partition has been computed with no errors
  bool IsPartitionComputed_;

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

  // singletons
  int NumSingletons_;
  vector<int> Singletons_;

  Teuchos::ParameterList List_;

  // required by certain partitioning schemes
  int nx_, ny_;
  int mx_, my_;
  int RootNode_;

}; // class GraphDecomposition

#endif
