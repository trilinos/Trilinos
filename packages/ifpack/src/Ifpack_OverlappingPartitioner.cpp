#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Ifpack_Partitioner.h"
#include "Ifpack_OverlappingPartitioner.h"
#include "Ifpack_Graph.h"

#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Teuchos_ParameterList.hpp"

//==============================================================================
Ifpack_OverlappingPartitioner::
Ifpack_OverlappingPartitioner(const Ifpack_Graph* Graph) :
  NumLocalParts_(0),
  Graph_(Graph),
  OverlappingLevel_(0),
  IsComputed_(false),
  verbose_(2),
  PrintMsg_("(Ifpack_Partitioner) "),
  ErrorMsg_("Ifpack_Partitioner ERROR "),
  NumMyNonDirichletRows_(0)
{
}

//==============================================================================
Ifpack_OverlappingPartitioner::~Ifpack_OverlappingPartitioner()
{
}

//==============================================================================
int Ifpack_OverlappingPartitioner::SetParameters(Teuchos::ParameterList& List)
{

  NumLocalParts_ = List.get("partitioner: local parts", NumLocalParts_);
  OverlappingLevel_ = List.get("partitioner: overlap", OverlappingLevel_);
  verbose_ = List.get("partitioner: print level", 2);

  SetPartitionParameters(List);

  return(0);
}

//==============================================================================
int Ifpack_OverlappingPartitioner::Compute()
{

  if (NumLocalParts_ < 1)
    IFPACK_CHK_ERR(-1); // incorrect value

  if (OverlappingLevel_ < 0)
    IFPACK_CHK_ERR(-1); // incorrect value

  // some output

  if ((verbose_ > 2) && (Comm().MyPID() == 0)) {
    cout << PrintMsg_ << "Number of local parts  = " << NumLocalParts_ << endl;
    cout << PrintMsg_ << "Number of global parts = " 
         << NumLocalParts_ * Comm().NumProc() << endl;
    cout << PrintMsg_ << "Amoung of overlap      = " << OverlappingLevel_ << endl;
  }
  // 1.- allocate memory 

  Partition_.resize(NumMyRows());
  Parts_.resize(NumLocalParts());

  // 2.- sanity checks on input graph
 
  if (Graph_->Filled() == false)
    IFPACK_CHK_ERR(-4); // need FillComplete() called

  if (Graph_->NumGlobalRows() != Graph_->NumGlobalCols())
    IFPACK_CHK_ERR(-3); // can partition square matrices only

  if (NumLocalParts_ < 1)
    IFPACK_CHK_ERR(-2); // value not valid
 
  // 3.- localize Dirichlet nodes. Mask_[i] will contain either
  //     -1 (node is Dirichlet), or the local ID (without
  //     counting the Dirichlet nodes).

  Mask_.resize(NumMyRows());

  NumMyNonDirichletRows_ = 0;
  int NumIndices;
  vector <int> Indices;
  Indices.resize(MaxNumEntries());

  for (int i = 0; i < NumMyRows() ; ++i) {  

    int ierr;
    ierr = Graph_->ExtractMyRowCopy(i, MaxNumEntries(), NumIndices, &Indices[0]);
    IFPACK_CHK_ERR(ierr);
    
    if (NumIndices <= 1) {
      Mask_[i] = -1;
    } else {
      Mask_[i] = NumMyNonDirichletRows_++;
    }
  }

  // 2.- perform non-overlapping partition
 
  IFPACK_CHK_ERR(ComputePartitions());

  // 3.- compute the partitions with overlapping
  
  IFPACK_CHK_ERR(ComputeOverlappingPartitions());

  // 4.- return to the user
 
  IsComputed_ = true;

  return(0);
}
int Ifpack_OverlappingPartitioner::ComputeOverlappingPartitions()
{

  // FIXME: the first part of this function should be elsewhere
  // start defining the subgraphs for no overlap

  vector<int> sizes;
  sizes.resize(NumLocalParts_);

  // 1.- compute how many rows are in each subgraph
  for (int i = 0 ; i < NumLocalParts_ ; ++i)
    sizes[i] = 0;

  for (int i = 0 ; i < NumMyRows() ; ++i) {
    if (Partition_[i] >= NumLocalParts_) {
      cout << "Partition[" << i << "] = "<< Partition_[i] 
	   << ", NumLocalParts = " << NumLocalParts_ << endl;
      IFPACK_CHK_ERR(-10);
    }
    // discard singletons
    if (Partition_[i] != -1)
      sizes[Partition_[i]]++;
  }

  // 2.- allocate space for each subgraph

  for (int i = 0 ; i < NumLocalParts_ ; ++i)
    Parts_[i].resize(sizes[i]);

  // 3.- cycle over all rows and populate the vectors

  for (int i = 0 ; i < NumLocalParts_ ; ++i)
    sizes[i] = 0;

  for (int i = 0 ; i < NumMyRows() ; ++i) {
    int part = Partition_[i];
    // discard singletons
    if (part == -1)
      continue;
    int count = sizes[part];
    Parts_[part][count] = i;
    sizes[part]++;
  }

  if (OverlappingLevel_ == 0)
    return(0);

  // wider overlap requires further computations
 
  for (int level = 1 ; level <= OverlappingLevel_ ; ++level) {

    vector<vector<int> > tmp;
    tmp.resize(NumLocalParts_);

    // cycle over all rows in the local graph (that is the overlapping
    // graph). For each row, all columns will belong to the subgraph of
    // row `i'.

    int MaxNumEntries = Graph_->MaxMyNumEntries();
    vector<int> Indices;
    Indices.resize(MaxNumEntries);

    for (int part = 0 ; part < NumLocalParts_ ; ++part) {

      for (int i = 0; i < Parts_[part].size() ; ++i) {  

	int LRID = Parts_[part][i];
	int NumIndices;
	int ierr = Graph_->ExtractMyRowCopy(LRID, MaxNumEntries, 
							NumIndices, &Indices[0]);
	IFPACK_CHK_ERR(ierr);

	for (int j = 0 ; j < NumIndices ; ++j) {

	  // use *local* indices
	  int col = Indices[j];
	  // has this column already been inserted?
	  vector<int>::iterator
	    where = find(tmp[part].begin(), tmp[part].end(), col);

	  if (where == tmp[part].end()) {
	    tmp[part].push_back(col);
	  }

	}
      }

    }

    // now I convert the STL vectors into Epetra_IntSerialDenseVectors.

    for (int i = 0 ; i < NumLocalParts_ ; ++i) {
      Parts_[i].resize(tmp[i].size());
      for (int j = 0 ; j < tmp[i].size() ; ++j)
	Parts_[i][j] = tmp[i][j];
    }

  }

  return(0);

}

//============================================================================
const int Ifpack_OverlappingPartitioner::NumMyRows() const
{
  return(Graph_->NumMyRows());
}

//============================================================================
const int Ifpack_OverlappingPartitioner::NumMyNonzeros() const
{
  return(Graph_->NumMyNonzeros());
}

//============================================================================
const int Ifpack_OverlappingPartitioner::NumMyNonDirichletRows() const
{
  return(NumMyNonDirichletRows_);
}

//============================================================================
const int Ifpack_OverlappingPartitioner::NumGlobalRows() const
{
  return(Graph_->NumGlobalRows());
}

//============================================================================
int Ifpack_OverlappingPartitioner::MaxNumEntries() const
{
  return(Graph_->MaxMyNumEntries());
}

//============================================================================
const Epetra_Comm& Ifpack_OverlappingPartitioner::Comm() const
{
  return(Graph_->Comm());
}
#endif // HAVE_IFPACK_TEUCHOS
