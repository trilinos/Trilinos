/* FIXME:
 * 
 * - singletons for all decomposition schemes
 * - add 3D partitions ??
 * - for ExtractMyRow, change Copy into View?
 * - timing
 * - other output
 *
 */
#include "Amesos_ConfigDefs.h"
#include "Amesos_Partitioner.h"

#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Teuchos_ParameterList.hpp"

// FIXME....
#define HAVE_AMESOS_METIS
#ifdef HAVE_AMESOS_METIS
typedef int idxtype;
extern "C" {
  void METIS_EstimateMemory(int *, idxtype *, idxtype *, int *, int *, int *);
  void METIS_PartGraphKway(int *, idxtype *, idxtype *, idxtype *, 
			   idxtype *, int *, int *, int *, int *, int *,
			   idxtype *);
  void METIS_PartGraphRecursive(int *, idxtype *, idxtype *, 
				idxtype *, idxtype *, int *, int *, int *, 
				int *, int *, idxtype *);

}
#endif

//==============================================================================
Amesos_Partitioner::
Amesos_Partitioner(const Epetra_CrsGraph* Graph,
		   const Epetra_CrsGraph* OverlappingGraph_,
		   Teuchos::ParameterList& List) :
  NumLocalParts_(0),
  Graph_(Graph),
  OverlappingGraph_(OverlappingGraph_),
  OverlappingMap_(0),
  OverlappingLevel_(0),
  IsPartitionComputed_(false),
  verbose_(2),
  PrintMsg_("(Amesos_Partitioner) "),
  ErrorMsg_("Amesos_Partitioner ERROR "),
  DecompositionType_("not-set"),
  NumSingletons_(0),
  nx_(-1),
  ny_(-1),
  mx_(-1),
  my_(-1),
  RootNode_(0),
  List_(List)
{
  SetParameters();
}

//==============================================================================
int Amesos_Partitioner::SetParameters()
{

  // get parameters from input list
  // FIXME: use isParameter() 
  NumLocalParts_ = List_.get("local parts", (int)1);
  OverlappingLevel_ = List_.get("overlap level", (int)0);
  DecompositionType_ = List_.get("partitioner", "METIS");
  verbose_ = List_.get("output level", 2);
      
  nx_ = List_.get("nx", nx_);
  ny_ = List_.get("ny", ny_);
  mx_ = List_.get("mx", mx_);
  my_ = List_.get("my", my_);
  RootNode_ = List_.get("root node", RootNode_);
  
  // sanity checks
 
  if (NumLocalParts_ < 1)
    AMESOS_CHK_ERR(-1); // incorrect value

  if (OverlappingLevel_ < 0)
    AMESOS_CHK_ERR(-1); // incorrect value

  // some output

  if ((verbose_ > 2) && (Comm().MyPID() == 0)) {
    cout << PrintMsg_ << "Number of local parts  = " << NumLocalParts_ << endl;
    cout << PrintMsg_ << "Number of global parts = " 
         << NumLocalParts_ * Comm().NumProc() << endl;
    cout << PrintMsg_ << "Amoung of overlap      = " << OverlappingLevel_ << endl;
    cout << PrintMsg_ << "Decomposition type     = " 
         << DecompositionType_ << endl;
  }
    
  return(0);

}

//==============================================================================
Amesos_Partitioner::~Amesos_Partitioner()
{

  IsPartitionComputed_ = false;

  if (OverlappingMap_)
    delete OverlappingMap_;

}

//==============================================================================
int Amesos_Partitioner::Compute()
{

  // 1.- allocate memory 

  Partition_.resize(NumMyRows());
  Parts_.resize(NumLocalParts());

  // 2.- sanity checks on input graph
 
  if (Graph_->Filled() == false)
    AMESOS_CHK_ERR(-4); // need FillComplete() called

  if (Graph_->NumGlobalRows() != Graph_->NumGlobalCols())
    AMESOS_CHK_ERR(-3); // can partition square matrices only

  if (NumMyRows() < NumLocalParts_)
    AMESOS_CHK_ERR(-2); // not enough local rows for this job
 
  // 2.- perform non-overlapping partition
 
  if (DecompositionType_ == "METIS") {
    AMESOS_CHK_ERR(ComputeMETISPartition());
  }
  else if (DecompositionType_ == "1D" ) {
    AMESOS_CHK_ERR(Compute1DPartition());
  }
  else if (DecompositionType_ == "2D" ) {
    AMESOS_CHK_ERR(Compute2DPartition());
  }
  else if (DecompositionType_ == "greedy" ) {
    AMESOS_CHK_ERR(ComputeGreedyPartition());
  }
  else
    AMESOS_CHK_ERR(-1); // value not valid

  // 3.- compute the partitions with overlapping
  
  AMESOS_CHK_ERR(ComputeOverlappingPartitions());

  // 4.- return to the user
 
  IsPartitionComputed_ = true;

  return(0);
}

//==============================================================================
int Amesos_Partitioner::Compute1DPartition()
{

  int mod = NumMyRows() / NumLocalParts_;
  for (int i = 0 ; i < NumMyRows() ; ++i) {
    Partition_[i] = i / mod;
    if (Partition_[i] >= NumLocalParts_)
      Partition_[i] = NumLocalParts_ - 1;
  }

  return(0);

}

//==============================================================================
int Amesos_Partitioner::Compute2DPartition()
{

  // for square matrices, automatically compute the nodes
  // along the X-axis and the Y-axis

  if ((nx_ == -1) || (ny_ == -1)) {
    nx_ = (int)sqrt((double)NumMyRows());
    ny_ = nx_;
  }
  if (nx_ * ny_ != NumMyRows()) {
    AMESOS_CHK_ERR(-1); // sqrt worked badly
  }

  // for square domains, automatically compute how
  // many subdomains will be created along the X-axis and the Y-axis
 
  if (mx_ == -1 || my_ == -1) {
    mx_ = (int)sqrt((double)NumLocalParts());
    my_ = mx_;
	
    if (mx_ * my_ != NumLocalParts()) {
      AMESOS_CHK_ERR(-1); // sqrt worked badly
    }

  } 

  // compute linear decomposition of global
  // how to divide the axis
      
  int nodes_x = nx_ / mx_;
  int nodes_y = ny_ / my_;
  int mody = (ny_ + (ny_ % my_)) / my_;
  int graph;
  int GID;

  for (int i = 0 ; i < NumMyRows() ; ++i) {

    GID = i;

    // 1.- compute the (X,Y) position on the grid
    int ix = GID % nx_;
    int iy = GID / ny_;
    // 2.- compute the domain to which the node belongs to
    int iix = ix / nodes_x;
    int iiy = iy / nodes_y;
    // 3.- compute the global numbering of subgraph (iix, iiy)
    graph = iix + iiy * mx_;
    // 4.- assign this number to node i
    Partition_[i] = graph;

  }

  return(0);

}

//==============================================================================
int Amesos_Partitioner::ComputeGreedyPartition()
{

  vector<int> ElementsPerPart;
  ElementsPerPart.resize(NumLocalParts());

  vector<int> count;
  count.resize(NumLocalParts());

  // define how many nodes have to be put on each part

  int div = NumMyRows() / NumLocalParts();
  int mod = NumMyRows() % NumLocalParts();

  for (int i = 0 ; i < NumLocalParts() ; ++i) {
    count[i] = 0;
    ElementsPerPart[i] = div;
    if (i < mod) ElementsPerPart[i]++;
  }

  for( int i=0 ; i<NumMyRows() ; ++i ) {
    Partition_[i] = -1;
  }

  int MaxNnzPerRow = MaxNumEntries();

  int CrsNumEntries;

  // start from row 0, assigned to domain 0
  Partition_[0] = 0;      
  int CurrentPart = 0;
  vector<int> Indices;
  Indices.resize(MaxNumEntries());
  
  bool ok = true;

  int RootNode = RootNode_;

  while( ok == true ) {

    ExtractMyRowCopy(RootNode, MaxNumEntries(),
		     CrsNumEntries, &Indices[0]);

    ok = false;

    for (int j = 0 ; j < CrsNumEntries ; ++j) {

      if (Indices[j] > NumMyRows()) continue;

      if (count[CurrentPart] == ElementsPerPart[CurrentPart]) {
	CurrentPart++;
      }

      if (Partition_[Indices[j]] == -1) {
	Partition_[Indices[j]] = CurrentPart;
	if( ok == false ) {
	  ok = true;
	  RootNode = Indices[j];
	}
	count[CurrentPart]++;
      }
    }

    // check if some -1 nodes are still available
    if (ok == false) {
      for (int j = 0 ; j < NumMyRows() ; ++j) {
	if (Partition_[j] == -1 ) {
	  RootNode = j;
	  ok = true;
	  break;
	}
      }
    }

  }

  return(0);
}

//==============================================================================
int Amesos_Partitioner::ComputeMETISPartition()
{

  int ierr;
  int nbytes;
  int nbytes_min;
  int nbytes_max;
  int edgecut;

  vector<char> DirichletRows;
  DirichletRows.resize(NumMyRows());

  vector<int> perm;
  perm.resize(NumMyRows());

  int NumMETISNonzeros = 0; // cannot use NumMyNonzeros() because
                        // I want to drop the Dirichlet rows
			// (rows with diagonal entry only)
  int NumMETISRows = 0;
  int NumIndices;
  vector<int> Indices;
  Indices.resize(MaxNumEntries());

  for (int i = 0; i < NumMyRows() ; ++i) {  

    ierr = ExtractMyRowCopy(i, MaxNumEntries(), NumIndices, &Indices[0]);
    AMESOS_CHK_ERR(ierr);
    
    if (NumIndices <= 1) {
      DirichletRows[i] = 'T';
      perm[i] = -1;
    } else {
      perm[i] = NumMETISRows++;
      DirichletRows[i] = 'F';
      NumMETISNonzeros += NumIndices;
    }
  }

  /* construct the CSR graph information of the LOCAL matrix
     using the get_row function */

  vector<idxtype> wgtflag;
  wgtflag.resize(4);

  vector<int> options;
  options.resize(4);
  
  int numflag;

  /* set parameters */
   
  wgtflag[0] = 0;    /* no weights */
  numflag    = 0;    /* C style */
  options[0] = 0;    /* default options */
   
  vector<idxtype> xadj;
  xadj.resize(NumMETISRows+1);

  vector<idxtype> adjncy;
  adjncy.resize(NumMETISNonzeros);
   
  int count = 0; 
  int count2 = 0; 
  xadj[0] = 0;
  
  for (int i = 0; i < NumMyRows() ; ++i) {

    if( DirichletRows[i] == 'F' ) {

      xadj[count2+1] = xadj[count2]; /* nonzeros in row i-1 */
    
      int NumIndices;
      ierr = ExtractMyRowCopy(i, MaxNumEntries(), NumIndices, &Indices[0]);
      AMESOS_CHK_ERR(ierr);

      /* need to avoid boundary nodes in METIS vectors. Skip them */
      /* (I am not pretty sure that rows with zero elements are   */
      /* well eated by METIS.) perm has been allocates of size    */
      /* Nrows, so columns corresponding to external nodes can not*/
      /* be given as input to perm                                */

      for (int j = 0 ; j < NumIndices ; ++j) {
	int jj = Indices[j];
	if (jj < NumMyRows()) {
	  if ((jj != i) && (perm[jj] != -1)) {
	    adjncy[count++] = perm[jj];
	    xadj[count2+1]++;
	  }
	}
      }
      count2++;
    }      
  }

  if ((count > NumMETISNonzeros) || (count2 != NumMETISRows)) {
    AMESOS_CHK_ERR(-11); // possible buffer overflow?
  }

  vector<idxtype> part;
  part.resize(NumMETISRows);

  vector<idxtype> NodesInSubgraph;
  NodesInSubgraph.resize(NumLocalParts_);

  /* ********************************************************************** */
  /* Before calling METIS, I verify that the two extreme situations are     */
  /* handled separately.                                                    */
  /* ********************************************************************** */
  
  int ok;

  if (NumLocalParts_ == 1) {

    if ((verbose_ > 8) && (Comm().MyPID() == 0)) {
      cout << PrintMsg_ << "The number of local parts is 1." << endl;
      cout << PrintMsg_ << "I put all the local rows in the same part," << endl;
      cout << PrintMsg_ << "METIS is not called." << endl;
    }

    for (int i = 0 ; i < NumMETISRows ; ++i) 
      part[i] = 0;
    
  } else if (NumLocalParts_ == NumMETISRows) {

    if ((verbose_ > 8) && (Comm().MyPID() == 0)) {
      cout << PrintMsg_ << "The requested number of local parts" << endl;
      cout << PrintMsg_ << "equal the number of non-Dirichlet rows." << endl;
      cout << PrintMsg_ << "I decompose the rows in a simple way," << endl;
      cout << PrintMsg_ << "METIS is not called." << endl;
    }
    for (int i = 0 ; i < NumMETISRows ; ++i) 
      part[i] = i;
  
  } else {

    ok = 0;

    while (ok == 0) {
      
      /* ****************************************************************** */
      /* Put -1 in part, so I can verify that METIS has filled each pos    */
      /* ****************************************************************** */

      for (int i = 0 ; i < NumMETISRows ; ++i) 
	part[i] = -1;
    
      /* ****************************************************************** */
      /* Estimate memory required by METIS. This memory will be dynamically */
      /* allocated inside; however this is a good estimation of how METIS   */
      /* will cost in terms of memory.                                      */
      /* Then, call METIS.                                                  */
      /* ****************************************************************** */

#ifdef HAVE_AMESOS_METIS
      if (NumLocalParts_ < 8) {

	int i = 1; /* optype in the METIS manual */
	numflag = 0;
	METIS_EstimateMemory(&NumMETISRows, &xadj[0], &adjncy[0], &numflag,
			     &i, &nbytes );
	
	METIS_PartGraphRecursive(&NumMETISRows, &xadj[0], &adjncy[0],
				 NULL, NULL,
				 &wgtflag[0], &numflag, &NumLocalParts_, 
				 &options[0], &edgecut, &part[0]);
      } else {

	int i = 2;
	numflag = 0;

	METIS_EstimateMemory(&NumMETISRows, &xadj[0], &adjncy[0], &numflag,
			     &i, &nbytes );
	
	METIS_PartGraphKway (&NumMETISRows, &xadj[0], &adjncy[0], NULL, 
			     NULL, &wgtflag[0], &numflag, 
			     &NumLocalParts_, &options[0],
			     &edgecut, &part[0]);
      }
#else
      for (int i = 0 ; i < NumMETISRows ; ++i) 
	part[i] = 0;
      NumLocalParts_ = 1;
#endif
      
      /* **************************************************************** */
      /* perform some checks. If aggregates with zero assigned nodes      */
      /* exist, then recall METIS, asking for a smaller number of sub     */
      /* graphs. This is the role of the `ok' variable.                   */
      /* Also, if the part vector contains some junk, recall METIS        */
      /* **************************************************************** */

      ok = 1;
      
      for (int i = 0 ; i < NumLocalParts_ ; ++i) 
	NodesInSubgraph[i] = 0;

      for (int i = 0 ; i < NumMETISRows ; ++i) {
	int j = part[i];
	if ((j < 0) || (j>= NumLocalParts_)) {
	  ok = 0;
	  break;
	} 
	else NodesInSubgraph[j]++;
      }
      
      for (int i = 0 ; i < NumLocalParts_ ; ++i) {
	if( NodesInSubgraph[i] == 0 ) {
	  ok = 0;
	  break;
	}
      }
      
      if (ok == 0) {
	cerr << "Specified number of subgraphs ("
	     << NumLocalParts_ << ") generates empty subgraphs." << endl;
	cerr << "Now I recall METIS with NumLocalParts_ = "
	     << NumLocalParts_ / 2 << "..." << endl;
	NumLocalParts_ = NumLocalParts_/2;
      }
      
      if (NumLocalParts_ == 0) {
	AMESOS_CHK_ERR(-20); // something went wrong
      }
      
      /* ************************************************************** */
      /* handle the case NumLocalParts_ = 1 separately. Do not recall METIS    */
      /* in this case, simply put everything to zero and continue       */
      /* ************************************************************** */
      
      if (NumLocalParts_ == 1) {
	for (int i = 0 ; i < NumMETISRows ; ++i) 
	  part[i] = 0;
	ok = 1;
      }
      
    } /* while( ok == 0 ) */
  
  } /* if( NumLocalParts_ == 1 ) */

  /* ********************************************************************** */
  /* Some fancy output for memory usage.                                    */
  /* ********************************************************************** */

  nbytes /= 1024;
  
  Comm().MaxAll(&nbytes,&nbytes_max,1);
  Comm().MinAll(&nbytes,&nbytes_min,1);

  if ((verbose_ > 8) && (Comm().MyPID() == 0)) {
   
    cout << PrintMsg_ << "Mininum estimated memory for METIS = "
         << nbytes_min << "Kb" << endl;
    cout << PrintMsg_ << "Maximum estimated memory for METIS = " 
         << nbytes_max << "Kb" << endl;
  }
  
  /* copy back part into aggr_index, and set to -1
     the aggr_index corresponding to ghost nodes */

  NumSingletons_ = 0;
  
  for (int i = 0 ; i < NumMyRows() ; ++i) {
    int j = perm[i];
    if( j != -1 ) 
      Partition_[i] = (int)part[j];
    else {
      NumSingletons_++;
      Partition_[i] = -7;
    }
  }

  assert (NumSingletons_ == NumMyRows() - NumMETISRows);
  return(0);

  // FIXME: some timing?
} /* ML_DecomposeGraph_with_METIS */

//============================================================================

int Amesos_Partitioner::ComputeOverlappingPartitions()
{

  // start defining the subgraphs for no overlap

  vector<int> sizes;
  sizes.resize(NumLocalParts_);

  // 1.- compute how many rows are in each subgraph
  for (int i = 0 ; i < NumLocalParts_ ; ++i)
    sizes[i] = 0;

  for (int i = 0 ; i < NumMyRows() ; ++i) {
    assert (Partition_[i] < NumLocalParts_);
    // discard singletons
    if (Partition_[i] != -7)
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
    if (part == -7)
      continue;
    int count = sizes[part];
    Parts_[part][count] = Graph_->GRID(i);
    sizes[part]++;
  }

  if (OverlappingLevel_ == 0)
    return(0);

  // wider overlap requires further computations
 
  for (int level = 1 ; level <= OverlappingLevel_ ; ++level) {

    vector<int> tmp[NumLocalParts_];

    // cycle over all rows in the local graph (that is the overlapping
    // graph). For each row, all columns will belong to the subgraph of
    // row `i'.

    for (int part = 0 ; part < NumLocalParts_ ; ++part) {

      for (int i = 0; i < Parts_[part].size() ; ++i) {  

	int row = Parts_[part][i];
	int LRID = OverlappingGraph()->LRID(row);
	int NumIndices;
	int* Indices;
	int ierr = OverlappingGraph()->ExtractMyRowView(LRID, NumIndices, Indices);
	AMESOS_CHK_ERR(ierr);

	for (int j = 0 ; j < NumIndices ; ++j) {

	  // need global ordering
	  int col = OverlappingGraph()->GCID(Indices[j]);
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
