#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Ifpack_Partitioner.h"
#include "Ifpack_OverlappingPartitioner.h"
#include "Ifpack_METISPartitioner.h"
#include "Ifpack_Graph.h"

#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Teuchos_ParameterList.hpp"

typedef int idxtype;
#ifdef HAVE_IFPACK_METIS
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
int Ifpack_METISPartitioner::ComputePartitions()
{

  int ierr;
  int nbytes = 0;
  int nbytes_min;
  int nbytes_max;
  int edgecut;

  int NumIndices;
  vector<int> Indices;
  Indices.resize(MaxNumEntries());

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
  xadj.resize(NumMyNonDirichletRows() + 1);

  vector<idxtype> adjncy;
  adjncy.resize(NumMyNonzeros());
   
  int count = 0; 
  int count2 = 0; 
  xadj[0] = 0;
  
  for (int i = 0; i < NumMyRows() ; ++i) {

    if (Mask_[i] != -1) {

      xadj[count2+1] = xadj[count2]; /* nonzeros in row i-1 */
    
      ierr = ExtractMyRowCopy(i, MaxNumEntries(), NumIndices, &Indices[0]);
      IFPACK_CHK_ERR(ierr);

      /* need to avoid boundary nodes in METIS vectors. Skip them */
      /* (I am not pretty sure that rows with zero elements are   */
      /* well eated by METIS.) perm has been allocates of size    */
      /* Nrows, so columns corresponding to external nodes can not*/
      /* be given as input to perm                                */

      for (int j = 0 ; j < NumIndices ; ++j) {
	int jj = Indices[j];
	if (jj < NumMyRows()) {
	  if ((jj != i) && (Mask_[jj] != -1)) {
	    adjncy[count++] = Mask_[jj];
	    xadj[count2+1]++;
	  }
	}
      }
      count2++;
    }      
  }

  vector<idxtype> part;
  part.resize(NumMyNonDirichletRows());

  vector<idxtype> NodesInSubgraph;
  NodesInSubgraph.resize(NumLocalParts_);

  // some cases can be handled separately
  
  int ok;

  if (NumLocalParts() == 1) {

    if ((verbose_ > 8) && (Comm().MyPID() == 0)) {
      cout << PrintMsg_ << "The number of local parts is 1." << endl;
      cout << PrintMsg_ << "I put all the local rows in the same part," << endl;
      cout << PrintMsg_ << "METIS is not called." << endl;
    }

    for (int i = 0 ; i < NumMyNonDirichletRows() ; ++i) 
      part[i] = 0;
    
  } else if (NumLocalParts() == NumMyNonDirichletRows()) {

    if ((verbose_ > 8) && (Comm().MyPID() == 0)) {
      cout << PrintMsg_ << "The requested number of local parts" << endl;
      cout << PrintMsg_ << "equal the number of non-Dirichlet rows." << endl;
      cout << PrintMsg_ << "I decompose the rows in a simple way," << endl;
      cout << PrintMsg_ << "METIS is not called." << endl;
    }
    for (int i = 0 ; i < NumMyNonDirichletRows() ; ++i) 
      part[i] = i;
  
  } else {

    ok = 0;

    while (ok == 0) {
      
      /* ****************************************************************** */
      /* Put -1 in part, so I can verify that METIS has filled each pos    */
      /* ****************************************************************** */

      for (int i = 0 ; i < NumMyNonDirichletRows() ; ++i) 
	part[i] = -1;
    
      /* ****************************************************************** */
      /* Estimate memory required by METIS. This memory will be dynamically */
      /* allocated inside; however this is a good estimation of how METIS   */
      /* will cost in terms of memory.                                      */
      /* Then, call METIS.                                                  */
      /* ****************************************************************** */

#ifdef HAVE_IFPACK_METIS
      if (NumLocalParts_ < 8) {

	int i = 1; /* optype in the METIS manual */
	numflag = 0;
	METIS_EstimateMemory(&NumMyNonDirichletRows_, &xadj[0], &adjncy[0], 
			     &numflag, &i, &nbytes );
	
	METIS_PartGraphRecursive(&NumMyNonDirichletRows_, &xadj[0], &adjncy[0],
				 NULL, NULL,
				 &wgtflag[0], &numflag, &NumLocalParts_, 
				 &options[0], &edgecut, &part[0]);
      } else {

	int i = 2;
	numflag = 0;

	METIS_EstimateMemory(&NumMyNonDirichletRows_, &xadj[0], &adjncy[0], 
			     &numflag, &i, &nbytes );
	
	METIS_PartGraphKway (&NumMyNonDirichletRows_, &xadj[0], &adjncy[0], 
			     NULL, 
			     NULL, &wgtflag[0], &numflag, 
			     &NumLocalParts_, &options[0],
			     &edgecut, &part[0]);
      }
#else
      if (Comm().MyPID() == 0) {
	cerr << ErrorMsg_ << "METIS was not linked; now I put all" << endl;
	cerr << ErrorMsg_ << "the local nodes in the same partition." << endl;
      }
      for (int i = 0 ; i < NumMyNonDirichletRows() ; ++i) 
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
      
      for (int i = 0 ; i < NumLocalParts() ; ++i) 
	NodesInSubgraph[i] = 0;

      for (int i = 0 ; i < NumMyNonDirichletRows() ; ++i) {
	int j = part[i];
	if ((j < 0) || (j>= NumLocalParts())) {
	  ok = 0;
	  break;
	} 
	else NodesInSubgraph[j]++;
      }
      
      for (int i = 0 ; i < NumLocalParts() ; ++i) {
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
      
      if (NumLocalParts() == 0) {
	IFPACK_CHK_ERR(-10); // something went wrong
      }
      
      if (NumLocalParts() == 1) {
	for (int i = 0 ; i < NumMyNonDirichletRows() ; ++i) 
	  part[i] = 0;
	ok = 1;
      }
      
    } /* while( ok == 0 ) */
  
  } /* if( NumLocalParts_ == 1 ) */

  // some fancy output for memory usage

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

  int NumSingletons = 0;
  
  for (int i = 0 ; i < NumMyRows() ; ++i) {
    int j = Mask_[i];
    if( j != -1 ) 
      Partition_[i] = part[j];
    else {
      NumSingletons++;
      Partition_[i] = -1;
    }
  }

  if (NumSingletons != NumMyRows() - NumMyNonDirichletRows())
    IFPACK_CHK_ERR(-10);

  return(0);

  // FIXME: some timing?
} 
#endif
