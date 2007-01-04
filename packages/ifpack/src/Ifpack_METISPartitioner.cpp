#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Partitioner.h"
#include "Ifpack_OverlappingPartitioner.h"
#include "Ifpack_METISPartitioner.h"
#include "Ifpack_Graph.h"
#include "Ifpack_Graph_Epetra_CrsGraph.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Map.h"
#include "Teuchos_ParameterList.hpp"

// may need to change this for wierd installations
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
// NOTE:
// - matrix is supposed to be localized, and passes through the
// singleton filter. This means that I do not have to look
// for Dirichlet nodes (singletons). Also, all rows and columns are 
// local.
int Ifpack_METISPartitioner::ComputePartitions()
{

  int ierr;
#ifdef HAVE_IFPACK_METIS
  int nbytes = 0;
  int edgecut;
#endif

  Epetra_CrsGraph* SymGraph = 0;
  Epetra_Map* SymMap = 0;
  Ifpack_Graph_Epetra_CrsGraph* SymIFPACKGraph = 0;
  Ifpack_Graph* IFPACKGraph = (Ifpack_Graph*)Graph_;

  int Length = 2 * MaxNumEntries();
  int NumIndices;
  vector<int> Indices;
  Indices.resize(Length);

  /* construct the CSR graph information of the LOCAL matrix
     using the get_row function */

  vector<idxtype> wgtflag;
  wgtflag.resize(4);

  vector<int> options;
  options.resize(4);
  
  int numflag;

  if (UseSymmetricGraph_) {

    // need to build a symmetric graph. 
    // I do this in two stages:
    // 1.- construct an Epetra_CrsMatrix, symmetric
    // 2.- convert the Epetra_CrsMatrix into METIS format
    SymMap = new Epetra_Map(NumMyRows(),0,Graph_->Comm());
    SymGraph = new Epetra_CrsGraph(Copy,*SymMap,0);

    for (int i = 0; i < NumMyRows() ; ++i) {

      ierr = Graph_->ExtractMyRowCopy(i, Length, NumIndices, 
				      &Indices[0]);
      IFPACK_CHK_ERR(ierr);

      for (int j = 0 ; j < NumIndices ; ++j) {
	int jj = Indices[j];
	if (jj != i) {
	  SymGraph->InsertGlobalIndices(i,1,&jj);
	  SymGraph->InsertGlobalIndices(jj,1,&i);
	}
      }      
    }
    IFPACK_CHK_ERR(SymGraph->FillComplete());
    SymIFPACKGraph = new Ifpack_Graph_Epetra_CrsGraph(SymGraph);
    IFPACKGraph = SymIFPACKGraph;
  }

  // now work on IFPACKGraph, that can be the symmetric or
  // the non-symmetric one

  /* set parameters */
   
  wgtflag[0] = 0;    /* no weights */
  numflag    = 0;    /* C style */
  options[0] = 0;    /* default options */
   
  vector<idxtype> xadj;
  xadj.resize(NumMyRows() + 1);

  vector<idxtype> adjncy;
  adjncy.resize(NumMyNonzeros());
   
  int count = 0; 
  int count2 = 0; 
  xadj[0] = 0;
  
  for (int i = 0; i < NumMyRows() ; ++i) {

    xadj[count2+1] = xadj[count2]; /* nonzeros in row i-1 */

    ierr = IFPACKGraph->ExtractMyRowCopy(i, Length, NumIndices, &Indices[0]);
    IFPACK_CHK_ERR(ierr);

    for (int j = 0 ; j < NumIndices ; ++j) {
      int jj = Indices[j];
      if (jj != i) {
	adjncy[count++] = jj;
	xadj[count2+1]++;
      }
    }
    count2++;
  }

  vector<idxtype> NodesInSubgraph;
  NodesInSubgraph.resize(NumLocalParts_);

  // some cases can be handled separately
  
  int ok;

  if (NumLocalParts() == 1) {

    for (int i = 0 ; i < NumMyRows() ; ++i) 
      Partition_[i] = 0;
    
  } else if (NumLocalParts() == NumMyRows()) {

    for (int i = 0 ; i < NumMyRows() ; ++i) 
      Partition_[i] = i;
  
  } else {

    ok = 0;

    // sometimes METIS creates less partitions than specified.
    // ok will check this problem, and recall metis, asking
    // for NumLocalParts_/2 partitions
    while (ok == 0) {
      
      for (int i = 0 ; i < NumMyRows() ; ++i) 
	Partition_[i] = -1;
    
#ifdef HAVE_IFPACK_METIS
      int j = NumMyRows();
      if (NumLocalParts_ < 8) {

	int i = 1; /* optype in the METIS manual */
	numflag = 0;
	METIS_EstimateMemory(&j, &xadj[0], &adjncy[0], 
			     &numflag, &i, &nbytes );
	
	METIS_PartGraphRecursive(&j, &xadj[0], &adjncy[0],
				 NULL, NULL,
				 &wgtflag[0], &numflag, &NumLocalParts_, 
				 &options[0], &edgecut, &Partition_[0]);
      } else {

	numflag = 0;
	
	METIS_PartGraphKway (&j, &xadj[0], &adjncy[0], 
			     NULL, 
			     NULL, &wgtflag[0], &numflag, 
			     &NumLocalParts_, &options[0],
			     &edgecut, &Partition_[0]);
      }
#else
      numflag = numflag * 2; // avoid warning for unused variable
      if (Graph_->Comm().MyPID() == 0) {
	cerr << "METIS was not linked; now I put all" << endl;
	cerr << "the local nodes in the same partition." << endl;
      }
      for (int i = 0 ; i < NumMyRows() ; ++i) 
	Partition_[i] = 0;
      NumLocalParts_ = 1;
#endif
      
      ok = 1;
      
      for (int i = 0 ; i < NumLocalParts() ; ++i) 
	NodesInSubgraph[i] = 0;

      for (int i = 0 ; i < NumMyRows() ; ++i) {
	int j = Partition_[i];
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
	for (int i = 0 ; i < NumMyRows() ; ++i) 
	  Partition_[i] = 0;
	ok = 1;
      }
      
    } /* while( ok == 0 ) */
  
  } /* if( NumLocalParts_ == 1 ) */

  if (SymGraph)
    delete SymGraph;
  if (SymMap)
    delete SymMap;
  if (SymIFPACKGraph)
    delete SymIFPACKGraph;

  return(0);
} 
