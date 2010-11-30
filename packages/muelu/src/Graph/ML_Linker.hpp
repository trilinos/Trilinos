#include "ml_aggregate.h"
#include "ml_epetra_utils.h"
extern MueLu_Graph  *MueLu_BuildGraph(ML_Operator *Amatrix, char *name);
extern int MueLu_DestroyGraph(MueLu_Graph *Graph);

/**********************************************************************************/
/* Function to execute new MueLu aggregation via old ML                          */
/* This function should go away soon as we should start executing new MueLu      */
/* aggregation inside MueLu.
/**********************************************************************************/
int ML_Aggregate_CoarsenUncoupled(ML_Aggregate *ml_ag, 
           ML_Operator *Amatrix, ML_Operator **Pmatrix, ML_Comm *comm)
{
   MueLu_Graph     *Graph;
   Graph = MueLu_BuildGraph(Amatrix,"ML_Uncoupled");

   if (Graph->EGraph->Comm().MyPID() == 0 && ml_ag->print_flag  < MueLu_PrintLevel())
       printf("ML_Aggregate_CoarsenUncoupled : \n");

   MueLu_AggOptions AggregateOptions;

   AggregateOptions.print_flag                 = ml_ag->print_flag;
   AggregateOptions.min_nodes_per_aggregate    = ml_ag->min_nodes_per_aggregate;
   AggregateOptions.max_neigh_already_selected = ml_ag->max_neigh_already_selected;
   AggregateOptions.ordering                   = ml_ag->ordering;
   AggregateOptions.phase3_agg_creation        = ml_ag->phase3_agg_creation;


   MueLu_Aggregate *Aggregates = NULL;

   Aggregates = MueLu_Aggregate_CoarsenUncoupled(&AggregateOptions,Graph);


   MueLu_AggregateLeftOvers(&AggregateOptions, Aggregates, "UC_CleanUp", Graph);

//#ifdef out
Epetra_IntVector Final( Aggregates->Vertex2AggId->Map() );
for (int i = 0; i < Aggregates->Vertex2AggId->Map().NumMyElements(); i++) 
   Final[i] = (*(Aggregates->Vertex2AggId))[i] + (*(Aggregates->ProcWinner))[i]*1000;
printf("finals\n");
cout << Final << endl; sleep(2);
//#endif

   MueLu_AggregateDestroy(Aggregates); 
   MueLu_DestroyGraph(Graph);
   return 0;
}

/**********************************************************************************/
/* Function to take an ML_Operator (which actually wraps an Epetra_CrsMatrix) and */
/* extract out the Epetra_CrsGraph. My guess is that this should be changed soon  */
/* so that the first argument is some MueLu API Operator.                        */
/**********************************************************************************/
MueLu_Graph *MueLu_BuildGraph(ML_Operator *Amatrix, char *name)
{
  MueLu_Graph *Graph;
  double *dtmp = NULL;
  Epetra_CrsMatrix *A;

  Graph = (MueLu_Graph *) malloc(sizeof(MueLu_Graph));
  Graph->EGraph     = NULL;
  Graph->name = NULL;
  Graph->name       = (char *) malloc(sizeof(char)*80); strcpy(Graph->name,name);
  Graph->NVertices  = Amatrix->invec_leng;

  if ( Amatrix->getrow->Nrows == 0) {
     Graph->VertexNeighbors    = NULL;
     Graph->VertexNeighborsPtr = NULL;
     Graph->NEdges             = 0;
  }
  else {
     Epetra_ML_GetCrsDataptrs(Amatrix, &dtmp, &(Graph->VertexNeighbors),&(Graph->VertexNeighborsPtr));
     if ( Graph->VertexNeighborsPtr == NULL) {
        printf("MueLu_BuildGraph: Only functions for an Epetra_CrsMatrix.\n");
        exit(1);
      }
      Graph->NEdges      = (Graph->VertexNeighborsPtr)[Amatrix->getrow->Nrows];
      Epetra_ML_GetCrsMatrix( Amatrix, (void **) &A );
      Graph->EGraph = &(A->Graph());
  }
  if (Graph->EGraph == NULL) Graph->NGhost = 0;
  else {
     Graph->NGhost = A->RowMatrixColMap().NumMyElements() - A->OperatorDomainMap().NumMyElements();
     if (Graph->NGhost < 0) Graph->NGhost = 0;
  }
  return Graph;
}
int MueLu_DestroyGraph(MueLu_Graph *Graph)
{
   if ( Graph != NULL) {
      if (Graph->name != NULL) free(Graph->name);
      free(Graph);
   }
   return 0;
}
