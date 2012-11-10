// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
// Unmaintained code snippet to use Muelu's aggregation algorithms in ML

#include "ml_aggregate.h"
#include "ml_epetra_utils.h"

extern MueLu_Graph  *MueLu_BuildGraph(ML_Matrix *Amatrix, char *name);
extern int MueLu_DestroyGraph(MueLu_Graph *graph);

/**********************************************************************************/
/* Function to execute new MueLu aggregation via old ML                          */
/* This function should go away soon as we should start executing new MueLu      */
/* aggregation inside MueLu.
/**********************************************************************************/
int ML_Aggregate_CoarsenUncoupled(ML_Aggregate *mlAggregates,
           ML_Matrix *Amatrix, ML_Matrix **Pmatrix, ML_Comm *comm)
{
   MueLu_Graph     *graph;
   graph = MueLu_BuildGraph(Amatrix,"ML_Uncoupled");

   if (graph->eGraph->Comm().MyPID() == 0 && mlAggregates->printFlag  < MueLu_PrintLevel())
       printf("ML_Aggregate_CoarsenUncoupled : \n");

   MueLu_AggOptions aggregateOptions;

   aggregateOptions.printFlag               = mlAggregates->print_flag;
   aggregateOptions.minNodesPerAggregate    = mlAggregates->min_nodes_per_aggregate;
   aggregateOptions.maxNeighAlreadySelected = mlAggregates->max_neigh_already_selected;
   aggregateOptions.ordering                = mlAggregates->ordering;
   aggregateOptions.phase3AggCreation       = mlAggregates->phase3_agg_creation;


   Aggregates *aggregates = NULL;

   aggregates = MueLu_Aggregate_CoarsenUncoupled(&aggregateOptions,graph);


   MueLu_AggregateLeftOvers(&aggregateOptions, aggregates, "UC_CleanUp", graph);

//#ifdef out
Epetra_IntVector Final( aggregates->vertex2AggId->Map() );
for (int i = 0; i < aggregates->vertex2AggId->Map().NumMyElements(); i++)
   Final[i] = (*(aggregates->vertex2AggId))[i] + (*(aggregates->procWinner))[i]*1000;
printf("finals\n");
cout << Final << endl; sleep(2);
//#endif

   MueLu_AggregateDestroy(aggregates);
   MueLu_DestroyGraph(graph);
   return 0;
}

/**********************************************************************************/
/* Function to take an ML_Matrix (which actually wraps an Epetra_CrsMatrix) and */
/* extract out the Epetra_CrsGraph. My guess is that this should be changed soon  */
/* so that the first argument is some MueLu API Matrix.                        */
/**********************************************************************************/
MueLu_Graph *MueLu_BuildGraph(ML_Matrix *Amatrix, char *name)
{
  MueLu_Graph *graph;
  double *dtmp = NULL;
  Epetra_CrsMatrix *A;

  graph = (MueLu_Graph *) malloc(sizeof(MueLu_Graph));
  graph->eGraph     = NULL;
  graph->name = NULL;
  graph->name       = (char *) malloc(sizeof(char)*80); strcpy(Graph->name,name);
  graph->nVertices  = Amatrix->invec_leng;

  if ( Amatrix->getrow->nrows == 0) {
     graph->vertexNeighbors    = NULL;
     graph->vertexNeighborsPtr = NULL;
     graph->nEdges             = 0;
  }
  else {
     Epetra_ML_GetCrsDataptrs(Amatrix, &dtmp, &(graph->vertexNeighbors),&(graph->vertexNeighborsPtr));
     if ( graph->vertexNeighborsPtr == NULL) {
        printf("MueLu_BuildGraph: Only functions for an Epetra_CrsMatrix.\n");
        exit(1);
      }
      graph->nEdges      = (graph->vertexNeighborsPtr)[Amatrix->getrow->Nrows];
      Epetra_ML_GetCrsMatrix( Amatrix, (void **) &A );
      graph->eGraph = &(A->graph());
  }
  if (graph->eGraph == NULL) graph->nGhost = 0;
  else {
     graph->nGhost = A->RowMatrixColMap().NumMyElements() - A->MatrixDomainMap().NumMyElements();
     if (graph->nGhost < 0) graph->nGhost = 0;
  }
  return graph;
}

int MueLu_DestroyGraph(MueLu_Graph *graph)
{
   if ( graph != NULL) {
      if (graph->name != NULL) free(graph->name);
      free(graph);
   }
   return 0;
}
