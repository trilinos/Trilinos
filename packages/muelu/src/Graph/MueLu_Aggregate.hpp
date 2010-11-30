#ifndef MUELU_AGGREGATE_HPP
#define MUELU_AGGREGATE_HPP

/***************************************************************************** 
   Structure holding aggregate information. Right now, NAggregates, IsRoot,
   Vertex2AggId, ProcWinner are populated.  This allows us to look at a node
   and determine the aggregate to which it has been assigned and the id of the 
   processor that owns this aggregte.  It is not so easy to determine vertices
   within the kth aggregate or the size of the kth aggregate. Thus, it might be
   useful to have a secondary structure which would be a rectangular CrsGraph 
   where rows (or vertices) correspond to aggregates and colunmns (or edges) 
   correspond to nodes.  While not strictly necessary, it might be convenient.
 *****************************************************************************/
typedef struct MueLoo_Aggregate_Struct
{
   char *name;
   int  NAggregates;              /* Number of aggregates on this processor  */

   Epetra_IntVector *Vertex2AggId;/* Vertex2AggId[k] gives a local id        */
                                  /* corresponding to the aggregate to which */
                                  /* local id k has been assigned.  While k  */
   Epetra_IntVector *ProcWinner;  /* is the local id on my processor (MyPID),*/
                                  /* Vertex2AggId[k] is the local id on the  */
                                  /* processor which actually owns the       */
                                  /* aggregate. This owning processor has id */
                                  /* given by ProcWinner[k].                 */

   bool *IsRoot;                  /* IsRoot[i] indicates whether vertex i  */
                                  /* is a root node.                       */

   Epetra_CrsGraph *Agg2Vertex;   /* Currently not used                    */
} MueLoo_Aggregate;


extern MueLoo_Aggregate *MueLoo_AggregateCreate(MueLoo_Graph *Graph, char *str);
extern int MueLoo_AggregateDestroy(MueLoo_Aggregate *agg);

// Constructors to create aggregates. This should really be replaced by an 
// aggregate class.
MueLoo_Aggregate *MueLoo_AggregateCreate(MueLoo_Graph *Graph, char *str)
{
   MueLoo_Aggregate *Aggregates;

   Aggregates = (MueLoo_Aggregate *) malloc(sizeof(MueLoo_Aggregate));
   Aggregates->name = (char *) malloc(sizeof(char)*80);
   strcpy(Aggregates->name,str);
   Aggregates->NAggregates  = 0;
   Aggregates->Agg2Vertex   = NULL;
   Aggregates->Vertex2AggId = new Epetra_IntVector(Graph->EGraph->ImportMap());
   Aggregates->Vertex2AggId->PutValue(MUELOO_UNAGGREGATED);
   Aggregates->ProcWinner = new Epetra_IntVector(Graph->EGraph->ImportMap());
   Aggregates->ProcWinner->PutValue(MUELOO_UNASSIGNED);
   Aggregates->IsRoot = new bool[Graph->EGraph->ImportMap().NumMyElements()];
   for (int i=0; i < Graph->EGraph->ImportMap().NumMyElements(); i++)
       Aggregates->IsRoot[i] = false;

   return Aggregates;
}
// Destructor for aggregates. This should really be replaced by an 
// aggregate class.
int MueLoo_AggregateDestroy(MueLoo_Aggregate *agg)
{
   if (agg != NULL) {
      if (agg->IsRoot       != NULL) delete [] (agg->IsRoot);
      if (agg->ProcWinner   != NULL) delete agg->ProcWinner;
      if (agg->Vertex2AggId != NULL) delete agg->Vertex2AggId;
      if (agg->name         != NULL) free(agg->name);
      if (agg->Agg2Vertex   != NULL) delete agg->Agg2Vertex;
      free(agg);
   }
   return 0;
}

#endif
