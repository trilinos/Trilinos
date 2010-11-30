#ifndef MUELU_GRAPH_HPP
#define MUELU_GRAPH_HPP

/******************************************************************************
   MueLoo representation of a graph. Some of this is redundant with an 
   Epetra_CrsGraph so we might want to clean up. In particular, things
   like VertexNeighbors, NVertices, NEdges, etc. are available somewhere 
   in EGraph (though perhaps protected).              
******************************************************************************/
typedef struct MueLoo_Graph_Struct
{
   char *name;
   int  NVertices, NEdges, NGhost;
   int  *VertexNeighborsPtr;  /* VertexNeighbors[VertexNeighborsPtr[i]:     */
   int  *VertexNeighbors;     /*                 VertexNeighborsPtr[i+1]-1] */
                              /* corresponds to vertices Adjacent to vertex */
                              /* i in graph                                 */ 
   const Epetra_CrsGraph *EGraph;
} MueLoo_Graph;

#endif
