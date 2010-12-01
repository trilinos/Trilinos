#ifndef MUELU_GRAPH_HPP
#define MUELU_GRAPH_HPP

/******************************************************************************
   MueLu representation of a graph. Some of this is redundant with an 
   Epetra_CrsGraph so we might want to clean up. In particular, things
   like VertexNeighbors, NVertices, NEdges, etc. are available somewhere 
   in EGraph (though perhaps protected).              
******************************************************************************/
typedef struct MueLu_Graph_Struct
{
   char *name;
   int  nVertices, nEdges, nGhost;
   int  *vertexNeighborsPtr;  /* VertexNeighbors[VertexNeighborsPtr[i]:     */
   int  *vertexNeighbors;     /*                 VertexNeighborsPtr[i+1]-1] */
                              /* corresponds to vertices Adjacent to vertex */
                              /* i in graph                                 */ 
   const Epetra_CrsGraph *eGraph;
} MueLu_Graph;

#endif
