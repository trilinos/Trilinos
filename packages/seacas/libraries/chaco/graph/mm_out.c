/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include <string.h>
#include "defs.h"
#include "structs.h"

/* Print out subgraph in matrix-market format. */
void      mm_out(graph, nvtxs, using_ewgts, tag, file_name)
struct vtx_data **graph;	/* graph data structure */
int       nvtxs;		/* number of vtxs in graph */
int       using_ewgts;		/* Are edges weighted? */
char     *tag;			/* message to include */
char     *file_name;		/* output file name if not null */
{
    FILE     *file;		/* output file */
    int       nedges;		/* number of edges in graph */
    int       i, j;		/* loop counter */
    int       DIAG;

DIAG = TRUE;

    if (file_name != NULL)
	file = fopen(file_name, "w");
    else
	file = stdout;

    /* Determine all the appropriate parameters. */
    nedges = 0;
    for (i = 1; i <= nvtxs; i++) {
	nedges += graph[i]->nedges - 1;
    }
if (DIAG) nedges += nvtxs;


    if (tag != NULL)
	fprintf(file, "%% graph_out: %s\n", tag);
    fprintf(file, " %d %d %d\n", nvtxs, nvtxs, nedges);
    for (i = 1; i <= nvtxs; i++) {
if (DIAG) {
 if (!using_ewgts)
  fprintf(file, "%d %d\n", i, i);
 else
  fprintf(file, "%d %d %.9f\n", i, i, 1.0);
}
	for (j = 1; j < graph[i]->nedges; j++) {
	        if (!using_ewgts)
	            fprintf(file, "%d %d\n", i, graph[i]->edges[j]);
	        else
	            fprintf(file, "%d %d %.9f\n", i, graph[i]->edges[j], graph[i]->ewgts[j]);
	}
    }

    if (file_name != NULL)
	fclose(file);
}
