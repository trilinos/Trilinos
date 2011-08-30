/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/



/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "ch_input_const.h"
#include "dr_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

int chaco_input_graph(
ZOLTAN_FILE *fin,			/* input file */
char     *inname,		/* name of input file */
int     **start,		/* start of edge list for each vertex */
int     **adjacency,		/* edge list data */
int      *nvtxs,		/* number of vertices in graph */
int      *vwgt_dim,		/* # of vertex weights per node */
float   **vweights,		/* vertex weight list data */
int      *ewgt_dim,		/* # of edge weights per edge */
float   **eweights 		/* edge weight list data */
)
{
    const char     *yo = "chaco_input_graph";
    int      *adjptr;		/* loops through adjacency data */
    float    *ewptr;		/* loops through edge weight data */
    int       narcs;		/* number of edges expected in graph */
    int       nedges;		/* twice number of edges really in graph */
    int       nedge;		/* loops through edges for each vertex */
    int       flag;		/* condition indicator */
    int       found_flag;	/* is vertex found in adjacency list? */
    int       skip_flag;	/* should this edge be ignored? */
    int       end_flag;		/* indicates end of line or file */
    int       vtx;		/* vertex in graph */
    int       line_num;		/* line number in input file */
    int       sum_edges;	/* total number of edges read so far */
    int       option = 0;	/* input option */
    int       using_ewgts;	/* are edge weights in input file? */
    int       using_vwgts;	/* are vertex weights in input file? */
    int       vtxnums;		/* are vertex numbers in input file? */
    int       vertex;		/* current vertex being read */
    int       new_vertex;	/* new vertex being read */
    float     weight;		/* weight being read */
    float     eweight;		/* edge weight being read */
    int       neighbor;		/* neighbor of current vertex */
    int       self_edge;	/* is a self edge encountered? */
    int       ignore_me;	/* is this edge being ignored? */
    int       ignored;		/* how many edges are ignored? */
    int       error_flag;	/* error reading input? */
    int       j;		/* loop counters */

    DEBUG_TRACE_START(0, yo);

    /* Read first line  of input (= nvtxs, narcs, option). */
    /* The (decimal) digits of the option variable mean: 1's digit not zero => input
       edge weights 10's digit not zero => input vertex weights 100's digit not zero
       => include vertex numbers */

    *start = NULL;
    *adjacency = NULL;
    *vweights = NULL;
    *eweights = NULL;

    error_flag = 0;
    line_num = 0;

    /* Read any leading comment lines */
    end_flag = 1;
    while (end_flag == 1) {
	*nvtxs = read_int(fin, &end_flag);
	++line_num;
    }
    if (*nvtxs <= 0) {
	printf("ERROR in graph file `%s':", inname);
	printf(" Invalid number of vertices (%d).\n", *nvtxs);
	ZOLTAN_FILE_close(fin);
	return(1);
    }

    narcs = read_int(fin, &end_flag);
    if (narcs < 0) {
	printf("ERROR in graph file `%s':", inname);
	printf(" Invalid number of expected edges (%d).\n", narcs);
	ZOLTAN_FILE_close(fin);
	return(1);
    }

    /*  Check if vertex or edge weights are used */
    if (!end_flag) {
	option = read_int(fin, &end_flag);
    }
    using_ewgts = option - 10 * (option / 10);
    option /= 10;
    using_vwgts = option - 10 * (option / 10);
    option /= 10;
    vtxnums = option - 10 * (option / 10);

    /* Get weight dimensions from Chaco option */
    (*vwgt_dim) = using_vwgts;
    (*ewgt_dim) = using_ewgts;

    /* Read weight dimensions if they are specified separately */
    if (!end_flag && using_vwgts==1){
       j = read_int(fin, &end_flag);
       if (!end_flag) (*vwgt_dim) = j;
    }
    if (!end_flag && using_ewgts==1){
       j = read_int(fin, &end_flag);
       if (!end_flag) (*ewgt_dim) = j;
    }

    /* Discard rest of line */
    while (!end_flag)
	j = read_int(fin, &end_flag);

    /* Allocate space for rows and columns. */
    *start = (int *) malloc((unsigned) (*nvtxs + 1) * sizeof(int));
    if (narcs != 0)
	*adjacency = (int *) malloc((unsigned) (2 * narcs + 1) * sizeof(int));
    else
	*adjacency = NULL;

    if (using_vwgts)
	*vweights = (float *) malloc((unsigned) (*nvtxs) * (*vwgt_dim) * sizeof(float));
    else
	*vweights = NULL;

    if (using_ewgts)
	*eweights = (float *)
                     malloc((unsigned) (2 * narcs + 1) * (*ewgt_dim) * sizeof(float));
    else
	*eweights = NULL;

    adjptr = *adjacency;
    ewptr = *eweights;
    self_edge = 0;
    ignored = 0;

    sum_edges = 0;
    nedges = 0;
    (*start)[0] = 0;
    vertex = 0;
    vtx = 0;
    new_vertex = TRUE;
    while ((using_vwgts || vtxnums || narcs) && end_flag != -1) {
	++line_num;

	/* If multiple input lines per vertex, read vertex number. */
	if (vtxnums) {
	    j = read_int(fin, &end_flag);
	    if (end_flag) {
		if (vertex == *nvtxs)
		    break;
		printf("ERROR in graph file `%s':", inname);
		printf(" no vertex number in line %d.\n", line_num);
		ZOLTAN_FILE_close(fin);
		return (1);
	    }
	    if (j != vertex && j != vertex + 1) {
		printf("ERROR in graph file `%s':", inname);
		printf(" out-of-order vertex number in line %d.\n", line_num);
		ZOLTAN_FILE_close(fin);
		return (1);
	    }
	    if (j != vertex) {
		new_vertex = TRUE;
		vertex = j;
	    }
	    else
		new_vertex = FALSE;
	}
	else
	    vertex = ++vtx;

	if (vertex > *nvtxs)
	    break;

	/* If vertices are weighted, read vertex weight. */
	if (using_vwgts && new_vertex) {
            for (j=0; j<(*vwgt_dim); j++){
	    	weight = read_val(fin, &end_flag);
	    	if (end_flag) {
			printf("ERROR in graph file `%s':", inname);
			printf(" not enough weights for vertex %d.\n", vertex);
			ZOLTAN_FILE_close(fin);
			return (1);
	    	}
	    	if ((weight < 0) && Debug_Chaco_Input) {
			printf("ERROR in graph file `%s':", inname);
			printf(" negative weight entered for vertex %d.\n", vertex);
			ZOLTAN_FILE_close(fin);
			return (1);
	    	}
	    	(*vweights)[(vertex-1)*(*vwgt_dim)+j] = weight;
	    }
	}

	nedge = 0;

	/* Read number of adjacent vertex. */
	neighbor = read_int(fin, &end_flag);

	while (!end_flag) {
	    skip_flag = FALSE;
	    ignore_me = FALSE;

            if (Debug_Chaco_Input){

    	    if (neighbor > *nvtxs) {
    		printf("ERROR in graph file `%s':", inname);
    		printf(" nvtxs=%d, but edge (%d,%d) was input.\n",
    		       *nvtxs, vertex, neighbor);
    		ZOLTAN_FILE_close(fin);
    		return (1);
    	    }
    	    if (neighbor < 0) {
    		printf("ERROR in graph file `%s':", inname);
    		printf(" negative vertex in edge (%d,%d).\n",
    		       vertex, neighbor);
    		ZOLTAN_FILE_close(fin);
    		return (1);
    	    }

    	    if (neighbor == vertex) {
    		if (!self_edge && Debug_Chaco_Input) {
    		    printf("WARNING: Self edge (%d, %d) being ignored.\n",
    			   vertex, vertex);
    		}
    		skip_flag = TRUE;
    		++self_edge;
    	    }

    	    /* Check if adjacency is repeated. */
    	    if (!skip_flag) {
    		found_flag = FALSE;
    		for (j = (*start)[vertex - 1]; !found_flag && j < sum_edges + nedge; j++) {
    		    if ((*adjacency)[j] == neighbor)
    			found_flag = TRUE;
    		}
    		if (found_flag) {
    		    printf("WARNING: Multiple occurences of edge (%d,%d) ignored\n",
    			vertex, neighbor);
    		    skip_flag = TRUE;
    		    if (!ignore_me) {
    			ignore_me = TRUE;
    			++ignored;
    		    }
    		}
    	    }
            }

	    if (using_ewgts) {	/* Read edge weight if it's being input. */
                for (j=0; j<(*ewgt_dim); j++){
		    eweight = read_val(fin, &end_flag);

		    if (end_flag) {
		        printf("ERROR in graph file `%s':", inname);
		        printf(" not enough weights for edge (%d,%d).\n", vertex, neighbor);
		        ZOLTAN_FILE_close(fin);
		        return (1);
		    }

		    if (eweight <= 0 && Debug_Chaco_Input) {
		        printf("WARNING: Bad weight entered for edge (%d,%d).  Edge ignored.\n",
			   vertex, neighbor);
		        skip_flag = TRUE;
		        if (!ignore_me) {
			    ignore_me = TRUE;
			    ++ignored;
		        }
		    }
		    else {
		        *ewptr++ = eweight;
		    }
		}
	    }

            if (Debug_Chaco_Input){

	        /* Check for edge only entered once. */
	        if (neighbor < vertex && !skip_flag) {
		    found_flag = FALSE;
		    for (j = (*start)[neighbor - 1]; !found_flag && j < (*start)[neighbor]; j++) {
		        if ((*adjacency)[j] == vertex)
			    found_flag = TRUE;
		    }
		    if (!found_flag) {
		        printf("ERROR in graph file `%s':", inname);
		        printf(" Edge (%d, %d) entered but not (%d, %d)\n",
			       vertex, neighbor, neighbor, vertex);
		        error_flag = 1;
		    }
	        }
	    }

	    /* Add edge to data structure. */
	    if (!skip_flag) {
		if (++nedges > 2*narcs) {
		    printf("ERROR in graph file `%s':", inname);
		    printf(" at least %d adjacencies entered, but nedges = %d\n",
			nedges, narcs);
		    ZOLTAN_FILE_close(fin);
		    return (1);
		}
		*adjptr++ = neighbor;
		nedge++;
	    }

	    /* Read number of next adjacent vertex. */
	    neighbor = read_int(fin, &end_flag);
	}

	sum_edges += nedge;
	(*start)[vertex] = sum_edges;
    }

    /* Make sure there's nothing else in file. */
    flag = FALSE;
    while (!flag && end_flag != -1) {
	read_int(fin, &end_flag);
	if (!end_flag)
	    flag = TRUE;
    }
    if (flag && Debug_Chaco_Input) {
	printf("WARNING: Possible error in graph file `%s'\n", inname);
	printf("         Data found after expected end of file\n");
    }

    (*start)[*nvtxs] = sum_edges;

    if (self_edge > 1 && Debug_Chaco_Input) {
	printf("WARNING: %d self edges were read and ignored.\n", self_edge);
    }

    if (vertex != 0) {		/* Normal file was read. */
        if (narcs) {
	    /* Make sure narcs was reasonable. */
	    if (nedges + 2 * self_edge != 2 * narcs &&
	        nedges + 2 * self_edge + ignored != 2 * narcs &&
	    	nedges + self_edge != 2 * narcs &&
		nedges + self_edge + ignored != 2 * narcs &&
		nedges != 2 * narcs &&
		nedges + ignored != 2 * narcs &&
		Debug_Chaco_Input) {
	        printf("WARNING: I expected %d edges entered twice, but I only count %d.\n",
	            narcs, nedges);
	    }
        }
        else { /* no edges, but did have vertex weights or vertex numbers */
  	    free(*start);
	    *start = NULL;
	    if (*adjacency != NULL)
	        free(*adjacency);
	    *adjacency = NULL;
	    if (*eweights != NULL)
	        free(*eweights);
            *eweights = NULL;
        }
    }

    else {
	/* Graph was empty */
	free(*start);
	if (*adjacency != NULL)
	    free(*adjacency);
	if (*vweights != NULL)
	    free(*vweights);
	if (*eweights != NULL)
	    free(*eweights);
	*start = NULL;
	*adjacency = NULL;
    }

    ZOLTAN_FILE_close(fin);

    DEBUG_TRACE_END(0, yo);
    return (error_flag);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
