/*
 * Copyright (c) 2014, Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Governement retains certain rights in this software.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>                      // for printf, fprintf, NULL, etc
#include "defs.h"                       // for TRUE, FALSE
#include "smalloc.h"                    // for smalloc, sfree


int       input_graph(FILE *fin, char *inname, int **start, int **adjacency, int *nvtxs, int **vweights, float **eweights)

/*
 *  FILE     *fin;	 input file
 *  char     *inname;	 name of input file
 *  int     **start;	 start of edge list for each vertex
 *  int     **adjacency; edge list data
 *  int      *nvtxs;	 number of vertices in graph 
 *  int     **vweights;	 vertex weight list data 
 *  float   **eweights;	 edge weight list data
 */
{
    extern FILE *Output_File;   /* output file or null */
    extern int CHECK_INPUT;	/* print warnings or not? */
    extern int DEBUG_INPUT;	/* echo that input file read successful? */
    extern int DEBUG_TRACE;	/* trace main execution path */
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
    int       option=0;		/* input option */
    int       using_ewgts;	/* are edge weights in input file? */
    int       using_vwgts;	/* are vertex weights in input file? */
    int       vtxnums;		/* are vertex numbers in input file? */
    int       vertex;		/* current vertex being read */
    int       new_vertex;	/* new vertex being read */
    int       weight;		/* weight being read */
    float     eweight;		/* edge weight being read */
    int       neighbor;		/* neighbor of current vertex */
    int       self_edge;	/* is a self edge encountered? */
    int       ignore_me;	/* is this edge being ignored? */
    int       ignored;		/* how many edges are ignored? */
    int       error_flag;	/* error reading input? */
    int       j;		/* loop counters */
    int       read_int();
    double    read_val();

    if (DEBUG_TRACE > 0) {
	printf("<Entering input_graph>\n");
    }

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
	if (Output_File != NULL) {
	    fprintf(Output_File, "ERROR in graph file `%s':", inname);
	    fprintf(Output_File, " Invalid number of vertices (%d).\n", *nvtxs);
	}
	fclose(fin);
	return(1);
    }

    narcs = read_int(fin, &end_flag);
    if (narcs < 0) {
	printf("ERROR in graph file `%s':", inname);
	printf(" Invalid number of expected edges (%d).\n", narcs);
	if (Output_File != NULL) {
	    fprintf(Output_File, "ERROR in graph file `%s':", inname);
	    fprintf(Output_File, " Invalid number of expected edges (%d).\n", narcs);
	}
	fclose(fin);
	return(1);
    }

    if (!end_flag) {
	option = read_int(fin, &end_flag);
    }
    while (!end_flag)
	j = read_int(fin, &end_flag);

    using_ewgts = option - 10 * (option / 10);
    option /= 10;
    using_vwgts = option - 10 * (option / 10);
    option /= 10;
    vtxnums = option - 10 * (option / 10);

    /* Allocate space for rows and columns. */
    *start = smalloc((*nvtxs + 1) * sizeof(int));
    if (narcs != 0)
	*adjacency = smalloc((2 * narcs + 1) * sizeof(int));
    else
	*adjacency = NULL;

    if (using_vwgts)
	*vweights = smalloc((*nvtxs) * sizeof(int));
    else
	*vweights = NULL;

    if (using_ewgts && narcs != 0)
	*eweights = smalloc((2 * narcs + 1) * sizeof(float));
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
    while ((narcs || using_vwgts || vtxnums) && end_flag != -1) {
	++line_num;

	/* If multiple input lines per vertex, read vertex number. */
	if (vtxnums) {
	    j = read_int(fin, &end_flag);
	    if (end_flag) {
		if (vertex == *nvtxs)
		    break;
		printf("ERROR in graph file `%s':", inname);
		printf(" no vertex number in line %d.\n", line_num);
		if (Output_File != NULL) {
		    fprintf(Output_File, "ERROR in graph file `%s':", inname);
		    fprintf(Output_File, " no vertex number in line %d.\n", line_num);
		}
		fclose(fin);
		return (1);
	    }
	    if (j != vertex && j != vertex + 1) {
		printf("ERROR in graph file `%s':", inname);
		printf(" out-of-order vertex number in line %d.\n", line_num);
		if (Output_File != NULL) {
		    fprintf(Output_File, "ERROR in graph file `%s':", inname);
		    fprintf(Output_File,
			" out-of-order vertex number in line %d.\n", line_num);
		}
		fclose(fin);
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
	    weight = read_int(fin, &end_flag);
	    if (end_flag) {
		printf("ERROR in graph file `%s':", inname);
		printf(" no weight for vertex %d.\n", vertex);
		if (Output_File != NULL) {
		    fprintf(Output_File, "ERROR in graph file `%s':", inname);
		    fprintf(Output_File, " no weight for vertex %d.\n", vertex);
		}
		fclose(fin);
		return (1);
	    }
	    if (weight <= 0) {
		printf("ERROR in graph file `%s':", inname);
		printf(" zero or negative weight entered for vertex %d.\n", vertex);
		if (Output_File != NULL) {
		    fprintf(Output_File, "ERROR in graph file `%s':", inname);
		    fprintf(Output_File,
			" zero or negative weight entered for vertex %d.\n", vertex);
		}
		fclose(fin);
		return (1);
	    }
	    (*vweights)[vertex-1] = weight;
	}

	nedge = 0;

	/* Read number of adjacent vertex. */
	neighbor = read_int(fin, &end_flag);

	while (!end_flag) {
	    skip_flag = FALSE;
	    ignore_me = FALSE;

	    if (neighbor > *nvtxs) {
		printf("ERROR in graph file `%s':", inname);
		printf(" nvtxs=%d, but edge (%d,%d) was input.\n",
		       *nvtxs, vertex, neighbor);
		if (Output_File != NULL) {
		    fprintf(Output_File, "ERROR in graph file `%s':", inname);
		    fprintf(Output_File,
			" nvtxs=%d, but edge (%d,%d) was input.\n",
			*nvtxs, vertex, neighbor);
		}
		fclose(fin);
		return (1);
	    }
	    if (neighbor <= 0) {
		printf("ERROR in graph file `%s':", inname);
		printf(" zero or negative vertex in edge (%d,%d).\n",
		       vertex, neighbor);
		if (Output_File != NULL) {
		    fprintf(Output_File, "ERROR in graph file `%s':", inname);
		    fprintf(Output_File,
		        " zero or negative vertex in edge (%d,%d).\n",
			vertex, neighbor);
		}
		fclose(fin);
		return (1);
	    }

	    if (neighbor == vertex) {
		if (!self_edge && CHECK_INPUT) {
		    printf("WARNING: Self edge (%d, %d) being ignored.\n",
			   vertex, vertex);
    		    if (Output_File != NULL) {
		        fprintf(Output_File,
			    "WARNING: Self edge (%d, %d) being ignored.\n", vertex, vertex);
		    }
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
		    if (Output_File != NULL) {
			fprintf(Output_File,
			    "WARNING: Multiple occurences of edge (%d,%d) ignored\n",
			    vertex, neighbor);
		    }
		    skip_flag = TRUE;
		    if (!ignore_me) {
			ignore_me = TRUE;
			++ignored;
		    }
		}
	    }

	    if (using_ewgts) {	/* Read edge weight if it's being input. */
		eweight = read_val(fin, &end_flag);

		if (end_flag) {
		    printf("ERROR in graph file `%s':", inname);
		    printf(" no weight for edge (%d,%d).\n", vertex, neighbor);
		    if (Output_File != NULL) {
		        fprintf(Output_File, "ERROR in graph file `%s':", inname);
		        fprintf(Output_File,
			    " no weight for edge (%d,%d).\n", vertex, neighbor);
		    }
		    fclose(fin);
		    return (1);
		}

		if (eweight <= 0 && CHECK_INPUT) {
		    printf("WARNING: Bad weight entered for edge (%d,%d).  Edge ignored.\n",
			   vertex, neighbor);
    		    if (Output_File != NULL) {
		        fprintf(Output_File,
			    "WARNING: Bad weight entered for edge (%d,%d).  Edge ignored.\n",
			    vertex, neighbor);
		    }
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
		    if (Output_File != NULL) {
		        fprintf(Output_File, "ERROR in graph file `%s':", inname);
		        fprintf(Output_File,
			    " Edge (%d, %d) entered but not (%d, %d)\n",
			    vertex, neighbor, neighbor, vertex);
		    }
		    error_flag = 1;
		}
	    }

	    /* Add edge to data structure. */
	    if (!skip_flag) {
		if (++nedges > 2*narcs) {
		    printf("ERROR in graph file `%s':", inname);
		    printf(" at least %d adjacencies entered, but nedges = %d\n",
			nedges, narcs);
		    if (Output_File != NULL) {
		        fprintf(Output_File, "ERROR in graph file `%s':", inname);
		        fprintf(Output_File,
			    " at least %d adjacencies entered, but nedges = %d\n",
			    nedges, narcs);
		    }
		    fclose(fin);
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
    if (flag && CHECK_INPUT) {
	printf("WARNING: Possible error in graph file `%s'\n", inname);
	printf("         Data found after expected end of file\n");
        if (Output_File != NULL) {
	    fprintf(Output_File,
		"WARNING: Possible error in graph file `%s'\n", inname);
	    fprintf(Output_File,
		"         Data found after expected end of file\n");
	}
    }

    (*start)[*nvtxs] = sum_edges;

    if (self_edge > 1 && CHECK_INPUT) {
	printf("WARNING: %d self edges were read and ignored.\n", self_edge);
        if (Output_File != NULL) {
	    fprintf(Output_File,
		"WARNING: %d self edges were read and ignored.\n", self_edge);
	}
    }

    if (vertex != 0) {		/* Normal file was read. */
	/* Make sure narcs was reasonable. */
	if (nedges + 2 * self_edge != 2 * narcs &&
	    nedges + 2 * self_edge + ignored != 2 * narcs &&
		nedges + self_edge != 2 * narcs && 
		nedges + self_edge + ignored != 2 * narcs && 
		nedges != 2 * narcs &&
		nedges + ignored != 2 * narcs &&
		CHECK_INPUT) {
	    printf("WARNING: I expected %d edges entered twice, but I only count %d.\n",
	        narcs, nedges);
    	    if (Output_File != NULL) {
	        fprintf(Output_File,
		    "WARNING: I expected %d edges entered twice, but I only count %d.\n",
	            narcs, nedges);
	    }
	}
    }

    else {
	/* Graph was empty => must be using inertial method. */
	sfree(*start);
	if (*adjacency != NULL)
	    sfree(*adjacency);
	if (*eweights != NULL)
	    sfree(*eweights);
	*start = NULL;
	*adjacency = NULL;
    }

    fclose(fin);

    if (DEBUG_INPUT > 0) {
	printf("Done reading graph file `%s'.\n", inname);
    }
    return (error_flag);
}
