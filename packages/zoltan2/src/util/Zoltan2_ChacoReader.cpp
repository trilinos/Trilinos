// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file ChacoReader.cpp
 *  \brief Reader for legacy zoltan1 files
 */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

namespace Zoltan2 {

#define LINE_LENGTH 200
static char chaco_line[LINE_LENGTH];	/* space to hold values */
static int chaco_offset = 0;		/* offset into line for next data */
static int chaco_break_pnt = LINE_LENGTH;	/* place in sequence to pause */
static int chaco_save_pnt;		/* place in sequence to save */
static void chaco_flush_line(FILE*);

/*! \brief read a double from a file
 *   Reader code was copied from zoltan1 test driver code.
 */

double    chaco_read_val(
  FILE* infile,		/* file to read value from */
  int      *end_flag 		/* 0 => OK, 1 => EOL, -1 => EOF */
)
{
    double    val;		/* return value */
    char     *ptr;		/* ptr to next string to read */
    char     *ptr2;		/* ptr to next string to read */
    int       length;		/* length of line to read */
    int       length_left;	/* length of line still around */
    int       white_seen;	/* have I detected white space yet? */
    int       done;		/* checking for end of scan */
    int       i;		/* loop counter */

    *end_flag = 0;

    if (chaco_offset == 0 || chaco_offset >= chaco_break_pnt) {
	if (chaco_offset >= chaco_break_pnt) { /* Copy rest of line back to beginning. */
	    length_left = LINE_LENGTH - chaco_save_pnt - 1;
	    ptr2 = chaco_line;
	    ptr = &chaco_line[chaco_save_pnt];
	    for (i=length_left; i; i--) *ptr2++ = *ptr++;
	    length = chaco_save_pnt + 1;
	}
	else {
	    length = LINE_LENGTH;
	    length_left = 0;
	}

	/* Now read next line, or next segment of current one. */
	ptr2 = fgets(&chaco_line[length_left], length, infile);

	if (ptr2 == (char *) NULL) {	/* We've hit end of file. */
	    *end_flag = -1;
	    return((double) 0.0);
	}

	if ((chaco_line[LINE_LENGTH - 2] != '\n') && (chaco_line[LINE_LENGTH - 2] != '\f')
	    && (strlen(chaco_line) == LINE_LENGTH - 1)){
	    /* Line too long.  Find last safe place in chaco_line. */
	    chaco_break_pnt = LINE_LENGTH - 1;
	    chaco_save_pnt = chaco_break_pnt;
	    white_seen = 0;
	    done = 0;
	    while (!done) {
		--chaco_break_pnt;
		if (chaco_line[chaco_break_pnt] != '\0') {
		    if (isspace((int)(chaco_line[chaco_break_pnt]))) {
			if (!white_seen) {
			    chaco_save_pnt = chaco_break_pnt + 1;
		            white_seen = 1;
			}
		    }
		    else if (white_seen) {
		        done= 1;
		    }
		}
	    }
	}
	else {
	    chaco_break_pnt = LINE_LENGTH;
	}

	chaco_offset = 0;
    }

    while (isspace((int)(chaco_line[chaco_offset])) && chaco_offset < LINE_LENGTH) chaco_offset++;
    if (chaco_line[chaco_offset] == '%' || chaco_line[chaco_offset] == '#') {
	*end_flag = 1;
	if (chaco_break_pnt < LINE_LENGTH) {
	    chaco_flush_line(infile);
	}
	return((double) 0.0);
    }

    ptr = &(chaco_line[chaco_offset]);
    val = strtod(ptr, &ptr2);

    if (ptr2 == ptr) {	/* End of input line. */
	chaco_offset = 0;
	*end_flag = 1;
	return((double) 0.0);
    }
    else {
	chaco_offset = (int) (ptr2 - chaco_line) / sizeof(char);
    }

    return(val);
}


/*! \brief read a int from a file
 *   Reader code was copied from zoltan1 test driver code.
 */

int     chaco_read_int(
FILE *infile,		/* file to read value from */
int      *end_flag 		/* 0 => OK, 1 => EOL, -1 => EOF */
)
{
    int       val;		/* return value */
    char     *ptr;		/* ptr to next string to read */
    char     *ptr2;		/* ptr to next string to read */
    int       length;		/* length of line to read */
    int       length_left;	/* length of line still around */
    int       white_seen;	/* have I detected white space yet? */
    int       done;		/* checking for end of scan */
    int       i;		/* loop counter */

    *end_flag = 0;

    if (chaco_offset == 0 || chaco_offset >= chaco_break_pnt) {
	if (chaco_offset >= chaco_break_pnt) { /* Copy rest of line back to beginning. */
	    length_left = LINE_LENGTH - chaco_save_pnt - 1;
	    ptr2 = chaco_line;
	    ptr = &chaco_line[chaco_save_pnt];
	    for (i=length_left; i; i--) *ptr2++ = *ptr++;
	    length = chaco_save_pnt + 1;
	}
	else {
	    length = LINE_LENGTH;
	    length_left = 0;
	}

	/* Now read next line, or next segment of current one. */
	ptr2 = fgets(&chaco_line[length_left], length, infile);

	if (ptr2 == (char *) NULL) {	/* We've hit end of file. */
	    *end_flag = -1;
	    return(0);
	}

	if ((chaco_line[LINE_LENGTH - 2] != '\n') && (chaco_line[LINE_LENGTH - 2] != '\f')
	    && (strlen(chaco_line) == LINE_LENGTH - 1)){
	    /* Line too long.  Find last safe place in line. */
	    chaco_break_pnt = LINE_LENGTH - 1;
	    chaco_save_pnt = chaco_break_pnt;
	    white_seen = 0;
	    done = 0;
	    while (!done) {
		--chaco_break_pnt;
		if (chaco_line[chaco_break_pnt] != '\0') {
		    if (isspace((int)(chaco_line[chaco_break_pnt]))) {
			if (!white_seen) {
			    chaco_save_pnt = chaco_break_pnt + 1;
		            white_seen = 1;
			}
		    }
		    else if (white_seen) {
		        done= 1;
		    }
		}
	    }
	}
	else {
	    chaco_break_pnt = LINE_LENGTH;
	}

	chaco_offset = 0;
    }

    while (isspace((int)(chaco_line[chaco_offset])) && chaco_offset < LINE_LENGTH) chaco_offset++;
    if (chaco_line[chaco_offset] == '%' || chaco_line[chaco_offset] == '#') {
	*end_flag = 1;
	if (chaco_break_pnt < LINE_LENGTH) {
	    chaco_flush_line(infile);
	}
	return(0);
    }

    ptr = &(chaco_line[chaco_offset]);
    val = (int) strtol(ptr, &ptr2, 10);

    if (ptr2 == ptr) {	/* End of input chaco_line. */
	chaco_offset = 0;
	*end_flag = 1;
	return(0);
    }
    else {
	chaco_offset = (int) (ptr2 - chaco_line) / sizeof(char);
    }

    return(val);
}


/*! \brief helper function for zoltan1 chaco reader code
 *   Reader code was copied from zoltan1 test driver code.
 */

void chaco_flush_line(
FILE *infile 		/* file to read value from */
)
{
    char      c;		/* character being read */

    c = fgetc(infile);
    while (c != '\n' && c != '\f')
	c = fgetc(infile);
}

/*! \brief read a Chaco graph file
 *   Reader code was copied from zoltan1 test driver code.
 *   File is closed when read is completed.
 *   \return 0 on success, 1 on failure
 */

int chaco_input_graph(
FILE *fin,			/* input file */
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
  int      *adjptr;		/* loops through adjacency data */
  float    *ewptr;		/* loops through edge weight data */
  int       narcs;		/* number of edges expected in graph */
  int       nedges;		/* twice number of edges really in graph */
  int       nedge;		/* loops through edges for each vertex */
  int       flag;		/* condition indicator */
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
	*nvtxs = chaco_read_int(fin, &end_flag);
	++line_num;
  }
  if (*nvtxs <= 0) {
	printf("ERROR in graph file `%s':", inname);
	printf(" Invalid number of vertices (%d).\n", *nvtxs);
	fclose(fin);
	return(1);
  }

  narcs = chaco_read_int(fin, &end_flag);
  if (narcs < 0) {
	printf("ERROR in graph file `%s':", inname);
	printf(" Invalid number of expected edges (%d).\n", narcs);
	fclose(fin);
	return(1);
  }

  /*  Check if vertex or edge weights are used */
  if (!end_flag) {
	option = chaco_read_int(fin, &end_flag);
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
     j = chaco_read_int(fin, &end_flag);
     if (!end_flag) (*vwgt_dim) = j;
  }
  if (!end_flag && using_ewgts==1){
     j = chaco_read_int(fin, &end_flag);
     if (!end_flag) (*ewgt_dim) = j;
  }

  /* Discard rest of line */
  while (!end_flag)
	j = chaco_read_int(fin, &end_flag);

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
  new_vertex = 1;
  while ((using_vwgts || vtxnums || narcs) && end_flag != -1) {
	++line_num;

	/* If multiple input lines per vertex, read vertex number. */
	if (vtxnums) {
	  j = chaco_read_int(fin, &end_flag);
	  if (end_flag) {
		if (vertex == *nvtxs)
		  break;
		printf("ERROR in graph file `%s':", inname);
		printf(" no vertex number in line %d.\n", line_num);
		fclose(fin);
		return (1);
	  }
	  if (j != vertex && j != vertex + 1) {
		printf("ERROR in graph file `%s':", inname);
		printf(" out-of-order vertex number in line %d.\n", line_num);
		fclose(fin);
		return (1);
	  }
	  if (j != vertex) {
		new_vertex = 1;
		vertex = j;
	  }
	  else
		new_vertex = 0;
	}
	else
	  vertex = ++vtx;

	if (vertex > *nvtxs)
	  break;

	/* If vertices are weighted, read vertex weight. */
	if (using_vwgts && new_vertex) {
          for (j=0; j<(*vwgt_dim); j++){
	  	weight = chaco_read_val(fin, &end_flag);
	  	if (end_flag) {
			printf("ERROR in graph file `%s':", inname);
			printf(" not enough weights for vertex %d.\n", vertex);
			fclose(fin);
			return (1);
	  	}
	  	(*vweights)[(vertex-1)*(*vwgt_dim)+j] = weight;
	  }
	}

	nedge = 0;

	/* Read number of adjacent vertex. */
	neighbor = chaco_read_int(fin, &end_flag);

	while (!end_flag) {
	  skip_flag = 0;
	  ignore_me = 0;

	  if (using_ewgts) {	/* Read edge weight if it's being input. */
              for (j=0; j<(*ewgt_dim); j++){
		  eweight = chaco_read_val(fin, &end_flag);

		  if (end_flag) {
		      printf("ERROR in graph file `%s':", inname);
		      printf(" not enough weights for edge (%d,%d).\n", vertex, neighbor);
		      fclose(fin);
		      return (1);
		  }

		  else {
		      *ewptr++ = eweight;
		  }
		}
	  }

	  /* Add edge to data structure. */
	  if (!skip_flag) {
		if (++nedges > 2*narcs) {
		  printf("ERROR in graph file `%s':", inname);
		  printf(" at least %d adjacencies entered, but nedges = %d\n",
			nedges, narcs);
		  fclose(fin);
		  return (1);
		}
		*adjptr++ = neighbor;
		nedge++;
	  }

	  /* Read number of next adjacent vertex. */
	  neighbor = chaco_read_int(fin, &end_flag);
	}

	sum_edges += nedge;
	(*start)[vertex] = sum_edges;
  }

  /* Make sure there's nothing else in file. */
  flag = 0;
  while (!flag && end_flag != -1) {
	chaco_read_int(fin, &end_flag);
	if (!end_flag)
	  flag = 1;
  }

  (*start)[*nvtxs] = sum_edges;

  if (vertex != 0) {		/* Normal file was read. */
      if (narcs) {
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

  fclose(fin);

  return (error_flag);
}


/*! \brief read a Chaco coordinates file
 *   Reader code was copied from zoltan1 test driver code.
 *   File is closed when read is completed.
 *   \return 0 on success, 1 on failure
 */

int chaco_input_geom(
FILE *fingeom,		/* geometry input file */
char     *geomname,		/* name of geometry file */
int       nvtxs,		/* number of coordinates to read */
int      *igeom,		/* dimensionality of geometry */
float   **x,         		/* coordinates of vertices */
float   **y,
float   **z
)
{
    float     xc, yc, zc =0;	/* first x, y, z coordinate */
    int       nread;		/* number of lines of coordinates read */
    int       flag;		/* any bad data at end of file? */
    int       line_num;		/* counts input lines in file */
    int       end_flag;		/* return conditional */
    int       ndims;		/* number of values in an input line */
    int       i=0;		/* loop counter */

    *x = *y = *z = NULL;
    line_num = 0;
    end_flag = 1;
    while (end_flag == 1) {
	xc = chaco_read_val(fingeom, &end_flag);
	++line_num;
    }

    if (end_flag == -1) {
	printf("No values found in geometry file `%s'\n", geomname);
	fclose(fingeom);
	return (1);
    }

    ndims = 1;
    yc = chaco_read_val(fingeom, &end_flag);
    if (end_flag == 0) {
	ndims = 2;
	zc = chaco_read_val(fingeom, &end_flag);
	if (end_flag == 0) {
	    ndims = 3;
	    chaco_read_val(fingeom, &end_flag);
	    if (!end_flag) {
		printf("Too many values on input line of geometry file `%s'\n",
		       geomname);

		printf(" Maximum dimensionality is 3\n");
		fclose(fingeom);
		return (1);
	    }
	}
    }

    *igeom = ndims;

    *x = (float *) malloc((unsigned) nvtxs * sizeof(float));
    (*x)[0] = xc;
    if (ndims > 1) {
	*y = (float *) malloc((unsigned) nvtxs * sizeof(float));
	(*y)[0] = yc;
    }
    if (ndims > 2) {
	*z = (float *) malloc((unsigned) nvtxs * sizeof(float));
	(*z)[0] = zc;
    }

    for (nread = 1; nread < nvtxs; nread++) {
	++line_num;
	if (ndims == 1) {
	    i = fscanf(fingeom, "%f", &((*x)[nread]));
	}
	else if (ndims == 2) {
	    i = fscanf(fingeom, "%f%f", &((*x)[nread]), &((*y)[nread]));
	}
	else if (ndims == 3) {
	    i = fscanf(fingeom, "%f%f%f", &((*x)[nread]), &((*y)[nread]),
		       &((*z)[nread]));
	}

	if (i == EOF) {
	    printf("Too few lines of values in geometry file; nvtxs=%d, but only %d read\n",
		   nvtxs, nread);
	    fclose(fingeom);
	    return (1);
	}
	else if (i != ndims) {
	    printf("Wrong number of values in line %d of geometry file `%s'\n",
		   line_num, geomname);
	    fclose(fingeom);
	    return (1);
	}
    }

    /* Check for spurious extra stuff in file. */
    flag = 0;
    end_flag = 0;
    while (!flag && end_flag != -1) {
	chaco_read_val(fingeom, &end_flag);
	if (!end_flag)
	    flag = 1;
    }

    fclose(fingeom);

    return (0);
}

}  // namespace Zoltan2
