/*
 * Copyright (c) 2005 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
 * retains certain rights in this software.
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
/*****************************************************************************
*
* testrd - read exodus file test.exo created by testwt
*
* author - Sandia National Laboratories
*          Larry A. Schoof - Original
*
*          
* environment - UNIX
*
* entry conditions - 
*   input parameters:
*       int     exoid                   exodus file id
*
* exit conditions - 
*
* revision history - 
*
*
*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "netcdf.h"
#include "exodusII.h"
/* #include "drmd.h" */


int main (int argc, char **argv)
{
  int exoid, num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets;
  int num_side_sets, error;
  int i, j, k, node_ctr;
  int *elem_map, *connect, *node_list, *node_ctr_list, *elem_list, *side_list;
  int *ids; 
  int *num_nodes_per_set = NULL;
  int *num_elem_per_set = NULL;
  int *num_df_per_set = NULL;
  int *node_ind = NULL;
  int *elem_ind = NULL;
  int *df_ind = NULL;
  int num_qa_rec, num_info;
  int num_glo_vars, num_nod_vars, num_ele_vars;
  int num_nset_vars, num_sset_vars;
  int *truth_tab;
  int num_time_steps;
  int *num_elem_in_block = NULL;
  int *num_nodes_per_elem = NULL;
  int *num_attr = NULL;
  int num_nodes_in_set, num_elem_in_set;
  int num_sides_in_set, num_df_in_set;
  int list_len, elem_list_len, node_list_len, df_list_len;
  int node_num, time_step, var_index, beg_time, end_time, elem_num;
  int CPU_word_size,IO_word_size;
  int num_props, prop_value, *prop_values;
  int idum;

  float time_value, *time_values, *var_values;
  float *x, *y, *z;
  float *attrib, *dist_fact;
  float version, fdum;

  char *coord_names[3], *qa_record[2][4], *info[3], *var_names[3];
  char *block_names[10], *nset_names[10], *sset_names[10];
  char *attrib_names[10];
  char name[MAX_STR_LENGTH+1];
  char title[MAX_LINE_LENGTH+1], elem_type[MAX_STR_LENGTH+1];
  char title_chk[MAX_LINE_LENGTH+1];
  char *cdum = 0;
  char *prop_names[3];

  CPU_word_size = 0;                    /* sizeof(float) */
  IO_word_size = 0;                     /* use what is stored in file */

  ex_opts (EX_VERBOSE | EX_ABORT );

  /* open EXODUS II files */
  exoid = ex_open ("test.exo",  /* filename path */
                   EX_READ,             /* access mode = READ */
                   &CPU_word_size,      /* CPU word size */
                   &IO_word_size,       /* IO word size */
                   &version);           /* ExodusII library version */

  printf ("\nafter ex_open\n");
  if (exoid < 0) exit(1);

  printf ("test.exo is an EXODUSII file; version %4.2f\n",
          version);
  /*   printf ("         CPU word size %1d\n",CPU_word_size);  */
  printf ("         I/O word size %1d\n",IO_word_size);
  ex_inquire(exoid,EX_INQ_API_VERS, &idum, &version, cdum);
  printf ("EXODUSII API; version %4.2f\n", version);

  ex_inquire(exoid,EX_INQ_LIB_VERS, &idum, &version, cdum);
  printf ("EXODUSII Library API; version %4.2f (%d)\n", version, idum);

  /* ncopts = NC_VERBOSE; */

  /* read database parameters */

  error = ex_get_init (exoid, title, &num_dim, &num_nodes, &num_elem,
                       &num_elem_blk, &num_node_sets, &num_side_sets);

  printf ("after ex_get_init, error = %3d\n", error);

  printf ("database parameters:\n");
  printf ("title =  '%s'\n",title);
  printf ("num_dim = %3d\n",num_dim);
  printf ("num_nodes = %3d\n",num_nodes);
  printf ("num_elem = %3d\n",num_elem);
  printf ("num_elem_blk = %3d\n",num_elem_blk);
  printf ("num_node_sets = %3d\n",num_node_sets);
  printf ("num_side_sets = %3d\n",num_side_sets);

  /* Check that ex_inquire gives same title */
  error = ex_inquire (exoid, EX_INQ_TITLE, &idum, &fdum, title_chk);
  if (strcmp(title, title_chk) != 0) {
    printf ("error in ex_inquire for EX_INQ_TITLE\n");
  }
  
  /* read nodal coordinates values and names from database */

  x = (float *) calloc(num_nodes, sizeof(float));
  y = (float *) calloc(num_nodes, sizeof(float));
  if (num_dim >= 3)
    z = (float *) calloc(num_nodes, sizeof(float));
  else
    z = 0;

  error = ex_get_coord (exoid, x, y, z);
  printf ("\nafter ex_get_coord, error = %3d\n", error);

  printf ("x coords = \n");
  for (i=0; i<num_nodes; i++)
    {
      printf ("%5.1f\n", x[i]);
    }

  printf ("y coords = \n");
  for (i=0; i<num_nodes; i++)
    {
      printf ("%5.1f\n", y[i]);
    }

  if (num_dim >= 3)
    {
      printf ("z coords = \n");
      for (i=0; i<num_nodes; i++)
        {
          printf ("%5.1f\n", z[i]);
        }
    }

  /*
    error = ex_get_1_coord (exoid, 2, x, y, z);
    printf ("\nafter ex_get_1_coord, error = %3d\n", error);

    printf ("x coord of node 2 = \n");
    printf ("%f \n", x[0]);

    printf ("y coord of node 2 = \n");
    printf ("%f \n", y[0]);
  */
  free (x);
  free (y);
  if (num_dim >= 3)
    free (z);


  for (i=0; i<num_dim; i++)
    {
      coord_names[i] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
    }

  error = ex_get_coord_names (exoid, coord_names);
  printf ("\nafter ex_get_coord_names, error = %3d\n", error);

  printf ("x coord name = '%s'\n", coord_names[0]);
  printf ("y coord name = '%s'\n", coord_names[1]);

  for (i=0; i<num_dim; i++)
    free(coord_names[i]);

  {
    int num_attrs = 0;
    error = ex_get_attr_param(exoid, EX_NODAL, 0, &num_attrs);
    printf ("num nodal attributes = %d\n", num_attrs);
    if (num_attrs > 0) {
      for (j=0; j<num_attrs; j++) {
	attrib_names[j] = (char *)calloc ((MAX_STR_LENGTH+1), sizeof(char));
      }
      error = ex_get_attr_names (exoid, EX_NODAL, 0, attrib_names);
      printf (" after ex_get_attr_names, error = %d\n", error);
      
      if (error == 0) {
	attrib = (float *) calloc(num_nodes,sizeof(float));
	for (j=0; j<num_attrs; j++) {
	  printf ("nodal attribute %d = '%s'\n", j, attrib_names[j]);
	  error = ex_get_one_attr(exoid, EX_NODAL, 0, j+1, attrib);
	  for (i=0; i < num_nodes; i++) {
	    printf ("%5.1f\n", attrib[i]);
	  }
	  free(attrib_names[j]);
	}
	free(attrib);
      }
    }
  }
  
  /* read element order map */

  elem_map = (int *) calloc(num_elem, sizeof(int));

  error = ex_get_map (exoid, elem_map);
  printf ("\nafter ex_get_map, error = %3d\n", error);

  for (i=0; i<num_elem; i++)
    {
      printf ("elem_map(%d) = %d \n", i, elem_map[i]);
    }

  free (elem_map);

  /* read element block parameters */

  if (num_elem_blk > 0) {
    ids = (int *) calloc(num_elem_blk, sizeof(int));
    num_elem_in_block = (int *) calloc(num_elem_blk, sizeof(int));
    num_nodes_per_elem = (int *) calloc(num_elem_blk, sizeof(int));
    num_attr = (int *) calloc(num_elem_blk, sizeof(int));
     
    error = ex_get_elem_blk_ids (exoid, ids);
    printf ("\nafter ex_get_elem_blk_ids, error = %3d\n", error);
     
    for (i=0; i<num_elem_blk; i++) {
      block_names[i] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
    }

    error = ex_get_names(exoid, EX_ELEM_BLOCK, block_names);
    printf ("\nafter ex_get_names, error = %3d\n", error);
    
    for (i=0; i<num_elem_blk; i++)
      {
	ex_get_name(exoid, EX_ELEM_BLOCK, ids[i], name);
	if (strcmp(name, block_names[i]) != 0) {
	  printf ("error in ex_get_name for block id %d\n", ids[i]);
	}
        error = ex_get_elem_block (exoid, ids[i], elem_type,
                                   &(num_elem_in_block[i]), 
                                   &(num_nodes_per_elem[i]), &(num_attr[i]));
        printf ("\nafter ex_get_elem_block, error = %d\n", error);
         
        printf ("element block id = %2d\n",ids[i]);
        printf ("element type = '%s'\n", elem_type);
        printf ("num_elem_in_block = %2d\n",num_elem_in_block[i]);
        printf ("num_nodes_per_elem = %2d\n",num_nodes_per_elem[i]);
        printf ("num_attr = %2d\n",num_attr[i]);
        printf ("name = '%s'\n",block_names[i]);
	free(block_names[i]);
      }
     
    /* read element block properties */
    error = ex_inquire (exoid, EX_INQ_EB_PROP, &num_props, &fdum, cdum);
    printf ("\nafter ex_inquire, error = %d\n", error);
    printf ("\nThere are %2d properties for each element block\n", num_props);
     
    for (i=0; i<num_props; i++)
      {
        prop_names[i] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
      }
     
    error = ex_get_prop_names(exoid,EX_ELEM_BLOCK,prop_names);
    printf ("after ex_get_prop_names, error = %d\n", error);
     
     
    for (i=1; i<num_props; i++) /* Prop 1 is id; skip that here */
      {
        for (j=0; j<num_elem_blk; j++)
          {
            error = ex_get_prop(exoid, EX_ELEM_BLOCK, ids[j], prop_names[i],
                                &prop_value);
            if (error == 0)
              printf ("element block %2d, property(%2d): '%s'= %5d\n",
                      j+1, i+1, prop_names[i], prop_value);
            else
              printf ("after ex_get_prop, error = %d\n", error);
          }
      }
     
    for (i=0; i<num_props; i++)
      free(prop_names[i]);
  }
   
  /* read element connectivity */

  for (i=0; i<num_elem_blk; i++)
    {
      if (num_elem_in_block[i] > 0) {
	connect = (int *) calloc((num_nodes_per_elem[i] * num_elem_in_block[i]), 
				 sizeof(int));
	
	error = ex_get_elem_conn (exoid, ids[i], connect);
	printf ("\nafter ex_get_elem_conn, error = %d\n", error);
	
	
	printf ("connect array for elem block %2d\n", ids[i]);
	
	for (j=0; j<num_nodes_per_elem[i]; j++)
	  {
	    printf ("%3d\n", connect[j]);
	  }
	/*
	  error = ex_get_1_elem_conn (exoid, 1, ids[i], connect);
	  printf ("\nafter ex_get_elem_conn, error = %d\n", error);
	  
	  printf ("node list for first element of element block %d \n ", ids[i]);
	  for (j=0; j<num_nodes_per_elem[i]; j++)
	  {
	  printf ("%d \n", connect[j]);
	  }
	*/
	free (connect);
      }
    }

  /* read element block attributes */

  for (i=0; i<num_elem_blk; i++)
    {
      if (num_elem_in_block[i] > 0) {
	for (j=0; j<num_attr[i]; j++)
	  attrib_names[j] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
	
	attrib = (float *) calloc(num_attr[i]*num_elem_in_block[i],sizeof(float));
	error = ex_get_elem_attr (exoid, ids[i], attrib);
	printf ("\n after ex_get_elem_attr, error = %d\n", error);
	
	if (error == 0) {
	  error = ex_get_elem_attr_names (exoid, ids[i], attrib_names);
	  printf (" after ex_get_elem_attr_names, error = %d\n", error);
	  
	  if (error == 0) {
	    printf ("element block %d attribute '%s' = %6.4f\n", ids[i], attrib_names[0], *attrib);
	  }
	}
	free (attrib);
	for (j=0; j<num_attr[i]; j++)
	  free (attrib_names[j]);
      }
    }
      
  if (num_elem_blk > 0) {
    free (ids);
    free (num_nodes_per_elem);
    free (num_attr);
  }

  /* read individual node sets */
  if (num_node_sets > 0) {
    ids = (int *) calloc(num_node_sets, sizeof(int));

    error = ex_get_node_set_ids (exoid, ids);
    printf ("\nafter ex_get_node_set_ids, error = %3d\n", error);

    for (i=0; i<num_node_sets; i++) {
      nset_names[i] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
    }

    error = ex_get_names(exoid, EX_NODE_SET, nset_names);
    printf ("\nafter ex_get_names, error = %3d\n", error);

    for (i=0; i<num_node_sets; i++)
      {
	ex_get_name(exoid, EX_NODE_SET, ids[i], name);
	if (strcmp(name, nset_names[i]) != 0) {
	  printf ("error in ex_get_name for nodeset id %d\n", ids[i]);
	}

        error = ex_get_node_set_param (exoid, ids[i], 
                                       &num_nodes_in_set, &num_df_in_set);
        printf ("\nafter ex_get_node_set_param, error = %3d\n", error);

        printf ("\nnode set %2d parameters: \n", ids[i]);
        printf ("num_nodes = %2d\n", num_nodes_in_set);
	printf ("name = '%s'\n", nset_names[i]);
	free(nset_names[i]);
        node_list = (int *) calloc(num_nodes_in_set, sizeof(int));
        dist_fact = (float *) calloc(num_nodes_in_set, sizeof(float));

        error = ex_get_node_set (exoid, ids[i], node_list);
        printf ("\nafter ex_get_node_set, error = %3d\n", error);

        if (num_df_in_set > 0)
          {
            error = ex_get_node_set_dist_fact (exoid, ids[i], dist_fact);
            printf ("\nafter ex_get_node_set_dist_fact, error = %3d\n", error);
          }

        printf ("\nnode list for node set %2d\n", ids[i]);

        for (j=0; j<num_nodes_in_set; j++)
          {
            printf ("%3d\n", node_list[j]);
          }

        if (num_df_in_set > 0)
          {
            printf ("dist factors for node set %2d\n", ids[i]);

            for (j=0; j<num_nodes_in_set; j++)
              {
                printf ("%5.2f\n", dist_fact[j]);
              }
          }
        else
          printf ("no dist factors for node set %2d\n", ids[i]);

        free (node_list);
        free (dist_fact);

	{
	  int num_attrs = 0;
	  error = ex_get_attr_param(exoid, EX_NODE_SET, ids[i], &num_attrs);
	  printf ("num nodeset attributes for nodeset %d = %d\n", ids[i], num_attrs);
	  if (num_attrs > 0) {
	    for (j=0; j<num_attrs; j++) {
	      attrib_names[j] = (char *)calloc ((MAX_STR_LENGTH+1), sizeof(char));
	    }
	    error = ex_get_attr_names (exoid, EX_NODE_SET, ids[i], attrib_names);
	    printf (" after ex_get_attr_names, error = %d\n", error);
	    
	    if (error == 0) {
	      attrib = (float *) calloc(num_nodes_in_set,sizeof(float));
	      for (j=0; j<num_attrs; j++) {
		printf ("nodeset attribute %d = '%s'\n", j, attrib_names[j]);
		error = ex_get_one_attr(exoid, EX_NODE_SET, ids[i], j+1, attrib);
		for (k=0; k < num_nodes_in_set; k++) {
		  printf ("%5.1f\n", attrib[k]);
		}
		free(attrib_names[j]);
	      }
	      free(attrib);
	    }
	  }
	}
      }
    free(ids);

    /* read node set properties */
    error = ex_inquire (exoid, EX_INQ_NS_PROP, &num_props, &fdum, cdum);
    printf ("\nafter ex_inquire, error = %d\n", error);
    printf ("\nThere are %2d properties for each node set\n", num_props);

    for (i=0; i<num_props; i++)
      {
        prop_names[i] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
      }
    prop_values = (int *) calloc (num_node_sets, sizeof(int));

    error = ex_get_prop_names(exoid,EX_NODE_SET,prop_names);
    printf ("after ex_get_prop_names, error = %d\n", error);


    for (i=0; i<num_props; i++)
      {
        error = ex_get_prop_array(exoid, EX_NODE_SET, prop_names[i],
                                  prop_values);
        if (error == 0)
          for (j=0; j<num_node_sets; j++)
            printf ("node set %2d, property(%2d): '%s'= %5d\n",
                    j+1, i+1, prop_names[i], prop_values[j]);
        else
          printf ("after ex_get_prop_array, error = %d\n", error);
      }
    for (i=0; i<num_props; i++)
      free(prop_names[i]);
    free(prop_values);

    /* read concatenated node sets; this produces the same information as
     * the above code which reads individual node sets
     */

    error = ex_inquire (exoid, EX_INQ_NODE_SETS, &num_node_sets, &fdum, cdum);
    printf ("\nafter ex_inquire, error = %3d\n",error);

    ids = (int *) calloc(num_node_sets, sizeof(int));
    num_nodes_per_set = (int *) calloc(num_node_sets, sizeof(int));
    num_df_per_set = (int *) calloc(num_node_sets, sizeof(int));
    node_ind = (int *) calloc(num_node_sets, sizeof(int));
    df_ind = (int *) calloc(num_node_sets, sizeof(int));

    error = ex_inquire (exoid, EX_INQ_NS_NODE_LEN, &list_len, &fdum, cdum);
    printf ("\nafter ex_inquire: EX_INQ_NS_NODE_LEN = %d, error = %3d\n",
            list_len, error);
    node_list = (int *) calloc(list_len, sizeof(int));

    error = ex_inquire (exoid, EX_INQ_NS_DF_LEN, &list_len, &fdum, cdum);
    printf ("\nafter ex_inquire: EX_INQ_NS_DF_LEN = %d, error = %3d\n",
            list_len, error);
    dist_fact = (float *) calloc(list_len, sizeof(float));

    error = ex_get_concat_node_sets (exoid,ids,num_nodes_per_set,num_df_per_set,
                                     node_ind, df_ind, node_list, dist_fact);
    printf ("\nafter ex_get_concat_node_sets, error = %3d\n", error);

    printf ("\nconcatenated node set info\n");

    printf ("ids = \n");
    for (i=0; i<num_node_sets; i++) printf ("%3d\n", ids[i]);

    printf ("num_nodes_per_set = \n");
    for (i=0; i<num_node_sets; i++) printf ("%3d\n", num_nodes_per_set[i]);

    printf ("node_ind = \n");
    for (i=0; i<num_node_sets; i++) printf ("%3d\n", node_ind[i]);

    printf ("node_list = \n");
    for (i=0; i<list_len; i++) printf ("%3d\n", node_list[i]);

    printf ("dist_fact = \n");
    for (i=0; i<list_len; i++) printf ("%5.3f\n", dist_fact[i]);

    free (ids);
    free (df_ind);
    free (node_ind);
    free (num_df_per_set);
    free (node_list);
    free (dist_fact);
  }

  /* read individual side sets */

  if (num_side_sets > 0) {
    ids = (int *) calloc(num_side_sets, sizeof(int));
    
    error = ex_get_side_set_ids (exoid, ids);
    printf ("\nafter ex_get_side_set_ids, error = %3d\n", error);
    
    for (i=0; i<num_side_sets; i++) {
      sset_names[i] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
    }

    error = ex_get_names(exoid, EX_SIDE_SET, sset_names);
    printf ("\nafter ex_get_names, error = %3d\n", error);

    for (i=0; i<num_side_sets; i++)
      {
	ex_get_name(exoid, EX_SIDE_SET, ids[i], name);
	if (strcmp(name, sset_names[i]) != 0) {
	  printf ("error in ex_get_name for sideset id %d\n", ids[i]);
	}

        error = ex_get_side_set_param (exoid, ids[i], &num_sides_in_set, 
                                       &num_df_in_set);
        printf ("\nafter ex_get_side_set_param, error = %3d\n", error);
        
        printf ("side set %2d parameters:\n",ids[i]);
	printf ("name = '%s'\n", sset_names[i]);
        printf ("num_sides = %3d\n",num_sides_in_set);
        printf ("num_dist_factors = %3d\n", num_df_in_set);
	free(sset_names[i]);
        
        
        /* Note: The # of elements is same as # of sides!  */
        num_elem_in_set = num_sides_in_set;
        elem_list = (int *) calloc(num_elem_in_set, sizeof(int));
        side_list = (int *) calloc(num_sides_in_set, sizeof(int));
        node_ctr_list = (int *) calloc(num_elem_in_set, sizeof(int));
        node_list = (int *) calloc(num_elem_in_set*21, sizeof(int));
        dist_fact = (float *) calloc(num_df_in_set, sizeof(float));
        
        error = ex_get_side_set (exoid, ids[i], elem_list, side_list);
        printf ("\nafter ex_get_side_set, error = %3d\n", error);
        
        error = ex_get_side_set_node_list (exoid, ids[i], node_ctr_list,
                                           node_list);
        printf ("\nafter ex_get_side_set_node_list, error = %3d\n", error);
        
        if (num_df_in_set > 0)
          {
            error = ex_get_side_set_dist_fact (exoid, ids[i], dist_fact);
            printf ("\nafter ex_get_side_set_dist_fact, error = %3d\n", error);
          }
        
        printf ("element list for side set %2d\n", ids[i]);
        for (j=0; j<num_elem_in_set; j++)
          {
            printf ("%3d\n", elem_list[j]);
          }
        
        printf ("side list for side set %2d\n", ids[i]);
        for (j=0; j<num_sides_in_set; j++)
          {
            printf ("%3d\n", side_list[j]);
          }
        
        node_ctr = 0;
        printf ("node list for side set %2d\n", ids[i]);
        for (k=0; k<num_elem_in_set; k++)
          {
            for (j=0; j<node_ctr_list[k]; j++)
              {
                printf ("%3d\n", node_list[node_ctr+j]);
              }
            node_ctr += node_ctr_list[k];
          }
        
        if (num_df_in_set > 0)
          {
            printf ("dist factors for side set %2d\n", ids[i]);
            
            for (j=0; j<num_df_in_set; j++)
              {
                printf ("%5.3f\n", dist_fact[j]);
              }
          }
        else
          printf ("no dist factors for side set %2d\n", ids[i]);
        
        free (elem_list);
        free (side_list);
        free (node_ctr_list);
        free (node_list);
        free (dist_fact);
        
      }
    
    /* read side set properties */
    error = ex_inquire (exoid, EX_INQ_SS_PROP, &num_props, &fdum, cdum);
    printf ("\nafter ex_inquire, error = %d\n", error);
    printf ("\nThere are %2d properties for each side set\n", num_props);
    
    for (i=0; i<num_props; i++)
      {
        prop_names[i] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
      }
    
    error = ex_get_prop_names(exoid,EX_SIDE_SET,prop_names);
    printf ("after ex_get_prop_names, error = %d\n", error);


    for (i=0; i<num_props; i++)
      {
        for (j=0; j<num_side_sets; j++)
          {
            error = ex_get_prop(exoid, EX_SIDE_SET, ids[j], prop_names[i],
                                &prop_value);
            if (error == 0)
              printf ("side set %2d, property(%2d): '%s'= %5d\n",
                      j+1, i+1, prop_names[i], prop_value);
            else
              printf ("after ex_get_prop, error = %d\n", error);
          }
      }
    for (i=0; i<num_props; i++)
      free(prop_names[i]);
    free (ids);

    error = ex_inquire (exoid, EX_INQ_SIDE_SETS, &num_side_sets, &fdum, cdum);
    printf ("\nafter ex_inquire: EX_INQ_SIDE_SETS = %d,  error = %d\n",
            num_side_sets, error);

    if (num_side_sets > 0)
      {
        error = ex_inquire(exoid, EX_INQ_SS_ELEM_LEN, &elem_list_len, &fdum, cdum);
        printf ("\nafter ex_inquire: EX_INQ_SS_ELEM_LEN = %d,  error = %d\n",
                elem_list_len, error);

        error = ex_inquire(exoid, EX_INQ_SS_NODE_LEN, &node_list_len, &fdum, cdum);
        printf ("\nafter ex_inquire: EX_INQ_SS_NODE_LEN = %d,  error = %d\n",
                node_list_len, error);

        error = ex_inquire(exoid, EX_INQ_SS_DF_LEN, &df_list_len, &fdum, cdum);
        printf ("\nafter ex_inquire: EX_INQ_SS_DF_LEN = %d,  error = %d\n",
                df_list_len, error);
      }

    /* read concatenated side sets; this produces the same information as
     * the above code which reads individual side sets
     */

    /* concatenated side set read */

    if (num_side_sets > 0) {
      ids = (int *) calloc(num_side_sets, sizeof(int));
      num_elem_per_set = (int *) calloc(num_side_sets, sizeof(int));
      num_df_per_set = (int *) calloc(num_side_sets, sizeof(int));
      elem_ind = (int *) calloc(num_side_sets, sizeof(int));
      df_ind = (int *) calloc(num_side_sets, sizeof(int));
      elem_list = (int *) calloc(elem_list_len, sizeof(int));
      side_list = (int *) calloc(elem_list_len, sizeof(int));
      dist_fact = (float *) calloc(df_list_len, sizeof(float));
     
      error = ex_get_concat_side_sets (exoid, ids, num_elem_per_set, 
                                       num_df_per_set, elem_ind, df_ind, 
                                       elem_list, side_list, dist_fact);
      printf ("\nafter ex_get_concat_side_sets, error = %3d\n", error);
     
      printf ("concatenated side set info\n");
     
      printf ("ids = \n");
      for (i=0; i<num_side_sets; i++) printf ("%3d\n", ids[i]);
     
      printf ("num_elem_per_set = \n");
      for (i=0; i<num_side_sets; i++) printf ("%3d\n", num_elem_per_set[i]);
     
      printf ("num_dist_per_set = \n");
      for (i=0; i<num_side_sets; i++) printf ("%3d\n", num_df_per_set[i]);
     
      printf ("elem_ind = \n");
      for (i=0; i<num_side_sets; i++) printf ("%3d\n", elem_ind[i]);
     
      printf ("dist_ind = \n");
      for (i=0; i<num_side_sets; i++) printf ("%3d\n", df_ind[i]);
     
      printf ("elem_list = \n");
      for (i=0; i<elem_list_len; i++) printf ("%3d\n", elem_list[i]);
     
      printf ("side_list = \n");
      for (i=0; i<elem_list_len; i++) printf ("%3d\n", side_list[i]);
     
      printf ("dist_fact = \n");
      for (i=0; i<df_list_len; i++) printf ("%5.3f\n", dist_fact[i]);
     
      free (ids);
      free (num_df_per_set);
      free (df_ind);
      free (elem_ind);
      free (elem_list);
      free (side_list);
      free (dist_fact);
    }
  }   
  /* end of concatenated side set read */

  /* read QA records */

  ex_inquire (exoid, EX_INQ_QA, &num_qa_rec, &fdum, cdum);

  for (i=0; i<num_qa_rec; i++)
    {
      for (j=0; j<4; j++)
        {
          qa_record[i][j] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
        }
    }

  error = ex_get_qa (exoid, qa_record); 
  printf ("\nafter ex_get_qa, error = %3d\n", error);

  printf ("QA records = \n");
  for (i=0; i<num_qa_rec; i++) 
    {
      for (j=0; j<4; j++)
        {
          printf (" '%s'\n", qa_record[i][j]);
          free(qa_record[i][j]);
        }
    }

  /* read information records */

  error = ex_inquire (exoid, EX_INQ_INFO, &num_info, &fdum, cdum);
  printf ("\nafter ex_inquire, error = %3d\n", error);

  for (i=0; i<num_info; i++)
    {
      info[i] = (char *) calloc ((MAX_LINE_LENGTH+1), sizeof(char));
    }

  error = ex_get_info (exoid, info); 
  printf ("\nafter ex_get_info, error = %3d\n", error);

  printf ("info records = \n");
  for (i=0; i<num_info; i++)
    {
      printf (" '%s'\n", info[i]);
      free(info[i]);
    }

  /* read global variables parameters and names */

  error = ex_get_var_param (exoid, "g", &num_glo_vars);
  printf ("\nafter ex_get_var_param, error = %3d\n", error);

  for (i=0; i<num_glo_vars; i++)
    {
      var_names[i] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
    }

  error = ex_get_var_names (exoid, "g", num_glo_vars, var_names);
  printf ("\nafter ex_get_var_names, error = %3d\n", error);

  printf ("There are %2d global variables; their names are :\n", 
          num_glo_vars);
  for (i=0; i<num_glo_vars; i++)
    {
      printf (" '%s'\n", var_names[i]);
      free(var_names[i]);
    }

  /* read nodal variables parameters and names */
  num_nod_vars = 0;
  if (num_nodes > 0) {
    error = ex_get_var_param (exoid, "n", &num_nod_vars);
    printf ("\nafter ex_get_var_param, error = %3d\n", error);

    for (i=0; i<num_nod_vars; i++)
      {
        var_names[i] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
      }

    error = ex_get_var_names (exoid, "n", num_nod_vars, var_names);
    printf ("\nafter ex_get_var_names, error = %3d\n", error);

    printf ("There are %2d nodal variables; their names are :\n", num_nod_vars);
    for (i=0; i<num_nod_vars; i++)
      {
        printf (" '%s'\n", var_names[i]);
        free(var_names[i]);
      }
  }

  /* read element variables parameters and names */

  num_ele_vars = 0;
  if (num_elem > 0) {
    error = ex_get_var_param (exoid, "e", &num_ele_vars);
    printf ("\nafter ex_get_var_param, error = %3d\n", error);
     
    for (i=0; i<num_ele_vars; i++)
      {
        var_names[i] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
      }
     
    error = ex_get_var_names (exoid, "e", num_ele_vars, var_names);
    printf ("\nafter ex_get_var_names, error = %3d\n", error);
     
    printf ("There are %2d element variables; their names are :\n", 
            num_ele_vars);
    for (i=0; i<num_ele_vars; i++)
      {
        printf (" '%s'\n", var_names[i]);
        free(var_names[i]);
      }

    /* read element variable truth table */

    if (num_ele_vars > 0) {
      truth_tab = (int *) calloc ((num_elem_blk*num_ele_vars), sizeof(int));
     
      error = ex_get_elem_var_tab (exoid, num_elem_blk, num_ele_vars, truth_tab);
      printf ("\nafter ex_get_elem_var_tab, error = %3d\n", error);
     
      printf ("This is the element variable truth table:\n");
     
      k = 0;
      for (i=0; i<num_elem_blk*num_ele_vars; i++)
        {
          printf ("%2d\n", truth_tab[k++]);
        }
      free (truth_tab);
    }
  }

  /* read nodeset variables parameters and names */

  num_nset_vars = 0;
  if (num_node_sets > 0) {
    error = ex_get_var_param (exoid, "m", &num_nset_vars);
    printf ("\nafter ex_get_var_param, error = %3d\n", error);
     
    if (num_nset_vars > 0) {
      for (i=0; i<num_nset_vars; i++)
	{
	  var_names[i] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
	}
     
      error = ex_get_var_names (exoid, "m", num_nset_vars, var_names);
      printf ("\nafter ex_get_var_names, error = %3d\n", error);
     
      printf ("There are %2d nodeset variables; their names are :\n", 
	      num_nset_vars);
      for (i=0; i<num_nset_vars; i++)
	{
	  printf (" '%s'\n", var_names[i]);
	  free(var_names[i]);
	}

      /* read nodeset variable truth table */

      if (num_nset_vars > 0) {
	truth_tab = (int *) calloc ((num_node_sets*num_nset_vars), sizeof(int));
     
	error = ex_get_nset_var_tab (exoid, num_node_sets, num_nset_vars, truth_tab);
	printf ("\nafter ex_get_nset_var_tab, error = %3d\n", error);
     
	printf ("This is the nodeset variable truth table:\n");
     
	k = 0;
	for (i=0; i<num_node_sets*num_nset_vars; i++)
	  {
	    printf ("%2d\n", truth_tab[k++]);
	  }
	free (truth_tab);
      }
    }
  }

  /* read sideset variables parameters and names */

  num_sset_vars = 0;
  if (num_side_sets > 0) {
    error = ex_get_var_param (exoid, "s", &num_sset_vars);
    printf ("\nafter ex_get_var_param, error = %3d\n", error);
     
    if (num_sset_vars > 0) {
      for (i=0; i<num_sset_vars; i++)
	{
	  var_names[i] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
	}
     
      error = ex_get_var_names (exoid, "s", num_sset_vars, var_names);
      printf ("\nafter ex_get_var_names, error = %3d\n", error);
     
      printf ("There are %2d sideset variables; their names are :\n", 
	      num_sset_vars);
      for (i=0; i<num_sset_vars; i++)
	{
	  printf (" '%s'\n", var_names[i]);
	  free(var_names[i]);
	}

      /* read sideset variable truth table */

      if (num_sset_vars > 0) {
	truth_tab = (int *) calloc ((num_side_sets*num_sset_vars), sizeof(int));
     
	error = ex_get_sset_var_tab (exoid, num_side_sets, num_sset_vars, truth_tab);
	printf ("\nafter ex_get_sset_var_tab, error = %3d\n", error);
     
	printf ("This is the sideset variable truth table:\n");
     
	k = 0;
	for (i=0; i<num_side_sets*num_sset_vars; i++)
	  {
	    printf ("%2d\n", truth_tab[k++]);
	  }
	free (truth_tab);
      }
    }
  }

  /* determine how many time steps are stored */

  error = ex_inquire (exoid, EX_INQ_TIME, &num_time_steps, &fdum, cdum);
  printf ("\nafter ex_inquire, error = %3d\n", error);
  printf ("There are %2d time steps in the database.\n", num_time_steps);

  /* read time value at one time step */

  time_step = 3;
  error = ex_get_time (exoid, time_step, &time_value);
  printf ("\nafter ex_get_time, error = %3d\n", error);

  printf ("time value at time step %2d = %5.3f\n", time_step, time_value);

  /* read time values at all time steps */

  time_values = (float *) calloc (num_time_steps, sizeof(float));

  error = ex_get_all_times (exoid, time_values);
  printf ("\nafter ex_get_all_times, error = %3d\n", error);

  printf ("time values at all time steps are:\n");
  for (i=0; i<num_time_steps; i++) printf ("%5.3f\n", time_values[i]);

  free (time_values);

  /* read all global variables at one time step */

  var_values = (float *) calloc (num_glo_vars, sizeof(float));

  error = ex_get_glob_vars (exoid, time_step, num_glo_vars, var_values);
  printf ("\nafter ex_get_glob_vars, error = %3d\n", error);

  printf ("global variable values at time step %2d\n", time_step);
  for (i=0; i<num_glo_vars; i++) printf ("%5.3f\n", var_values[i]);

  free (var_values); 

  /* read a single global variable through time */

  var_index = 1;
  beg_time = 1;
  end_time = -1;

  var_values = (float *) calloc (num_time_steps, sizeof(float));

  error = ex_get_glob_var_time (exoid, var_index, beg_time, end_time, 
                                var_values);
  printf ("\nafter ex_get_glob_var_time, error = %3d\n", error);

  printf ("global variable %2d values through time:\n", var_index);
  for (i=0; i<num_time_steps; i++) printf ("%5.3f\n", var_values[i]);

  free (var_values); 

  /* read a nodal variable at one time step */

  if (num_nodes > 0) {
    var_values = (float *) calloc (num_nodes, sizeof(float));

    error = ex_get_nodal_var (exoid, time_step, var_index, num_nodes, 
                              var_values);
    printf ("\nafter ex_get_nodal_var, error = %3d\n", error);

    printf ("nodal variable %2d values at time step %2d\n", var_index, 
            time_step);
    for (i=0; i<num_nodes; i++) printf ("%5.3f\n", var_values[i]);

    free (var_values); 

    /* read a nodal variable through time */

    var_values = (float *) calloc (num_time_steps, sizeof(float));

    node_num = 1;
    error = ex_get_nodal_var_time (exoid, var_index, node_num, beg_time, 
                                   end_time, var_values);
    printf ("\nafter ex_get_nodal_var_time, error = %3d\n", error);

    printf ("nodal variable %2d values for node %2d through time:\n", var_index,
            node_num);
    for (i=0; i<num_time_steps; i++) printf ("%5.3f\n", var_values[i]);

    free (var_values); 
  }
  /* read an element variable at one time step */

  if (num_elem_blk > 0) {
    ids = (int *) calloc(num_elem_blk, sizeof(int));
     
    error = ex_get_elem_blk_ids (exoid, ids);
    printf ("\n after ex_get_elem_blk_ids, error = %3d\n", error);
     
    for (i=0; i<num_elem_blk; i++)
      {
	if (num_elem_in_block[i] > 0) {
	  var_values = (float *) calloc (num_elem_in_block[i], sizeof(float));
         
	  error = ex_get_elem_var (exoid, time_step, var_index, ids[i], 
				   num_elem_in_block[i], var_values);
	  printf ("\nafter ex_get_elem_var, error = %3d\n", error);
         
	  if (!error)
	    {
	      printf 
		("element variable %2d values of element block %2d at time step %2d\n",
		 var_index, ids[i], time_step);
	      for (j=0; j<num_elem_in_block[i]; j++) 
		printf ("%5.3f\n", var_values[j]);
	    }
         
	  free (var_values);
	}
      }
    free (num_elem_in_block);
    free(ids);
  }
  /* read an element variable through time */

  if (num_ele_vars > 0) {
    var_values = (float *) calloc (num_time_steps, sizeof(float));
     
    var_index = 2;
    elem_num = 2;
    error = ex_get_elem_var_time (exoid, var_index, elem_num, beg_time, 
                                  end_time, var_values);
    printf ("\nafter ex_get_elem_var_time, error = %3d\n", error);
     
    printf ("element variable %2d values for element %2d through time:\n", 
            var_index, elem_num);
    for (i=0; i<num_time_steps; i++) printf ("%5.3f\n", var_values[i]);
     
    free (var_values);
  }
   
  /* read a sideset variable at one time step */

  if (num_sset_vars > 0) {
    ids = (int *) calloc(num_side_sets, sizeof(int));
     
    error = ex_get_side_set_ids (exoid, ids);
    printf ("\n after ex_get_side_set_ids, error = %3d\n", error);
     
    for (i=0; i<num_side_sets; i++)
      {
        var_values = (float *) calloc (num_elem_per_set[i], sizeof(float));
         
        error = ex_get_sset_var (exoid, time_step, var_index, ids[i], 
                                 num_elem_per_set[i], var_values);
        printf ("\nafter ex_get_sset_var, error = %3d\n", error);
         
        if (!error)
          {
            printf 
              ("sideset variable %2d values of sideset %2d at time step %2d\n",
               var_index, ids[i], time_step);
            for (j=0; j<num_elem_per_set[i]; j++) 
              printf ("%5.3f\n", var_values[j]);
          }
         
        free (var_values); 
      }
    free (num_elem_per_set);
    free(ids);
  }

  /* read a nodeset variable at one time step */

  if (num_nset_vars > 0) {
    ids = (int *) calloc(num_node_sets, sizeof(int));
     
    error = ex_get_node_set_ids (exoid, ids);
    printf ("\n after ex_get_node_set_ids, error = %3d\n", error);
     
    for (i=0; i<num_node_sets; i++)
      {
        var_values = (float *) calloc (num_nodes_per_set[i], sizeof(float));
         
        error = ex_get_nset_var (exoid, time_step, var_index, ids[i], 
                                 num_nodes_per_set[i], var_values);
        printf ("\nafter ex_get_nset_var, error = %3d\n", error);
         
        if (!error)
          {
            printf 
              ("nodeset variable %2d values of nodeset %2d at time step %2d\n",
               var_index, ids[i], time_step);
            for (j=0; j<num_nodes_per_set[i]; j++) 
              printf ("%5.3f\n", var_values[j]);
          }
         
        free (var_values); 
      }
    free(ids);
  }
  if (num_node_sets > 0)
    free (num_nodes_per_set);

  error = ex_close (exoid);
  printf ("\nafter ex_close, error = %3d\n", error);
  return 0;
}
