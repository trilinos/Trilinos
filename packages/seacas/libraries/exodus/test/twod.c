/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "exodusII.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
  int         num_glo_vars  = 10;
  int         num_nod_vars  = 2;
  int         CPU_word_size = 8;
  int         IO_word_size  = 8;
  const char *title         = "This is a 2D mesh example with tri, quad, beam, truss, circle";
  int         ebids[]       = {100, 200, 300, 400, 500};
  int         num_dim       = 2;
  int         num_nodes     = 13;
  int         num_elem      = 20;
  int         num_elem_blk  = 5;
  int         num_node_sets = 2;
  int         num_side_sets = 2;

  /* create EXODUS II file */
  int exoid = ex_create("twod.e",       /* filename path */
                        EX_CLOBBER,     /* create mode */
                        &CPU_word_size, /* CPU float word size in bytes */
                        &IO_word_size); /* I/O float word size in bytes */

  ex_opts(EX_VERBOSE);

  /* initialize file with parameters */
  ex_put_init(exoid, title, num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets,
              num_side_sets);

  /* write nodal coordinates values and names to database */
  {
    double x[13], y[13];
    x[0]  = 0.0;
    y[0]  = 0.0;
    x[1]  = -0.5;
    y[1]  = -0.5;
    x[2]  = 0.5;
    y[2]  = -0.5;
    x[3]  = 0.5;
    y[3]  = 0.5;
    x[4]  = -0.5;
    y[4]  = 0.5;
    x[5]  = -1.0;
    y[5]  = -1.0;
    x[6]  = 1.0;
    y[6]  = -1.0;
    x[7]  = 1.0;
    y[7]  = 1.0;
    x[8]  = -1.0;
    y[8]  = 1.0;
    x[9]  = -2.0;
    y[9]  = 0.0;
    x[10] = 0.0;
    y[10] = -2.0;
    x[11] = 2.0;
    y[11] = 0.0;
    x[12] = 0.0;
    y[12] = 2.0;

    ex_put_coord(exoid, x, y, 0);
  }

  {
    const char *coord_names[] = {"xcoor", "ycoor"};
    ex_put_coord_names(exoid, (char **)coord_names);
  }

  {
    int node_map[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130};
    ex_put_node_num_map(exoid, node_map);
  }

  /* write element order map */
  {
    int elem_map[] = {11,  21,  31,  41,  52,  62,  72,  82,  93,  103,
                      113, 123, 133, 143, 153, 163, 174, 184, 194, 204};
    ex_put_elem_num_map(exoid, elem_map);
  }

  /* write element block parameters */
  {
    const char *block_names[]        = {"Triangles", "Quadrilaterals", "", "Trusses", "Circles"};
    int         num_elem_in_block[]  = {4, 4, 4, 4, 4};
    int         num_nodes_per_elem[] = {3, 4, 2, 2, 1};

    ex_put_block(exoid, EX_ELEM_BLOCK, ebids[0], "triangle", num_elem_in_block[0],
                 num_nodes_per_elem[0], 0, 0, 0);
    ex_put_block(exoid, EX_ELEM_BLOCK, ebids[1], "quad", num_elem_in_block[1],
                 num_nodes_per_elem[1], 0, 0, 0);
    ex_put_block(exoid, EX_ELEM_BLOCK, ebids[2], "beam", num_elem_in_block[2],
                 num_nodes_per_elem[2], 0, 0, 3);
    ex_put_block(exoid, EX_ELEM_BLOCK, ebids[3], "truss", num_elem_in_block[3],
                 num_nodes_per_elem[3], 0, 0, 1);
    ex_put_block(exoid, EX_ELEM_BLOCK, ebids[4], "circle", num_elem_in_block[4],
                 num_nodes_per_elem[4], 0, 0, 2);

    /* Write element block names */
    ex_put_names(exoid, EX_ELEM_BLOCK, (char **)block_names);
  }

  /* write element connectivity */
  {
    int conn_t[] = {2, 3, 1, 3, 4, 1, 4, 5, 1, 5, 2, 1};
    int conn_q[] = {6, 7, 3, 2, 7, 8, 4, 3, 8, 9, 5, 4, 9, 6, 2, 5};
    int conn_B[] = {11, 7, 8, 13, 13, 9, 6, 11};
    int conn_T[] = {10, 6, 9, 10, 7, 12, 12, 8};
    int conn_c[] = {6, 7, 8, 9};

    ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[0], conn_t, NULL, NULL);
    ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[1], conn_q, NULL, NULL);
    ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[2], conn_B, NULL, NULL);
    ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[3], conn_T, NULL, NULL);
    ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[4], conn_c, NULL, NULL);
  }

  /* write element block attributes */
  {
    const char *attn_T[] = {"Area"};
    double      attr_T[] = {1.0, 1.1, 1.2, 1.3};

    const char *attn_B[] = {"A", "I", "J"};
    double attr_B[] = {1.0, 100.0, 200.0, 1.1, 100.1, 200.1, 1.2, 100.2, 200.2, 1.3, 100.3, 200.3};

    const char *attn_c[] = {"Radius", "A"};
    double      attr_c[] = {1.0, 3.14, 1.1, 4.14, 1.2, 5.14, 1.3, 6.14};

    ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[2], attr_B);
    ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[3], attr_T);
    ex_put_attr(exoid, EX_ELEM_BLOCK, ebids[4], attr_c);

    ex_put_attr_names(exoid, EX_ELEM_BLOCK, ebids[2], (char **)attn_B);
    ex_put_attr_names(exoid, EX_ELEM_BLOCK, ebids[3], (char **)attn_T);
    ex_put_attr_names(exoid, EX_ELEM_BLOCK, ebids[4], (char **)attn_c);
  }

  /* write individual node sets */
  {
    int         num_nodes_in_nset[] = {5, 8};
    int         nsids[]             = {20, 22};
    int         nod1[]              = {5, 4, 3, 2, 1};
    int         nod2[]              = {6, 7, 8, 9, 2, 3, 4, 5};
    const char *nset_names[]        = {"Triangle_Nodes", "Quadrilateral_Nodes"};

    ex_put_set_param(exoid, EX_NODE_SET, nsids[0], num_nodes_in_nset[0], 0);
    ex_put_set_param(exoid, EX_NODE_SET, nsids[1], num_nodes_in_nset[1], 0);

    ex_put_set(exoid, EX_NODE_SET, nsids[0], nod1, 0);
    ex_put_set(exoid, EX_NODE_SET, nsids[1], nod2, 0);
    ex_put_names(exoid, EX_NODE_SET, (char **)nset_names);
  }

  {
    /* write individual side sets */
    int num_face_in_sset[] = {4, 4};
    int ssids[]            = {100, 200};
    int ss1el[]            = {1, 2, 3, 4};
    int ss1si[]            = {1, 1, 1, 1};

    int         ss2el[]      = {5, 7, 6, 8};
    int         ss2si[]      = {1, 1, 1, 1};
    const char *sset_names[] = {"A", "B"};

    ex_put_set_param(exoid, EX_SIDE_SET, ssids[0], num_face_in_sset[0], 0);
    ex_put_set_param(exoid, EX_SIDE_SET, ssids[1], num_face_in_sset[1], 0);

    ex_put_set(exoid, EX_SIDE_SET, ssids[0], ss1el, ss1si);
    ex_put_set(exoid, EX_SIDE_SET, ssids[1], ss2el, ss2si);
    ex_put_names(exoid, EX_SIDE_SET, (char **)sset_names);
  }

  /* write results variables parameters and names */
  {
    const char *gvarn[] = {"g_01", "g_02", "g_03", "g_04", "g_05",
                           "g_06", "g_07", "g_08", "g_09", "g_10"};
    ex_put_variable_param(exoid, EX_GLOBAL, num_glo_vars);
    ex_put_variable_names(exoid, EX_GLOBAL, num_glo_vars, (char **)gvarn);
  }

  {
    const char *nvarn[] = {"disp_x", "disp_y"};
    ex_put_variable_param(exoid, EX_NODAL, num_nod_vars);
    ex_put_variable_names(exoid, EX_NODAL, num_nod_vars, (char **)nvarn);
  }

#if 0
   num_ele_vars = 3;
   /*              0        1         2         3   */
   /*              12345678901234567890123456789012 */
   var_names[0] = "this_variable_name_is_short";
   var_names[1] = "this_variable_name_is_just_right";
   var_names[2] = "this_variable_name_is_tooooo_long";

   ex_put_variable_param (exoid, EX_ELEM_BLOCK, num_ele_vars);
   printf ("after ex_put_variable_param, %d\n", error);
   if (error) {
     ex_close (exoid);
     exit(-1);
   }

   ex_put_variable_names (exoid, EX_ELEM_BLOCK, num_ele_vars, var_names);
   printf ("after ex_put_variable_names, %d\n", error);
   if (error) {
     ex_close (exoid);
     exit(-1);
   }

   {
     num_nset_vars = 3;

     var_names[0] = "ns_var0";
     var_names[1] = "ns_var1";
     var_names[2] = "ns_var2";

     ex_put_variable_param (exoid, "m", num_nset_vars);
     printf ("after ex_put_variable_param, %d\n", error);
     if (error) {
       ex_close (exoid);
       exit(-1);
     }

     ex_put_variable_names (exoid, "m", num_nset_vars, var_names);
     printf ("after ex_put_variable_names, %d\n", error);
     if (error) {
       ex_close (exoid);
       exit(-1);
     }
   }

   {
     num_sset_vars = 3;

     var_names[0] = "ss_var0";
     var_names[1] = "ss_var1";
     var_names[2] = "ss_var2";

     ex_put_variable_param (exoid, EX_SIDE_SET, num_sset_vars);
     printf ("after ex_put_variable_param, %d\n", error);
     if (error) {
       ex_close (exoid);
       exit(-1);
     }

     ex_put_variable_names (exoid, EX_SIDE_SET, num_sset_vars, var_names);
     printf ("after ex_put_variable_names, %d\n", error);
     if (error) {
       ex_close (exoid);
       exit(-1);
     }
   }
#endif

  /* for each time step, write the analysis results;
   * the code below fills the arrays glob_var_vals,
   * nodal_var_vals, and elem_var_vals with values for debugging purposes;
   * obviously the analysis code will populate these arrays
   */

  {
    int i, j, k;
    int whole_time_step = 1;
    int num_time_steps  = 10;

    double gvar[10];
    double nvar[20];

    for (i = 0; i < num_time_steps; i++) {
      double time_value = (double)(i) / 100.;

      ex_put_time(exoid, whole_time_step, &time_value);

      for (j = 0; j < num_glo_vars; j++) {
        gvar[j] = (double)(j + 2) * time_value;
      }
      ex_put_var(exoid, whole_time_step, EX_GLOBAL, 1, 1, num_glo_vars, gvar);

      /* write nodal variables */
      for (k = 0; k < num_nod_vars; k++) {
        for (j = 0; j < num_nodes; j++) {
          nvar[j] = (double)k + ((double)(j + 1) * time_value);
        }

        ex_put_nodal_var(exoid, whole_time_step, k + 1, num_nodes, nvar);
      }

#if 0
/* write element variables */

     for (k=1; k<=num_ele_vars; k++)
     {
       for (j=0; j<num_elem_blk; j++)
       {
         for (m=0; m<num_elem_in_block[j]; m++)
         {
           elem_var_vals[m] = (float)(k+1) + (float)(j+2) +
                              ((float)(m+1)*time_value);
           /* printf("elem_var_vals[%d]: %f\n",m,elem_var_vals[m]); */
         }
         ex_put_elem_var (exoid, whole_time_step, k, ebids[j],
                                  num_elem_in_block[j], elem_var_vals);
         printf ("after ex_put_elem_var, %d\n", error);
         if (error) {
           ex_close (exoid);
           exit(-1);
         }
       }
     }

/* write sideset variables */

     for (k=1; k<=num_sset_vars; k++)
     {
       for (j=0; j<num_side_sets; j++)
       {
         for (m=0; m<num_face_in_sset[j]; m++)
         {
           sset_var_vals[m] = (float)(k+2) + (float)(j+3) +
                              ((float)(m+1)*time_value);
           /* printf("sset_var_vals[%d]: %f\n",m,sset_var_vals[m]); */
         }
         ex_put_sset_var (exoid, whole_time_step, k, ssids[j],
                                  num_face_in_sset[j], sset_var_vals);
         printf ("after ex_put_sset_var, %d\n", error);
         if (error) {
           ex_close (exoid);
           exit(-1);
         }
       }
     }

/* write nodeset variables */

     for (k=1; k<=num_nset_vars; k++)
     {
       for (j=0; j<num_node_sets; j++)
       {
         for (m=0; m<num_nodes_in_nset[j]; m++)
         {
           nset_var_vals[m] = (float)(k+3) + (float)(j+4) +
                              ((float)(m+1)*time_value);
           /* printf("nset_var_vals[%d]: %f\n",m,nset_var_vals[m]); */
         }
         ex_put_nset_var (exoid, whole_time_step, k, nsids[j],
                                  num_nodes_in_nset[j], nset_var_vals);
         printf ("after ex_put_nset_var, %d\n", error);
         if (error) {
           ex_close (exoid);
           exit(-1);
         }
       }
     }
#endif

      whole_time_step++;
    }
  }
  ex_close(exoid);
  return 0;
}
