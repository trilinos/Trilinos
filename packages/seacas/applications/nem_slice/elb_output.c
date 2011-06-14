/*
 * Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <exodusII.h>
#include <ne_nemesisI.h>

#include "elb_const.h"
#include "elb_allo_const.h"
#include "elb_output_const.h"
#include "elb_err_const.h"
#include "elb_loadbal_const.h"
#include "elb_elem_const.h"
#include "elb_util_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function outputs a load balance file using the ExodusII and NemesisI
 * API.
 *****************************************************************************/
int write_nemesis(char *nemI_out_file,
                  MACHINE_PTR machine,
                  PROB_INFO_PTR problem,
                  MESH_INFO_PTR mesh,
                  LB_INFO_PTR lb,
                  SPHERE_INFO_PTR sphere)
{
  int     cnt, cnt2, cnt3, proc, exoid, cpu_ws, io_ws;
  char    title[MAX_LINE_LENGTH+1], method1[MAX_LINE_LENGTH+1];
  char    method2[MAX_LINE_LENGTH+1], *info[3];

  /* For QA record */
  time_t  time_val;
  char   *ct_ptr, tm_date[30];
  char    qa_date[15], qa_time[10], qa_name[MAX_STR_LENGTH];
  char    qa_vers[10];
  char  **lqa_record;

  int     nsize, i2;
  int    *n_cmap_nodes, *n_cmap_procs;

  int    *num_nmap_cnts, *num_emap_cnts;
  int    *node_proc_ptr, *node_cmap_ids_cc, *node_cmap_cnts_cc;
  int    *elem_proc_ptr, *elem_cmap_ids_cc, *elem_cmap_cnts_cc;

  /* Function prototypes */
  int  cmp_ints(int *, int *);
  int  (*func_ptr)();
/*-----------------------------Execution Begins------------------------------*/

  func_ptr = &cmp_ints;
  cpu_ws = sizeof(float);
  io_ws  = sizeof(float);

  printf("Outputting load balance to file %s\n", nemI_out_file);

  /* Create the load balance file */
  /* NOTE: Typically, opening with EX_SHARE is a bad thing for speed.
     However, on lustre and panasas filesystems, it turns out that the
     EX_SHARE works very good (orders of magnitude faster) for the
     nemesis files.  Until we figure out how best to determine the
     underlying filesystem type, we just always open it EX_SHARE
  */
  if((exoid=ex_create(nemI_out_file, EX_CLOBBER|EX_SHARE, &cpu_ws, &io_ws)) < 0)
  {
    Gen_Error(0, "fatal: failed to create Nemesis file");
    return 0;
  }

  /* Set the error reporting value */
  if(error_lev > 1)
    ex_opts(EX_VERBOSE | EX_DEBUG);

  /* Create the title */
  if(problem->type == NODAL)
    strcpy(method1, "nodal");
  else
    strcpy(method1, "elemental");

  sprintf(title, "nem_slice %s load balance file", method1);

  strcpy(method1, "method1: ");
  strcpy(method2, "method2: ");

  switch(lb->type)
  {
  case MULTIKL:
    strcat(method1, "Multilevel-KL decomposition");
    strcat(method2, "With Kernighan-Lin refinement");
    break;
  case SPECTRAL:
    strcat(method1, "Spectral decomposition");
    break;
  case INERTIAL:
    strcat(method1, "Inertial decomposition");
    break;
  case ZPINCH:
    strcat(method1, "ZPINCH decomposition");
    break;
  case BRICK:
    strcat(method1, "BRICK decomposition");
    break;
  case ZOLTAN_RCB:
    strcat(method1, "RCB decomposition");
    break;
  case ZOLTAN_RIB:
    strcat(method1, "RIB decomposition");
    break;
  case ZOLTAN_HSFC:
    strcat(method1, "HSFC decomposition");
    break;
  case LINEAR:
    strcat(method1, "Linear decomposition");
    break;
  case RANDOM:
    strcat(method1, "Random decomposition");
    break;
  case SCATTERED:
    strcat(method1, "Scattered decomposition");
    break;
  }

  if(lb->refine == KL_REFINE && lb->type != MULTIKL)
    strcat(method2, "with Kernighan-Lin refinement");
  else if(lb->type != MULTIKL)
    strcat(method2, "no refinement");

  switch(lb->num_sects)
  {
  case 1:
    strcat(method1, " via bisection");
    break;
  case 2:
    strcat(method1, " via quadrasection");
    break;
  case 3:
    strcat(method1, " via octasection");
    break;
  }

  info[0] = title;
  info[1] = method1;
  info[2] = method2;

  /* Do some sorting */
  for(proc=0; proc < machine->num_procs; proc++)
  {

    /* Sort node maps */
    qsort(lb->int_nodes[proc], lb->num_int_nodes[proc],
          sizeof(int), func_ptr);
    if(problem->type == NODAL)
    {
      sort2_int_int(lb->num_ext_nodes[proc], (lb->ext_nodes[proc]) - 1,
                    (lb->ext_procs[proc]) - 1);
    }

    /* Sort element maps */
    qsort(lb->int_elems[proc], lb->num_int_elems[proc],
          sizeof(int), func_ptr);
  }

  /* Output the info records */
  if(ex_put_info(exoid, 3, info) < 0)
    Gen_Error(0, "warning: output of info records failed");

  /* Generate a QA record for the utility */
  time_val = time(NULL);
  ct_ptr   = asctime(localtime(&time_val));
  strcpy(tm_date, ct_ptr);

  /* Break string with null characters */
  tm_date[3]  = '\0';
  tm_date[7]  = '\0';
  tm_date[10] = '\0';
  tm_date[19] = '\0';

  sprintf(qa_date, "%s %s %s", &tm_date[8], &tm_date[4], &tm_date[20]);
  sprintf(qa_time, "%s", &tm_date[11]);
  strcpy(qa_name, UTIL_NAME);
  strcpy(qa_vers, ELB_VERSION);

  if(qa_date[strlen(qa_date)-1] == '\n')
    qa_date[strlen(qa_date)-1] = '\0';

  lqa_record = (char **)array_alloc(1, 4, sizeof(char *));
  for(i2=0; i2 < 4; i2++)
    lqa_record[i2] = (char *)array_alloc(1, MAX_STR_LENGTH+1, sizeof(char));

  strcpy(lqa_record[0], qa_name);
  strcpy(lqa_record[1], qa_vers);
  strcpy(lqa_record[2], qa_date);
  strcpy(lqa_record[3], qa_time);

  printf("QA Record:\n");
  for(i2=0; i2 < 4; i2++) {
    printf("\t%s\n", lqa_record[i2]);
  }

  if(ex_put_qa(exoid, 1, (char *(*)[]) &lqa_record[0]) < 0) {
    Gen_Error(0, "fatal: unable to output QA records");
    ex_close(exoid);
    return 0;
  }

  /* free up memory */
  for(i2=0; i2 < 4; i2++)
    free(lqa_record[i2]);

  free(lqa_record);

  /* Output the the initial Nemesis global information */
  if(ne_put_init_global(exoid, mesh->num_nodes, mesh->num_elems,
                        mesh->num_el_blks, 0, 0) < 0)
  {
    Gen_Error(0, "fatal: failed to output initial Nemesis parameters");
    ex_close(exoid);
    return 0;
  }

  /* Set up dummy arrays for ouput */
  if(problem->type == NODAL)
  {
    lb->num_bor_elems = malloc(machine->num_procs * sizeof(int));
    if(!(lb->num_bor_elems))
    {
      Gen_Error(0, "fatal: insufficient memory");
      ex_close(exoid);
      return 0;
    }
    for(cnt=0; cnt < machine->num_procs; lb->num_bor_elems[cnt++] = 0);

    num_nmap_cnts = malloc(2 * machine->num_procs * sizeof(int));
    if(!num_nmap_cnts)
    {
      Gen_Error(0, "fatal: insufficient memory");
      ex_close(exoid);
      return 0;
    }
    num_emap_cnts = num_nmap_cnts + machine->num_procs;

    /* need to check and make sure that there really are comm maps */
    for(cnt=0; cnt < machine->num_procs; cnt++) {
      if (lb->num_bor_nodes[cnt] > 0) num_nmap_cnts[cnt] = 1;
      else                            num_nmap_cnts[cnt] = 0;
    }
    for(cnt=0; cnt < machine->num_procs; num_emap_cnts[cnt++] = 0);

  }
  else	/* Elemental load balance */
  {
    lb->num_ext_nodes = malloc(machine->num_procs * sizeof(int));
    if(!(lb->num_ext_nodes))
    {
      Gen_Error(0, "fatal: insufficient memory");
      ex_close(exoid);
      return 0;
    }
    for(cnt=0; cnt < machine->num_procs; lb->num_ext_nodes[cnt++] = 0);

    num_nmap_cnts = malloc(2 * machine->num_procs * sizeof(int));
    if(!num_nmap_cnts)
    {
      Gen_Error(0, "fatal: insufficient memory");
      ex_close(exoid);
      return 0;
    }
    num_emap_cnts = num_nmap_cnts + machine->num_procs;

    if(((problem->num_vertices)-(sphere->num)) > 0)
    {
      /* need to check and make sure that there really are comm maps */
      for(cnt=0; cnt < machine->num_procs; cnt++) {
        if (lb->num_bor_nodes[cnt] > 0) num_nmap_cnts[cnt] = 1;
        else                            num_nmap_cnts[cnt] = 0;
      }
      for(cnt=0; cnt < machine->num_procs; cnt++) {
        if (lb->num_bor_elems[cnt] > 0) num_emap_cnts[cnt] = 1;
        else                            num_emap_cnts[cnt] = 0;
      }
    }
    else
    {
      for(cnt=0; cnt < machine->num_procs; num_nmap_cnts[cnt++] = 0);
      for(cnt=0; cnt < machine->num_procs; num_emap_cnts[cnt++] = 0);
    }
  }

  if(ne_put_init_info(exoid, machine->num_procs, machine->num_procs, "s") < 0)
  {
    Gen_Error(0, "fatal: unable to output init info");
    ex_close(exoid);
    return 0;
  }

  if(ne_put_loadbal_param_cc(exoid, lb->num_int_nodes, lb->num_bor_nodes,
                             lb->num_ext_nodes, lb->num_int_elems,
                             lb->num_bor_elems, num_nmap_cnts,
                             num_emap_cnts) < 0)
  {
    Gen_Error(0, "fatal: unable to output load-balance parameters");
    ex_close(exoid);
    return 0;
  }
  free(num_nmap_cnts);

  if(problem->type == NODAL)		/* Nodal load balance output */
  {
    /* Free unused memory */
    free(lb->num_bor_elems);

    /* Set up for the concatenated communication map parameters */
    node_proc_ptr = malloc(((3*machine->num_procs)+1)*sizeof(int));
    if(!node_proc_ptr)
    {
      Gen_Error(0, "fatal: insufficient memory");
      ex_close(exoid);
      return 0;
    }
    node_cmap_ids_cc  = node_proc_ptr + machine->num_procs + 1;
    node_cmap_cnts_cc = node_cmap_ids_cc + machine->num_procs;

    node_proc_ptr[0] = 0;
    for(proc=0; proc < machine->num_procs; proc++)
    {
      node_proc_ptr[proc+1]   = node_proc_ptr[proc] + 1;
      node_cmap_cnts_cc[proc] = lb->num_ext_nodes[proc];
      node_cmap_ids_cc[proc]  = 1;
    }

    /* Output the communication map parameters */
    if(ne_put_cmap_params_cc(exoid, node_cmap_ids_cc, node_cmap_cnts_cc,
                             node_proc_ptr, NULL, NULL, NULL) < 0)
    {
      Gen_Error(0, "fatal: unable to output communication map parameters");
      ex_close(exoid);
      return 0;
    }

    /* Output the node and element maps */
    for(proc=0; proc < machine->num_procs; proc++)
    {
      /* Output the nodal map */
      if(ne_put_node_map(exoid, lb->int_nodes[proc], lb->bor_nodes[proc],
                         lb->ext_nodes[proc], proc) < 0)
      {
        Gen_Error(0, "fatal: failed to output node map");
        ex_close(exoid);
        return 0;
      }

      /* Output the elemental map */
      if(ne_put_elem_map(exoid, lb->int_elems[proc], NULL, proc) < 0)
      {
        Gen_Error(0, "fatal: failed to output element map");
        ex_close(exoid);
        return 0;
      }

      /*
       * Reorder the nodal communication maps so that they are ordered
       * by processor and then by global ID.
       */

      /* This is a 2-key sort */
      qsort2(lb->ext_procs[proc], lb->ext_nodes[proc], lb->num_ext_nodes[proc]);

      /* Output the nodal communication map */
      if(ne_put_node_cmap(exoid, 1, lb->ext_nodes[proc],
                          lb->ext_procs[proc], proc) < 0)
      {
        Gen_Error(0, "fatal: failed to output nodal communication map");
        ex_close(exoid);
        return 0;
      }

    } /* End "for(proc=0; proc < machine->num_procs; proc++)" */

    free(node_proc_ptr);
  }
  else if(problem->type == ELEMENTAL)	/* Elemental load balance output */
  {
    /* Free unused memory */
    free(lb->num_ext_nodes);

    node_proc_ptr = malloc(((3*machine->num_procs)+1) * sizeof(int));
    if(!node_proc_ptr)
    {
      Gen_Error(0, "fatal: insufficient memory");
      ex_close(exoid);
      return 0;
    }
    node_cmap_ids_cc  = node_proc_ptr + machine->num_procs + 1;
    node_cmap_cnts_cc = node_cmap_ids_cc + machine->num_procs;

    node_proc_ptr[0] = 0;
    for(proc=0; proc < machine->num_procs; proc++)
    {
      node_proc_ptr[proc+1]   = node_proc_ptr[proc] + 1;

      node_cmap_cnts_cc[proc] = 0;
      for(cnt=0; cnt < lb->num_bor_nodes[proc]; cnt++)
        node_cmap_cnts_cc[proc] += lb->born_proc_cnts[proc][cnt];

      node_cmap_ids_cc[proc]  = 1;
    }

    elem_proc_ptr = malloc(((3*machine->num_procs)+1) * sizeof(int));
    if(!elem_proc_ptr)
    {
      Gen_Error(0, "fata: insufficient memory");
      ex_close(exoid);
      return 0;
    }
    elem_cmap_ids_cc  = elem_proc_ptr + machine->num_procs + 1;
    elem_cmap_cnts_cc = elem_cmap_ids_cc + machine->num_procs;

    elem_proc_ptr[0] = 0;
    for(proc=0; proc < machine->num_procs; proc++)
    {
      elem_proc_ptr[proc+1]   = elem_proc_ptr[proc] + 1;
      elem_cmap_cnts_cc[proc] = lb->e_cmap_size[proc];
      elem_cmap_ids_cc[proc]  = 1;
    }

    /* Output the communication map parameters */
    if(ne_put_cmap_params_cc(exoid, node_cmap_ids_cc, node_cmap_cnts_cc,
                             node_proc_ptr, elem_cmap_ids_cc,
                             elem_cmap_cnts_cc, elem_proc_ptr) < 0)
    {
      Gen_Error(0, "fatal: unable to output communication map parameters");
      ex_close(exoid);
      return 0;
    }
    free(elem_proc_ptr);
    free(node_proc_ptr);

    /* Output the node and element maps */
    for(proc=0; proc < machine->num_procs; proc++)
    {
      /* Output the nodal map */
      if(ne_put_node_map(exoid, lb->int_nodes[proc], lb->bor_nodes[proc],
                         NULL, proc) < 0)
      {
        Gen_Error(0, "fatal: failed to output node map");
        ex_close(exoid);
        return 0;
      }

      /* Output the elemental map */
      if(ne_put_elem_map(exoid, lb->int_elems[proc], lb->bor_elems[proc],
                         proc) < 0)
      {
        Gen_Error(0, "fatal: failed to output element map");
        ex_close(exoid);
        return 0;
      }

      /*
       * Build a nodal communication map from the list of border nodes
       * and their associated processors and side IDs.
       */
      nsize = 0;
      for(cnt=0; cnt < lb->num_bor_nodes[proc]; cnt++)
        nsize += lb->born_proc_cnts[proc][cnt];

      if (nsize > 0) {
        n_cmap_nodes = malloc(2*nsize*sizeof(int));
        if(!n_cmap_nodes)
        {
          Gen_Error(0, "fatal: insufficient memory");
          ex_close(exoid);
          return 0;
        }
        n_cmap_procs = n_cmap_nodes + nsize;

        cnt3 = 0;
        for(cnt=0; cnt < lb->num_bor_nodes[proc]; cnt++)
        {
          for(cnt2=0; cnt2 < lb->born_proc_cnts[proc][cnt]; cnt2++)
          {
            n_cmap_nodes[cnt3]   = lb->bor_nodes[proc][cnt];
            n_cmap_procs[cnt3++] = lb->born_procs[proc][cnt][cnt2];
          }
        }

        /*
         * Reorder the nodal communication maps so that they are ordered
         * by processor and then by global ID.
         */
	/* This is a 2-key sort */
	qsort2(n_cmap_procs, n_cmap_nodes, cnt3);

        /* Output the nodal communication map */
        if(ne_put_node_cmap(exoid, 1, n_cmap_nodes, n_cmap_procs, proc) < 0)
        {
          Gen_Error(0, "fatal: unable to output nodal communication map");
          ex_close(exoid);
          return 0;
        }

        free(n_cmap_nodes);

      } /* End "if (nsize > 0)" */

      /* Output the elemental communication map */
      if(lb->e_cmap_size[proc] > 0)
      {
        if(ne_put_elem_cmap(exoid, 1, lb->e_cmap_elems[proc],
                            lb->e_cmap_sides[proc],
                            lb->e_cmap_procs[proc], proc) < 0)
        {
          Gen_Error(0, "fatal: unable to output elemental communication map");
          ex_close(exoid);
          return 0;
        }
      }

    } /* End "for(proc=0; proc < machine->num_procs; proc++)" */

  }

  /* Close the Nemesis file */
  ex_close(exoid);

  return 1;

} /*------------------------End write_nemesis()------------------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function outputs an ExodusII file for the purposes of visualizing
 * the load balance.
 *****************************************************************************/
int write_vis(char *nemI_out_file,
              char *exoII_inp_file,
              MACHINE_PTR machine,
              PROB_INFO_PTR prob,
              MESH_INFO_PTR mesh,
              LB_INFO_PTR lb)
{
  int    exid_vis, exid_inp, cpu_ws=0, io_ws=0, acc_vis;
  int    icpu_ws=0, iio_ws=0;

  char  *cptr, vis_file_name[MAX_FNL+1];
  char   title[MAX_LINE_LENGTH+1];
  char   *coord_names[] = {"X", "Y", "Z"};

  float *xptr, *yptr, *zptr;

  char **elem_name=NULL;
  char  *var_names[] = {"proc"};
  int    ncnt, ecnt, pcnt, nsize, nsize_old, max_np_elem, nnodes, proc;
  int    vis_nelem_blks, pos, old_pos, bcnt, ccnt;
  int   *el_blk_ids=NULL, *el_cnt_blk, *node_pel_blk, *el_ptr;
  int   *nattr_el_blk, *elem_block, *tmp_connect, *vis_el_blk_ptr;
  int   *elem_map;
  float  vers, time_val, *node_vals;
/*-----------------------------Execution Begins------------------------------*/

  /* Generate the file name for the visualization file */
  strcpy(vis_file_name, nemI_out_file);
  cptr = strrchr(vis_file_name, '.');
  strcpy(cptr, "-vis");
  strcat(vis_file_name, ".exoII");

  /* Generate the title for the file */
  strcpy(title, UTIL_NAME);
  strcat(title, " ");
  strcat(title, ELB_VERSION);
  strcat(title, " load balance visualization file");

  /*
   * If the vis technique is to be by element block then calculate the
   * number of element blocks.
   */
  if(prob->type == ELEMENTAL)
    vis_nelem_blks = machine->num_procs;
  else
    vis_nelem_blks = machine->num_procs + 1;

  /* Create the ExodusII file */
  fprintf(stdout, "Outputting load balance visualization file %s\n",
          vis_file_name);
  if((exid_vis=ex_create(vis_file_name, EX_CLOBBER, &cpu_ws, &io_ws)) < 0)
  {
    Gen_Error(0, "fatal: unable to create visualization output file");
    return 0;
  }


  /*
   * Open the original input ExodusII file, read the values for the
   * element blocks and output them to the visualization file.
   */
  if((exid_inp=ex_open(exoII_inp_file, EX_READ, &icpu_ws, &iio_ws, &vers)) < 0)
  {
    Gen_Error(0, "fatal: unable to open input ExodusII file");
    ex_close(exid_vis);
    return 0;
  }

  el_blk_ids = malloc((4*(mesh->num_el_blks)+2*mesh->num_elems+
                       vis_nelem_blks+1)
                      *sizeof(int));
  elem_name  = array_alloc(2, mesh->num_el_blks, MAX_STR_LENGTH+1,
                           sizeof(char));
  if(!el_blk_ids || !elem_name)
  {
    Gen_Error(0, "fatal: insufficient memory");
    ex_close(exid_vis);
    return 0;
  }
  el_cnt_blk     = el_blk_ids + mesh->num_el_blks;
  node_pel_blk   = el_cnt_blk + mesh->num_el_blks;
  nattr_el_blk   = node_pel_blk + mesh->num_el_blks;
  elem_block     = nattr_el_blk + mesh->num_el_blks;
  vis_el_blk_ptr = elem_block + mesh->num_elems;
  elem_map       = vis_el_blk_ptr + (vis_nelem_blks+1);

  if(ex_get_elem_blk_ids(exid_inp, el_blk_ids) < 0)
  {
    Gen_Error(0, "fatal: unable to get element block IDs");
    ex_close(exid_vis);
    ex_close(exid_inp);
    return 0;
  }

  acc_vis = ELB_TRUE;
  max_np_elem = 0;
  nsize = 0;

  /*
   * Find out if the mesh consists of mixed elements. If not then
   * element blocks will be used to visualize the partitioning. Otherwise
   * nodal results will be used.
   */
  for(ecnt=0; ecnt < mesh->num_el_blks; ecnt++)
  {
    if(ex_get_elem_block(exid_inp, el_blk_ids[ecnt], elem_name[ecnt],
                         &el_cnt_blk[ecnt], &node_pel_blk[ecnt],
                         &nattr_el_blk[ecnt]) < 0)
    {
      Gen_Error(0, "fatal: unable to get element block parameters");
      ex_close(exid_vis);
      ex_close(exid_inp);
      return 0;
    }

    if(node_pel_blk[ecnt] > max_np_elem)
      nsize += el_cnt_blk[ecnt]*node_pel_blk[ecnt];

    if(strcmp(elem_name[0], elem_name[ecnt]) == 0 && acc_vis != ELB_FALSE)
    {
      if(node_pel_blk[0] == node_pel_blk[ecnt])
        acc_vis = ELB_TRUE;
      else
        acc_vis = ELB_FALSE;
    }
    else
      acc_vis = ELB_FALSE;
  }

  /*
   * For clearer, more accurate, element block visualization of the
   * partitioning.
   */
  if(acc_vis == ELB_TRUE)
  {

    /* Output the initial information */
    if(ex_put_init(exid_vis, title, mesh->num_dims, mesh->num_nodes,
                   mesh->num_elems, vis_nelem_blks, 0, 0) < 0)
    {
      Gen_Error(0, "fatal: unable to output initial params to vis file");
      ex_close(exid_vis);
      return 0;
    }

    /* Output the nodal coordinates */
    xptr = yptr = zptr = NULL;
    switch(mesh->num_dims)
    {
    case 3:
      zptr = (mesh->coords) + 2*mesh->num_nodes;
      /* FALLTHRU */
    case 2:
      yptr = (mesh->coords) + mesh->num_nodes;
      /* FALLTHRU */
    case 1:
      xptr = mesh->coords;
    }
    if(ex_put_coord(exid_vis, xptr, yptr, zptr) < 0)
    {
      Gen_Error(0, "fatal: unable to output coords to vis file");
      ex_close(exid_vis);
      return 0;
    }
    if(ex_put_coord_names(exid_vis, coord_names) < 0)
    {
      Gen_Error(0, "fatal: unable to output coordinate names");
      ex_close(exid_vis);
      return 0;
    }

    tmp_connect = malloc(nsize*sizeof(int));
    if(!tmp_connect)
    {
      Gen_Error(0, "fatal: insufficient memory");
      ex_close(exid_inp);
      ex_close(exid_vis);
      return 0;
    }
    for(ecnt=0; ecnt < mesh->num_elems; ecnt++)
    {
      elem_map[ecnt] = ecnt+1;
      if(prob->type == ELEMENTAL)
        elem_block[ecnt] = lb->vertex2proc[ecnt];
      else
      {
        proc   = lb->vertex2proc[mesh->connect[ecnt][0]];
        nnodes = get_elem_info(NNODES, mesh->elem_type[ecnt]);
        elem_block[ecnt] = proc;
        for(ncnt=1; ncnt < nnodes; ncnt++)
        {
          if(lb->vertex2proc[mesh->connect[ecnt][ncnt]] != proc)
          {
            elem_block[ecnt] = machine->num_procs;
            break;
          }
        }
      }
    }

    ccnt = 0;
    for(bcnt=0; bcnt < vis_nelem_blks; bcnt++)
    {
      vis_el_blk_ptr[bcnt] = ccnt;
      pos = 0;
      old_pos = 0;
      el_ptr = elem_block;
      ecnt   = mesh->num_elems;
      while(pos != -1)
      {
        pos = in_list(bcnt, ecnt, el_ptr);
        if(pos != -1)
        {
          old_pos += pos + 1;
          ecnt     = mesh->num_elems - old_pos;
          el_ptr   = elem_block + old_pos;
          nnodes = get_elem_info(NNODES, mesh->elem_type[old_pos-1]);
          for(ncnt=0; ncnt < nnodes; ncnt++)
            tmp_connect[ccnt++] = mesh->connect[old_pos-1][ncnt] + 1;
        }
      }
    }
    vis_el_blk_ptr[vis_nelem_blks] = ccnt;

    /* Output the element map */
    if(ex_put_map(exid_vis, elem_map) < 0)
    {
      Gen_Error(0, "fatal: unable to output element number map");
      ex_close(exid_vis);
      return 0;
    }

    /* Output the visualization element blocks */
    for(bcnt=0; bcnt < vis_nelem_blks; bcnt++)
    {
      /*
       * Note this assumes all the blocks contain the same type
       * element.
       */
      ecnt = (vis_el_blk_ptr[bcnt+1]-vis_el_blk_ptr[bcnt])/node_pel_blk[0];
      if(ex_put_elem_block(exid_vis, bcnt+1, elem_name[0],
                           ecnt, node_pel_blk[0], 0) < 0)
      {
        Gen_Error(0, "fatal: unable to output element block params");
        ex_close(exid_vis);
        return 0;
      }

      /* Output the connectivity */
      if(ex_put_elem_conn(exid_vis, bcnt+1,
                          &tmp_connect[vis_el_blk_ptr[bcnt]]) < 0)
      {
        Gen_Error(0, "fatal: unable to output element connectivity");
        ex_close(exid_vis);
        return 0;
      }
    }

    /* Free some unused memory */
    if(tmp_connect)
      free(tmp_connect);
    if(el_blk_ids)
      free(el_blk_ids);
    if(elem_name)
      free(elem_name);

  }
  else	/* For nodal results visualization of the partioning. */
  {
    /* Output the initial information */
    if(ex_put_init(exid_vis, title, mesh->num_dims, mesh->num_nodes,
                   mesh->num_elems, mesh->num_el_blks, 0, 0) < 0)
    {
      Gen_Error(0, "fatal: unable to output initial params to vis file");
      ex_close(exid_vis);
      return 0;
    }

    /* Output the nodal coordinates */
    xptr = yptr = zptr = NULL;
    switch(mesh->num_dims)
    {
    case 3:
      zptr = (mesh->coords) + 2*mesh->num_nodes;
      /* FALLTHRU */
    case 2:
      yptr = (mesh->coords) + mesh->num_nodes;
      /* FALLTHRU */
    case 1:
      xptr = mesh->coords;
    }
    if(ex_put_coord(exid_vis, xptr, yptr, zptr) < 0)
    {
      Gen_Error(0, "fatal: unable to output coords to vis file");
      ex_close(exid_vis);
      return 0;
    }
    if(ex_put_coord_names(exid_vis, coord_names) < 0)
    {
      Gen_Error(0, "fatal: unable to output coordinate names");
      ex_close(exid_vis);
      return 0;
    }

    nsize_old = 0;
    for(ecnt=0; ecnt < mesh->num_el_blks; ecnt++)
    {
      nsize = el_cnt_blk[ecnt] * node_pel_blk[ecnt] * sizeof(int);
      if(nsize > nsize_old)
      {
        if(nsize_old == 0)
          tmp_connect = malloc(nsize);
        else
        {
          tmp_connect = realloc(tmp_connect, nsize);
          nsize_old = nsize;
        }

        if(!tmp_connect)
        {
          Gen_Error(0, "fatal: insufficient memory");
          ex_close(exid_vis);
          ex_close(exid_inp);
          return 0;
        }
      }

      if(ex_get_elem_conn(exid_inp, el_blk_ids[ecnt], tmp_connect) < 0)
      {
        Gen_Error(0, "fatal: unable to get element connectivity");
        ex_close(exid_vis);
        ex_close(exid_inp);
        return 0;
      }

      if(ex_put_elem_block(exid_vis, el_blk_ids[ecnt], elem_name[ecnt],
                           el_cnt_blk[ecnt], node_pel_blk[ecnt],
                           nattr_el_blk[ecnt]) < 0)
      {
        Gen_Error(0, "fatal: unable to output element block parameters");
        ex_close(exid_vis);
        ex_close(exid_inp);
        return 0;
      }

      if(ex_put_elem_conn(exid_vis, el_blk_ids[ecnt], tmp_connect) < 0)
      {
        Gen_Error(0, "fatal: unable to output element connectivity");
        ex_close(exid_vis);
        ex_close(exid_inp);
        return 0;
      }

    }


    /* Free some memory */
    if(tmp_connect)
      free(tmp_connect);

    if(el_blk_ids)
      free(el_blk_ids);

    if(elem_name)
      free(elem_name);

    /* Allocate memory for the nodal values */
    node_vals = malloc(mesh->num_nodes * sizeof(float));
    if(!node_vals)
    {
      Gen_Error(0, "fatal: insufficient memory");
      ex_close(exid_vis);
      return 0;
    }

    /* Set up the file for nodal results */
    time_val = 0.0;
    if(ex_put_time(exid_vis, 1, &time_val) < 0)
    {
      Gen_Error(0, "fatal: unable to output time to vis file");
      ex_close(exid_vis);
      return 0;
    }
    if(ex_put_var_param(exid_vis, "n", 1) < 0)
    {
      Gen_Error(0, "fatal: unable to output var params to vis file");
      ex_close(exid_vis);
      return 0;
    }
    if(ex_put_var_names(exid_vis, "n", 1, var_names) < 0)
    {
      Gen_Error(0, "fatal: unable to output variable name");
      ex_close(exid_vis);
      return 0;
    }

    /* Do some problem specific assignment */
    if(prob->type == NODAL)
    {
      for(ncnt=0; ncnt < mesh->num_nodes; ncnt++)
        node_vals[ncnt] = lb->vertex2proc[ncnt];

      for(pcnt=0; pcnt < machine->num_procs; pcnt++)
      {
        for(ncnt=0; ncnt < lb->num_bor_nodes[pcnt]; ncnt++)
          node_vals[lb->bor_nodes[pcnt][ncnt]] = machine->num_procs + 1;
      }

    }
    else if(prob->type == ELEMENTAL)
    {
      for(pcnt=0; pcnt < machine->num_procs; pcnt++)
      {
        for(ncnt=0; ncnt < lb->num_int_nodes[pcnt]; ncnt++)
          node_vals[lb->int_nodes[pcnt][ncnt]] = pcnt;
        for(ncnt=0; ncnt < lb->num_bor_nodes[pcnt]; ncnt++)
          node_vals[lb->bor_nodes[pcnt][ncnt]] = machine->num_procs;
      }
    }

    /* Output the nodal variables */
    if(ex_put_nodal_var(exid_vis, 1, 1, mesh->num_nodes, node_vals) < 0)
    {
      Gen_Error(0, "fatal: unable to output nodal variables");
      ex_close(exid_vis);
      return 0;
    }

    /* Free unused memory */
    free(node_vals);
  }

  /* Close the visualization file */
  ex_close(exid_vis);

  /* Close the input ExodusII file */
  ex_close(exid_inp);

  return 1;

} /*---------------------------End write_vis()-------------------------------*/

/*****************************************************************************/
/* This function is used by the qsort() to compare two integer values
 *****************************************************************************/
int cmp_ints(int *pval1, int *pval2)
{

  if(*pval1 < *pval2)return -1;
  if(*pval1 == *pval2)return 0;

  return 1;

}
