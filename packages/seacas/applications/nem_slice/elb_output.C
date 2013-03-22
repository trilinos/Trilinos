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

#include <exodusII.h>                   // for ex_close, MAX_LINE_LENGTH, etc

#include <stddef.h>                     // for size_t
#include <stdio.h>                      // for NULL, printf, sprintf, etc
#include <stdlib.h>                     // for free, malloc, etc
#include <string.h>                     // for strcat, strcpy, strlen, etc
#include <time.h>                       // for asctime, localtime, time, etc

#include "elb_output.h"
#include "elb_allo.h"             // for array_alloc
#include "elb.h"                  // for LB_Description<INT>, etc
#include "elb_elem.h"             // for get_elem_info, E_Type, etc
#include "elb_err.h"              // for Gen_Error, error_lev
#include "elb_util.h"             // for qsort2, in_list, etc

namespace {
  const std::string remove_extension(const std::string &filename)
  {
    // Strip off the extension
    size_t ind = filename.find_last_of('.', filename.size());
    if (ind != std::string::npos)
      return filename.substr(0,ind);
    else
      return filename;
  }
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function outputs a load balance file using the ExodusII and NemesisI
 * API.
 *****************************************************************************/
template int write_nemesis(std::string &nemI_out_file, Machine_Description* machine,
			   Problem_Description* problem, Mesh_Description<int>* mesh,
			   LB_Description<int>* lb, Sphere_Info* sphere);
template int write_nemesis(std::string &nemI_out_file, Machine_Description* machine,
			   Problem_Description* problem, Mesh_Description<int64_t>* mesh,
			   LB_Description<int64_t>* lb, Sphere_Info* sphere);


template <typename INT>
int write_nemesis(std::string &nemI_out_file,
                  Machine_Description* machine,
                  Problem_Description* problem,
                  Mesh_Description<INT>* mesh,
                  LB_Description<INT>* lb,
                  Sphere_Info* sphere)
{
  int     exoid;
  char    title[MAX_LINE_LENGTH+1], method1[MAX_LINE_LENGTH+1];
  char    method2[MAX_LINE_LENGTH+1];

  int cpu_ws = sizeof(float);
  int io_ws  = sizeof(float);

  printf("Outputting load balance to file %s\n", nemI_out_file.c_str());

  /* Create the load balance file */
  /* Attempt to create a netcdf4-format file; if it fails, then assume
     that the netcdf library does not support that mode and fall back
     to classic netcdf3 format.  If that fails, issue an error and
     return failure.
  */
  int mode3 = EX_CLOBBER;
  int mode4 = mode3|EX_NETCDF4|EX_NOCLASSIC|problem->int64db|problem->int64api;

  ex_opts(EX_DEFAULT); // Eliminate misleading error if the first ex_create fails, but the second succeeds.
  if((exoid=ex_create(nemI_out_file.c_str(), mode4, &cpu_ws, &io_ws)) < 0) {
    /* If int64api or int64db non-zero, then netcdf-4 format is required, so
       fail now...
    */
    if (problem->int64db|problem->int64api) {
      Gen_Error(0, "fatal: failed to create Nemesis netcdf-4 file");
      return 0;
    }
    if((exoid=ex_create(nemI_out_file.c_str(), mode3, &cpu_ws, &io_ws)) < 0) {
      Gen_Error(0, "fatal: failed to create Nemesis file");
      return 0;
    }
  }

  /* Set the error reporting value */
  if (error_lev > 1)
    ex_opts(EX_VERBOSE | EX_DEBUG);
  else
    ex_opts(EX_VERBOSE);

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

  /* Do some sorting */
  for(int proc=0; proc < machine->num_procs; proc++) {

    /* Sort node maps */
    gds_qsort(TOPTR(lb->int_nodes[proc]), lb->int_nodes[proc].size());
    if(problem->type == NODAL) {
      sort2(lb->ext_nodes[proc].size(), TOPTR(lb->ext_nodes[proc]) - 1,
	    TOPTR(lb->ext_procs[proc]) - 1);
    }

    /* Sort element maps */
    gds_qsort(TOPTR(lb->int_elems[proc]), lb->int_elems[proc].size());
  }

  /* Output the info records */
  char *info[3];
  info[0] = title;
  info[1] = method1;
  info[2] = method2;

  if(ex_put_info(exoid, 3, info) < 0)
    Gen_Error(0, "warning: output of info records failed");

  /* Generate a QA record for the utility */
  time_t time_val = time(NULL);
  char *ct_ptr   = asctime(localtime(&time_val));
  char tm_date[30];
  strcpy(tm_date, ct_ptr);

  /* Break string with null characters */
  tm_date[3]  = '\0';
  tm_date[7]  = '\0';
  tm_date[10] = '\0';
  tm_date[19] = '\0';

  char    qa_date[15], qa_time[10], qa_name[MAX_STR_LENGTH];
  char    qa_vers[10];

  sprintf(qa_date, "%s %s %s", &tm_date[8], &tm_date[4], &tm_date[20]);
  sprintf(qa_time, "%s", &tm_date[11]);
  strcpy(qa_name, UTIL_NAME);
  strcpy(qa_vers, ELB_VERSION);

  if(qa_date[strlen(qa_date)-1] == '\n')
    qa_date[strlen(qa_date)-1] = '\0';

  char **lqa_record = (char **)array_alloc(1, 4, sizeof(char *));
  for(int i2=0; i2 < 4; i2++)
    lqa_record[i2] = (char *)array_alloc(1, MAX_STR_LENGTH+1, sizeof(char));

  strcpy(lqa_record[0], qa_name);
  strcpy(lqa_record[1], qa_vers);
  strcpy(lqa_record[2], qa_date);
  strcpy(lqa_record[3], qa_time);

  printf("QA Record:\n");
  for(int i2=0; i2 < 4; i2++) {
    printf("\t%s\n", lqa_record[i2]);
  }

  if(ex_put_qa(exoid, 1, (char *(*)[4]) &lqa_record[0]) < 0) {
    Gen_Error(0, "fatal: unable to output QA records");
    ex_close(exoid);
    return 0;
  }

  /* free up memory */
  for(int i2=0; i2 < 4; i2++)
    free(lqa_record[i2]);

  free(lqa_record);

  /* Output the the initial Nemesis global information */
  if(ex_put_init_global(exoid, mesh->num_nodes, mesh->num_elems,
                        mesh->num_el_blks, 0, 0) < 0) {
    Gen_Error(0, "fatal: failed to output initial Nemesis parameters");
    ex_close(exoid);
    return 0;
  }
  
  /* Set up dummy arrays for ouput */
  std::vector<INT> num_nmap_cnts(machine->num_procs);
  std::vector<INT> num_emap_cnts(machine->num_procs);
  
  if(problem->type == NODAL) {
    /* need to check and make sure that there really are comm maps */
    for(int cnt=0; cnt < machine->num_procs; cnt++) {
      if (!lb->bor_nodes[cnt].empty())
	num_nmap_cnts[cnt] = 1;
    }
  }
  else {	/* Elemental load balance */
    if(((problem->num_vertices)-(sphere->num)) > 0) {
      /* need to check and make sure that there really are comm maps */
      for(int cnt=0; cnt < machine->num_procs; cnt++) {
        if (!lb->bor_nodes[cnt].empty()) num_nmap_cnts[cnt] = 1;
      }
      for(int cnt=0; cnt < machine->num_procs; cnt++) {
        if (!lb->bor_elems[cnt].empty()) num_emap_cnts[cnt] = 1;
      }
    }
  }

  if(ex_put_init_info(exoid, machine->num_procs, machine->num_procs, (char*)"s") < 0) {
    Gen_Error(0, "fatal: unable to output init info");
    ex_close(exoid);
    return 0;
  }

  // Need to create 5 arrays with the sizes of lb->int_nodes[i].size()...
  {
    std::vector<INT> ins(machine->num_procs);
    std::vector<INT> bns(machine->num_procs);
    std::vector<INT> ens(machine->num_procs);
    std::vector<INT> ies(machine->num_procs);
    std::vector<INT> bes(machine->num_procs);

    for (int iproc = 0; iproc < machine->num_procs; iproc++) {
      ins[iproc] = lb->int_nodes[iproc].size();
      bns[iproc] = lb->bor_nodes[iproc].size();
      ens[iproc] = lb->ext_nodes[iproc].size();
      ies[iproc] = lb->int_elems[iproc].size();
      bes[iproc] = lb->bor_elems[iproc].size();
    }

    if(ex_put_loadbal_param_cc(exoid,
			       TOPTR(ins), TOPTR(bns), TOPTR(ens),
			       TOPTR(ies), TOPTR(bes), TOPTR(num_nmap_cnts),
			       TOPTR(num_emap_cnts)) < 0)
      {
	Gen_Error(0, "fatal: unable to output load-balance parameters");
	ex_close(exoid);
	return 0;
      }
  }

  if(problem->type == NODAL)		/* Nodal load balance output */
    {
      /* Set up for the concatenated communication map parameters */
      INT *node_proc_ptr = (INT*)malloc(((3*machine->num_procs)+1)*sizeof(INT));
      if(!node_proc_ptr)
	{
	  Gen_Error(0, "fatal: insufficient memory");
	  ex_close(exoid);
	  return 0;
	}
      INT *node_cmap_ids_cc  = node_proc_ptr + machine->num_procs + 1;
      INT *node_cmap_cnts_cc = node_cmap_ids_cc + machine->num_procs;

      node_proc_ptr[0] = 0;
      for(int proc=0; proc < machine->num_procs; proc++) {
	node_proc_ptr[proc+1]   = node_proc_ptr[proc] + 1;
	node_cmap_cnts_cc[proc] = lb->ext_nodes[proc].size();
	node_cmap_ids_cc[proc]  = 1;
      }

      /* Output the communication map parameters */
      if(ex_put_cmap_params_cc(exoid, node_cmap_ids_cc, node_cmap_cnts_cc,
			       node_proc_ptr, NULL, NULL, NULL) < 0)
	{
	  Gen_Error(0, "fatal: unable to output communication map parameters");
	  ex_close(exoid);
	  return 0;
	}

      /* Output the node and element maps */
      for(int proc=0; proc < machine->num_procs; proc++) {
	/* Output the nodal map */
	if(ex_put_processor_node_maps(exoid,
				      TOPTR(lb->int_nodes[proc]),
				      TOPTR(lb->bor_nodes[proc]),
				      TOPTR(lb->ext_nodes[proc]), proc) < 0)
	  {
	    Gen_Error(0, "fatal: failed to output node map");
	    ex_close(exoid);
	    return 0;
	  }

	/* Output the elemental map */
	if(ex_put_processor_elem_maps(exoid, TOPTR(lb->int_elems[proc]), NULL, proc) < 0)
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
	qsort2(TOPTR(lb->ext_procs[proc]), TOPTR(lb->ext_nodes[proc]), lb->ext_nodes[proc].size());

	/* Output the nodal communication map */
	if(ex_put_node_cmap(exoid, 1,
			    TOPTR(lb->ext_nodes[proc]),
			    TOPTR(lb->ext_procs[proc]), proc) < 0)
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
      INT *node_proc_ptr = (INT*)malloc(((3*machine->num_procs)+1) * sizeof(INT));
      if(!node_proc_ptr)
	{
	  Gen_Error(0, "fatal: insufficient memory");
	  ex_close(exoid);
	  return 0;
	}
      INT *node_cmap_ids_cc  = node_proc_ptr + machine->num_procs + 1;
      INT *node_cmap_cnts_cc = node_cmap_ids_cc + machine->num_procs;

      node_proc_ptr[0] = 0;
      for(int proc=0; proc < machine->num_procs; proc++) {
	node_proc_ptr[proc+1]   = node_proc_ptr[proc] + 1;

	node_cmap_cnts_cc[proc] = 0;
	for(size_t cnt=0; cnt < lb->bor_nodes[proc].size(); cnt++)
	  node_cmap_cnts_cc[proc] += lb->born_procs[proc][cnt].size();

	node_cmap_ids_cc[proc]  = 1;
      }

      INT *elem_proc_ptr = (INT*)malloc(((3*machine->num_procs)+1) * sizeof(INT));
      if(!elem_proc_ptr)
	{
	  Gen_Error(0, "fata: insufficient memory");
	  ex_close(exoid);
	  return 0;
	}
      INT *elem_cmap_ids_cc  = elem_proc_ptr + machine->num_procs + 1;
      INT *elem_cmap_cnts_cc = elem_cmap_ids_cc + machine->num_procs;

      elem_proc_ptr[0] = 0;
      for(int proc=0; proc < machine->num_procs; proc++) {
	elem_proc_ptr[proc+1]   = elem_proc_ptr[proc] + 1;
	elem_cmap_cnts_cc[proc] = lb->e_cmap_elems[proc].size();
	elem_cmap_ids_cc[proc]  = 1;
      }

      /* Output the communication map parameters */
      if(ex_put_cmap_params_cc(exoid, node_cmap_ids_cc, node_cmap_cnts_cc,
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
      for(int proc=0; proc < machine->num_procs; proc++)
	{
	  /* Output the nodal map */
	  if(ex_put_processor_node_maps(exoid,
					TOPTR(lb->int_nodes[proc]),
					TOPTR(lb->bor_nodes[proc]),
					NULL, proc) < 0)
	    {
	      Gen_Error(0, "fatal: failed to output node map");
	      ex_close(exoid);
	      return 0;
	    }

	  /* Output the elemental map */
	  if(ex_put_processor_elem_maps(exoid,
					TOPTR(lb->int_elems[proc]),
					TOPTR(lb->bor_elems[proc]),
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
	  size_t nsize = 0;
	  for(size_t cnt=0; cnt < lb->bor_nodes[proc].size(); cnt++)
	    nsize += lb->born_procs[proc][cnt].size();

	  if (nsize > 0) {
	    INT *n_cmap_nodes = (INT*)malloc(2*nsize*sizeof(INT));
	    if(!n_cmap_nodes)
	      {
		Gen_Error(0, "fatal: insufficient memory");
		ex_close(exoid);
		return 0;
	      }
	    INT *n_cmap_procs = n_cmap_nodes + nsize;

	    size_t cnt3 = 0;
	    for(size_t cnt=0; cnt < lb->bor_nodes[proc].size(); cnt++) {
	      for(size_t cnt2=0; cnt2 < lb->born_procs[proc][cnt].size(); cnt2++) {
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
	    if(ex_put_node_cmap(exoid, 1, n_cmap_nodes, n_cmap_procs, proc) < 0)
	      {
		Gen_Error(0, "fatal: unable to output nodal communication map");
		ex_close(exoid);
		return 0;
	      }

	    free(n_cmap_nodes);

	  } /* End "if (nsize > 0)" */

	    /* Output the elemental communication map */
	  if(!lb->e_cmap_elems[proc].empty()) {
	    if(ex_put_elem_cmap(exoid, 1,
				TOPTR(lb->e_cmap_elems[proc]),
				TOPTR(lb->e_cmap_sides[proc]),
				TOPTR(lb->e_cmap_procs[proc]), proc) < 0)
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
  template int write_vis(std::string &nemI_out_file, std::string &exoII_inp_file,
			 Machine_Description* machine, Problem_Description* prob,
			 Mesh_Description<int>* mesh, LB_Description<int>* lb);
  template int write_vis(std::string &nemI_out_file, std::string &exoII_inp_file,
			 Machine_Description* machine, Problem_Description* prob,
			 Mesh_Description<int64_t>* mesh, LB_Description<int64_t>* lb);


  template <typename INT>
    int write_vis(std::string &nemI_out_file,
		  std::string &exoII_inp_file,
		  Machine_Description* machine,
		  Problem_Description* prob,
		  Mesh_Description<INT>* mesh,
		  LB_Description<INT>* lb)
    {
      int    exid_vis, exid_inp, acc_vis;

      char  title[MAX_LINE_LENGTH+1];
      const char   *coord_names[] = {"X", "Y", "Z"};

      /*-----------------------------Execution Begins------------------------------*/

      /* Generate the file name for the visualization file */
      std::string vis_file_name = remove_extension(nemI_out_file);
      vis_file_name += "-vis.exoII";

      /* Generate the title for the file */
      strcpy(title, UTIL_NAME);
      strcat(title, " ");
      strcat(title, ELB_VERSION);
      strcat(title, " load balance visualization file");

      /*
       * If the vis technique is to be by element block then calculate the
       * number of element blocks.
       */
      int    vis_nelem_blks;
      if(prob->type == ELEMENTAL)
	vis_nelem_blks = machine->num_procs;
      else
	vis_nelem_blks = machine->num_procs + 1;

      /* Create the ExodusII file */
      fprintf(stdout, "Outputting load balance visualization file %s\n",
	      vis_file_name.c_str());
      int cpu_ws = 0;
      int io_ws = 0;
      int mode = EX_CLOBBER;
      if (prob->int64db|prob->int64api) {
	mode |= EX_NETCDF4|EX_NOCLASSIC|prob->int64db|prob->int64api;
      }
      if((exid_vis=ex_create(vis_file_name.c_str(), mode, &cpu_ws, &io_ws)) < 0)
	{
	  Gen_Error(0, "fatal: unable to create visualization output file");
	  return 0;
	}

      /*
       * Open the original input ExodusII file, read the values for the
       * element blocks and output them to the visualization file.
       */
      int icpu_ws=0;
      int iio_ws=0;
      float vers=0.0;
      mode = EX_READ | prob->int64api;
      if((exid_inp=ex_open(exoII_inp_file.c_str(), mode, &icpu_ws, &iio_ws, &vers)) < 0)
	{
	  Gen_Error(0, "fatal: unable to open input ExodusII file");
	  ex_close(exid_vis);
	  return 0;
	}

      INT *el_blk_ids = (INT*)malloc((4*(mesh->num_el_blks)+2*mesh->num_elems+
				      vis_nelem_blks+1)
				     *sizeof(INT));
      char **elem_type  = (char**)array_alloc(2, mesh->num_el_blks, MAX_STR_LENGTH+1,
					      sizeof(char));
      if(!el_blk_ids || !elem_type)
	{
	  Gen_Error(0, "fatal: insufficient memory");
	  ex_close(exid_vis);
	  return 0;
	}
      INT *el_cnt_blk     = el_blk_ids + mesh->num_el_blks;
      INT *node_pel_blk   = el_cnt_blk + mesh->num_el_blks;
      INT *nattr_el_blk   = node_pel_blk + mesh->num_el_blks;
      INT *elem_block     = nattr_el_blk + mesh->num_el_blks;
      INT *vis_el_blk_ptr = elem_block + mesh->num_elems;
      INT *elem_map       = vis_el_blk_ptr + (vis_nelem_blks+1);

      if(ex_get_elem_blk_ids(exid_inp, el_blk_ids) < 0)
	{
	  Gen_Error(0, "fatal: unable to get element block IDs");
	  ex_close(exid_vis);
	  ex_close(exid_inp);
	  return 0;
	}

      acc_vis = ELB_TRUE;
      size_t nsize = 0;

      /*
       * Find out if the mesh consists of mixed elements. If not then
       * element blocks will be used to visualize the partitioning. Otherwise
       * nodal results will be used.
       */
      for(size_t ecnt=0; ecnt < mesh->num_el_blks; ecnt++)
	{
	  if(ex_get_elem_block(exid_inp, el_blk_ids[ecnt], elem_type[ecnt],
			       &el_cnt_blk[ecnt], &node_pel_blk[ecnt],
			       &nattr_el_blk[ecnt]) < 0)
	    {
	      Gen_Error(0, "fatal: unable to get element block parameters");
	      ex_close(exid_vis);
	      ex_close(exid_inp);
	      return 0;
	    }

	  nsize += el_cnt_blk[ecnt]*node_pel_blk[ecnt];

	  if(strcmp(elem_type[0], elem_type[ecnt]) == 0 && acc_vis != ELB_FALSE)
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
	  float *xptr = NULL;
	  float *yptr = NULL;
	  float *zptr = NULL;
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

	  if(ex_put_coord_names(exid_vis, (char**)coord_names) < 0) {
	    Gen_Error(0, "fatal: unable to output coordinate names");
	    ex_close(exid_vis);
	    return 0;
	  }

	  INT* tmp_connect = (INT*)malloc(nsize*sizeof(INT));
	  if(!tmp_connect)
	    {
	      Gen_Error(0, "fatal: insufficient memory");
	      ex_close(exid_inp);
	      ex_close(exid_vis);
	      return 0;
	    }
	  for(size_t ecnt=0; ecnt < mesh->num_elems; ecnt++)
	    {
	      elem_map[ecnt] = ecnt+1;
	      if(prob->type == ELEMENTAL)
		elem_block[ecnt] = lb->vertex2proc[ecnt];
	      else
		{
		  int proc   = lb->vertex2proc[mesh->connect[ecnt][0]];
		  int nnodes = get_elem_info(NNODES, mesh->elem_type[ecnt]);
		  elem_block[ecnt] = proc;
		  for(int ncnt=1; ncnt < nnodes; ncnt++) {
		    if(lb->vertex2proc[mesh->connect[ecnt][ncnt]] != proc) {
		      elem_block[ecnt] = machine->num_procs;
		      break;
		    }
		  }
		}
	    }

	  int ccnt = 0;
	  for(INT bcnt=0; bcnt < vis_nelem_blks; bcnt++) {
	      vis_el_blk_ptr[bcnt] = ccnt;
	      int pos = 0;
	      int old_pos = 0;
	      INT* el_ptr = elem_block;
	      size_t ecnt   = mesh->num_elems;
	      while(pos != -1)
		{
		  pos = in_list(bcnt, ecnt, el_ptr);
		  if(pos != -1)
		    {
		      old_pos += pos + 1;
		      ecnt     = mesh->num_elems - old_pos;
		      el_ptr   = elem_block + old_pos;
		      int nnodes = get_elem_info(NNODES, mesh->elem_type[old_pos-1]);
		      for(int ncnt=0; ncnt < nnodes; ncnt++)
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
	  for(int bcnt=0; bcnt < vis_nelem_blks; bcnt++)
	    {
	      /*
	       * Note this assumes all the blocks contain the same type
	       * element.
	       */
	      int ecnt = (vis_el_blk_ptr[bcnt+1]-vis_el_blk_ptr[bcnt])/node_pel_blk[0];
	      if(ex_put_elem_block(exid_vis, bcnt+1, elem_type[0],
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
	  if(elem_type)
	    free(elem_type);

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
	  float *xptr = NULL;
	  float *yptr = NULL;
	  float *zptr = NULL;
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
	  if(ex_put_coord_names(exid_vis, (char**)coord_names) < 0)
	    {
	      Gen_Error(0, "fatal: unable to output coordinate names");
	      ex_close(exid_vis);
	      return 0;
	    }

	  size_t nsize_old = 0;
	  INT *tmp_connect = NULL;
	  for(size_t ecnt=0; ecnt < mesh->num_el_blks; ecnt++)
	    {
	      nsize = el_cnt_blk[ecnt] * node_pel_blk[ecnt] * sizeof(INT);
	      if(nsize > nsize_old)
		{
		  if(nsize_old == 0)
		    tmp_connect = (INT*)malloc(nsize);
		  else
		    {
		      tmp_connect = (INT*)realloc(tmp_connect, nsize);
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

	      if(ex_put_elem_block(exid_vis, el_blk_ids[ecnt], elem_type[ecnt],
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

	  if(elem_type)
	    free(elem_type);

	  /* Allocate memory for the nodal values */
	  float *node_vals = (float*)malloc(mesh->num_nodes * sizeof(float));
	  if(!node_vals)
	    {
	      Gen_Error(0, "fatal: insufficient memory");
	      ex_close(exid_vis);
	      return 0;
	    }

	  /* Set up the file for nodal results */
	  float time_val = 0.0;
	  if(ex_put_time(exid_vis, 1, &time_val) < 0)
	    {
	      Gen_Error(0, "fatal: unable to output time to vis file");
	      ex_close(exid_vis);
	      return 0;
	    }
	  if(ex_put_variable_param(exid_vis, EX_NODAL, 1) < 0)
	    {
	      Gen_Error(0, "fatal: unable to output var params to vis file");
	      ex_close(exid_vis);
	      return 0;
	    }

	  const char  *var_names[] = {"proc"};
	  if(ex_put_variable_names(exid_vis, EX_NODAL, 1, (char**)var_names) < 0)
	    {
	      Gen_Error(0, "fatal: unable to output variable name");
	      ex_close(exid_vis);
	      return 0;
	    }

	  /* Do some problem specific assignment */
	  if(prob->type == NODAL)
	    {
	      for(size_t ncnt=0; ncnt < mesh->num_nodes; ncnt++)
		node_vals[ncnt] = lb->vertex2proc[ncnt];

	      for(int pcnt=0; pcnt < machine->num_procs; pcnt++)
		{
		  for(size_t ncnt=0; ncnt < lb->bor_nodes[pcnt].size(); ncnt++)
		    node_vals[lb->bor_nodes[pcnt][ncnt]] = machine->num_procs + 1;
		}

	    }
	  else if(prob->type == ELEMENTAL)
	    {
	      for(int pcnt=0; pcnt < machine->num_procs; pcnt++)
		{
		  for(size_t ncnt=0; ncnt < lb->int_nodes[pcnt].size(); ncnt++)
		    node_vals[lb->int_nodes[pcnt][ncnt]] = pcnt;

		  for(size_t ncnt=0; ncnt < lb->bor_nodes[pcnt].size(); ncnt++)
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

