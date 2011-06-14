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

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Functions contained in this file:
 *	main()
 *	print_input()
 *	print_usage()
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef SGI10K
#include <sys/resource.h>
#endif

#include "elb_const.h"
#include "elb_err_const.h"
#include "elb_inp_const.h"
#include "elb_exo_const.h"
#include "elb_allo_const.h"
#include "elb_graph_const.h"
#include "elb_loadbal_const.h"
#include "elb_output_const.h"

extern void add_to_log(const char *name);
/* Local function prototypes */
void print_usage(void);

char* elem_names[NULL_EL] = {"SPHERE", "BAR2", "BAR3", "QUAD4", "QUAD8", "QUAD9",
			     "SHELL4", "SHELL8", "TRI3", "TRI6", "TSHELL3", "TSHELL6", "HEX8",
			     "HEX20", "HEX27", "HEXSHELL", "TET4", "TET10", "TET8", "WEDGE6", "WEDGE15",
			     "WEDGE16", "PYRAMID5", "PYRAMID13", "SHELL2", "SHELL3"};


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Routine which reads in a EXODUS II mesh database and load balances the
 * mesh by either a nodal or element assignment.
 *----------------------------------------------------------------------------
 * Functions called:
 *	cmd_line_arg_parse()
 *	print_usage()
 *	read_cmd_file()
 *	check_inp_specs()
 *	print_input()
 *	read_exo_weights()
 *	read_mesh_params()
 *      read_mesh()
 *	generate_gaph()
 *	generate_loadbal()
 *	write_vis()
 *	generate_maps()
 *	write_nemesis()
 *----------------------------------------------------------------------------
 * Variable Index:
 *	machine       - structure containing information about the machine
 *			for which the load balance is to be generated.
 *	lb	      - structure containing information about what type
 *			of load balance is to be performed.
 *	problem	      - structure containing information about the problem
 *			type. Currently whether the problem is a nodal or
 *			elemental based decomposition.
 *	solver	      - structure containing various parameters for the
 *			eigensolver used by Chaco.
 *	weight	      - structure used to store parameters about how to
 *			weight the resulting graph before it is fed into
 *			Chaco.
 *	graph	      - structure containing the graph of the problem.
 *	mesh	      - structure containing a description of the FEM mesh.
 *	exoII_inp_file - name of the ExodusII file containing the problem
 *			 geometry.
 *	ascii_inp_file - name of the input command file.
 *	nemI_out_file  - name of the output NemesisI file.
 *
 *****************************************************************************/
int main (int argc, char *argv[])
{
  /* Machine dependent variables */
#ifdef SGI10K
  struct rlimit64 rlp64;
#endif

  /* Local variables */
  int    cnt, cnt1;
  char   exoII_inp_file[MAX_FNL+1], ascii_inp_file[MAX_FNL+1];
  char   nemI_out_file[MAX_FNL+1];
  double time1, time2, start_time, end_time;

  MACHINE_TYPE machine;
  LB_INFO      lb;
  PROB_INFO    problem;
  SOLVE_INFO   solver;
  WEIGHT_INFO  weight;
  MESH_INFO    mesh;
  SPHERE_INFO  sphere;
  GRAPH_INFO   graph;

  /* Local function prototypes */
  void    print_input(MACHINE_PTR,
                      LB_INFO_PTR,
                      PROB_INFO_PTR,
                      SOLVE_INFO_PTR,
                      WEIGHT_INFO_PTR);

/*-----------------------------Execution Begins------------------------------*/

  start_time = get_time();

#ifdef SGI10K
  printf("Setting up for large memory usage on 64bit SGI...");
  fflush(stdout);
  rlp64.rlim_cur = RLIM64_INFINITY;
  rlp64.rlim_max = RLIM64_INFINITY;

  if(setrlimit64(RLIMIT_VMEM, &rlp64) != 0)
    fprintf(stderr, "Warning: failed to set virt memory usage to maximum\n");

  if(setrlimit64(RLIMIT_DATA, &rlp64) != 0)
    fprintf(stderr, "Warning: failed to set data memory usage to maximum\n");

  printf("done\n");
#endif

  /* Initialization */
  graph.adj       = NULL;
  graph.sur_elem  = NULL;
  graph.max_nsur  = 0;
  graph.alloc_cnt = NULL;
  graph.nsur_elem = NULL;
  graph.start     = NULL;

  machine.num_boxes = -1;
  machine.type      = -1;
  machine.num_dims  = -1;
  machine.dim[0]    = -1;

  lb.type       = -1;
  lb.refine     = -1;
  lb.num_sects  = -1;
  lb.cnctd_dom  = -1;
  lb.outfile    = -1;

  solver.tolerance = -1.0;
  solver.rqi_flag  = -1;
  solver.vmax  = -1;

  problem.type             = -1;
  problem.read_coords      = -1;
  problem.coarse_flag      = -1;
  problem.alloc_graph      = -1;
  problem.vis_out          = -1;
  problem.skip_checks      = -1;
  problem.face_adj         = -1;
  problem.partial_adj      =  0;
  problem.global_mech      = -1;
  problem.local_mech       = -1;
  problem.find_cnt_domains = -1;
  problem.mech_add_procs   = -1;
  problem.dsd_add_procs    = -1;
  problem.no_sph           = -1;
  problem.groups           = NULL;
  problem.group_no         = NULL;
  problem.num_groups       = -1;

  weight.type            = -1;
  weight.ow_read         = 0;
  weight.exo_filename[0] = '\0';
  weight.exo_varname[0]  = '\0';
  weight.exo_tindx       = -1;
  weight.exo_vindx       = -1;
  weight.num_ebw         = 0;
  weight.elemblk         = NULL;
  weight.elemblk_wgt     = NULL;
  weight.ow              = NULL;
  weight.vertices        = NULL;
  weight.edges           = NULL;

  sphere.num             = 0;
  sphere.adjust          = NULL;
  sphere.begin           = NULL;
  sphere.end             = NULL;

  mesh.title[0]          = '\0';

  /* Initialize file names */
  ascii_inp_file[0] = '\0';
  exoII_inp_file[0] = '\0';
  nemI_out_file[0]  = '\0';

  /* Initialize to just reporting the error */
  error_lev      = 1;

  /* check if the user just wants to know the version */
  for(cnt=0; cnt < argc; cnt++) {
    if(strcmp(argv[cnt],"-V") == 0) {
      printf("%s version %s\n", UTIL_NAME, ELB_VERSION);
      exit(0);
    }
  }

  /* Parse the command line */
  if(!cmd_line_arg_parse(argc, argv,
                         exoII_inp_file,
                         ascii_inp_file,
                         nemI_out_file,
                         &machine,
                         &lb,
                         &problem,
                         &solver,
                         &weight))
  {
    fprintf(stderr, "error parsing command line\n");
    error_report();
    exit(1);
  }

  /*
   *If the information is to be taken from an ASCII input file then
   * read that file.
   */
  if(ascii_inp_file[0] != '\0')
  {
    if(!read_cmd_file(ascii_inp_file, exoII_inp_file, nemI_out_file,
                      &machine, &lb, &problem, &solver, &weight))
    {
      fprintf(stderr, "error parsing command file\n");
      error_report();
      exit(1);
    }
  }

  /* make sure that this type is set */
  if (weight.type < 0) weight.type = NO_WEIGHT;

  /*
   * Perform at least some rudimentary error checks on the user
   * specified input.
   */
  if(!check_inp_specs(exoII_inp_file, nemI_out_file, &machine, &lb,
                      &problem, &solver, &weight))
  {
    fprintf(stderr, "Error in user specified parameters\n");
    error_report();
    exit(1);
  }

  /* Output the parameters for the run to the screen */
  print_input(&machine, &lb, &problem, &solver, &weight);

  /* Read in mesh parameters */
  time1 = get_time();
  if(!read_mesh_params(exoII_inp_file, &problem, &mesh, &sphere))
  {
    fprintf(stderr, "Error reading mesh parameters\n");
    error_report();
    exit(1);
  }
  time2 = get_time();
  printf("Time to read mesh parameters: %fs\n", time2-time1);

  /* Check for error conditions */
  if((mesh.num_nodes)/(machine.num_procs) < 1)
  {
    Gen_Error(0, "fatal: problem divided among too many processors");
    error_report();
    exit(1);
  }
  else if ((problem.type == ELEMENTAL) &&
            ((mesh.num_elems)/(machine.num_procs) < 1))
  {
    Gen_Error(0, "fatal: problem divided among too many processors");
    error_report();
    exit(1);
  }

  /* If vertex weighting is turned on, prepare weight struct */
  if ((weight.type & READ_EXO) || (weight.type & EL_BLK))
  {
    time1 = get_time();
    if(!init_weight_struct(&problem, &mesh, &weight))
    {
      fprintf(stderr, "Error during initialization of weight struct\n");
      error_report();
      exit(1);
    }
    time2 = get_time();
    printf("Time to initialize weight struct: %fs\n", time2-time1);
  }

  /* If desired, read in the weighting factors from the ExodusII file */
  if(weight.type & READ_EXO)
  {
    time1 = get_time();
    if(!read_exo_weights(&problem, &weight))
    {
      fprintf(stderr, "Error during read of ExodusII weights\n");
      error_report();
      exit(1);
    }
    time2 = get_time();
    printf("Time to read ExodusII weights: %fs\n", time2-time1);
  }

  /* Initialize various parameters */
  if(lb.type == INERTIAL || lb.type == ZPINCH ||
     lb.type == BRICK || lb.type == ZOLTAN_RCB ||
     lb.type == ZOLTAN_RIB || lb.type == ZOLTAN_HSFC ||
     problem.vis_out == ELB_TRUE)
    problem.read_coords = ELB_TRUE;
  else
    problem.read_coords = ELB_FALSE;

  if(lb.type != SPECTRAL)
    problem.coarse_flag = ELB_FALSE;
  else
    problem.coarse_flag = ELB_TRUE;

  if(lb.refine == KL_REFINE)
    problem.alloc_graph = ELB_TRUE;
  else if(lb.type == SPECTRAL)
    problem.alloc_graph = ELB_TRUE;
  else
    problem.alloc_graph = ELB_FALSE;

  /* Allocate necessary memory */
  if(problem.type == NODAL)
    problem.num_vertices = mesh.num_nodes;
  else if(problem.type == ELEMENTAL)
    problem.num_vertices = (mesh.num_elems - sphere.num);

  if(problem.read_coords == ELB_TRUE)
  {
    mesh.coords = malloc((mesh.num_dims)*(mesh.num_nodes)*sizeof(float));
    if(!(mesh.coords))
    {
      Gen_Error(0, "fatal: insufficient memory for coordinates");
      error_report();
      exit(1);
    }
  }
  else
    mesh.coords = NULL;

  mesh.elem_type = (E_Type *)array_alloc(1, mesh.num_elems, sizeof(E_Type));
  mesh.connect   = (int **)array_alloc(2, mesh.num_elems,
                                       mesh.max_np_elem, sizeof(int));
  if(!(mesh.elem_type) || !(mesh.connect))
  {
    Gen_Error(0, "fatal: insufficient memory");
    error_report();
    exit(1);
  }

  /* Read the mesh */
  time1 = get_time();
  if(!read_mesh(exoII_inp_file, &problem, &mesh, &weight))
  {
    Gen_Error(0, "fatal: could not read mesh");
    error_report();
    exit(1);
  }
  time2 = get_time();
  printf("Time to read mesh: %fs\n", time2-time1);

  /* free unneeded memory */
  if (weight.type & READ_EXO)
    if (weight.ow) free (weight.ow);

  if ((weight.type & EL_BLK) && weight.num_ebw > 0) {
    if (weight.elemblk)     free (weight.elemblk);
    if (weight.elemblk_wgt) free (weight.elemblk_wgt);
  }

  /* Generate the graph for the mesh */
  time1 = get_time();
  if(!generate_graph(&problem, &mesh, &graph, &weight, &sphere))
  {
    Gen_Error(0, "fatal: could not generate graph");
    error_report();
    exit(1);
  }
  time2 = get_time();
  printf("Time to generate graph: %fs\n", time2-time1);

  /* Generate load balance */
  time1 = get_time();
  if(!generate_loadbal(&machine, &problem, &mesh, &lb, &solver, &graph,
                       &weight, &sphere, argc, argv))
  {
    Gen_Error(0, "fatal: could not generate load balance");
    error_report();
    exit(1);
  }
  time2 = get_time();
  printf("Time to generate load balance: %fs\n", time2-time1);

  /* free up memory */
  if (sphere.adjust) free(sphere.adjust);


#ifdef PRINT_VERT
  for (cnt=0; cnt < problem.num_vertices; cnt++)
    printf("element = %i, proc = %i\n", cnt, lb.vertex2proc[cnt]);
#endif

  /*
   * NOTE: in Chaco, if FREE_GRAPH is set to 1, the following arrays
   * are freed: graph.start
   *		graph.adj
   *		weight.vertices
   *		weight.edges
   *
   * Need to take into account special case where the mesh contains
   * only spheres. In this case Chaco is not called, and the above
   * arrays are not freed
   */

  if (sphere.num >= mesh.num_elems) {
    if (graph.start)     free (graph.start);
    if (graph.adj)       free (graph.adj);
    if (weight.vertices) free (weight.vertices);
    if (weight.edges)    free (weight.edges);
  }


  /* Generate the load balance maps */
  time1 = get_time();
  if(!generate_maps(&machine, &problem, &mesh, &lb, &graph))
  {
    Gen_Error(0, "fatal: could not generate load-balance maps");
    error_report();
    exit(1);
  }
  time2 = get_time();
  printf("Time to generate load-balance maps: %fs\n", time2-time1);

  /* Output the visualization file */
  if(problem.vis_out == ELB_TRUE)
  {
    time1 = get_time();
    if(!write_vis(nemI_out_file, exoII_inp_file, &machine, &problem,
                  &mesh, &lb))
    {
      Gen_Error(0, "warning: unable to output visualization file");
    }
    time2 = get_time();
    printf("Time to output visualization file: %fs\n", time2-time1);
  }

  /* Free up memory no longer needed */
  if(problem.read_coords == ELB_TRUE)
    free(mesh.coords);

  free(mesh.elem_type);
  free(mesh.connect);

  if(graph.sur_elem)
  {
    for(cnt=0; cnt < mesh.num_nodes; cnt++)
      free(graph.sur_elem[cnt]);

    free(graph.sur_elem);

    free(graph.nsur_elem);
  }

  /* Output a Nemesis load balance file */
  time1 = get_time();
  if(!write_nemesis(nemI_out_file, &machine, &problem, &mesh, &lb, &sphere))
  {
    Gen_Error(0, "fatal: could not output Nemesis file");
    error_report();
    exit(1);
  }
  time2 = get_time() - time1;
  printf("Time to write Nemesis file: %fs\n", time2);

  /* Free up unused memory for leak checking */
  for(cnt=0; cnt < machine.num_procs; cnt++)
  {
    free(lb.int_nodes[cnt]);
    free(lb.bor_nodes[cnt]);
    if(problem.type == NODAL)
    {
      free(lb.ext_nodes[cnt]);
      free(lb.ext_procs[cnt]);
    }

    free(lb.int_elems[cnt]);
    if(problem.type == ELEMENTAL)
      free(lb.bor_elems[cnt]);
  }

  free(lb.int_nodes);

  if(problem.type == NODAL)
  {
    free(lb.ext_nodes);
    free(lb.ext_procs);
  }

  free(lb.int_elems);

  if(problem.type == ELEMENTAL)
  {
    free(lb.bor_elems);
    for(cnt=0; cnt < machine.num_procs; cnt++)
    {
      for(cnt1=0; cnt1 < lb.num_bor_nodes[cnt]; cnt1++)
        free(lb.born_procs[cnt][cnt1]);

      if (lb.born_procs[cnt])     free(lb.born_procs[cnt]);
      if (lb.born_proc_cnts[cnt]) free(lb.born_proc_cnts[cnt]);
      free(lb.ext_procs[cnt]);
      free(lb.e_cmap_elems[cnt]);
      free(lb.e_cmap_sides[cnt]);
      free(lb.e_cmap_procs[cnt]);
      free(lb.e_cmap_neigh[cnt]);

    }
    free(lb.e_cmap_elems);
    free(lb.e_cmap_size);
    free(lb.born_procs);
    free(lb.born_proc_cnts);
    free(lb.ext_procs);

  }

  free(lb.num_int_nodes);
  free(lb.bor_nodes);
  free(lb.vertex2proc);

  /* Report any non-fatal errors that may have occured */
  error_report();

  /* Get ending time */
  end_time = get_time();
  printf("The entire load balance took %fs\n", end_time-start_time);
  add_to_log(argv[0]);
  return 0;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function prints out parameters as read from the command line and/or
 * the input ASCII load-balance file.
 *---------------------------------------------------------------------------*/
void print_input(
  MACHINE_PTR machine,
  LB_INFO_PTR lb,
  PROB_INFO_PTR prob,
  SOLVE_INFO_PTR solver,
  WEIGHT_INFO_PTR weight
  )
{
  int cnt;

/*-----------------------------Execution Begins------------------------------*/

  printf("%s version %s\n", UTIL_NAME, ELB_VERSION);

  printf("Performing ");
  switch(prob->type)
  {
  case NODAL:
    printf("a nodal ");
    break;

  case ELEMENTAL:
    printf("an elemental ");
    break;
  }

  printf("load balance with the following parameters...\n");

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                         MACHINE PARAMETERS                                */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  printf("Machine Description\n");
  if (machine->num_boxes > 1)
  {
    printf("\tCluster of %d boxes\n", machine->num_boxes);
    printf("\tarchitecture of each box: ");
  }
  else
    printf("\tarchitecture: ");
  switch(machine->type)
  {
  case HCUBE:
    printf("hypercube\n");
    break;

  case MESH:
    printf("mesh\n");
    break;
  }
  if (machine->num_boxes > 1)
    printf("\tdimension(s) of each box: ");
  else
    printf("\tdimension(s): ");
  switch(machine->type)
  {
  case HCUBE:
    printf("%d\n", machine->dim[0]);
    break;

  case MESH:
    for(cnt=0; cnt < (machine->num_dims)-1; cnt++)
      printf("%dx", machine->dim[cnt]);

    printf("%d\n", machine->dim[(machine->num_dims)-1]);
    break;
  }
  printf("\ttotal number of processors: %d\n", machine->num_procs);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                       LOAD BALANCE PARAMETERS                             */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  printf("Load Balance Parameters\n");
  switch(lb->type)
  {
  case MULTIKL:
    printf("\ttype: multilevel\n");
    printf("\tnumber of sections: %d\n", lb->num_sects);
    break;

  case SPECTRAL:
    printf("\ttype: spectral\n");
    printf("\tnumber of sections: %d\n", lb->num_sects);
    break;

  case INERTIAL:
    printf("\ttype: inertial\n");
    break;

  case ZPINCH:
    printf("\ttype: zpinch\n");
    break;

  case BRICK:
    printf("\ttype: brick\n");
    break;

  case ZOLTAN_RCB:
    printf("\ttype: rcb\n");
    break;

  case ZOLTAN_RIB:
    printf("\ttype: rib\n");
    break;

  case ZOLTAN_HSFC:
    printf("\ttype: hsfc\n");
    break;

  case LINEAR:
    printf("\ttype: linear\n");
    break;

  case RANDOM:
    printf("\ttype: random\n");
    break;

  case SCATTERED:
    printf("\ttype: scattered\n");
    break;

  case INFILE:
    printf("\ttype: input from file\n");
    printf("\tfile name: %s\n", lb->file);
    break;
  }
  printf("\trefinement: ");
  switch(lb->refine)
  {
  case KL_REFINE:
    printf("Kernighan-Lin\n");
    break;

  case NO_REFINE:
    printf("none\n");
    break;
  }
  if (lb->cnctd_dom) {
    printf("\tConnected Domain enforced\n");
  }
  if (lb->outfile) {
    printf("\toutput assignment vector\n");
    printf("\tfile name: %s\n", lb->file);
  }

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                         EIGENSOLVER PARAMETERS                            */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if(lb->type == MULTIKL || lb->type == SPECTRAL)
  {
    printf("Eigensolver Parameters\n");
    printf("\teignsolver tolerance: %f\n", solver->tolerance);
    if(solver->rqi_flag == USE_RQI)
    {
      printf("\tusing RQI/Symmlq eigensolver\n");
      printf("\tnumber of vertices to coarsen down to: %d\n",
             solver->vmax);
    }
  }

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                          WEIGHTING PARAMETERS                             */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  printf("Weighting Parameters\n");
  if(weight->type == NO_WEIGHT)
    printf("\tno weighting\n");
  else if(weight->type & READ_EXO)
  {
    printf("\tweights from: ExodusII file\n");
    printf("\ttime index: %d\n", weight->exo_tindx);
    printf("\tvariable index: %d\n", weight->exo_vindx);
    printf("\tvariable name: %s\n", weight->exo_varname);
  }
  else if(weight->type & EL_BLK)
  {
    printf("\tElement Block weights specified\n");
    for (cnt=0; cnt < weight->num_ebw; cnt++)
      printf("\telement block: %d, weight: %d\n", weight->elemblk[cnt],
             weight->elemblk_wgt[cnt]);
  }
  else if(weight->type & EDGE_WGT)
    printf("\tedge weights turned on\n");

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                          WEIGHTING PARAMETERS                             */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  printf("Miscellaneous Options\n");
  if(prob->face_adj == 1)
    printf("\tusing face definition of adjacency\n");
  if(prob->global_mech == 1)
    printf("\tlooking for global mechanisms\n");
  if(prob->local_mech == 1)
    printf("\tlooking for local mechanisms\n");
  if(prob->find_cnt_domains == 1)
    printf("\tidentifying the number of disconnected element blocks on a subdomain\n");
  if(prob->mech_add_procs == 1)
    printf("\tincreasing the number of processors if mechanisms are found\n");
  if(prob->dsd_add_procs == 1)
    printf("\tincreasing the number of processors if disconnected sudomains are found\n");
  if(prob->no_sph == 1)
    printf("\tSPHERES are being treated as concentrated mass - connectivity exists\n");
  if (prob->skip_checks == 1)
    printf("\tWARNING: side id error checks turned off\n");
  if (prob->groups != NULL)
  {
    printf("\telement block groups defined\n");
    printf("\tgroup string: \"%s\"\n", prob->groups);
  }
  return;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function prints out usage information for the utility.
 *---------------------------------------------------------------------------*/
void print_usage(void)
{

  printf("\n");
  printf("usage:\t%s [-h] [<-n|-e> -o <output file>", UTIL_NAME);
  printf(" -m <machine description>\n");
  printf("\t -l <load bal description> -s <eigen solver specs>\n");
  printf("\t -w <weighting options> -g <group list> -f]\n");
  printf("\t [-a <ascii file>] exoII_file\n\n");
  printf(" -n\t\tperform a nodal based load balance\n");
  printf(" -e\t\tperform an elemental based load balance\n");
  printf(" -o NemI file\toutput NemesisI load-balance file name\n");
  printf(" -m mach desc\tdescription of the machine to load balance for\n");
  printf(" -l LB data\tload balance specifications\n");
  printf(" -s Eigen specs\tEigen solver specifications\n");
  printf(" -w weighting\tweighting specs for load balance\n");
  printf(" -g groupings\tgrouping specifications for load balance\n");
  printf(" -f\t\tuse face definition of adjacency\n");
  printf(" -p\t\tuse partial definition of adjacency: \n");
  printf("   \t\trequire only 3 matching quad face nodes\n");
  printf("\t\t  OR\n");
  printf(" -h\t\tusage information\n");
  printf("\t\t  OR\n");
  printf(" -a ascii file\tget info from ascii input file name\n");

  return;

} /*----------End print_usage()-------------*/
