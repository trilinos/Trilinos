/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
#ifndef lint
static char *cvs_dr_exoII_io = "$Id$";
#endif

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
/* Author(s):  Matthew M. St.John (9226)                                    */
/*--------------------------------------------------------------------------*/
/* Revision History:                                                        */
/*                                                                          */
/*    30 March 1999:    Date of creation                                    */
/*--------------------------------------------------------------------------*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>

#include "exodusII.h"
#include "ne_nemesisI.h"

#include "dr_const.h"
#include "dr_input_const.h"
#include "dr_elem_const.h"
#include "dr_util_const.h"
#include "dr_par_util_const.h"
#include "dr_err_const.h"
#include "dr_output_const.h"

#define LIST_ALLOC 10

static int read_elem_info(int, int, PROB_INFO_PTR, ELEM_INFO *);
static int find_surnd_elem(ELEM_INFO *, int **, int *, int *);
static int find_adjacency(int, ELEM_INFO *, int **, int *, int);
static int read_comm_map_info(int, int, PROB_INFO_PTR, ELEM_INFO *);

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int read_exoII_mesh(int Proc,
                    int Num_Proc,
                    PROB_INFO_PTR prob,
                    PARIO_INFO_PTR pio_info,
                    ELEM_INFO **elements)
{
#ifdef LB_NO_NEMESIS
  Gen_Error(0, "Fatal:  Nemesis requested but not linked with driver.");
  return 0;

#else /* !LB_NO_NEMESIS */
  /* Local declarations. */
  char   par_nem_fname[FILENAME_MAX+1], title[MAX_LINE_LENGTH+1];
  char   cmesg[256];

  float  ver;

  int    i, pexoid, cpu_ws = 0, io_ws = 0;
#ifdef DEBUG_EXO
  int    j, k, elem;
#endif
  FILE  *fdtmp;

/***************************** BEGIN EXECUTION ******************************/

  /* since this is a test driver, set error reporting in exodus */
  ex_opts(EX_VERBOSE | EX_DEBUG);

  /* generate the parallel filename for this processor */
  gen_par_filename(pio_info->pexo_fname, par_nem_fname, pio_info, Proc,
                   Num_Proc);

  /* 
   * check whether parallel file exists.  do the check with fopen 
   * as ex_open coredumps on the paragon when files do not exist.
   */

  if ((fdtmp = fopen(par_nem_fname, "r")) == NULL) {
    sprintf(cmesg,"fatal: parallel Exodus II file %s does not exist",
            par_nem_fname);
    Gen_Error(0, cmesg);
    return 0;
  }
  else
    fclose(fdtmp);

  /*
   * now open the existing parallel file using Exodus calls.
   */

  if ((pexoid = ex_open(par_nem_fname, EX_READ, &cpu_ws, &io_ws,
                        &ver)) < 0) {
    sprintf(cmesg,"fatal: could not open parallel Exodus II file %s",
            par_nem_fname);
    Gen_Error(0, cmesg);
    return 0;
  }

  /* and get initial information */
  if (ex_get_init(pexoid, title, &(Mesh.num_dims),
                  &(Mesh.num_nodes), &(Mesh.num_elems),
                  &(Mesh.num_el_blks), &(Mesh.num_node_sets),
                  &(Mesh.num_side_sets)) < 0) {
    Gen_Error(0, "fatal: Error returned from ex_get_init");
    return 0;
  }


  /* alocate some memory for the element blocks */
  Mesh.eb_ids = (int *) malloc (4 * Mesh.num_el_blks * sizeof(int));
  if (!Mesh.eb_ids) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }
  Mesh.eb_cnts = Mesh.eb_ids + Mesh.num_el_blks;
  Mesh.eb_nnodes = Mesh.eb_cnts + Mesh.num_el_blks;
  Mesh.eb_nattrs = Mesh.eb_nnodes + Mesh.num_el_blks;

  Mesh.eb_names = (char **) malloc (Mesh.num_el_blks * sizeof(char *));
  if (!Mesh.eb_names) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }

  if (ex_get_elem_blk_ids(pexoid, Mesh.eb_ids) < 0) {
    Gen_Error(0, "fatal: Error returned from ex_get_elem_blk_ids");
    return 0;
  }

  /* get the element block information */
  for (i = 0; i < Mesh.num_el_blks; i++) {

    /* allocate space for name */
    Mesh.eb_names[i] = (char *) malloc((MAX_STR_LENGTH+1) * sizeof(char));
    if (!Mesh.eb_names[i]) {
      Gen_Error(0, "fatal: insufficient memory");
      return 0;
    }

    if (ex_get_elem_block(pexoid, Mesh.eb_ids[i], Mesh.eb_names[i],
                          &(Mesh.eb_cnts[i]), &(Mesh.eb_nnodes[i]),
                          &(Mesh.eb_nattrs[i])) < 0) {
      Gen_Error(0, "fatal: Error returned from ex_get_elem_block");
      return 0;
    }

  }

  /*
   * allocate memory for the elements
   * allocate a little extra for element migration latter
   */
  Mesh.elem_array_len = Mesh.num_elems + 5;
  *elements = (ELEM_INFO_PTR) malloc (Mesh.elem_array_len * sizeof(ELEM_INFO));
  if (!(*elements)) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }

  /*
   * intialize all of the element structs as unused by
   * setting the globalID to -1
   */
  for (i = 0; i < Mesh.elem_array_len; i++) {
    (*elements)[i].globalID = -1;
    (*elements)[i].coord = NULL;
    (*elements)[i].connect = NULL;
    (*elements)[i].adj = NULL;
    (*elements)[i].adj_proc = NULL;
    (*elements)[i].edge_wgt = NULL;
  }

  /* read the information for the individual elements */
  if (!read_elem_info(pexoid, Proc, prob, *elements)) {
    Gen_Error(0, "fatal: Error returned from read_elem_info");
    return 0;
  }

  /* read the communication information */
  if (!read_comm_map_info(pexoid, Proc, prob, *elements)) {
    Gen_Error(0, "fatal: Error returned from read_comm_map_info");
    return 0;
  }

  /* Close the parallel file */
  if(ex_close (pexoid) < 0) {
    Gen_Error(0, "fatal: Error returned from ex_close");
    return 0;
  }

#ifdef DEBUG_EXO

  /* print out the distributed mesh */
  print_distributed_mesh(Proc, Num_Proc, prob, *elements);

#endif

  return 1;

#endif /* LB_NO_NEMESIS */
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

#ifndef LB_NO_NEMESIS

static int read_elem_info(int pexoid, int Proc, PROB_INFO_PTR prob,
                          ELEM_INFO elements[])
{
  /* Local declarations. */
  int    iblk, ielem, inode, lnode, cnode, iplace, len;
  int    max_nsur;
  int   *nmap, *emap, *connect;
  int  **sur_elem, *nsurnd;

  float *xptr = NULL, *yptr = NULL, *zptr = NULL;
/***************************** BEGIN EXECUTION ******************************/


  /* allocate memory for the global number maps */
  nmap = (int *) malloc ((Mesh.num_nodes + Mesh.num_elems) * sizeof(int));
  if (!nmap) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }
  emap = nmap + Mesh.num_nodes;


  /*
   * get the global maps
   */
  if (ex_get_elem_num_map(pexoid, emap) < 0) {
    Gen_Error(0, "fatal: Error returned from ex_get_elem_num_map");
    return 0;
  }

  if (ex_get_node_num_map(pexoid, nmap) < 0) {
    Gen_Error(0, "fatal: Error returned from ex_get_node_num_map");
    return 0;
  }

  /* allocate memory for the coordinates */
  if (prob->read_coord) {
    xptr = (float *) malloc (Mesh.num_dims * Mesh.num_nodes * sizeof(float));
    if (!xptr) {
      Gen_Error(0, "fatal: insufficient memory");
      return 0;
    }
    switch (Mesh.num_dims) {
      case 3:
        zptr = xptr + 2 * Mesh.num_nodes;
        /* FALLTHRU */
      case 2:
        yptr = xptr + Mesh.num_nodes;
    }

    if (ex_get_coord(pexoid, xptr, yptr, zptr) < 0) {
      Gen_Error(0, "fatal: Error returned from ex_get_coord");
      return 0;
    }
  }

  /*
   * figure out which element block needs
   * the most space for its connect table
   */
  len = 0;
  for (iblk = 0; iblk < Mesh.num_el_blks; iblk++)
    if ((iplace = Mesh.eb_cnts[iblk] * Mesh.eb_nnodes[iblk]) > len)
      len = iplace;

  connect = (int *) malloc (len * sizeof(int));
  if (!connect) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }

  /***************************************************************************/
  /* Fill the Connect table, Coordinates, Global Ids for each element        */
  /***************************************************************************/
  iplace = 0;
  for (iblk = 0; iblk < Mesh.num_el_blks; iblk++) {

    if (Mesh.eb_cnts[iblk] > 0) {
      if (ex_get_elem_conn(pexoid, Mesh.eb_ids[iblk], connect) < 0) {
        Gen_Error(0, "fatal: Error returned from ex_get_elem_conn");
        return 0;
      }

      cnode = 0;
      for (ielem = 0; ielem < Mesh.eb_cnts[iblk]; ielem++) {
        /* set some fields in the element structure */
        elements[iplace].border = 0;
        elements[iplace].globalID = emap[iplace];
        elements[iplace].elem_blk = iblk;
        elements[iplace].nadj = 0;
        /* weights are 1 for now */
        elements[iplace].cpu_wgt = 1.0;
        elements[iplace].mem_wgt = 1.0;

        /* allocate space for the connect list and the coordinates */
        elements[iplace].connect = (int *) malloc(Mesh.eb_nnodes[iblk] *
                                                  sizeof(int));
        if (!(elements[iplace].connect)) {
          Gen_Error(0, "fatal: insufficient memory");
          return 0;
        }

        if (prob->read_coord) {
          elements[iplace].coord = (float **) malloc(Mesh.eb_nnodes[iblk] *
                                                     sizeof(float *));
          if (!(elements[iplace].coord)) {
            Gen_Error(0, "fatal: insufficient memory");
            return 0;
          }
        }
        else
          elements[iplace].coord = NULL;

        /* save the connect table as local numbers for the moment */
        for (inode = 0; inode < Mesh.eb_nnodes[iblk]; inode++) {
          lnode = connect[cnode] - 1;
          elements[iplace].connect[inode] = lnode;
          cnode++;

          if (prob->read_coord) {
            elements[iplace].coord[inode] = (float *) malloc(Mesh.num_dims *
                                                             sizeof(float));
            if (!(elements[iplace].coord[inode])) {
              Gen_Error(0, "fatal: insufficient memory");
              return 0;
            }

            switch (Mesh.num_dims) {
              case 3:
                elements[iplace].coord[inode][2] = zptr[lnode];
                /* FALLTHRU */
              case 2:
                elements[iplace].coord[inode][1] = yptr[lnode];
                /* FALLTHRU */
              case 1:
                elements[iplace].coord[inode][0] = xptr[lnode];
            }
          }
        } /* End: "for (inode = 0; inode < Mesh.eb_nnodes[iblk]; inode++)" */

        iplace++;

      } /* End: "for (ielem = 0; ielem < Mesh.eb_cnts[iblk]; ielem++)" */
    } /* End: "if (Mesh.eb_cnts[iblk] > 0)" */
  } /* End: "for (iblk = 0; iblk < Mesh.num_el_blks; iblk++)" */

  /* free some memory */
  free(connect);
  if (prob->read_coord) free(xptr);

  if (prob->gen_graph) {
    /*************************************************************************/
    /* Find the adjacency list for each element                              */
    /*	Part one: find the surrounding elements for each node                */
    /*************************************************************************/
    sur_elem = (int **) malloc(Mesh.num_nodes * sizeof(int *));
    if (!sur_elem) {
      Gen_Error(0, "fatal: insufficient memory");
      return 0;
    }
    nsurnd = (int *) malloc(Mesh.num_nodes * sizeof(int));
    if (!nsurnd) {
      Gen_Error(0, "fatal: insufficient memory");
      return 0;
    }

    if (!find_surnd_elem(elements, sur_elem, nsurnd, &max_nsur)) {
      Gen_Error(0, "fatal: Error returned from find_surnd_elems");
      return 0;
    }

    /*************************************************************************/
    /*	Part two: Find the adjacencies on this processor                     */
    /*		and get the edge weights                                     */ 
    /*************************************************************************/
    if (!find_adjacency(Proc, elements, sur_elem, nsurnd, max_nsur)) {
      Gen_Error(0, "fatal: Error returned from find_adjacency");
      return 0;
    }

    /*
     * convert the node numbers in the connect lists to Global IDs
     * since they will be much easier to work with
     */
    for (ielem = 0; ielem < Mesh.num_elems; ielem++) {
      iblk = elements[ielem].elem_blk;
      for (inode = 0; inode < Mesh.eb_nnodes[iblk]; inode++) {
        elements[ielem].connect[inode] = nmap[elements[ielem].connect[inode]];
      }
    }

    for (inode = 0; inode < Mesh.num_nodes; inode++) free(sur_elem[inode]);
    free(sur_elem);
    free(nsurnd);
  }

  free(nmap);

  return 1;
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static int find_surnd_elem(ELEM_INFO elements[], int **sur_elem, int *nsurnd,
                           int *max_nsur)
{
  /* Local declarations. */
  int     ielem, inode, lnode;
  int    *alloc_cnt, *tmp_ptr;
/***************************** BEGIN EXECUTION ******************************/

  alloc_cnt = (int *) malloc(Mesh.num_nodes * sizeof(int));
  if (!alloc_cnt) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }

  /* Allocate rows of the structure */
  for(inode=0; inode < Mesh.num_nodes; inode++)
  {
    sur_elem[inode] = (int *) malloc(LIST_ALLOC*sizeof(int));
    if(!(sur_elem[inode])) {
      Gen_Error(0, "fatal: insufficient memory");
      return 0;
    }
    alloc_cnt[inode] = LIST_ALLOC;
    nsurnd[inode] = 0;
  }

  /* Find the surrounding elements for each node in the mesh */
  for(ielem=0; ielem < Mesh.num_elems; ielem++) {
    for(inode=0; inode < Mesh.eb_nnodes[elements[ielem].elem_blk]; inode++) {
      lnode = elements[ielem].connect[inode];

      /*
       * in the case of degenerate elements, where a node can be
       * entered into the connect table twice, need to check to
       * make sure that this element is not already listed as
       * surrounding this node
       */
      if (nsurnd[lnode] > 0 &&
            ielem == sur_elem[lnode][(nsurnd[inode]) - 1])
        continue;

      nsurnd[lnode]++;

      /* keep track of the largest number of surrounding elements */
      if (nsurnd[lnode] > *max_nsur) *max_nsur = nsurnd[lnode];

      /* Check to see if this rows size should be increased */
      if(nsurnd[lnode] > alloc_cnt[lnode]) {
        alloc_cnt[lnode] += LIST_ALLOC;
        tmp_ptr = (int *) realloc(sur_elem[lnode],
                                  (alloc_cnt[lnode]) * sizeof(int));
        if(!tmp_ptr) {
          Gen_Error(0, "fatal: insufficient memory");
          return 0;
        }
        sur_elem[lnode] = tmp_ptr;
      }

      /* Add the element to the list */
      sur_elem[lnode][nsurnd[lnode]-1] = ielem;
    }

  } /* End "for(ielem=0; ielem < Mesh.num_elems; ielem++)" */

  free (alloc_cnt);

  return 1;
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static int find_adjacency(int Proc, ELEM_INFO elements[],
                          int **sur_elem, int *nsurnd, int max_nsur)
{
  /* Local declarations. */
  int     i, iblk, nsides, ielem, nscnt, inode, entry;
  int     side_cnt, nnodes, sid;
  int     side_nodes[MAX_SIDE_NODES], mirror_nodes[MAX_SIDE_NODES];
  int    *hold_elem, *pt_list, nhold, nelem;

  E_Type *eb_etype;
/***************************** BEGIN EXECUTION ******************************/
  /*
   * Use face definition of adjacencies. So, one elements that are
   * connected by an entire face will be considered adjacent. This
   * is temporary, and will be expanded. This will make determining
   * off processor adjacencies much easier.
   */

  /* allocate memory and determine the element type for each element block */
  eb_etype = (E_Type *) malloc (Mesh.num_el_blks * sizeof(E_Type));
  if(!eb_etype) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }

  for (iblk = 0; iblk < Mesh.num_el_blks; iblk++) {
    if (Mesh.eb_cnts[iblk] > 0) {
      if ((eb_etype[iblk] =  get_elem_type(Mesh.eb_names[iblk],
                                           Mesh.eb_nnodes[iblk],
                                           Mesh.num_dims)) < 0) {
        Gen_Error(0, "fatal: could not get element type");
        return 0;
      }
    }
    else eb_etype[iblk] = NULL_EL;
  }

  /* allocate space to hold info about surounding elements */
  pt_list = (int *) malloc(2 * max_nsur * sizeof(int));
  if(!pt_list) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }
  hold_elem = pt_list + max_nsur;

  for (ielem = 0; ielem < Mesh.num_elems; ielem++) {
    elements[ielem].adj = (int *) malloc(LIST_ALLOC*sizeof(int));
    elements[ielem].adj_proc = (int *) malloc(LIST_ALLOC*sizeof(int));
    elements[ielem].edge_wgt = (float *) malloc(LIST_ALLOC*sizeof(float));
    if(!(elements[ielem].adj) || !(elements[ielem].edge_wgt) ||
       !(elements[ielem].adj_proc)) {
      Gen_Error(0, "fatal: insufficient memory");
      return 0;
    }
    /* NOTE: nadj set in read_elem_info in case graph not generated */
    elements[ielem].adj_len = LIST_ALLOC;

    iblk = elements[ielem].elem_blk;

    /* exclude circle and sphere elements from graph */
    if (Mesh.eb_nnodes[iblk] > 1) {

      if ((nsides = get_elem_info(NSIDES, eb_etype[iblk], 0)) < 0) {
        Gen_Error(0, "fatal: could not get element information");
        return 0;
      }

      /* check each side of this element */
      for (nscnt = 0; nscnt < nsides; nscnt++) {

        /* get the list of nodes on this side set */
        side_cnt = ss_to_node_list(eb_etype[iblk], elements[ielem].connect,
                                   (nscnt+1), side_nodes);

        /*
         * now I need to determine how many side set nodes I
         * need to use to determine if there is an element
         * connected to this side.
         *
         * 2-D - need two nodes, so find one intersection
         * 3-D - need three nodes, so find two intersections
         * NOTE: must check to make sure that this number is not
         *       larger than the number of nodes on the sides (ie - SHELL).
         */
        nnodes = Mesh.num_dims;
        if (side_cnt < nnodes)   nnodes = side_cnt;
        nnodes--;      /* decrement to find the number of intersections  */

        nelem = 0;     /* reset this in case no intersections are needed */

        /* copy the first array into temp storage */
        nhold = nsurnd[side_nodes[0]];
        for (i = 0; i < nhold; i++)
          hold_elem[i] = sur_elem[side_nodes[0]][i];

        for (inode = 0; inode < nnodes; inode++) {
          nelem = find_inter(hold_elem, sur_elem[side_nodes[(inode+1)]],
                             nhold, nsurnd[side_nodes[(inode+1)]], 2, pt_list);

          if (nelem < 2) break;
          else {
            nhold = nelem;
            for (i = 0; i < nelem; i++)
              hold_elem[i] = hold_elem[pt_list[i]];
          }
        }

        /*
         * if there is an element on this side of ielem, then there
         * will be at least two elements in the intersection (one
         * will be ielem)
         */
        if (nelem > 1) {

          /*
           * now go through and check each element in the list
           * to see if it is different than ielem.
           */
          for(i=0; i < nelem; i++) {

            entry = hold_elem[i];

            if(entry != ielem) {
              /*
               * get the side id of entry. Make sure that ielem is
               * trying to communicate to a valid side of elem
               */
              side_cnt = get_ss_mirror(eb_etype[iblk], side_nodes, (nscnt+1),
                                       mirror_nodes);

              /*
               * in order to get the correct side order for elem,
               * get the mirror of the side of ielem
               */
              sid = get_side_id(eb_etype[elements[entry].elem_blk],
                                elements[entry].connect,
                                side_cnt, mirror_nodes);
              if (sid > 0) {
                (elements[ielem].nadj)++;
                if(elements[ielem].nadj > elements[ielem].adj_len) {
                  elements[ielem].adj_len += LIST_ALLOC;
                  elements[ielem].adj = (int *) realloc(elements[ielem].adj,
                                        elements[ielem].adj_len * sizeof(int));
                  elements[ielem].adj_proc = (int *) realloc(
                                        elements[ielem].adj_proc,
                                        elements[ielem].adj_len * sizeof(int));
                  elements[ielem].edge_wgt = (float *) realloc(
                                      elements[ielem].edge_wgt,
                                      elements[ielem].adj_len * sizeof(float));
                  if(!(elements[ielem].adj) || !(elements[ielem].adj_proc) ||
                     !(elements[ielem].edge_wgt)) {
                    Gen_Error(0, "fatal: insufficient memory");
                    return 0;
                  }
                }
                elements[ielem].adj[(elements[ielem].nadj)-1] = entry;
                elements[ielem].adj_proc[(elements[ielem].nadj)-1] = Proc;

                /*
                 * the edge weight is the number of nodes in the
                 * connecting face
                 */
                elements[ielem].edge_wgt[(elements[ielem].nadj)-1] = 
                  (float) get_elem_info(NSNODES, eb_etype[iblk], (nscnt+1));

              } /* End: "if (sid > 0)" */
              else if (sid < 0) {
                Gen_Error(0, "fatal: could not find side id");
                return 0;
              }
            } /* End: "if(ielem != entry)" */
          } /* End: "for(i=0; i < nelem; i++)" */
        } /* End: "if (nelem > 1)" */
      } /* End: "for (nscnt = 0; ...)" */
    } /* End: "if (nnode > 1)" */
  } /* End: "for (ielem=0; ...)" */

  free(pt_list);
  free(eb_etype);

  return 1;
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static int read_comm_map_info(int pexoid, int Proc, PROB_INFO_PTR prob,
                              ELEM_INFO elements[])
{
  /* Local declarations. */
  int  ielem, imap, loc_elem, iblk, max_len, offset, index;
  int  nnodei, nnodeb, nnodee, nelemi, nelemb, nncmap, necmap;
  int *int_elem, *bor_elem, *ecmap_cnt, *ecmap_id;
  int *elem_ids, *side_ids, *proc_ids, *neigh_ids;

  E_Type etype;

  MPI_Status status;
/***************************** BEGIN EXECUTION ******************************/

  /*
   * Currently, the communication map information is only used to
   * determine the graph for this mesh. So, if this method does not
   * need a graph generated, return from this function. In the future
   * communication map information will probably be needed for other
   * things.
   */
  if (!prob->gen_graph) return 1;

  if (ne_get_loadbal_param(pexoid, &nnodei, &nnodeb, &nnodee,
                           &nelemi, &nelemb, &nncmap, &necmap, Proc) < 0) {
    Gen_Error(0, "fatal: Error returned from ne_get_loadbal_param");
    return 0;
  }

  /*
   * get the list of the border elements in order to set
   * the border flag in the element structures
   */
  int_elem = (int *) malloc ((nelemi + nelemb) * sizeof(int));
  if (!int_elem) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }
  bor_elem = int_elem + nelemi;

  if (ne_get_elem_map(pexoid, int_elem, bor_elem, Proc) < 0) {
    Gen_Error(0, "fatal: Error returned from ne_get_elem_map");
    return 0;
  }

  for (ielem = 0; ielem < nelemb; ielem++) {
    elements[ielem-1].border = 1;
  }

  free(int_elem);

  /*
   * For now, only get the elemental communication maps,
   * since, in the driver, elements are only considered
   * adjacent if they share a face (same definition used
   * in element communication maps). Eventually, the ability
   * to consider elements that are connected by any nodes
   * adjacent will have to be added. When that happens,
   * the nodal communication maps will be needed.
   */
  ecmap_cnt = (int *) malloc (2 * necmap * sizeof(int));
  if (!ecmap_cnt) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }
  ecmap_id = ecmap_cnt + necmap;

  if (ne_get_cmap_params(pexoid, NULL, NULL, ecmap_id, ecmap_cnt, Proc) < 0) {
    Gen_Error(0, "fatal: Error returned from ne_get_cmap_params");
    return 0;
  }

  max_len = 0;
  for (imap = 0; imap < necmap; imap++)
    max_len += ecmap_cnt[imap];

  elem_ids = (int *) malloc(4 * max_len * sizeof(int));
  if (!elem_ids) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }
  side_ids = elem_ids + max_len;
  proc_ids = side_ids + max_len;
  neigh_ids = proc_ids + max_len;

  offset = 0;
  for (imap = 0; imap < necmap; imap++) {

    if(ne_get_elem_cmap(pexoid, ecmap_id[imap], &(elem_ids[offset]),
                        &(side_ids[offset]), &(proc_ids[offset]), Proc) < 0) {
      Gen_Error(0, "fatal: Error returned from ne_get_elem_cmap");
      return 0;
    }
    offset += ecmap_cnt[imap];
  } /* End: "for (imap = 0; imap < necmap; imap++)" */

  /*
   * convert the element ids to global ids to send to
   * the neighboring processor. Store them in the proc_ids
   * array, since the proc_id information is not used for
   * anything here.
   */
  for (ielem = 0; ielem < max_len; ielem++)
    proc_ids[ielem] = elements[elem_ids[ielem]-1].globalID;

#ifdef DEBUG_ALL
  printf("\nCommunication maps for Proc %d\n", Proc);
  printf("Number of maps: %d\n", necmap);
  printf("Map Counts:");
  for (imap = 0; imap < necmap; imap++)
    printf(" %d", ecmap_cnt[imap]);
  printf("\n");
  printf("Map Ids:");
  for (imap = 0; imap < necmap; imap++)
    printf(" %d", ecmap_id[imap]);
  printf("\n");
  printf("elem side globalID\n");
  for (ielem = 0; ielem < max_len; ielem++)
    printf("%d   %d   %d\n", elem_ids[ielem], side_ids[ielem], proc_ids[ielem]);
#endif

  /*
   * Now communicate with other processor to get global IDs
   * for the adjacent elements in this communication map.
   *
   * parallel nemesis trick...
   * each communication map is only for a single neigboring
   * processor, and the communication map id is the processor
   * number that it is for
   */
  offset = 0;
  for (imap = 0; imap < necmap; imap++) {

    /*
     * handshake with processor and wait until it is ready
     * to talk
     */
    MPI_Send(NULL, 0, MPI_INT, ecmap_id[imap], 0, MPI_COMM_WORLD);
    MPI_Recv(NULL, 0, MPI_INT, ecmap_id[imap], 0, MPI_COMM_WORLD, &status);

    /* now send list of global element ids to the processor for this map */
    MPI_Send(&(proc_ids[offset]), ecmap_cnt[imap], MPI_INT, ecmap_id[imap], 0,
             MPI_COMM_WORLD);
    MPI_Recv(&(neigh_ids[offset]), ecmap_cnt[imap], MPI_INT, ecmap_id[imap],
             0, MPI_COMM_WORLD, &status);
    offset += ecmap_cnt[imap];
  }

  /* now process all of the element ids that have been received */
  offset = 0;
  for (imap = 0; imap < necmap; imap++) {
    for (ielem = 0; ielem < ecmap_cnt[imap]; ielem++) {
      index = ielem + offset;
      /* translate from element id in the communication map to local elem id */
      loc_elem = elem_ids[index] - 1;
      iblk = elements[loc_elem].elem_blk;
      etype = get_elem_type(Mesh.eb_names[iblk], Mesh.eb_nnodes[iblk],
                            Mesh.num_dims);

      (elements[loc_elem].nadj)++;
      if(elements[loc_elem].nadj > elements[loc_elem].adj_len) {
        elements[loc_elem].adj_len += LIST_ALLOC;
        elements[loc_elem].adj = (int *) realloc(elements[loc_elem].adj,
                                     elements[loc_elem].adj_len * sizeof(int));
        elements[loc_elem].adj_proc = (int *) realloc(
                                     elements[loc_elem].adj_proc,
                                     elements[loc_elem].adj_len * sizeof(int));
        elements[loc_elem].edge_wgt = (float *) realloc(
                                   elements[loc_elem].edge_wgt,
                                   elements[loc_elem].adj_len * sizeof(float));
        if(!(elements[loc_elem].adj) || !(elements[loc_elem].adj_proc) ||
           !(elements[loc_elem].edge_wgt)) {
          Gen_Error(0, "fatal: insufficient memory");
          return 0;
        }
      }

      elements[loc_elem].adj[elements[loc_elem].nadj-1] = neigh_ids[index];
      elements[loc_elem].adj_proc[elements[loc_elem].nadj-1] = ecmap_id[imap];
      elements[loc_elem].edge_wgt[elements[loc_elem].nadj-1] =
                     (float) get_elem_info(NSNODES, etype, side_ids[index]);

    } /* End: "for (ielem = 0; ielem < ecmap_cnt[imap]; ielem++)" */
    offset += ecmap_cnt[imap];
  } /* End: "for for (imap = 0; imap < necmap; imap++)" */
  
  free (ecmap_cnt);
  free (elem_ids);

  return 1;
}
#endif /* !LB_NO_NEMESIS */
