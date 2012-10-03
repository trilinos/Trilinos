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
 *	generate_graph()
 *	find_surnd_elems()
 *	find_adjacency()
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include <sstream>
#include <assert.h>                     // for assert
#include <stdio.h>                      // for sprintf, printf, NULL
#include <stdlib.h>                     // for realloc, malloc, free
#include <string.h>                     // for strcat, strcpy
#include "elb_graph.h"
#include "elb.h"                  // for Graph_Description<INT>, etc
#include "elb_elem.h"             // for E_Type, get_elem_info, etc
#include "elb_err.h"              // for Gen_Error
#include "elb_util.h"             // for in_list, find_inter

extern int is_hex(E_Type etype);
extern int is_tet(E_Type etype);
extern int is_3d_element(E_Type etype);

/* Local function prototypes */
namespace {
  template <typename INT>
  int  find_surnd_elems(Mesh_Description<INT>*, Graph_Description<INT>*);

  template <typename INT>
  int  find_adjacency(Problem_Description*, Mesh_Description<INT>*, Graph_Description<INT>*,
		      Weight_Description<INT>*, Sphere_Info*);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function generate_graph() begins:
 *----------------------------------------------------------------------------
 * This function does the work to generate the graph from the FE mesh.
 *****************************************************************************/
template int generate_graph(Problem_Description* problem,
			    Mesh_Description<int>* mesh,
			    Graph_Description<int>* graph,
			    Weight_Description<int>* weight,
			    Sphere_Info* sphere);

template int generate_graph(Problem_Description* problem,
			    Mesh_Description<int64_t>* mesh,
			    Graph_Description<int64_t>* graph,
			    Weight_Description<int64_t>* weight,
			    Sphere_Info* sphere);

template <typename INT>
int generate_graph(Problem_Description* problem,
                   Mesh_Description<INT>* mesh,
                   Graph_Description<INT>* graph,
                   Weight_Description<INT>* weight,
                   Sphere_Info* sphere)
{
  double time1 = get_time();
  /* Find the elements surrounding a node */
  if(!find_surnd_elems(mesh, graph)) {
    Gen_Error(0, "fatal: could not find surrounding elements");
    return 0;
  }
  double time2 = get_time();
  printf("Time to find surrounding elements: %fs\n", time2-time1);

  /* Find the adjacency, if required */
  if(problem->alloc_graph == ELB_TRUE) {
    if(!find_adjacency(problem, mesh, graph, weight, sphere)) {
      Gen_Error(0, "fatal: could not find adjacency");
      return 0;
    }
    time1 = get_time();
    printf("Time to find the adjacency: %fs\n", time1-time2);
  }
  return 1;
}

namespace {
  /*****************************************************************************/
  /*****************************************************************************/
  /*****************************************************************************/
  /* Function find_surnd_elems() begins:
   *----------------------------------------------------------------------------
   * This function finds the elements surrounding a given FEM node. In other
   * words, this function generates a list of elements containing a given
   * FEM node.
   *****************************************************************************/
  template <typename INT>
  int find_surnd_elems(Mesh_Description<INT>* mesh,
		       Graph_Description<INT>* graph)
  {
    graph->sur_elem.resize(mesh->num_nodes);

    /* Find the surrounding elements for each node in the mesh */
    for(size_t ecnt=0; ecnt < mesh->num_elems; ecnt++) {
      int nnodes = get_elem_info(NNODES, mesh->elem_type[ecnt]);
      for(int ncnt=0; ncnt < nnodes; ncnt++) {
	INT node = mesh->connect[ecnt][ncnt];

	/*
	 * in the case of degenerate elements, where a node can be
	 * entered into the connect table twice, need to check to
	 * make sure that this element is not already listed as
	 * surrounding this node
	 */
	if (graph->sur_elem[node].empty() ||
	    ecnt != (size_t)graph->sur_elem[node][graph->sur_elem[node].size()-1]) {
	  /* Add the element to the list */
	  graph->sur_elem[node].push_back(ecnt);
	}
      }
    } /* End "for(ecnt=0; ecnt < mesh->num_elems; ecnt++)" */

    for(size_t ncnt=0; ncnt < mesh->num_nodes; ncnt++) {
      if(graph->sur_elem[ncnt].empty()) {
	printf("WARNING: Node = %lu has no elements\n", ncnt+1);
      } else {
	size_t nsur = graph->sur_elem[ncnt].size();
	if (nsur > graph->max_nsur)
	  graph->max_nsur = nsur;
      }
    }
    return 1;
  }

  /*****************************************************************************/
  /*****************************************************************************/
  /*****************************************************************************/
  /* Function find_adjacency() begins:
   *----------------------------------------------------------------------------
   * This function finds adjacency (or graph) of the problem.
   *****************************************************************************/
  template <typename INT>
  int find_adjacency(Problem_Description* problem,
		     Mesh_Description<INT>* mesh,
		     Graph_Description<INT>* graph,
		     Weight_Description<INT>* weight,
		     Sphere_Info* sphere)
  {
    int     iret;

    size_t nelem = 0;
    size_t nhold = 0;
    int sid = 0;
    
    INT *pt_list=NULL;
    INT *hold_elem=NULL;
    INT side_nodes[MAX_SIDE_NODES];
    INT mirror_nodes[MAX_SIDE_NODES];

    static int count = 0;

    int     hflag1, hflag2, tflag1, tflag2;
    /*-----------------------------Execution Begins------------------------------*/

    /* Allocate memory necessary for the adjacency */
    graph->start.resize(problem->num_vertices+1);
    
    for (int i=0; i < MAX_SIDE_NODES; i++) {
      side_nodes[i]=-999;
      mirror_nodes[i]=-999;
    }
  
    /* Find the adjacency for a nodal based decomposition */
    if(problem->type == NODAL) {
      graph->nadj = 0;
      for(size_t ncnt=0; ncnt < mesh->num_nodes; ncnt++) {
	graph->start[ncnt] = graph->nadj;
	for(size_t ecnt=0; ecnt < graph->sur_elem[ncnt].size(); ecnt++) {
	  size_t elem   = graph->sur_elem[ncnt][ecnt];
	  int nnodes = get_elem_info(NNODES, mesh->elem_type[elem]);
	  for(int i=0; i < nnodes; i++) {
	    INT entry = mesh->connect[elem][i];
	    
	    if(ncnt != (size_t)entry &&
	       in_list(entry,
		       graph->adj.size()-graph->start[ncnt],
		       &graph->adj[graph->start[ncnt]]) < 0) {
	      graph->adj.push_back(entry);
	    }
	  }
	} /* End "for(ecnt=0; ecnt < graph->nsur_elem[ncnt]; ecnt++)" */
      } /* End "for(ncnt=0; ncnt < mesh->num_nodes; ncnt++)" */
    }
    /* Find the adjacency for a elemental based decomposition */
    else {

      /* for face adjacencies, need to allocate some memory */
      if (problem->face_adj) {
	/* allocate space to hold info about surounding elements */
	pt_list   = (INT*)malloc(sizeof(INT)*graph->max_nsur);
	hold_elem = (INT*)malloc(sizeof(INT)*graph->max_nsur);
	if(!(pt_list) || !(hold_elem)) {
	  Gen_Error(0, "fatal: insufficient memory");
	  return 0;
	}
      }
      graph->nadj = 0;
      size_t cnt = 0;

      /* tmp_element used to speed up the in_list calc */
      int *tmp_element = (int*)malloc(sizeof(int)*mesh->num_elems);
      if(!tmp_element) {
	Gen_Error(0, "fatal: insufficient memory");
	return 0;
      }
      for(size_t ecnt=0; ecnt < mesh->num_elems; ecnt++)
	tmp_element[ecnt] = -1;

    
      /* Cycle through the elements */
      E_Type etype_last = NULL_EL;
      E_Type etype = NULL_EL;

      int element_3d = 0;
      int nnodes = mesh->num_dims; 
      int nsides = 0;
    
      for(size_t ecnt=0; ecnt < mesh->num_elems; ecnt++) {
	etype = mesh->elem_type[ecnt];
	if (etype != etype_last) {
	  etype_last = etype;
	  element_3d = is_3d_element(mesh->elem_type[ecnt]);
	  if (problem->face_adj == 0) {
	    nnodes = get_elem_info(NNODES, etype);
	  }
	  nsides = get_elem_info(NSIDES, etype);
	}
	  
	if(etype != SPHERE || (etype == SPHERE && problem->no_sph == 1)) {
	  graph->start[cnt] = graph->nadj;

	  /*
	   * now have to decide how to determine adjacency
	   * !face_adj - any element that connects to any node in this
	   *             element is an adjacent element
	   * face_adj - a) for 3D elements only those that share an 
	   *               entire face with this element are considered 
	   *               adjacent
	   *            b) do not connect 1D/2D elements to 3D elements
	   *            c) 1D and 2D elements can connect to each other
	   */

	  /* If not forcing face adjaceny */
	  if (problem->face_adj == 0) {

	    /* ncnt = 0,...,7 for hex */
	    for(int ncnt=0; ncnt < nnodes; ncnt++) {
	      /* node is the node number 'ncnt' of element 'ecnt' */
	      size_t node = mesh->connect[ecnt][ncnt];

	      /* i varies from 0 -> # of elements touching 'node' */
	      for(size_t i=0; i < graph->sur_elem[node].size(); i++){

		/* 'entry' is element # i touching node 'node' */
		INT entry = graph->sur_elem[node][i];

		/* make sure we're not checking if the element
		   is connected to itself */
		if(ecnt != (size_t)entry && mesh->elem_type[entry] != SPHERE) {
		  /* If tmp_element[entry] != ecnt, then entry is not in list... */
		  if ((size_t)tmp_element[entry] != ecnt) {
#if 0
		    assert(in_list(entry, (graph->nadj)-(graph->start[cnt]),
				   (graph->adj)+(graph->start[cnt])) < 0);
#endif
		    tmp_element[entry] = ecnt;
		    (graph->nadj)++;
		    graph->adj.push_back(entry);
		    if (weight->type & EDGE_WGT)
		      weight->edges.push_back(1.0);
		  }
		  else if (weight->type & EDGE_WGT) {
		    iret=in_list(entry, (graph->nadj)-(graph->start[cnt]),
				 &graph->adj[graph->start[cnt]]);
		    assert(iret >= 0);
		    weight->edges[iret+(graph->start[cnt])] += 1.0;
		  }
		}
	      }
	    } /* End "for(ncnt=0; ...)" */
	  }  /* End: "if (problem->face_adj == 0)" */
        
	  /* So if this is a 3-d element and we're forcing face
	   * adjacency, if it gets to this else below 
	   *
	   * if this element is 1d/2d allow connections to 1d and 2d 
	   * elements but not to 3d elements
	   *
	   */

	  else {
	    if(element_3d) {
	      /* need to check for hex's or tet's */

	      /*
	       * If the first element is a hex or tet, set flags
	       * hflag1/tflag1 to 1
	       */
	      hflag1 = is_hex(etype);
	      tflag1 = is_tet(etype);

	      /* check each side of this element */
	      for (int nscnt = 0; nscnt < nsides; nscnt++) {
		/* get the list of nodes on this side set */
		int side_cnt = ss_to_node_list(etype, mesh->connect[ecnt], (nscnt+1),
					       side_nodes);

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

		nnodes = mesh->num_dims;

		/*
		 * In case the number of nodes on this side are less 
		 * than the minimum number, set nnodes to side_cnt,
		 * i.e., if a 3-D mesh contains a bar, then nnodes=3, 
		 * and side_cnt = 2
		 */

		if (nnodes > side_cnt)   
		  nnodes = side_cnt;

		nnodes--;    /* decrement to find the number of intersections  */

		nelem = 0;     /* reset this in case no intersections are needed */

		/* copy the first array into temp storage */

#if 0
		/* nhold is the number of elements touching node 0 on
		   the side of this element */
		size_t nhold = graph->sur_elem[side_nodes[0]].size();

		/* Now that we have the number of elements touching
		   side 0, get their element ids and store them in hold_elem */
		for (size_t ncnt = 0; ncnt < nhold; ncnt++)
		  hold_elem[ncnt] = graph->sur_elem[side_nodes[0]][ncnt];
#endif
		/*
		 * need to handle hex's differently because of
		 * the tet/hex combination
		 */
              

		if (!hflag1) {
		  /* Get the number of elements ( and their ids )
		     that touch node (ncnt+1) and see if any elements touch 
		     both node 0 and node (ncnt+1), and if so, return to nelem
		     the number of elements touching both nodes and their 
		     indices in pt_list.  When ncnt != 0, hold_elem and nhold 
		     change */
		  nhold = graph->sur_elem[side_nodes[0]].size();
		  for (size_t ncnt = 0; ncnt < nhold; ncnt++)
		    hold_elem[ncnt] = graph->sur_elem[side_nodes[0]][ncnt];

		  for (int ncnt = 0; ncnt < nnodes; ncnt++) {
		    nelem = find_inter(hold_elem,
				       &graph->sur_elem[side_nodes[(ncnt+1)]][0],
				       nhold,
				       graph->sur_elem[side_nodes[(ncnt+1)]].size(),
				       pt_list);
		    
		    /*  If less than 2 ( 0 or 1 ) elements only
			touch nodes 0 and ncnt+1 then try next side node, i.e., 
			repeat loop ncnt */
		    if (nelem < 2)
		      break;
		    else {
		      nhold = nelem;
		      for (size_t i = 0; i < nelem; i++)
			hold_elem[i] = hold_elem[pt_list[i]];
		    }
		  }
		}

		/* If this element is a hex type */
		else {

		  /*
		   * To handle hex's, check opposite corners. First check
		   * 1 and 3 and then 2 and 4. Only an element connected
		   * to this face will be connected to both corners. If there
		   * are tet's connected to this face, both will show up in
		   * one of the intersections (nothing will show up in the
		   * other intersection).
		   */

		  /* See if hexes share nodes 0 and nodes (ncnt+2) */
		  int inode = 0;
		  for (int ncnt = 0; ncnt < nnodes; ncnt++) {
		    nelem = find_inter(&graph->sur_elem[side_nodes[inode]][0],
				       &graph->sur_elem[side_nodes[(ncnt+2)]][0],
				       graph->sur_elem[side_nodes[inode]].size(),
				       graph->sur_elem[side_nodes[(ncnt+2)]].size(),
				       pt_list);

		    /*
		     * If there are multiple elements in the intersection, then
		     * they must share the face, since the intersection is between
		     * the corner nodes. No element could connect with both of
		     * those nodes without being connected elsewhere.
		     */
		    if (nelem > 1) {

		      /* Then get the correct elements out of the hold array */
		      nhold = nelem;
		      for (size_t i = 0; i < nelem; i++)
			hold_elem[i] = graph->sur_elem[side_nodes[inode]][pt_list[i]];
		      break;
		    }
		    else {
		      /*
		       * if there aren't multiple elements in the intersection,
		       * then check the opposite corners (1 & 3)
		       */
		      inode = 1;
		    }
		  }
		} /* "if (!hflag)" */

		/*
		 * if there is an element on this side of ecnt, then there
		 * will be at least two elements in the intersection (one
		 * will be ecnt)
		 */
		if (nelem > 1) {

		  /*
		   * now go through and check each element in the list
		   * to see if it is different than ecnt.
		   */

		  for(size_t i=0; i < nelem; i++) {
		    size_t entry = hold_elem[i];

		    if(ecnt != entry) {

		      /*
		       * Need to verify that this side of ecnt is actually
		       * connected to a face of entry. The problem case is
		       * when an entire face of a shell (one of the ends)
		       * is connected to only an edge of a quad/tet
		       */

		      E_Type etype2 = mesh->elem_type[entry];

		      /* make sure this is a 3d element*/

		      if (is_3d_element(etype2)) {

			/* need to check for hex's */
			hflag2 = is_hex(etype2);

			/* TET10 cannnot connect to a HEX */
			tflag2 = is_tet(etype2);

			/* check here for tet/hex combinations */
			if ((tflag1 && hflag2) || (hflag1 && tflag2)) {
			  /*
			   * have to call a special function to get the side id
			   * in these cases. In both cases, the number of side
			   * nodes for the element will not be consistent with
			   * side_cnt, and:
			   *
			   * TET/HEX - side_nodes only contains three of the
			   *           the side nodes of the hex.
			   *
			   * HEX/TET - Have to check that this tet shares a side
			   *           with the hex.
			   */
			  sid = get_side_id_hex_tet(mesh->elem_type[entry],
						    mesh->connect[entry],
						    side_cnt, side_nodes);
			}
			else {
			  /*
			   * get the side id of elem. Make sure that ecnt is
			   * trying to communicate to a valid side of elem
			   */

			  side_cnt = get_ss_mirror(etype, side_nodes, (nscnt+1),
						   mirror_nodes);

			  /*
			   * small kludge to handle 6 node faces butted up against
			   * 4 node faces
			   */

			  /* if this element 1 is a hexshell, then only
			     require 4 of the 6 nodes to match between elements 
			     1 and 2 */
			  if (etype == HEXSHELL && side_cnt == 6) side_cnt = 4;

			  /* side_cnt is the number of nodes on the face
			     of a particular element.  This number is passed
			     to get_side_id and the error with two hexes
			     only sharing 3 nodes is in get_side_id 
			     Additional comments can be found there */
		      
			  /*
			   * in order to get the correct side order for elem,
			   * get the mirror of the side of ecnt
			   */

			  /* Based on elements intersecting, get the side
			     of element 1 that is connected to the element in the list
			     which it intersects with.  The two elements must have
			     (originally) side_cnt nodes in common */

			  sid = get_side_id(mesh->elem_type[entry],
						   mesh->connect[entry],
						   side_cnt, mirror_nodes,
						   problem->skip_checks,
						   problem->partial_adj);

			  /* printf("sid = %d\n",sid); */
			}

			if (sid > 0) {
			  (graph->nadj)++;
			  graph->adj.push_back(entry);
			  if (weight->type & EDGE_WGT) {
			    /*
			     * the edge weight is the number of nodes in the
			     * connecting face
			     */
			    weight->edges.push_back(side_cnt);

			    /*
			     * have to put a kluge in here for the
			     * tet/hex problem
			     */
			    if (hflag1 && tflag2)
			      (weight->edges[weight->edges.size()-1])--;
			  }
			}
			else if ((sid < 0) && (!problem->skip_checks)) {
			  /*
			   * too many errors with bad meshes, print out
			   * more information here for diagnostics
			   */
			  char tmpstr[80];
			  char   cmesg[256];
			  sprintf(cmesg,
				  "Error returned while getting side id for communication map.");
			  Gen_Error(0, cmesg);
			  sprintf(cmesg, "Element 1: %lu", (ecnt+1));
			  Gen_Error(0, cmesg);
			  nnodes = get_elem_info(NNODES, etype);
			  strcpy(cmesg, "connect table:");
			  for (int ii = 0; ii < nnodes; ii++) {
			    sprintf(tmpstr, " %lu", (size_t)(mesh->connect[ecnt][ii]+1));
			    strcat(cmesg, tmpstr);
			  }
			  Gen_Error(0, cmesg);
			  sprintf(cmesg, "side id: %d", (nscnt+1));
			  Gen_Error(0, cmesg);
			  strcpy(cmesg, "side nodes:");
			  for (int ii = 0; ii < side_cnt; ii++) {
			    sprintf(tmpstr, " %lu", (size_t)(side_nodes[ii]+1));
			    strcat(cmesg, tmpstr);
			  }
			  Gen_Error(0, cmesg);
			  sprintf(cmesg, "Element 2: %lu", (entry+1));
			  Gen_Error(0, cmesg);
			  nnodes = get_elem_info(NNODES, etype2);
			  strcpy(cmesg, "connect table:");
			  for (int ii = 0; ii < nnodes; ii++) {
			    sprintf(tmpstr, " %lu", (size_t)(mesh->connect[entry][ii]+1));
			    strcat(cmesg, tmpstr);
			  }
			  Gen_Error(0, cmesg);
			  count++;
			  printf("Now we have %d bad element connections.\n",count);
			} /* End "if (sid > 0)" */


		      } /* End: "if(ecnt != entry)" */
		    }

		  } /* End: "for(i=0; i < nelem; i++)" */

		} /* End: "if (nelem > 1)" */

	      } /* End: "for (nscnt = 0; nscnt < nsides; nscnt++)" */

	    } /* End: "if(element_3d)" */

	    else {

	      /* this is either a 2d or 1d element. Only allow attachments to other
	       * 1d or 2d elements
	       */

	      nnodes = get_elem_info(NNODES, mesh->elem_type[ecnt]);

	      for(int ncnt=0; ncnt < nnodes; ncnt++) {
		/* node is the node number 'ncnt' of element 'ecnt' */
		size_t node = mesh->connect[ecnt][ncnt];

		/* i varies from 0 -> # of elements touching 'node' */
		for(size_t i=0; i < graph->sur_elem[node].size(); i++) {

		  /* 'entry' is element # i touching node 'node' */
		  INT entry = graph->sur_elem[node][i];

		  /* make sure we're not checking if the element
		     is connected to itself */
		  if(ecnt != (size_t)entry) {

		    /* now make sure that the entry is not a 3d element */
		    if (!is_3d_element(mesh->elem_type[entry])) {
		      
		      if((iret=in_list(entry, graph->adj.size()-graph->start[cnt],
				       &graph->adj[graph->start[cnt]])) < 0) {

			(graph->nadj)++;
			graph->adj.push_back(entry);
			if (weight->type & EDGE_WGT)
			  weight->edges.push_back(1.0);
		      }
		      else if (weight->type & EDGE_WGT)
			weight->edges[iret+(graph->start[cnt])] += 1.0;
		    }
		  } /* End: if(ecnt != entry) */
		} /* for(i=0; i < graph->nsur_elem[node]; i++) */
	      } /* End "for(ncnt=0; ...)" */
	    } /* End: "else" (if !element_3d) */
	  } /* End: "else" (if face_adj != 0) */

	  cnt++;

	} /* End "if(etype != SPHERE)" */
      } /* End "for(ecnt=0; ecnt < mesh->num_elems; ecnt++)" */

      if (problem->face_adj) {
	free(hold_elem);
	free(pt_list);
      }
      free(tmp_element);
    }

    graph->start[problem->num_vertices] = graph->adj.size();
    graph->nadj = graph->adj.size();

    if ((size_t)graph->start[problem->num_vertices] != graph->nadj) {
      // Possibly an integer overflow... Output error message and stop.
      std::ostringstream errmsg;
      errmsg << "fatal: Graph adjacency edge count (" << graph->nadj << ") exceeds chaco 32-bit integer range.\n";
      Gen_Error(0, errmsg.str().c_str());
      return 0;
    }

    /* Adjust for a mesh with spheres */
    if(problem->type == ELEMENTAL && sphere->num) {
      /* Decrement adjacency entries */
      for(size_t cnt1=0; cnt1 < graph->adj.size(); cnt1++) {
	for(size_t ecnt=0; ecnt < mesh->num_el_blks; ecnt++) {
	  if(graph->adj[cnt1] >= sphere->begin[ecnt] &&
	     graph->adj[cnt1] < sphere->end[ecnt]) {
	    graph->adj[cnt1] -= sphere->adjust[ecnt];
	    break;
	  }
	}
      }
    }
    return 1;
  }
}
