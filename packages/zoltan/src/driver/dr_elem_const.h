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
 * $Name$
 *====================================================================*/

#ifndef _DR_ELM_CONST_H
#define _DR_ELM_CONST_H
#ifndef lint
static char *cvs_elemch_id = "$Id$";
#endif

/* Define element types */
typedef enum {E_TYPE_ERROR=-1, SPHERE, BAR1, BAR2, QUAD1, S_QUAD2, QUAD2,
              SHELL1, SHELL2, TRI1, TRI2, TSHELL1, TSHELL2, HEX1,
              S_HEX2, HEX2, TET1, TET2, WEDGE1, WEDGE2,
              HEXSHELL, NULL_EL} E_Type;

extern
E_Type get_elem_type(
  const char *elem_name,	/* ExodusII element name */
  const int   num_nodes,	/* Number of nodes in the element */
  const int   num_dim		/* Number of dimensions of the mesh */
);

extern
int get_elem_info(
  const int info,		/* The requested information */
  const E_Type elem_type,	/* The element type */
  const int sid			/* side id (to get number of side nodes) */
);

extern
int get_side_id(
  const E_Type  etype,		/* The element type */
  const int *conn,		/* The element connectivity */
  const int  nsnodes,		/* The number of side nodes */
  int  side_nodes[]		/* The list of side node IDs */
);

extern
int ss_to_node_list(
  const E_Type  etype,		/* The element type */
  const int *connect,		/* The element connectivity */
  int  side_num,		/* The element side number */
  int  ss_node_list[]		/* The list of side node IDs */
);

extern
int get_ss_mirror(
  const E_Type etype,		/* The element type */
  const int *ss_node_list,	/* The list of side node IDs */
  int side_num,			/* The element side number */
  int mirror_node_list[]	/* The list of the mirror side node IDs */
);


/* Define element info requests */
#define NNODES		0
#define NQUAD		1
#define NDIM		2
#define NQUAD_SURF	3
#define NSNODES		4
#define NSIDES		5

/* Define for the maximum number of nodes on an element side/face */
#define MAX_SIDE_NODES	9
/*
 * Define for the maximum number of sides (and hence communications
 * entries) that an element can have
 */
#define MAX_ELEM_SIDES	6

#endif /* _DR_ELM_CONST_H */
