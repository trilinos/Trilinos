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
static char *cvs_parmetis_const_h_id = "$Id$";
#endif

/* Data structures used in parmetis interface routines */

struct LB_vtx_list {
  int length;     /* Length of (remainder of ) list */
  LB_GID my_gid;     /* Global id of local vtx */
  int my_gno;     /* Global number of local vtx */
  LB_GID nbor_gid;   /* Global id of off-proc vtx */
  int * adj;      /* Pointer to adjcny array */
  struct LB_vtx_list * next;
};

struct LB_hash_node {
  LB_GID gid;  /* Global id */
  int gno;  /* Global number */
  struct LB_hash_node * next;
};

/* Function prototypes */

void LB_ParMETIS_Part(LB *lb, int *num_imp, LB_GID** imp_gids,
                       LB_LID** imp_lids, int **imp_procs);
int LB_hashf(LB_GID key, int n);
int LB_hash_lookup (struct LB_hash_node **hashtab, int n, LB_GID key);


/* ParMETIS data types. Must be consistent with the definitions
   in ParMetis/ParMETISlib.
*/

/* Undefine the following #define in order to use short int as the idxtype */
#define IDXTYPE_INT

/* Indexes are as long as integers for now */
#ifdef IDXTYPE_INT
typedef int idxtype;
#define IDX_DATATYPE    MPI_INT
#else
typedef short idxtype;
#define IDX_DATATYPE    MPI_SHORT
#endif

