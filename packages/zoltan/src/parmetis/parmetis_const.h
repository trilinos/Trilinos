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

extern int LB_Set_ParMetis_Param(char *, char *);
extern int LB_hashf(LB_GID, int);
extern int LB_hash_lookup (struct LB_hash_node **, LB_GID, int);


/* ParMETIS data types and definitions. */

/* Undefine the following #define in order to use short as the idxtype.
 * NB: Make sure these defs are consistent with those in your 
 * ParMetis installation !
*/
#define IDXTYPE_INT

#ifdef IDXTYPE_INT
typedef int idxtype;
#define IDX_DATATYPE    MPI_INT
#else
typedef short idxtype;
#define IDX_DATATYPE    MPI_SHORT
#endif

/* ParMetis 2.0 function prototypes */
extern void ParMETIS_PartKway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
extern void ParMETIS_PartGeomKway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, int *, idxtype *, MPI_Comm *);
extern void ParMETIS_PartGeom(idxtype *, int *, float *, idxtype *, MPI_Comm *);
extern void ParMETIS_RefineKway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
extern void ParMETIS_RepartLDiffusion(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
extern void ParMETIS_RepartGDiffusion(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
extern void ParMETIS_RepartRemap(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
extern void ParMETIS_RepartMLRemap(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);

