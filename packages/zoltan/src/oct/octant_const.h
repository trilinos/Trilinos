#include "lb_const.h"
#include <math.h>

typedef struct Region_Node* pRegion;   /* typedef for a pointer to a region  */
typedef struct Region_Node {           /* region = area in 3-space           */
  struct Region_Node *next;            /* pointer to next region in list     */
  double Coord[3];                     /* centroid location of region        */
  double Weight;                       /* weight of Region - default is 1    */
  LB_TAG Tag;                          /* Tag containing IDs for the object  */
  int    attached;                     /* flag to see if region was attached */
} Region;

typedef struct _Octant* pOctant;     /* typedef for a pointer to an octant   */
typedef struct _Octant {             /* octant tree node that has 8 children */
  pRegion list;                      /* list of regions associated to octant */
  struct _Octant *child[8];          /* array of octant's children           */
  struct _Octant *parent;            /* parent of the octant                 */
  double min[3];                     /* minimum bounds of an octant          */
  double max[3];                     /* max bounds of an octant              */
  int ppid;                          /* parent pid, -1 mean a local root     */
  int id;                            /* octant's id number                   */
  int which;                         /* which child of parent                */
  int numChild;                      /* number of children, 0 == terminal    */
  float cost;                        /* cost of the octant                   */
  int npid;                          /* where to migrate to                  */
  int cpid[8];                       /* the pid of the children              */
  /* int orientation;*/              /* octant traversal orientation         */
} Octant;

typedef struct RL_Node *pRList;         /* typedef for a pointer to rootlist */
typedef struct RL_Node {                /* an entry in the local root list   */
  struct RL_Node *next;                 /* pointer to next node in the list  */
  pOctant oct;                          /* pointer to the root octant        */
} RList;

typedef struct {
  int npid;
  double min[3];
  double max[3];
} Map;

extern pRList OCT_rootlist;          /* list of all the local roots          */
extern int OCT_count;                /* count of all local octants           */
extern int OCT_localpid;             /* the processor id                     */
extern int idcount;                  /* count for id's, help with uniqueness */
extern double gmin[3];
extern double gmax[3];
extern int dimension;
extern int GRAY;
extern int HILBERT;
/* extern Map *array; */

#define vector_set(r,a)     \
        ( (r)[0]=(a)[0], (r)[1]=(a)[1], (r)[2]=(a)[2] )
#define vector_eq(a,b)      \
        ( (a)[0]==(b)[0] && (a)[1]==(b)[1] && (a)[2]==(b)[2] )
#define vector_lt(a,b)      \
        ( (a)[0]<(b)[0] && (a)[1]<(b)[1] && (a)[2]<(b)[2] )
#define vector_gt(a,b)      \
        ( (a)[0]>(b)[0] && (a)[1]>(b)[1] && (a)[2]>(b)[2] )
#define vector_add(r,a,b)     \
        ( (r)[0]=(a)[0]+(b)[0],     \
          (r)[1]=(a)[1]+(b)[1],     \
          (r)[2]=(a)[2]+(b)[2] )
#define vector_set_comp(r,x,y,z) \
        ( (r)[0]=(x), (r)[1]=(y), (r)[2]=(z) )
#define vector_cmult(r,c,a)       \
        ( (r)[0]=c*(a)[0],        \
          (r)[1]=c*(a)[1],        \
          (r)[2]=c*(a)[2] )
#define vector_divc(r,a,c)        \
        ( (r)[0]=(a)[0]/c,        \
          (r)[1]=(a)[1]/c,        \
          (r)[2]=(a)[2]/c )
#define vector_dist(a,b) \
        ( sqrt( ((a)[0]-(b)[0])*((a)[0]-(b)[0]) + \
                ((a)[1]-(b)[1])*((a)[1]-(b)[1]) + \
                ((a)[2]-(b)[2])*((a)[2]-(b)[2]) ) )

extern void    POC_init(int pid, int dim);
extern pOctant POC_new();
extern void    POC_free(pOctant oct);
#ifdef LGG_MIGOCT
extern void    POC_setID(pOctant oct, int id);
#endif /* LGG_MIGOCT */
extern int     POC_id(pOctant oct);
extern void    POC_setparent(pOctant oct, pOctant parent, int ppid);
extern void    POC_setchildnum(pOctant oct, int childnum);
extern int     POC_childnum(pOctant oct);
extern void    POC_setchild(pOctant oct, int i, pOctant child);
extern void    POC_setchildren(pOctant oct, pOctant children[8], int cpids[8]);
extern void    POC_setbounds(pOctant oct, double min[3], double max[3]);
extern void    POC_bounds(pOctant oct, double min[3], double max[3]);
extern pOctant POC_parent(pOctant oct);
extern pOctant POC_child(pOctant oct, int i);
extern int     POC_children(pOctant oct, pOctant children[8]);
extern int     POC_isTerminal(pOctant oct);
extern pRegion POC_regionlist(pOctant oct);
extern void    POC_addRegion(pOctant oct, pRegion region);
extern void    POC_remRegion(pOctant oct, pRegion region);
extern void    POC_clearRegions(pOctant oct);
extern int     POC_nRegions(pOctant oct);
extern pRList  POC_localroots();
extern void    POC_modify_cost(pOctant oct, float cost);
extern void    POC_modify_newpid(pOctant oct, int newpid);
extern float   POC_data_cost(pOctant oct);
extern int     POC_data_newpid(pOctant oct);
extern int     POC_nlocal(pOctant oct);
extern int     POC_nOctants(void);
extern void    POC_origin_volume(pOctant oct,double origin[3],double *volume);
extern void    POC_printResults();
extern pOctant POC_malloc();
extern void    POC_DfsTraversal(pOctant oct);
extern void    POC_printRegionInfo(pOctant oct);
extern pOctant POC_nextDfs(pOctant octant);
extern int     POC_local(pOctant octant, int i);
extern void    POC_setCpid(pOctant octant, int i, int cpid); 
extern int     POC_delTree(pOctant root);
/* extern pOctant POC_nextDfs(pOcant octant);
 * extern void    POC_setOrientation(pOctant octant, int orientation);
 * extern int     POC_getOrientation(pOctant octant);
 */
