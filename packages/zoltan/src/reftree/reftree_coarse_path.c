#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*****************************************************************************/
/* This file contains the functions that define a path, along with in        */
/* and out vertices, through the elements of the coarse grid                 */
/*****************************************************************************/

#include <stdio.h>
#include "zz_const.h"
#include "reftree.h"
#include "hsfc_hilbert_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* Prototypes for functions internal to this file */

static void free_lists(int all_triangles);
static int init_stack(ZZ *zz);
static int push(int elem, int list, ZZ *zz);
static int pop(int list);
static int set_neigh(ZOLTAN_ID_PTR vertices, int *num_vert,
                     int all_triangles, ZZ *zz);
static int add_neigh_pair(int v,int e1,int e2, ZZ *zz);
static int add_to_to_add(int element, ZZ *zz);
static int initial_cycle(ZZ *zz);
static int add_around_vertices(ZZ *zz);
static int add_around_vertex(int vert, ZZ *zz);
static int get_element_to_add();
static int add_to_cycle(int element);
static int add_to_cycle_by_element(int elementD, int indexB1, int indexB2);
static int try_adding_all(int *ierr, ZZ *zz);
static int element_swap(int *ierr, ZZ *zz);
static int element_swap_recur(int element, int *ierr, ZZ *zz);
static int hanging_nodes(int *ierr, ZZ *zz);
static int vertex_swap(int *ierr, ZZ *zz);
static int broken_link(int *ierr, ZZ *zz);
static int sfc_coarse_grid_path(int nobj,int *num_vert, ZOLTAN_ID_PTR vertices,
                          ZOLTAN_ID_PTR in_vertex, ZOLTAN_ID_PTR out_vertex,
                          int *order, ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
                          char *initpath_method, int all_triangles, ZZ *zz);
static int find_inout(int elem, int prev, int prevprev,
                      ZOLTAN_ID_PTR in_vertex, ZOLTAN_ID_PTR out_vertex,
                      ZOLTAN_ID_PTR vertices, int *num_vert, int *first_vert,
                      ZZ *zz);
static double InvSierpinski2d(ZZ *zz, double *coord);
double sier2d(double x, double y, int state, int level, double addf,
              double addp, double peakx, double peaky);
static void sort_index(int n, double ra[], int indx[]);


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* Variables global to this file */

/* num_neigh[i][j] - number of neighbors of i that share j vertices */
static int **num_neigh;
/* neigh[i][j][k] - the kth neighbor that shares j vertices with element i */
static int ***neigh;
/* neigh_dim[i][j] - the dimension (of the 3rd index) of neigh[i][j] */
static int **neigh_dim;
/* shared_vert[i][j][k][l] - the l'th shared vertex of the kth neighbor ... */
static int ****shared_vert;

/* num_elem[i] - number of elements that contain vertex i */
static int *num_elem;
/* element_list[i][j] - the jth element that contains vertex i */
static int **element_list;
/* element_list_dim - the dimension (of the first index) of element_list */
static int element_list_dim;
/* elem_dim[i] - the dimension (of the second index) of element_list[i] */
static int *elem_dim;

/* to_add[j][] - an array used as a stack of elements that share j vertices
                 with an element on the path, hence is a candidate to add.
                 to_add[0] is used as a stack of vertices to examine when
                 using add_around_vertices. */
static int **to_add;
/* to_add_dim[j] - current size of to_add[j][] */
static int *to_add_dim;
/* to_add_ptr[j] - pointer top of stack j */
static int *to_add_ptr;

/* in this file the vertices are integers from 0 to nvert-1, but in the
   input/output arguments they are gids.  vertex_gid maps from int to gid */
ZOLTAN_ID_PTR vertex_gid;

/* variables for each element which define the path */
static int *prev, *next, *in, *out, *onpath;

/* visited marks which elements have already been visited in element_swap */
static int *visited;

static int num_obj, path_length;

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void free_lists(int all_triangles)
{

/*
 * Free the space in the lists
 */

   int i, j, k;

   for (i=0; i<num_obj; i++) {
      for (j=1; j<=MAXVERT; j++) {
         for (k=0; k<neigh_dim[i][j]; k++) {
            ZOLTAN_FREE(&(shared_vert[i][j][k]));
         }
         ZOLTAN_FREE(&(shared_vert[i][j]));
         ZOLTAN_FREE(&(neigh[i][j]));
      }
      ZOLTAN_FREE(&(shared_vert[i]));
      ZOLTAN_FREE(&(neigh[i]));
      ZOLTAN_FREE(&(neigh_dim[i]));
      ZOLTAN_FREE(&(num_neigh[i]));
   }
   for (j=0; j<=MAXVERT; j++) {
      ZOLTAN_FREE(&(to_add[j]));
   }
   ZOLTAN_FREE(&shared_vert);
   ZOLTAN_FREE(&neigh);
   ZOLTAN_FREE(&neigh_dim);
   ZOLTAN_FREE(&num_neigh);
   ZOLTAN_FREE(&onpath);
   ZOLTAN_FREE(&prev);
   ZOLTAN_FREE(&next);
   ZOLTAN_FREE(&in);
   ZOLTAN_FREE(&out);
   ZOLTAN_FREE(&to_add);
   ZOLTAN_FREE(&to_add_dim);
   ZOLTAN_FREE(&to_add_ptr);
   ZOLTAN_FREE(&vertex_gid);

   if (all_triangles) {
      for (j=0; j<element_list_dim; j++) ZOLTAN_FREE(&(element_list[j]));
      ZOLTAN_FREE(&element_list);
      ZOLTAN_FREE(&num_elem);
      ZOLTAN_FREE(&elem_dim);
   }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int init_stack(ZZ *zz)
{

/*
 * allocate memory for the to_add stacks and initialize them to empty
 */

char *yo = "init_stack";
int i, j;

   to_add = (int **) ZOLTAN_MALLOC(sizeof(int *)*(MAXVERT+1));
   to_add_dim = (int *) ZOLTAN_MALLOC(sizeof(int)*(MAXVERT+1));
   to_add_ptr = (int *) ZOLTAN_MALLOC(sizeof(int)*(MAXVERT+1));
   if (to_add == NULL || to_add_dim == NULL || to_add_ptr == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      Zoltan_Multifree(__FILE__, __LINE__, 3, &to_add,
                                              &to_add_dim,
                                              &to_add_ptr);
      return(ZOLTAN_MEMERR);
   }

   for (i=0; i<=MAXVERT; i++) {
      to_add_ptr[i] = -1;
      to_add_dim[i] = 100; /* arbitrary initial size */
      to_add[i] = (int *) ZOLTAN_MALLOC(sizeof(int)*to_add_dim[i]);
      if (to_add[i] == NULL) {
         ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
         for (j=0; j<i; j++) ZOLTAN_FREE(&(to_add[j]));
         Zoltan_Multifree(__FILE__, __LINE__, 3, &to_add,
                                                 &to_add_dim,
                                                 &to_add_ptr);
         return(ZOLTAN_MEMERR);
      }
   }
   return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int push(int elem, int list, ZZ *zz)
{

/*
 * Pushes element elem onto stack number list
 */

char *yo = "push";

/* make sure there's enough memory */

   if (to_add_ptr[list] >= to_add_dim[list]-1) {
      to_add_dim[list] = 2*to_add_dim[list];
      to_add[list] = (int *) ZOLTAN_REALLOC(to_add[list],
                                            sizeof(int *)*to_add_dim[list]);
      if (to_add[list] == NULL) {
         ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
         return(ZOLTAN_MEMERR);
      }
   }

/* push */

   to_add_ptr[list]++;
   to_add[list][to_add_ptr[list]] = elem;

   return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int pop(int list)
{

/*
 * Returns the top of stack number list, or -1 if it is empty
 */

   if (to_add_ptr[list] < 0) {
      return(-1);
   }
   else {
      to_add_ptr[list]--;
      return(to_add[list][to_add_ptr[list]+1]);
   }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int set_neigh(ZOLTAN_ID_PTR vertices, int *num_vert,
                     int all_triangles, ZZ *zz)
{

/* 
 * determine the neighbor relationships and list of shared vertices
 */

char *yo = "set_neigh";
struct Zoltan_Reftree_inthash_node **hashtable;
int **temp_element_list;
int i, j, k, l, nvert, vert, index, element, vert_count, ierr;

/* 
 * first create a list of elements for each vertex
 */

/* Allocate memory for the element lists.  Start by assuming there are more
   elements than vertices, and no more than 8 elements around each vertex,
   and increase memory later if needed */

   element_list = (int **) ZOLTAN_MALLOC(sizeof(int *)*num_obj);
   num_elem = (int *) ZOLTAN_MALLOC(sizeof(int)*num_obj);
   elem_dim = (int *) ZOLTAN_MALLOC(sizeof(int)*num_obj);
   vertex_gid = ZOLTAN_MALLOC_GID_ARRAY(zz,num_obj);
   if (element_list == NULL || num_elem == NULL || elem_dim == NULL ||
       vertex_gid == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      Zoltan_Multifree(__FILE__, __LINE__, 3, &element_list,
                                              &num_elem,
                                              &elem_dim);
      return(ZOLTAN_MEMERR);
   }
   element_list_dim = num_obj;
   for (i=0; i<num_obj; i++) {
      num_elem[i] = 0;
      elem_dim[i] = 8;
      element_list[i] = (int *) ZOLTAN_MALLOC(sizeof(int)*8);
      if (element_list[i] == NULL) {
         ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
         for (j=0; j<i; j++) ZOLTAN_FREE(&(element_list[j]));
         Zoltan_Multifree(__FILE__, __LINE__, 3, &element_list,
                                                 &num_elem,
                                                 &elem_dim);
         return(ZOLTAN_MEMERR);
      }
   }

/* create the hash table for mapping vertex GIDs into array indices */

   hashtable = (struct Zoltan_Reftree_inthash_node **)
               ZOLTAN_MALLOC(sizeof(struct Zoltan_Reftree_inthash_node *)
                             *DEFAULT_HASH_TABLE_SIZE);
   if (hashtable == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      for (j=0; j<i; j++) ZOLTAN_FREE(&(element_list[j]));
      Zoltan_Multifree(__FILE__, __LINE__, 3, &element_list,
                                              &num_elem,
                                              &elem_dim);
      return(ZOLTAN_MEMERR);
   }
   for (i=0; i<DEFAULT_HASH_TABLE_SIZE; i++)
      hashtable[i] = (struct Zoltan_Reftree_inthash_node *)NULL;

/*
 * go through the elements, putting them on the list of each of their
 * vertices.  And while we're going through, see if they're all triangles.
 */

   nvert = 0;
   index = -1;

   for (element=0; element<num_obj; element++) {

/* go through the vertices of this element */

      for (vert_count=0; vert_count<num_vert[element]; vert_count++) {
         index++;
         vert = Zoltan_Reftree_inthash_lookup(zz,hashtable,
                       &(vertices[index*zz->Num_GID]),DEFAULT_HASH_TABLE_SIZE);
         if (vert == -1) {
            Zoltan_Reftree_IntHash_Insert(zz,&(vertices[index*zz->Num_GID]),
                                      nvert,hashtable,DEFAULT_HASH_TABLE_SIZE);
            vert = nvert;
            nvert++;
            ZOLTAN_SET_GID(zz, &(vertex_gid[vert*zz->Num_GID]),
                               &(vertices[index*zz->Num_GID]));
         }

/* see if we need more memory for the element_list's */

         if (nvert >= element_list_dim) {
            element_list_dim *= 2;
            temp_element_list = (int **) ZOLTAN_REALLOC(element_list,
                                                sizeof(int *)*element_list_dim);
            if (temp_element_list != NULL) element_list = temp_element_list;
            num_elem = (int *) ZOLTAN_REALLOC(num_elem,
                                              sizeof(int)*element_list_dim);
            elem_dim = (int *) ZOLTAN_REALLOC(elem_dim,
                                              sizeof(int)*element_list_dim);
            vertex_gid = ZOLTAN_REALLOC_GID_ARRAY(zz,vertex_gid,element_list_dim);
            if (temp_element_list == NULL || num_elem == NULL ||
                elem_dim == NULL || vertex_gid == NULL) {
               ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
               for (j=0; j<element_list_dim/2; j++) ZOLTAN_FREE(&(element_list[j]));
               Zoltan_Multifree(__FILE__, __LINE__, 3, &element_list,
                                                       &num_elem,
                                                       &elem_dim);
               return(ZOLTAN_MEMERR);
            }
            for (i=element_list_dim/2; i<element_list_dim; i++) {
               num_elem[i] = 0;
               elem_dim[i] = 8;
               element_list[i] = (int *) ZOLTAN_MALLOC(sizeof(int)*8);
               if (element_list[i] == NULL) {
                  ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
                  for (j=0; j<i; j++) ZOLTAN_FREE(&(element_list[j]));
                  Zoltan_Multifree(__FILE__, __LINE__, 3, &element_list,
                                                          &num_elem,
                                                          &elem_dim);
                  return(ZOLTAN_MEMERR);
               }
            }
         }
         if (num_elem[vert] > elem_dim[vert]-1) {
            elem_dim[vert] *= 2;
            element_list[vert] = (int *) ZOLTAN_REALLOC(element_list[vert],
                                            sizeof(int)*elem_dim[vert]);
            if (element_list[vert] == NULL) {
               ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
               for (j=0; j<element_list_dim; j++) ZOLTAN_FREE(&(element_list[j]));
               Zoltan_Multifree(__FILE__, __LINE__, 3, &element_list,
                                                       &num_elem,
                                                       &elem_dim);
               return(ZOLTAN_MEMERR);
            }
         }

/* add the element to the list of this vertex */

         element_list[vert][num_elem[vert]] = element;
         num_elem[vert]++;

      } /* next vert */
   } /* next element */

/* done with the hash table */

   Zoltan_Reftree_Clear_IntHash_Table(hashtable,DEFAULT_HASH_TABLE_SIZE);
   ZOLTAN_FREE(&hashtable);

/*
 * initial memory allocations for neighbor and shared vertices lists,
 * and the in and out vertices
 */

   onpath = (int *) ZOLTAN_MALLOC(sizeof(int)*num_obj);
   prev = (int *) ZOLTAN_MALLOC(sizeof(int)*num_obj);
   next = (int *) ZOLTAN_MALLOC(sizeof(int)*num_obj);
   in = (int *) ZOLTAN_MALLOC(sizeof(int)*num_obj);
   out = (int *) ZOLTAN_MALLOC(sizeof(int)*num_obj);
   if (onpath == NULL || prev == NULL || next == NULL || in == NULL || out == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      for (j=0; j<element_list_dim; j++) ZOLTAN_FREE(&(element_list[j]));
      Zoltan_Multifree(__FILE__, __LINE__, 3, &element_list,
                                              &num_elem,
                                              &elem_dim);
      return(ZOLTAN_MEMERR);
   }
   for (i=0; i<num_obj; i++) onpath[i] = FALSE;

   num_neigh = (int **) ZOLTAN_MALLOC(sizeof(int *)*num_obj);
   neigh_dim = (int **) ZOLTAN_MALLOC(sizeof(int *)*num_obj);
   neigh = (int ***) ZOLTAN_MALLOC(sizeof(int **)*num_obj);
   shared_vert = (int ****) ZOLTAN_MALLOC(sizeof(int ***)*num_obj);
   if (num_neigh == NULL || neigh_dim == NULL || neigh == NULL || shared_vert == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      for (j=0; j<element_list_dim; j++) ZOLTAN_FREE(&(element_list[j]));
      Zoltan_Multifree(__FILE__, __LINE__, 3, &element_list,
                                              &num_elem,
                                              &elem_dim);
      return(ZOLTAN_MEMERR);
   }

   for (i=0; i<num_obj; i++) {
      num_neigh[i] = (int *) ZOLTAN_MALLOC(sizeof(int)*(MAXVERT+1));
      neigh_dim[i] = (int *) ZOLTAN_MALLOC(sizeof(int)*(MAXVERT+1));
      neigh[i] = (int **) ZOLTAN_MALLOC(sizeof(int *)*(MAXVERT+1));
      shared_vert[i] = (int ***) ZOLTAN_MALLOC(sizeof(int **)*(MAXVERT+1));
      if (num_neigh[i] == NULL || neigh_dim[i] == NULL || neigh[i] == NULL || shared_vert[i] == NULL) {
         ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
         for (j=0; j<element_list_dim; j++) ZOLTAN_FREE(&(element_list[j]));
         Zoltan_Multifree(__FILE__, __LINE__, 3, &element_list,
                                                 &num_elem,
                                                 &elem_dim);
         return(ZOLTAN_MEMERR);
      }

      for (j=0; j<MAXVERT+1; j++) {
         num_neigh[i][j] = 0;
         if (j != 0) {
            neigh_dim[i][j] = 4;
            neigh[i][j] = (int *) ZOLTAN_MALLOC(sizeof(int)*neigh_dim[i][j]);
            shared_vert[i][j] = (int **) ZOLTAN_MALLOC(sizeof(int *)*neigh_dim[i][j]);
            if (neigh[i][j] == NULL || shared_vert[i][j] == NULL) {
               ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
               for (k=0; k<element_list_dim; k++) ZOLTAN_FREE(&(element_list[k]));
               Zoltan_Multifree(__FILE__, __LINE__, 3, &element_list,
                                                       &num_elem,
                                                       &elem_dim);
               return(ZOLTAN_MEMERR);
            }
         }
         else {
            neigh_dim[i][j] = 0;
            neigh[i][j] = NULL;
            shared_vert[i][j] = NULL;
         }

         for (k=0; k<neigh_dim[i][j]; k++) {
            shared_vert[i][j][k] = (int *) ZOLTAN_MALLOC(sizeof(int)*j);
            if (shared_vert[i][j][k] == NULL) {
               ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
               for (l=0; l<element_list_dim; l++) ZOLTAN_FREE(&(element_list[l]));
               Zoltan_Multifree(__FILE__, __LINE__, 3, &element_list,
                                                       &num_elem,
                                                       &elem_dim);
               return(ZOLTAN_MEMERR);
            }
         }
      }
   }

/*
 * now go through the vertices noting the neighbor relationship of all
 * pairs of elements that contain that vertex
 */

   for (i = 0; i < nvert; i++) {
      for (j = 0; j < num_elem[i]-1; j++) {
         for (k=j+1; k < num_elem[i]; k++) {
            ierr = add_neigh_pair(i,element_list[i][j],element_list[i][k],zz);
            if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
               ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                                  "Error returned from add_neigh_pair.");
               return(ierr);
            }
         }
      }
   }

/* free memory for element lists, unless all elements are triangles */

   if (!all_triangles) {
      for (j=0; j<element_list_dim; j++) ZOLTAN_FREE(&(element_list[j]));
      Zoltan_Multifree(__FILE__, __LINE__, 3, &element_list,
                                              &num_elem,
                                              &elem_dim);
   }

   return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int add_neigh_pair(int v,int e1,int e2, ZZ *zz)
{

/* 
 * adds the elements e1 and e2, which share vertex v, to the neigh
 * and shared_vert lists
 */

char *yo = "add_neigh_pair";
int nshare, index, i, j, k;

/* 
 * see if they are already listed as neighbors
 */

   nshare = 0;
   for (j=1; j<=MAXVERT && nshare==0; j++) {
      for (k=0; k<num_neigh[e1][j] && nshare==0; k++) {
         if (neigh[e1][j][k] == e2) {
            nshare = j;
            index = k;
         }
      }
   }

/*
 * add e2 to the list for e1 with nshare+1 shared vertices
 */

/* see if we need more memory */

   if (neigh_dim[e1][nshare+1] <= num_neigh[e1][nshare+1]) {
      neigh_dim[e1][nshare+1] *= 2;
      neigh[e1][nshare+1] = (int *) ZOLTAN_REALLOC(neigh[e1][nshare+1],
                                           sizeof(int)*neigh_dim[e1][nshare+1]);
      shared_vert[e1][nshare+1] = (int **) ZOLTAN_REALLOC(shared_vert[e1][nshare+1],
                                      sizeof(int *)*neigh_dim[e1][nshare+1]);
      if (neigh[e1][nshare+1] == NULL || shared_vert[e1][nshare+1] == NULL) {
         ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
         return(ZOLTAN_MEMERR);
      }
      for (k=neigh_dim[e1][nshare+1]/2; k<neigh_dim[e1][nshare+1]; k++) {
         shared_vert[e1][nshare+1][k] = (int *) ZOLTAN_MALLOC(sizeof(int)*(nshare+1));
         if (shared_vert[e1][nshare+1][k] == NULL) {
            ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
            return(ZOLTAN_MEMERR);
         }
      }
   }

/* add e2 to the list */

   neigh[e1][nshare+1][num_neigh[e1][nshare+1]] = e2;

/* copy the shared vertices from the previous list (if nshare!=0) and add v */

   for (i=0; i<nshare; i++) shared_vert[e1][nshare+1][num_neigh[e1][nshare+1]][i] = shared_vert[e1][nshare][index][i];
   shared_vert[e1][nshare+1][num_neigh[e1][nshare+1]][nshare] = v;
   num_neigh[e1][nshare+1]++;

/* remove e2 from the previous list by moving the end of list to e2's spot */

   if (nshare != 0) {
      num_neigh[e1][nshare]--;
      neigh[e1][nshare][index] = neigh[e1][nshare][num_neigh[e1][nshare]];
      for (i=0; i<nshare; i++) shared_vert[e1][nshare][index][i] = shared_vert[e1][nshare][num_neigh[e1][nshare]][i];
   }

/* 
 * add e1 to the list for e2
 */

/* find e1 in e2's list, if they previously shared vertices */

   if (nshare != 0) {
      index = -1;
      for (k=0; k<num_neigh[e2][nshare]; k++) {
         if (neigh[e2][nshare][k] == e1) index = k;
      }
      if (index == -1) {
         ZOLTAN_PRINT_ERROR(zz->Proc, yo, "e1 missing from e2's neighbor list");
         return(ZOLTAN_FATAL);
      }
   }

/* see if we need more memory */

   if (neigh_dim[e2][nshare+1] <= num_neigh[e2][nshare+1]) {
      neigh_dim[e2][nshare+1] *= 2;
      neigh[e2][nshare+1] = (int *) ZOLTAN_REALLOC(neigh[e2][nshare+1],
                                           sizeof(int)*neigh_dim[e2][nshare+1]);
      shared_vert[e2][nshare+1] = (int **) ZOLTAN_REALLOC(shared_vert[e2][nshare+1],
                                      sizeof(int *)*neigh_dim[e2][nshare+1]);
      if (neigh[e2][nshare+1] == NULL || shared_vert[e2][nshare+1] == NULL) {
         ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
         return(ZOLTAN_MEMERR);
      }
      for (k=neigh_dim[e2][nshare+1]/2; k<neigh_dim[e2][nshare+1]; k++) {
         shared_vert[e2][nshare+1][k] = (int *) ZOLTAN_MALLOC(sizeof(int)*(nshare+1));
         if (shared_vert[e2][nshare+1][k] == NULL) {
            ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
            return(ZOLTAN_MEMERR);
         }
      }
   }

/* add e1 to the list */

   neigh[e2][nshare+1][num_neigh[e2][nshare+1]] = e1;

/* copy the shared vertices from the previous list (if nshare!=0) and add v */

   for (i=0; i<nshare; i++) shared_vert[e2][nshare+1][num_neigh[e2][nshare+1]][i] = shared_vert[e2][nshare][index][i];
   shared_vert[e2][nshare+1][num_neigh[e2][nshare+1]][nshare] = v;
   num_neigh[e2][nshare+1]++;

/* remove e1 from the previous list by moving the end of list to e1's spot */

   if (nshare != 0) {
      num_neigh[e2][nshare]--;
      neigh[e2][nshare][index] = neigh[e2][nshare][num_neigh[e2][nshare]];
      for (i=0; i<nshare; i++) shared_vert[e2][nshare][index][i] = shared_vert[e2][nshare][num_neigh[e2][nshare]][i];
   }

   return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int add_to_to_add(int element, ZZ *zz)
{

/* 
 * add all the neighbors of element that are not on the cycle to the
 * lists of elements that can be added to the cycle
 */

int num_share, nay, ierr;

   for (num_share=1; num_share<=MAXVERT; num_share++) {
      for (nay=0; nay<num_neigh[element][num_share]; nay++) {
         if (!onpath[neigh[element][num_share][nay]]) {
            ierr = push(neigh[element][num_share][nay],num_share,zz);
            if (ierr == ZOLTAN_FATAL || ierr == ZOLTAN_MEMERR) return(ierr);
         }
      }
   }
   return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int initial_cycle(ZZ *zz)
{

/*
 * create the initial cycle of two elements
 */

char *yo = "initial_cycle";
int elementA, elementB, j, verta, vertb, ierr;

/*
 * find two elements that share at least two vertices
 */

   elementA = 0;
   elementB = -1;
   while (elementB == -1) {
      j = MAXVERT;
      while (j > 1 && num_neigh[elementA][j] == 0) j--;
      if (j==1) {
         elementA++;
         if (elementA >= num_obj) {
            ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                               "couldn't find two elements that share a side");
            return(ZOLTAN_FATAL);
         }
      }
      else {
         elementB = neigh[elementA][j][0];
      }
   }
   verta = shared_vert[elementA][j][0][0];
   vertb = shared_vert[elementA][j][0][1];

/*
 * create the cycle aAbBa
 */

   prev[elementA] = elementB;
   in  [elementA] = verta;
   onpath[elementA] = TRUE;
   out [elementA] = vertb;
   next[elementA] = elementB;
   prev[elementB] = elementA;
   in  [elementB] = vertb;
   onpath[elementB] = TRUE;
   out [elementB] = verta;
   next[elementB] = elementA;
   path_length = 2;

/* put all the neighbors onto a to_add list */

   ierr = add_to_to_add(elementA,zz);
   if (ierr == ZOLTAN_FATAL || ierr == ZOLTAN_MEMERR) return(ierr);
   ierr = add_to_to_add(elementB,zz);
   if (ierr == ZOLTAN_FATAL || ierr == ZOLTAN_MEMERR) return(ierr);

/* put the in and out vertices onto the vertex to_add (to add around) list */

   ierr = push(verta, 0, zz);
   if (ierr == ZOLTAN_FATAL || ierr == ZOLTAN_MEMERR) return(ierr);
   ierr = push(vertb, 0, zz);
   if (ierr == ZOLTAN_FATAL || ierr == ZOLTAN_MEMERR) return(ierr);

   return(ZOLTAN_OK);

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int add_around_vertices(ZZ *zz)
{

/*
 * This function adds elements to the path by adding several elements around
 * vertices on the path at once.  See add_around_vertex for details.
 */

   int vert, element, i, ierr;

/* repeat until there are no more vertices to consider */

   vert = pop(0);
   while (vert != -1) {
      ierr = add_around_vertex(vert,zz);
      if (ierr == ZOLTAN_FATAL || ierr == ZOLTAN_MEMERR) return(ierr);
      vert = pop(0);
   }

/* if that didn't get all the elements, put elements that are not on the path
   but are adjacent to an element on the path onto the to_add lists */

   if (path_length != num_obj) {
      element = 0;
      i = 0;
      while ((element != 0 && i <= num_obj) || i == 0) {
         ierr = add_to_to_add(element,zz);
         if (ierr == ZOLTAN_FATAL || ierr == ZOLTAN_MEMERR) return(ierr);
         element = next[element];
         i++;
      }
   }

   return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int add_around_vertex(int vert, ZZ *zz)
{

/*
 * This function adds a contiguous set of elements around vertex vert to
 * the path.  vert must be on the path, there must be more than one
 * element of vert that is not on the path, and they must share a vertex
 * other than vert.  If the path currently contains A vert B and there are
 * elements C, D, ..., E which share vertices c, d, ... which are not vert,
 * then the new path is A vert C c D d ... E vert B
 *
 * This comes from the proof that a path exists for triangles.  I don't know
 * what it will do with other element shapes, especially in 3D where the
 * notion of going around a vertex is not well defined.
 */

int elementA, elementB, elementC, elementD, elementE;
int vertexc, vertexd;
int i, j, k, l, num_share, pass, success, ierr;

/* find an element with vert as the out vertex */

  elementA = -1;
  for (i=0; i<num_elem[vert] && elementA == -1; i++) {
    if (onpath[element_list[vert][i]] && out[element_list[vert][i]] == vert) {
       elementA = element_list[vert][i];
    }
  }

  if (elementA == -1) { /* vert is not on the path */
    return(ZOLTAN_OK);
  }

  elementB = next[elementA];

/* look at each element of this vertex that is not on the path, until
   some elements are added */

  success = FALSE;

  for (pass=0; pass<2 && !success; pass++) {
   for (i=0; i<num_elem[vert] && !success; i++) {
    elementC = element_list[vert][i];
    if (!onpath[elementC]) {

/* on the first pass, require that C share at least 2 vertices with A */
/* TEMP if this ever gets used for 3D, want C and A to share a face */

      if (pass == 0) {
        for (num_share=2; num_share<=MAXVERT && !success; num_share++) {
          for (j=0; j<num_neigh[elementA][num_share] && !success; j++) {
            if (neigh[elementA][num_share][j] == elementC) success = TRUE;
          }
        }
      }
      else {
        success = TRUE;
      }

      if (success) {
       success = FALSE;

/* look at the elements around vert to see if there is one that shares at
   least two vertices with C and is not on the path, and identify a
   shared vertex that is not vert. */

       for (j=0; j<num_elem[vert] && !success; j++) {
        elementD = element_list[vert][j];
        if (!onpath[elementD] && elementD != elementC) {
          for (num_share=2; num_share <= MAXVERT && !success; num_share++) {
            for (k=0; k<num_neigh[elementC][num_share] && !success; k++) {
              if (neigh[elementC][num_share][k] == elementD) {
                success = TRUE;
                if (shared_vert[elementC][num_share][k][0] == vert) {
                  vertexc = shared_vert[elementC][num_share][k][1];
                }
                else {
                  vertexc = shared_vert[elementC][num_share][k][0];
                }
              }
            }
          }
        }
      }

/* if we found one, then add C and D to the path, add vert back to the
   vertices to be considered (in case we didn't get all the elements
   around it), and add vertexc to the vertices to be considered. */

       if (success) {
        next[elementA] = elementC;
        prev[elementC] = elementA;
        in  [elementC] = vert;
        onpath[elementC] = TRUE;
        out [elementC] = vertexc;
        next[elementC] = elementD;
        prev[elementD] = elementC;
        in  [elementD] = vertexc;
        onpath[elementD] = TRUE;
        out [elementD] = vert;
        next[elementD] = elementB;
        prev[elementB] = elementD;
        path_length += 2;
        ierr = push(vert,0,zz);
        if (ierr == ZOLTAN_FATAL || ierr == ZOLTAN_MEMERR) return(ierr);
        ierr = push(vertexc,0,zz);
        if (ierr == ZOLTAN_FATAL || ierr == ZOLTAN_MEMERR) return(ierr);
       }
     }
    }
   }
  }

/* look for additional elements that can be added */

  while (success) {

/* look at the elements around vert to see if there is one that shares at
   least two vertices with D, one of which is not vert or vertexc,  and is
   not on the path, and identify a shared vertex that is not vert or vertexc */

    success = FALSE;
    for (j=0; j<num_elem[vert] && !success; j++) {
      elementE = element_list[vert][j];
      if (!onpath[elementE]) {
        for (num_share=2; num_share<=MAXVERT && !success; num_share++) {
          for (k=0; k<num_neigh[elementD][num_share] && !success; k++) {
            if (neigh[elementD][num_share][k] == elementE) {
              for (l=0; l<num_share && !success; l++) {
                vertexd = shared_vert[elementD][num_share][k][l];
                if (vertexd != vert && vertexd != vertexc) success = TRUE;
              }
            }
          }
        }
      }
    }

/* if we found one, add it to the path and assign it as elementD to look more,
   and add vertexd to the vertices to be considered. */

    if (success) {

      out [elementD] = vertexd;
      next[elementD] = elementE;
      prev[elementE] = elementD;
      in  [elementE] = vertexd;
      onpath[elementE] = TRUE;
      out [elementE] = vert;
      next[elementE] = elementB;
      prev[elementB] = elementE;
      path_length += 1;

      vertexc = vertexd;
      elementD = elementE;

      ierr = push(vertexd,0,zz);
      if (ierr == ZOLTAN_FATAL || ierr == ZOLTAN_MEMERR) return(ierr);
    }
  }
  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int get_element_to_add()
{

/*
 * returns pop(to_add[j]) for largest j where to_add is not empty
 */

int j, element;

   j = MAXVERT;
   element = -1;
   while (j > 0 && element == -1) {
      element = pop(j);
      j--;
   }
   return(element);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int add_to_cycle(int element)
{

/*
 * Adds element to the cycle, if it goes in next to one of its neighbors
 */

int success, num_share, nay;

/*
 * can't add it if it is already on the cycle
 */

   if (onpath[element]) return(FALSE);

/*
 * go through the neighbors until element is added to the cycle or no
 * more neighbors
 */

   success = FALSE;
   num_share = MAXVERT;
   while (num_share > 0 && !success) {
      nay = 0;
      while (nay < num_neigh[element][num_share] && !success) {
         if (onpath[neigh[element][num_share][nay]]) {
            success = add_to_cycle_by_element(element,num_share,nay);
         }
         nay++;
      }
      num_share--;
   }
   return(success);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int add_to_cycle_by_element(int elementD, int indexB1, int indexB2)
{

/*
 * Adds elementD to the cycle next to elementB (given by
 * neigh[elementD][indexB1][indexB2]) if it meets one of the
 * criteria for easily slipping it in.
 */

int elementA, elementB, elementC;
int verta, vertb, vertc, vertd, verte, vertf;
/* indexA1, indexB1 and indexC1 also tell how many vertices are shared with D */
int indexA1, indexA2, indexC1, indexC2;
int j, k, ivert1, ivert2, D_has_b_and_c, success;

/*
 * The nearby segment of the cycle is aAbBcCd where lower case
 * letters are vertices and upper case letters are elements
 */

   elementB = neigh[elementD][indexB1][indexB2];
   elementA = prev[elementB];
   elementC = next[elementB];
   verta = in [elementA];
   vertb = in [elementB];
   vertc = out[elementB];
   vertd = out[elementC];

/* see if A and C share vertices with D, and get indices if they do */

   indexA1 = 0; indexA2 = 0; indexC1 = 0; indexC2 = 0;
   for (j=1; j<=MAXVERT; j++) {
      for (k=0; k<num_neigh[elementD][j]; k++) {
         if (neigh[elementD][j][k] == elementA) {
            indexA1 = j;
            indexA2 = k;
         }
         if (neigh[elementD][j][k] == elementC) {
            indexC1 = j;
            indexC2 = k;
         }
      }
   }

   success = FALSE;

/*
 * If D contains b and shares another vertex e with B, e!=c
 * then new cycle is aAbDeBcCd
 */

   for (ivert1=0; ivert1<indexB1 && !success; ivert1++) {
      if (shared_vert[elementD][indexB1][indexB2][ivert1] == vertb) {
         for (ivert2=0; ivert2<indexB1 && !success; ivert2++) {
            verte = shared_vert[elementD][indexB1][indexB2][ivert2];
            if (!success && verte != vertb && verte != vertc) {
               next[elementA] = elementD;
               prev[elementD] = elementA;
               in  [elementD] = vertb;
               onpath[elementD] = TRUE;
               out [elementD] = verte;
               next[elementD] = elementB;
               prev[elementB] = elementD;
               in  [elementB] = verte;
               path_length++;
               success = TRUE;
            }
         }
      }
   }

/*
 * If D contains c and shares another e with B, e!=b
 * then new cycle is aAbBeDcCd
 */

   for (ivert1=0; ivert1<indexB1 && !success; ivert1++) {
      if (shared_vert[elementD][indexB1][indexB2][ivert1] == vertc) {
         for (ivert2=0; ivert2<indexB1 && !success; ivert2++) {
            verte = shared_vert[elementD][indexB1][indexB2][ivert2];
            if (!success && verte != vertb && verte != vertc) {
               out [elementB] = verte;
               next[elementB] = elementD;
               prev[elementD] = elementB;
               in  [elementD] = verte;
               onpath[elementD] = TRUE;
               out [elementD] = vertc;
               next[elementD] = elementC;
               prev[elementC] = elementD;
               path_length++;
               success = TRUE;
            }
         }
      }
   }

/*
 * If D contains both b and c, and A and B share e with e!=a and e!= b
 * then new cycle is aAeBbDcCd
 */

   D_has_b_and_c = FALSE;
   for (ivert1=0; ivert1<indexB1 && !success && !D_has_b_and_c; ivert1++) {
      if (shared_vert[elementD][indexB1][indexB2][ivert1] == vertb) {
         for (ivert2=0; ivert2<indexB1 && !D_has_b_and_c; ivert2++) {
            if (shared_vert[elementD][indexB1][indexB2][ivert1] == vertc) {
               D_has_b_and_c = TRUE;
            }
         }
      }
   }

   if (D_has_b_and_c) {
      for (j=MAXVERT; j>1 && !success; j--) {
         for (k=0; k<num_neigh[elementB][j] && !success; k++) {
            if (neigh[elementB][j][k] == elementA) {
               for (ivert1=0; ivert1<j && !success; ivert1++) {
                  verte = shared_vert[elementB][j][k][ivert1];
                  if (!success && verte != verta && verte != vertb) {
                     out [elementA] = verte;
                     in  [elementB] = verte;
                     out [elementB] = vertb;
                     next[elementB] = elementD;
                     prev[elementD] = elementB;
                     in  [elementD] = vertb;
                     onpath[elementD] = TRUE;
                     out [elementD] = vertc;
                     next[elementD] = elementC;
                     path_length++;
                     success = TRUE;
                  }
               }
            }
         }
      }
   }

/*
 * If D contains both b and c, and A and C share e with e!=d and e!= c
 * then new cycle is aAbDcBeCd
 */

   if (D_has_b_and_c) {
      for (j=MAXVERT; j>1 && !success; j--) {
         for (k=0; k<num_neigh[elementB][j] && !success; k++) {
            if (neigh[elementB][j][k] == elementC) {
               for (ivert1=0; ivert1<j && !success; ivert1++) {
                  verte = shared_vert[elementB][j][k][ivert1];
                  if (!success && verte != vertd && verte != vertc) {
                     prev[elementD] = elementA;
                     in  [elementD] = vertb;
                     onpath[elementD] = TRUE;
                     out [elementD] = vertc;
                     next[elementD] = elementB;
                     prev[elementB] = elementD;
                     in  [elementB] = vertc;
                     out [elementB] = verte;
                     in  [elementC] = verte;
                     path_length++;
                     success = TRUE;
                  }
               }
            }
         }
      }
   }

/*
 * If D shares an e != c with B and shares f with A with f!=e and f!=a
 * then new cycle is aAfDeBcCd
 */

   for (ivert1=0; ivert1<indexB1 && !success; ivert1++) {
      verte = shared_vert[elementD][indexB1][indexB2][ivert1];
      if (verte != vertc) {
         for (ivert2=0; ivert2<indexA1 && !success; ivert2++) {
            vertf = shared_vert[elementD][indexA1][indexA2][ivert2];
            if (!success && vertf != verte && vertf != verta) {
               out [elementA] = vertf;
               next[elementA] = elementD;
               prev[elementD] = elementA;
               in  [elementD] = vertf;
               onpath[elementD] = TRUE;
               out [elementD] = verte;
               next[elementD] = elementB;
               prev[elementB] = elementD;
               in  [elementB] = verte;
               path_length++;
               success = TRUE;
            }
         }
      }
   }

/*
 * If D shares an e != b with B and shares f with C with f!=e and f!=d
 * then new cycle is aAbBeDfCd
 */

   for (ivert1=0; ivert1<indexB1 && !success; ivert1++) {
      verte = shared_vert[elementD][indexB1][indexB2][ivert1];
      if (verte != vertb) {
         for (ivert2=0; ivert2<indexC1 && !success; ivert2++) {
            vertf = shared_vert[elementD][indexC1][indexC2][ivert2];
            if (!success && vertf != verte && vertf != vertd) {
               out [elementB] = verte;
               next[elementB] = elementD;
               prev[elementD] = elementB;
               in  [elementD] = verte;
               onpath[elementD] = TRUE;
               out [elementD] = vertf;
               next[elementD] = elementC;
               prev[elementC] = elementD;
               in  [elementC] = vertf;
               path_length++;
               success = TRUE;
            }
         }
      }
   }

   return(success);

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int try_adding_all(int *ierr, ZZ *zz)
{

/*
 * Try adding every element that is not already on the path, until we find one
 */

int success, element;

   success = FALSE;
   for (element=0; element < num_obj && !success; element++) {
      if (!onpath[element]) {
         success = add_to_cycle(element);
         if (success) {
           *ierr = add_to_to_add(element,zz);
           if (*ierr == ZOLTAN_FATAL || *ierr == ZOLTAN_MEMERR) return(success);
         }
      }
   }
   return(success);
}

static int element_swap(int *ierr, ZZ *zz)
{

/* When an element shares both the in and out vertices of an element on the
   path, this element can replace the other element on the path, and it may
   be possible to add that element back into the path.  This function attempts
   all possibilities of such swapping of which element is on the path.  For
   tetrahedra this should never be necessary unless the grid contains a local
   cut edge.  For triangles, this procedure is guaranteed to succeed if the
   grid contains at least one interior vertex unless the grid contains a
   local cut point.  For other element shapes, it might work.  Note we only
   have to consider neighbors that share exactly 2 vertices because with only
   1 shared vertex it can't have both in and out, and with more than two it
   would be already added to the path under condition 1 in add_to_cycle. */

   char *yo = "element_swap";
   int i, success, element;

   visited = (int *) ZOLTAN_MALLOC(sizeof(int)*num_obj);
   if (visited == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      *ierr = ZOLTAN_MEMERR;
      return(FALSE);
   }

   for (i=0; i<num_obj; i++) visited[i] = FALSE;

   success = FALSE;
   for (element=0; element<num_obj && !success; element++) {
      if (!onpath[element]) {
         visited[element] = TRUE;
         success = element_swap_recur(element,ierr,zz);
         if (*ierr == ZOLTAN_FATAL || *ierr == ZOLTAN_MEMERR) {
            ZOLTAN_FREE(&visited);
            return(success);
         }
         visited[element] = FALSE;
      } /* endif !onpath */
   } /* end for each element */
   element--;

   if (success && !(*ierr == ZOLTAN_FATAL || *ierr == ZOLTAN_MEMERR)) {
      *ierr = add_to_to_add(element,zz);
   }

   ZOLTAN_FREE(&visited);
   return(success);
}

static int element_swap_recur(int element, int *ierr, ZZ *zz)
{

/* Recursively swap elements that share in and out vertices.  See the comments
   in element_swap.  */

/* arbitrarily limit the number of swapable neighbors for simplicity. */
/* 18 is the number of neighbors in a uniform hexahedron grid */
#define MAXSWAP 18

   int swapable[MAXSWAP];
   int i, k, num_swapable, swap, success;

   success = FALSE;

/* create a list of all neighbors that:
   1) are on the path
   2) share exactly two vertices
   3) those vertices are the in and out vertices of that neighbor
   4) have not already been visited in this branch of the recursion
*/

   num_swapable = 0;
   for (k=0; k<num_neigh[element][2] && num_swapable<MAXSWAP; k++) {
      if (onpath[neigh[element][2][k]] && !visited[neigh[element][2][k]]) {
         if ((shared_vert[element][2][k][1] == in [neigh[element][2][k]] ||
              shared_vert[element][2][k][1] == out[neigh[element][2][k]]) &&
             (shared_vert[element][2][k][2] == in [neigh[element][2][k]] ||
              shared_vert[element][2][k][2] == out[neigh[element][2][k]])) {
            swapable[num_swapable] = neigh[element][2][k];
            num_swapable++;
         }
      }
   }

/* for each swapable neighbor, swap it and try to add it back or
   recursively try it's swapable neighbors */

   for (i=0; i<num_swapable && !success; i++) {
      swap = swapable[i];
      visited[swap] = TRUE;

/* swap that element off the path */

      prev[element] = prev[swap];
      in  [element] = in  [swap];
      onpath[element] = TRUE;
      out [element] = out [swap];
      next[element] = next[swap];
      onpath[swap] = FALSE;

/* try to add it to the path */

      success = add_to_cycle(swap);

/* if it can't be added, recursively try it's neighbors */

      if (!success) {
         success = element_swap_recur(swap,ierr,zz);
         if (*ierr == ZOLTAN_FATAL || *ierr == ZOLTAN_MEMERR) return(success);
      }

/* if that didn't work, swap it back on the path and move on to the next one */

      if (!success) {
         prev[swap] = prev[element];
         in  [swap] = in  [element];
         onpath[swap] = TRUE;
         out [swap] = out [element];
         next[swap] = next[element];
         onpath[element] = FALSE;
      }

/* if something was added, put the swapped element's neighbors on the
   to_add list */

      else {
         *ierr = add_to_to_add(swap,zz);
      }

      visited[swap] = FALSE;

   } /* end for each swapable neighbor */

   return(success);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int hanging_nodes(int *ierr,ZZ *zz)
{

/*
 * This function tries two path insertions which deal with hanging nodes
 * where one element is on the path and two are not, and the vertices
 * are in the right configuration.  (Some other hanging nodes with two
 * elements on the path and one not on the path are handled by case 6
 * in insert_by_element, i.e. aAbBeDfCd).  These cases are not included
 * in insert_by_element for two reasons: 1) we don't want to apply this
 * case to triangles and tetrahedra, and 2) this adds two elements to
 * the path and we need to apply add_to_to_add to both of them.
 */

int success;
int j, k, j2, k2, j3, k3;
int elementB, elementC, elementD;
int vertb, vertc, verte, vertf, vert1, vert2, vert3;

/*
 * go through all the elements, D, until one is added
 */

 success = FALSE;
 for (elementD=0; elementD < num_obj && !success; elementD++) {
  if (!onpath[elementD]) {

/* look for a neighbor, B, of D that is on the path */

   for (j=MAXVERT; j>0 && !success; j--) {
    for (k=0; k<num_neigh[elementD][j] && !success; k++) {
     elementB = neigh[elementD][j][k];
     if (onpath[elementB]) {

/* the path currently contains bBc */

      vertb = in [elementB];
      vertc = out[elementB];

/*
 * if D contains vertex b ...
 */

      for (vert1=0; vert1<j && !success; vert1++) {
       if (shared_vert[elementD][j][k][vert1] == vertb) {

/* look for an element, C, not on the path, which shares a vertex, e, with
   D, e!=b and a vertex, f, with B, f!=e and f!=c */

        for (j2=MAXVERT; j2>0 && !success; j2--) {
         for (k2=0; k2<num_neigh[elementD][j2] && !success; k2++) {
          elementC = neigh[elementD][j2][k2];
          if (!onpath[elementC]) {
           for (vert2=0; vert2<j2 && !success; vert2++) {
            verte = shared_vert[elementD][j2][k2][vert2];
            if (verte != vertb) {
             for (j3=MAXVERT; j3>0 && !success; j3--) {
              for (k3=0; k3<num_neigh[elementC][j3] && !success; k3++) {
               if (neigh[elementC][j3][k3] == elementB) {
                for (vert3=0; vert3<j3 && !success; vert3++) {
                 vertf = shared_vert[elementC][j3][k3][vert3];
                 if (vertf != verte && vertf != vertc) {

/* found one, the new path is bDeCfBc */

                  next[prev[elementB]] = elementD;
                  prev[elementD] = prev[elementB];
                  in  [elementD] = vertb;
                  onpath[elementD] = TRUE;
                  out [elementD] = verte;
                  next[elementD] = elementC;
                  prev[elementC] = elementD;
                  in  [elementC] = verte;
                  onpath[elementC] = TRUE;
                  out [elementC] = vertf;
                  next[elementC] = elementB;
                  prev[elementB] = elementC;
                  in  [elementB] = vertf;
                  path_length += 2;
                  success = TRUE;
                  *ierr = add_to_to_add(elementD,zz);
                  if (*ierr == ZOLTAN_FATAL || *ierr == ZOLTAN_MEMERR)
                     return(success);
                  *ierr = add_to_to_add(elementC,zz);
                  if (*ierr == ZOLTAN_FATAL || *ierr == ZOLTAN_MEMERR)
                     return(success);
                 }
                }
               }
              }
             }
            }
           }
          }
         }
        }
       }
      }

/*
 * if that didn't work and D contains vertex c ...
 */

      for (vert1=0; vert1<j && !success; vert1++) {
       if (shared_vert[elementD][j][k][vert1] == vertc) {

/* look for an element, C, not on the path, which shares a vertex, e, with
   D, e!=c and a vertex, f, with B, f!=e and f!=b */

        for (j2=MAXVERT; j2>0 && !success; j2--) {
         for (k2=0; k2<num_neigh[elementD][j2] && !success; k2++) {
          elementC = neigh[elementD][j2][k2];
          if (!onpath[elementC]) {
           for (vert2=0; vert2<j2 && !success; vert2++) {
            verte = shared_vert[elementD][j2][k2][vert2];
            if (verte != vertc) {
             for (j3=MAXVERT; j3>0 && !success; j3--) {
              for (k3=0; k3<num_neigh[elementC][j3] && !success; k3++) {
               if (neigh[elementC][j3][k3] == elementB) {
                for (vert3=0; vert3<j3 && !success; vert3++) {
                 vertf = shared_vert[elementC][j3][k3][vert3];
                 if (vertf != verte && vertf != vertb) {

/* found one, the new path is bBfCeDc */

                  prev[next[elementB]] = elementD;
                  next[elementD] = next[elementB];
                  out [elementD] = vertc;
                  onpath[elementD] = TRUE;
                  in  [elementD] = verte;
                  prev[elementD] = elementC;
                  next[elementC] = elementD;
                  out [elementC] = verte;
                  onpath[elementC] = TRUE;
                  in  [elementC] = vertf;
                  prev[elementC] = elementB;
                  next[elementB] = elementC;
                  out [elementB] = vertf;
                  path_length += 2;
                  success = TRUE;
                  *ierr = add_to_to_add(elementD,zz);
                  if (*ierr == ZOLTAN_FATAL || *ierr == ZOLTAN_MEMERR)
                     return(success);
                  *ierr = add_to_to_add(elementC,zz);
                  if (*ierr == ZOLTAN_FATAL || *ierr == ZOLTAN_MEMERR)
                     return(success);
                 }
                }
               }
              }
             }
            }
           }
          }
         }
        }
       }
      }
     }
    }
   }
  }
 }
 return(success);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int vertex_swap(int *ierr, ZZ *zz)
{

/*
 * look for a vertex on the path that can be replaced by a different vertex
 * between the same two elements such that the change enables a neighbor
 * element to be added to the path
 */

int success;
int element, neighbor, count, j, k, j2, k2, found_prev, vert1, vert2, v;

/*
 * start with the first element that is on the cycle and follow the cycle until
 * a swap is made or we have gone all the way around
 */

   success = FALSE;
   element = 0;
   while (!onpath[element] && element<num_obj) element++;
   if (!onpath[element]) return(success);

   for (count=0; count<path_length && !success; count++) {

/*
 * See if this element it has a neighbor not on the path.
 * Need not consider neighbors with only one shared vertex.
 */

      for (j=MAXVERT; j>1 && !success; j--) {
         for (k=0; k<num_neigh[element][j] && !success; k++) {
            neighbor = neigh[element][j][k];
            if (!onpath[neighbor]) {

/*
 * Find the previous element on the path in element's neighbors
 * (don't have to look at it if they only share one vertex)
 */

               found_prev = FALSE;
               for (j2=MAXVERT; j2>1 && !found_prev; j2--) {
                  for (k2=0; k2<num_neigh[element][j2] && !found_prev; k2++) {
                     if (neigh[element][j2][k2] == prev[element]) {
                        found_prev = TRUE;

/*
 * either all the shared vertices between element and neighbor are among the
 * in and out of element, or none of them are.  See which case this is.
 */

                        if (shared_vert[element][j][k][0] ==  in[element] ||
                            shared_vert[element][j][k][0] == out[element]) {

/*
 * They are ins and outs.
 *
 * Look for a vertex shared by element and prev[element] which is not
 * in[element], in[prev[element] or out[element], and make that the
 * new in/out between them if you find one.
 */

                           for (vert1=0; vert1<j2 && !success; vert1++) {
                              v = shared_vert[element][j2][k2][vert1];
                              if (v != in[element] && v != out[element] &&
                                  v != in[prev[element]]) {
                                 in[element] = v;
                                 out[prev[element]] = v;
                                 success = TRUE;
                                 *ierr = add_to_to_add(element,zz);
                                 if (*ierr == ZOLTAN_FATAL ||
                                     *ierr == ZOLTAN_MEMERR) return(success);
                                 *ierr = add_to_to_add(prev[element],zz);
                                 if (*ierr == ZOLTAN_FATAL ||
                                     *ierr == ZOLTAN_MEMERR) return(success);
                              }
                           }
                        }
                        else {

/*
 * The shared vertices are not the in and out vertices of element.
 *
 * Look for a vertex shared by all three that is not out[element] or
 * in[prev[element]], and make that the in/out if found
 */

                           for (vert1=0; vert1<j && !success; vert1++) {
                              v = shared_vert[element][j][k][vert1];
                              if (v != out[element] && v != in[prev[element]]) {
                                 for (vert2=0; vert2<j2 && !success; vert2++) {
                                    if (shared_vert[element][j2][k2][vert2] == v) {
                                       in[element] = v;
                                       out[prev[element]] = v;
                                       success = TRUE;
                                       *ierr = add_to_to_add(element,zz);
                                       if (*ierr == ZOLTAN_FATAL ||
                                           *ierr == ZOLTAN_MEMERR) return(success);
                                       *ierr = add_to_to_add(prev[element],zz);
                                       if (*ierr == ZOLTAN_FATAL ||
                                           *ierr == ZOLTAN_MEMERR) return(success);
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }

      element = next[element];
   }

   return(success);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int broken_link(int *ierr, ZZ *zz)
{

/*
 * If all else fails, then put in a broken link so that we can at least
 * compute partitions, although they might not be connected
 */

int success;
int elementB, elementC, elementD;
int vertb, verte, vertf, vert1;
int j, k;
int success1, j1, k1, l1;

/*
 * Find an element, D, not on the path that is next to an element, B,
 * on the path, bBcCd.  Let D share vertex e!=b with B and let f be any
 * other vertex of D.  The new path is bBeDfcCd.  Note the broken link
 * where f!=c.
 */

/* go through all the elements, D, until one is added */

 success = FALSE;
 for (elementD=0; elementD < num_obj && !success; elementD++) {
  if (!onpath[elementD]) {

/* look for a neighbor, B, of D that is on the path */

   for (j=MAXVERT; j>0 && !success; j--) {
    for (k=0; k<num_neigh[elementD][j] && !success; k++) {
     elementB = neigh[elementD][j][k];
     if (onpath[elementB]) {

/* the path currently contains bBcC */

      vertb = in [elementB];
      /* vertc = out[elementB];  KDD Commented out to remove set-but-never-used
                                     compiler warning */
      elementC = next[elementB];

/* find a vertex, e, in D and B that is not b */

      for (vert1=0; vert1<j && !success; vert1++) {
       verte = shared_vert[elementD][j][k][vert1];
       if (verte != vertb) {

/* and a vertex, f, in D that is not e */

        success1 = FALSE;
        for (j1=MAXVERT; j1>0 && !success1; j1--) {
         for (k1=0; k1<num_neigh[elementD][j1] && !success1; k1++) {
          for (l1=0; l1<j1 && !success1; l1++) {
           vertf = shared_vert[elementD][j1][k1][l1];
           if (vertf != verte) {
            success1 = TRUE;
           }
          }
         }
        }

/* and set the path to bBeDfcCd */

        if (success1) {
         out [elementB] = verte;
         next[elementB] = elementD;
         prev[elementD] = elementB;
         in  [elementD] = verte;
         onpath[elementD] = TRUE;
         out [elementD] = vertf;
         next[elementD] = elementC;
         prev[elementC] = elementD;
         path_length++;
         success = TRUE;
         *ierr = add_to_to_add(elementD,zz);
         if (*ierr == ZOLTAN_FATAL || *ierr == ZOLTAN_MEMERR) return(success);
         *ierr = add_to_to_add(elementB,zz);
         if (*ierr == ZOLTAN_FATAL || *ierr == ZOLTAN_MEMERR) return(success);
        }
       }
      }
     }
    }
   }
  }
 }
 return(success);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int sfc_coarse_grid_path(int nobj, int *num_vert, ZOLTAN_ID_PTR vertices,
                         ZOLTAN_ID_PTR in_vertex, ZOLTAN_ID_PTR out_vertex,
                         int *order, ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
                         char *initgrid_method, int all_triangles, ZZ *zz)
{

/*
 * This routine uses either the Hilbert or Sierpinski Space Filling Curve to
 * find a path through the coarse grid.  This path is not necessarily connected
 * in non-convex domains.  Even if it is connected, the path through refinements
 * of the * grid may not be connected if the in/out vertices can't be assigned
 * to match up.
 * This uses the Hilbert routines to map from 2D and 3D to 1D in
 * hsfc/hsfc_hilbert.c.  It has its own Sierpinski routine, which only
 * applies to 2D problems.
 */

  char *yo = "sfc_coarse_grid_path";
  int ierr, num_geom, i, elem, prev, prevprev;
  int *ind, *first_vert;
  double *sfccoord, *coords;
  double xmin,xmax,ymin,ymax,zmin,zmax;

  ZOLTAN_TRACE_ENTER(zz, yo);
  ierr = ZOLTAN_OK;

/*
 * Find out if this is a 2D or 3D problem
 */

  if (zz->Get_Num_Geom == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Must register ZOLTAN_NUM_GEOM_FN.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ZOLTAN_FATAL);
  }
  num_geom = zz->Get_Num_Geom(zz->Get_Num_Geom_Data, &ierr);
  if (ierr) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                   "Error returned from user function Get_Num_Geom.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ierr); 
  }
  if (!(num_geom==2 || num_geom==3)) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Geometry must be either 2D or 3D.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ZOLTAN_FATAL);
  }

/*
 * For each element, determine its SFC mapping from its coordinates to 1D
 */

/* verify the geometry function has be registered */

  if (zz->Get_Geom == NULL && zz->Get_Geom_Multi == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                      "Must register ZOLTAN_GEOM_FN or ZOLTAN_GEOM_MULTI_FN.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ZOLTAN_FATAL);
  }

/* allocate space for the element coordinates and the SFC coordinate */

  coords = (double *) ZOLTAN_MALLOC(nobj*num_geom*sizeof(double));
  sfccoord = (double *) ZOLTAN_MALLOC(nobj*sizeof(double));
  if (sfccoord == NULL || coords == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ZOLTAN_MEMERR);
  }

/* get the coordinates of each element and find max and min */

  xmin =  1.0e50;
  xmax = -1.0e50;
  ymin =  1.0e50;
  ymax = -1.0e50;
  zmin =  1.0e50;
  zmax = -1.0e50;
  if (zz->Get_Geom) {
    for (i=0; i<nobj; i++) {

      zz->Get_Geom(zz->Get_Geom_Data, zz->Num_GID, zz->Num_LID, 
                   &(gids[zz->Num_GID*i]), &(lids[zz->Num_GID*i]),
                   &(coords[num_geom*i]),&ierr);
      if (ierr) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                           "Error returned from user function Get_Geom.");
        Zoltan_Multifree(__FILE__,__LINE__, 2, &coords, &sfccoord);
        ZOLTAN_TRACE_EXIT(zz, yo);
        return(ierr);
      }
    }
  }
  else {
    /* Use MULTI Function */
    zz->Get_Geom_Multi(zz->Get_Geom_Multi_Data, zz->Num_GID, zz->Num_LID, nobj,
                       gids, lids, num_geom, coords, &ierr);
    if (ierr) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                         "Error returned from user function Get_Geom_Multi.");
      Zoltan_Multifree(__FILE__,__LINE__, 2, &coords, &sfccoord);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return(ierr);
    }
  }

  for (i = 0; i < nobj; i++) {
    if (coords[num_geom*i  ] < xmin) xmin = coords[num_geom*i  ];
    if (coords[num_geom*i  ] > xmax) xmax = coords[num_geom*i  ];
    if (coords[num_geom*i+1] < ymin) ymin = coords[num_geom*i+1];
    if (coords[num_geom*i+1] > ymax) ymax = coords[num_geom*i+1];
    if (num_geom == 3) {
      if (coords[num_geom*i+1] < zmin) zmin = coords[num_geom*i+2];
      if (coords[num_geom*i+2] > zmax) zmax = coords[num_geom*i+2];
    }
  }

/* scale the coordinates to [0,1] and compute the SFC mapping */

  for (i=0; i<nobj; i++) {
    coords[num_geom*i  ] = (coords[num_geom*i  ]-xmin)/(xmax-xmin);
    coords[num_geom*i+1] = (coords[num_geom*i+1]-ymin)/(ymax-ymin);
    if ((strcmp(initgrid_method, "SIERPINSKI") == 0) ||
        (strcmp(initgrid_method, "REFTREE_DEFAULT") == 0 && all_triangles)) {
      if (num_geom == 2) {
        sfccoord[i] = InvSierpinski2d(zz,&(coords[num_geom*i]));
      } else {
        ierr = ZOLTAN_FATAL;
        ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                           "Sierpinski SFC only applies to 2D problems.");
        Zoltan_Multifree(__FILE__,__LINE__, 2, &coords, &sfccoord);
        ZOLTAN_TRACE_EXIT(zz, yo);
        return(ierr);
      }
    } else { /* Hilbert */
      if (num_geom == 2) {
        sfccoord[i] = Zoltan_HSFC_InvHilbert2d(zz,&(coords[num_geom*i]));
      } else {
        coords[num_geom*i+2] = (coords[num_geom*i+2]-zmin)/(zmax-zmin);
        sfccoord[i] = Zoltan_HSFC_InvHilbert3d(zz,&(coords[num_geom*i]));
      }
    }
  }

/*
 * sort the Hilbert coordinates to get the order of the elements
 */

  ind = (int *) ZOLTAN_MALLOC(nobj*sizeof(int));
  if (ind == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    Zoltan_Multifree(__FILE__,__LINE__, 2, &coords, &sfccoord);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ZOLTAN_MEMERR);
  }

  for (i=0; i<nobj; i++) ind[i] = i;
  sort_index(nobj, sfccoord, ind);

/*
 * Determine where each element's vertices start in vertices
 */

  first_vert = (int *) ZOLTAN_MALLOC(nobj*sizeof(int));
  if (first_vert == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    Zoltan_Multifree(__FILE__,__LINE__, 3, &coords, &sfccoord, &ind);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ZOLTAN_MEMERR);
  }

  first_vert[0] = 0;
  for (i=1; i<nobj; i++) first_vert[i] = first_vert[i-1] + num_vert[i-1];

/* 
 * pass through the elements in order setting order and looking for
 * in/out vertices
 */

  
  elem = ind[0];
  order[elem] = 0;
  ZOLTAN_SET_GID(zz,&(in_vertex[zz->Num_GID*elem]),
                    &(vertices[zz->Num_GID*first_vert[elem]]));

  for (i=1; i<nobj; i++) {
    elem = ind[i];
    order[elem] = i;
    prev = ind[i-1];
    if (i==1) {
      prevprev = -1;
    } else {
      prevprev = ind[i-2];
    }

    ierr = find_inout(elem, prev, prevprev, in_vertex, out_vertex,
                      vertices, num_vert, first_vert, zz);
    if (ierr) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                     "Error returned from find_inout.");
      Zoltan_Multifree(__FILE__,__LINE__, 4, &coords,&sfccoord,&ind,&first_vert);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return(ierr); 
    }
  }

  if (ZOLTAN_EQ_GID(zz,&(in_vertex[zz->Num_GID*elem]),
                       &( vertices[zz->Num_GID*first_vert[elem]]))) {
    ZOLTAN_SET_GID(zz,&(out_vertex[zz->Num_GID*elem]),
                      &(  vertices[zz->Num_GID*(first_vert[elem]+1)]));
  } else {
    ZOLTAN_SET_GID(zz,&(out_vertex[zz->Num_GID*elem]),
                      &(  vertices[zz->Num_GID*first_vert[elem]]));
  }

  Zoltan_Multifree(__FILE__,__LINE__, 4, &coords, &sfccoord, &ind, &first_vert);
  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int find_inout(int elem, int prev, int prevprev,
               ZOLTAN_ID_PTR in_vertex, ZOLTAN_ID_PTR out_vertex,
               ZOLTAN_ID_PTR vertices, int *num_vert, int *first_vert,
               ZZ *zz)
{

/*
 * This routine finds and sets an in-vertex for elem and out-vertex for prev.
 * Usually they will be the same vertex, but if necessary they will be
 * different which will cause the Hamiltonian path to be disconnected.
 * If necessary and possible, this routine will change the in-vertex of prev
 * and out-vertex of prevprev to make in and out be the same.
 */

  char *yo = "find_inout";
  int looking, i, j, ngid, shared;

  ZOLTAN_TRACE_ENTER(zz, yo);
  ngid = zz->Num_GID;

/*
 * first look for a shared vertex of elem and prev that is not the
 * in-vertex of prev
 */

  looking = TRUE;
  for (i=0; i<num_vert[elem] && looking; i++) {
    for (j=0; j<num_vert[prev] && looking; j++) {
      if (ZOLTAN_EQ_GID(zz,&(vertices[ngid*(first_vert[prev]+j)]),
                           &(vertices[ngid*(first_vert[elem]+i)])) &&
          !ZOLTAN_EQ_GID(zz,&( vertices[ngid*(first_vert[prev]+j)]),
                            &(in_vertex[ngid*prev]))) {
        ZOLTAN_SET_GID(zz,&(in_vertex[ngid*elem]),
                          &( vertices[ngid*(first_vert[elem]+i)]));
        ZOLTAN_SET_GID(zz,&(out_vertex[ngid*prev]),
                          &(vertices[ngid*(first_vert[elem]+i)]));
        looking = FALSE;
      }
    }
  }

/*
 * If that failed then the in-vertex of prev must be shared, or else there
 * are no shared vertices between elem and prev.
 */

  if (looking) {
    shared = FALSE;
    for (i=0; i<num_vert[elem] && !shared; i++) {
      if (ZOLTAN_EQ_GID(zz,&( vertices[ngid*(first_vert[elem]+i)]),
                           &(in_vertex[ngid*prev]))) {
        shared = TRUE;
      }
    }

/*
 * If the in-vertex of prev is not shared, then prev and elem are not adjacent,
 * so pick any vertices for the in and out.
 */

    if (!shared) {
      ZOLTAN_SET_GID(zz,&(in_vertex[ngid*elem]),
                        &( vertices[ngid*(first_vert[elem])]));
      if (ZOLTAN_EQ_GID(zz,&(in_vertex[ngid*prev]),
                           &( vertices[ngid*first_vert[prev]]))) {
        ZOLTAN_SET_GID(zz,&(out_vertex[ngid*prev]),
                          &(  vertices[ngid*(first_vert[prev]+1)]));
      } else {
        ZOLTAN_SET_GID(zz,&(out_vertex[ngid*prev]),
                          &(  vertices[ngid*first_vert[prev]]));
      }
      looking = FALSE;
    }
  }

/*
 * If the in-vertex of prev is shared, then it must be the only shared vertex
 * between elem and prev (otherwise the first approach would have worked).
 * Try to change the in-vertex of prev and use it as the out-vertex of prev
 * and in-vertex of elem.
 */

/*
 * First, if prev is the beginning of the path then set in-vertex to be
 * any other vertex.
 */

  if (looking) {
    if (prevprev == -1) {
      ZOLTAN_SET_GID(zz,&( in_vertex[ngid*elem]),
                        &( in_vertex[ngid*prev]));
      ZOLTAN_SET_GID(zz,&(out_vertex[ngid*prev]),
                        &( in_vertex[ngid*prev]));
      if (ZOLTAN_EQ_GID(zz,&(in_vertex[ngid*prev]),
                           &( vertices[ngid*first_vert[prev]]))) {
        ZOLTAN_SET_GID(zz,&(in_vertex[ngid*prev]),
                          &( vertices[ngid*(first_vert[prev]+1)]));
      } else {
        ZOLTAN_SET_GID(zz,&(in_vertex[ngid*prev]),
                          &( vertices[ngid*first_vert[prev]]));
      }
      looking = FALSE;
    }
  }

/*
 * Second, if prev and prevprev are not adjacent or contain a broken link,
 * i.e., out_vertex of prevprev is not the in_vertex of prev, then change
 * the in_vertex of prev to any other vertex.
 */

  if (looking) {
    if (!ZOLTAN_EQ_GID(zz,&( in_vertex[ngid*prev]),
                          &(out_vertex[ngid*prevprev]))) {
      ZOLTAN_SET_GID(zz,&( in_vertex[ngid*elem]),
                        &( in_vertex[ngid*prev]));
      ZOLTAN_SET_GID(zz,&(out_vertex[ngid*prev]),
                        &( in_vertex[ngid*prev]));
      if (ZOLTAN_EQ_GID(zz,&(in_vertex[ngid*prev]),
                           &( vertices[ngid*first_vert[prev]]))) {
        ZOLTAN_SET_GID(zz,&(in_vertex[ngid*prev]),
                          &( vertices[ngid*(first_vert[prev]+1)]));
      } else {
        ZOLTAN_SET_GID(zz,&(in_vertex[ngid*prev]),
                          &( vertices[ngid*first_vert[prev]]));
      }
      looking = FALSE;
    }
  }

/*
 * Third and final, look for a shared vertex between prev and prevprev that
 * is not the in-vertex of either prev or prevprev and use that as the new
 * out-vertex of prevprev and in-vertex of prev.
 */

  if (looking) {
    for (i=0; i<num_vert[prev] && looking; i++) {
      if (!ZOLTAN_EQ_GID(zz,&(in_vertex[ngid*prev]),
                            &( vertices[ngid*(first_vert[prev]+i)]))) {
        for (j=0; j<num_vert[prevprev] && looking; j++) {
          if (ZOLTAN_EQ_GID(zz,&(vertices[ngid*(first_vert[prevprev]+j)]),
                               &(vertices[ngid*(first_vert[prev]+i)])) &&
              !ZOLTAN_EQ_GID(zz,&( vertices[ngid*(first_vert[prevprev]+j)]),
                                &(in_vertex[ngid*prevprev]))) {
            ZOLTAN_SET_GID(zz,&( in_vertex[ngid*elem]),
                              &( in_vertex[ngid*prev]));
            ZOLTAN_SET_GID(zz,&(out_vertex[ngid*prev]),
                              &( in_vertex[ngid*prev]));
            ZOLTAN_SET_GID(zz,&(in_vertex[ngid*prev]),
                              &( vertices[ngid*(first_vert[prev]+i)]));
            ZOLTAN_SET_GID(zz,&(out_vertex[ngid*prevprev]),
                              &(  vertices[ngid*(first_vert[prev]+i)]));
            looking = FALSE;
          }
        }
      }
    }
  }

/*
 * If that failed, give up and put in a broken link by using the in-vertex
 * of prev (the only shared vertex between elem and prev) as the in-vertex
 * of elem and any other vertex of prev as the out-vertex of prev
 */

  if (looking) {
    ZOLTAN_SET_GID(zz,&( in_vertex[ngid*elem]),
                      &( in_vertex[ngid*prev]));
    if (ZOLTAN_EQ_GID(zz,&(in_vertex[ngid*prev]),
                         &( vertices[ngid*first_vert[prev]]))) {
      ZOLTAN_SET_GID(zz,&(out_vertex[ngid*prev]),
                        &(  vertices[ngid*(first_vert[prev]+1)]));
    } else {
      ZOLTAN_SET_GID(zz,&(out_vertex[ngid*prev]),
                        &(  vertices[ngid*first_vert[prev]]));
    }
  }

  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double InvSierpinski2d(ZZ *zz, double *coord)
{

/*
 * Given x,y coordinates in [0,1]x[0,1], returns the Sierpinski key [0,1]
 */

  char *yo = "InvSierpinski2d";

   /* sanity check for input arguments */
  if ((coord[0] < 0.0) || (coord[0] > 1.0) || (coord[1] < 0.0) ||
      (coord[1] > 1.0))
     ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Spatial Coordinates out of range.");

/* begin recursion that computes key */

  return sier2d(coord[0], coord[1], 0, 0, (double) 0.5, (double) 0.5,
                (double) 0.0, (double) 0.0);

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double sier2d(double x, double y, int state, int level, double addf,
              double addp, double peakx, double peaky)
{

static int MAXLEV = 24;

/*
 * Recursively compute the Sierpinski mapping from the unit square to the
 * unit interval
 */

  if (level >= MAXLEV) return (double) 0;

  switch (state) {

  case 0:
    if (y > x)
      return sier2d(x,y,1,level+1,addf/2,addp,(double) 0.0, (double) 1.0);
    else
      return sier2d(x,y,2,level+1,addf/2,addp,(double) 1.0, (double) 0.0)+addf;

  case 1:
    if ( y-peaky < -(x-peakx) )
      return sier2d(x,y,5,level+1,addf/2,addp,peakx+addp,peaky-addp);
    else
      return sier2d(x,y,6,level+1,addf/2,addp,peakx+addp,peaky-addp) + addf;

  case 2:
    if ( y-peaky >= -(x-peakx) )
      return sier2d(x,y,7,level+1,addf/2,addp,peakx-addp,peaky+addp);
    else
      return sier2d(x,y,8,level+1,addf/2,addp,peakx-addp,peaky+addp) + addf;

  case 3:
    if ( y-peaky <= x-peakx )
      return sier2d(x,y,8,level+1,addf/2,addp,peakx+addp,peaky+addp);
    else
      return sier2d(x,y,5,level+1,addf/2,addp,peakx+addp,peaky+addp) + addf;

  case 4:
    if ( y-peaky > x-peakx )
      return sier2d(x,y,6,level+1,addf/2,addp,peakx-addp,peaky-addp);
    else
      return sier2d(x,y,7,level+1,addf/2,addp,peakx-addp,peaky-addp) + addf;

  case 5:
    if ( y < peaky )
      return sier2d(x,y,1,level+1,addf/2,addp/2,peakx-addp,peaky);
    else
      return sier2d(x,y,3,level+1,addf/2,addp/2,peakx-addp,peaky) + addf;

  case 6:
    if ( x < peakx )
      return sier2d(x,y,4,level+1,addf/2,addp/2,peakx,peaky+addp);
    else
      return sier2d(x,y,1,level+1,addf/2,addp/2,peakx,peaky+addp) + addf;

  case 7:
    if ( y >= peaky )
      return sier2d(x,y,2,level+1,addf/2,addp/2,peakx+addp,peaky);
    else
      return sier2d(x,y,4,level+1,addf/2,addp/2,peakx+addp,peaky) + addf;

  case 8:
    if ( x >= peakx )
      return sier2d(x,y,3,level+1,addf/2,addp/2,peakx,peaky-addp);
    else
      return sier2d(x,y,2,level+1,addf/2,addp/2,peakx,peaky-addp) + addf;

  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void sort_index(int n, double ra[], int indx[])
  
/*
*       Numerical Recipies in C source code
*       modified to have first argument an integer array
*       WFM 8/4/2004 I copied this routine from driver/dr_util.c and
*       modified it to take ra as a double.
*
*       Sorts the array ra[0,..,(n-1)] in ascending numerical order using
*       heapsort algorithm.
*       Array ra is not reorganized.  An index array indx is built that
*       gives the new order.
*
*/

{
  int   l, j, ir, i;
  int   irra;
  double rra;

  /*
   *  No need to sort if one or fewer items.
   */
  if (n <= 1) return;

  l=n >> 1;
  ir=n-1;
  for (;;) {
    if (l > 0) {
      rra=ra[indx[--l]];
      irra = indx[l];
    }
    else {
      rra=ra[indx[ir]];
      irra=indx[ir];

      indx[ir]=indx[0];
      if (--ir == 0) {
        indx[0]=irra;
        return;
      }
    }
    i=l; 
    j=(l << 1)+1;
    while (j <= ir) {
      if (j < ir &&
          (ra[indx[j]] <  ra[indx[j+1]]))
        ++j;
      if (rra <  ra[indx[j]]) {
        indx[i] = indx[j];
        j += (i=j)+1;
      }
      else j=ir+1;
    }
    indx[i]=irra;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_Reftree_Coarse_Grid_Path(int nobj, int *num_vert, 
                               ZOLTAN_ID_PTR vertices, ZOLTAN_ID_PTR in_vertex,
                               ZOLTAN_ID_PTR out_vertex, int *order,
                               ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
                               char *initpath_method, ZZ *zz)
{

/*
 * Find a cycle through all the elements of the coarse grid which passes
 * from element to element through vertices, not going in and out of an
 * element through the same vertex.
 */

char *yo = "Zoltan_Reftree_Coarse_Grid_Path";
int all_triangles = 0, i, element, success, ierr;

   ZOLTAN_TRACE_ENTER(zz, yo);
   num_obj = nobj;
   ierr = ZOLTAN_OK;

/* check for a non-triangle */

   all_triangles = TRUE;
   for (element=0; element<num_obj; element++) {
      if (num_vert[element] != 3) all_triangles = FALSE;
   }

/* see if an space filling curve method was requested */

   if (strcmp(initpath_method, "HILBERT") == 0 ||
       strcmp(initpath_method, "SIERPINSKI") == 0 ||
       strcmp(initpath_method, "REFTREE_DEFAULT") == 0 ) {
      ierr = sfc_coarse_grid_path(nobj, num_vert, vertices, in_vertex,
                                  out_vertex, order, gids, lids, 
                                  initpath_method, all_triangles, zz);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ierr);
   }

/* otherwise, use the method for connected partitions */

/* initialize the stacks of elements that can be added to the path */

   ierr = init_stack(zz);
   if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                         "Error returned from init_stack.");
      free_lists(all_triangles);
      return(ierr);
   }

/* determine the neighbor relationships and list of shared vertices */

   ierr = set_neigh(vertices, num_vert, all_triangles, zz);
   if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                         "Error returned from set_neigh.");
      free_lists(all_triangles);
      return(ierr);
   }

/* create the initial cycle */

   ierr = initial_cycle(zz);
   if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                         "Error returned from initial_cycle.");
      free_lists(all_triangles);
      return(ierr);
   }

/* if the grid is all triangles, add elements around vertices first */

   if (all_triangles) {
      ierr = add_around_vertices(zz);
      if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
         ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                            "Error returned from add_around_vertices.");
         free_lists(all_triangles);
         return(ierr);
      }
   }

/* repeat until all elements are on the path */

   while (path_length < num_obj) {

/* while there are elements on the to_add list, get one, add it if you can,
   and put its neighbors that are not on the path onto the to_add list */

      while ( (element = get_element_to_add() ) != -1 && path_length < num_obj) {
         success = add_to_cycle(element);
         if (success) {
            ierr = add_to_to_add(element,zz);
            if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
               ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                                  "Error returned from add_to_to_add.");
               free_lists(all_triangles);
               return(ierr);
            }
         }
      }

/* if all elements are on the path, then done */

      success = (path_length == num_obj);

/* if not done, then first try adding all elements that are not on the
   cycle to make sure I didn't miss any */

      if (!success) {
         success = try_adding_all(&ierr,zz);
         if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
            ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                               "Error returned from try_adding_all.");
            free_lists(all_triangles);
            return(ierr);
         }
      }

/* if the elements are triangles, then swap which elements are on the
   path to start it up again */

      if (!success && all_triangles) {
         success = element_swap(&ierr, zz);
         if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
            ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                               "Error returned from element_swap.");
            free_lists(all_triangles);
            return(ierr);
         }
      }

/* look for a vertex in the path that can be replaced to start it up again */

      if (!success) {
         success = vertex_swap(&ierr, zz);
         if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
            ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                               "Error returned from vertex_swap.");
            free_lists(all_triangles);
            return(ierr);
         }
      }

/* look for hanging node type relationships that can get it started again */

      if (!success) {
         success = hanging_nodes(&ierr,zz);
         if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
            ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                               "Error returned from hanging_nodes.");
            free_lists(all_triangles);
            return(ierr);
         }
      }

/* for non-triangles, try swapping which elements are on the path */

      if (!success && !all_triangles) {
         success = element_swap(&ierr, zz);
         if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
            ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                               "Error returned from element_swap.");
            free_lists(all_triangles);
            return(ierr);
         }
      }


/* give up and put in a broken link */

      if (!success) {
         success = broken_link(&ierr, zz);
         if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
            ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                               "Error returned from broken_link.");
            free_lists(all_triangles);
            return(ierr);
         }
      }

/* broken_link should always succeed, so if success is still false, something
   went wrong */

      if (!success) {
         ZOLTAN_PRINT_ERROR(zz->Proc, yo, "failed to insert broken link");
         ZOLTAN_TRACE_EXIT(zz, yo);
         free_lists(all_triangles);
         return(ZOLTAN_FATAL);
      }

   }

/* set the return values */

   element = 0;
   i = 0;
   while ((element != 0 && i <= num_obj) || i == 0) {
      order[element] = i;
      ZOLTAN_SET_GID(zz,&(in_vertex[element*zz->Num_GID]),
                        &(vertex_gid[in[element]*zz->Num_GID]));
      ZOLTAN_SET_GID(zz,&(out_vertex[element*zz->Num_GID]),
                        &(vertex_gid[out[element]*zz->Num_GID]));
      element = next[element];
      i++;
   }

   free_lists(all_triangles);
   ZOLTAN_TRACE_EXIT(zz, yo);
   return(ZOLTAN_OK);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
