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
static char *cvs_gr_hash_c = "$Id$";
#endif

#include "all_const.h"
#include "gr_const.h"
#include "gr_hash_const.h"
#include "gr_list_const.h"
#include "id_util_const.h"
#include "vx_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  Functions implementing the graph as hash tables.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 * Prototypes
 */

static GRAPH_ID_FN find_vertex_hash_graph;
static GRAPH_VERTEX_FN add_vertex_hash_graph;
static GRAPH_VERTEX_FN delete_vertex_hash_graph;
static GRAPH_FREE_FN free_hash_graph;
static GRAPH_ITERATOR_FN first_vertex_hash_graph;
static GRAPH_ITERATOR_FN next_vertex_hash_graph;

static HASH_TABLE *new_hash_table();
static void free_hash_table(HASH_TABLE **);

static HASH_FN hash_function_one;

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_initialize_hash_graph(GRAPH *graph)
{
/*
 *  Function that initializes function pointers for graphs implemented as 
 *  hash tables.
 */

  graph->Find_Vertex     = find_vertex_hash_graph;
  graph->Add_Vertex      = add_vertex_hash_graph;
  graph->Delete_Vertex   = delete_vertex_hash_graph;
  graph->Free_Graph      = free_hash_graph;
  graph->First_Vertex    = first_vertex_hash_graph;
  graph->Next_Vertex     = next_vertex_hash_graph;
  graph->Graph_Data      = (void *) new_hash_table();
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static HASH_TABLE *new_hash_table()
{
/*
 *  Routine to allocate and initialize a new hash table.
 */

HASH_TABLE *table;
int i, num_buckets;
  
  table = (HASH_TABLE *) LB_MALLOC(sizeof(HASH_TABLE));

  num_buckets = table->Num_Buckets = NUM_HASH_BUCKETS;
  for (i = 0; i < num_buckets; i++) 
    table->Bucket[i] = NULL;
  table->Hash_Fn = hash_function_one;

  return(table);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void free_hash_table(HASH_TABLE **table) 
{
int i; 
int num_buckets = (*table)->Num_Buckets;

  for (i = 0; i < num_buckets; i++) 
    LB_free_bucket(&((*table)->Bucket[i]));

  LB_FREE(table);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static VERTEX *find_vertex_hash_graph(GRAPH *graph, ID *id)
{
/*
 *  Function to find a vertex in a graph.  Returns a pointer to the vertex.
 */

int hash_value;
HASH_TABLE *table = (HASH_TABLE *) (graph->Graph_Data);
  
  hash_value = table->Hash_Fn(id);
  return(LB_search_bucket(table->Bucket[hash_value], id));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void add_vertex_hash_graph(GRAPH *graph, VERTEX *vertex)
{
/*
 *  Function to add a vertex to the graph.
 */

int hash_value;
ID *id = &(vertex->Id);
HASH_TABLE *table = (HASH_TABLE *) (graph->Graph_Data);
  
  hash_value = table->Hash_Fn(id);
  LB_add_to_bucket(&(table->Bucket[hash_value]), vertex);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void delete_vertex_hash_graph(GRAPH *graph, VERTEX *vertex)
{
/*
 *  Function to remove a vertex from the graph.
 */

int hash_value;
ID *id = &(vertex->Id);
HASH_TABLE *table = (HASH_TABLE *) (graph->Graph_Data);

  hash_value = table->Hash_Fn(id);
  LB_remove_from_bucket(&(table->Bucket[hash_value]), vertex);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void free_hash_graph(GRAPH **graph)
{
/*
 * Function that frees all storage associated with a graph
 * (including all vertices of the graph).
 */

  free_hash_table((HASH_TABLE **) (&(*graph)->Graph_Data));
  LB_FREE(graph);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static VERTEX *first_vertex_hash_graph(GRAPH *graph, 
                                       LOOP_CONTROL *loop_control)
{
/*
 *  Initializes the loop control variable;
 *  returns the first entry in the hash table for the graph.
 */
HASH_TABLE_LOOP *loop;

  *loop_control = (LOOP_CONTROL) LB_MALLOC(sizeof(HASH_TABLE_LOOP));
  loop = (HASH_TABLE_LOOP *) (*loop_control);

  /*  Initialize LOOP_CONTROL */

  loop->Bucket = -1;
  loop->Entry = NULL;

  /*  Return the first entry in the hash table */

  return(next_vertex_hash_graph(graph, loop_control));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static VERTEX *next_vertex_hash_graph(GRAPH *graph, 
                                      LOOP_CONTROL *loop_control)
{
HASH_TABLE_LOOP *loop = (HASH_TABLE_LOOP *) (*loop_control);
HASH_TABLE *table = (HASH_TABLE *) (graph->Graph_Data);
LIST_ENTRY *next;
VERTEX *vertex;
int i;

  if (loop->Entry != NULL && (next = loop->Entry->Next) != NULL) {
    /*
     *  Return next value in the bucket.
     */
    vertex = next->Vertex;
    loop->Entry = next;
  }
  else {
    /* 
     * Need to search for the next non-empty bucket.
     */

    i = loop->Bucket + 1;
    while (i < table->Num_Buckets && table->Bucket[i] == NULL) i++;

    if (i < table->Num_Buckets) {

      /*
       *  A new non-empty bucket was found.  Return its first entry.
       */
      loop->Entry = table->Bucket[i];
      vertex = loop->Entry->Vertex;
      loop->Bucket = i;
    }
    else {
      /*
       *  No non-empty buckets were found; return NULL, and free LOOP_CONTROL.
       */

      vertex = NULL;
      LB_FREE(loop_control);
    }
  }

  return(vertex);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int hash_function_one(ID *id)
{
  return(id->Number % NUM_HASH_BUCKETS);
}
