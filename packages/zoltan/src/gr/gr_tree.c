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
static char *cvs_gr_tree_c = "$Id$";
#endif

#include "all_const.h"
#include "gr_const.h"
#include "gr_tree_const.h"
#include "id_util_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  Functions implementing the graph as AVL trees.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 * Prototypes
 */

static GRAPH_ID_FN find_vertex_tree_graph;
static GRAPH_VERTEX_FN add_vertex_tree_graph;
static GRAPH_VERTEX_FN delete_vertex_tree_graph;
static GRAPH_ITERATOR_FN first_vertex_tree_graph;
static GRAPH_ITERATOR_FN next_vertex_tree_graph;
static GRAPH_FREE_FN free_tree_graph;
static VERTEX *find_entry(TREE *, ID *, ID_INT_FN *);
static void free_tree(TREE **p_root);

/*
 *  Prototypes for static routines that manipulate the AVL trees:
 */

static void single_LL_rotation(TREE **, TREE *, TREE *);
static void double_LR_rotation(TREE **, TREE *, TREE *);
static void single_RR_rotation(TREE **, TREE *, TREE *);
static void double_RL_rotation(TREE **, TREE *, TREE *);
static BOOLEAN avl_insert(TREE **, TREE *, VERTEX *, ID_INT_FN *fn, 
                          BOOLEAN *);
static void avl_del(TREE **, BOOLEAN *);
static void *avl_delete(TREE **, VERTEX *, ID_INT_FN *fn, BOOLEAN *);
static void avl_balance1(TREE **, BOOLEAN *);
static void avl_balance2(TREE **, BOOLEAN *);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_initialize_tree_graph(GRAPH *graph)
{
/*
 *  Function that initializes function pointers for graphs implemented as 
 *  AVL trees.
 */

  graph->Find_Vertex     = find_vertex_tree_graph;
  graph->Add_Vertex      = add_vertex_tree_graph;
  graph->Delete_Vertex   = delete_vertex_tree_graph;
  graph->First_Vertex    = first_vertex_tree_graph;
  graph->Next_Vertex     = next_vertex_tree_graph;
  graph->Free_Graph      = free_tree_graph;
  graph->Graph_Data      = NULL;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static VERTEX *find_vertex_tree_graph(GRAPH *graph, ID *id)
{
/*
 *  Function to find a vertex in a graph.  Returns a pointer to the vertex.
 */

  return(find_entry((TREE *) (graph->Graph_Data), id, BL_ID_Util.Compare));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static VERTEX *find_entry(TREE *p_root, ID *id, ID_INT_FN *fn)
{
int comparison;
VERTEX *vertex;

/*
 *  Recursively search through the AVL tree for the vertex with Id == id.
 */

  if (p_root) {
    comparison = fn(id, &(p_root->Vertex->Id));
    if (comparison < 0)
      vertex = find_entry(p_root->Left, id, fn);
    else if (comparison > 0)
      vertex = find_entry(p_root->Right, id, fn);
    else
      vertex = p_root->Vertex;
  }

  return(vertex);

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void add_vertex_tree_graph(GRAPH *graph, VERTEX *vertex)
{
/*
 *  Function to add a vertex to the graph.
 */

BOOLEAN h = FALSE;
TREE *p_root = (TREE *) graph->Graph_Data;
TREE *p_parent = (p_root ? p_root->Parent : NULL);
void print_tree_graph(TREE *p_root);

  avl_insert((TREE **) (&(graph->Graph_Data)), p_parent, vertex, 
             BL_ID_Util.Compare, &h);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void delete_vertex_tree_graph(GRAPH *graph, VERTEX *vertex)
{
/*
 *  Function to remove a vertex from the graph.
 */

BOOLEAN h = FALSE;

  avl_delete((TREE **) (&(graph->Graph_Data)), vertex, BL_ID_Util.Compare, &h);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void free_tree_graph(GRAPH **graph)
{
/* 
 * Function that frees all storage associated with a graph
 * (including all vertices of the graph).
 */

  free_tree((TREE **) (&((*graph)->Graph_Data)));
  LB_FREE(graph);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void free_tree(TREE **p_root)
{
  if (*p_root) {
    free_tree(&((*p_root)->Left));
    free_tree(&((*p_root)->Right));
    LB_free_vertex(&((*p_root)->Vertex));
    LB_FREE(p_root);
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static VERTEX *first_vertex_tree_graph(GRAPH *graph, 
                                       LOOP_CONTROL *loop_control)
{
/*
 *  Initializes the loop control variable;
 *  returns the first entry in the AVL tree for the graph.
 */
TREE_LOOP *loop;

  /*
   * Allocate space for loop control structure at pointer loop_control.
   * Assign loop to same address.
   */

  *loop_control = (LOOP_CONTROL) LB_MALLOC(sizeof(TREE_LOOP));
  loop = (TREE_LOOP *) (*loop_control);

  /* Initialize LOOP_CONTROL */

  loop->Current_Node = (TREE *) (graph->Graph_Data);
  loop->Go_Left = TRUE;
  loop->Go_Mid = FALSE;
  loop->Go_Right = FALSE;
  loop->Go_Up = FALSE;

  /* Return the first tree node of the tree */

  return(next_vertex_tree_graph(graph, loop_control));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static VERTEX *next_vertex_tree_graph(GRAPH *graph, LOOP_CONTROL *loop_control)
{
/*
 *  Returns the next loop entry in the dynamic data structure.
 */

TREE_LOOP *loop = (TREE_LOOP *) (*loop_control);
TREE *node = loop->Current_Node;

  /*
   * This routine loops through the entries of the tree using an inorder
   * traversal of the tree (i.e., recursively process left subtree first;
   * process root node next; then process right subtree).  The bits of
   * the LOOP_CONTROL "loop" indicate which in direction of the tree
   * to proceed.  For non-NULL trees, values are returned only from the
   * Go_Mid section of the code.
   *
   * If a NULL ptr is returned, the looping is over; the loop control pointer
   * "loop_control" is freed.
   */

  if (node == NULL) {
    /*
     *  Function has been called with a NULL tree.  Return NULL.
     */
    LB_FREE(loop_control);
    return(NULL);
  }

  while (1) {
    if (loop->Go_Left) {
      /*
       *  Move node ptr to the left-most entry in this subtree.
       *  When the left-most entry is reached, return its data through
       *  the Go_Mid section of code.
       */

      loop->Go_Left = FALSE;
      loop->Go_Mid = TRUE;
      while (node->Left != NULL) {
        node = node->Left;
      }
    }
    else if (loop->Go_Mid) {
      /*
       *  Set the Current_Node pointer to node so that, on the
       *  next loop iteration, the routine can continue traversing the
       *  tree starting at Current_Node.  Set the direction bits
       *  so that, in the next loop iteration, the right subtree of node
       *  is traversed.
       *  Return a pointer to the data pointed to by node.
       */
      loop->Go_Mid = FALSE;
      loop->Go_Right = TRUE;
      loop->Current_Node = node;
      if (node != NULL)
        return(node->Vertex);
      else {
        LB_FREE(loop_control);
        return(NULL);
      }
    }
    else if (loop->Go_Right) {
      /*
       *  If node has a right subtree, set node to the root of
       *  that subtree and perform inorder traversal of that subtree (i.e.,
       *  set Go_Left bit to process its left subtree as if recursively).
       *  If node does not have a right subtree, node and all its
       *  children have been processed, so move up the tree to other
       *  unprocessed entries.
       */
      loop->Go_Right = FALSE;
      if (node->Right != NULL) {
        node = node->Right;
        loop->Go_Left = TRUE;
      }
      else {
        loop->Go_Up = TRUE;
      }
    }
    else /* loop->Go_Up */ {
      /*
       *  Move up the tree to an unprocessed parent.  All parents of
       *  right subtrees have already been processed (since the traversal
       *  is inorder), so move up until the node is the left child
       *  of its parent.  Then process the parent in Go_Mid.
       */
      loop->Go_Up = FALSE;
      loop->Go_Mid = TRUE;
      while ((node->Parent != NULL) &&
             (node == node->Parent->Right)) {
        node = node->Parent;
      }
      node = node->Parent;
    }
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/***************************  AVL ROUTINES ***********************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void single_LL_rotation(TREE **p_root, TREE *p1, TREE *p2)
{
/* Perform single LL rotation for avl rebalance */

  (*p_root)->Left = p2;
  p1->Right = *p_root;
  p1->Parent = (*p_root)->Parent;
  (*p_root)->Parent = p1;
  if (p2 != NULL)
    p2->Parent = *p_root;
  *p_root = p1;
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void double_LR_rotation(TREE **p_root, TREE *p1, TREE *p2)
{
/* Perform double LR rotation for avl rebalance */
TREE *p3, *p4;

  p3 = p2->Left;
  p4 = p2->Right;
  p1->Right = p3;
  p2->Left = p1;
  (*p_root)->Left = p4;
  p2->Right = *p_root;
  p2->Parent = (*p_root)->Parent;
  (*p_root)->Parent = p2;
  p1->Parent = p2;
  if (p3 != NULL)
    p3->Parent = p1;
  if (p4 != NULL)
    p4->Parent = *p_root;
  *p_root = p2;
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void single_RR_rotation(TREE **p_root, TREE *p1, TREE *p2)
{
/* perform single RR rotation for avl rebalance */
  (*p_root)->Right = p1->Left;
  p1->Left = *p_root;
  p1->Parent = (*p_root)->Parent;
  (*p_root)->Parent = p1;
  if (p2 != NULL)
    p2->Parent = *p_root;
  *p_root = p1;
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void double_RL_rotation(TREE **p_root, TREE *p1, TREE *p2)
{
/* perform double RL rotation for avl rebalance */
TREE *p3, *p4;

  p3 = p2->Right;
  p4 = p2->Left;
  p1->Left = p3;
  p2->Right = p1;
  (*p_root)->Right = p4;
  p2->Left = *p_root;
  p2->Parent = (*p_root)->Parent;
  (*p_root)->Parent = p2;
  p1->Parent = p2;
  if (p3 != NULL)
    p3->Parent = p1;
  if (p4 != NULL)
    p4->Parent = *p_root;
  *p_root = p2;
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static BOOLEAN avl_insert(TREE **p_root, TREE *p_parent, 
                          VERTEX *p_data, ID_INT_FN *fn, BOOLEAN *h)
{
TREE *p1, *p2, *p3, *p4;
BOOLEAN rc;
int comparison;

  if (*p_root == NULL) {
    /* element is not in the tree; insert it */
    *p_root = (TREE *)LB_MALLOC(sizeof(TREE));
    if (*p_root == NULL) {
      fprintf(stderr, "Out of memory in avl_insert\n");
      exit(-1);
    }
    (*p_root)->Vertex = p_data;
    (*p_root)->Left = NULL;
    (*p_root)->Right = NULL;
    (*p_root)->Bal = 0;
    (*p_root)->Parent = p_parent;

    *h = TRUE;
    rc = TRUE;
  } 
  else {
    comparison = fn(&(p_data->Id), &((*p_root)->Vertex->Id));
    if (comparison < 0) {
      rc = avl_insert(&((*p_root)->Left), (*p_root), p_data, fn, h);
      if (*h) {
        /* Left branch has grown higher */
        switch ((*p_root)->Bal) {
        case 1:
          (*p_root)->Bal = 0;
          *h = FALSE;
          break;

        case 0:
          (*p_root)->Bal = -1;
          break;

        case -1:
          p1 = (*p_root)->Left;
          p2 = p1->Right;
          if (p1->Bal == -1) {
            /* single LL rotation */
            (*p_root)->Bal = 0;
            single_LL_rotation(p_root, p1, p2);
          } 
          else {
            /* double LR rotation */
            if (p2->Bal == -1) {
              (*p_root)->Bal = 1;
            } else {
              (*p_root)->Bal = 0;
            }
            if (p2->Bal == 1) {
              p1->Bal = -1;
            } else {
              p1->Bal = 0;
            }
            double_LR_rotation(p_root, p1, p2);
          }

          (*p_root)->Bal = 0;
          *h = FALSE;
          break;
        }
      }
    } 
    else if (comparison > 0) {
      rc = avl_insert(&((*p_root)->Right), (*p_root), p_data, fn, h);
      if (*h) {
        /* right branch has grown higher */
        switch ((*p_root)->Bal) {
        case -1:
          (*p_root)->Bal = 0;
          *h = FALSE;
          break;

        case 0:
          (*p_root)->Bal = 1;
          break;

        case 1:
          p1 = (*p_root)->Right;
          p2 = p1->Left;
          if (p1->Bal == 1) {
            /* single RR rotation */
            (*p_root)->Bal = 0;
            single_RR_rotation(p_root, p1, p2);
          } else {
            /* double RL rotation */
            if (p2->Bal == 1) {
              (*p_root)->Bal = -1;
            } else {
              (*p_root)->Bal = 0;
            }
            if (p2->Bal == -1) {
              p1->Bal = 1;
            } else {
              p1->Bal = 0;
            }
            double_RL_rotation(p_root, p1, p2);
          }

          (*p_root)->Bal = 0;
          *h = FALSE;
          break;
        }
      }
    } else {
      /* duplicate entry; don't insert */
      *h = FALSE;
      rc = FALSE;
    }
  }

  return(rc);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void avl_balance1(TREE **p_root, BOOLEAN *h)
{
TREE *p1, *p2;
int b1, b2;


  /* h is true if left branch has become less high */
  switch ((*p_root)->Bal) {
  case -1:
    (*p_root)->Bal = 0;
    break;

  case 0:
    (*p_root)->Bal = 1;
    *h = FALSE;
    break;

  case 1:
    p1 = (*p_root)->Right;
    p2 = p1->Left;
    b1 = p1->Bal;
    if (b1 >= 0) {
      /* single RR rotation */
      if (b1 == 0) {
        (*p_root)->Bal = 1;
        p1->Bal = -1;
        *h = FALSE;
      } else {
        (*p_root)->Bal = 0;
        p1->Bal = 0;
      }
      single_RR_rotation(p_root, p1, p2);
    } else {
      /* double RL rotation */
      b2 = p2->Bal;
      if (b2 == 1) {
        (*p_root)->Bal = -1;
      } else {
        (*p_root)->Bal = 0;
      }
      if (b2 == -1) {
        p1->Bal = 1;
      } else {
        p1->Bal = 0;
      }
      p2->Bal = 0;
      double_RL_rotation(p_root, p1, p2);
    }
    break;
  }
}

/******************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void avl_balance2(TREE **p_root, BOOLEAN *h)
{
TREE *p1, *p2;
int b1, b2;


  /* h is true if right branch has become less high */
  switch ((*p_root)->Bal) {
  case 1:
    (*p_root)->Bal = 0;
    break;

  case 0:
    (*p_root)->Bal = -1;
    *h = FALSE;
    break;

  case -1:
    p1 = (*p_root)->Left;
    p2 = p1->Right;
    b1 = p1->Bal;
    if (b1 <= 0) {
      /* single LL rotation */
      if (b1 == 0) {
        (*p_root)->Bal = -1;
        p1->Bal = 1;
        *h = FALSE;
      } else {
        (*p_root)->Bal = 0;
        p1->Bal = 0;
      }
      single_LL_rotation(p_root, p1, p2);
    } else {
      /* double LR rotation */
      b2 = p2->Bal;
      if (b2 == -1) {
        (*p_root)->Bal = 1;
      } else {
        (*p_root)->Bal = 0;
      }
      if (b2 == 1) {
        p1->Bal = -1;
      } else {
        p1->Bal = 0;
      }
      p2->Bal = 0;
      double_LR_rotation(p_root, p1, p2);
    }
    break;
  }
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static TREE *q;

static void avl_del(TREE **p_root, BOOLEAN *h)
{
TREE *x;

  if ((*p_root)->Right != NULL) {
    avl_del(&((*p_root)->Right), h);
    if (*h) {
      avl_balance2(p_root, h);
    }
  } else {
    *h = TRUE;
    x = *p_root;
    q->Vertex = x->Vertex;
    *p_root = x->Left;
    if (*p_root != NULL)
      (*p_root)->Parent = q;
    LB_FREE(&x);
  }
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void *avl_delete(TREE **p_root, VERTEX *p_data, 
                        ID_INT_FN *fn, BOOLEAN *h)
{
void *rc;
int comparison;

  if (*p_root) {
    comparison = fn(&(p_data->Id), &((*p_root)->Vertex->Id));
    if (comparison < 0) {
      rc = avl_delete(&((*p_root)->Left), p_data, fn, h);
      if (*h) {
        avl_balance1(p_root, h);
      }
    } else if (comparison > 0) {
      rc = avl_delete(&((*p_root)->Right), p_data, fn, h);
      if (*h) {
        avl_balance2(p_root, h);
      }
    } else {
      rc = (*p_root)->Vertex;
      q = *p_root;
      if (q->Right == NULL) {
        *p_root = q->Left;
        if (*p_root != NULL)
          (*p_root)->Parent = q->Parent;
        *h = TRUE;
        LB_FREE(&q);
      } else if (q->Left == NULL) {
        *p_root = q->Right;
        if (*p_root != NULL)
          (*p_root)->Parent = q->Parent;
        *h = TRUE;
        LB_FREE(&q);
      } else {
        avl_del(&((*p_root)->Left), h);
        if (*h) {
          avl_balance1(p_root, h);
        }
      }
    }
  } else {
    *h = FALSE;
    rc = NULL;
  }

  return(rc);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/************************  DEBUGGING ROUTINES *******************************/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

#include "vx_const.h"

static void print_recursively(TREE *p_root, int level)
{
int i;

  if (p_root) {

    print_recursively(p_root->Right, level+1);

    for (i = 0; i < 8 * level; i++) printf(" ");
    BL_ID_Util.Print_ID(&(p_root->Vertex->Id));
    printf("[");
    if (p_root->Parent)
      BL_ID_Util.Print_ID(&( p_root->Parent->Vertex->Id));
    printf("]\n");

    print_recursively(p_root->Left, level+1);
  }
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void print_tree_graph(TREE *p_root)
{
  /*
   * Code for DYNAMIC_STRUCT == TREE_STRUCT
   */
  printf("-----------------------------------------\n");
  print_recursively(p_root, 0);
  printf("-----------------------------------------\n");
}

