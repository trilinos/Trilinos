/*
 * Copyright(C) 1999-2025 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*
 * test-tet15-sideset - side set node lists of 14- and 15-node tetrahedra
 *
 * A 14- or 15-node tetrahedron has seven nodes per side: three vertices,
 * three mid-edge nodes and one mid-face node.  ex_get_side_set_node_list
 * reported seven nodes per side for these elements but filled only six,
 * leaving the seventh entry of every side uninitialized.  This test writes
 * one element of each type with a side set over all four sides and checks
 * every entry of the node list against the element connectivity.
 */

#include "exodusII.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#undef NDEBUG
#include <assert.h>

/* Side node table of a 14- or 15-node tetrahedron: vertices, mid-edge nodes,
 * mid-face node, for sides 1 to 4 (1-based node numbers). */
static const int tet_side_nodes[4][7] = {{1, 2, 4, 5, 9, 8, 14},
                                         {2, 3, 4, 6, 10, 9, 12},
                                         {1, 4, 3, 8, 10, 7, 13},
                                         {1, 3, 2, 7, 6, 5, 11}};

static int test_tet(const char *filename, const char *elem_type, int num_nodes_per_elem)
{
  int CPU_word_size = sizeof(double);
  int IO_word_size  = sizeof(double);

  int exoid = ex_create(filename, EX_CLOBBER, &CPU_word_size, &IO_word_size);
  assert(exoid >= 0);

  ex_init_params par = {.title         = "tetra side sets",
                        .num_dim       = 3,
                        .num_nodes     = num_nodes_per_elem,
                        .num_elem      = 1,
                        .num_elem_blk  = 1,
                        .num_node_sets = 0,
                        .num_side_sets = 1};
  assert(ex_put_init_ext(exoid, &par) == EX_NOERR);

  double x[15], y[15], z[15];
  for (int i = 0; i < num_nodes_per_elem; i++) {
    x[i] = i;
    y[i] = 2 * i;
    z[i] = 3 * i;
  }
  assert(ex_put_coord(exoid, x, y, z) == EX_NOERR);

  assert(ex_put_block(exoid, EX_ELEM_BLOCK, 10, elem_type, 1, num_nodes_per_elem, 0, 0, 0) ==
         EX_NOERR);

  /* Connectivity of the one element: local node i is global node i. */
  int connect[15];
  for (int i = 0; i < num_nodes_per_elem; i++) {
    connect[i] = i + 1;
  }
  assert(ex_put_conn(exoid, EX_ELEM_BLOCK, 10, connect, NULL, NULL) == EX_NOERR);

  int elem_list[4] = {1, 1, 1, 1};
  int side_list[4] = {1, 2, 3, 4};
  assert(ex_put_set_param(exoid, EX_SIDE_SET, 20, 4, 0) == EX_NOERR);
  assert(ex_put_set(exoid, EX_SIDE_SET, 20, elem_list, side_list) == EX_NOERR);
  assert(ex_close(exoid) == EX_NOERR);

  /* Read back. */
  CPU_word_size = 0;
  IO_word_size  = 0;
  float version;
  exoid = ex_open(filename, EX_READ, &CPU_word_size, &IO_word_size, &version);
  assert(exoid >= 0);

  int node_list_len = 0;
  assert(ex_get_side_set_node_list_len(exoid, 20, &node_list_len) == EX_NOERR);
  assert(node_list_len == 4 * 7);

  int node_cnt_list[4];
  int node_list[4 * 7];
  for (int i = 0; i < 4 * 7; i++) {
    node_list[i] = -1; /* a sentinel that no connectivity entry equals */
  }
  assert(ex_get_side_set_node_list(exoid, 20, node_cnt_list, node_list) == EX_NOERR);

  int errors = 0;
  int pos    = 0;
  for (int s = 0; s < 4; s++) {
    if (node_cnt_list[s] != 7) {
      fprintf(stderr, "%s side %d: %d nodes per side, expected 7\n", elem_type, s + 1,
              node_cnt_list[s]);
      errors++;
    }
    for (int k = 0; k < 7; k++, pos++) {
      int expected = connect[tet_side_nodes[s][k] - 1];
      if (node_list[pos] != expected) {
        fprintf(stderr, "%s side %d node %d: got %d, expected %d\n", elem_type, s + 1, k + 1,
                node_list[pos], expected);
        errors++;
      }
    }
  }
  assert(ex_close(exoid) == EX_NOERR);
  unlink(filename);
  return errors;
}

int main(int argc, char **argv)
{
  ex_opts(EX_VERBOSE);
  int errors = 0;
  errors += test_tet("test-tet14-sideset.exo", "TETRA14", 14);
  errors += test_tet("test-tet15-sideset.exo", "TETRA15", 15);
  if (errors == 0) {
    printf("test-tet15-sideset: all side set node lists match the connectivity\n");
  }
  return errors == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
