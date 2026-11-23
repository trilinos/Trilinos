#include "exodusII.h"
#include <stdio.h>

#define EXCHECK(funcall)                                                                           \
  do {                                                                                             \
    if ((funcall) != EX_NOERR) {                                                                   \
      fprintf(stderr, "%s failed\n", #funcall);                                                    \
      return 1;                                                                                    \
    }                                                                                              \
  } while (0)

int main(void)
{
  int cpu_word_size = 8;
  int io_word_size  = 8;

  int exoid = ex_create("side-set-tri4.exo", EX_CLOBBER, &cpu_word_size, &io_word_size);
  if (exoid < 0) {
    return 1;
  }

  EXCHECK(ex_put_init(exoid, "TRI4 shell faces", 3, 4, 1, 1, 0, 1));
  double x[4] = {0, 3, 0, 1};
  double y[4] = {0, 0, 3, 1};
  double z[4] = {0, 0, 0, 0};
  EXCHECK(ex_put_coord(exoid, x, y, z));
  EXCHECK(ex_put_block(exoid, EX_ELEM_BLOCK, 1, "TRI4", 1, 4, 0, 0, 0));
  int connect[4] = {1, 2, 3, 4};
  EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, 1, connect, NULL, NULL));

  /* The two faces have opposite corner orderings and share the center node. */
  int elements[2] = {1, 1};
  int sides[2]    = {1, 2};
  EXCHECK(ex_put_set_param(exoid, EX_SIDE_SET, 1, 2, 0));
  EXCHECK(ex_put_set(exoid, EX_SIDE_SET, 1, elements, sides));
  EXCHECK(ex_close(exoid));

  float version;
  io_word_size = 0;
  exoid        = ex_open("side-set-tri4.exo", EX_READ, &cpu_word_size, &io_word_size, &version);
  if (exoid < 0) {
    return 1;
  }

  int length    = 0;
  int counts[2] = {0, 0};
  int nodes[8]  = {-1, -1, -1, -1, -1, -1, -1, -1};
  EXCHECK(ex_get_side_set_node_list_len(exoid, 1, &length));
  EXCHECK(ex_get_side_set_node_list(exoid, 1, counts, nodes));
  EXCHECK(ex_close(exoid));

  printf("TRI4 faces: length=%d counts=[%d,%d] nodes=[%d,%d,%d,%d,%d,%d,%d,%d]\n", length,
         counts[0], counts[1], nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6],
         nodes[7]);
  if (length != 8 || counts[0] != 4 || counts[1] != 4) {
    return 1;
  }
  int expected[8] = {1, 2, 3, 4, 3, 2, 1, 4};
  for (int i = 0; i < 8; i++) {
    if (nodes[i] != expected[i]) {
      fprintf(stderr, "node %d: expected %d, got %d\n", i, expected[i], nodes[i]);
      return 1;
    }
  }
  return 0;
}
