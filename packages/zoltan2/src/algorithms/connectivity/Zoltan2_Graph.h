#ifndef __graph_h__
#define __graph_h__

struct graph {
  int n;
  unsigned m;
  int* out_array;
  unsigned* out_degree_list;
  int max_degree_vert;
  double avg_out_degree;
};

#define out_degree(g, n) (g->out_degree_list[n+1] - g->out_degree_list[n])
#define out_vertices(g, n) (&g->out_array[g->out_degree_list[n]])

#endif
