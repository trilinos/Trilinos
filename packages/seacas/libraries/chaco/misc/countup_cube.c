/*
 * Copyright(C) 1999-2021 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "smalloc.h" // for sfree, smalloc
#include "structs.h" // for vtx_data
#include <stdio.h>   // for fprintf, printf, FILE, NULL
#include <stdlib.h>

/* Print metrics of partition quality. */

void countup_cube(struct vtx_data **graph,      /* graph data structure */
                  int               nvtxs,      /* number of vtxs in graph */
                  int              *assignment, /* set number of each vtx (length nvtxs+1) */
                  int               ndims,      /* number of cuts at each level */
                  int               ndims_tot,  /* total number of divisions of graph */
                  int               print_lev,  /* level of output */
                  FILE             *outfile,    /* output file if not NULL */
                  int               using_ewgts /* are edge weights being used? */
)
{
  int print2file = (outfile != NULL);

  int     nsets   = (1 << ndims_tot);
  double *cutsize = smalloc(nsets * sizeof(double));
  double *hopsize = smalloc(nsets * sizeof(double));
  int    *setsize = smalloc(nsets * sizeof(int));

  int *setseen  = smalloc(nsets * sizeof(int));
  int *startptr = smalloc((nsets + 1) * sizeof(int));
  int *inorder  = smalloc(nvtxs * sizeof(int));
  for (int j = 0; j < nsets; j++) {
    setsize[j] = 0;
  }
  for (int i = 1; i <= nvtxs; i++) {
    ++setsize[assignment[i]];
  }

  /* Modify setsize to become index into vertex list. */
  for (int j = 1; j < nsets; j++) {
    setsize[j] += setsize[j - 1];
  }
  for (int j = nsets - 1; j > 0; j--) {
    startptr[j] = setsize[j] = setsize[j - 1];
  }
  startptr[0] = setsize[0] = 0;
  startptr[nsets]          = nvtxs;
  for (int i = 1; i <= nvtxs; i++) {
    int set               = assignment[i];
    inorder[setsize[set]] = i;
    setsize[set]++;
  }

  int start_dims;
  int level;
  if (abs(print_lev) > 1) { /* Print data from all levels of recursion. */
    start_dims = ndims;
    level      = 0;
  }
  else { /* Only print data from final level. */
    start_dims = ndims_tot;
    level      = (ndims_tot + ndims - 1) / ndims - 1;
  }
  int k = start_dims;
  while (k <= ndims_tot) {
    level++;
    nsets = (1 << k);
    for (int j = 0; j < nsets; j++) {
      cutsize[j] = 0;
      hopsize[j] = 0;
      setsize[j] = 0;
    }
    int mask = 0;
    for (int j = 0; j < k; j++) {
      mask = (mask << 1) + 1;
    }

    for (int i = 1; i <= nvtxs; i++) {
      int set = assignment[i] & mask;
      setsize[set] += graph[i]->vwgt;
      for (int j = 1; j < graph[i]->nedges; j++) {
        int neighbor = graph[i]->edges[j];
        int set2     = assignment[neighbor] & mask;
        if (set != set2) {
          double ewgt = 1;
          if (using_ewgts) {
            ewgt = graph[i]->ewgts[j];
          }
          cutsize[set] += ewgt;
          int bits = set ^ set2;
          for (int l = bits; l; l >>= 1) {
            if (l & 1) {
              hopsize[set] += ewgt;
            }
          }
        }
      }
    }

    int tot_size = 0;
    int max_size = 0;
    for (int set = 0; set < nsets; set++) {
      tot_size += setsize[set];
      if (setsize[set] > max_size) {
        max_size = setsize[set];
      }
    }

    int min_size = max_size;
    for (int set = 0; set < nsets; set++) {
      if (setsize[set] < min_size) {
        min_size = setsize[set];
      }
    }

    double ncuts           = 0;
    double nhops           = 0;
    double total_bdyvtxs   = 0;
    int    total_neighbors = 0;
    double bdyvtx_hops_tot = 0;
    double bdyvtx_hops_max = 0;
    double bdyvtx_hops_min = 0;
    double maxcuts         = 0;
    double mincuts         = 0;
    double maxhops         = 0;
    double minhops         = 0;
    int    total_internal  = 0;
    int    min_internal    = max_size;
    int    max_internal    = 0;
    double maxbdy          = 0;
    double minbdy          = 0;
    int    maxneighbors    = 0;
    int    minneighbors    = 0;

    printf("\nAfter level %d  (nsets = %d):\n", level, nsets);
    if (print2file) {
      fprintf(outfile, "\nAfter level %d  (nsets = %d):\n", level, nsets);
    }
    if (print_lev < 0) {
      printf("    set    size      cuts       hops   bndy_vtxs    adj_sets\n");
      if (print2file) {
        fprintf(outfile, "    set    size      cuts       hops   bndy_vtxs    adj_sets\n");
      }
    }
    for (int set = 0; set < nsets; set++) {
      int internal = setsize[set];
      for (int i = 0; i < nsets; i++) {
        setseen[i] = 0;
      }
      /* Compute number of set neighbors, and number of vtxs on boundary. */
      /* Loop through multiple assignments defining current set. */
      int bdyvtxs     = 0;
      int bdyvtx_hops = 0;
      for (int l = 0; l < (1 << (ndims_tot - k)); l++) {
        int set2 = (l << k) + set;
        for (int i = startptr[set2]; i < startptr[set2 + 1]; i++) {
          int onbdy = 0;
          int vtx   = inorder[i];
          for (int j = 1; j < graph[vtx]->nedges; j++) {
            int neighbor = graph[vtx]->edges[j];
            int set3     = assignment[neighbor] & mask;
            if (set3 != set) { /* Is vtx on boundary? */
              /* Has this neighboring set been seen already? */
              if (setseen[set3] >= 0) {
                int bits = set ^ set3;
                for (int ll = bits; ll; ll >>= 1) {
                  if (ll & 1) {
                    ++bdyvtx_hops;
                  }
                }
                ++onbdy;
                setseen[set3] = -setseen[set3] - 1;
              }
            }
          }
          /* Now reset all the setseen values to be positive. */
          if (onbdy != 0) {
            for (int j = 1; j < graph[vtx]->nedges; j++) {
              int neighbor = graph[vtx]->edges[j];
              int set3     = assignment[neighbor] & mask;
              if (setseen[set3] < 0) {
                setseen[set3] = -setseen[set3];
              }
            }
            internal -= graph[vtx]->vwgt;
          }
          bdyvtxs += onbdy;
        }
      }

      total_internal += internal;
      bdyvtx_hops_tot += bdyvtx_hops;
      if (bdyvtx_hops > bdyvtx_hops_max) {
        bdyvtx_hops_max = bdyvtx_hops;
      }
      if (set == 0 || bdyvtx_hops < bdyvtx_hops_min) {
        bdyvtx_hops_min = bdyvtx_hops;
      }
      if (internal > max_internal) {
        max_internal = internal;
      }
      if (set == 0 || internal < min_internal) {
        min_internal = internal;
      }

      /* Now count up the number of neighboring sets. */
      int neighbor_sets = 0;
      for (int i = 0; i < nsets; i++) {
        if (setseen[i] != 0) {
          ++neighbor_sets;
        }
      }

      if (print_lev < 0) {
        printf(" %5d    %5d    %6g     %6g   %6d      %6d\n", set, setsize[set], cutsize[set],
               hopsize[set], bdyvtxs, neighbor_sets);
        if (print2file) {
          fprintf(outfile, " %5d    %5d    %6g     %6g   %6d      %6d\n", set, setsize[set],
                  cutsize[set], hopsize[set], bdyvtxs, neighbor_sets);
        }
      }
      if (cutsize[set] > maxcuts) {
        maxcuts = cutsize[set];
      }
      if (set == 0 || cutsize[set] < mincuts) {
        mincuts = cutsize[set];
      }
      if (hopsize[set] > maxhops) {
        maxhops = hopsize[set];
      }
      if (set == 0 || hopsize[set] < minhops) {
        minhops = hopsize[set];
      }
      if (bdyvtxs > maxbdy) {
        maxbdy = bdyvtxs;
      }
      if (set == 0 || bdyvtxs < minbdy) {
        minbdy = bdyvtxs;
      }
      if (neighbor_sets > maxneighbors) {
        maxneighbors = neighbor_sets;
      }
      if (set == 0 || neighbor_sets < minneighbors) {
        minneighbors = neighbor_sets;
      }
      ncuts += cutsize[set];
      nhops += hopsize[set];
      total_bdyvtxs += bdyvtxs;
      total_neighbors += neighbor_sets;
    }
    ncuts /= 2;
    nhops /= 2;

    printf("\n");
    printf("                            Total      Max/Set      Min/Set\n");
    printf("                            -----      -------      -------\n");
    printf("Set Size:             %11d  %11d  %11d\n", tot_size, max_size, min_size);
    printf("Edge Cuts:            %11g  %11g  %11g\n", ncuts, maxcuts, mincuts);
    printf("Hypercube Hops:       %11g  %11g  %11g\n", nhops, maxhops, minhops);
    printf("Boundary Vertices:    %11g  %11g  %11g\n", total_bdyvtxs, maxbdy, minbdy);
    printf("Boundary Vertex Hops: %11g  %11g  %11g\n", bdyvtx_hops_tot, bdyvtx_hops_max,
           bdyvtx_hops_min);
    printf("Adjacent Sets:        %11d  %11d  %11d\n", total_neighbors, maxneighbors, minneighbors);
    printf("Internal Vertices:    %11d  %11d  %11d\n\n", total_internal, max_internal,
           min_internal);

    if (print2file) {
      fprintf(outfile, "\n");
      fprintf(outfile, "                            Total      Max/Set      Min/Set\n");
      fprintf(outfile, "                            -----      -------      -------\n");
      fprintf(outfile, "Set Size:             %11d  %11d  %11d\n", tot_size, max_size, min_size);
      fprintf(outfile, "Edge Cuts:            %11g  %11g  %11g\n", ncuts, maxcuts, mincuts);
      fprintf(outfile, "Hypercube Hops:       %11g  %11g  %11g\n", nhops, maxhops, minhops);
      fprintf(outfile, "Boundary Vertices:    %11g  %11g  %11g\n", total_bdyvtxs, maxbdy, minbdy);
      fprintf(outfile, "Boundary Vertex Hops: %11g  %11g  %11g\n", bdyvtx_hops_tot, bdyvtx_hops_max,
              bdyvtx_hops_min);
      fprintf(outfile, "Adjacent Sets:        %11d  %11d  %11d\n", total_neighbors, maxneighbors,
              minneighbors);
      fprintf(outfile, "Internal Vertices:    %11d  %11d  %11d\n\n", total_internal, max_internal,
              min_internal);
    }

    if (k == ndims_tot) {
      k++;
    }
    else {
      k += ndims;
      if (k > ndims_tot) {
        k = ndims_tot;
      }
    }
  }

  sfree(cutsize);
  sfree(hopsize);
  sfree(setsize);
  sfree(setseen);
  sfree(startptr);
  sfree(inorder);
}
