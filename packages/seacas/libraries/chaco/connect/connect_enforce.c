/*
 * Copyright(C) 1999-2020, 2023 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "defs.h"    // for TRUE
#include "smalloc.h" // for sfree, smalloc
#include "structs.h" // for vtx_data, heap
#include <stdio.h>   // for NULL, printf

/* Move vertices between domains to ensure each domain is connected. */
/* Note: This will likely result in load imbalance. */

void connect_enforce(struct vtx_data **graph,       /* data structure for graph */
                     int               nvtxs,       /* number of vertices in full graph */
                     int               using_ewgts, /* are edge weights being used? */
                     int              *assignment,  /* set number of each vtx (length n) */
                     double           *goal,        /* desired sizes for each set */
                     int               nsets_tot,   /* total number sets created */
                     int              *total_move,  /* total number of vertices moved */
                     int              *max_move     /* largest connected component moved */
)
{
  struct vtx_data **subgraph;           /* data structure for domain graph */
  int               subnvtxs;           /* number of vertices in a domain */
  int               subnedges;          /* number of edges in a domain */
  struct heap      *heap;               /* data structure for sorting set sizes */
  int              *heap_map;           /* pointers from sets to heap locations */
  int              *list_ptrs;          /* header of vtx list for each domain */
  int              *setlists;           /* linked list of vtxs for each domain */
  int              *vtxlist;            /* space for breadth first search list */
  int              *comp_lists;         /* list of vtxs in each connected comp */
  int              *clist_ptrs;         /* pointers to heads of comp_lists */
  int              *subsets;            /* list of active domains (all of them) */
  int              *subsets2;           /* list of active domains (all of them) */
  double           *set_size;           /* weighted sizes of different partitions */
  double            size;               /* size of subset being moved to new domain */
  int              *bndy_list;          /* list of domains adjacent to component */
  double           *bndy_size;          /* size of these boundaries */
  double           *comp_size;          /* sizes of different connected components */
  double            comp_size_max;      /* size of largest connected component */
  int               comp_max_index = 0; /* which component is largest? */
  int              *glob2loc;           /* global to domain renumbering */
  int              *loc2glob;           /* domain to global renumbering */
  int              *degree;             /* number of neighbors of a vertex */
  int              *comp_flag;          /* component number for each vtx */
  double            ewgt;               /* edge weight */
  int               nbndy;              /* number of sets adjacent to component */
  int               domain;             /* which subdomain I'm working on */
  int               new_domain;         /* subdomain to move some vertices to */
  double            max_bndy;           /* max connectivity to other domain */
  int               ncomps;             /* number of connected components */
  int               change;             /* how many vertices change set? */
  int               max_change;         /* largest subset moved together */
  int               vtx;                /* vertex in a connected component */
  int               set;                /* set a neighboring vertex is in */
  int               i, j, k, l;         /* loop counters */

  change     = 0;
  max_change = 0;

  /* Allocate space & initialize some values. */

  set_size  = smalloc(nsets_tot * sizeof(double));
  bndy_size = smalloc(nsets_tot * sizeof(double));
  bndy_list = smalloc(nsets_tot * sizeof(int));

  setlists  = smalloc((nvtxs + 1) * sizeof(int));
  list_ptrs = smalloc(nsets_tot * sizeof(int));

  glob2loc = smalloc((nvtxs + 1) * sizeof(int));
  loc2glob = smalloc((nvtxs + 1) * sizeof(int));
  subsets  = smalloc(nsets_tot * sizeof(int));
  heap     = (struct heap *)smalloc((nsets_tot + 1) * sizeof(struct heap));
  heap_map = smalloc(nsets_tot * sizeof(int));

  for (i = 0; i < nsets_tot; i++) {
    set_size[i]  = 0;
    bndy_list[i] = 0;
    bndy_size[i] = 0;
    subsets[i]   = i;
  }

  for (i = 1; i <= nvtxs; i++) {
    set_size[assignment[i]] += graph[i]->vwgt;
  }

  for (i = 0; i < nsets_tot; i++) {
    heap[i + 1].tag = i;
    heap[i + 1].val = set_size[i] - goal[i];
  }

  make_setlists(setlists, list_ptrs, nsets_tot, subsets, assignment, NULL, nvtxs, TRUE);

  heap_build(heap, nsets_tot, heap_map);

  for (i = 0; i < nsets_tot; i++) {
    /* Find largest remaining set to work on next */
    size = heap_extract_max(heap, nsets_tot - i, &domain, heap_map);

    /* Construct subdomain graph. */
    subnvtxs = make_maps(setlists, list_ptrs, domain, glob2loc, loc2glob);
    if (subnvtxs > 1) {

      subgraph = (struct vtx_data **)smalloc((subnvtxs + 1) * sizeof(struct vtx_data *));
      degree   = smalloc((subnvtxs + 1) * sizeof(int));

      make_subgraph(graph, subgraph, subnvtxs, &subnedges, assignment, domain, glob2loc, loc2glob,
                    degree, using_ewgts);

      /* Find connected components. */
      comp_flag = smalloc((subnvtxs + 1) * sizeof(int));
      vtxlist   = smalloc(subnvtxs * sizeof(int));
      ncomps    = find_comps(subgraph, subnvtxs, comp_flag, vtxlist);
      sfree(vtxlist);

      /* Restore original graph */
      remake_graph(subgraph, subnvtxs, loc2glob, degree, using_ewgts);
      sfree(degree);
      sfree(subgraph);

      if (ncomps > 1) {

        /* Figure out sizes of components */
        comp_size = smalloc(ncomps * sizeof(double));
        for (j = 0; j < ncomps; j++) {
          comp_size[j] = 0;
        }
        for (j = 1; j <= subnvtxs; j++) {
          comp_size[comp_flag[j]] += graph[loc2glob[j]]->vwgt;
        }
        comp_size_max = 0;
        for (j = 0; j < ncomps; j++) {
          if (comp_size[j] > comp_size_max) {
            comp_size_max  = comp_size[j];
            comp_max_index = j;
          }
        }
        for (j = 0; j < ncomps; j++) {
          if (j != comp_max_index) {
            change += comp_size[j];
            if (comp_size[j] > max_change) {
              max_change = comp_size[j];
            }
          }
        }
        sfree(comp_size);

        /* Make data structures for traversing components */
        comp_lists = smalloc((subnvtxs + 1) * sizeof(int));
        clist_ptrs = smalloc(ncomps * sizeof(int));
        if (ncomps > nsets_tot) {
          subsets2 = smalloc(ncomps * sizeof(int));
          for (j = 0; j < ncomps; j++) {
            subsets2[j] = j;
          }
        }
        else {
          subsets2 = subsets;
        }
        make_setlists(comp_lists, clist_ptrs, ncomps, subsets2, comp_flag, NULL, subnvtxs, TRUE);
        if (ncomps > nsets_tot) {
          sfree(subsets2);
        }

        /* Move all but the largest component. */
        ewgt = 1;
        for (j = 0; j < ncomps; j++) {
          if (j != comp_max_index) {

            /* Figure out to which other domain it is most connected. */
            nbndy = 0;
            k     = clist_ptrs[j];
            while (k != 0) {
              vtx = loc2glob[k];
              for (l = 1; l <= graph[vtx]->nedges; l++) {
                set = assignment[graph[vtx]->edges[l]];
                if (set != domain) {
                  if (bndy_size[set] == 0) {
                    bndy_list[nbndy++] = set;
                  }
                  if (using_ewgts) {
                    ewgt = graph[vtx]->ewgts[l];
                  }
                  bndy_size[set] += ewgt;
                }
              }

              k = comp_lists[k];
            }

            /* Select a new domain. */
            /* Instead of just big boundary, penalize too-large sets. */
            /* Could be more aggressive to improve balance. */
            max_bndy   = 0;
            new_domain = -1;
            for (k = 0; k < nbndy; k++) {
              l = bndy_list[k];
              if (bndy_size[l] * goal[l] / (set_size[l] + 1) > max_bndy) {
                new_domain = bndy_list[k];
                max_bndy   = bndy_size[l] * goal[l] / (set_size[l] + 1);
              }
            }
            if (new_domain == -1) {
              printf("Error in connect_enforce: new_domain = -1.  Disconnected graph?\n");
              new_domain = domain;
            }

            /* Clear bndy_size array */
            for (k = 0; k < nbndy; k++) {
              bndy_size[bndy_list[k]] = 0;
            }

            k = clist_ptrs[j];

            size = 0;

            while (k != 0) {
              vtx             = loc2glob[k];
              assignment[vtx] = new_domain;
              size += graph[vtx]->vwgt;

              /* Finally, update setlists and list_ptrs */
              /* Note: current domain setlist now bad, but not used
                 again */
              setlists[vtx]         = list_ptrs[new_domain];
              list_ptrs[new_domain] = vtx;

              k = comp_lists[k];
            }
            /*
            printf("Updating set %d (from %d) to size %g\n", new_domain, domain,
            set_size[new_domain] + size - goal[new_domain]);
            */
            if (heap_map[new_domain] > 0) {
              set_size[new_domain] += size;
              heap_update_val(heap, heap_map[new_domain], set_size[new_domain] - goal[new_domain],
                              heap_map);
            }
          }
        }

        sfree(clist_ptrs);
        sfree(comp_lists);
      }
      sfree(comp_flag);
    }
  }

  sfree(heap_map);
  sfree(heap);
  sfree(subsets);
  sfree(loc2glob);
  sfree(glob2loc);
  sfree(list_ptrs);
  sfree(setlists);
  sfree(bndy_list);
  sfree(bndy_size);
  sfree(set_size);

  *total_move = change;
  *max_move   = max_change;
}
