// Copyright(C) 1999-2020, 2022, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifdef SEACAS_HAVE_MPI
#include "mpi.h"
#endif
#include "gtest/gtest.h"

#include "ElementPartition.h"

#include "exodusII.h" // for ex_opts, ex_int64_status, etc

#include <cstdint> // for int64_t
#include <cstdio>  // for stderr, etc
#include <cstdlib> // for exit
#include <cstring>
#include <unistd.h> // for getopt, optarg, optind

#include <Ionit_Initializer.h>
#include <Ioss_CodeTypes.h>
#include <Ioss_Enumerate.h>
#include <Ioss_SmartAssert.h>
#include <Ioss_SubSystem.h>
#include <Ioss_Transform.h>
#include <Ioss_Utils.h>
#include <tokenize.h>

#include <cstdarg> // for va_end, va_arg, va_list, etc
#include <cstddef> // for size_t
#include <cstdio>  // for stderr
#include <cstdlib> // for exit, malloc
#include <fmt/format.h>
#include <fmt/ostream.h>

#include "DataUtils.h"
#include "ElementUtils.h"

namespace utest_util {

  template ElementPartition<int>::ElementPartition(IossMesh                      *mesh,
                                                   const std::vector<EntityProc> &procAssign,
                                                   const int                      nProc);
  template ElementPartition<int64_t>::ElementPartition(IossMesh                      *mesh,
                                                       const std::vector<EntityProc> &procAssign,
                                                       const int                      nProc);

  template <typename INT>
  ElementPartition<INT>::ElementPartition(IossMesh *mesh, const std::vector<EntityProc> &procAssign,
                                          const int nProc)
      : Partition(mesh, procAssign, nProc)
  {
    verify_input_partition();
    generate_partition_maps();
  }

  template <typename INT> size_t ElementPartition<INT>::get_num_local_entities() const
  {
    return m_mesh->get_num_local_elements();
  }

  template <typename INT> int ElementPartition<INT>::get_processor(size_t localId) const
  {
    IossElementData elem = m_mesh->get_local_element(localId);
    EXPECT_TRUE(elem.is_valid());

    auto lowerBound = std::lower_bound(m_entityPartition.begin(), m_entityPartition.end(), elem.id,
                                       EntityProcLess());
    EXPECT_FALSE(lowerBound == m_entityPartition.end() || lowerBound->id != elem.id)
        << "element " << elem.id << " was not assigned a processor";

    return lowerBound->proc;
  }

  template <typename INT> int ElementPartition<INT>::get_node_conn_index(size_t elem, INT node)
  {
    IossElementData elemData = m_mesh->get_local_element(elem);
    int             nnodes   = elemData.topology->number_nodes();

    for (int i = 0; i < nnodes; i++) {
      if (elemData.localConnectivity[i] == node) {
        return i + 1;
      }
    }
    return -1;
  }

  template <typename INT>
  void ElementPartition<INT>::find_beam_internal_and_border_elements(
      const std::vector<std::vector<INT>> &sur_elem, const int max_nsur)
  {
    size_t numNodes = m_mesh->get_num_local_nodes();
    size_t numElems = m_mesh->get_num_local_elements();

    std::vector<int> categorized(numElems, 0);
    for (size_t ecnt = 0; ecnt < numElems; ecnt++) {
      IossElementData elemData = m_mesh->get_local_element(ecnt);
      int             nsides   = elemData.topology->number_boundaries();
      assert(nsides == 2);

      int proce = get_processor(ecnt);
      /* check each side of this element */
      for (int nscnt = 0; nscnt < nsides; nscnt++) {

        /* get the node on this element side (should only be one)*/
        std::vector<utest_util::EntityId> side_nodes;

        get_local_side_nodes(elemData, nscnt + 1 /* 1 based */, side_nodes);
        assert(side_nodes.size() == 1);

        size_t nhold = sur_elem[side_nodes[0]].size();
        assert(nhold == 1 || nhold == 2);
        if (nhold == 2) {
          // 2 elements connected to this node -- `ecnt` and the other one...
          for (size_t ncnt = 0; ncnt < nhold; ncnt++) {
            size_t elem = sur_elem[side_nodes[0]][ncnt];
            if (elem == ecnt) {
              continue;
            }
            int proc2 = get_processor(elem);
            if (proce == proc2) {
              continue;
            }
            else {
              // Processors of the two elements are different, so we are
              // at a processor boundary...
              if (!categorized[ecnt]) {
                bor_elems[proce].push_back(ecnt);
                categorized[ecnt] = 1;
              }
              auto nodeConnIndex = get_node_conn_index(elem, side_nodes[0]);
              e_cmap_elems[proc2].push_back(elem);
              e_cmap_sides[proc2].push_back(nodeConnIndex);
              e_cmap_procs[proc2].push_back(proce);
              e_cmap_neigh[proc2].push_back(ecnt);
            }
          }
        }
      }
      if (!categorized[ecnt]) {
        int_elems[proce].push_back(ecnt);
        categorized[ecnt] = 1;
      }
    }
  }

  template <typename INT>
  void ElementPartition<INT>::categorize_elements(const std::vector<std::vector<INT>> &sur_elem,
                                                  const int                            max_nsur)
  {
    int mesh_dim = m_mesh->get_spatial_dimension();

    if (mesh_dim == 1) {
      find_beam_internal_and_border_elements(sur_elem, max_nsur);
      return;
    }

    /* allocate space to hold info about surrounding elements */
    std::vector<INT> pt_list(max_nsur, 0);
    std::vector<INT> hold_elem(max_nsur, 0);

    size_t numNodes = m_mesh->get_num_local_nodes();
    size_t numElems = m_mesh->get_num_local_elements();

    /* Find the internal and border elements */
    std::vector<utest_util::EntityId> side_nodes;
    for (size_t ecnt = 0; ecnt < numElems; ecnt++) {

      IossElementData elemData = m_mesh->get_local_element(ecnt);
      int             nsides   = elemData.topology->number_boundaries();
      int             proc     = get_processor(ecnt);
      assert(proc < m_numProcs);

      bool        internal = true;
      int         flag     = 0;
      ElementType etype    = get_elem_type(elemData.topology, m_mesh->get_spatial_dimension());
      int         dim1     = elemData.topology->spatial_dimension();

      /* need to check for hex's or tet's */
      bool hflag1 = is_hex(elemData.topology);

      /* a TET10 cannot connect to a HEX */
      bool tflag1 = is_tet(elemData.topology);

      /* check each side of this element */
      for (int nscnt = 0; nscnt < nsides; nscnt++) {

        /* get the list of nodes on this element side */
        get_local_side_nodes(elemData, nscnt + 1 /* 1 based */, side_nodes);
        int side_cnt = side_nodes.size();

        /*
         * now determine how many side set nodes are needed to
         * determine if there is an element connected to this side.
         *
         * 2-D - need two nodes, so find one intersection
         * 3-D - need three nodes, so find two intersections
         * NOTE: must check to make sure that this number is not
         *       larger than the number of nodes on the sides (ie - SHELL).
         */

        int nnodes = mesh_dim;
        if (side_cnt < nnodes) {
          nnodes = side_cnt;
        }
        nnodes--; /* decrement to find the number of intersections needed */

        int64_t nelem = 0; /* reset this in case no intersections are needed */

        /*
         * need to handle hex's differently because of
         * the tet/hex combination
         */

        if (!hflag1) { /* Not a hex */

          /* ignore degenerate bars */

          if (!((elemData.topology->is_alias("bar2") ||
                 elemData.topology->is_alias("shell_line_2")) &&
                side_nodes[0] == side_nodes[1])) {

            size_t nhold = sur_elem[side_nodes[0]].size();
            for (size_t ncnt = 0; ncnt < nhold; ncnt++) {
              hold_elem[ncnt] = sur_elem[side_nodes[0]][ncnt];
            }

            for (int ncnt = 0; ncnt < nnodes; ncnt++) {
              /* Find elements connected to both node '0' and node 'ncnt+1' */
              nelem =
                  find_intersection(Data(hold_elem), Data(sur_elem[side_nodes[(ncnt + 1)]]), nhold,
                                    sur_elem[side_nodes[(ncnt + 1)]].size(), Data(pt_list));

              if (nelem < 2) {
                break;
              }

              nhold = nelem;
              for (int ncnt2 = 0; ncnt2 < nelem; ncnt2++) {
                hold_elem[ncnt2] = hold_elem[pt_list[ncnt2]];
              }
            }
          }
          else {
            fmt::print("WARNING: Element = {} is a DEGENERATE BAR\n", ecnt + 1);
          }
        }
        else { /* Is a hex */
          /*
           * Hex faces are fairly complicated now. There are two
           * exceptions to the standard case:
           *   1. it is a degenerate hex (mimics a wedge). this has
           *      two special faces, the first is a triangle, and the
           *      second is a 2d line
           *   2. two tets are connected to this hex face
           */

          /* first need to check for a degenerate element */
          int dflag = 0;
          if (side_nodes[0] == side_nodes[1] || side_nodes[0] == side_nodes[3]) {
            dflag++;
          }
          if (side_nodes[2] == side_nodes[1] || side_nodes[2] == side_nodes[3]) {
            dflag++;
          }

          /*
           * if both flags are set, then this face is the 2d line,
           * and should be ignored with respect to elemental
           * communication maps
           */
          if (dflag == 2) {
            nelem = 1;
          }
          else {
            /*
             * In order to check for two tets connected to this face,
             * check the intersection of opposite corners of this face.
             * Both tets should show up in the intersection of one of the
             * sets of opposite corners (nothing should show up in the
             * other).
             */

            /*
             * Initial check is side nodes 0 and 2 which are
             * diagonally opposite
             */
            int    inode = 0;
            int    node  = 2;
            size_t nhold = 0;
            for (int ncnt = 0; ncnt < nnodes; ncnt++) {
              /* Find elements connected to both node 'inode' and node 'node' */
              nelem = find_intersection(Data(sur_elem[side_nodes[inode]]),
                                        Data(sur_elem[side_nodes[node]]),
                                        sur_elem[side_nodes[inode]].size(),
                                        sur_elem[side_nodes[node]].size(), Data(pt_list));

              if (nelem > 1) {
                if (ncnt == 0) {
                  nhold = nelem;
                  for (int ncnt2 = 0; ncnt2 < nelem; ncnt2++) {
                    hold_elem[ncnt2] = sur_elem[side_nodes[inode]][pt_list[ncnt2]];
                  }

                  if (dflag) {
                    /*
                     * in this case, need to get an intersection with
                     * another (unique) point since nodes 0 and 2
                     * may represent an edge and not the diagonal
                     */
                    if (side_nodes[1] != side_nodes[0] && side_nodes[1] != side_nodes[2]) {
                      node = 1;
                    }
                    else {
                      node = 3;
                    }
                  }
                  else {
                    /*
                     * in the non-degenerate case, if an element is connected
                     * to two opposite nodes, then it must share a face.
                     */
                    break;
                  }
                }
                else {
                  /* This is the second or later time through this
                         loop and each time through, there have been two
                         or more elements that are connected to 'node'
                         (which changes) and 'inode'.  We want to make
                         sure that the elements matched this time through
                         were also in the list the first time through so
                         that the elements contain nodes 0 1 2 of the face
                         and not just 0 1 and 0 2...
                         So, this time, only put an element in the list if
                         it was in the list before.
                   */
                  for (size_t ncnt2 = 0; ncnt2 < nhold; ncnt2++) {
                    hold_elem[ncnt2] = -hold_elem[ncnt2];
                  }
                  for (int ncnt3 = 0; ncnt3 < nelem; ncnt3++) {
                    for (size_t ncnt2 = 0; ncnt2 < nhold; ncnt2++) {
                      if (-hold_elem[ncnt2] == sur_elem[side_nodes[inode]][pt_list[ncnt3]]) {
                        hold_elem[ncnt2] = sur_elem[side_nodes[inode]][pt_list[ncnt3]];
                        break;
                      }
                    }
                  }
                  /* Now, go through list and cull out element < 0 */
                  size_t ncnt3 = 0;
                  for (size_t ncnt2 = 0; ncnt2 < nhold; ncnt2++) {
                    if (hold_elem[ncnt2] >= 0) {
                      hold_elem[ncnt3] = hold_elem[ncnt2];
                      ncnt3++;
                    }
                  }
                  nelem = ncnt3;
                  if (!dflag && nelem > 2) {
                    fmt::print(
                        stderr,
                        "Possible corrupted mesh detected at element {}, strange connectivity.\n",
                        ecnt);
                  }
                }
              }
              else { /* nelem == 1 or 0 */
                if (!dflag) {
                  nhold = sur_elem[side_nodes[1]].size();
                  for (size_t ncnt2 = 0; ncnt2 < nhold; ncnt2++) {
                    hold_elem[ncnt2] = sur_elem[side_nodes[1]][ncnt2];
                  }
                }
                inode = 1;
                node  = 3; /* The node diagonally opposite node 1 */
              }
            }
          }
        } /* "if (!hflag1)" */

        /*
         * if there is an element on this side of ecnt, then there
         * will be at least two elements in the intersection (one
         * will be ecnt)
         */
        if (nelem > 1) {

          /*
           * now go through and check each element in the list to see
           * if it on a different processor than ecnt.  Don't need to
           * worry about ecnt (which is in the list) since it is on
           * the same processor as itself.  Note that due to filtering
           * done above, we are guaranteed to either have an element
           * on a different processor or elem==ecnt.
           */
          for (int ncnt = 0; ncnt < nelem; ncnt++) {

            INT elem  = hold_elem[ncnt];
            int proc2 = get_processor(elem);
            assert(proc2 < m_numProcs);

            if (proc != proc2) {
              IossElementData elemData2 = m_mesh->get_local_element(elem);

              int dim2 = elemData2.topology->spatial_dimension();
              int diff = abs(dim1 - dim2);

              /*
               * hex's to shells - ok
               * shells to bar - ok
               * hex to bar - BAD since a BAR will see a HEX but a HEX will not
               *              see a BAR
               */

              if (diff < 2) {

                /* need to check for hex's */
                bool hflag2 = is_hex(elemData2.topology);
                bool tflag2 = is_tet(elemData2.topology);

                /* check here for tet/hex combinations */
                int sid = 0;
                if ((tflag1 && hflag2) || (hflag1 && tflag2)) {
                  /*
                   * have to call a special function to get the side id
                   * in these cases. In both cases, the number of side
                   * nodes for the element will not be consistent with
                   * side_cnt, and:
                   *
                   * TET/HEX - side_nodes only contains three of the
                   *           the side nodes of the hex.
                   *
                   * HEX/TET - Have to check that this tet shares a side
                   *           with the hex.
                   */
                  sid = get_side_id_hex_tet(elemData2.topology, Data(elemData2.localConnectivity),
                                            side_cnt, Data(side_nodes));
                }
                else {
                  /*
                   * get the side id of elem. Make sure that ecnt is
                   * trying to communicate to a valid side of elem
                   */
                  std::vector<utest_util::EntityId> mirror_nodes;

                  get_local_side_mirror_nodes(elemData, nscnt + 1 /* 1 based */, mirror_nodes);
                  side_cnt = side_nodes.size();

                  /*
                   * small kludge to handle 6 node faces butted up against
                   * 4 node faces
                   */
                  if (etype == ElementType::HEXSHELL && side_cnt == 6) {
                    side_cnt = 4;
                  }

                  /*
                   * in order to get the correct side order for elem,
                   * get the mirror of the side of ecnt
                   */
                  sid = get_side_id(elemData2, m_mesh->get_spatial_dimension(), side_cnt,
                                    Data(mirror_nodes));
                }

                if (sid > 0) {
                  /* Element is a border element */
                  internal = false;
                  if (!flag) {
                    flag = 1;
                    bor_elems[proc].push_back(ecnt);
                  }

                  /* now put ecnt into proc2's communications map */
                  e_cmap_elems[proc2].push_back(elem);
                  e_cmap_sides[proc2].push_back(sid);
                  e_cmap_procs[proc2].push_back(proc);
                  e_cmap_neigh[proc2].push_back(ecnt);
                }
                else if (sid < 0) {
                  /*
                   * too many errors with bad meshes, print out
                   * more information here for diagnostics
                   */
                  std::string cmesg =
                      fmt::format("Error returned while getting side id for communication map.\n");
                  cmesg += fmt::format("Element 1: {}\n", (ecnt + 1));

                  nnodes = elemData.topology->number_nodes();
                  cmesg += "connect table:";
                  for (int i = 0; i < nnodes; i++) {
                    std::string tmpstr =
                        fmt::format(" {}", (size_t)(elemData.localConnectivity[i] + 1));
                    cmesg += tmpstr;
                  }
                  cmesg += "\n";
                  cmesg += fmt::format("side id: {}\n", static_cast<size_t>(nscnt + 1));
                  cmesg += "side nodes:";
                  for (int i = 0; i < side_cnt; i++) {
                    std::string tmpstr = fmt::format(" {}", (size_t)(side_nodes[i] + 1));
                    cmesg += tmpstr;
                  }
                  cmesg += "\n";
                  cmesg += fmt::format("Element 2: {}\n", (size_t)(elem + 1));

                  nnodes = elemData2.topology->number_nodes();
                  cmesg += "connect table:";
                  for (int i = 0; i < nnodes; i++) {
                    std::string tmpstr =
                        fmt::format(" {}", (size_t)(elemData2.localConnectivity[i] + 1));
                    cmesg += tmpstr;
                  }
                  cmesg += "\n";
                  EXPECT_TRUE(false) << cmesg; /* and get out of here */

                } /* End "if sid < 0" */
              } /* End "if (sid > 0)" */
            } /* End "if (proc != proc2)" */
          } /* End "for (ncnt = 0; ncnt < nelem; ncnt++)" */
        } /* End "if (nelem > 1)" */
      } /* End "for (nscnt = 0; nscnt < nsides; nscnt++)" */

      if (internal) {
        int_elems[proc].push_back(ecnt);
      }
    } /* End "for(ecnt=0; ecnt < m_mesh->num_elems; ecnt++)" */
  }

  template <typename INT>
  void ElementPartition<INT>::categorize_nodes(const std::vector<std::vector<INT>> &sur_elem,
                                               const int                            max_nsur)
  {
    size_t numNodes = m_mesh->get_num_local_nodes();

    for (size_t ncnt = 0; ncnt < numNodes; ncnt++) {
      bool internal = true;
      int  proc     = 0;

      /* If a node is not connected to any elements (nsur_elem[ncnt] == 0),
         then it will be assigned to processor 0 and treated as internal.
      */

      if (!sur_elem[ncnt].empty()) {
        size_t elem = sur_elem[ncnt][0];
        proc        = get_processor(elem);
        assert(proc < m_numProcs);
        int flag = 0;
        for (size_t ecnt = 1; ecnt < sur_elem[ncnt].size(); ecnt++) {
          int proc2 = get_processor(sur_elem[ncnt][ecnt]);
          assert(proc2 < m_numProcs);
          /* check if the processor for any two surrounding elems differ */
          if (proc != proc2) {
            /* ncnt is a border node  of proc */
            internal = false;
            /* first, I have to deal with node being border for proc */
            if (!flag) {
              flag = 1; /* only want to do this once */
              bor_nodes[proc].push_back(ncnt);
            }

            /*
             * now I have to put ncnt in the border list for proc2
             * I need to check to make sure that this node has not
             * already been added to this list. If it has, then it
             * is in the last position in the array
             */
            if ((bor_nodes[proc2].empty()) || ((INT)ncnt != bor_nodes[proc2].back())) {
              bor_nodes[proc2].push_back(ncnt);
            }
          } /* if (proc != proc2) */
        } /* for(ecnt=1; ecnt < nsur_elem[ncnt]; ecnt++) */
      } /* if(nsur_elem[ncnt]) */

      if (internal) {
        /*
         * NOTE: if all of the processors above were the same, then
         * the one held in proc is the correct one
         */
        int_nodes[proc].push_back(ncnt);
      }
    } /* for(ncnt=0; ncnt < machine->num_nodes; ncnt++) */
  }

  template <typename INT>
  void ElementPartition<INT>::assign_border_node_processors(
      const std::vector<std::vector<INT>> &sur_elem, const int max_nsur)
  {
    size_t numNodes = m_mesh->get_num_local_nodes();

    /* Allocate memory for the border node processor IDs */
    for (int proc = 0; proc < m_numProcs; proc++) {
      if (!bor_nodes[proc].empty()) {
        born_procs[proc].resize(bor_nodes[proc].size());
      }
    }

    /* Now find the processor(s) associated with each border node */
    for (int pcnt = 0; pcnt < m_numProcs; pcnt++) {
      for (size_t ncnt = 0; ncnt < bor_nodes[pcnt].size(); ncnt++) {
        size_t node = bor_nodes[pcnt][ncnt];

        for (size_t ecnt = 0; ecnt < sur_elem[node].size(); ecnt++) {
          size_t elem = sur_elem[node][ecnt];
          int    proc = get_processor(elem);
          assert(proc < m_numProcs);
          if (proc != pcnt) {
            if (in_list(proc, born_procs[pcnt][ncnt]) < 0) {
              born_procs[pcnt][ncnt].push_back(proc);
            }
          } /* End "if(proc != pcnt)" */
        } /* End "for(ecnt=0; ecnt < graph->nsur_elems[node]; ecnt++)" */
      } /* End "for(ncnt=0; ncnt < num_bor_nodes[pcnt]; ncnt++)" */
    } /* End "for(pcnt=0; pcnt < machine->num_procs; pcnt++)" */
  }

  template <typename INT> void ElementPartition<INT>::order_element_communication_maps()
  {
    /* Order the element communication maps by processor */
    for (int pcnt = 0; pcnt < m_numProcs; pcnt++) {
      /* Note that this sort is multi-key */
      qsort4(Data(e_cmap_procs[pcnt]),   /* 1st key */
             Data(e_cmap_elems[pcnt]),   /* 2nd key */
             Data(e_cmap_neigh[pcnt]),   /* 3rd key */
             Data(e_cmap_sides[pcnt]),   /* 4th key */
             e_cmap_procs[pcnt].size()); /* Size */
    }
    /*
     * At this point, each processors arrays are sorted on three keys:
     * [processor, element, neighbor]
     */

    /*
     * Now order the elemental communication maps so that they are
     * consistent between processors.
     */
    for (INT pcnt = 1; pcnt < m_numProcs; pcnt++) {
      int  save_fv1 = 0;
      auto size = static_cast<INT>(e_cmap_procs[pcnt].size()); /* Define shortcuts size and procs */
      std::vector<INT> procs = e_cmap_procs[pcnt];

      INT fv1 = -1;
      INT lv1 = -1;

      for (int pcnt2 = 0; pcnt2 < pcnt; pcnt2++) {

        /*
         * Find the first and last entries for processor "pcnt2" in
         * the list of processor "pcnt".
         */

        /* c_cmap_procs[pcnt] is sorted based on processor.
         * From point 'save_fv1' search for 'pcnt2'
         * If not found, value is -1; else search for !pcnt2 from
         * that point forward.
         */
        INT i = save_fv1;
        while (i < size && procs[i] < pcnt2) {
          i++;
        }
        if (i >= size || procs[i] != pcnt2) {
          fv1 = -1;
          lv1 = -1;
        }
        else {
          fv1 = i;
          assert(procs[i] == pcnt2);
          for (lv1 = fv1; lv1 < size; lv1++) {
            if (procs[lv1] != pcnt2) {
              lv1 = lv1 - 1;
              assert(procs[lv1] == pcnt2);
              break;
            }
          }
        }

        if (lv1 >= size) {
          lv1 = size - 1;
        }

        if (lv1 != -1) {
          save_fv1 = lv1 + 1;
        }

#if 0
        /* Old method -- can use for verification by uncommenting this if block  */
        {
          int tst_fv1, tst_lv1;
          find_first_last(pcnt2, size, procs, &tst_fv1, &tst_lv1);
          assert(tst_fv1 == fv1);
          assert(tst_lv1 == lv1);
        }
#endif

        if (fv1 >= 0) {
          /* Sort based on neighbor element */
          sort3(lv1 - fv1 + 1, (&e_cmap_neigh[pcnt][fv1]), (&e_cmap_elems[pcnt][fv1]),
                (&e_cmap_sides[pcnt][fv1]));
          /*
           * Find the first and last entries for processor "pcnt" in
           * the list of processor "pcnt2".
           */
          INT fv2 = -1;
          INT lv2 = -1;
          find_first_last(pcnt, e_cmap_procs[pcnt2].size(), &e_cmap_procs[pcnt2][0], &fv2, &lv2);
#if 1
          if (lv2 - fv2 != lv1 - fv1) {
            fmt::print(stderr, "{}: {} to {}\n", static_cast<size_t>(pcnt2), (size_t)fv1,
                       (size_t)lv1);
            for (i = fv1; i <= lv1; i++) {
              fmt::print(stderr, "{}: {}\t{}\t{}\t{}\n", static_cast<size_t>(i),
                         (size_t)e_cmap_elems[pcnt][i], (size_t)e_cmap_neigh[pcnt][i],
                         (size_t)e_cmap_procs[pcnt][i], (size_t)e_cmap_sides[pcnt][i]);
            }
            fmt::print(stderr, "{}: {} to {}\n", (size_t)pcnt, (size_t)fv2, (size_t)lv2);
            for (i = fv2; i <= lv2; i++) {
              fmt::print(stderr, "{}: {}\t{}\t{}\t{}\n", static_cast<size_t>(i),
                         (size_t)e_cmap_elems[pcnt2][i], (size_t)e_cmap_neigh[pcnt2][i],
                         (size_t)e_cmap_procs[pcnt2][i], (size_t)e_cmap_sides[pcnt2][i]);
            }
          }
#endif
          assert(lv2 - fv2 == lv1 - fv1);

          /* Sort based on element -- This will then match order of
           * the fv1->lv1 arrays.
           */
          sort3(lv2 - fv2 + 1, (&e_cmap_elems[pcnt2][fv2]), (&e_cmap_neigh[pcnt2][fv2]),
                (&e_cmap_sides[pcnt2][fv2]));

        } /* End "if(fv1 >= 0)" */
      } /* End "for(pcnt2=0; pcnt2 < pcnt; pcnt2++)" */
    } /* End "for(pcnt=0; pcnt < machine->num_procs; pcnt++)" */
  }

  template <typename INT> void ElementPartition<INT>::generate_partition_maps()
  {
    /* Allocate memory */
    int_nodes.resize(m_numProcs);
    bor_nodes.resize(m_numProcs);
    ext_nodes.resize(m_numProcs);

    int_elems.resize(m_numProcs);
    bor_elems.resize(m_numProcs);

    ext_procs.resize(m_numProcs);
    born_procs.resize(m_numProcs);

    e_cmap_elems.resize(m_numProcs);
    e_cmap_sides.resize(m_numProcs);
    e_cmap_procs.resize(m_numProcs);
    e_cmap_neigh.resize(m_numProcs);

    int                           max_nsur{0};
    std::vector<std::vector<INT>> sur_elem;

    fill_node_element_connectivity(m_mesh, sur_elem, max_nsur);

    categorize_elements(sur_elem, max_nsur);

    categorize_nodes(sur_elem, max_nsur);

    assign_border_node_processors(sur_elem, max_nsur);

    order_element_communication_maps();

    for (int proc = 0; proc < m_numProcs; proc++) {
      /* Sort node maps */
      gds_qsort(Data(int_nodes[proc]), int_nodes[proc].size());

      /* Sort element maps */
      gds_qsort(Data(int_elems[proc]), int_elems[proc].size());
    }
  }

  template void ElementPartition<int>::write_nemesis_data(int exoid) const;
  template void ElementPartition<int64_t>::write_nemesis_data(int exoid) const;

  template <typename INT> void ElementPartition<INT>::write_nemesis_data(int exoid) const
  {
    int    status        = 0;
    size_t numNodes      = m_mesh->get_num_local_nodes();
    size_t numElems      = m_mesh->get_num_local_elements();
    size_t numElemBlocks = m_mesh->get_num_global_element_blocks();

    /* Output the the initial Nemesis global information */
    status = ex_put_init_global(exoid, numNodes, numElems, numElemBlocks, 0, 0);
    ASSERT_FALSE(status < 0) << "fatal: failed to output initial Nemesis parameters";

    {
      auto ebs = m_mesh->get_region()->get_element_blocks();

      std::vector<INT> eb_ids(ebs.size(), 0);
      std::vector<INT> eb_cnts(ebs.size(), 0);

      for (size_t i = 0; i < ebs.size(); i++) {
        eb_ids[i]  = ebs[i]->get_property("id").get_int();
        eb_cnts[i] = ebs[i]->entity_count();
      }
      ex_put_eb_info_global(exoid, Data(eb_ids), Data(eb_cnts));
    }

    /* Set up dummy arrays for output */
    std::vector<INT> num_nmap_cnts(m_numProcs, 0);
    std::vector<INT> num_emap_cnts(m_numProcs, 0);

    /* need to check and make sure that there really are comm maps */
    for (int cnt = 0; cnt < m_numProcs; cnt++) {
      if (!bor_nodes[cnt].empty()) {
        num_nmap_cnts[cnt] = 1;
      }
    }
    for (int cnt = 0; cnt < m_numProcs; cnt++) {
      if (!bor_elems[cnt].empty()) {
        num_emap_cnts[cnt] = 1;
      }
    }

    status = ex_put_init_info(exoid, m_numProcs, m_numProcs, (char *)"s");
    ASSERT_FALSE(status < 0) << "fatal: unable to output init info";

    // Need to create 5 arrays with the sizes of int_nodes[i].size()...
    {
      std::vector<INT> ins(m_numProcs, 0);
      std::vector<INT> bns(m_numProcs, 0);
      std::vector<INT> ens(m_numProcs, 0);
      std::vector<INT> ies(m_numProcs, 0);
      std::vector<INT> bes(m_numProcs, 0);

      for (int iproc = 0; iproc < m_numProcs; iproc++) {
        ins[iproc] = int_nodes[iproc].size();
        bns[iproc] = bor_nodes[iproc].size();
        ens[iproc] = ext_nodes[iproc].size();
        ies[iproc] = int_elems[iproc].size();
        bes[iproc] = bor_elems[iproc].size();
      }

      status = ex_put_loadbal_param_cc(exoid, Data(ins), Data(bns), Data(ens), Data(ies), Data(bes),
                                       Data(num_nmap_cnts), Data(num_emap_cnts));
      ASSERT_FALSE(status < 0) << "fatal: unable to output load-balance parameters";
    }

    std::vector<INT> node_proc_ptr(m_numProcs + 1, 0);
    std::vector<INT> node_cmap_ids_cc(m_numProcs, 0);
    std::vector<INT> node_cmap_cnts_cc(m_numProcs, 0);

    node_proc_ptr[0] = 0;
    for (int proc = 0; proc < m_numProcs; proc++) {
      node_proc_ptr[proc + 1] = node_proc_ptr[proc] + 1;

      node_cmap_cnts_cc[proc] = 0;
      for (size_t cnt = 0; cnt < bor_nodes[proc].size(); cnt++) {
        node_cmap_cnts_cc[proc] += born_procs[proc][cnt].size();
      }

      node_cmap_ids_cc[proc] = 1;
    }

    std::vector<INT> elem_proc_ptr(m_numProcs + 1, 0);
    std::vector<INT> elem_cmap_ids_cc(m_numProcs, 0);
    std::vector<INT> elem_cmap_cnts_cc(m_numProcs, 0);

    elem_proc_ptr[0] = 0;
    for (int proc = 0; proc < m_numProcs; proc++) {
      elem_proc_ptr[proc + 1] = elem_proc_ptr[proc] + 1;
      elem_cmap_cnts_cc[proc] = e_cmap_elems[proc].size();
      elem_cmap_ids_cc[proc]  = 1;
    }

    /* Output the communication map parameters */
    status = ex_put_cmap_params_cc(exoid, Data(node_cmap_ids_cc), Data(node_cmap_cnts_cc),
                                   Data(node_proc_ptr), Data(elem_cmap_ids_cc),
                                   Data(elem_cmap_cnts_cc), Data(elem_proc_ptr));
    ASSERT_FALSE(status < 0) << "fatal: unable to output communication map parameters";

    /* Output the node and element maps */
    for (int proc = 0; proc < m_numProcs; proc++) {
      /* Output the nodal map */
      status = ex_put_processor_node_maps(exoid, Data(int_nodes[proc]), Data(bor_nodes[proc]),
                                          nullptr, proc);
      ASSERT_FALSE(status < 0) << "fatal: failed to output node map for proc: " << proc;

      /* Output the elemental map */
      status =
          ex_put_processor_elem_maps(exoid, Data(int_elems[proc]), Data(bor_elems[proc]), proc);
      ASSERT_FALSE(status < 0) << "fatal: failed to output element map for proc: " << proc;

      /*
       * Build a nodal communication map from the list of border nodes
       * and their associated processors and side IDs.
       */
      size_t nsize = 0;
      for (size_t cnt = 0; cnt < bor_nodes[proc].size(); cnt++) {
        nsize += born_procs[proc][cnt].size();
      }

      if (nsize > 0) {
        std::vector<INT> n_cmap_nodes(nsize, 0);
        std::vector<INT> n_cmap_procs(nsize, 0);

        size_t cnt3 = 0;
        for (size_t cnt = 0; cnt < bor_nodes[proc].size(); cnt++) {
          for (size_t cnt2 = 0; cnt2 < born_procs[proc][cnt].size(); cnt2++) {
            n_cmap_nodes[cnt3]   = bor_nodes[proc][cnt];
            n_cmap_procs[cnt3++] = born_procs[proc][cnt][cnt2];
          }
        }

        /*
         * Reorder the nodal communication maps so that they are ordered
         * by processor and then by global ID.
         */
        /* This is a 2-key sort */
        qsort2(Data(n_cmap_procs), Data(n_cmap_nodes), cnt3);

        /* Output the nodal communication map */
        status = ex_put_node_cmap(exoid, 1, Data(n_cmap_nodes), Data(n_cmap_procs), proc);
        ASSERT_FALSE(status < 0) << "fatal: unable to output nodal communication map for proc: "
                                 << proc;
      } /* End "if (nsize > 0)" */

      /* Output the elemental communication map */
      if (!e_cmap_elems[proc].empty()) {
        status = ex_put_elem_cmap(exoid, 1, Data(e_cmap_elems[proc]), Data(e_cmap_sides[proc]),
                                  Data(e_cmap_procs[proc]), proc);
        ASSERT_FALSE(status < 0) << "fatal: unable to output elemental communication map for proc: "
                                 << proc;
      }

    } /* End "for(proc=0; proc < m_numProcs; proc++)" */
  }

} // namespace utest_util
