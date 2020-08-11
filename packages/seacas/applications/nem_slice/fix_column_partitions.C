/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "elb.h" // for LB_Description<INT>, etc
#include "elb_elem.h"
#include "elb_err.h"
#include "fix_column_partitions.h"
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

namespace {
  // Opposite side IDs in a hex according to Exodus II convention

  int hex_opp_side[6] = {3, 4, 1, 2, 6, 5};

  /*! @brief Given an element and a side, find the adjacent element to that side
    @param cur_elem  Current element under consideration
    @param etype     Element type
    @param side_id   Side across which we want to find an adjacent element
    @param nadj      Number of elements adjacent to cur_elem (from graph description)
    @param adj       Pointer to elements adjacent to cur_elem (from graph description)
    @param global_connect Global connectivity array
    @param adj_elem  ID of adjacent element (-1 if not found)
    @param adj_side  Local ID of common side in adjacent element (0 if adj_elem not found)
  */

  template <typename INT>
  void find_adjacent_element(INT cur_elem, E_Type etype, int side_id, int nadj, INT const *adj,
                             Mesh_Description<INT> const *const mesh, INT *adj_elem, int *adj_side)
  {
    *adj_elem = -1;

    // Get nodes of this side (face) of the element

    int nsnodes = is_hex(etype) ? 4 : 3;
    INT side_nodes[9]; // SHELL9 has 9 nodes on a face.

    INT *elnodes = mesh->connect[cur_elem];
    ss_to_node_list(etype, elnodes, side_id, side_nodes);

    // How would these side's nodes be if they were viewed from the
    // adjacent element

    INT side_nodes_flipped[9]; // Realistically: 4, max possible: 9
    get_ss_mirror(etype, side_nodes, side_id, side_nodes_flipped);

    for (int i = 0; i < nadj; i++) {
      INT    adj_elem_id = adj[i] - 1; // Adjacency graph entries start from 1
      E_Type etype2      = mesh->elem_type[adj_elem_id];

      // Does this side occurs in the adjacent element?

      INT *elnodes2 = mesh->connect[adj_elem_id];

      // options to keep 'get_side_id' to not flip out if side nodes are
      // not found in element

      int skip_check  = 2;
      int partial_adj = 1;
      *adj_side =
          get_side_id(etype2, elnodes2, nsnodes, side_nodes_flipped, skip_check, partial_adj);

      if (*adj_side > 0) {
        *adj_elem = adj_elem_id;
        return;
      }
    }
  }
} // namespace

/*! @brief If the mesh is columnar, ensure each column is fully in on partition
  @param lb    Load balancing or partitioning information (may be modified)
  @param mesh  Description of the mesh
  @param graph Description of the adjacency graph

  **** ASSUMES COLUMNS ARE STRICTLY VERTICAL, i.e., NO LATERAL FACE ****
  **** HAS A Z-COMPONENT IN ITS NORMAL ****

*/

template int fix_column_partitions(LB_Description<int> *lb, Mesh_Description<int> const *const mesh,
                                   Graph_Description<int> const *const graph);
template int fix_column_partitions(LB_Description<int64_t> *               lb,
                                   Mesh_Description<int64_t> const *const  mesh,
                                   Graph_Description<int64_t> const *const graph);

template <typename INT>
int fix_column_partitions(LB_Description<INT> *lb, Mesh_Description<INT> const *const mesh,
                          Graph_Description<INT> const *const graph)
{
  int nmoved = 0;
  int nel    = mesh->num_elems;
  int nnod   = mesh->num_nodes;

  // a flag to indicate if a particular element has been processed
  std::vector<bool> processed_flag(nel, false);

  // Go through elements and attempt to discover a column of elements
  // that contain it. Then check if the column is all on one partition
  // - if not, fix it

  for (int i = 0; i < nel; i++) {
    if (processed_flag[i])
      continue;

    E_Type etype = mesh->elem_type[i];

    // Only hexes and wedges can be stacked in columns
    if (!is_hex(etype) && !is_wedge(etype)) {
      continue;
    }

    // retrieve the faces of this element - faces are described by the
    // local numbering of nodes with respect to the element node list

    INT *elnodes  = mesh->connect[i];
    int  nelnodes = get_elem_info(NNODES, etype);

    float elcoord[27][3];
    for (int j = 0; j < nelnodes; j++) {
      for (int d = 0; d < 3; d++) {
        elcoord[j][d] = mesh->coords[elnodes[j] + d * nnod];
      }
    }

    int top_side0 = 0;
    int bot_side0 = 0;
    int nelfaces  = get_elem_info(NSIDES, etype);

    // Find top and bottom faces by eliminating lateral faces under
    // the assumption that lateral face normals have no Z component

    int count = 0;
    for (int j = 0; j < nelfaces; j++) {
      INT fnodes[9]; // Should only need 4, but ss_to_node_list can potentially access 9 (SHELL9).

      int nfn = 4;
      if (is_wedge(etype)) {
        if (j < 3) {
          // lateral faces of wedge according to Exodus II - cannot be
          // up/down faces in a column
          continue;
        }
        nfn = 3;
      }

      // Nodes of the side/face
      ss_to_node_list(etype, mesh->connect[i], j + 1, fnodes);

      // Translate global IDs of side nodes to local IDs in element
      int fnodes_loc[9];
      for (int k = 0; k < nfn; k++) {
        bool found = false;
        for (int k2 = 0; k2 < nelnodes; k2++) {
          if (fnodes[k] == elnodes[k2]) {
            fnodes_loc[k] = k2;
            found         = true;
            break;
          }
        }
        if (!found)
          Gen_Error(0, "FATAL: side/face node not found in element node list?");
      }

      double normal[3] = {0.0, 0.0, 0.0};
      if (nfn == 3) {
        double v0[3];
        double v1[3];
        for (int d = 0; d < 3; d++) {
          v0[d] = elcoord[fnodes_loc[1]][d] - elcoord[fnodes_loc[0]][d];
          v1[d] = elcoord[fnodes_loc[nfn - 1]][d] - elcoord[fnodes_loc[0]][d];
        }

        // cross product to get normal corner
        for (int d = 0; d < 3; d++) {
          normal[d] = v0[(d + 1) % 3] * v1[(d + 2) % 3] - v0[(d + 2) % 3] * v1[(d + 1) % 3];
        }
      }
      else {
        for (int k = 0; k < nfn; k++) {
          double v0[3];
          double v1[3];
          for (int d = 0; d < 3; d++) {
            v0[d] = elcoord[fnodes_loc[(k + 1) % nfn]][d] - elcoord[fnodes_loc[k]][d];
            v1[d] = elcoord[fnodes_loc[(k - 1 + nfn) % nfn]][d] - elcoord[fnodes_loc[k]][d];
          }

          // cross product to get normal at corner - add to face normal
          for (int d = 0; d < 3; d++) {
            normal[d] += v0[(d + 1) % 3] * v1[(d + 2) % 3] - v0[(d + 2) % 3] * v1[(d + 1) % 3];
          }
        }
        double len = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
        if (len > 1.0e-24) { // Don't normalize nearly zero-length vectors
          for (double &d : normal) {
            d /= len;
          }
        }
      }
      if (fabs(normal[2]) > 1.0e-12) { // non-zero
        if (normal[2] > 0.0) {
          top_side0 = j + 1; // side id counting starts from 1
        }
        else {
          bot_side0 = j + 1;
        }
        count++;
      }
    } // for (j = 0; j < nelfaces; j++)
    if (count > 2) {
#ifdef DEBUG
      Gen_Error(1, "WARNING: more than two faces with non-zero Z-components of normal.");
#endif
      Gen_Error(1, "WARNING: Mesh may not be strictly columnar. Initial partitioning unchanged.");
      return 0;
    }
    if (count < 2) {
#ifdef DEBUG
      Gen_Error(1, "WARNING: could not find up and down faces for element.");
#endif
      Gen_Error(1, "WARNING: Mesh may not be strictly columnar. Initial partitioning unchanged.");
      return 0;
    }

    // Found top and bottom sides/faces of current element. Build a
    // lists of elements stacked above it and below it.

    std::vector<INT> above_list;
    std::vector<INT> below_list;

    INT  cur_elem      = i;
    INT  adj_elem      = -1;
    int  adj_side      = -1;
    int  bot_side      = bot_side0;
    int  top_side      = top_side0;
    bool upsearch_done = false;
    while (!upsearch_done) {

      int        nadj = graph->start[cur_elem + 1] - graph->start[cur_elem];
      INT const *adj  = graph->adj.data() + graph->start[cur_elem];
      find_adjacent_element(cur_elem, etype, top_side, nadj, adj, mesh, &adj_elem, &adj_side);
      if (adj_elem == -1) {
        upsearch_done = true;
      }
      else {
        above_list.push_back(adj_elem);
        cur_elem = adj_elem;
        bot_side = adj_side;

        if (is_hex(etype)) {
          top_side = hex_opp_side[bot_side - 1];
        }
        else { // wedge
          if (bot_side != 4 && bot_side != 5) {
            Gen_Error(0, "FATAL: Expected bottom side in wedge to be side 4 or 5");
            return 0;
          }
          top_side = (bot_side == 4) ? 5 : 4;
        }
        if (processed_flag[adj_elem]) {
          Gen_Error(0, "FATAL: repeated column elements");
          return 0;
        }
        processed_flag[adj_elem] = true;
      }
    } // while (!upsearch_done)
    int nabove = above_list.size();

    cur_elem             = i;
    bot_side             = bot_side0;
    top_side             = top_side0;
    bool downsearch_done = false;
    while (!downsearch_done) {

      int        nadj = graph->start[cur_elem + 1] - graph->start[cur_elem];
      INT const *adj  = graph->adj.data() + graph->start[cur_elem];
      find_adjacent_element(cur_elem, etype, bot_side, nadj, adj, mesh, &adj_elem, &adj_side);
      if (adj_elem == -1) {
        downsearch_done = true;
      }
      else {
        below_list.push_back(adj_elem);
        cur_elem = adj_elem;
        top_side = adj_side;

        if (is_hex(etype)) {
          bot_side = hex_opp_side[top_side - 1];
        }
        else { // wedge
          if (top_side != 4 && top_side != 5) {
            Gen_Error(0, "FATAL: Expected top side in wedge to be side 4 or 5");
            return 0;
          }
          bot_side = (top_side == 4) ? 5 : 4;
        }
        if (processed_flag[adj_elem]) {
          Gen_Error(0, "FATAL: repeated column elements");
          return 0;
        }
        processed_flag[adj_elem] = true;
      }
    } // while (!upsearch_done)
    int nbelow = below_list.size();

    // Build list of elements in column from top to bottom
    // Code below assumes we are NOT compiling with C++11 standard enabled
    // If we do, we can simplify the code quite a bit

    std::vector<INT> colelems;
    colelems.reserve(nabove + nbelow + 1);
    typename std::vector<INT>::reverse_iterator rit = above_list.rbegin();
    while (rit != above_list.rend()) {
      colelems.push_back(*rit);
      ++rit;
    }
    colelems.push_back(i);

    typename std::vector<INT>::iterator it = below_list.begin();
    while (it != below_list.end()) {
      colelems.push_back(*it);
      ++it;
    }

    // Make all the other elements in the column be on the same
    // processor as the one that majority of the elements are on.
    // To do that make a unique list of processors
    std::map<int, int> procid_elcount;

    it = colelems.begin();
    while (it != colelems.end()) {
      INT elem2  = *it;
      int procid = lb->vertex2proc[elem2];

      // Try to insert 'procid' with element count of 1 as a new entry
      // If it fails, the processor is already in the map; then just
      // increment the element count

      std::pair<std::map<int, int>::iterator, bool> status =
          procid_elcount.insert(std::pair<int, int>(procid, 1));
      if (status.second == false) { // procid already in map; could not insert
        std::map<int, int>::iterator itmap = status.first;
        (itmap->second)++;
      }
      ++it;
    }

    // Which processor has a dominant presence in this column?
    int                          max_procid = -1;
    int                          max_elems  = 0;
    std::map<int, int>::iterator itmap      = procid_elcount.begin();
    while (itmap != procid_elcount.end()) {
      if (itmap->second > max_elems) {
        max_procid = itmap->first;
        max_elems  = itmap->second;
      }
      ++itmap;
    }

    // Switch all elements in the column to the dominant processor
    it = colelems.begin();
    while (it != colelems.end()) {
      INT elem2 = *it;
      if (lb->vertex2proc[elem2] != max_procid) {
#ifdef DEBUG
        fmt::print(" Reassigning element {} from proc {} to {}\n", elem2, lb->vertex2proc[elem2],
                   max_procid);
#endif
        lb->vertex2proc[elem2] = max_procid;
        ++nmoved;
      }
      ++it;
    }

  } // for (int i = 0; i < nel; i++)

  return nmoved;
} // fix_column_partitions.C
