// Copyright(C) 1999-2020, 2022, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifdef SEACAS_HAVE_MPI
#include "mpi.h"
#endif
#include "gtest/gtest.h"

#include "ElementUtils.h"

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

namespace utest_util {

  bool is_hex(const Ioss::ElementTopology *topo)
  {
    return topo->shape() == Ioss::ElementShape::HEX;
  }

  bool is_hex(ElementType etype)
  {
    return etype == ElementType::HEX8 || etype == ElementType::HEX27 ||
           etype == ElementType::HEX20 || etype == ElementType::HEXSHELL;
  }

  bool is_tet(const Ioss::ElementTopology *topo)
  {
    return topo->shape() == Ioss::ElementShape::TET;
  }

  bool is_tet(ElementType etype)
  {
    return etype == ElementType::TET4 || etype == ElementType::TET10 ||
           etype == ElementType::TET8 || etype == ElementType::TET14 || etype == ElementType::TET15;
  }

  bool is_wedge(const Ioss::ElementTopology *topo)
  {
    return topo->shape() == Ioss::ElementShape::WEDGE;
  }

  bool is_wedge(ElementType etype)
  {
    return etype == ElementType::WEDGE6 || etype == ElementType::WEDGE15 ||
           etype == ElementType::WEDGE16 || etype == ElementType::WEDGE20 ||
           etype == ElementType::WEDGE21;
  }

  bool is_pyramid(const Ioss::ElementTopology *topo)
  {
    return topo->shape() == Ioss::ElementShape::PYRAMID;
  }

  bool is_pyramid(ElementType etype)
  {
    return etype == ElementType::PYRAMID5 || etype == ElementType::PYRAMID13 ||
           etype == ElementType::PYRAMID14 || etype == ElementType::PYRAMID18 ||
           etype == ElementType::PYRAMID19;
  }

  bool is_3d_element(const Ioss::ElementTopology *topo)
  {
    return is_hex(topo) || is_tet(topo) || is_wedge(topo) || is_pyramid(topo);
  }

  bool is_3d_element(ElementType etype)
  {
    return is_hex(etype) || is_tet(etype) || is_wedge(etype) || is_pyramid(etype);
  }

  ElementType get_elem_type(const char *elem_name, const int num_nodes, const int num_dim)
  {
    ElementType answer = ElementType::NULL_EL;
    switch (elem_name[0]) {
    case 'h':
    case 'H':
      if (strncasecmp(elem_name, "HEX", 3) == 0) {
        switch (num_nodes) {
        case 8: answer = ElementType::HEX8; break;
        case 12: answer = ElementType::HEXSHELL; break;
        case 16: answer = ElementType::HEX16; break;
        case 20: answer = ElementType::HEX20; break;
        case 27: answer = ElementType::HEX27; break;
        default: EXPECT_TRUE(false) << "fatal: unsupported HEX element"; exit(1);
        }
      }
      break;

    case 'c':
    case 'C':
      if (strncasecmp(elem_name, "CIRCLE", 6) == 0) {
        answer = ElementType::SPHERE;
      }
      break;

    case 's':
    case 'S':
      if (strncasecmp(elem_name, "SPHERE", 6) == 0) {
        answer = ElementType::SPHERE;
      }
      else if (strncasecmp(elem_name, "SHELL", 5) == 0) {
        switch (num_nodes) {
        case 2:
          if (num_dim == 2) {
            answer = ElementType::SHELL2;
          }
          else {
            EXPECT_TRUE(false) << "fatal: unsupported SHELL element";
            exit(1);
          }
          break;
        case 3:
          if (num_dim == 2) {
            answer = ElementType::SHELL3;
          }
          else {
            EXPECT_TRUE(false) << "fatal: unsupported SHELL element";
            exit(1);
          }
          break;
        case 4: answer = ElementType::SHELL4; break;
        case 8: answer = ElementType::SHELL8; break;
        case 9: answer = ElementType::SHELL9; break;
        default: EXPECT_TRUE(false) << "fatal: unsupported SHELL element"; exit(1);
        }
      }
      break;

    case 'b':
    case 'B':
    case 't':
    case 'T':
    case 'r':
    case 'R':
      if (strncasecmp(elem_name, "BEAM", 4) == 0 || strncasecmp(elem_name, "TRUSS", 5) == 0 ||
          strncasecmp(elem_name, "ROD", 3) == 0 || strncasecmp(elem_name, "BAR", 3) == 0) {
        switch (num_nodes) {
        case 2: answer = num_dim == 1 ? ElementType::BAR1D2 : ElementType::BAR2; break;
        case 3: answer = num_dim == 1 ? ElementType::BAR1D3 : ElementType::BAR3; break;
        default: EXPECT_TRUE(false) << "fatal: unsupported BAR/BEAM/TRUSS element"; exit(1);
        }
      }
      else if (strncasecmp(elem_name, "TRI", 3) == 0) {
        switch (num_nodes) {
        case 3:
          if (num_dim == 2) {
            answer = ElementType::TRI3;
          }
          else {
            answer = ElementType::TSHELL3;
          }
          break;
        case 4:
          if (num_dim == 2) {
            answer = ElementType::TRI4;
          }
          else {
            answer = ElementType::TSHELL4;
          }
          break;
        case 6:
          if (num_dim == 2) {
            answer = ElementType::TRI6;
          }
          else {
            answer = ElementType::TSHELL6;
          }
          break;
        case 7:
          if (num_dim == 2) {
            answer = ElementType::TRI7;
          }
          else {
            answer = ElementType::TSHELL7;
          }
          break;
        default: EXPECT_TRUE(false) << "fatal: unsupported TRI element"; exit(1);
        }
      }
      else if (strncasecmp(elem_name, "TET", 3) == 0) {
        switch (num_nodes) {
        case 4: answer = ElementType::TET4; break;
        case 8: answer = ElementType::TET8; break;
        case 10: answer = ElementType::TET10; break;
        case 14: answer = ElementType::TET14; break;
        case 15: answer = ElementType::TET15; break;
        default: EXPECT_TRUE(false) << "fatal: unsupported TET element"; exit(1);
        }
      }
      break;

    case 'q':
    case 'Q':
      if (strncasecmp(elem_name, "QUAD", 4) == 0) {
        switch (num_nodes) {
        case 4:
          if (num_dim == 2) {
            answer = ElementType::QUAD4;
          }
          else {
            answer = ElementType::SHELL4;
          }
          break;
        case 8:
          if (num_dim == 2) {
            answer = ElementType::QUAD8;
          }
          else {
            answer = ElementType::SHELL8;
          }
          break;
        case 9:
          if (num_dim == 2) {
            answer = ElementType::QUAD9;
          }
          else {
            answer = ElementType::SHELL9;
          }
          break;
        default: EXPECT_TRUE(false) << "fatal: unsupported QUAD element"; exit(1);
        }
      }
      break;

    case 'w':
    case 'W':
      if (strncasecmp(elem_name, "WEDGE", 5) == 0) {
        switch (num_nodes) {
        case 6: answer = ElementType::WEDGE6; break;
        case 12: answer = ElementType::WEDGE12; break;
        case 15: answer = ElementType::WEDGE15; break;
        case 16: answer = ElementType::WEDGE16; break;
        case 20: answer = ElementType::WEDGE20; break;
        case 21: answer = ElementType::WEDGE21; break;
        default: EXPECT_TRUE(false) << "fatal: unsupported WEDGE element"; exit(1);
        }
      }
      break;

    case 'p':
    case 'P':
      if (strncasecmp(elem_name, "PYR", 3) == 0) {
        switch (num_nodes) {
        case 5: answer = ElementType::PYRAMID5; break;
        case 13: answer = ElementType::PYRAMID13; break;
        case 14: answer = ElementType::PYRAMID14; break;
        case 18: answer = ElementType::PYRAMID18; break;
        case 19: answer = ElementType::PYRAMID19; break;
        default: EXPECT_TRUE(false) << "fatal: unsupported PYRAMID element"; exit(1);
        }
      }
      break;

    default: break;
    }

    if (answer == ElementType::NULL_EL) {
      std::string errstr;
      errstr = fmt::format("fatal: unsupported element type '{}' read", elem_name);
      EXPECT_TRUE(false) << errstr;
      exit(1);
    }

    return answer;

  } /*---------------------------End get_elem_type()---------------------------*/

  ElementType get_elem_type(const Ioss::ElementTopology *topo, int num_dim)
  {
    return get_elem_type(topo->name().c_str(), topo->number_nodes(), num_dim);
  }

  /*
   *----------------------------------------------------------------------------
   * This function returns the list of nodes in a side of an element given
   * the element type, and the side id. It also returns the number of nodes
   * in that side.
   *****************************************************************************/
  void get_local_side_nodes(const IossElementData &elem, int side, std::vector<EntityId> &sideNodes)
  {
    Ioss::IntVector indices = elem.topology->boundary_connectivity(side);

    sideNodes.clear();
    for (auto i : indices) {
      sideNodes.push_back(elem.localConnectivity[i]);
    }
  }

  /*
   *----------------------------------------------------------------------------
   * This function returns the node list for the mirror of the list
   * given. This will be the node list of a face that is connected
   * to this element on this face i.e the first negative permutation
   *****************************************************************************/
  void get_local_side_mirror_nodes(const IossElementData &elem, int side,
                                   std::vector<EntityId> &mirrorNodes)
  {
    // Get positive permutation
    std::vector<EntityId> sideNodes;
    get_local_side_nodes(elem, side, sideNodes);

    auto sideTopo = elem.topology->boundary_type(side);
    auto sidePerm = sideTopo->permutation();

    unsigned firstNegativePermutation = sidePerm->num_positive_permutations();
    if (firstNegativePermutation >= sidePerm->num_permutations()) {
      firstNegativePermutation = 0;
    }

    std::vector<Ioss::Ordinal> permutationIndex =
        sidePerm->permutation_indices(static_cast<Ioss::Permutation>(firstNegativePermutation));
    mirrorNodes.clear();
    for (size_t i = 0; i < sideNodes.size(); i++) {
      mirrorNodes.push_back(sideNodes[permutationIndex[i]]);
    }
  }

  template int get_side_id(const IossElementData &elemData, const int num_dim, int nsnodes,
                           const int side_nodes[]);
  template int get_side_id(const IossElementData &elemData, const int num_dim, int nsnodes,
                           const int64_t side_nodes[]);

  template <typename INT>
  int get_side_id(const IossElementData &elemData,
                  const int              num_dim, // Mesh dimension
                  int                    nsnodes, /* The number of side nodes */
                  const INT              side_nodes[])         /* The list of side node IDs */
  {
    const Ioss::ElementTopology *topo    = elemData.topology;
    const EntityId              *connect = Data(elemData.localConnectivity);

    int count = 0;
    /*  min_match for hex elements means that min_match+1 nodes
        on a face of a hex must match to return the side of the
        hex on which the nodes exist, i.e., if 3/4 nodes of a hex
        match, then it might be considered connected.  If it is
        connected, then min_match states if 3/4 nodes is considered a face of
        a hex or not */

    /*  Default for min_match is 3, 2 is trial and error stuff */

    const int min_match = 3;
    /* const int min_match = 2; */

    /* check if this is a degenerate face */
    int                dup = 0;
    std::array<int, 9> location{};
    for (int i = 0; i < (nsnodes - 1); i++) {
      for (int j = (i + 1); j < nsnodes; j++) {
        if (side_nodes[i] == side_nodes[j]) {
          location[dup++] = i; /* location of duplicated node */
        }
      }
    }

    int nnodes = topo->number_nodes();

    /* Find all of the side nodes in the connect table */
    int num = 0;
    for (int i = 0; i < nnodes; i++) {
      for (int j = 0; j < nsnodes; j++) {
        if (connect[i] == side_nodes[j]) {
          num++;
          break;
        }
      }
      if (num == nsnodes) {
        break;
      }
    }

    /* I commented out the conditional statement causing the
       error if 2 hexes only share 3 out of 4 nodes.  I replaced
       this with what is seen below.  It works, but only for
       this particular case */

    /* the following ifdef is used to determine face adjacency
       old way:  numnodes on face must match on both elements
       new way:  only 3/4 of hex nodes have to match to be face adjacent */

    if (((num < nsnodes - 1) && (num >= 2)) || (num != nsnodes)) {
      EXPECT_TRUE(false) << "fatal: not all side nodes in connect table for element";
      return -1;
    }

    ElementType etype = get_elem_type(topo, num_dim);

    /* Find the side ID */
    switch (etype) {
    case ElementType::BAR1D2:
    case ElementType::BAR1D3:
      /* SIDE 1 */
      if (side_nodes[0] == connect[0]) {
        return 1;
      }
      if (side_nodes[0] == connect[1]) {
        return 2;
      }
      break;

    case ElementType::BAR2:
    case ElementType::BAR3:
    case ElementType::SHELL2:
    case ElementType::SHELL3:
      /* SIDE 1 */
      if (side_nodes[0] == connect[0] && side_nodes[1] == connect[1]) {
        return 1;
      }
      break;
    case ElementType::QUAD4:
    case ElementType::QUAD8:
    case ElementType::QUAD9:
      /* SIDE 1 */
      if (side_nodes[0] == connect[0] && side_nodes[1] == connect[1]) {
        return 1;
      }

      /* SIDE 2 */
      if (side_nodes[0] == connect[1] && side_nodes[1] == connect[2]) {
        return 2;
      }

      /* SIDE 3 */
      if (side_nodes[0] == connect[2] && side_nodes[1] == connect[3]) {
        return 3;
      }

      /* SIDE 4 */
      if (side_nodes[0] == connect[3] && side_nodes[1] == connect[0]) {
        return 4;
      }

      break;

    case ElementType::TRI3:
    case ElementType::TRI4:
    case ElementType::TRI6:
    case ElementType::TRI7:
      /* SIDE 1 */
      if (side_nodes[0] == connect[0] && side_nodes[1] == connect[1]) {
        return 1;
      }

      /* SIDE 2 */
      if (side_nodes[0] == connect[1] && side_nodes[1] == connect[2]) {
        return 2;
      }

      /* SIDE 3 */
      if (side_nodes[0] == connect[2] && side_nodes[1] == connect[0]) {
        return 3;
      }

      break;

    case ElementType::TET4:
    case ElementType::TET10:
    case ElementType::TET14:
    case ElementType::TET15:
    case ElementType::TET8:
      /* check the # of side nodes */
      if (nsnodes < 3) {
        return 0;
      }

      /* SIDE 1, 3, or 4 */
      if ((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 3] == connect[1] && side_nodes[(2 + num) % 3] == connect[3]) {
          return 1;
        }
        if (side_nodes[(1 + num) % 3] == connect[3] && side_nodes[(2 + num) % 3] == connect[2]) {
          return 3;
        }
        if (side_nodes[(1 + num) % 3] == connect[2] && side_nodes[(2 + num) % 3] == connect[1]) {
          return 4;
        }
      }

      /* SIDE 2 */
      if ((num = in_list(connect[1], nsnodes, side_nodes)) >= 0) {
        if (side_nodes[(1 + num) % 3] == connect[2] && side_nodes[(2 + num) % 3] == connect[3]) {
          return 2;
        }
      }

      break;

    case ElementType::HEX8:
    case ElementType::HEX16:
    case ElementType::HEX20:
    case ElementType::HEX27:
    case ElementType::HEXSHELL: /* this should be the same as a HEX element */
      /* check the # of side nodes */
      if (nsnodes < 4) {
        return 0;
      }

      /* SIDE 1 */
      if ((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
        count = 0;
        count += numbermatch(side_nodes, 1, num, 4, connect[1]);
        count += numbermatch(side_nodes, 2, num, 4, connect[5]);
        count += numbermatch(side_nodes, 3, num, 4, connect[4]);
        if (count >= min_match) {
          return 1;
        }

        /* if this is the duplicated node, then find the next occurrence */
        if (dup) {
          for (int i = 0; i < dup; i++) {
            if (connect[0] == side_nodes[location[i]]) {
              num = in_list(connect[0], (nsnodes - num), &(side_nodes[num + 1])) + location[i] + 1;
              count = 0;
              count += numbermatch(side_nodes, 1, num, 4, connect[1]);
              count += numbermatch(side_nodes, 2, num, 4, connect[5]);
              count += numbermatch(side_nodes, 3, num, 4, connect[4]);
              if (count >= min_match) {
                return 1;
              }
            }
          }
        }
      }

      /* SIDE 2 */
      if ((num = in_list(connect[1], nsnodes, side_nodes)) >= 0) {
        count = 0;
        count += numbermatch(side_nodes, 1, num, 4, connect[2]);
        count += numbermatch(side_nodes, 2, num, 4, connect[6]);
        count += numbermatch(side_nodes, 3, num, 4, connect[5]);
        if (count >= min_match) {
          return 2;
        }

        /* if this is the duplicated node, then find the next occurrence */
        if (dup) {
          for (int i = 0; i < dup; i++) {
            if (connect[1] == side_nodes[location[i]]) {
              num = in_list(connect[1], (nsnodes - num), &(side_nodes[num + 1])) + location[i] + 1;
              count = 0;
              count += numbermatch(side_nodes, 1, num, 4, connect[2]);
              count += numbermatch(side_nodes, 2, num, 4, connect[6]);
              count += numbermatch(side_nodes, 3, num, 4, connect[5]);
              if (count >= min_match) {
                return 2;
              }
            }
          }
        }
      }

      /* SIDE 3 */
      if ((num = in_list(connect[2], nsnodes, side_nodes)) >= 0) {
        count = 0;
        count += numbermatch(side_nodes, 1, num, 4, connect[3]);
        count += numbermatch(side_nodes, 2, num, 4, connect[7]);
        count += numbermatch(side_nodes, 3, num, 4, connect[6]);
        if (count >= min_match) {
          return 3;
        }

        /* if this is the duplicated node, then find the next occurrence */
        if (dup) {
          for (int i = 0; i < dup; i++) {
            if (connect[2] == side_nodes[location[i]]) {
              num = in_list(connect[2], (nsnodes - num), &(side_nodes[num + 1])) + location[i] + 1;
              count = 0;
              count += numbermatch(side_nodes, 1, num, 4, connect[3]);
              count += numbermatch(side_nodes, 2, num, 4, connect[7]);
              count += numbermatch(side_nodes, 3, num, 4, connect[6]);
              if (count >= min_match) {
                return 3;
              }
            }
          }
        }
      }

      /* SIDE 4 */
      if ((num = in_list(connect[3], nsnodes, side_nodes)) >= 0) {
        count = 0;
        count += numbermatch(side_nodes, 1, num, 4, connect[0]);
        count += numbermatch(side_nodes, 2, num, 4, connect[4]);
        count += numbermatch(side_nodes, 3, num, 4, connect[7]);
        if (count >= min_match) {
          return 4;
        }

        /* if this is the duplicated node, then find the next occurrence */
        if (dup) {
          for (int i = 0; i < dup; i++) {
            if (connect[3] == side_nodes[location[i]]) {
              num = in_list(connect[3], (nsnodes - num), &(side_nodes[num + 1])) + location[i] + 1;
              count = 0;
              count += numbermatch(side_nodes, 1, num, 4, connect[0]);
              count += numbermatch(side_nodes, 2, num, 4, connect[4]);
              count += numbermatch(side_nodes, 3, num, 4, connect[7]);
              if (count >= min_match) {
                return 4;
              }
            }
          }
        }
      }

      /* SIDE 5 */
      if ((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
        count = 0;
        count += numbermatch(side_nodes, 1, num, 4, connect[3]);
        count += numbermatch(side_nodes, 2, num, 4, connect[2]);
        count += numbermatch(side_nodes, 3, num, 4, connect[1]);
        if (count >= min_match) {
          return 5;
        }

        /* if this is the duplicated node, then find the next occurrence */
        if (dup) {
          for (int i = 0; i < dup; i++) {
            if (connect[0] == side_nodes[location[i]]) {
              num = in_list(connect[0], (nsnodes - num), &(side_nodes[num + 1])) + location[i] + 1;
              count = 0;
              count += numbermatch(side_nodes, 1, num, 4, connect[3]);
              count += numbermatch(side_nodes, 2, num, 4, connect[2]);
              count += numbermatch(side_nodes, 3, num, 4, connect[1]);
              if (count >= min_match) {
                return 5;
              }
            }
          }
        }
      }

      /* SIDE 6 */
      if ((num = in_list(connect[4], nsnodes, side_nodes)) >= 0) {
        count = 0;
        count += numbermatch(side_nodes, 1, num, 4, connect[5]);
        count += numbermatch(side_nodes, 2, num, 4, connect[6]);
        count += numbermatch(side_nodes, 3, num, 4, connect[7]);
        if (count >= min_match) {
          return 6;
        }

        /* if this is the duplicated node, then find the next occurrence */
        if (dup) {
          for (int i = 0; i < dup; i++) {
            if (connect[4] == side_nodes[location[i]]) {
              num = in_list(connect[4], (nsnodes - num), &(side_nodes[num + 1])) + location[i] + 1;
              count = 0;
              count += numbermatch(side_nodes, 1, num, 4, connect[5]);
              count += numbermatch(side_nodes, 2, num, 4, connect[6]);
              count += numbermatch(side_nodes, 3, num, 4, connect[7]);
              if (count >= min_match) {
                return 6;
              }
            }
          }
        }
      }

      break;

    case ElementType::TSHELL3:
    case ElementType::TSHELL4:
    case ElementType::TSHELL6:
    case ElementType::TSHELL7:

      /* 2D sides */
      if (nsnodes == 2 || (etype == ElementType::TSHELL6 && nsnodes == 3) ||
          (etype == ElementType::TSHELL7 && nsnodes == 3)) {
        /* SIDE 3 */
        if (side_nodes[0] == connect[0] && side_nodes[1] == connect[1]) {
          return 3;
        }

        /* SIDE 4 */
        if (side_nodes[0] == connect[1] && side_nodes[1] == connect[2]) {
          return 4;
        }

        /* SIDE 5 */
        if (side_nodes[0] == connect[2] && side_nodes[1] == connect[0]) {
          return 5;
        }
      }

      /* 3D faces */
      else if (nsnodes == 3 || nsnodes == 4 || nsnodes == 6 || nsnodes == 7) {
        if ((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
          if (side_nodes[(1 + num) % 3] == connect[1] && side_nodes[(2 + num) % 3] == connect[2]) {
            return 1;
          }
          if (side_nodes[(1 + num) % 3] == connect[2] && side_nodes[(2 + num) % 3] == connect[1]) {
            return 2;
          }
        }
      }

      break;

    case ElementType::SHELL4:
    case ElementType::SHELL8:
    case ElementType::SHELL9:

      /* 2D sides */
      if (nsnodes == 2 || nsnodes == 3) {
        /* SIDE 3 */
        if (side_nodes[0] == connect[0] && side_nodes[1] == connect[1]) {
          return 3;
        }

        /* SIDE 4 */
        if (side_nodes[0] == connect[1] && side_nodes[1] == connect[2]) {
          return 4;
        }

        /* SIDE 5 */
        if (side_nodes[0] == connect[2] && side_nodes[1] == connect[3]) {
          return 5;
        }

        /* SIDE 6 */
        if (side_nodes[0] == connect[3] && side_nodes[1] == connect[0]) {
          return 6;
        }
      }

      /* 3D faces */
      else if (nsnodes == 4 || nsnodes == 8 || nsnodes == 9) {

        /* SIDE 1 */
        if ((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
          if (side_nodes[(1 + num) % 4] == connect[1] && side_nodes[(2 + num) % 4] == connect[2] &&
              side_nodes[(3 + num) % 4] == connect[3]) {
            return 1;
          }
          if (side_nodes[(1 + num) % 4] == connect[3] && side_nodes[(2 + num) % 4] == connect[2] &&
              side_nodes[(3 + num) % 4] == connect[1]) {
            return 2;
          }
        }
      }

      break;

    case ElementType::WEDGE6:
    case ElementType::WEDGE12:
    case ElementType::WEDGE15:
    case ElementType::WEDGE16:
    case ElementType::WEDGE20:
    case ElementType::WEDGE21:

      /* quad sides */
      if (nsnodes == 4 || nsnodes == 8 || nsnodes == 9) {
        if ((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
          if (side_nodes[(1 + num) % 4] == connect[1] && side_nodes[(2 + num) % 4] == connect[4] &&
              side_nodes[(3 + num) % 4] == connect[3]) {
            return 1;
          }
          if (side_nodes[(1 + num) % 4] == connect[3] && side_nodes[(2 + num) % 4] == connect[5] &&
              side_nodes[(3 + num) % 4] == connect[2]) {
            return 3;
          }
        }

        if ((num = in_list(connect[1], nsnodes, side_nodes)) >= 0) {
          if (side_nodes[(1 + num) % 4] == connect[2] && side_nodes[(2 + num) % 4] == connect[5] &&
              side_nodes[(3 + num) % 4] == connect[4]) {
            return 2;
          }
        }
      }

      /* triangle sides */
      else if (nsnodes == 3 || nsnodes == 6 || nsnodes == 7) {
        if ((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
          if (side_nodes[(1 + num) % 3] == connect[2] && side_nodes[(2 + num) % 3] == connect[1]) {
            return 4;
          }
        }

        if ((num = in_list(connect[3], nsnodes, side_nodes)) >= 0) {
          if (side_nodes[(1 + num) % 3] == connect[4] && side_nodes[(2 + num) % 3] == connect[5]) {
            return 5;
          }
        }
      }

      break;

    case ElementType::PYRAMID5:
    case ElementType::PYRAMID13:
    case ElementType::PYRAMID14:
    case ElementType::PYRAMID18:
    case ElementType::PYRAMID19:
      /* triangular sides */
      if (nsnodes == 3 || nsnodes == 6 || nsnodes == 7) {
        /* SIDE 1 */
        if ((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
          if (side_nodes[(1 + num) % 3] == connect[1] && side_nodes[(2 + num) % 3] == connect[4]) {
            return 1;
          }
          if (side_nodes[(1 + num) % 3] == connect[4] && side_nodes[(2 + num) % 3] == connect[3]) {
            return 4;
          }
        }

        /* SIDE 2 */
        if ((num = in_list(connect[1], nsnodes, side_nodes)) >= 0) {
          if (side_nodes[(1 + num) % 3] == connect[2] && side_nodes[(2 + num) % 3] == connect[4]) {
            return 2;
          }
        }

        /* SIDE 3 */
        if ((num = in_list(connect[2], nsnodes, side_nodes)) >= 0) {
          if (side_nodes[(1 + num) % 3] == connect[3] && side_nodes[(2 + num) % 3] == connect[4]) {
            return 3;
          }
        }
      }

      else if (nsnodes == 4 || nsnodes == 8 || nsnodes == 9) {
        /* SIDE 5 */
        if ((num = in_list(connect[0], nsnodes, side_nodes)) >= 0) {
          if (side_nodes[(1 + num) % 4] == connect[3] && side_nodes[(2 + num) % 4] == connect[2] &&
              side_nodes[(3 + num) % 4] == connect[1]) {
            return 5;
          }
        }
      }

      break;

    case ElementType::SPHERE: break;

    default: {
      std::string err_buff;
      err_buff = fmt::format("fatal: unknown element shape {} in function {}",
                             static_cast<int>(topo->shape()), __func__);
      EXPECT_TRUE(false) << err_buff;
      exit(1);
    }
    } /* End "switch()" */

    return 0;
  }

  template int get_side_id_hex_tet(const Ioss::ElementTopology *topo, const int *connect,
                                   int nsnodes, const int side_nodes[]);
  template int get_side_id_hex_tet(const Ioss::ElementTopology *topo, const int64_t *connect,
                                   int nsnodes, const int64_t side_nodes[]);

  template <typename INT>
  int get_side_id_hex_tet(const Ioss::ElementTopology *topo,    /* The element type */
                          const INT                   *connect, /* The element connectivity */
                          int                          nsnodes, /* The number of side nodes */
                          const INT                    side_nodes[])               /* The list of side node IDs */
  {
    std::vector<int> loc_node_ids;

    int nnodes = topo->number_nodes();

    /* Find the local node numbers for nodes forming the side */
    for (int i1 = 0; i1 < nnodes; i1++) {
      for (int i2 = 0; i2 < nsnodes; i2++) {
        if (connect[i1] == side_nodes[i2]) {
          loc_node_ids.push_back(i1 + 1);
          break;
        }
      }
      if (loc_node_ids.size() == nsnodes) {
        break;
      }
    }

    int lcnt = loc_node_ids.size();

    switch (topo->shape()) {
    case Ioss::ElementShape::TET: {
      auto il1 = in_list(1, lcnt, Data(loc_node_ids)) >= 0;
      auto il2 = in_list(2, lcnt, Data(loc_node_ids)) >= 0;
      auto il3 = in_list(3, lcnt, Data(loc_node_ids)) >= 0;
      auto il4 = in_list(4, lcnt, Data(loc_node_ids)) >= 0;

      if (il1 && il2 && il4) {
        return 1;
      }

      if (il2 && il3 && il4) {
        return 2;
      }

      if (il1 && il3 && il4) {
        return 3;
      }

      if (il1 && il2 && il3) {
        return 4;
      }
    } break;

    case Ioss::ElementShape::HEX: {
      auto il1 = in_list(1, lcnt, Data(loc_node_ids)) >= 0 ? 1 : 0;
      auto il2 = in_list(2, lcnt, Data(loc_node_ids)) >= 0 ? 1 : 0;
      auto il3 = in_list(3, lcnt, Data(loc_node_ids)) >= 0 ? 1 : 0;
      auto il4 = in_list(4, lcnt, Data(loc_node_ids)) >= 0 ? 1 : 0;
      auto il5 = in_list(5, lcnt, Data(loc_node_ids)) >= 0 ? 1 : 0;
      auto il6 = in_list(6, lcnt, Data(loc_node_ids)) >= 0 ? 1 : 0;
      auto il7 = in_list(7, lcnt, Data(loc_node_ids)) >= 0 ? 1 : 0;
      auto il8 = in_list(8, lcnt, Data(loc_node_ids)) >= 0 ? 1 : 0;

      if (il1 + il2 + il5 + il6 > 2) {
        return 1;
      }

      if (il2 + il3 + il6 + il7 > 2) {
        return 2;
      }

      if (il3 + il4 + il7 + il8 > 2) {
        return 3;
      }

      if (il1 + il4 + il5 + il8 > 2) {
        return 4;
      }

      if (il1 + il2 + il3 + il4 > 2) {
        return 5;
      }

      if (il5 + il6 + il7 + il8 > 2) {
        return 6;
      }
    } break;

    default: {
      std::string err_buff;
      err_buff = fmt::format("fatal: unknown element type {} in function {}",
                             static_cast<int>(topo->shape()), __func__);
      EXPECT_TRUE(false) << err_buff;
      exit(1);
    }

    } /* End "switch()" */

    return 0;
  } /*-------------------------End get_side_id_hex()---------------------------*/
} // namespace utest_util
