// Copyright(C) 1999-2020, 2022, 2023, 2024, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include <gtest/gtest.h>
#ifdef SEACAS_HAVE_MPI
#include <mpi.h>
#endif
#include <string>
#include <unordered_map>

#include <memory>
#include <string>
#if defined(_WIN32) && !defined(__MINGW32__)
#include <string.h>
#define strcasecmp  _stricmp
#define strncasecmp _strnicmp
#else
#include <strings.h>
#endif
#include <ostream>
#include <vector>

#include <cstdarg> // for va_end, va_arg, va_list, etc
#include <cstddef> // for size_t
#include <cstdio>  // for stderr
#include <cstdlib> // for exit, malloc
#include <fmt/format.h>
#include <fmt/ostream.h>

#include "Ioss_ElementTopology.h"

#include "IossMeshTypes.h"

namespace utest_util {

  /* Define element types */
  enum class ElementType {
    SPHERE,
    BAR2,
    BAR3,
    QUAD4,
    QUAD8,
    QUAD9,
    SHELL4,
    SHELL8,
    SHELL9,
    TRI3,
    TRI4,
    TRI6,
    TRI7,
    TSHELL3,
    TSHELL4,
    TSHELL6,
    TSHELL7,
    HEX8,
    HEX16,
    HEX20,
    HEX27,
    HEXSHELL,
    TET4,
    TET10,
    TET8,
    TET14,
    TET15,
    WEDGE6,
    WEDGE12,
    WEDGE15,
    WEDGE16,
    WEDGE20,
    WEDGE21,
    PYRAMID5,
    PYRAMID13,
    PYRAMID14,
    PYRAMID18,
    PYRAMID19,
    SHELL2,
    SHELL3,
    BAR1D2,
    BAR1D3,
    NULL_EL
  };

  bool is_hex(const Ioss::ElementTopology *topo);

  bool is_hex(ElementType etype);

  bool is_tet(const Ioss::ElementTopology *topo);

  bool is_tet(ElementType etype);

  bool is_wedge(const Ioss::ElementTopology *topo);

  bool is_wedge(ElementType etype);

  bool is_pyramid(const Ioss::ElementTopology *topo);

  bool is_pyramid(ElementType etype);

  bool is_3d_element(const Ioss::ElementTopology *topo);

  bool is_3d_element(ElementType etype);

  ElementType get_elem_type(const char *elem_name, const int num_nodes, const int num_dim);

  ElementType get_elem_type(const Ioss::ElementTopology *topo, int num_dim);

  /*
   *----------------------------------------------------------------------------
   * This function returns the list of nodes in a side of an element given
   * the element type, and the side id. It also returns the number of nodes
   * in that side.
   *****************************************************************************/
  void get_local_side_nodes(const IossElementData &elem, int side,
                            std::vector<EntityId> &sideNodes);

  /*
   *----------------------------------------------------------------------------
   * This function returns the node list for the mirror of the list
   * given. This will be the node list of a face that is connected
   * to this element on this face i.e the first negative permutation
   *****************************************************************************/
  void get_local_side_mirror_nodes(const IossElementData &elem, int side,
                                   std::vector<EntityId> &mirrorNodes);

  template <typename INT>
  int get_side_id(const IossElementData &elemData,
                  const int              num_dim, // Mesh dimension
                  int                    nsnodes, /* The number of side nodes */
                  const INT              side_nodes[]);        /* The list of side node IDs */

  template <typename INT>
  int get_side_id_hex_tet(const Ioss::ElementTopology *topo,    /* The element type */
                          const INT                   *connect, /* The element connectivity */
                          int                          nsnodes, /* The number of side nodes */
                          const INT side_nodes[]);              /* The list of side node IDs */

} // namespace utest_util
