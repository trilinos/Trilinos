/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "elb_elem.h"
#include "elb_err.h"  // for error_report, Gen_Error
#include "elb_util.h" // for in_list
#include <cstddef>    // for size_t
#include <cstdlib>    // for exit
#include <cstring>    // for strncasecmp
#include <fmt/ostream.h>
#include <vector> // for vector

/*****************************************************************************/
/*****************************************************************************/
namespace {
  template <typename INT>
  inline int numbermatch(INT *sidenodes, size_t i, size_t j, size_t k, size_t value)
  {
    if ((size_t)sidenodes[(i + j) % k] == value) {
      return 1;
    }
    return 0;
  }
} // namespace

/*****************************************************************************/
/* Function get_elem_type() begins:
 *----------------------------------------------------------------------------
 * This function returns the type of element based on the ExodusII element
 * string and number of nodes.
 *
 * Need the number of dimensions in order to distinguish between
 * TRI elements in a 2d mesh from TRI elements in a 3d mesh.
 *****************************************************************************/
const char *elem_name_from_enum(const E_Type elem_type)
{
  static const char *elem_names[NULL_EL] = {
      "SPHERE",    "BAR2",      "BAR3",      "QUAD4",   "QUAD8",   "QUAD9",    "SHELL4",
      "SHELL8",    "SHELL9",    "TRI3",      "TRI4",    "TRI6",    "TRI7",     "TSHELL3",
      "TSHELL4",   "TSHELL6",   "TSHELL7",   "HEX8",    "HEX16",   "HEX20",    "HEX27",
      "HEXSHELL",  "TET4",      "TET10",     "TET8",    "TET14",   "TET15",    "WEDGE6",
      "WEDGE12",   "WEDGE15",   "WEDGE16",   "WEDGE20", "WEDGE21", "PYRAMID5", "PYRAMID13",
      "PYRAMID14", "PYRAMID18", "PYRAMID19", "SHELL2",  "SHELL3"};
  return elem_names[elem_type];
}

E_Type get_elem_type(const char *elem_name, const int num_nodes, const int num_dim)
{

  E_Type answer = NULL_EL;
  switch (elem_name[0]) {
  case 'h':
  case 'H':
    if (strncasecmp(elem_name, "HEX", 3) == 0) {
      switch (num_nodes) {
      case 8: answer = HEX8; break;
      case 12: answer = HEXSHELL; break;
      case 16: answer = HEX16; break;
      case 20: answer = HEX20; break;
      case 27: answer = HEX27; break;
      default:
        Gen_Error(0, "fatal: unsupported HEX element");
        error_report();
        exit(1);
      }
    }
    break;

  case 'c':
  case 'C':
    if (strncasecmp(elem_name, "CIRCLE", 6) == 0) {
      answer = SPHERE;
    }
    break;

  case 's':
  case 'S':
    if (strncasecmp(elem_name, "SPHERE", 6) == 0) {
      answer = SPHERE;
    }
    else if (strncasecmp(elem_name, "SHELL", 5) == 0) {
      switch (num_nodes) {
      case 2:
        if (num_dim == 2) {
          answer = SHELL2;
        }
        else {
          Gen_Error(0, "fatal: unsupported SHELL element");
          error_report();
          exit(1);
        }
        break;
      case 3:
        if (num_dim == 2) {
          answer = SHELL3;
        }
        else {
          Gen_Error(0, "fatal: unsupported SHELL element");
          error_report();
          exit(1);
        }
        break;
      case 4: answer = SHELL4; break;
      case 8: answer = SHELL8; break;
      case 9: answer = SHELL9; break;
      default:
        Gen_Error(0, "fatal: unsupported SHELL element");
        error_report();
        exit(1);
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
      case 2: answer = BAR2; break;
      case 3: answer = BAR3; break;
      default:
        Gen_Error(0, "fatal: unsupported BAR/BEAM/TRUSS element");
        error_report();
        exit(1);
      }
    }
    else if (strncasecmp(elem_name, "TRI", 3) == 0) {
      switch (num_nodes) {
      case 3:
        if (num_dim == 2) {
          answer = TRI3;
        }
        else {
          answer = TSHELL3;
        }
        break;
      case 4:
        if (num_dim == 2) {
          answer = TRI4;
        }
        else {
          answer = TSHELL4;
        }
        break;
      case 6:
        if (num_dim == 2) {
          answer = TRI6;
        }
        else {
          answer = TSHELL6;
        }
        break;
      case 7:
        if (num_dim == 2) {
          answer = TRI7;
        }
        else {
          answer = TSHELL7;
        }
        break;
      default:
        Gen_Error(0, "fatal: unsupported TRI element");
        error_report();
        exit(1);
      }
    }
    else if (strncasecmp(elem_name, "TET", 3) == 0) {
      switch (num_nodes) {
      case 4: answer = TET4; break;
      case 8: answer = TET8; break;
      case 10: answer = TET10; break;
      case 14: answer = TET14; break;
      case 15: answer = TET15; break;
      default:
        Gen_Error(0, "fatal: unsupported TET element");
        error_report();
        exit(1);
      }
    }
    break;

  case 'q':
  case 'Q':
    if (strncasecmp(elem_name, "QUAD", 4) == 0) {
      switch (num_nodes) {
      case 4:
        if (num_dim == 2) {
          answer = QUAD4;
        }
        else {
          answer = SHELL4;
        }
        break;
      case 8:
        if (num_dim == 2) {
          answer = QUAD8;
        }
        else {
          answer = SHELL8;
        }
        break;
      case 9:
        if (num_dim == 2) {
          answer = QUAD9;
        }
        else {
          answer = SHELL9;
        }
        break;
      default:
        Gen_Error(0, "fatal: unsupported QUAD element");
        error_report();
        exit(1);
      }
    }
    break;

  case 'w':
  case 'W':
    if (strncasecmp(elem_name, "WEDGE", 5) == 0) {
      switch (num_nodes) {
      case 6: answer = WEDGE6; break;
      case 12: answer = WEDGE12; break;
      case 15: answer = WEDGE15; break;
      case 16: answer = WEDGE16; break;
      case 20: answer = WEDGE20; break;
      case 21: answer = WEDGE21; break;
      default:
        Gen_Error(0, "fatal: unsupported WEDGE element");
        error_report();
        exit(1);
      }
    }
    break;

  case 'p':
  case 'P':
    if (strncasecmp(elem_name, "PYR", 3) == 0) {
      switch (num_nodes) {
      case 5: answer = PYRAMID5; break;
      case 13: answer = PYRAMID13; break;
      case 14: answer = PYRAMID14; break;
      case 18: answer = PYRAMID18; break;
      case 19: answer = PYRAMID19; break;
      default:
        Gen_Error(0, "fatal: unsupported PYRAMID element");
        error_report();
        exit(1);
      }
    }
    break;

  default: break;
  }

  if (answer == NULL_EL) {
    std::string errstr;
    errstr = fmt::format("fatal: unknown element type '{}' read", elem_name);
    Gen_Error(0, errstr);
    error_report();
    exit(1);
  }

  return answer;

} /*---------------------------End get_elem_type()---------------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Convenience functions for code readability
 *****************************************************************************/
int is_hex(E_Type etype)
{
  return static_cast<int>(etype == HEX8 || etype == HEX27 || etype == HEX20 || etype == HEXSHELL);
}

int is_tet(E_Type etype)
{
  return static_cast<int>(etype == TET4 || etype == TET10 || etype == TET8 || etype == TET14 ||
                          etype == TET15);
}

int is_wedge(E_Type etype)
{
  return static_cast<int>(etype == WEDGE6 || etype == WEDGE15 || etype == WEDGE16 ||
                          etype == WEDGE20 || etype == WEDGE21);
}

int is_pyramid(E_Type etype)
{
  return static_cast<int>(etype == PYRAMID5 || etype == PYRAMID13 || etype == PYRAMID14 ||
                          etype == PYRAMID18 || etype == PYRAMID19);
}

int is_3d_element(E_Type etype)
{
  return static_cast<int>((is_hex(etype) != 0) || (is_tet(etype) != 0) || (is_wedge(etype) != 0) ||
                          (is_pyramid(etype) != 0));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function get_elem_info() begins:
 *----------------------------------------------------------------------------
 * This function returns various information about the input element type.
 *****************************************************************************/
int get_elem_info(const int req, const E_Type etype)
{

  int answer = 0;

  switch (etype) /* Switch over the element type */
  {
  case BAR2:
    switch (req) {
    case NNODES: answer = 2; break;
    case NSIDE_NODES: answer = 2; break;
    case NSIDES: answer = 1; break;
    case NDIM: /* number of physical dimensions */ answer = 1; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case SHELL2:
    switch (req) {
    case NNODES: answer = 2; break;
    case NSIDE_NODES: answer = 2; break;
    case NSIDES: answer = 1; break;
    case NDIM: /* number of physical dimensions */ answer = 1; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case SHELL3:
    switch (req) {
    case NNODES: answer = 3; break;
    case NSIDE_NODES: answer = 2; break;
    case NSIDES: answer = 1; break;
    case NDIM: /* number of physical dimensions */ answer = 1; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case BAR3:
    switch (req) {
    case NNODES: answer = 3; break;
    case NSIDE_NODES: answer = 2; break;
    case NSIDES: answer = 1; break;
    case NDIM: /* number of physical dimensions */ answer = 1; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case SPHERE:
    switch (req) {
    case NNODES: answer = 1; break;
    case NSIDE_NODES: answer = 0; break;
    case NSIDES: answer = 0; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    }
    break;

  case QUAD4:    /* First order quad */
    switch (req) /* select type of information required*/
    {
    case NNODES: /* number of nodes */ answer = 4; break;
    case NDIM: /* number of physical dimensions */ answer = 2; break;
    case NSIDE_NODES: answer = 2; break;
    case NSIDES: answer = 4; break;
    default:
      Gen_Error(0, "fatal:unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case QUAD8:    /* 2nd order serendipity quad */
    switch (req) /* select type of information required */
    {
    case NNODES: /* number of nodes */ answer = 8; break;
    case NDIM: /* number of physical dimensions */ answer = 2; break;
    case NSIDE_NODES: answer = 3; break;
    case NSIDES: answer = 4; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case QUAD9:    /* biquadratic quadrilateral */
    switch (req) /* select type of information required */
    {
    case NNODES: /* number of nodes */ answer = 9; break;
    case NDIM: /* number of physical dimensions */ answer = 2; break;
    case NSIDE_NODES: answer = 3; break;
    case NSIDES: answer = 4; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  /* NOTE: cannot determine NSIDE_NODES for SHELL element */
  case SHELL4:
    switch (req) {
    case NNODES: answer = 4; break;
    case NSIDES: answer = 6; break;
    case NDIM: /* number of physical dimensions */ answer = 2; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case SHELL8:
    switch (req) {
    case NNODES: answer = 8; break;
    case NSIDES: answer = 6; break;
    case NDIM: /* number of physical dimensions */ answer = 2; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case SHELL9:
    switch (req) {
    case NNODES: answer = 9; break;
    case NSIDES: answer = 6; break;
    case NDIM: /* number of physical dimensions */ answer = 2; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case TRI3:
    switch (req) /* select type of information required */
    {
    case NNODES: /* number of nodes */ answer = 3; break;
    case NDIM: /* number of physical dimensions */ answer = 2; break;
    case NSIDE_NODES: answer = 2; break;
    case NSIDES: answer = 3; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case TRI4:
    switch (req) /* select type of information required */
    {
    case NNODES: /* number of nodes */ answer = 4; break;
    case NDIM: /* number of physical dimensions */ answer = 2; break;
    case NSIDE_NODES: answer = 2; break;
    case NSIDES: answer = 3; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case TRI6:
    switch (req) /* select type of information required */
    {
    case NNODES: /* number of nodes */ answer = 6; break;
    case NDIM: /* number of physical dimensions */ answer = 2; break;
    case NSIDE_NODES: answer = 3; break;
    case NSIDES: answer = 3; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case TRI7:
    switch (req) /* select type of information required */
    {
    case NNODES: /* number of nodes */ answer = 7; break;
    case NDIM: /* number of physical dimensions */ answer = 2; break;
    case NSIDE_NODES: answer = 3; break;
    case NSIDES: answer = 3; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  /* NOTE: cannot determine NSIDE_NODES for TSHELL element */
  case TSHELL3:
    switch (req) /* select type of information required */
    {
    case NNODES: /* number of nodes */ answer = 3; break;
    case NDIM: /* number of physical dimensions */ answer = 2; break;
    case NSIDES: answer = 5; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case TSHELL4:
    switch (req) /* select type of information required */
    {
    case NNODES: /* number of nodes */ answer = 4; break;
    case NDIM: /* number of physical dimensions */ answer = 2; break;
    case NSIDES: answer = 5; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case TSHELL6:
    switch (req) /* select type of information required */
    {
    case NNODES: /* number of nodes */ answer = 6; break;
    case NDIM: /* number of physical dimensions */ answer = 2; break;
    case NSIDES: answer = 5; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case TSHELL7:
    switch (req) /* select type of information required */
    {
    case NNODES: /* number of nodes */ answer = 7; break;
    case NDIM: /* number of physical dimensions */ answer = 2; break;
    case NSIDES: answer = 5; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case HEX8:     /* trilinear hexahedron */
    switch (req) /* select type of information required */
    {
    case NNODES: /* number of nodes */ answer = 8; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    case NSIDE_NODES: answer = 4; break;
    case NSIDES: answer = 6; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case HEX16: /* localization element NSNODES is not consistent... */
    switch (req) {
    case NNODES: /* number of nodes */ answer = 16; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    case NSIDES: answer = 6; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case HEX20:    /* serendipity triquadratic hexahedron */
    switch (req) /* select type of information required */
    {
    case NNODES: /* number of nodes */ answer = 20; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    case NSIDE_NODES: answer = 8; break;
    case NSIDES: answer = 6; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case HEX27:    /* triquadratic hexahedron */
    switch (req) /* select type of information required*/
    {
    case NNODES: /* number of nodes */ answer = 27; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    case NSIDE_NODES: answer = 9; break;
    case NSIDES: answer = 6; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  /* NOTE: cannot determine NSIDE_NODES for HEXSHELL element */
  case HEXSHELL:
    switch (req) {
    case NNODES: answer = 12; break;
    case NSIDES: answer = 6; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case TET4:     /* trilinear tetrahedron */
    switch (req) /* select type of information required*/
    {
    case NNODES: /* number of nodes */ answer = 4; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    case NSIDE_NODES: answer = 3; break;
    case NSIDES: answer = 4; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case TET10:    /* triquadradic tetrahedron */
    switch (req) /* select type of information required */
    {
    case NNODES: /* number of nodes */ answer = 10; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    case NSIDE_NODES: answer = 6; break;
    case NSIDES: answer = 4; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case TET14:
    switch (req) /* select type of information required */
    {
    case NNODES: /* number of nodes */ answer = 14; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    case NSIDE_NODES: answer = 7; break;
    case NSIDES: answer = 4; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case TET15:
    switch (req) /* select type of information required */
    {
    case NNODES: /* number of nodes */ answer = 15; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    case NSIDE_NODES: answer = 7; break;
    case NSIDES: answer = 4; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case TET8:     /* 8-node (midface nodes) tetrahedron */
    switch (req) /* select type of information required */
    {
    case NNODES: /* number of nodes */ answer = 8; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    case NSIDE_NODES: answer = 4; break;
    case NSIDES: answer = 4; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  /* NOTE: cannot determine NSIDE_NODES for WEDGE elements */
  case WEDGE6:
    switch (req) {
    case NNODES: answer = 6; break;
    case NSIDES: answer = 5; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case WEDGE12:
    switch (req) {
    case NNODES: answer = 12; break;
    case NSIDES: answer = 5; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case WEDGE15:
    switch (req) {
    case NNODES: answer = 15; break;
    case NSIDES: answer = 5; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case WEDGE16:
    switch (req) {
    case NNODES: answer = 16; break;
    case NSIDES: answer = 5; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case WEDGE20:
    switch (req) {
    case NNODES: answer = 20; break;
    case NSIDES: answer = 5; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case WEDGE21:
    switch (req) {
    case NNODES: answer = 21; break;
    case NSIDES: answer = 5; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  /* NOTE: cannot determine NSIDE_NODES for PYRAMID element */
  case PYRAMID5:
    switch (req) {
    case NNODES: answer = 5; break;
    case NSIDES: answer = 5; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case PYRAMID13:
    switch (req) {
    case NNODES: answer = 13; break;
    case NSIDES: answer = 5; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case PYRAMID14:
    switch (req) {
    case NNODES: answer = 14; break;
    case NSIDES: answer = 5; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case PYRAMID18:
    switch (req) {
    case NNODES: answer = 18; break;
    case NSIDES: answer = 5; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  case PYRAMID19:
    switch (req) {
    case NNODES: answer = 19; break;
    case NSIDES: answer = 5; break;
    case NDIM: /* number of physical dimensions */ answer = 3; break;
    default:
      Gen_Error(0, "fatal: unknown quantity");
      error_report();
      exit(1);
    }
    break;

  default:
    Gen_Error(0, "fatal: unknown or unimplemented element type");
    error_report();
    exit(1);
  }

  return answer;

} /*---------------------------End get_elem_info()---------------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function get_side_id() begins:
 *----------------------------------------------------------------------------
 * This function returns the Side ID (as used in ExodusII) given a list of
 * nodes on that side.
 *
 * Changed so that it is now order dependent, but independent of starting
 * node for 3-D sides. On 2-D sides (lines), the starting node is important.
 *
 * Now supoports degenerate faces in HEX elements.
 *****************************************************************************/
template int get_side_id(const E_Type etype, const int *connect, const int nsnodes,
                         int side_nodes[], const int skip_check, const int partial_adj);
template int get_side_id(const E_Type etype, const int64_t *connect, const int nsnodes,
                         int64_t side_nodes[], const int skip_check, const int partial_adj);

template <typename INT>
int get_side_id(const E_Type etype, const INT *connect, const int nsnodes, INT side_nodes[],
                const int skip_check, const int partial_adj)
{
  int dup;
  int location[9];
  int count;
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
  dup = 0;
  for (int i = 0; i < (nsnodes - 1); i++) {
    for (int j = (i + 1); j < nsnodes; j++) {
      if (side_nodes[i] == side_nodes[j]) {
        location[dup++] = i; /* location of duplicated node */
      }
    }
  }

  int nnodes = get_elem_info(NNODES, etype);

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

  if (((partial_adj == 1) && (num < nsnodes - 1) && (num >= 2)) ||
      ((partial_adj != 1) && (num != nsnodes))) {
    if (skip_check) {
      if (skip_check == 1) /* print only if skip_check is 1 (not > 1) */
        Gen_Error(0, "warning: not all side nodes in connect table for element");
    }
    else {
      Gen_Error(0, "fatal: not all side nodes in connect table for element");
      return -1;
    }
  }

  if ((partial_adj == 1) && (num != nsnodes)) {
    return 0;
  }

  /* Find the side ID */
  switch (etype) {
  case BAR2:
  case BAR3:
  case SHELL2:
  case SHELL3:
    /* SIDE 1 */
    if (side_nodes[0] == connect[0] && side_nodes[1] == connect[1]) {
      return 1;
    }
    break;
  case QUAD4:
  case QUAD8:
  case QUAD9:
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

  case TRI3:
  case TRI4:
  case TRI6:
  case TRI7:
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

  case TET4:
  case TET10:
  case TET14:
  case TET15:
  case TET8:
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

  case HEX8:
  case HEX16:
  case HEX20:
  case HEX27:
  case HEXSHELL: /* this should be the same as a HEX element */
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
            num   = in_list(connect[0], (nsnodes - num), &(side_nodes[num + 1])) + location[i] + 1;
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
            num   = in_list(connect[1], (nsnodes - num), &(side_nodes[num + 1])) + location[i] + 1;
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
            num   = in_list(connect[2], (nsnodes - num), &(side_nodes[num + 1])) + location[i] + 1;
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
            num   = in_list(connect[3], (nsnodes - num), &(side_nodes[num + 1])) + location[i] + 1;
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
            num   = in_list(connect[0], (nsnodes - num), &(side_nodes[num + 1])) + location[i] + 1;
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
            num   = in_list(connect[4], (nsnodes - num), &(side_nodes[num + 1])) + location[i] + 1;
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

  case SHELL4:
  case SHELL8:
  case SHELL9:

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

  case WEDGE6:
  case WEDGE12:
  case WEDGE15:
  case WEDGE16:
  case WEDGE20:
  case WEDGE21:

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

  case TSHELL3:
  case TSHELL4:
  case TSHELL6:
  case TSHELL7:

    /* 2D sides */
    if (nsnodes == 2 || (etype == TSHELL6 && nsnodes == 3) || (etype == TSHELL7 && nsnodes == 3)) {
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

  case PYRAMID5:
  case PYRAMID13:
  case PYRAMID14:
  case PYRAMID18:
  case PYRAMID19:
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

  case SPHERE: break;

  default: {
    std::string err_buff;
    err_buff = fmt::format("fatal: unknown element type {} in function {}", static_cast<int>(etype),
                           __func__);
    Gen_Error(0, err_buff);
    error_report();
    exit(1);
  }
  } /* End "switch(etype)" */

  return 0;
} /*---------------------------End get_side_id()-----------------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function get_side_id_hex_tet() begins:
 *----------------------------------------------------------------------------
 * This function returns the Side ID (as used in ExodusII) of a HEX element
 * given a list of at least 3 nodes on that side. This function works on
 * the fact that any three nodes define a side of a hex. When a tet is
 * connected to a side of a hex, there are only three nodes connecting
 * the two. In this case a side id can be found.
 *****************************************************************************/
template int get_side_id_hex_tet(const E_Type etype, const int *connect, int nsnodes,
                                 const int side_nodes[]);
template int get_side_id_hex_tet(const E_Type etype, const int64_t *connect, int nsnodes,
                                 const int64_t side_nodes[]);

template <typename INT>
int get_side_id_hex_tet(const E_Type etype,     /* The element type */
                        const INT *  connect,   /* The element connectivity */
                        int          nsnodes,   /* The number of side nodes */
                        const INT    side_nodes[]) /* The list of side node IDs */
{
  int              nnodes;
  int              lcnt;
  int              i1;
  int              i2;
  std::vector<int> loc_node_ids(MAX_SIDE_NODES);

  nnodes = get_elem_info(NNODES, etype);

  /* Find the local node numbers for nodes forming the side */
  lcnt = 0;
  for (i1 = 0; i1 < nnodes; i1++) {
    for (i2 = 0; i2 < nsnodes; i2++) {
      if (connect[i1] == side_nodes[i2]) {
        loc_node_ids[lcnt++] = i1 + 1;
        break;
      }
    }
    if (lcnt == nsnodes) {
      break;
    }
  }

  switch (etype) {
  case TET4:
  case TET10:
  case TET8:
  case TET14:
  case TET15: {
    auto il1 = in_list(1, lcnt, loc_node_ids.data()) >= 0;
    auto il2 = in_list(2, lcnt, loc_node_ids.data()) >= 0;
    auto il3 = in_list(3, lcnt, loc_node_ids.data()) >= 0;
    auto il4 = in_list(4, lcnt, loc_node_ids.data()) >= 0;

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

  case HEX8:
  case HEX16:
  case HEX20:
  case HEX27: {
    auto il1 = in_list(1, lcnt, loc_node_ids.data()) >= 0 ? 1 : 0;
    auto il2 = in_list(2, lcnt, loc_node_ids.data()) >= 0 ? 1 : 0;
    auto il3 = in_list(3, lcnt, loc_node_ids.data()) >= 0 ? 1 : 0;
    auto il4 = in_list(4, lcnt, loc_node_ids.data()) >= 0 ? 1 : 0;
    auto il5 = in_list(5, lcnt, loc_node_ids.data()) >= 0 ? 1 : 0;
    auto il6 = in_list(6, lcnt, loc_node_ids.data()) >= 0 ? 1 : 0;
    auto il7 = in_list(7, lcnt, loc_node_ids.data()) >= 0 ? 1 : 0;
    auto il8 = in_list(8, lcnt, loc_node_ids.data()) >= 0 ? 1 : 0;

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
    err_buff = fmt::format("fatal: unknown element type {} in function {}", static_cast<int>(etype),
                           __func__);
    Gen_Error(0, err_buff);
    error_report();
    exit(1);
  }

  } /* End "switch(etype)" */

  return 0;
} /*-------------------------End get_side_id_hex()---------------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function ss_to_node_list() begins:
 *----------------------------------------------------------------------------
 * This function returns the list of nodes in a side of an element given
 * the element type, and the side id. It also returns the number of nodes
 * in that side.
 *****************************************************************************/
template int ss_to_node_list(const E_Type etype, const int *connect, int side_num,
                             int ss_node_list[]);
template int ss_to_node_list(const E_Type etype, const int64_t *connect, int side_num,
                             int64_t ss_node_list[]);

template <typename INT>
int ss_to_node_list(const E_Type etype,    /* The element type */
                    const INT *  connect,  /* The element connectivity */
                    int          side_num, /* The element side number */
                    INT          ss_node_list[])    /* The list of side node IDs */
{
  int i = 0;

  /*
   * This function returns a list of global node numbers forming a
   * side set.
   */

  /* triangle */
  static int tri_table[3][3] = {
      {1, 2, 4}, // side 1
      {2, 3, 5}, // side 2
      {3, 1, 6}  // side 3
  };

  /* tshell */
  static int tshell_table[2][7] = {
      {1, 2, 3, 4, 5, 6, 7}, // side 1
      {1, 3, 2, 6, 5, 4, 7}  // side 2
  };

  /* quad */
  static int quad_table[4][3] = {
      {1, 2, 5}, // side 1
      {2, 3, 6}, // side 2
      {3, 4, 7}, // side 3
      {4, 1, 8}  // side 4
  };

  /* shell */
  static int shell_table[2][9] = {
      {1, 2, 3, 4, 5, 6, 7, 8, 9}, // side 1
      {1, 4, 3, 2, 8, 7, 6, 5, 9}  // side 2
  };

  /* tetra */
  static int tetra_table[4][7] = {
      {1, 2, 4, 5, 9, 8, 14},  /* Side 1 nodes */
      {2, 3, 4, 6, 10, 9, 12}, /* Side 2 nodes */
      {1, 4, 3, 8, 10, 7, 13}, /* Side 3 nodes */
      {1, 3, 2, 7, 6, 5, 11}   /* Side 4 nodes */
  };

  /* wedge */
  /* wedge 6 or 7 */
  static int wedge6_table[5][4] = {
      {1, 2, 5, 4}, /* Side 1 nodes -- quad     */
      {2, 3, 6, 5}, /* Side 2 nodes -- quad     */
      {1, 4, 6, 3}, /* Side 3 nodes -- quad     */
      {1, 3, 2, 0}, /* Side 4 nodes -- triangle */
      {4, 5, 6, 0}  /* Side 5 nodes -- triangle */
  };

  /* wedge 12 */
  static int wedge12_table[5][6] = {
      {1, 2, 5, 4, 7, 10},  /* Side 1 nodes -- quad4     */
      {2, 3, 6, 5, 8, 11},  /* Side 2 nodes -- quad4     */
      {3, 1, 4, 6, 9, 12},  /* Side 3 nodes -- quad4     */
      {1, 3, 2, 9, 8, 7},   /* Side 4 nodes -- triangle6 */
      {4, 5, 6, 10, 11, 12} /* Side 5 nodes -- triangle6 */
  };

  /* wedge 15 or 16 */
  static int wedge15_table[5][8] = {
      {1, 2, 5, 4, 7, 11, 13, 10}, /* Side 1 nodes -- quad     */
      {2, 3, 6, 5, 8, 12, 14, 11}, /* Side 2 nodes -- quad     */
      {1, 4, 6, 3, 10, 15, 12, 9}, /* Side 3 nodes -- quad     */
      {1, 3, 2, 9, 8, 7, 0, 0},    /* Side 4 nodes -- triangle */
      {4, 5, 6, 13, 14, 15, 0, 0}  /* Side 5 nodes -- triangle */
  };

  /* wedge 20 */
  static int wedge20_table[5][9] = {
      {1, 2, 5, 4, 7, 11, 13, 10, 20}, /* Side 1 nodes -- quad     */
      {2, 3, 6, 5, 8, 12, 14, 11, 18}, /* Side 2 nodes -- quad     */
      {1, 4, 6, 3, 10, 15, 12, 9, 19}, /* Side 3 nodes -- quad     */
      {1, 3, 2, 9, 8, 7, 16, 0, 0},    /* Side 4 nodes -- triangle */
      {4, 5, 6, 13, 14, 15, 17, 0, 0}  /* Side 5 nodes -- triangle */
  };

  /* wedge 21 */
  static int wedge21_table[5][9] = {
      {1, 2, 5, 4, 7, 11, 13, 10, 21}, /* Side 1 nodes -- quad     */
      {2, 3, 6, 5, 8, 12, 14, 11, 19}, /* Side 2 nodes -- quad     */
      {1, 4, 6, 3, 10, 15, 12, 9, 20}, /* Side 3 nodes -- quad     */
      {1, 3, 2, 9, 8, 7, 17, 0, 0},    /* Side 4 nodes -- triangle */
      {4, 5, 6, 13, 14, 15, 18, 0, 0}  /* Side 5 nodes -- triangle */
  };

  /* hex */
  static int hex_table[6][9] = {
      {1, 2, 6, 5, 9, 14, 17, 13, 26},  /* side 1 */
      {2, 3, 7, 6, 10, 15, 18, 14, 25}, /* side 2 */
      {3, 4, 8, 7, 11, 16, 19, 15, 27}, /* side 3 */
      {1, 5, 8, 4, 13, 20, 16, 12, 24}, /* side 4 */
      {1, 4, 3, 2, 12, 11, 10, 9, 22},  /* side 5 */
      {5, 6, 7, 8, 17, 18, 19, 20, 23}  /* side 6 */
  };

  /* hex16 */
  static int hex16_table[6][9] = {
      {1, 2, 6, 5, 9, 13, 0, 0},   /* side 1 */
      {2, 3, 7, 6, 10, 14, 0, 0},  /* side 2 */
      {3, 4, 8, 7, 11, 15, 0, 0},  /* side 3 */
      {4, 1, 5, 8, 4, 12, 16, 0},  /* side 4 */
      {1, 4, 3, 2, 12, 11, 10, 9}, /* side 5 */
      {5, 6, 7, 8, 13, 14, 15, 16} /* side 6 */
  };

  /* hexshell */
  static int hexshell_table[6][6] = {
      {1, 2, 6, 5, 10, 9},  // side 1
      {2, 3, 7, 6, 11, 10}, // side 2
      {3, 4, 8, 7, 12, 11}, // side 3
      {4, 1, 5, 8, 9, 12},  // side 4
      {1, 4, 3, 2, 0, 0},   // side 5
      {5, 6, 7, 8, 0, 0}    // side 6
  };

  /* pyramid */
  static int pyramid_table[5][9] = {
      {1, 2, 5, 6, 11, 10, 15, 0, 0}, // side 1 (tri)
      {2, 3, 5, 7, 12, 11, 16, 0, 0}, // side 2 (tri)
      {3, 4, 5, 8, 13, 12, 17, 0, 0}, // side 3 (tri)
      {4, 1, 5, 9, 10, 13, 18, 0, 0}, // side 4 (tri)
      {1, 4, 3, 2, 9, 8, 7, 6, 14}    // side 5 (quad)
  };

  static int bar_table[1][3] = {{1, 2, 3}};

  /* Locally decrement side_num */
  side_num--;

  /* Switch over the element type. */
  switch (etype) {
  case BAR2:
  case SHELL2:
    /* Bar1 has 1 side */
    for (i = 0; i < 2; i++) {
      ss_node_list[i] = connect[(bar_table[side_num][i] - 1)];
    }
    break;

  case BAR3:
  case SHELL3:
    /* Bar has 1 side */
    for (i = 0; i < 3; i++) {
      ss_node_list[i] = connect[(bar_table[side_num][i] - 1)];
    }
    break;

  case QUAD4:
    for (i = 0; i < 2; i++) {
      ss_node_list[i] = connect[(quad_table[side_num][i] - 1)];
    }
    break;

  case QUAD8:
  case QUAD9:
    for (i = 0; i < 3; i++) {
      ss_node_list[i] = connect[(quad_table[side_num][i] - 1)];
    }
    break;

  case SHELL4:
    switch (side_num) {
    case 0:
    case 1:
      for (i = 0; i < 4; i++) {
        ss_node_list[i] = connect[(shell_table[side_num][i] - 1)];
      }
      break;

    default:
      /*
       * sides 3, 4, 5, & 6 correspond to sides 1, 2, 3, & 4
       * of the quad element.
       */
      for (i = 0; i < 2; i++) {
        ss_node_list[i] = connect[(quad_table[(side_num - 2)][i] - 1)];
      }
      break;
    }
    break;

  case SHELL8:
    switch (side_num) {
    case 0:
    case 1:
      for (i = 0; i < 8; i++) {
        ss_node_list[i] = connect[(shell_table[side_num][i] - 1)];
      }
      break;

    default:
      /*
       * sides 3, 4, 5, & 6 correspond to sides 1, 2, 3, & 4
       * of the quad element.
       */
      for (i = 0; i < 3; i++) {
        ss_node_list[i] = connect[(quad_table[(side_num - 2)][i] - 1)];
      }
      break;
    }
    break;

  case SHELL9:
    switch (side_num) {
    case 0:
    case 1:
      for (i = 0; i < 9; i++) {
        ss_node_list[i] = connect[(shell_table[side_num][i] - 1)];
      }
      break;

    default:
      /*
       * sides 3, 4, 5, & 6 correspond to sides 1, 2, 3, & 4
       * of the quad element.
       */
      for (i = 0; i < 3; i++) {
        ss_node_list[i] = connect[(quad_table[(side_num - 2)][i] - 1)];
      }
      break;
    }
    break;

  case TRI3:
  case TRI4:
    for (i = 0; i < 2; i++) {
      ss_node_list[i] = connect[(tri_table[side_num][i] - 1)];
    }
    break;

  case TRI6:
  case TRI7:
    for (i = 0; i < 3; i++) {
      ss_node_list[i] = connect[(tri_table[side_num][i] - 1)];
    }
    break;

  case TSHELL3:
  case TSHELL4:
    switch (side_num) {
    case 0:
    case 1:
      for (i = 0; i < 3; i++) {
        ss_node_list[i] = connect[(tshell_table[side_num][i] - 1)];
      }
      break;

    default:
      /*
       * sides 3, 4 & 5 correspond to sides 1, 2 & 3
       * of the tri element.
       */
      for (i = 0; i < 2; i++) {
        ss_node_list[i] = connect[(tri_table[(side_num - 2)][i] - 1)];
      }
      break;
    }
    break;

  case TSHELL6:
  case TSHELL7:
    switch (side_num) {
    case 0:
    case 1:
      for (i = 0; i < 6; i++) {
        ss_node_list[i] = connect[(tshell_table[side_num][i] - 1)];
      }
      break;

    default:
      /*
       * sides 3, 4 & 5 correspond to sides 1, 2 & 3
       * of the tri element.
       */
      for (i = 0; i < 3; i++) {
        ss_node_list[i] = connect[(tri_table[(side_num - 2)][i] - 1)];
      }
      break;
    }
    break;

  case HEX8:
    for (i = 0; i < 4; i++) {
      ss_node_list[i] = connect[(hex_table[side_num][i] - 1)];
    }
    break;

  case HEX16:
    switch (side_num) {
    case 4:
    case 5:
      for (i = 0; i < 8; i++) {
        ss_node_list[i] = connect[(hex16_table[side_num][i] - 1)];
      }
      break;

    default:
      for (i = 0; i < 6; i++) {
        ss_node_list[i] = connect[(hex16_table[side_num][i] - 1)];
      }
      break;
    }
    break;

  case HEX20:
    for (i = 0; i < 8; i++) {
      ss_node_list[i] = connect[(hex_table[side_num][i] - 1)];
    }
    break;

  case HEX27:
    for (i = 0; i < 9; i++) {
      ss_node_list[i] = connect[(hex_table[side_num][i] - 1)];
    }
    break;

  case TET4:
    for (i = 0; i < 3; i++) {
      ss_node_list[i] = connect[(tetra_table[side_num][i] - 1)];
    }
    break;

  case TET10:
    for (i = 0; i < 6; i++) {
      ss_node_list[i] = connect[(tetra_table[side_num][i] - 1)];
    }
    break;

  case TET14:
  case TET15:
    for (i = 0; i < 7; i++) {
      ss_node_list[i] = connect[(tetra_table[side_num][i] - 1)];
    }
    break;

  case TET8:
    for (i = 0; i < 4; i++) {
      ss_node_list[i] = connect[(tetra_table[side_num][i] - 1)];
    }
    break;

  case WEDGE6:
    switch (side_num) {
    case 3:
    case 4:
      for (i = 0; i < 3; i++) {
        ss_node_list[i] = connect[(wedge6_table[side_num][i] - 1)];
      }
      break;

    default:
      for (i = 0; i < 4; i++) {
        ss_node_list[i] = connect[(wedge6_table[side_num][i] - 1)];
      }
      break;
    }
    break;

  case WEDGE12:
    for (i = 0; i < 6; i++) {
      ss_node_list[i] = connect[(wedge12_table[side_num][i] - 1)];
    }
    break;

  case WEDGE15:
  case WEDGE16:
    switch (side_num) {
    case 3:
    case 4:
      for (i = 0; i < 6; i++) {
        ss_node_list[i] = connect[(wedge15_table[side_num][i] - 1)];
      }
      break;

    default:
      for (i = 0; i < 8; i++) {
        ss_node_list[i] = connect[(wedge15_table[side_num][i] - 1)];
      }
      break;
    }
    break;

  case WEDGE20:
    switch (side_num) {
    case 3:
    case 4:
      for (i = 0; i < 7; i++) {
        ss_node_list[i] = connect[(wedge20_table[side_num][i] - 1)];
      }
      break;

    default:
      for (i = 0; i < 9; i++) {
        ss_node_list[i] = connect[(wedge20_table[side_num][i] - 1)];
      }
      break;
    }
    break;

  case WEDGE21:
    switch (side_num) {
    case 3:
    case 4:
      for (i = 0; i < 7; i++) {
        ss_node_list[i] = connect[(wedge21_table[side_num][i] - 1)];
      }
      break;

    default:
      for (i = 0; i < 9; i++) {
        ss_node_list[i] = connect[(wedge21_table[side_num][i] - 1)];
      }
      break;
    }
    break;

  case HEXSHELL:
    switch (side_num) {
    case 4:
    case 5:
      for (i = 0; i < 4; i++) {
        ss_node_list[i] = connect[(hexshell_table[side_num][i] - 1)];
      }
      break;

    default:
      for (i = 0; i < 6; i++) {
        ss_node_list[i] = connect[(hexshell_table[side_num][i] - 1)];
      }
      break;
    }
    break;

  case PYRAMID5:
    switch (side_num) {
    case 4:
      for (i = 0; i < 4; i++) {
        ss_node_list[i] = connect[(pyramid_table[side_num][i] - 1)];
      }
      break;

    default:
      for (i = 0; i < 3; i++) {
        ss_node_list[i] = connect[(pyramid_table[side_num][i] - 1)];
      }
      break;
    }
    break;

  case PYRAMID13:
    switch (side_num) {
    case 4:
      for (i = 0; i < 8; i++) {
        ss_node_list[i] = connect[(pyramid_table[side_num][i] - 1)];
      }
      break;

    default:
      for (i = 0; i < 6; i++) {
        ss_node_list[i] = connect[(pyramid_table[side_num][i] - 1)];
      }
      break;
    }
    break;

  case PYRAMID14:
    switch (side_num) {
    case 4:
      for (i = 0; i < 9; i++) {
        ss_node_list[i] = connect[(pyramid_table[side_num][i] - 1)];
      }
      break;

    default:
      for (i = 0; i < 6; i++) {
        ss_node_list[i] = connect[(pyramid_table[side_num][i] - 1)];
      }
      break;
    }
    break;

  case PYRAMID18:
  case PYRAMID19: /* Pyramid18 with mid-volume node */
    switch (side_num) {
    case 4:
      for (i = 0; i < 9; i++) {
        ss_node_list[i] = connect[(pyramid_table[side_num][i] - 1)];
      }
      break;

    default:
      for (i = 0; i < 7; i++) {
        ss_node_list[i] = connect[(pyramid_table[side_num][i] - 1)];
      }
      break;
    }
    break;

  case SPHERE: /* SHPERE's have no side sets */
  case NULL_EL: i = 0; break;

  } /* End "switch (etype)" */

  /* the variable "i" should be the number of positions that I filled */
  return (i);

} /*-------------------------End ss_to_node_list()---------------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function get_ss_mirror() begins:
 *----------------------------------------------------------------------------
 * This function returns the node list for the mirror of the list
 * given. This will be the node list of a face that is connected
 * to this element on this face.
 *****************************************************************************/
template int get_ss_mirror(const E_Type etype, const int *ss_node_list, int side_num,
                           int mirror_node_list[]);
template int get_ss_mirror(const E_Type etype, const int64_t *ss_node_list, int side_num,
                           int64_t mirror_node_list[]);

template <typename INT>
int get_ss_mirror(const E_Type etype,             /* The element type */
                  const INT *  ss_node_list,      /* The list of side node IDs */
                  int          side_num,          /* The element side number */
                  INT          mirror_node_list[] /* The list of the mirror side node IDs */
)
{
  int i = 0;

  /*
   * the following arrays are the conversion from the side to
   * an opposing face
   */

  /* line (1-d) */
  static int line_table[3] = {1, 0, 2};

  /* square (2-d) */
  static int sqr_table[9] = {0, 3, 2, 1, 7, 6, 5, 4, 8};

  /* square hexshell (2-d) */
  static int hs_table[6] = {0, 3, 2, 1, 5, 4};

  /* triangle (2-d) */
  static int tri_table[7] = {1, 0, 2, 3, 5, 4, 6};

  /***************************** execution begins ******************************/

  /* Switch over the element type. */
  switch (etype) {
  case BAR2:
  case SHELL2:
    for (i = 0; i < 2; i++) {
      mirror_node_list[i] = ss_node_list[line_table[i]];
    }
    break;
  case BAR3:
  case SHELL3:
    for (i = 0; i < 3; i++) {
      mirror_node_list[i] = ss_node_list[line_table[i]];
    }
    break;
  case QUAD4:
    for (i = 0; i < 2; i++) {
      mirror_node_list[i] = ss_node_list[line_table[i]];
    }
    break;

  case QUAD8:
  case QUAD9:
    for (i = 0; i < 3; i++) {
      mirror_node_list[i] = ss_node_list[line_table[i]];
    }
    break;

  case SHELL4:
    switch (side_num) {
    case 1:
    case 2:
      for (i = 0; i < 4; i++) {
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      }
      break;

    default:
      for (i = 0; i < 2; i++) {
        mirror_node_list[i] = ss_node_list[line_table[i]];
      }
      break;
    }
    break;

  case SHELL8:
    switch (side_num) {
    case 1:
    case 2:
      for (i = 0; i < 8; i++) {
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      }
      break;

    default:
      for (i = 0; i < 3; i++) {
        mirror_node_list[i] = ss_node_list[line_table[i]];
      }
      break;
    }
    break;

  case SHELL9:
    switch (side_num) {
    case 1:
    case 2:
      for (i = 0; i < 9; i++) {
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      }
      break;

    default:
      for (i = 0; i < 3; i++) {
        mirror_node_list[i] = ss_node_list[line_table[i]];
      }
      break;
    }
    break;

  case TRI3:
  case TRI4:
    for (i = 0; i < 2; i++) {
      mirror_node_list[i] = ss_node_list[line_table[i]];
    }
    break;

  case TRI6:
  case TRI7:
    for (i = 0; i < 3; i++) {
      mirror_node_list[i] = ss_node_list[line_table[i]];
    }
    break;

  case TSHELL3:
    switch (side_num) {
    case 1:
    case 2:
      for (i = 0; i < 3; i++) {
        mirror_node_list[i] = ss_node_list[tri_table[i]];
      }
      break;

    default:
      for (i = 0; i < 2; i++) {
        mirror_node_list[i] = ss_node_list[line_table[i]];
      }
      break;
    }
    break;

  case TSHELL4:
    switch (side_num) {
    case 1:
    case 2:
      for (i = 0; i < 4; i++) {
        mirror_node_list[i] = ss_node_list[tri_table[i]];
      }
      break;

    default:
      for (i = 0; i < 2; i++) {
        mirror_node_list[i] = ss_node_list[line_table[i]];
      }
      break;
    }
    break;

  case TSHELL6:
    switch (side_num) {
    case 1:
    case 2:
      for (i = 0; i < 6; i++) {
        mirror_node_list[i] = ss_node_list[tri_table[i]];
      }
      break;

    default:
      for (i = 0; i < 3; i++) {
        mirror_node_list[i] = ss_node_list[line_table[i]];
      }
      break;
    }
    break;

  case TSHELL7:
    switch (side_num) {
    case 1:
    case 2:
      for (i = 0; i < 7; i++) {
        mirror_node_list[i] = ss_node_list[tri_table[i]];
      }
      break;

    default:
      for (i = 0; i < 3; i++) {
        mirror_node_list[i] = ss_node_list[line_table[i]];
      }
      break;
    }
    break;

  case HEX8:
    for (i = 0; i < 4; i++) {
      mirror_node_list[i] = ss_node_list[sqr_table[i]];
    }
    break;

  case HEX16:
    switch (side_num) {
    case 4:
    case 5:
      for (i = 0; i < 8; i++) {
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      }
      break;

    default:
      for (i = 0; i < 6; i++) {
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      }
    }
    break;

  case HEX27:
    for (i = 0; i < 9; i++) {
      mirror_node_list[i] = ss_node_list[sqr_table[i]];
    }
    break;

  case HEX20:
    for (i = 0; i < 8; i++) {
      mirror_node_list[i] = ss_node_list[sqr_table[i]];
    }
    break;

  case TET4:
    for (i = 0; i < 3; i++) {
      mirror_node_list[i] = ss_node_list[tri_table[i]];
    }
    break;

  case TET8:
    for (i = 0; i < 4; i++) {
      mirror_node_list[i] = ss_node_list[tri_table[i]];
    }
    break;

  case TET10:
    for (i = 0; i < 6; i++) {
      mirror_node_list[i] = ss_node_list[tri_table[i]];
    }
    break;

  case TET14:
  case TET15:
    for (i = 0; i < 7; i++) {
      mirror_node_list[i] = ss_node_list[tri_table[i]];
    }
    break;

  case WEDGE6:
    switch (side_num) {
    case 4:
    case 5:
      for (i = 0; i < 3; i++) {
        mirror_node_list[i] = ss_node_list[tri_table[i]];
      }
      break;

    default:
      for (i = 0; i < 4; i++) {
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      }
      break;
    }
    break;

  case WEDGE12:
    for (i = 0; i < 6; i++) {
      mirror_node_list[i] = ss_node_list[tri_table[i]];
    }
    break;

  case WEDGE15:
  case WEDGE16:
    switch (side_num) {
    case 4:
    case 5:
      for (i = 0; i < 6; i++) {
        mirror_node_list[i] = ss_node_list[tri_table[i]];
      }
      break;

    default:
      for (i = 0; i < 8; i++) {
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      }
      break;
    }
    break;

  case WEDGE20:
  case WEDGE21:
    switch (side_num) {
    case 4:
    case 5:
      for (i = 0; i < 7; i++) {
        mirror_node_list[i] = ss_node_list[tri_table[i]];
      }
      break;

    default:
      for (i = 0; i < 9; i++) {
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      }
      break;
    }
    break;

  case HEXSHELL:
    switch (side_num) {
    case 5:
    case 6:
      for (i = 0; i < 4; i++) {
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      }
      break;

    default:
      for (i = 0; i < 6; i++) {
        mirror_node_list[i] = ss_node_list[hs_table[i]];
      }
      break;
    }
    break;

  case PYRAMID5:
    switch (side_num) {
    case 5:
      for (i = 0; i < 4; i++) {
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      }
      break;

    default:
      for (i = 0; i < 3; i++) {
        mirror_node_list[i] = ss_node_list[tri_table[i]];
      }
      break;
    }
    break;

  case PYRAMID13:
    switch (side_num) {
    case 5:
      for (i = 0; i < 8; i++) {
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      }
      break;

    default:
      for (i = 0; i < 6; i++) {
        mirror_node_list[i] = ss_node_list[tri_table[i]];
      }
      break;
    }
    break;

  case PYRAMID14:
    switch (side_num) {
    case 5:
      for (i = 0; i < 9; i++) {
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      }
      break;

    default:
      for (i = 0; i < 6; i++) {
        mirror_node_list[i] = ss_node_list[tri_table[i]];
      }
      break;
    }
    break;

  case PYRAMID18:
  case PYRAMID19:
    switch (side_num) {
    case 5:
      for (i = 0; i < 9; i++) {
        mirror_node_list[i] = ss_node_list[sqr_table[i]];
      }
      break;

    default:
      for (i = 0; i < 7; i++) {
        mirror_node_list[i] = ss_node_list[tri_table[i]];
      }
      break;
    }
    break;

  case SPHERE: /* SHPERE's have no side sets */
  case NULL_EL: i = 0; break;

  } /* End "switch (etype)" */

  /* the variable "i" should be the number of positions that I filled */
  return (i);

} /*-------------------------Ed get_ss_mirror()---------------------------*/
