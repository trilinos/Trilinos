/*
 * Copyright(C) 1999-2020, 2022, 2024 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

/*
 * Define element types.
 */

/* 1-d elements */
enum Elements {
  BAR2,
  BAR3,
  SHELL2,
  SHELL3,

  /* 2-d elements */
  QUAD4,
  QUAD8,
  QUAD9,
  TRI3,
  TRI4,
  TRI6,
  TRI7,

  /* 3-d elements */
  HEX8,
  HEX16,
  HEX20,
  HEX27,
  TET4,
  TET10,
  TET8,
  TET14,
  TET15,
  SHELL4,
  SHELL8,
  SHELL9,
  SPHERE,
  WEDGE6,
  WEDGE12,
  WEDGE15,
  WEDGE16,
  WEDGE20,
  WEDGE21,
  HEXSHELL,
  TSHELL3,
  TSHELL4,
  TSHELL6,
  TSHELL7,
  PYRAMID5,
  PYRAMID13,
  PYRAMID14,
  PYRAMID18,
  PYRAMID19
};
/* define element data "request for information" types */

enum ElementRequest { NNODES = 1, NDIM = 3, NINTERP = 5, NN_SIDE = 6 };

/******************************* PROTOTYPES FOR el_elm_info.c ****************/

extern int elem_info(int info, int ielem_type, int supp);

extern int get_type(const char string[], int nodes, int num_dim);
