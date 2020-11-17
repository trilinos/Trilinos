/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/*
 * Define element types.
 */

/* 1-d elements */

#define BAR2 0
#define BAR3 1
#define SHELL2 2
#define SHELL3 3

/* 2-d elements */

#define QUAD4 14
#define QUAD8 18
#define QUAD9 19
#define TRI3 23
#define TRI4 24
#define TRI6 26
#define TRI7 27

/* 3-d elements */

#define HEX8 108
#define HEX16 116
#define HEX20 120
#define HEX27 127
#define TET4 204
#define TET10 210
#define TET8 208
#define TET14 214
#define TET15 215
#define SHELL4 304
#define SHELL8 308
#define SHELL9 309
#define SPHERE 401
#define WEDGE6 506
#define WEDGE12 512
#define WEDGE15 515
#define WEDGE16 516
#define WEDGE20 520
#define WEDGE21 521
#define HEXSHELL 608
#define TSHELL3 703
#define TSHELL4 704
#define TSHELL6 706
#define TSHELL7 707
#define PYRAMID5 805
#define PYRAMID13 813
#define PYRAMID14 814
#define PYRAMID18 818
#define PYRAMID19 819

/* define element data "request for information" types */

#define NNODES 1
#define NDIM 3
#define NINTERP 5
#define NN_SIDE 6

/******************************* PROTOTYPES FOR el_elm_info.c ****************/

extern int elem_info(int info, int ielem_type, int supp);

extern int get_type(const char string[], int nodes, int num_dim);
