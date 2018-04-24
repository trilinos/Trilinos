/*
 * Copyright (C) 2009-2017 National Technology & Engineering Solutions of
 * Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
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

extern int get_type(char string[], int nodes, int num_dim);
