/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/* stdtyp.h - standard types, from Plum, "C Programming Guidelines", p 17
 * adapted to Sun3 Unix by Pat McGee, C6
 */
#ifndef STDTYP_H
#define STDTYP_H

typedef char          tiny;  /* 8+ bit signed number */
typedef unsigned char utiny; /* 8+ bit unsigned number */
typedef char          tbits; /* 8+ bit thing for bit manipulation*/
typedef short         bits;  /* 16+ bit thing for bit manipulation */
typedef long          lbits; /* 32+ bit thing for bit manipulation */

typedef char anything; /* type to hold a pointer to anything */

#endif
