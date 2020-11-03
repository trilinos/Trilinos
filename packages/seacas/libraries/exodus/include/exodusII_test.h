/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef __exodusII_test_h
#define __exodusII_test_h

#include "exodusII.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

int cCreateEdgeFace(int, char *[]);
int cReadEdgeFace(int, char *[]);

#ifdef __cplusplus
}
#endif /* __cplusplus */

inline int CreateEdgeFace(int argc, char *argv[]) { return cCreateEdgeFace(argc, argv); }
inline int ReadEdgeFace(int argc, char *argv[]) { return cReadEdgeFace(argc, argv); }

#endif /* __exodusII_test_h */
