/*
 * Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/* cdrcom.h - external structure is used to hook up with fortran
 *            common block /cdrcom/
 */
#pragma once

extern struct cdr
{
  int KWRTFL;
  int KRDFL;
  int KOUTFL;
  int KINFL;
  int KWRDSZ;
  int KBYTEL;
  int KCPW;
  int KBAUD;
  int KCOMTP;
} cdrcom;
