/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/* Switchable timer routine  */
double lanc_seconds(void)
{
  extern int LANCZOS_TIME; /* perform detailed timing on Lanczos_SO? */
  double     seconds();

  if (LANCZOS_TIME) {
    return (seconds());
  }

  return (0);
}
