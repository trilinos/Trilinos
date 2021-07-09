/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "defs.h" // for FALSE, TRUE
#include <math.h>
#include <stdio.h> // for printf

/* Determine whether to pause in Lanczos */
int lanpause(int      j,         /* current step */
             int      lastpause, /* when last paused */
             int      interval,  /* interval between pauses */
             double **q,         /* the Lanczos vectors */
             int      n,         /* length of Lanczos vectors */
             int *    pausemode, /* which pausing criterion to use */
             int      version,   /* which version of sel. orth. we are using */
             double   beta       /* current off-diagonal value */
)
{
  extern int    DEBUG_EVECS;    /* debugging level for eigen computation */
  extern double DOUBLE_EPSILON; /* machine precision */
  double        paige_dot;      /* q[j]^T q[1] */
  double        paigetol;       /* pause if paigedot > paigetol */
  double        dot();          /* standard dot product */
  void          checkorth();

  /* Check orthogonality of last Lanczos vector against previous ones */
  if (DEBUG_EVECS > 3) {
    checkorth(q, n, j);
  }

  /* periodic reorthogonalization */
  if (version == 1 || version == 2) {
    if ((j - lastpause) == interval || beta < 1000 * DOUBLE_EPSILON) {
      return (TRUE);
    }

    return (FALSE);
  }

  /* Run until orthogonality with first Lanczos vector deteriorates, then switch
     switch to periodic reorthog. */
  if (version == 3) {
    paigetol = 1.0e-3;
    if (*pausemode == 1) {
      paige_dot = fabs(dot(q[1], 1, n, q[j]));
      if ((paige_dot > paigetol && j > 1) || beta < 1000 * DOUBLE_EPSILON) {
        if (DEBUG_EVECS > 1) {
          printf("  Pausing on step %3d with Paige prod. = %g\n", j, paige_dot);
        }
        *pausemode = 2;
        return (TRUE);
      }

      return (FALSE);
    }
    if (*pausemode == 2) {
      if ((j - lastpause) == interval || beta < 1000 * DOUBLE_EPSILON) {
        return (TRUE);
      }

      return (FALSE);
    }
  }

  /* shouldn't ever get this far, but alint really wants a return value */
  return (FALSE);
}

int lanpause_float(int     j,         /* current step */
                   int     lastpause, /* when last paused */
                   int     interval,  /* interval between pauses */
                   float **q,         /* the Lanczos vectors */
                   int     n,         /* length of Lanczos vectors */
                   int *   pausemode, /* which pausing criterion to use */
                   int     version,   /* which version of sel. orth. we are using */
                   double  beta       /* current off-diagonal value */
)
{
  extern int    DEBUG_EVECS;    /* debugging level for eigen computation */
  extern double DOUBLE_EPSILON; /* machine precision */
  double        paige_dot;      /* q[j]^T q[1] */
  double        paigetol;       /* pause if paigedot > paigetol */
  double        dot_float();    /* standard dot product */
  void          checkorth_float();

  /* Check orthogonality of last Lanczos vector against previous ones */
  if (DEBUG_EVECS > 3) {
    checkorth_float(q, n, j);
  }

  /* periodic reorthogonalization */
  if (version == 1 || version == 2) {
    if ((j - lastpause) == interval || beta < 1000 * DOUBLE_EPSILON) {
      return (TRUE);
    }

    return (FALSE);
  }

  /* Run until orthogonality with first Lanczos vector deteriorates, then switch
     switch to periodic reorthog. */
  if (version == 3) {
    paigetol = 1.0e-3;
    if (*pausemode == 1) {
      paige_dot = fabs(dot_float(q[1], 1, n, q[j]));
      if ((paige_dot > paigetol && j > 1) || beta < 1000 * DOUBLE_EPSILON) {
        if (DEBUG_EVECS > 1) {
          printf("  Pausing on step %3d with Paige prod. = %g\n", j, paige_dot);
        }
        *pausemode = 2;
        return (TRUE);
      }

      return (FALSE);
    }
    if (*pausemode == 2) {
      if ((j - lastpause) == interval || beta < 1000 * DOUBLE_EPSILON) {
        return (TRUE);
      }

      return (FALSE);
    }
  }

  /* shouldn't ever get this far, but alint really wants a return value */
  return (FALSE);
}
