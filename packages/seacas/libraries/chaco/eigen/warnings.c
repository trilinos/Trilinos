/*
 * Copyright (c) 2005-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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

#include "defs.h" // for FALSE, TRUE
#include "structs.h"
#include <stdio.h> // for fprintf, NULL, FILE, stdout

/* Post various warnings about computation  */
void warnings(double *          workn,    /* work vector (1..n) */
              struct vtx_data **A,        /* graph */
              double **         y,        /* eigenvectors */
              int               n,        /* number of vtxs */
              double *          lambda,   /* ritz approximation to eigenvals of A */
              double *          vwsqrt,   /* square roots of vertex weights */
              double *          Ares,     /* how well Lanczos calc. eigpair lambda,y */
              double *          bound,    /* on ritz pair approximations to eig pairs of A */
              int *             index,    /* the Ritz index of an eigenpair */
              int               d,        /* problem dimension = number of eigvecs to find */
              int               j,        /* number of Lanczos iterations used */
              int               maxj,     /* maximum number of Lanczos iterations */
              double            Sres_max, /* Max value of Sres */
              double            eigtol,   /* tolerance on eigenvectors */
              double *          u,        /* Lanczos vector; here used as workspace */
              double            Anorm,    /* Gershgorin bound on eigenvalue */
              FILE *            out_file  /* output file */
)
{
  extern int    DEBUG_EVECS;              /* print debugging output? */
  extern int    WARNING_EVECS;            /* print warning messages? */
  extern double WARNING_ORTHTOL;          /* Warning: modest loss of orthogonality */
  extern double WARNING_MISTOL;           /* Warning: serious loss of orthogonality */
  extern double SRESTOL;                  /* limit on relative residual tol for evec of T */
  extern int    LANCZOS_CONVERGENCE_MODE; /* type of Lanczos convergence test */
  extern int    SRES_SWITCHES;            /* # switches to backup routine for computing s */
  int           warning1 = 0;             /* warning1 cond. (eigtol not achieved) true? */
  int           warning2 = 0;             /* warning2 cond. (premature orth. loss) true? */
  int           warning3 = 0;             /* warning3 cond. (suspected misconvergence) true? */
  int           i;                        /* loop index */
  int           hosed;                    /* flag for serious Lanczos problems */
  int           pass;                     /* which time through we are on */
  FILE *        outfile = NULL;           /* set to output file or stdout */
  double        checkeig();               /* calculate residual of eigenvector of A */
  void          doubleout_file();         /* print a double precision number */
  void          bail();                   /* our exit routine */

  hosed = FALSE;
  for (pass = 1; pass <= 2; pass++) {

    if (pass == 1) {
      outfile = stdout;
    }
    if (pass == 2) {
      if (out_file != NULL) {
        outfile = out_file;
      }
      else if (hosed) {
        bail(NULL, 1);
      }
      else {
        return;
      }
    }

    if (DEBUG_EVECS > 0 || WARNING_EVECS > 0) {
      if (LANCZOS_CONVERGENCE_MODE == 1) {
        fprintf(outfile, "Note about warnings: in partition convergence monitoring mode.\n");
      }
      for (i = 1; i <= d; i++) {
        Ares[i] = checkeig(workn, A, y[i], n, lambda[i], vwsqrt, u);
      }
    }

    if (DEBUG_EVECS > 0) {
      if (pass == 1) {
        fprintf(outfile, "Lanczos itns. = %d\n", j);
      }
      fprintf(outfile,
              "          lambda                Ares est.              Ares          index\n");
      for (i = 1; i <= d; i++) {
        fprintf(outfile, "%2d.", i);
        doubleout_file(outfile, lambda[i], 1);
        doubleout_file(outfile, bound[i], 1);
        doubleout_file(outfile, Ares[i], 1);
        fprintf(outfile, "   %3d\n", index[i]);
      }
      fprintf(outfile, "\n");
    }

    if (WARNING_EVECS > 0) {
      warning1 = FALSE;
      warning2 = FALSE;
      warning3 = FALSE;
      for (i = 1; i <= d; i++) {
        if (Ares[i] > eigtol) {
          warning1 = TRUE;
        }
        if (Ares[i] > WARNING_ORTHTOL * bound[i] && Ares[i] > .01 * eigtol) {
          warning2 = TRUE;
        }
        if (Ares[i] > WARNING_MISTOL * bound[i] && Ares[i] > .01 * eigtol) {
          warning3 = TRUE;
        }
      }
      if (j == maxj) {
        fprintf(outfile, "WARNING: Maximum number of Lanczos iterations reached.\n");
      }
      if (warning2 && !warning3) {
        fprintf(outfile, "WARNING: Minor loss of orthogonality (Ares/est. > %g).\n",
                WARNING_ORTHTOL);
      }
      if (warning3) {
        fprintf(outfile, "WARNING: Substantial loss of orthogonality (Ares/est. > %g).\n",
                WARNING_MISTOL);
      }
      if (warning1) {
        fprintf(outfile, "WARNING: Eigen pair tolerance (%g) not achieved.\n", eigtol);
      }
    }

    if (WARNING_EVECS > 1) {
      if (warning1 || warning2 || warning3) {
        if (DEBUG_EVECS <= 0) {
          fprintf(outfile,
                  "          lambda                Ares est.              Ares          index\n");
          for (i = 1; i <= d; i++) {
            fprintf(outfile, "%2d.", i);
            doubleout_file(outfile, lambda[i], 1);
            doubleout_file(outfile, bound[i], 1);
            doubleout_file(outfile, Ares[i], 1);
            fprintf(outfile, "   %3d\n", index[i]);
          }
        }
        /* otherwise gets printed above */
      }
    }

    if (warning1 || warning2 || warning3 || WARNING_EVECS > 2) {
      if (Sres_max > SRESTOL) {
        fprintf(outfile, "WARNING: Maximum eigen residual of T (%g) exceeds SRESTOL.\n", Sres_max);
      }
    }

    if (WARNING_EVECS > 2) {
      if (SRES_SWITCHES > 0) {
        fprintf(outfile, "WARNING: Switched routine for computing evec of T %d times.\n",
                SRES_SWITCHES);
        SRES_SWITCHES = 0;
      }
    }

    /* Put the best face on things ... */
    for (i = 1; i <= d; i++) {
      if (lambda[i] < 0 || lambda[i] > Anorm + eigtol) {
        hosed = TRUE;
      }
    }
    if (hosed) {
      fprintf(outfile, "ERROR: Sorry, out-of-bounds eigenvalue indicates serious breakdown.\n");
      fprintf(outfile, "       Try different parameters or another eigensolver.\n");
      if (pass == 2) {
        bail(NULL, 1);
      }
    }

  } /* Pass loop */
}
