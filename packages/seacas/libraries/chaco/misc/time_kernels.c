/*
 * Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "structs.h"
#include <math.h>  // for sqrt
#include <stdio.h> // for printf, NULL

static double checkvec();

/* Benchmark certain kernel operations */
void time_kernels(struct vtx_data **A,     /* matrix/graph being analyzed */
                  int               n,     /* number of rows/columns in matrix */
                  double           *vwsqrt /* square roots of vertex weights */
)
{
  extern int    DEBUG_PERTURB; /* debug flag for matrix perturbation */
  extern int    PERTURB;       /* randomly perturb to break symmetry? */
  extern int    NPERTURB;      /* number of edges to perturb */
  extern int    DEBUG_TRACE;   /* trace main execution path */
  extern double PERTURB_MAX;   /* maximum size of perturbation */
  int           i, beg, end;
  double       *dvec1, *dvec2, *dvec3;
  float        *svec1, *svec2, *svec3, *vwsqrt_float;
  double        norm_dvec, norm_svec;
  double        dot_dvec, dot_svec;
  double        time, time_dvec, time_svec;
  double        diff;
  double        factor, fac;
  float         factor_float, fac_float;
  int           loops;
  double        min_time, target_time;

  double *mkvec(int nl, int nh);
  float  *mkvec_float(int nl, int nh);
  void    frvec(double *v, int nl), frvec_float(float *v, int nl);
  void    vecran();
  double  ch_norm(double *vec, int beg, int end), dot(double *vec1, int beg, int end, double *vec2);
  double  norm_float(), dot_float(float *vec1, int beg, int end, float *vec2);
  double  seconds(void);
  void    scadd(), scadd_float(),
      update(double *vec1, int beg, int end, double *vec2, double fac, double *vec3),
      update_float(float *vec1, int beg, int end, float *vec2, float fac, float *vec3);
  void splarax(), splarax_float();
  void perturb_init(), perturb_clear();

  if (DEBUG_TRACE > 0) {
    printf("<Entering time_kernels>\n");
  }

  beg = 1;
  end = n;

  dvec1 = mkvec(beg, end);
  dvec2 = mkvec(beg, end);
  dvec3 = mkvec(beg - 1, end);
  svec1 = mkvec_float(beg, end);
  svec2 = mkvec_float(beg, end);
  svec3 = mkvec_float(beg - 1, end);

  if (vwsqrt == NULL) {
    vwsqrt_float = NULL;
  }
  else {
    vwsqrt_float = mkvec_float(beg - 1, end);
    for (i = beg - 1; i <= end; i++) {
      vwsqrt_float[i] = vwsqrt[i];
    }
  }

  vecran(dvec1, beg, end);
  vecran(dvec2, beg, end);
  vecran(dvec3, beg, end);
  for (i = beg; i <= end; i++) {
    svec1[i] = dvec1[i];
    svec2[i] = dvec2[i];
    svec3[i] = dvec3[i];
  }

  /* Set number of loops so that ch_norm(double *vec, int beg, int end) takes about one second. This
     should insulate against inaccurate timings on faster machines. */

  loops       = 1;
  time_dvec   = 0;
  min_time    = 0.5;
  target_time = 1.0;
  while (time_dvec < min_time) {
    time = seconds();
    for (i = loops; i; i--) {
      norm_dvec = ch_norm(dvec1, beg, end);
    }
    time_dvec = seconds() - time;
    if (time_dvec < min_time) {
      loops = 10 * loops;
    }
  }
  loops = (target_time / time_dvec) * loops;
  if (loops < 1) {
    loops = 1;
  }

  printf("                Kernel benchmarking\n");
  printf("Time (in seconds) for %d loops of each operation:\n\n", loops);

  printf("Routine      Double     Float      Discrepancy      Description\n");
  printf("-------      ------     -----      -----------      -----------\n");

  /* Norm operation */
  time = seconds();
  for (i = loops; i; i--) {
    norm_dvec = ch_norm(dvec1, beg, end);
  }
  time_dvec = seconds() - time;

  time = seconds();
  for (i = loops; i; i--) {
    norm_svec = norm_float(svec1, beg, end);
  }
  time_svec = seconds() - time;

  diff = norm_dvec - norm_svec;
  printf("norm        %6.2f    %6.2f    %14.5e", time_dvec, time_svec, diff);
  printf("      2 norm\n");

  /* Dot operation */
  time = seconds();
  for (i = loops; i; i--) {
    dot_dvec = dot(dvec1, beg, end, dvec2);
  }
  time_dvec = seconds() - time;

  time = seconds();
  for (i = loops; i; i--) {
    dot_svec = dot_float(svec1, beg, end, svec2);
  }
  time_svec = seconds() - time;

  diff = dot_dvec - dot_svec;
  printf("dot         %6.2f    %6.2f    %14.5e", time_dvec, time_svec, diff);
  printf("      scalar product\n");

  /* Scadd operation */
  factor       = 1.01;
  factor_float = factor;

  fac  = factor;
  time = seconds();
  for (i = loops; i; i--) {
    scadd(dvec1, beg, end, fac, dvec2);
    fac = -fac; /* to keep things in scale */
  }
  time_dvec = seconds() - time;

  fac_float = factor_float;
  time      = seconds();
  for (i = loops; i; i--) {
    scadd_float(svec1, beg, end, fac_float, svec2);
    fac_float = -fac_float; /* to keep things in scale */
  }
  time_svec = seconds() - time;

  diff = checkvec(dvec1, beg, end, svec1);
  printf("scadd       %6.2f    %6.2f    %14.5e", time_dvec, time_svec, diff);
  printf("      vec1 <- vec1 + alpha*vec2\n");

  /* Update operation */
  time = seconds();
  for (i = loops; i; i--) {
    update(dvec1, beg, end, dvec2, factor, dvec3);
  }
  time_dvec = seconds() - time;

  time = seconds();
  for (i = loops; i; i--) {
    update_float(svec1, beg, end, svec2, factor_float, svec3);
  }
  time_svec = seconds() - time;

  diff = checkvec(dvec1, beg, end, svec1);
  printf("update      %6.2f    %6.2f    %14.2g", time_dvec, time_svec, diff);
  printf("      vec1 <- vec2 + alpha*vec3\n");

  /* splarax operation */
  if (PERTURB) {
    if (NPERTURB > 0 && PERTURB_MAX > 0.0) {
      perturb_init(n);
      if (DEBUG_PERTURB > 0) {
        printf("Matrix being perturbed with scale %e\n", PERTURB_MAX);
      }
    }
    else if (DEBUG_PERTURB > 0) {
      printf("Matrix not being perturbed\n");
    }
  }

  time = seconds();
  for (i = loops; i; i--) {
    splarax(dvec1, A, n, dvec2, vwsqrt, dvec3);
  }
  time_dvec = seconds() - time;

  time = seconds();
  for (i = loops; i; i--) {
    splarax_float(svec1, A, n, svec2, vwsqrt_float, svec3);
  }

  time_svec = seconds() - time;

  diff = checkvec(dvec1, beg, end, svec1);
  printf("splarax     %6.2f    %6.2f    %14.5e", time_dvec, time_svec, diff);
  printf("      sparse matrix vector multiply\n");

  if (PERTURB && NPERTURB > 0 && PERTURB_MAX > 0.0) {
    perturb_clear();
  }
  printf("\n");

  /* Free memory */
  frvec(dvec1, 1);
  frvec(dvec2, 1);
  frvec(dvec3, 0);
  frvec_float(svec1, 1);
  frvec_float(svec2, 1);
  frvec_float(svec3, 0);
  if (vwsqrt_float != NULL) {
    frvec_float(vwsqrt_float, beg - 1);
  }
}

/* Compute norm of difference between a double and float vector. */
static double checkvec(double *dvec, int beg, int end, float *svec)
{
  double sum, diff;
  int    i;

  sum = 0;
  for (i = beg; i <= end; i++) {
    diff = dvec[i] - svec[i];
    sum += diff * diff;
  }
  return (sqrt(sum));
}
