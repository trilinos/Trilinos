/*******************************************************************************
 * Copyright 1999, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


/*
 * This code is a wrapper layer around AZTEC.  It takes the real and imaginary
 * parts of a complex valued linear system and forms an equivalent real system,
 * calls AZTEC to solve the real system and returns the solution.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "az_aztec.h"
#include "az_ifpack.h"

void AZ_ifpack_iterate(double *x, double *b,
			int *options, double *params,
			double *status, int *proc_config,
			AZ_MATRIX *Amat )

/*******************************************************************************

  Author:          Mike Heroux, SNL, 9222
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  x:               On input, contains the initial guess, On output 
                   contains the solution to 
                   the linear system.

  b:               Right hand side of linear system.

  options:         Determines specific solution method and other parameters.

  params:          Drop tolerance and convergence tolerance info.

  status:          On output, indicates termination status:
                    0:  terminated normally.
                   -1:  maximum number of iterations taken without achieving
                        convergence.
                   -2:  Breakdown. The algorithm can not proceed due to
                        numerical difficulties (usually a divide by zero).
                   -3:  Internal residual differs from the computed residual due
                        to a significant loss of precision.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

 Amat:             The matrix operator, stored as an AZ_MATRIX structure.

  Internal Parameter list:
  ========================

 x:                Komplex version of initial guess and solution.
 b:                Komplex version of RHS.
 Prec:             Preconditioner stored as an AZ_PRECOND structure.

  Overview
  ========

*******************************************************************************/


{

AZ_PRECOND *Prec;    /* Structure representing entire preconditioner.        */
                    /*                                                      */

   AZ_ifpack_prec_create (x, b, options,  params, proc_config, Amat, &Prec);

   /* solve linear system using Aztec. */

   AZ_iterate(x, b, options, params, status, proc_config, Amat, Prec, NULL);

   AZ_ifpack_prec_destroy (options,  params, proc_config, Amat, Prec);

/* AZ_ifpack_iterate*/
}
