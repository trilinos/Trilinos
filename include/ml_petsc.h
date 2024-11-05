/*!
 * \file ml_petsc.h
 *
 * \brief ML wrappers for PETSc data stuctures.
 *
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/*#############################################################################
# CVS File Information
#    Current revision: $Revision$
#    Last modified:    $Date$
#    Modified by:      $Author$
#############################################################################*/

#ifndef ML_PETSC_H
#define ML_PETSC_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#ifdef HAVE_ML_PETSC

#include "petscksp.h"

/*wrap PETSc's pointers to preconditioner and solver structures*/
typedef PC ML_PetscPC;
typedef KSP ML_PetscKSP;

#endif /*ifdef HAVE_ML_PETSC*/

#endif /* define ML_PETSC_H */
