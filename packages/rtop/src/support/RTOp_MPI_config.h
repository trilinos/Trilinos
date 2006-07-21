/*
// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

/* If the macro RTOp_USE_MPI is defined, then these */
/* declarations will be MPI compatible.  If not then */
/* dummy MPI declarations will be used. */
/* */

#ifndef RTOP_MPI_CONFIG_H
#define RTOP_MPI_CONFIG_H

#include <RTOp_ConfigDefs.hpp>  /* This C++ header also has a proper C mode */

#ifdef HAVE_MPI
#define RTOp_USE_MPI  /* This macro is used in several places so we must keep it */
#endif

#ifdef RTOp_USE_MPI
#include "mpi.h"       /* Use real MPI declarations */
#else
/* #warning causes errors on Atlantis
#warning "Compiling in support for dummy MPI, real MPI will not be available!" */
#include "RTOp_mpi.h"  /* Use dummy MPI declarations */
#endif

/*
#ifdef __cplusplus
extern "C" {
#endif
*/

/** \file RTOp_config.h Platform dependent configuration options for RTOp.
 *
 * These typedefs and macros can be adjusted to the specific requirements of the platform.
 * For example, <tt>long double</tt> could be used instead of \c double if greater
 * precsion is needed.
 *
 * Also included are a few macros for MPI interoperability.  Also included in this default
 * header file is is \c RTOp_mpi.h which contains dummy MPI declarations (which are defined
 * in \c RTOp_mpi.c) for a subset of the MPI functions that are correct for the serial
 * case.
 */
/*@{ */

typedef MPI_Datatype            RTOp_Datatype;   /*< Compatible with MPI_Datatype? */
#define RTOpMPI_VALUE_TYPE      MPI_DOUBLE       /*< (MPI only) Compatible with fortran DOUBLE PRECISION? */
#define RTOpMPI_INDEX_TYPE      MPI_INT          /*< (MPI only) Compatible with fortran INTEGER? */
#define RTOpMPI_CHAR_TYPE       MPI_CHAR         /*< (MPI only) Compatible with fortran CHARACTER? */

/* The maxinum number of characters in a name of a reduction/transformation operator class */
#define RTOp_MAX_REDUCT_TRANS_OP_CLASS_NAME  50

/*@} */

/*
#ifdef __cplusplus
}
#endif
*/

#endif /* RTOP_MPI_CONFIG_H */
