/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __CH_INPUT_CONST_H
#define __CH_INPUT_CONST_H

#include <mpi.h>
#include "dr_const.h"
#include "dr_input_const.h"

extern int chaco_input_graph(FILE *, char *, int **, int **, int *, 
           float **, float **);

extern int chaco_input_geom(FILE *, char *, int, int *, float **, float **, 
                            float **);

extern int chaco_dist_graph(MPI_Comm, PARIO_INFO_PTR, 
                            int, int *, int *, int **, int **, 
                            float **, float **, 
                            int *, float **, float **, float **);

extern double read_val(FILE *, int *);
extern int read_int(FILE *, int *);

#ifndef TRUE
#define FALSE (0)
#define TRUE  (1)
#endif /* !TRUE */

#endif
