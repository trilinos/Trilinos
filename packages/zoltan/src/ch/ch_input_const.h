/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

#ifndef __CH_INPUT_CONST_H
#define __CH_INPUT_CONST_H

#ifndef lint
static char *cvs_ch_input_const_h = "$Id$";
#endif

extern int chaco_input_graph(FILE *, char *, int **, int **, int *, int **, 
                             float **);

extern int chaco_input_geom(FILE *, char *, int, int *, float **, float **, 
                            float **);

extern int DEBUG_TRACE;
extern int CHECK_INPUT;
extern int DEBUG_INPUT;
#endif
