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

/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1986 Thomas L. Quarles
**********/

/*
 */

#ifndef IFERRMSGS
#define IFERRMSGS


    /* common error message descriptions */

#define E_PAUSE -1      /* pausing on demand */

#define OK 0

#define E_PANIC 1       /* vague internal error for "can't get here" cases */
#define E_EXISTS 2      /* warning/error - attempt to create duplicate */
                        /* instance or model. Old one reused instead */
#define E_NODEV 3       /* attempt to modify a non-existant instance */
#define E_NOMOD 4       /* attempt to modify a non-existant model */
#define E_NOANAL 5      /* attempt to modify a non-existant analysis */
#define E_NOTERM 6      /* attempt to bind to a non-existant terminal */
#define E_BADPARM 7     /* attempt to specify a non-existant parameter */
#define E_NOMEM 8       /* insufficient memory available - VERY FATAL */
#define E_NODECON 9     /* warning/error - node already connected, old */
                        /* connection replaced */
#define E_UNSUPP 10     /* the specified operation is unsupported by the */
                        /* simulator */
#define E_PARMVAL 11    /* the parameter value specified is illegal */
#define E_NOTEMPTY 12   /* deleted still referenced item. */
#define E_NOCHANGE 13   /* simulator can't tolerate any more topology changes */
#define E_NOTFOUND 14   /* simulator can't find something it was looking for */
#define E_BAD_DOMAIN 15 /* output interface begin/end domain calls mismatched */
#define E_NOCIRCUIT 16  /* no circuit found - likely premature .end */
#define E_NOTGIVEN 17   /* parameter value was not given explicitly */
#define E_MULTIPARM 18  /* Multiple specification of parameter */
#define E_BADEXP   19   /* Too many math error in math expression (for asrc) */
#define E_DEV_OVERFLOW 20
#define E_NOGROUND 21   /* no path to ground for node(s) */
#define E_SUBCIRCUIT 22 /* Error in subcircuit expansion or parameter handling */
#define E_NOGLOBAL 23   /* global variable not found on rerun command */


#define E_PRIVATE 100   /* messages above this number are private to */
                        /* the simulator and MUST be accompanied by */
                        /* a proper setting of errMsg */
                        /* this constant should be added to all such messages */
                        /* to ensure error free operation if it must be */
                        /* changed in the future */

extern char *errMsg;    /* descriptive message about what went wrong */
                        /* MUST be malloc()'d - front end will free() */
                        /* this should be a detailed message,and is assumed */
                        /* malloc()'d so that you will feel free to add */
                        /* lots of descriptive information with sprintf*/

extern char *errRtn;    /* name of the routine declaring error */
                        /* should not be malloc()'d, will not be free()'d */
                        /* This should be a simple constant in your routine */
                        /* and thus can be set correctly even if we run out */
                        /* of memory */

#endif /*IFERRMSGS*/
